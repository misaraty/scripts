import os, sys, re, random, shutil, subprocess, glob

os.chdir(os.path.split(os.path.realpath(__file__))[0])

from pathlib import Path
from functools import lru_cache
import math
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from rdkit import Chem
import selfies as sf

os.environ["MPLBACKEND"] = "Agg"
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

DATA_CSV = "D.csv"
BASE_DIR = Path("runs/mol_chemprop_multi")
CP_DIR = BASE_DIR / "models_chemprop"
OUT_DIR = BASE_DIR / "outputs"
SPLIT_DIR = CP_DIR / "cv_splits"
for p in [CP_DIR, OUT_DIR, SPLIT_DIR]:
    p.mkdir(parents=True, exist_ok=True)

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(SEED)
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

TARGET_D = 20

DM_DUPLICATE_MODE = "mean"

EXCLUDED_ELEMENTS = []

SIZE_FILTER_QUANTILES = (0.01, 0.99)

NUM_GEN = 2048
ENRICH_ROUNDS = 5

EPOCHS_GEN = 30
LATENT_DIM = 384
TEMPERATURE = 0.9
TOP_K = 100
TOP_P = 0.9
REPEAT_PENALTY = 1.15

MAX_LEN = 0
BATCH_SIZE_GEN = 0
SAMPLE_BATCH_SIZE = 0

EPOCHS_VAE = 30
LATENT_DIM_VAE = 128
BETA_VAE = 1.0
KL_ANNEAL = True

CHEMPROP_FOLDS = 5
CHEMPROP_EPOCHS = 200
CHEMPROP_BATCH_SIZE = 2

PAD_ID = 0

DM_CANONICAL_SET = set()
DM_HEAVY_MIN = None
DM_HEAVY_MAX = None
DM_ELEMENTS = set()
RESOLVED_DIRECTION = "max"
SELFIES_LENGTH_STATS = {}


def _bin(cmd):
    p = shutil.which(cmd)
    if p:
        return p
    if os.name == "nt":
        cand = os.path.join(sys.prefix, "Scripts", cmd + ".exe")
    else:
        cand = os.path.join(sys.prefix, "bin", cmd)
    if not os.path.exists(cand):
        raise FileNotFoundError(f"Cannot find {cmd} at {cand}")
    return cand


CHEMPROP_TRAIN = _bin("chemprop_train")
CHEMPROP_PREDICT = _bin("chemprop_predict")


def run(cmd, cwd=None):
    print(">>", " ".join(map(str, cmd)))
    subprocess.check_call(cmd, cwd=cwd)


def _as_tensor(batch):
    if isinstance(batch, (list, tuple)):
        batch = batch[0]
    if isinstance(batch, np.ndarray):
        batch = torch.from_numpy(batch)
    return batch


def normalize_element_symbol(symbol):
    s = str(symbol).strip()
    if not s:
        return None
    s = s[0].upper() + s[1:].lower()
    try:
        if Chem.GetPeriodicTable().GetAtomicNumber(s) <= 0:
            return None
    except Exception:
        return None
    return s


EXCLUDED_ELEMENT_SET = {
    x for x in (normalize_element_symbol(v) for v in EXCLUDED_ELEMENTS) if x is not None
}


@lru_cache(maxsize=300000)
def _canonicalize_smiles_cached(s):
    try:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        return None


def canonicalize_smiles(smiles):
    if smiles is None:
        return None
    s = str(smiles).strip()
    if not s or s.lower() == "nan":
        return None
    return _canonicalize_smiles_cached(s)


@lru_cache(maxsize=300000)
def _mol_properties_cached(cs):
    try:
        mol = Chem.MolFromSmiles(cs)
        if mol is None:
            return None
        neutral = Chem.rdmolops.GetFormalCharge(mol) == 0
        heavy = int(mol.GetNumHeavyAtoms())
        elements = tuple(
            sorted({a.GetSymbol() for a in mol.GetAtoms() if a.GetAtomicNum() > 1})
        )
        return neutral, heavy, elements
    except Exception:
        return None


def is_neutral_canonical(cs):
    p = _mol_properties_cached(cs)
    return p is not None and p[0]


def heavy_atom_count_canonical(cs):
    p = _mol_properties_cached(cs)
    return None if p is None else p[1]


def elements_canonical(cs):
    p = _mol_properties_cached(cs)
    return set() if p is None else set(p[2])


def contains_excluded_elements_canonical(cs):
    if not EXCLUDED_ELEMENT_SET:
        return False
    return bool(elements_canonical(cs) & EXCLUDED_ELEMENT_SET)


def within_dm_size_range_canonical(cs):
    if DM_HEAVY_MIN is None or DM_HEAVY_MAX is None:
        return True
    n = heavy_atom_count_canonical(cs)
    return n is not None and DM_HEAVY_MIN <= n <= DM_HEAVY_MAX


def resolve_target_direction(values):
    global RESOLVED_DIRECTION
    med = float(np.median(np.asarray(values, dtype=float)))
    RESOLVED_DIRECTION = "max" if float(TARGET_D) >= med else "min"
    print(
        f"[TARGET] TARGET_D={TARGET_D}, training median={med:.6f}, "
        f"direction={RESOLVED_DIRECTION}"
    )
    return RESOLVED_DIRECTION


def _aggregate_duplicate_group(g, mode):
    if mode == "first":
        return float(g.iloc[0])
    if mode == "mean":
        return float(g.mean())
    if mode == "median":
        return float(g.median())
    if mode == "max":
        return float(g.max())
    if mode == "min":
        return float(g.min())
    raise ValueError("DM_DUPLICATE_MODE must be first/mean/median/max/min")


def _align_up(value, multiple=8):
    value = int(value)
    return int(math.ceil(value / float(multiple)) * multiple)


def _auto_sequence_batch_sizes(max_len):
    if DEVICE != "cuda":
        train_bs, sample_bs = 16, 32
    else:
        try:
            gb = torch.cuda.get_device_properties(0).total_memory / (1024**3)
        except Exception:
            gb = 12.0
        if gb >= 28:
            train_bs, sample_bs = 64, 256
        elif gb >= 20:
            train_bs, sample_bs = 48, 192
        elif gb >= 14:
            train_bs, sample_bs = 32, 128
        elif gb >= 10:
            train_bs, sample_bs = 16, 64
        else:
            train_bs, sample_bs = 8, 32
        if max_len > 256:
            train_bs = max(8, train_bs // 2)
            sample_bs = max(32, sample_bs // 2)
        if max_len > 384:
            train_bs = max(4, train_bs // 2)
            sample_bs = max(16, sample_bs // 2)
    return int(train_bs), int(sample_bs)


def _auto_seed_count(n_candidates):
    n_candidates = int(n_candidates)
    if n_candidates <= 0:
        return 0
    k = int(math.ceil(n_candidates * 0.10))
    k = max(32, min(256, k))
    return min(n_candidates, k)


def initialize_sequence_config(df):
    global MAX_LEN, BATCH_SIZE_GEN, SAMPLE_BATCH_SIZE, SELFIES_LENGTH_STATS

    lengths = []
    for cs in df["SMILES"].tolist():
        try:
            ss = sf.encoder(cs)
            lengths.append(len(list(sf.split_selfies(ss))))
        except Exception:
            pass
    if not lengths:
        raise RuntimeError("Cannot determine SELFIES lengths from Dm.csv")

    ls = pd.Series(lengths, dtype=float)
    max_content = int(max(lengths))
    MAX_LEN = _align_up(max_content + 2, 8)
    BATCH_SIZE_GEN, SAMPLE_BATCH_SIZE = _auto_sequence_batch_sizes(MAX_LEN)
    SELFIES_LENGTH_STATS = {
        "p50": float(ls.quantile(0.50)),
        "p95": float(ls.quantile(0.95)),
        "p99": float(ls.quantile(0.99)),
        "max": max_content,
    }

    print(
        f"[AUTO SELFIES] P50={SELFIES_LENGTH_STATS['p50']:.1f} "
        f"P95={SELFIES_LENGTH_STATS['p95']:.1f} "
        f"P99={SELFIES_LENGTH_STATS['p99']:.1f} max={max_content} MAX_LEN={MAX_LEN}"
    )
    print(
        f"[AUTO HARDWARE] BATCH_SIZE_GEN={BATCH_SIZE_GEN} "
        f"SAMPLE_BATCH_SIZE={SAMPLE_BATCH_SIZE} device={DEVICE}"
    )

    cfg_path = OUT_DIR / "auto_config_sequence.txt"
    with open(cfg_path, "w", encoding="utf-8") as f:
        f.write(f"unique_molecules={len(df)}\n")
        f.write(f"elements={','.join(sorted(DM_ELEMENTS))}\n")
        f.write(f"selfies_p50={SELFIES_LENGTH_STATS['p50']:.6f}\n")
        f.write(f"selfies_p95={SELFIES_LENGTH_STATS['p95']:.6f}\n")
        f.write(f"selfies_p99={SELFIES_LENGTH_STATS['p99']:.6f}\n")
        f.write(f"selfies_max_content={max_content}\n")
        f.write(f"MAX_LEN={MAX_LEN}\n")
        f.write(f"BATCH_SIZE_GEN={BATCH_SIZE_GEN}\n")
        f.write(f"SAMPLE_BATCH_SIZE={SAMPLE_BATCH_SIZE}\n")
        f.write(f"heavy_atom_min={DM_HEAVY_MIN}\n")
        f.write(f"heavy_atom_max={DM_HEAVY_MAX}\n")
        f.write(f"TARGET_D={TARGET_D}\n")
        f.write(f"direction={RESOLVED_DIRECTION}\n")


def load_dm_dataframe(path):
    global DM_CANONICAL_SET, DM_HEAVY_MIN, DM_HEAVY_MAX, DM_ELEMENTS

    df = pd.read_csv(path)
    smiles_col = None
    d_col = None
    for c in df.columns:
        if smiles_col is None and re.search("smiles", str(c), re.I):
            smiles_col = c
        if d_col is None and re.search(
            r"(^|[^A-Za-z])D([^A-Za-z]|$)|velocity|label|target", str(c), re.I
        ):
            d_col = c
    if smiles_col is None or d_col is None:
        raise ValueError("Cannot find SMILES/D column.")

    raw = (
        df[[smiles_col, d_col]]
        .rename(columns={smiles_col: "SMILES", d_col: "D"})
        .dropna()
        .reset_index(drop=True)
    )
    raw["D"] = pd.to_numeric(raw["D"], errors="coerce")
    raw = raw.dropna(subset=["D"]).reset_index(drop=True)

    rows = []
    invalid = []
    for i, row in raw.iterrows():
        cs = canonicalize_smiles(row["SMILES"])
        if cs is None:
            invalid.append((i, str(row["SMILES"])))
            continue
        rows.append((cs, float(row["D"])))

    if invalid:
        bad_path = OUT_DIR / "invalid_smiles_in_Dm.csv"
        pd.DataFrame(invalid, columns=["row_index", "SMILES"]).to_csv(
            bad_path, index=False
        )
        print(f"[WARN] invalid SMILES filtered: {len(invalid)} -> {bad_path}")

    valid = pd.DataFrame(rows, columns=["SMILES", "D"])
    before = len(valid)

    grouped = []
    for cs, g in valid.groupby("SMILES", sort=False):
        grouped.append((cs, _aggregate_duplicate_group(g["D"], DM_DUPLICATE_MODE)))
    out = pd.DataFrame(grouped, columns=["SMILES", "D"]).reset_index(drop=True)

    removed = before - len(out)
    DM_CANONICAL_SET = set(out["SMILES"])

    heavy = []
    elems = set()
    for cs in out["SMILES"]:
        n = heavy_atom_count_canonical(cs)
        if n is not None:
            heavy.append(n)
        elems.update(elements_canonical(cs))

    DM_ELEMENTS = elems

    if SIZE_FILTER_QUANTILES is None:
        DM_HEAVY_MIN = None
        DM_HEAVY_MAX = None
    else:
        qlo, qhi = SIZE_FILTER_QUANTILES
        if not (0 <= qlo < qhi <= 1):
            raise ValueError("SIZE_FILTER_QUANTILES must satisfy 0 <= low < high <= 1")
        if heavy:
            DM_HEAVY_MIN = max(1, int(math.floor(np.quantile(heavy, qlo))))
            DM_HEAVY_MAX = int(math.ceil(np.quantile(heavy, qhi)))

    resolve_target_direction(out["D"].values)

    dedup_path = OUT_DIR / "Dm_deduplicated.csv"
    out.to_csv(dedup_path, index=False)

    print(
        f"[Dm] input={len(raw)} invalid={len(invalid)} duplicate_mode={DM_DUPLICATE_MODE} "
        f"removed_internal={removed} kept={len(out)}"
    )
    print(f"[Dm] elements={sorted(DM_ELEMENTS)}")
    if DM_HEAVY_MIN is not None:
        print(
            f"[SIZE] heavy-atom range learned from Dm "
            f"{SIZE_FILTER_QUANTILES}: {DM_HEAVY_MIN}..{DM_HEAVY_MAX}"
        )
    if EXCLUDED_ELEMENT_SET:
        print(f"[FILTER] excluded generated elements={sorted(EXCLUDED_ELEMENT_SET)}")
    else:
        print("[FILTER] EXCLUDED_ELEMENTS is empty: element filtering disabled")

    initialize_sequence_config(out)
    return out


def finalize_generated(smiles_list, tag="Generated", exclude_dm=True):
    kept = []
    seen = set()
    stats = {
        "invalid": 0,
        "nonneutral": 0,
        "excluded_element": 0,
        "size": 0,
        "internal_duplicate": 0,
        "dm_duplicate": 0,
    }

    dm_set = DM_CANONICAL_SET if exclude_dm else set()

    for s in smiles_list:
        cs = canonicalize_smiles(s)
        if cs is None:
            stats["invalid"] += 1
            continue
        if not is_neutral_canonical(cs):
            stats["nonneutral"] += 1
            continue
        if contains_excluded_elements_canonical(cs):
            stats["excluded_element"] += 1
            continue
        if not within_dm_size_range_canonical(cs):
            stats["size"] += 1
            continue
        if cs in seen:
            stats["internal_duplicate"] += 1
            continue
        if cs in dm_set:
            stats["dm_duplicate"] += 1
            continue
        seen.add(cs)
        kept.append(cs)

    print(
        f"[FILTER:{tag}] input={len(smiles_list)} kept={len(kept)} "
        f"invalid={stats['invalid']} nonneutral={stats['nonneutral']} "
        f"element={stats['excluded_element']} size={stats['size']} "
        f"dup={stats['internal_duplicate']} vsDm={stats['dm_duplicate']}"
    )
    return kept


def simple_kfold_indices(n_splits, n_samples, seed=SEED):
    rng = np.random.RandomState(seed)
    idx = np.arange(n_samples)
    rng.shuffle(idx)
    folds = []
    sizes = [n_samples // n_splits] * n_splits
    for i in range(n_samples % n_splits):
        sizes[i] += 1
    start = 0
    for size in sizes:
        folds.append(idx[start : start + size])
        start += size
    return folds


def make_fixed_cv_splits(df, n_folds=CHEMPROP_FOLDS):
    folds = simple_kfold_indices(n_folds, len(df), seed=SEED)
    assignment = np.full(len(df), -1, dtype=int)
    for k, ids in enumerate(folds):
        assignment[ids] = k

    manifest = df.copy()
    manifest["test_fold"] = assignment
    manifest.to_csv(SPLIT_DIR / "cv_assignments.csv", index=False)

    split_meta = []
    all_idx = np.arange(len(df))
    for k in range(n_folds):
        test_idx = np.asarray(folds[k], dtype=int)
        val_idx = np.asarray(folds[(k + 1) % n_folds], dtype=int)
        holdout = set(test_idx.tolist()) | set(val_idx.tolist())
        train_idx = np.asarray([i for i in all_idx if i not in holdout], dtype=int)

        train_df = df.iloc[train_idx].reset_index(drop=True)
        val_df = df.iloc[val_idx].reset_index(drop=True)
        test_df = df.iloc[test_idx].reset_index(drop=True)

        train_path = SPLIT_DIR / f"fold_{k}_train.csv"
        val_path = SPLIT_DIR / f"fold_{k}_val.csv"
        test_path = SPLIT_DIR / f"fold_{k}_test.csv"

        train_df.rename(columns={"SMILES": "smiles"}).to_csv(train_path, index=False)
        val_df.rename(columns={"SMILES": "smiles"}).to_csv(val_path, index=False)
        test_df.rename(columns={"SMILES": "smiles"}).to_csv(test_path, index=False)

        split_meta.append(
            {
                "fold": k,
                "train_path": train_path,
                "val_path": val_path,
                "test_path": test_path,
                "test_df": test_df,
            }
        )
    return split_meta


def find_checkpoint(save_dir):
    patterns = [
        save_dir / "model_0" / "model.pt",
        save_dir / "model.pt",
    ]
    for p in patterns:
        if p.exists():
            return str(p)
    cands = glob.glob(str((save_dir / "**" / "*.pt").as_posix()), recursive=True)
    return cands[0] if cands else None


def mae(y_true, y_pred):
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    return float(np.mean(np.abs(y_true - y_pred)))


def rmse(y_true, y_pred):
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    return float(np.sqrt(np.mean((y_true - y_pred) ** 2)))


def r2(y_true, y_pred):
    y_true = np.asarray(y_true, float)
    y_pred = np.asarray(y_pred, float)
    sse = np.sum((y_true - y_pred) ** 2)
    sst = np.sum((y_true - np.mean(y_true)) ** 2)
    return float(1.0 - sse / sst) if sst > 0 else float("nan")


def pick_pred_column(df):
    lower_cols = {str(c).lower(): c for c in df.columns}
    for key in ["d_pred", "pred", "prediction", "y_pred", "d"]:
        if key in lower_cols:
            return lower_cols[key]
    num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    return num_cols[0] if num_cols else None


def chemprop_kfold_train(df):
    split_meta = make_fixed_cv_splits(df, CHEMPROP_FOLDS)

    for info in split_meta:
        k = info["fold"]
        fold_dir = CP_DIR / f"fold_{k}"
        if fold_dir.exists():
            shutil.rmtree(fold_dir)
        fold_dir.mkdir(parents=True, exist_ok=True)

        cmd = [
            CHEMPROP_TRAIN,
            "--data_path",
            str(info["train_path"]),
            "--separate_val_path",
            str(info["val_path"]),
            "--separate_test_path",
            str(info["test_path"]),
            "--dataset_type",
            "regression",
            "--target_columns",
            "D",
            "--save_dir",
            str(fold_dir),
            "--epochs",
            str(CHEMPROP_EPOCHS),
            "--batch_size",
            str(CHEMPROP_BATCH_SIZE),
            "--metric",
            "rmse",
            "--seed",
            str(SEED + k),
            "--quiet",
        ]
        run(cmd)

    kfold_eval_and_artifacts(split_meta, OUT_DIR)
    return CP_DIR


def kfold_eval_and_artifacts(split_meta, outs_dir):
    cols = {}
    per_fold = []
    colors = [plt.get_cmap("tab10")(i) for i in range(len(split_meta))]
    fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
    xlim = [float("inf"), -float("inf")]
    ylim = [float("inf"), -float("inf")]

    for info in split_meta:
        k = info["fold"]
        test_df = info["test_df"]
        test_input = SPLIT_DIR / f"fold_{k}_eval_input.csv"
        test_df[["SMILES"]].rename(columns={"SMILES": "smiles"}).to_csv(
            test_input, index=False
        )

        ck = find_checkpoint(CP_DIR / f"fold_{k}")
        if not ck:
            raise FileNotFoundError(f"No checkpoint for fold {k}")

        pred_csv = outs_dir / f"preds_test_fold{k}.csv"
        run(
            [
                CHEMPROP_PREDICT,
                "--test_path",
                str(test_input),
                "--checkpoint_path",
                ck,
                "--preds_path",
                str(pred_csv),
            ]
        )

        pr = pd.read_csv(pred_csv)
        pcol = pick_pred_column(pr)
        if pcol is None:
            raise ValueError(f"No numeric prediction column found in {pred_csv}")

        y_true = test_df["D"].astype(float).values
        y_pred = pr[pcol].astype(float).values
        if len(y_true) != len(y_pred):
            raise ValueError(f"Fold {k}: prediction length mismatch")

        mma = mae(y_true, y_pred)
        mmr = rmse(y_true, y_pred)
        rr = r2(y_true, y_pred)
        per_fold.append((k, mma, mmr, rr, len(y_true)))

        cols[f"fold{k}_true"] = y_true.tolist()
        cols[f"fold{k}_pred"] = y_pred.tolist()

        ax.scatter(
            y_true,
            y_pred,
            s=16,
            alpha=0.7,
            edgecolors="none",
            label=f"Fold {k}",
            color=colors[k],
        )
        xlim = [min(xlim[0], np.min(y_true)), max(xlim[1], np.max(y_true))]
        ylim = [min(ylim[0], np.min(y_pred)), max(ylim[1], np.max(y_pred))]

    low = min(xlim[0], ylim[0])
    high = max(xlim[1], ylim[1])
    ax.plot([low, high], [low, high], ls="--", lw=1, color="black")
    ax.set_xlabel("True D")
    ax.set_ylabel("Predicted D")
    ax.set_title("Fixed 5-fold held-out predictions")
    ax.legend(frameon=False, markerscale=1.2, fontsize=8)
    ax.set_xlim(low, high)
    ax.set_ylim(low, high)
    fig.tight_layout()

    fig_path = outs_dir / "kfold_val_pred_300dpi.jpg"
    fig.savefig(fig_path, dpi=300)
    plt.close(fig)

    max_len = max(len(v) for v in cols.values())
    dat = pd.DataFrame({k: v + [np.nan] * (max_len - len(v)) for k, v in cols.items()})
    dat_path = outs_dir / "kfold_val_preds.dat"
    dat.to_csv(dat_path, sep="\t", index=False)

    maes = [p[1] for p in per_fold]
    rmses = [p[2] for p in per_fold]
    r2s = [p[3] for p in per_fold]
    mlog = outs_dir / "verbose.log"
    with open(mlog, "w", encoding="utf-8") as f:
        f.write("Fixed 5-fold held-out evaluation\n")
        for k, mma, mmr, rr, n in per_fold:
            f.write(
                f"\tFold {k}: n={n}, MAE={mma:.6f}, " f"RMSE={mmr:.6f}, R2={rr:.6f}\n"
            )
        f.write(f"Overall MAE={np.mean(maes):.6f} +/- {np.std(maes, ddof=1):.6f}\n")
        f.write(f"Overall RMSE={np.mean(rmses):.6f} +/- {np.std(rmses, ddof=1):.6f}\n")
        f.write(f"Overall R2={np.nanmean(r2s):.6f} +/- {np.nanstd(r2s, ddof=1):.6f}\n")

    print(f"[KFold] metrics -> {mlog}")
    print(f"[KFold] data -> {dat_path}")
    print(f"[KFold] figure -> {fig_path}")


def chemprop_predict_paths(ckpt_dir, gen_csv, tag):
    out_csv = OUT_DIR / f"generated_with_pred_{tag}.csv"
    hits_csv = OUT_DIR / f"generated_hits_{tag}.csv"

    if not os.path.exists(gen_csv) or os.path.getsize(gen_csv) == 0:
        return None, None
    try:
        _input_df = pd.read_csv(gen_csv)
        if len(_input_df) == 0:
            print(f"[WARN:{tag}] empty generated CSV, prediction skipped")
            return None, None
    except Exception:
        return None, None

    run(
        [
            CHEMPROP_PREDICT,
            "--test_path",
            str(gen_csv),
            "--checkpoint_dir",
            str(ckpt_dir),
            "--preds_path",
            str(out_csv),
        ]
    )

    if not out_csv.exists() or out_csv.stat().st_size == 0:
        return None, None

    df = pd.read_csv(out_csv)
    if "smiles" not in df.columns:
        if "SMILES" in df.columns:
            df = df.rename(columns={"SMILES": "smiles"})
        elif "smiles_0" in df.columns:
            df = df.rename(columns={"smiles_0": "smiles"})

    pcol = pick_pred_column(df)
    if pcol is None:
        return str(out_csv), None

    df["D_pred"] = pd.to_numeric(df[pcol], errors="coerce")
    df = df.dropna(subset=["D_pred"]).copy()
    ascending = RESOLVED_DIRECTION == "min"
    df = df.sort_values("D_pred", ascending=ascending).reset_index(drop=True)

    if RESOLVED_DIRECTION == "max":
        hits = df[df["D_pred"] >= float(TARGET_D)].copy()
    else:
        hits = df[df["D_pred"] <= float(TARGET_D)].copy()

    df.to_csv(out_csv, index=False)
    hits.to_csv(hits_csv, index=False)

    print(
        f"[OK:{tag}] pred={len(df)} hits={len(hits)} "
        f"target={TARGET_D} direction={RESOLVED_DIRECTION}"
    )
    return str(out_csv), str(hits_csv)


def extract_top_seeds(pred_csv):
    if pred_csv is None or not os.path.exists(pred_csv):
        return []
    df = pd.read_csv(pred_csv)
    scol = (
        "smiles"
        if "smiles" in df.columns
        else ("SMILES" if "SMILES" in df.columns else None)
    )
    if scol is None:
        return []
    pcol = "D_pred" if "D_pred" in df.columns else pick_pred_column(df)
    if pcol is None:
        return []
    ascending = RESOLVED_DIRECTION == "min"
    df = df.sort_values(pcol, ascending=ascending)
    topk = _auto_seed_count(len(df))
    print(f"[AUTO SEEDS] candidates={len(df)} selected={topk}")
    return df[scol].dropna().astype(str).head(topk).tolist()


def enrich_by_predictions(ckpt, gen_csv, tag):
    preds, hits = chemprop_predict_paths(ckpt, gen_csv, tag)
    seeds = extract_top_seeds(preds) if preds else []
    return preds, seeds


def merge_and_dedup(csv_list, out_name):
    all_sm = []
    for p in csv_list:
        if not p or not os.path.exists(p):
            continue
        try:
            df = pd.read_csv(p)
            col = (
                "smiles"
                if "smiles" in df.columns
                else ("SMILES" if "SMILES" in df.columns else None)
            )
            if col is not None:
                all_sm.extend(df[col].dropna().astype(str).tolist())
        except Exception as e:
            print(f"[WARN] merge failed for {p}: {e}")

    uniq = finalize_generated(all_sm, tag=f"MERGE:{out_name}", exclude_dm=True)
    out = OUT_DIR / out_name
    pd.DataFrame({"smiles": uniq}).to_csv(out, index=False)
    return out


def selfies_tokenize(smiles_list):
    selfies_list = []
    bad = 0
    for s in smiles_list:
        try:
            selfies_list.append(sf.encoder(s))
        except Exception:
            bad += 1
    if bad:
        print(f"[WARN] SELFIES encoding failed: {bad}")
    if not selfies_list:
        raise ValueError("No valid SELFIES after filtering/encoding.")

    tokens = []
    for s in selfies_list:
        tokens.extend(sf.split_selfies(s))

    vocab = sorted(set(tokens))
    stoi = {t: i + 1 for i, t in enumerate(vocab)}
    itos = {i: t for t, i in stoi.items()}
    bos_id = len(vocab) + 1
    eos_id = len(vocab) + 2
    itos[bos_id] = "<BOS>"
    itos[eos_id] = "<EOS>"
    return stoi, itos, vocab, bos_id, eos_id


def selfies_to_tensor(smiles_list, stoi, max_len, bos_id, eos_id):
    seqs = []
    bad = 0
    for sm in smiles_list:
        try:
            s = sf.encoder(sm)
            toks = list(sf.split_selfies(s))
            ids = [stoi[t] for t in toks if t in stoi]
        except Exception:
            bad += 1
            continue

        ids = ids[: max_len - 1]
        ids_inp = [bos_id] + ids
        ids_tgt = ids + [eos_id]

        ids_inp += [PAD_ID] * (max_len - len(ids_inp))
        ids_tgt += [PAD_ID] * (max_len - len(ids_tgt))
        seqs.append((ids_inp, ids_tgt))

    if bad:
        print(f"[WARN] selfies_to_tensor skipped {bad}")
    if not seqs:
        raise ValueError("No sequences created in selfies_to_tensor.")

    x = torch.tensor([a for a, _ in seqs], dtype=torch.long)
    y = torch.tensor([b for _, b in seqs], dtype=torch.long)
    return x, y


def prepare_shared_selfies_data(df):
    smiles = [canonicalize_smiles(s) for s in df["SMILES"].tolist()]
    smiles = [s for s in smiles if s is not None]
    stoi, itos, vocab, bos_id, eos_id = selfies_tokenize(smiles)
    x, y = selfies_to_tensor(smiles, stoi, MAX_LEN, bos_id, eos_id)
    print(
        f"[SELFIES] shared vocabulary={len(vocab)} training_sequences={len(x)} "
        f"automatically learned from Dm"
    )
    return {
        "smiles": smiles,
        "stoi": stoi,
        "itos": itos,
        "vocab": vocab,
        "bos_id": bos_id,
        "eos_id": eos_id,
        "x": x,
        "y": y,
    }


class MolRNN(nn.Module):
    def __init__(self, vocab_size, hidden=LATENT_DIM, extra_tokens=2):
        super().__init__()
        self.embed = nn.Embedding(
            vocab_size + 1 + extra_tokens, hidden, padding_idx=PAD_ID
        )
        self.gru = nn.GRU(hidden, hidden, batch_first=True)
        self.fc = nn.Linear(hidden, vocab_size + 1 + extra_tokens)

    def forward(self, x):
        emb = self.embed(x)
        out, _ = self.gru(emb)
        return self.fc(out)


class PositionalEncoding(nn.Module):
    def __init__(self, d_model, max_len=2048):
        super().__init__()
        pe = torch.zeros(max_len, d_model)
        position = torch.arange(0, max_len, dtype=torch.float).unsqueeze(1)
        div_term = torch.exp(
            torch.arange(0, d_model, 2).float() * (-np.log(10000.0) / d_model)
        )
        pe[:, 0::2] = torch.sin(position * div_term)
        pe[:, 1::2] = torch.cos(position * div_term)
        self.register_buffer("pe", pe.unsqueeze(0))

    def forward(self, x):
        return x + self.pe[:, : x.size(1), :]


class MolTransformer(nn.Module):
    def __init__(
        self,
        vocab_size,
        d_model=LATENT_DIM,
        nhead=8,
        num_layers=8,
        dim_ff=None,
        dropout=0.1,
        extra_tokens=2,
    ):
        super().__init__()
        dim_ff = dim_ff or 4 * d_model
        self.embed = nn.Embedding(
            vocab_size + 1 + extra_tokens, d_model, padding_idx=PAD_ID
        )
        self.pos = PositionalEncoding(d_model)
        enc_layer = nn.TransformerEncoderLayer(
            d_model=d_model,
            nhead=nhead,
            dim_feedforward=dim_ff,
            dropout=dropout,
            batch_first=True,
            activation="gelu",
        )
        self.encoder = nn.TransformerEncoder(enc_layer, num_layers=num_layers)
        self.lm_head = nn.Linear(d_model, vocab_size + 1 + extra_tokens)

    def forward(self, x, attn_mask=None):
        h = self.pos(self.embed(x))
        if attn_mask is None:
            T = x.size(1)
            attn_mask = torch.triu(torch.ones(T, T, device=x.device), diagonal=1).bool()
        pad_mask = x.eq(PAD_ID)
        out = self.encoder(h, mask=attn_mask, src_key_padding_mask=pad_mask)
        return self.lm_head(out)


class MolVAE(nn.Module):
    def __init__(
        self, vocab_size, hidden=LATENT_DIM, latent_dim=LATENT_DIM_VAE, extra_tokens=2
    ):
        super().__init__()
        self.vocab_size = vocab_size + 1 + extra_tokens
        self.hidden = hidden
        self.latent_dim = latent_dim
        self.embed = nn.Embedding(self.vocab_size, hidden, padding_idx=PAD_ID)
        self.enc_rnn = nn.GRU(hidden, hidden, batch_first=True)
        self.fc_mu = nn.Linear(hidden, latent_dim)
        self.fc_logvar = nn.Linear(hidden, latent_dim)
        self.dec_in = nn.Linear(latent_dim, hidden)
        self.dec_rnn = nn.GRU(hidden, hidden, batch_first=True)
        self.fc_out = nn.Linear(hidden, self.vocab_size)

    def encode(self, x):
        lengths = x.ne(PAD_ID).sum(dim=1).clamp(min=1).cpu()
        emb = self.embed(x)
        packed = torch.nn.utils.rnn.pack_padded_sequence(
            emb, lengths, batch_first=True, enforce_sorted=False
        )
        _, h = self.enc_rnn(packed)
        h = h[-1]
        return self.fc_mu(h), self.fc_logvar(h)

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z, x_dec):
        h0 = self.dec_in(z).unsqueeze(0)
        emb = self.embed(x_dec)
        out, _ = self.dec_rnn(emb, h0)
        return self.fc_out(out)

    def forward(self, x_enc, x_dec):
        mu, logvar = self.encode(x_enc)
        z = self.reparameterize(mu, logvar)
        logits = self.decode(z, x_dec)
        return logits, mu, logvar


def _make_loader(x, y):
    ds = torch.utils.data.TensorDataset(x, y)
    return torch.utils.data.DataLoader(
        ds,
        batch_size=BATCH_SIZE_GEN,
        shuffle=True,
        drop_last=False,
        pin_memory=(DEVICE == "cuda"),
    )


def train_rnn_model(shared):
    x, y = shared["x"], shared["y"]
    model = MolRNN(len(shared["vocab"])).to(DEVICE)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    loader = _make_loader(x, y)

    for epoch in range(EPOCHS_GEN):
        model.train()
        losses = []
        for xb, yb in loader:
            xb = xb.to(DEVICE, non_blocking=True)
            yb = yb.to(DEVICE, non_blocking=True)
            out = model(xb)
            loss = F.cross_entropy(
                out.reshape(-1, out.shape[-1]),
                yb.reshape(-1),
                ignore_index=PAD_ID,
                label_smoothing=0.05,
            )
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            losses.append(loss.item())
        print(f"[RNN] epoch {epoch+1}/{EPOCHS_GEN} loss={np.mean(losses):.4f}")

    torch.save(model.state_dict(), OUT_DIR / "molrnn.pt")
    return model


def train_transformer_model(shared):
    x, y = shared["x"], shared["y"]
    model = MolTransformer(len(shared["vocab"])).to(DEVICE)
    opt = torch.optim.AdamW(
        model.parameters(), lr=2e-4, betas=(0.9, 0.98), weight_decay=0.01
    )
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=EPOCHS_GEN)
    loader = _make_loader(x, y)

    for epoch in range(EPOCHS_GEN):
        model.train()
        losses = []
        for xb, yb in loader:
            xb = xb.to(DEVICE, non_blocking=True)
            yb = yb.to(DEVICE, non_blocking=True)
            Tlen = xb.size(1)
            attn_mask = torch.triu(
                torch.ones(Tlen, Tlen, device=DEVICE), diagonal=1
            ).bool()
            logits = model(xb, attn_mask=attn_mask)
            loss = F.cross_entropy(
                logits.reshape(-1, logits.size(-1)),
                yb.reshape(-1),
                ignore_index=PAD_ID,
                label_smoothing=0.05,
            )
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            losses.append(loss.item())
        sched.step()
        print(f"[Transformer] epoch {epoch+1}/{EPOCHS_GEN} loss={np.mean(losses):.4f}")

    torch.save(model.state_dict(), OUT_DIR / "moltransformer.pt")
    return model


def kl_divergence(mu, logvar):
    return -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1).mean()


def train_vae_model(shared):
    x_enc, x_tgt = shared["x"], shared["y"]
    model = MolVAE(len(shared["vocab"])).to(DEVICE)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    loader = _make_loader(x_enc, x_tgt)

    for epoch in range(EPOCHS_VAE):
        model.train()
        losses = []
        for xb_enc, yb in loader:
            xb_enc = xb_enc.to(DEVICE, non_blocking=True)
            yb = yb.to(DEVICE, non_blocking=True)
            logits, mu, logvar = model(xb_enc, xb_enc)
            rec_loss = F.cross_entropy(
                logits.reshape(-1, logits.size(-1)),
                yb.reshape(-1),
                ignore_index=PAD_ID,
                label_smoothing=0.05,
            )
            kl = kl_divergence(mu, logvar)
            coef = (
                BETA_VAE * min(1.0, float(epoch + 1) / float(EPOCHS_VAE))
                if KL_ANNEAL
                else BETA_VAE
            )
            loss = rec_loss + coef * kl
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            losses.append(loss.item())
        print(f"[VAE] epoch {epoch+1}/{EPOCHS_VAE} loss={np.mean(losses):.4f}")

    torch.save(model.state_dict(), OUT_DIR / "molvae.pt")
    return model


def top_k_top_p_filtering(logits, top_k=0, top_p=0.0, min_tokens_to_keep=1):
    B, V = logits.shape
    if top_k and top_k > 0:
        k = int(min(max(top_k, 1), V - min_tokens_to_keep))
        if k > 0:
            kth = torch.topk(logits, k, dim=-1).values[:, -1].unsqueeze(-1)
            logits = torch.where(
                logits < kth,
                torch.full_like(logits, float("-inf")),
                logits,
            )

    if 0.0 < top_p < 1.0:
        sorted_logits, sorted_idx = torch.sort(logits, descending=True, dim=-1)
        probs = F.softmax(sorted_logits, dim=-1)
        cum = torch.cumsum(probs, dim=-1)
        mask = cum > top_p
        mask[..., 1:] = mask[..., :-1].clone()
        mask[..., 0] = False
        sorted_logits = torch.where(
            mask, torch.full_like(sorted_logits, float("-inf")), sorted_logits
        )
        new_logits = torch.full_like(logits, float("-inf"))
        new_logits.scatter_(1, sorted_idx, sorted_logits)
        logits = new_logits

    all_neg_inf = torch.isinf(logits).all(dim=-1, keepdim=True)
    if all_neg_inf.any():
        logits = torch.where(all_neg_inf, torch.zeros_like(logits), logits)
    return logits


def decode_selfies_sequence(seq, itos, bos_id, eos_id):
    out = []
    for t in seq:
        t = int(t)
        if t == bos_id or t == PAD_ID:
            continue
        if t == eos_id:
            break
        out.append(itos.get(t, ""))

    selfies = "".join(out)
    if not selfies:
        return None

    try:
        sm = sf.decoder(selfies)
        return canonicalize_smiles(sm)
    except Exception:
        return None


def select_seed_smiles(seed_smiles, num_gen, round_id):
    if not seed_smiles:
        return []
    uniq = []
    seen = set()
    for s in seed_smiles:
        cs = canonicalize_smiles(s)
        if cs is not None and cs not in seen:
            seen.add(cs)
            uniq.append(cs)
    k = min(_auto_seed_count(len(uniq)), int(num_gen))
    if k <= 0:
        return []
    rng = random.Random(SEED + round_id)
    return uniq if len(uniq) <= k else rng.sample(uniq, k)


def smiles_to_prefix_ids(smiles, stoi, bos_id, max_len):
    try:
        ss = sf.encoder(smiles)
        ids = [stoi[t] for t in sf.split_selfies(ss) if t in stoi]
        ids = ids[: max_len - 1]
        return [bos_id] + ids
    except Exception:
        return None


def build_forced_prefix_tensor(seed_smiles, stoi, bos_id, num_gen, max_len):
    forced = torch.full((num_gen, max_len), -1, dtype=torch.long, device=DEVICE)
    forced[:, 0] = bos_id
    for i, sm in enumerate(seed_smiles[:num_gen]):
        ids = smiles_to_prefix_ids(sm, stoi, bos_id, max_len)
        if ids is None:
            continue
        forced[i, : len(ids)] = torch.tensor(ids, dtype=torch.long, device=DEVICE)
    return forced


def apply_repeat_penalty(logits, seq, bos_id, repeat_penalty):
    if repeat_penalty is None or repeat_penalty == 1.0:
        return logits
    penal = torch.zeros_like(logits, dtype=torch.bool)
    tok = seq.clamp(min=0)
    penal.scatter_(1, tok, torch.ones_like(tok, dtype=torch.bool))
    penal[:, PAD_ID] = False
    penal[:, bos_id] = False
    return logits / (penal.float() * (repeat_penalty - 1.0) + 1.0)


def choose_next_token(logits, seq, bos_id, eos_id, temp, top_k, top_p, repeat_penalty):
    logits = logits / max(1e-6, float(temp))
    logits = apply_repeat_penalty(logits, seq, bos_id, repeat_penalty)

    logits[:, PAD_ID] = float("-inf")
    logits[:, bos_id] = float("-inf")

    V = logits.size(-1)
    dyn_k = min(int(top_k), max(1, V - 2)) if top_k else 0
    logits = top_k_top_p_filtering(logits, top_k=dyn_k, top_p=top_p)
    probs = F.softmax(logits, dim=-1)
    return torch.multinomial(probs, 1).squeeze(1)


def _finalize_sequences(seqs, itos, bos_id, eos_id, tag):
    smiles = []
    for row in seqs:
        sm = decode_selfies_sequence(row, itos, bos_id, eos_id)
        if sm is not None:
            smiles.append(sm)
    return finalize_generated(smiles, tag=tag, exclude_dm=True)


def _sample_rnn_batch(
    model,
    stoi,
    itos,
    bos_id,
    eos_id,
    batch_n,
    seed_smiles,
    temp,
    top_k,
    top_p,
    repeat_penalty,
):
    forced = build_forced_prefix_tensor(seed_smiles, stoi, bos_id, batch_n, MAX_LEN)
    seq = torch.full((batch_n, 1), bos_id, dtype=torch.long, device=DEVICE)
    finished = torch.zeros(batch_n, dtype=torch.bool, device=DEVICE)

    model.eval()
    with torch.no_grad():
        while seq.size(1) < MAX_LEN:
            pos = seq.size(1)
            logits = model(seq)[:, -1, :]
            nxt = choose_next_token(
                logits, seq, bos_id, eos_id, temp, top_k, top_p, repeat_penalty
            )

            forced_tok = forced[:, pos]
            force_mask = forced_tok.ge(0) & (~finished)
            nxt = torch.where(force_mask, forced_tok, nxt)
            nxt = torch.where(finished, torch.full_like(nxt, PAD_ID), nxt)

            newly_finished = (~finished) & (nxt == eos_id)
            seq = torch.cat([seq, nxt.unsqueeze(1)], dim=1)
            finished |= newly_finished
            if finished.all():
                break

    return seq.detach().cpu().numpy()


def _sample_transformer_batch(
    model,
    stoi,
    itos,
    bos_id,
    eos_id,
    batch_n,
    seed_smiles,
    temp,
    top_k,
    top_p,
    repeat_penalty,
):
    forced = build_forced_prefix_tensor(seed_smiles, stoi, bos_id, batch_n, MAX_LEN)
    seq = torch.full((batch_n, 1), bos_id, dtype=torch.long, device=DEVICE)
    finished = torch.zeros(batch_n, dtype=torch.bool, device=DEVICE)

    model.eval()
    with torch.no_grad():
        while seq.size(1) < MAX_LEN:
            pos = seq.size(1)
            Tlen = seq.size(1)
            attn_mask = torch.triu(
                torch.ones(Tlen, Tlen, device=DEVICE), diagonal=1
            ).bool()
            logits = model(seq, attn_mask=attn_mask)[:, -1, :]
            nxt = choose_next_token(
                logits, seq, bos_id, eos_id, temp, top_k, top_p, repeat_penalty
            )

            forced_tok = forced[:, pos]
            force_mask = forced_tok.ge(0) & (~finished)
            nxt = torch.where(force_mask, forced_tok, nxt)
            nxt = torch.where(finished, torch.full_like(nxt, PAD_ID), nxt)

            newly_finished = (~finished) & (nxt == eos_id)
            seq = torch.cat([seq, nxt.unsqueeze(1)], dim=1)
            finished |= newly_finished
            if finished.all():
                break

    return seq.detach().cpu().numpy()


def seed_tensor_for_vae(seed_smiles, shared):
    if not seed_smiles:
        return None
    x, _ = selfies_to_tensor(
        seed_smiles,
        shared["stoi"],
        MAX_LEN,
        shared["bos_id"],
        shared["eos_id"],
    )
    return x


def _sample_vae_batch(
    model, shared, batch_n, seed_smiles, temp, top_k, top_p, repeat_penalty
):
    stoi = shared["stoi"]
    bos_id = shared["bos_id"]
    eos_id = shared["eos_id"]

    forced = build_forced_prefix_tensor(seed_smiles, stoi, bos_id, batch_n, MAX_LEN)

    model.eval()
    with torch.no_grad():
        k = min(len(seed_smiles), batch_n)
        if k > 0:
            xseed = seed_tensor_for_vae(seed_smiles[:k], shared).to(DEVICE)
            mu, logvar = model.encode(xseed)
            z_seed = model.reparameterize(mu, logvar)
        else:
            z_seed = torch.empty((0, model.latent_dim), device=DEVICE)

        if k < batch_n:
            z_rand = torch.randn(batch_n - k, model.latent_dim, device=DEVICE)
            z = torch.cat([z_seed, z_rand], dim=0)
        else:
            z = z_seed

        seq = torch.full((batch_n, 1), bos_id, dtype=torch.long, device=DEVICE)
        finished = torch.zeros(batch_n, dtype=torch.bool, device=DEVICE)

        while seq.size(1) < MAX_LEN:
            pos = seq.size(1)
            logits = model.decode(z, seq)[:, -1, :]
            nxt = choose_next_token(
                logits, seq, bos_id, eos_id, temp, top_k, top_p, repeat_penalty
            )

            forced_tok = forced[:, pos]
            force_mask = forced_tok.ge(0) & (~finished)
            nxt = torch.where(force_mask, forced_tok, nxt)
            nxt = torch.where(finished, torch.full_like(nxt, PAD_ID), nxt)

            newly_finished = (~finished) & (nxt == eos_id)
            seq = torch.cat([seq, nxt.unsqueeze(1)], dim=1)
            finished |= newly_finished
            if finished.all():
                break

    return seq.detach().cpu().numpy()


def _sample_in_batches(
    sampler,
    model,
    shared,
    num_gen,
    seed_smiles,
    round_id,
    tag,
    temp=TEMPERATURE,
    top_k=TOP_K,
    top_p=TOP_P,
    repeat_penalty=REPEAT_PENALTY,
):
    selected = select_seed_smiles(seed_smiles, num_gen, round_id)
    all_sequences = []
    generated = 0
    seed_offset = 0

    while generated < num_gen:
        n = min(SAMPLE_BATCH_SIZE, num_gen - generated)
        batch_seeds = selected[seed_offset : seed_offset + n]
        seed_offset += len(batch_seeds)

        if sampler == "rnn":
            seq = _sample_rnn_batch(
                model,
                shared["stoi"],
                shared["itos"],
                shared["bos_id"],
                shared["eos_id"],
                n,
                batch_seeds,
                temp,
                top_k,
                top_p,
                repeat_penalty,
            )
        elif sampler == "transformer":
            seq = _sample_transformer_batch(
                model,
                shared["stoi"],
                shared["itos"],
                shared["bos_id"],
                shared["eos_id"],
                n,
                batch_seeds,
                temp,
                top_k,
                top_p,
                repeat_penalty,
            )
        elif sampler == "vae":
            seq = _sample_vae_batch(
                model, shared, n, batch_seeds, temp, top_k, top_p, repeat_penalty
            )
        else:
            raise ValueError(sampler)

        all_sequences.extend(seq)
        generated += n

    smiles = _finalize_sequences(
        all_sequences,
        shared["itos"],
        shared["bos_id"],
        shared["eos_id"],
        tag=tag,
    )
    return smiles


def sample_rnn(model, shared, seed_smiles=None, round_id=0):
    smiles = _sample_in_batches(
        "rnn", model, shared, NUM_GEN, seed_smiles, round_id, tag=f"RNN_R{round_id}"
    )
    out_csv = OUT_DIR / f"generated_smiles_RNN_R{round_id}.csv"
    pd.DataFrame({"smiles": smiles}).to_csv(out_csv, index=False)
    return out_csv


def sample_transformer(model, shared, seed_smiles=None, round_id=0):
    smiles = _sample_in_batches(
        "transformer",
        model,
        shared,
        NUM_GEN,
        seed_smiles,
        round_id,
        tag=f"Transformer_R{round_id}",
    )
    out_csv = OUT_DIR / f"generated_smiles_Transformer_R{round_id}.csv"
    pd.DataFrame({"smiles": smiles}).to_csv(out_csv, index=False)
    return out_csv


def sample_vae(model, shared, seed_smiles=None, round_id=0):
    smiles = _sample_in_batches(
        "vae", model, shared, NUM_GEN, seed_smiles, round_id, tag=f"VAE_R{round_id}"
    )
    out_csv = OUT_DIR / f"generated_smiles_VAE_R{round_id}.csv"
    pd.DataFrame({"smiles": smiles}).to_csv(out_csv, index=False)
    return out_csv


def main():
    df = load_dm_dataframe(DATA_CSV)

    ckpt_root = chemprop_kfold_train(df)

    shared = prepare_shared_selfies_data(df)
    rnn = train_rnn_model(shared)
    transformer = train_transformer_model(shared)
    vae = train_vae_model(shared)

    gen_rnn = sample_rnn(rnn, shared, seed_smiles=None, round_id=0)
    gen_tf = sample_transformer(transformer, shared, seed_smiles=None, round_id=0)
    gen_vae = sample_vae(vae, shared, seed_smiles=None, round_id=0)

    enrich_by_predictions(ckpt_root, gen_rnn, "RNN_R0")
    enrich_by_predictions(ckpt_root, gen_tf, "Transformer_R0")
    enrich_by_predictions(ckpt_root, gen_vae, "VAE_R0")

    merged_prev = merge_and_dedup(
        [gen_rnn, gen_tf, gen_vae], "generated_merged_round0.csv"
    )
    merged_pred, _ = chemprop_predict_paths(ckpt_root, merged_prev, "MERGED_R0")
    seeds = extract_top_seeds(merged_pred)

    for r in range(1, ENRICH_ROUNDS + 1):
        gen_rnn = sample_rnn(rnn, shared, seed_smiles=seeds, round_id=r)
        gen_tf = sample_transformer(transformer, shared, seed_smiles=seeds, round_id=r)
        gen_vae = sample_vae(vae, shared, seed_smiles=seeds, round_id=r)

        enrich_by_predictions(ckpt_root, gen_rnn, f"RNN_R{r}")
        enrich_by_predictions(ckpt_root, gen_tf, f"Transformer_R{r}")
        enrich_by_predictions(ckpt_root, gen_vae, f"VAE_R{r}")

        merged_current = merge_and_dedup(
            [merged_prev, gen_rnn, gen_tf, gen_vae], f"generated_merged_round{r}.csv"
        )
        merged_pred, _ = chemprop_predict_paths(
            ckpt_root, merged_current, f"MERGED_R{r}"
        )
        seeds = extract_top_seeds(merged_pred)
        merged_prev = merged_current

    print("Done.")


if __name__ == "__main__":
    main()

import os, sys, re, random, shutil, subprocess, glob, math

os.chdir(os.path.split(os.path.realpath(__file__))[0])
from pathlib import Path
from functools import lru_cache
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from rdkit import Chem
from rdkit.Chem import AllChem, BRICS
import selfies as sf

os.environ["MPLBACKEND"] = "Agg"
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import warnings

warnings.filterwarnings("ignore")

SEED = 42
random.seed(SEED)
np.random.seed(SEED)
torch.manual_seed(SEED)
if torch.cuda.is_available():
    torch.cuda.manual_seed_all(SEED)
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

DATA_CSV = "D.csv"
BASE_DIR = Path("runs/mol_chemprop_multi")
CP_DIR = BASE_DIR / "models_chemprop"
OUT_DIR = BASE_DIR / "outputs"
for p in [CP_DIR, OUT_DIR]:
    p.mkdir(parents=True, exist_ok=True)

TARGET_D = 20
DM_DUPLICATE_MODE = "mean"
EXCLUDED_ELEMENTS = []
SIZE_FILTER_QUANTILES = (0.01, 0.99)

NUM_GEN = 1024
ENRICH_ROUNDS = 5

T_STEPS = 200
EPOCHS = 30
LR = 1e-3

SELFIES_MIN_CONTENT_TOK = 4
SELFIES_HIDDEN = 384
SELFIES_TOPK = None
SELFIES_TOPP = 0.95

BATCH_SIZE = 0
MIN_VALID = 0
SELFIES_MAX_TOK = 0
SELFIES_PREFIX_MAX = 0

GRAPH_HIDDEN = 512
GRAPH_SIGMA_MIN = 0.01
GRAPH_SIGMA_MAX = 1.0
GRAPH_STEPS = 60

GEODIFF_HIDDEN = 512
GEODIFF_SIGMA_MIN = 0.01
GEODIFF_SIGMA_MAX = 1.0
GEODIFF_STEPS = 60

DM_CANONICAL_SET = set()
DM_STRUCTURE_SET = set()
DM_HEAVY_ATOM_MIN = None
DM_HEAVY_ATOM_MAX = None
TARGET_DIRECTION_RESOLVED = None
GRAPH_ATOMS = []
GRAPH_TYPES = []
GRAPH_ATOM_TO_INDEX = {}
GRAPH_NODE_DIM = 0
GRAPH_N_MAX = 0
GEODIFF_ATOMS = []
GEODIFF_TYPES = []
GEODIFF_ATOM_TO_INDEX = {}
GEODIFF_NODE_DIM = 0
GEODIFF_N_MAX = 0
EXCLUDED_ELEMENT_SET = set()

BOND_NONE = 0
BOND_SINGLE = 1
BOND_DOUBLE = 2
BOND_TRIPLE = 3
BOND_AROMATIC = 4
BOND_CODE_MAX = 4.0


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
        z = Chem.GetPeriodicTable().GetAtomicNumber(s)
    except Exception:
        z = 0
    if not z:
        raise ValueError(f"Unknown element symbol in EXCLUDED_ELEMENTS: {symbol}")
    return s


def initialize_excluded_elements():
    global EXCLUDED_ELEMENT_SET
    EXCLUDED_ELEMENT_SET = {
        normalize_element_symbol(x) for x in EXCLUDED_ELEMENTS if str(x).strip()
    }
    if EXCLUDED_ELEMENT_SET:
        print("[ELEMENT FILTER] excluded =", sorted(EXCLUDED_ELEMENT_SET))
    else:
        print("[ELEMENT FILTER] disabled: EXCLUDED_ELEMENTS is empty")


@lru_cache(maxsize=300000)
def _canonicalize_smiles_cached(s):
    try:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    except Exception:
        return None


@lru_cache(maxsize=300000)
def _structure_key_cached(cs):
    try:
        mol = Chem.MolFromSmiles(cs)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True, isomericSmiles=False)
    except Exception:
        return None


def structure_key(smiles):
    cs = canonicalize_smiles(smiles)
    return None if cs is None else _structure_key_cached(cs)


def canonicalize_smiles(smiles):
    if smiles is None:
        return None
    s = str(smiles).strip()
    if not s or s.lower() == "nan":
        return None
    return _canonicalize_smiles_cached(s)


@lru_cache(maxsize=300000)
def _is_neutral_canonical_cached(cs):
    try:
        mol = Chem.MolFromSmiles(cs)
        return mol is not None and Chem.rdmolops.GetFormalCharge(mol) == 0
    except Exception:
        return False


def is_neutral_smiles(smiles):
    cs = canonicalize_smiles(smiles)
    return False if cs is None else _is_neutral_canonical_cached(cs)


@lru_cache(maxsize=300000)
def _contains_excluded_elements_cached(cs):
    if not EXCLUDED_ELEMENT_SET:
        return False
    mol = Chem.MolFromSmiles(cs)
    if mol is None:
        return False
    return any(a.GetSymbol() in EXCLUDED_ELEMENT_SET for a in mol.GetAtoms())


def contains_excluded_elements(smiles):
    if not EXCLUDED_ELEMENT_SET:
        return False
    cs = canonicalize_smiles(smiles)
    return False if cs is None else _contains_excluded_elements_cached(cs)


@lru_cache(maxsize=300000)
def _heavy_atom_count_cached(cs):
    mol = Chem.MolFromSmiles(cs)
    return -1 if mol is None else int(mol.GetNumHeavyAtoms())


def _aggregate_duplicate_values(group, mode):
    vals = pd.to_numeric(group["D"], errors="coerce").dropna()
    if len(vals) == 0:
        return np.nan
    mode = str(mode).strip().lower()
    if mode == "first":
        return float(vals.iloc[0])
    if mode == "mean":
        return float(vals.mean())
    if mode == "median":
        return float(vals.median())
    if mode == "max":
        return float(vals.max())
    if mode == "min":
        return float(vals.min())
    raise ValueError("DM_DUPLICATE_MODE must be one of: first, mean, median, max, min")


def load_dm_dataframe(path):
    global DM_CANONICAL_SET, DM_STRUCTURE_SET
    df = pd.read_csv(path)
    smiles_col = None
    d_col = None
    for c in df.columns:
        if smiles_col is None and re.search("smiles", str(c), re.I):
            smiles_col = c
        if d_col is None and re.search(r"\bD\b|velocity|target|label", str(c), re.I):
            d_col = c
    if smiles_col is None or d_col is None:
        raise ValueError("Missing SMILES or target column")

    raw = (
        df[[smiles_col, d_col]]
        .rename(columns={smiles_col: "SMILES", d_col: "D"})
        .dropna()
        .reset_index(drop=True)
    )
    raw["D"] = pd.to_numeric(raw["D"], errors="coerce")
    raw = raw.dropna(subset=["D"]).reset_index(drop=True)

    valid = []
    bad_rows = []
    for i, row in raw.iterrows():
        cs = canonicalize_smiles(row["SMILES"])
        if cs is None:
            bad_rows.append((i, str(row["SMILES"])))
            continue
        valid.append((cs, float(row["D"])))

    if bad_rows:
        bad_path = OUT_DIR / "invalid_smiles_in_Dm.csv"
        pd.DataFrame(bad_rows, columns=["row_index", "SMILES"]).to_csv(
            bad_path, index=False
        )
        print(f"[WARN] invalid SMILES filtered: {len(bad_rows)} -> {bad_path}")

    vdf = pd.DataFrame(valid, columns=["SMILES", "D"])
    before = len(vdf)
    records = []
    for cs, g in vdf.groupby("SMILES", sort=False):
        records.append((cs, _aggregate_duplicate_values(g, DM_DUPLICATE_MODE), len(g)))
    out = pd.DataFrame(records, columns=["SMILES", "D", "DuplicateCount"])
    out = out.dropna(subset=["D"]).reset_index(drop=True)
    removed = before - len(out)
    DM_CANONICAL_SET = set(out["SMILES"].tolist())
    DM_STRUCTURE_SET = {structure_key(x) for x in out["SMILES"].tolist()}
    DM_STRUCTURE_SET.discard(None)
    print(
        f"[DEDUP:Dm] input={len(raw)} invalid={len(bad_rows)} "
        f"removed_internal={removed} kept={len(out)} mode={DM_DUPLICATE_MODE}"
    )
    return out


def resolve_target_direction(df):
    global TARGET_DIRECTION_RESOLVED
    median_d = float(pd.to_numeric(df["D"], errors="coerce").median())
    TARGET_DIRECTION_RESOLVED = "max" if float(TARGET_D) >= median_d else "min"
    print(
        f"[TARGET] threshold={TARGET_D} training_median={median_d:.6f} "
        f"direction={TARGET_DIRECTION_RESOLVED}"
    )
    return TARGET_DIRECTION_RESOLVED


def _align_up(value, multiple=8):
    value = int(value)
    return int(math.ceil(value / float(multiple)) * multiple)


def _auto_batch_size(selfies_max_len, graph_n_max):
    if DEVICE != "cuda":
        return 16
    try:
        gb = torch.cuda.get_device_properties(0).total_memory / (1024**3)
    except Exception:
        gb = 12.0
    if gb >= 28:
        batch = 64
    elif gb >= 20:
        batch = 48
    elif gb >= 14:
        batch = 32
    elif gb >= 10:
        batch = 16
    else:
        batch = 8
    if graph_n_max > 192 or selfies_max_len > 256:
        batch = max(8, batch // 2)
    if graph_n_max > 256 or selfies_max_len > 384:
        batch = max(4, batch // 2)
    return int(batch)


def _auto_seed_count(n_candidates):
    n_candidates = int(n_candidates)
    if n_candidates <= 0:
        return 0
    k = int(math.ceil(n_candidates * 0.10))
    k = max(32, min(256, k))
    return min(n_candidates, k)


def initialize_chemical_space(df):
    global GRAPH_ATOMS, GRAPH_TYPES, GRAPH_ATOM_TO_INDEX, GRAPH_NODE_DIM, GRAPH_N_MAX
    global GEODIFF_ATOMS, GEODIFF_TYPES, GEODIFF_ATOM_TO_INDEX, GEODIFF_NODE_DIM, GEODIFF_N_MAX
    global DM_HEAVY_ATOM_MIN, DM_HEAVY_ATOM_MAX
    global SELFIES_MAX_TOK, SELFIES_PREFIX_MAX, MIN_VALID, BATCH_SIZE

    pt = Chem.GetPeriodicTable()
    elements = set()
    heavy_counts = []
    selfies_lengths = []

    for cs in df["SMILES"].tolist():
        mol = Chem.MolFromSmiles(cs)
        if mol is None:
            continue
        heavy_counts.append(int(mol.GetNumHeavyAtoms()))
        for atom in mol.GetAtoms():
            sym = atom.GetSymbol()
            if sym != "H":
                elements.add(sym)
        try:
            ss = sf.encoder(cs)
            selfies_lengths.append(len(list(sf.split_selfies(ss))))
        except Exception:
            pass

    if not elements or not heavy_counts or not selfies_lengths:
        raise RuntimeError("Cannot initialize chemical space from Dm.csv")

    GRAPH_ATOMS = sorted(elements, key=lambda x: pt.GetAtomicNumber(x))
    GRAPH_TYPES = ["<PAD>"] + GRAPH_ATOMS
    GRAPH_ATOM_TO_INDEX = {s: i for i, s in enumerate(GRAPH_TYPES)}
    GRAPH_NODE_DIM = len(GRAPH_TYPES)

    GEODIFF_ATOMS = GRAPH_ATOMS.copy()
    GEODIFF_TYPES = ["<PAD>"] + GEODIFF_ATOMS
    GEODIFF_ATOM_TO_INDEX = {s: i for i, s in enumerate(GEODIFF_TYPES)}
    GEODIFF_NODE_DIM = len(GEODIFF_TYPES)

    max_heavy = int(max(heavy_counts))
    GRAPH_N_MAX = _align_up(max_heavy + 2, 8)
    GEODIFF_N_MAX = GRAPH_N_MAX

    qlow, qhigh = SIZE_FILTER_QUANTILES
    if not (0 <= qlow < qhigh <= 1):
        raise ValueError("SIZE_FILTER_QUANTILES must satisfy 0 <= low < high <= 1")
    hs = pd.Series(heavy_counts, dtype=float)
    DM_HEAVY_ATOM_MIN = max(1, int(math.floor(float(hs.quantile(qlow)))))
    DM_HEAVY_ATOM_MAX = int(math.ceil(float(hs.quantile(qhigh))))

    ls = pd.Series(selfies_lengths, dtype=float)
    max_content = int(max(selfies_lengths))
    SELFIES_MAX_TOK = _align_up(max_content + 2, 8)
    median_content = float(ls.median())
    SELFIES_PREFIX_MAX = int(round(median_content * 0.25))
    SELFIES_PREFIX_MAX = max(10, min(40, SELFIES_PREFIX_MAX))
    SELFIES_PREFIX_MAX = min(SELFIES_PREFIX_MAX, max(2, SELFIES_MAX_TOK - 2))

    MIN_VALID = min(NUM_GEN, max(32, int(math.ceil(NUM_GEN * 0.125))))
    BATCH_SIZE = _auto_batch_size(SELFIES_MAX_TOK, GRAPH_N_MAX)

    print("[CHEM SPACE] Graph/Geo elements =", GRAPH_ATOMS)
    print(f"[CHEM SPACE] node capacity = {GRAPH_N_MAX}")
    print(
        f"[SIZE FILTER] heavy atoms = {DM_HEAVY_ATOM_MIN}..{DM_HEAVY_ATOM_MAX} "
        f"from quantiles {SIZE_FILTER_QUANTILES}"
    )
    print(
        f"[AUTO SELFIES] P50={ls.quantile(0.50):.1f} P95={ls.quantile(0.95):.1f} "
        f"P99={ls.quantile(0.99):.1f} max={max_content} "
        f"SELFIES_MAX_TOK={SELFIES_MAX_TOK} PREFIX_MAX={SELFIES_PREFIX_MAX}"
    )
    print(f"[AUTO GENERATION] MIN_VALID={MIN_VALID} from NUM_GEN={NUM_GEN}")
    print(f"[AUTO HARDWARE] BATCH_SIZE={BATCH_SIZE} device={DEVICE}")

    cfg_path = OUT_DIR / "auto_config_diffusion.txt"
    with open(cfg_path, "w", encoding="utf-8") as f:
        f.write(f"unique_molecules={len(df)}\n")
        f.write(f"elements={','.join(GRAPH_ATOMS)}\n")
        f.write(f"selfies_p50={ls.quantile(0.50):.6f}\n")
        f.write(f"selfies_p95={ls.quantile(0.95):.6f}\n")
        f.write(f"selfies_p99={ls.quantile(0.99):.6f}\n")
        f.write(f"selfies_max_content={max_content}\n")
        f.write(f"SELFIES_MAX_TOK={SELFIES_MAX_TOK}\n")
        f.write(f"SELFIES_PREFIX_MAX={SELFIES_PREFIX_MAX}\n")
        f.write(f"GRAPH_N_MAX={GRAPH_N_MAX}\n")
        f.write(f"GEODIFF_N_MAX={GEODIFF_N_MAX}\n")
        f.write(f"heavy_atom_min={DM_HEAVY_ATOM_MIN}\n")
        f.write(f"heavy_atom_max={DM_HEAVY_ATOM_MAX}\n")
        f.write(f"BATCH_SIZE={BATCH_SIZE}\n")
        f.write(f"MIN_VALID={MIN_VALID}\n")
        f.write(f"TARGET_D={TARGET_D}\n")
        f.write(f"direction={TARGET_DIRECTION_RESOLVED}\n")


def _candidate_passes_size_filter(cs):
    if DM_HEAVY_ATOM_MIN is None or DM_HEAVY_ATOM_MAX is None:
        return True
    n = _heavy_atom_count_cached(cs)
    return DM_HEAVY_ATOM_MIN <= n <= DM_HEAVY_ATOM_MAX


def _extend_unique_generated(target, seen, candidates, exclude_dm=True):
    added = 0
    duplicate_count = 0
    dm_count = 0
    invalid_count = 0
    nonneutral_count = 0
    element_count = 0
    size_count = 0
    dm_set = DM_STRUCTURE_SET if exclude_dm else set()

    for s in candidates:
        cs = canonicalize_smiles(s)
        if cs is None:
            invalid_count += 1
            continue
        skey = _structure_key_cached(cs)
        if skey is None:
            invalid_count += 1
            continue
        if not _is_neutral_canonical_cached(cs):
            nonneutral_count += 1
            continue
        if EXCLUDED_ELEMENT_SET and _contains_excluded_elements_cached(cs):
            element_count += 1
            continue
        if not _candidate_passes_size_filter(cs):
            size_count += 1
            continue
        if skey in seen:
            duplicate_count += 1
            continue
        if skey in dm_set:
            dm_count += 1
            continue
        seen.add(skey)
        target.append(cs)
        added += 1

    return (
        added,
        duplicate_count,
        dm_count,
        invalid_count,
        nonneutral_count,
        element_count,
        size_count,
    )


def canonicalize_and_filter_neutral(smiles_list, exclude_dm=False):
    out = []
    seen = set()
    _extend_unique_generated(out, seen, smiles_list, exclude_dm=exclude_dm)
    return out


def finalize_generated(smiles_list, tag="Generated"):
    out = []
    seen = set()
    added, dup, dmdup, invalid, nonneutral, elem, size = _extend_unique_generated(
        out, seen, smiles_list, exclude_dm=True
    )
    print(
        f"[FILTER:{tag}] input={len(smiles_list)} kept={added} "
        f"dup={dup} vs_Dm={dmdup} invalid={invalid} nonneutral={nonneutral} "
        f"excluded_element={elem} size={size}"
    )
    return out


def make_fixed_cv_splits(n_samples, n_folds=5, seed=SEED):
    rng = np.random.RandomState(seed)
    idx = np.arange(n_samples)
    rng.shuffle(idx)
    folds = np.array_split(idx, n_folds)
    splits = []
    n_val = max(1, int(round(0.10 * n_samples)))
    for k in range(n_folds):
        test_idx = np.array(folds[k], dtype=int)
        rest = np.concatenate(
            [np.array(folds[j], dtype=int) for j in range(n_folds) if j != k]
        )
        rng_k = np.random.RandomState(seed + 1000 + k)
        rng_k.shuffle(rest)
        val_idx = rest[: min(n_val, max(1, len(rest) - 1))]
        train_idx = rest[len(val_idx) :]
        splits.append((np.sort(train_idx), np.sort(val_idx), np.sort(test_idx)))
    return splits


def find_fold_checkpoint(fold_dir):
    patterns = [
        fold_dir / "model_0" / "model.pt",
        fold_dir / "model.pt",
    ]
    for p in patterns:
        if p.exists():
            return str(p)
    cands = glob.glob(str((fold_dir / "**" / "*.pt").as_posix()), recursive=True)
    return cands[0] if cands else None


def mae(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    return float(np.mean(np.abs(y_true - y_pred)))


def rmse(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    return float(np.sqrt(np.mean((y_true - y_pred) ** 2)))


def r2(y_true, y_pred):
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
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


def _to_chemprop_csv(df, path, include_target=True):
    if include_target:
        df[["SMILES", "D"]].rename(columns={"SMILES": "smiles"}).to_csv(
            path, index=False
        )
    else:
        df[["SMILES"]].rename(columns={"SMILES": "smiles"}).to_csv(path, index=False)


def kfold_eval_and_artifacts(df, split_records, save_dir, outs_dir):
    cols = {}
    per_fold = []
    colors = [plt.get_cmap("tab10")(i) for i in range(len(split_records))]
    fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
    xlim = [float("inf"), -float("inf")]
    ylim = [float("inf"), -float("inf")]

    for k, rec in enumerate(split_records):
        test_df = df.iloc[rec["test_idx"]].reset_index(drop=True)
        test_input = outs_dir / f"cv_test_input_fold{k}.csv"
        pred_csv = outs_dir / f"preds_test_fold{k}.csv"
        _to_chemprop_csv(test_df, test_input, include_target=False)
        ck = find_fold_checkpoint(save_dir / f"fold_{k}")
        if ck is None:
            raise FileNotFoundError(f"No checkpoint found for fold {k}")
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
        pred_col = pick_pred_column(pr)
        if pred_col is None:
            raise ValueError(f"No numeric prediction column found in {pred_csv}")
        y_true = test_df["D"].astype(float).to_numpy()
        y_pred = pr[pred_col].astype(float).to_numpy()
        if len(y_pred) != len(y_true):
            raise ValueError(f"Prediction length mismatch in fold {k}")
        m_mae = mae(y_true, y_pred)
        m_rmse = rmse(y_true, y_pred)
        m_r2 = r2(y_true, y_pred)
        per_fold.append((k, m_mae, m_rmse, m_r2, len(y_true)))
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
    ax.set_title("Fixed 5-fold test predictions")
    ax.legend(frameon=False, markerscale=1.2, fontsize=8)
    ax.set_xlim(low, high)
    ax.set_ylim(low, high)
    fig.tight_layout()
    fig_path = outs_dir / "kfold_test_pred_300dpi.jpg"
    fig.savefig(fig_path, dpi=300)
    plt.close(fig)

    max_len = max(len(v) for v in cols.values())
    dat = pd.DataFrame({k: v + [np.nan] * (max_len - len(v)) for k, v in cols.items()})
    dat_path = outs_dir / "kfold_test_preds.dat"
    dat.to_csv(dat_path, sep="\t", index=False)

    maes = [p[1] for p in per_fold]
    rmses = [p[2] for p in per_fold]
    r2s = [p[3] for p in per_fold]
    mlog = outs_dir / "verbose.log"
    with open(mlog, "w", encoding="utf-8") as f:
        f.write("Fixed 5-fold cross validation\n")
        for k, mma, mmr, rrr, n in per_fold:
            f.write(f"\tFold {k}: n={n}, MAE={mma:.6f}, RMSE={mmr:.6f}, R2={rrr:.6f}\n")
        f.write(f"Overall MAE = {np.mean(maes):.6f} +/- {np.std(maes, ddof=1):.6f}\n")
        f.write(
            f"Overall RMSE = {np.mean(rmses):.6f} +/- {np.std(rmses, ddof=1):.6f}\n"
        )
        f.write(f"Overall R2 = {np.mean(r2s):.6f} +/- {np.std(r2s, ddof=1):.6f}\n")
    print(f"[KFold] metrics log -> {mlog}")
    print(f"[KFold] dat -> {dat_path}")
    print(f"[KFold] figure -> {fig_path}")


def chemprop_kfold_train(df, n_folds=5):
    if CP_DIR.exists():
        shutil.rmtree(CP_DIR)
    CP_DIR.mkdir(parents=True, exist_ok=True)

    splits = make_fixed_cv_splits(len(df), n_folds=n_folds, seed=SEED)
    split_records = []
    for k, (train_idx, val_idx, test_idx) in enumerate(splits):
        fold_dir = CP_DIR / f"fold_{k}"
        fold_dir.mkdir(parents=True, exist_ok=True)
        train_df = df.iloc[train_idx].reset_index(drop=True)
        val_df = df.iloc[val_idx].reset_index(drop=True)
        test_df = df.iloc[test_idx].reset_index(drop=True)
        train_csv = fold_dir / "train.csv"
        val_csv = fold_dir / "val.csv"
        test_csv = fold_dir / "test.csv"
        _to_chemprop_csv(train_df, train_csv, True)
        _to_chemprop_csv(val_df, val_csv, True)
        _to_chemprop_csv(test_df, test_csv, True)
        split_records.append(
            {
                "train_idx": train_idx,
                "val_idx": val_idx,
                "test_idx": test_idx,
            }
        )
        print(
            f"[CV fold {k}] train={len(train_df)} val={len(val_df)} test={len(test_df)}"
        )
        run(
            [
                CHEMPROP_TRAIN,
                "--data_path",
                str(train_csv),
                "--separate_val_path",
                str(val_csv),
                "--separate_test_path",
                str(test_csv),
                "--dataset_type",
                "regression",
                "--target_columns",
                "D",
                "--save_dir",
                str(fold_dir),
                "--epochs",
                "200",
                "--batch_size",
                "2",
                "--metric",
                "rmse",
                "--seed",
                str(SEED + k),
                "--quiet",
            ]
        )

    try:
        kfold_eval_and_artifacts(df, split_records, CP_DIR, OUT_DIR)
    except Exception as e:
        print("[WARN] kfold post-analysis failed:", e)
    return CP_DIR


def target_hit_mask(values):
    vals = pd.to_numeric(values, errors="coerce")
    if TARGET_DIRECTION_RESOLVED == "min":
        return vals <= float(TARGET_D)
    return vals >= float(TARGET_D)


def chemprop_predict_paths(ckpt_dir, gen_csv, tag):
    out_csv = OUT_DIR / f"generated_with_pred_{tag}.csv"
    hits_csv = OUT_DIR / f"generated_hits_{tag}.csv"
    gen_csv = str(gen_csv)
    if (not os.path.exists(gen_csv)) or os.path.getsize(gen_csv) == 0:
        return None, None
    try:
        run(
            [
                CHEMPROP_PREDICT,
                "--test_path",
                gen_csv,
                "--checkpoint_dir",
                str(ckpt_dir),
                "--preds_path",
                str(out_csv),
            ]
        )
    except subprocess.CalledProcessError:
        print("[WARN] chemprop_predict failed")
        return None, None
    if not out_csv.exists() or out_csv.stat().st_size == 0:
        return None, None
    df = pd.read_csv(out_csv)
    if "smiles" not in df.columns:
        if "SMILES" in df.columns:
            df = df.rename(columns={"SMILES": "smiles"})
        elif "smiles_0" in df.columns:
            df = df.rename(columns={"smiles_0": "smiles"})
    pred_col = pick_pred_column(df)
    if pred_col is None:
        return str(out_csv), None
    df["D_pred"] = pd.to_numeric(df[pred_col], errors="coerce")
    ascending = TARGET_DIRECTION_RESOLVED == "min"
    df = df.sort_values("D_pred", ascending=ascending, na_position="last").reset_index(
        drop=True
    )
    df.to_csv(out_csv, index=False)
    hits = df[target_hit_mask(df["D_pred"])].copy()
    hits.to_csv(hits_csv, index=False)
    print(
        f"[OK:{tag}] pred={len(df)} hits={len(hits)} "
        f"direction={TARGET_DIRECTION_RESOLVED} threshold={TARGET_D}"
    )
    return str(out_csv), str(hits_csv)


def to_selfies(sm):
    try:
        cs = canonicalize_smiles(sm)
        return None if cs is None else sf.encoder(cs)
    except Exception:
        return None


def from_selfies(s):
    try:
        sm = sf.decoder(s)
        return canonicalize_smiles(sm)
    except Exception:
        return None


def build_selfies_vocab(selfies_list):
    sym = set()
    for s in selfies_list:
        try:
            sym.update(sf.split_selfies(s))
        except Exception:
            pass
    itoks = ["<PAD>", "<BOS>", "<EOS>", "<MASK>"] + sorted(sym)
    stoi = {t: i for i, t in enumerate(itoks)}
    itos = {i: t for t, i in stoi.items()}
    return stoi, itos


def selfies_to_ids(s, stoi, max_len):
    bos, eos, pad, mask = stoi["<BOS>"], stoi["<EOS>"], stoi["<PAD>"], stoi["<MASK>"]
    toks = [bos] + [stoi.get(t, mask) for t in sf.split_selfies(s)] + [eos]
    toks = toks[:max_len]
    if toks[-1] != eos:
        toks[-1] = eos
    x = np.full(max_len, pad, dtype=np.int64)
    x[: len(toks)] = toks
    return x


def selfies_to_prefix_ids(selfies_str, stoi, max_len):
    bos = stoi["<BOS>"]
    mask = stoi["<MASK>"]
    toks = [bos] + [stoi.get(t, mask) for t in sf.split_selfies(selfies_str)]
    toks = toks[: min(max_len, SELFIES_PREFIX_MAX)]
    if len(toks) < 2:
        return None
    return torch.tensor(toks, dtype=torch.long)


def ids_to_selfies(ids, itos):
    toks = []
    for i in ids:
        tok = itos.get(int(i), "")
        if tok == "<EOS>":
            break
        if tok in ["<PAD>", "<BOS>", "<MASK>"]:
            continue
        toks.append(tok)
    return "".join(toks)


class TimestepEmbed(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(128, dim), nn.SiLU(), nn.Linear(dim, dim))

    def forward(self, t):
        device = t.device
        half = 64
        freqs = torch.exp(torch.arange(half, device=device) * (-np.log(10000.0) / half))
        args = t.float().unsqueeze(1) * freqs.unsqueeze(0)
        pos = torch.cat([torch.sin(args), torch.cos(args)], dim=1)
        return self.fc(pos)


class SelfiesDenoiser(nn.Module):
    def __init__(self, vocab, hidden, max_len):
        super().__init__()
        self.emb = nn.Embedding(vocab, hidden, padding_idx=0)
        self.pos = nn.Embedding(max_len, hidden)
        self.time = TimestepEmbed(hidden)
        self.gru = nn.GRU(hidden, hidden, batch_first=True, bidirectional=True)
        self.fc = nn.Linear(2 * hidden, vocab)

    def forward(self, x, t):
        b, l = x.shape
        h = self.emb(x) + self.pos(
            torch.arange(l, device=x.device).unsqueeze(0).expand(b, l)
        )
        te = self.time(t).unsqueeze(1).expand_as(h)
        y, _ = self.gru(h + te)
        return self.fc(y)


class SelfiesDiffusion:
    def __init__(self, stoi, itos, hidden, max_len, T=200):
        self.stoi = stoi
        self.itos = itos
        self.vocab = len(stoi)
        self.max_len = max_len
        self.model = SelfiesDenoiser(self.vocab, hidden, max_len).to(DEVICE)
        self.opt = torch.optim.Adam(self.model.parameters(), lr=LR)
        self.T = T
        self.pad = stoi["<PAD>"]
        self.mask = stoi["<MASK>"]
        self.eos = stoi["<EOS>"]
        self.bos = stoi["<BOS>"]

    def mask_ratio(self, t):
        return (0.05 + 0.90 * (t.float() / float(self.T))).clamp(0.05, 0.95)

    def corrupt(self, x0, t):
        b, l = x0.shape
        ratio = self.mask_ratio(t).view(-1, 1)
        valid = (x0 != self.pad) & (x0 != self.bos)
        mask_positions = (torch.rand(b, l, device=x0.device) < ratio) & valid
        for i in range(b):
            if valid[i].any() and not mask_positions[i].any():
                pos = torch.where(valid[i])[0]
                j = pos[torch.randint(0, len(pos), (1,), device=x0.device)]
                mask_positions[i, j] = True
        xt = x0.clone()
        xt[mask_positions] = self.mask
        return xt, mask_positions

    def train_epoch(self, loader):
        self.model.train()
        losses = []
        for batch in loader:
            x0 = _as_tensor(batch).to(DEVICE)
            t = torch.randint(1, self.T + 1, (x0.size(0),), device=DEVICE)
            xt, mask_positions = self.corrupt(x0, t)
            logits = self.model(xt, t)
            labels = x0.clone()
            labels[~mask_positions] = -100
            loss = F.cross_entropy(
                logits.reshape(-1, logits.size(-1)),
                labels.reshape(-1),
                ignore_index=-100,
            )
            self.opt.zero_grad()
            loss.backward()
            self.opt.step()
            losses.append(float(loss.item()))
        return float(np.mean(losses)) if losses else float("nan")

    def _sample_from_probs(self, probs, topk=None, topp=0.95):
        shape = probs.shape
        flat = probs.reshape(-1, probs.size(-1))
        if topk is not None:
            k = min(int(topk), flat.size(-1))
            sp, si = torch.topk(flat, k, dim=-1)
            sp = sp / (sp.sum(-1, keepdim=True) + 1e-12)
            choice = torch.multinomial(sp, 1)
            idx = si.gather(-1, choice).squeeze(-1)
            conf = sp.gather(-1, choice).squeeze(-1)
        elif topp is not None and topp < 1.0:
            sp, si = torch.sort(flat, dim=-1, descending=True)
            cum = sp.cumsum(-1)
            remove = cum > topp
            remove[:, 1:] = remove[:, :-1].clone()
            remove[:, 0] = False
            sp[remove] = 0.0
            sp = sp / (sp.sum(-1, keepdim=True) + 1e-12)
            choice = torch.multinomial(sp, 1)
            idx = si.gather(-1, choice).squeeze(-1)
            conf = sp.gather(-1, choice).squeeze(-1)
        else:
            choice = torch.multinomial(flat, 1)
            idx = choice.squeeze(-1)
            conf = flat.gather(-1, choice).squeeze(-1)
        return idx.view(shape[:-1]), conf.view(shape[:-1])

    @torch.no_grad()
    def sample(self, n, seed_ids=None, topk=SELFIES_TOPK, topp=SELFIES_TOPP):
        self.model.eval()
        x = torch.full((n, self.max_len), self.mask, device=DEVICE, dtype=torch.long)
        x[:, 0] = self.bos
        fixed = torch.zeros_like(x, dtype=torch.bool)
        fixed[:, 0] = True

        if seed_ids is not None:
            k = min(seed_ids.size(0), n)
            for i in range(k):
                pref = seed_ids[i].to(DEVICE)
                pad_pos = torch.where(pref == self.pad)[0]
                raw_len = int(pad_pos[0].item()) if len(pad_pos) else int(pref.numel())
                L = min(raw_len, self.max_len)
                x[i, :L] = pref[:L]
                fixed[i, :L] = True
                x[i, 0] = self.bos

        mutable = ~fixed
        mutable[:, 0] = False
        for step, tval in enumerate(range(self.T, 0, -1), start=1):
            masked = (x == self.mask) & mutable
            if not masked.any():
                break
            t = torch.full((n,), tval, device=DEVICE, dtype=torch.long)
            logits = self.model(x, t)
            probs = F.softmax(logits, dim=-1)
            probs[..., self.pad] = 0.0
            probs[..., self.bos] = 0.0
            probs[..., self.mask] = 0.0
            min_eos_pos = min(self.max_len - 1, SELFIES_MIN_CONTENT_TOK + 1)
            if min_eos_pos > 1:
                probs[:, 1:min_eos_pos, self.eos] = 0.0
            probs = probs / (probs.sum(dim=-1, keepdim=True) + 1e-12)
            proposal, confidence = self._sample_from_probs(probs, topk=topk, topp=topp)

            for i in range(n):
                pos = torch.where(masked[i])[0]
                if len(pos) == 0:
                    continue
                desired_remaining = int(
                    round((tval - 1) / float(self.T) * int(mutable[i].sum().item()))
                )
                current_remaining = len(pos)
                commit_n = max(1, current_remaining - desired_remaining)
                commit_n = min(commit_n, current_remaining)
                conf_i = confidence[i, pos]
                chosen = pos[torch.topk(conf_i, k=commit_n, largest=True).indices]
                x[i, chosen] = proposal[i, chosen]

        remaining = (x == self.mask) & mutable
        if remaining.any():
            t = torch.ones((n,), device=DEVICE, dtype=torch.long)
            probs = F.softmax(self.model(x, t), dim=-1)
            probs[..., self.pad] = 0.0
            probs[..., self.bos] = 0.0
            probs[..., self.mask] = 0.0
            probs = probs / (probs.sum(dim=-1, keepdim=True) + 1e-12)
            proposal, _ = self._sample_from_probs(probs, topk=topk, topp=topp)
            x[remaining] = proposal[remaining]

        out = []
        for row in x.tolist():
            s = ids_to_selfies(row, self.itos)
            sm = from_selfies(s)
            if sm is not None:
                out.append(sm)
        return canonicalize_and_filter_neutral(out, exclude_dm=False)


def prepare_selfies_training_data(df):
    selfies_list = [to_selfies(s) for s in df["SMILES"].tolist()]
    selfies_list = [s for s in selfies_list if s]
    stoi, itos = build_selfies_vocab(selfies_list)
    X = np.stack(
        [selfies_to_ids(s, stoi, SELFIES_MAX_TOK) for s in selfies_list], axis=0
    )
    return stoi, itos, torch.tensor(X, dtype=torch.long)


def train_selfies_generator(df):
    stoi, itos, X = prepare_selfies_training_data(df)
    loader = torch.utils.data.DataLoader(
        torch.utils.data.TensorDataset(X),
        batch_size=BATCH_SIZE,
        shuffle=True,
        drop_last=False,
        pin_memory=(DEVICE == "cuda"),
    )
    model = SelfiesDiffusion(stoi, itos, SELFIES_HIDDEN, SELFIES_MAX_TOK, T=T_STEPS)
    for ep in range(EPOCHS):
        loss = model.train_epoch(loader)
        print(f"[SELFIES-Diffusion] epoch {ep+1}/{EPOCHS} loss={loss:.4f}")
    return {"model": model, "stoi": stoi, "itos": itos}


def build_selfies_seed_ids(seed_smiles, stoi):
    if not seed_smiles:
        return None
    seed_pref = []
    seen = set()
    for sm in seed_smiles:
        cs = canonicalize_smiles(sm)
        if cs is None or cs in seen:
            continue
        seen.add(cs)
        ss = to_selfies(cs)
        if not ss:
            continue
        pref = selfies_to_prefix_ids(ss, stoi, SELFIES_MAX_TOK)
        if pref is not None:
            seed_pref.append(pref)
    if not seed_pref:
        return None
    return torch.nn.utils.rnn.pad_sequence(
        seed_pref, batch_first=True, padding_value=stoi["<PAD>"]
    )


def selfies_pipeline(bundle, df, seed_smiles=None, round_id=0):
    seed_ids = build_selfies_seed_ids(seed_smiles, bundle["stoi"])
    all_valid = []
    seen = set()
    attempts = 0
    while len(all_valid) < MIN_VALID and attempts < 10:
        batch = bundle["model"].sample(NUM_GEN, seed_ids=seed_ids)
        _extend_unique_generated(all_valid, seen, batch, exclude_dm=True)
        attempts += 1
    if len(all_valid) < MIN_VALID:
        all_valid = brics_augment(df, need=MIN_VALID, base=all_valid)
    all_valid = finalize_generated(all_valid, f"SELFIES_R{round_id}")
    out_csv = OUT_DIR / f"generated_selfies_diffusion_R{round_id}.csv"
    pd.DataFrame({"smiles": all_valid}).to_csv(out_csv, index=False)
    print(f"[SELFIES-Diffusion] valid={len(all_valid)} -> {out_csv}")
    return out_csv


def bond_to_code(bond):
    bt = bond.GetBondType()
    if bt == Chem.BondType.SINGLE:
        return BOND_SINGLE
    if bt == Chem.BondType.DOUBLE:
        return BOND_DOUBLE
    if bt == Chem.BondType.TRIPLE:
        return BOND_TRIPLE
    if bt == Chem.BondType.AROMATIC:
        return BOND_AROMATIC
    return None


def mol_to_graph(smiles):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None, None
    atoms = list(m.GetAtoms())
    if len(atoms) > GRAPH_N_MAX:
        return None, None
    X = np.zeros((GRAPH_N_MAX, GRAPH_NODE_DIM), dtype=np.float32)
    X[:, 0] = 1.0
    for i, atom in enumerate(atoms):
        sym = atom.GetSymbol()
        idx = GRAPH_ATOM_TO_INDEX.get(sym)
        if idx is None:
            return None, None
        X[i, :] = 0.0
        X[i, idx] = 1.0
    A = np.zeros((GRAPH_N_MAX, GRAPH_N_MAX), dtype=np.float32)
    for bond in m.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        code = bond_to_code(bond)
        if code is None:
            return None, None
        A[i, j] = code / BOND_CODE_MAX
        A[j, i] = code / BOND_CODE_MAX
    return X, A


def graph_to_mol(X, A):
    try:
        node_idx = np.argmax(X, axis=1)
        active_positions = [i for i, t in enumerate(node_idx) if t != 0]
        if len(active_positions) < 2:
            return None
        rw = Chem.RWMol()
        pos_to_new = {}
        for pos in active_positions:
            sym = GRAPH_TYPES[int(node_idx[pos])]
            if sym == "<PAD>" or sym not in GRAPH_ATOM_TO_INDEX:
                continue
            pos_to_new[pos] = rw.AddAtom(Chem.Atom(sym))
        for ii, i in enumerate(active_positions):
            if i not in pos_to_new:
                continue
            for j in active_positions[ii + 1 :]:
                if j not in pos_to_new:
                    continue
                code = int(np.rint(np.clip(A[i, j], 0.0, 1.0) * BOND_CODE_MAX))
                if code == BOND_NONE:
                    continue
                if code == BOND_SINGLE:
                    bt = Chem.BondType.SINGLE
                elif code == BOND_DOUBLE:
                    bt = Chem.BondType.DOUBLE
                elif code == BOND_TRIPLE:
                    bt = Chem.BondType.TRIPLE
                elif code == BOND_AROMATIC:
                    bt = Chem.BondType.AROMATIC
                else:
                    continue
                try:
                    rw.AddBond(pos_to_new[i], pos_to_new[j], bt)
                    if bt == Chem.BondType.AROMATIC:
                        rw.GetAtomWithIdx(pos_to_new[i]).SetIsAromatic(True)
                        rw.GetAtomWithIdx(pos_to_new[j]).SetIsAromatic(True)
                except Exception:
                    pass
        m = rw.GetMol()
        Chem.SanitizeMol(m)
        sm = Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)
        return sm if Chem.MolFromSmiles(sm) is not None else None
    except Exception:
        return None


class TimeMLP(nn.Module):
    def __init__(self, dim):
        super().__init__()
        self.fc = nn.Sequential(nn.Linear(128, dim), nn.SiLU(), nn.Linear(dim, dim))

    def forward(self, t):
        device = t.device
        half = 64
        freqs = torch.exp(torch.arange(half, device=device) * (-np.log(10000.0) / half))
        args = t.float().unsqueeze(1) * freqs.unsqueeze(0)
        pos = torch.cat([torch.sin(args), torch.cos(args)], dim=1)
        return self.fc(pos)


class GraphDenoiser(nn.Module):
    def __init__(self, n_max, node_dim, hidden):
        super().__init__()
        self.n_max = n_max
        self.node_dim = node_dim
        self.time = TimeMLP(hidden)
        in_dim = n_max * node_dim + n_max * n_max
        self.pre = nn.Linear(in_dim, hidden)
        self.gnn = nn.Sequential(nn.SiLU(), nn.Linear(hidden, hidden), nn.SiLU())
        self.out_x = nn.Linear(hidden, n_max * node_dim)
        self.out_a = nn.Linear(hidden, n_max * n_max)

    def forward(self, X, A, t):
        b = X.size(0)
        h = torch.cat([X.reshape(b, -1), A.reshape(b, -1)], dim=1)
        h = self.pre(h)
        h = self.gnn(h + self.time(t))
        ex = self.out_x(h).view(b, self.n_max, self.node_dim)
        ea = self.out_a(h).view(b, self.n_max, self.n_max)
        ea = (ea + ea.transpose(1, 2)) / 2.0
        return ex, ea


class GraphDiffusion:
    def __init__(self, n_max, node_dim):
        self.n_max = n_max
        self.node_dim = node_dim
        self.model = GraphDenoiser(n_max, node_dim, GRAPH_HIDDEN).to(DEVICE)
        self.opt = torch.optim.Adam(self.model.parameters(), lr=LR)
        self.T = T_STEPS

    def sigma(self, t):
        return GRAPH_SIGMA_MIN + (GRAPH_SIGMA_MAX - GRAPH_SIGMA_MIN) * (
            t.float() / float(self.T)
        )

    def train_epoch(self, loader):
        self.model.train()
        losses = []
        for X0, A0 in loader:
            X0 = X0.to(DEVICE, non_blocking=True)
            A0 = A0.to(DEVICE, non_blocking=True)
            t = torch.randint(1, self.T + 1, (X0.size(0),), device=DEVICE)
            sigma = self.sigma(t).view(-1, 1, 1)
            Xt = X0 + sigma * torch.randn_like(X0)
            At = A0 + sigma * torch.randn_like(A0)
            ex, ea = self.model(Xt, At, t)
            loss = ((ex - X0) ** 2).mean() + ((ea - A0) ** 2).mean()
            self.opt.zero_grad()
            loss.backward()
            self.opt.step()
            losses.append(float(loss.item()))
        return float(np.mean(losses)) if losses else float("nan")

    @torch.no_grad()
    def sample(self, n, steps=GRAPH_STEPS, seed_graphs=None):
        self.model.eval()
        t_values = np.unique(np.linspace(self.T, 1, steps, dtype=int))[::-1]
        sigma_T = float(GRAPH_SIGMA_MAX)
        X = torch.randn(n, self.n_max, self.node_dim, device=DEVICE) * sigma_T
        A = torch.randn(n, self.n_max, self.n_max, device=DEVICE) * sigma_T
        A = (A + A.transpose(1, 2)) / 2.0

        if seed_graphs is not None:
            Xs, As = seed_graphs
            k = min(Xs.shape[0], n)
            X[:k] = Xs[:k].to(DEVICE) + sigma_T * torch.randn_like(X[:k])
            A[:k] = As[:k].to(DEVICE) + sigma_T * torch.randn_like(A[:k])

        for idx, tval in enumerate(t_values):
            tt = torch.full((n,), int(tval), device=DEVICE, dtype=torch.long)
            x0_hat, a0_hat = self.model(X, A, tt)
            a0_hat = (a0_hat + a0_hat.transpose(1, 2)) / 2.0
            if idx == len(t_values) - 1:
                X = x0_hat
                A = a0_hat
            else:
                next_t = int(t_values[idx + 1])
                next_sigma = float(
                    GRAPH_SIGMA_MIN
                    + (GRAPH_SIGMA_MAX - GRAPH_SIGMA_MIN) * (next_t / float(self.T))
                )
                X = x0_hat + next_sigma * torch.randn_like(X)
                A = a0_hat + next_sigma * torch.randn_like(A)
                A = (A + A.transpose(1, 2)) / 2.0

        X = F.softmax(X, dim=-1).cpu().numpy()
        A = torch.clamp(A, 0.0, 1.0).cpu().numpy()
        out = []
        for i in range(n):
            sm = graph_to_mol(X[i], A[i])
            if sm:
                out.append(sm)
        return list(dict.fromkeys(out))


def prepare_graph_dataset(df):
    Xs, As = [], []
    skipped = 0
    for sm in df["SMILES"].tolist():
        x, a = mol_to_graph(sm)
        if x is None:
            skipped += 1
            continue
        Xs.append(x)
        As.append(a)
    if not Xs:
        raise RuntimeError("Empty Graph dataset")
    if skipped:
        print(f"[Graph] skipped unsupported molecules: {skipped}")
    X = torch.tensor(np.stack(Xs), dtype=torch.float32)
    A = torch.tensor(np.stack(As), dtype=torch.float32)
    loader = torch.utils.data.DataLoader(
        torch.utils.data.TensorDataset(X, A),
        batch_size=BATCH_SIZE,
        shuffle=True,
        drop_last=False,
        pin_memory=(DEVICE == "cuda"),
    )
    return loader


def train_graph_generator(df):
    loader = prepare_graph_dataset(df)
    gd = GraphDiffusion(GRAPH_N_MAX, GRAPH_NODE_DIM)
    for ep in range(EPOCHS):
        loss = gd.train_epoch(loader)
        print(f"[Graph-Diffusion] epoch {ep+1}/{EPOCHS} loss={loss:.6f}")
    return gd


def build_graph_seeds(seed_smiles):
    if not seed_smiles:
        return None
    Xs, As, seen = [], [], set()
    for s in seed_smiles:
        cs = canonicalize_smiles(s)
        if cs is None or cs in seen:
            continue
        seen.add(cs)
        x, a = mol_to_graph(cs)
        if x is not None:
            Xs.append(x)
            As.append(a)
    if not Xs:
        return None
    return torch.tensor(np.stack(Xs), dtype=torch.float32), torch.tensor(
        np.stack(As), dtype=torch.float32
    )


def graph_pipeline(model, df, seed_smiles=None, round_id=0):
    seed_graphs = build_graph_seeds(seed_smiles)
    all_valid, seen = [], set()
    attempts = 0
    while len(all_valid) < MIN_VALID and attempts < 20:
        batch = model.sample(NUM_GEN, seed_graphs=seed_graphs)
        _extend_unique_generated(all_valid, seen, batch, exclude_dm=True)
        attempts += 1
    if len(all_valid) < MIN_VALID:
        all_valid = brics_augment(df, need=MIN_VALID, base=all_valid)
    all_valid = finalize_generated(all_valid, f"GRAPH_R{round_id}")
    out_csv = OUT_DIR / f"generated_graph_diffusion_R{round_id}.csv"
    pd.DataFrame({"smiles": all_valid}).to_csv(out_csv, index=False)
    print(f"[Graph-Diffusion] valid={len(all_valid)} -> {out_csv}")
    return out_csv


def _build_geodiff_feature(sm, max_iters=100):
    m = Chem.MolFromSmiles(sm)
    if m is None or m.GetNumHeavyAtoms() > GEODIFF_N_MAX:
        return None
    mh = Chem.AddHs(m)
    try:
        params = AllChem.ETKDGv3()
        params.randomSeed = SEED
        params.useRandomCoords = True
        params.maxIterations = 1000
        status = AllChem.EmbedMolecule(mh, params)
        if status != 0:
            status = AllChem.EmbedMolecule(
                mh, randomSeed=SEED, useRandomCoords=True, maxAttempts=1000
            )
        if status != 0:
            return None
        try:
            AllChem.UFFOptimizeMolecule(mh, maxIters=max_iters)
        except Exception:
            pass
        m3 = Chem.RemoveHs(mh)
    except Exception:
        return None
    conf = m3.GetConformer()
    atoms = list(m3.GetAtoms())
    n = len(atoms)
    C = np.zeros((GEODIFF_N_MAX, 3), dtype=np.float32)
    T = np.zeros((GEODIFF_N_MAX, GEODIFF_NODE_DIM), dtype=np.float32)
    T[:, 0] = 1.0
    coords = []
    for i, atom in enumerate(atoms):
        sym = atom.GetSymbol()
        idx = GEODIFF_ATOM_TO_INDEX.get(sym)
        if idx is None:
            return None
        p = conf.GetAtomPosition(i)
        coords.append([p.x, p.y, p.z])
        T[i, :] = 0.0
        T[i, idx] = 1.0
    coords = np.asarray(coords, dtype=np.float32)
    coords -= coords.mean(axis=0, keepdims=True)
    C[:n] = coords
    return C, T


class GeoDenoiser(nn.Module):
    def __init__(self, n_max, type_dim, hidden):
        super().__init__()
        self.n_max = n_max
        self.type_dim = type_dim
        self.time = TimeMLP(hidden)
        in_dim = n_max * 3 + n_max * type_dim
        self.pre = nn.Linear(in_dim, hidden)
        self.encoder = nn.Sequential(nn.SiLU(), nn.Linear(hidden, hidden), nn.SiLU())
        self.out_c = nn.Linear(hidden, n_max * 3)
        self.out_t = nn.Linear(hidden, n_max * type_dim)

    def forward(self, coords, types, t):
        b = coords.size(0)
        h = torch.cat([coords.reshape(b, -1), types.reshape(b, -1)], dim=1)
        h = self.pre(h)
        h = self.encoder(h + self.time(t))
        c0 = self.out_c(h).view(b, self.n_max, 3)
        t0 = self.out_t(h).view(b, self.n_max, self.type_dim)
        return c0, t0


class GeoDiffusion:
    def __init__(self, n_max, type_dim):
        self.n_max = n_max
        self.type_dim = type_dim
        self.model = GeoDenoiser(n_max, type_dim, GEODIFF_HIDDEN).to(DEVICE)
        self.opt = torch.optim.Adam(self.model.parameters(), lr=LR)
        self.T = T_STEPS

    def sigma(self, t):
        return GEODIFF_SIGMA_MIN + (GEODIFF_SIGMA_MAX - GEODIFF_SIGMA_MIN) * (
            t.float() / float(self.T)
        )

    def train_epoch(self, loader):
        self.model.train()
        losses = []
        for C0, T0 in loader:
            C0 = C0.to(DEVICE, non_blocking=True)
            T0 = T0.to(DEVICE, non_blocking=True)
            t = torch.randint(1, self.T + 1, (C0.size(0),), device=DEVICE)
            sigma = self.sigma(t).view(-1, 1, 1)
            Ct = C0 + sigma * torch.randn_like(C0)
            Tt = T0 + sigma * torch.randn_like(T0)
            c0_hat, t0_hat = self.model(Ct, Tt, t)
            loss = ((c0_hat - C0) ** 2).mean() + ((t0_hat - T0) ** 2).mean()
            self.opt.zero_grad()
            loss.backward()
            self.opt.step()
            losses.append(float(loss.item()))
        return float(np.mean(losses)) if losses else float("nan")

    @torch.no_grad()
    def sample(self, n, steps=GEODIFF_STEPS, seed=None):
        self.model.eval()
        t_values = np.unique(np.linspace(self.T, 1, steps, dtype=int))[::-1]
        sigma_T = float(GEODIFF_SIGMA_MAX)
        C = torch.randn(n, self.n_max, 3, device=DEVICE) * sigma_T
        T = torch.randn(n, self.n_max, self.type_dim, device=DEVICE) * sigma_T

        if seed is not None:
            Cs, Ts = seed
            k = min(Cs.shape[0], n)
            C[:k] = Cs[:k].to(DEVICE) + sigma_T * torch.randn_like(C[:k])
            T[:k] = Ts[:k].to(DEVICE) + sigma_T * torch.randn_like(T[:k])

        for idx, tval in enumerate(t_values):
            tt = torch.full((n,), int(tval), device=DEVICE, dtype=torch.long)
            c0_hat, t0_hat = self.model(C, T, tt)
            if idx == len(t_values) - 1:
                C = c0_hat
                T = t0_hat
            else:
                next_t = int(t_values[idx + 1])
                next_sigma = float(
                    GEODIFF_SIGMA_MIN
                    + (GEODIFF_SIGMA_MAX - GEODIFF_SIGMA_MIN) * (next_t / float(self.T))
                )
                C = c0_hat + next_sigma * torch.randn_like(C)
                T = t0_hat + next_sigma * torch.randn_like(T)
        return C, F.softmax(T, dim=-1)


def coords_to_smiles(C, T):
    try:
        types = np.argmax(T, axis=1)
        active = [i for i, t in enumerate(types) if int(t) != 0]
        if len(active) < 2:
            return None
        rw = Chem.RWMol()
        pos_to_new = {}
        atomic_nums = {}
        for pos in active:
            sym = GEODIFF_TYPES[int(types[pos])]
            if sym == "<PAD>" or sym not in GEODIFF_ATOM_TO_INDEX:
                continue
            idx = rw.AddAtom(Chem.Atom(sym))
            pos_to_new[pos] = idx
            atomic_nums[pos] = Chem.GetPeriodicTable().GetAtomicNumber(sym)
        if len(pos_to_new) < 2:
            return None

        pt = Chem.GetPeriodicTable()
        for ii, i in enumerate(active):
            if i not in pos_to_new:
                continue
            for j in active[ii + 1 :]:
                if j not in pos_to_new:
                    continue
                d = float(np.linalg.norm(C[i] - C[j]))
                ri = float(pt.GetRcovalent(atomic_nums[i]))
                rj = float(pt.GetRcovalent(atomic_nums[j]))
                cutoff = 1.25 * (ri + rj)
                if d <= 0.4 or d > cutoff:
                    continue
                ratio = d / max(ri + rj, 1e-6)
                if ratio < 0.72:
                    bt = Chem.BondType.TRIPLE
                elif ratio < 0.86:
                    bt = Chem.BondType.DOUBLE
                else:
                    bt = Chem.BondType.SINGLE
                try:
                    rw.AddBond(pos_to_new[i], pos_to_new[j], bt)
                except Exception:
                    pass
        m = rw.GetMol()
        Chem.SanitizeMol(m)
        sm = Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)
        return sm if Chem.MolFromSmiles(sm) is not None else None
    except Exception:
        return None


def prepare_geodiff_dataset(df):
    Cs, Ts = [], []
    skipped = 0
    for sm in df["SMILES"].tolist():
        item = _build_geodiff_feature(sm, 100)
        if item is None:
            skipped += 1
            continue
        C, T = item
        Cs.append(C)
        Ts.append(T)
    if not Cs:
        raise RuntimeError("Empty GeoDiff dataset")
    if skipped:
        print(f"[GeoDiff] skipped unsupported molecules: {skipped}")
    C = torch.tensor(np.stack(Cs), dtype=torch.float32)
    T = torch.tensor(np.stack(Ts), dtype=torch.float32)
    loader = torch.utils.data.DataLoader(
        torch.utils.data.TensorDataset(C, T),
        batch_size=BATCH_SIZE,
        shuffle=True,
        drop_last=False,
        pin_memory=(DEVICE == "cuda"),
    )
    return loader


def train_geodiff_generator(df):
    loader = prepare_geodiff_dataset(df)
    gd = GeoDiffusion(GEODIFF_N_MAX, GEODIFF_NODE_DIM)
    for ep in range(EPOCHS):
        loss = gd.train_epoch(loader)
        print(f"[GeoDiff] epoch {ep+1}/{EPOCHS} loss={loss:.6f}")
    return gd


def build_geodiff_seeds(seed_smiles):
    if not seed_smiles:
        return None
    Cs, Ts, seen = [], [], set()
    for s in seed_smiles:
        cs = canonicalize_smiles(s)
        if cs is None or cs in seen:
            continue
        seen.add(cs)
        item = _build_geodiff_feature(cs, 50)
        if item is None:
            continue
        C, T = item
        Cs.append(C)
        Ts.append(T)
    if not Cs:
        return None
    return torch.tensor(np.stack(Cs), dtype=torch.float32), torch.tensor(
        np.stack(Ts), dtype=torch.float32
    )


def geodiff_pipeline(model, df, seed_smiles=None, round_id=0):
    seed = build_geodiff_seeds(seed_smiles)
    all_valid, seen = [], set()
    tries = 0
    while len(all_valid) < MIN_VALID and tries < 30:
        C, T = model.sample(NUM_GEN, seed=seed)
        C = C.cpu().numpy()
        T = T.cpu().numpy()
        batch = []
        for i in range(C.shape[0]):
            sm = coords_to_smiles(C[i], T[i])
            if sm:
                batch.append(sm)
        _extend_unique_generated(all_valid, seen, batch, exclude_dm=True)
        tries += 1
    if len(all_valid) < MIN_VALID:
        all_valid = brics_augment(df, need=MIN_VALID, base=all_valid)
    all_valid = finalize_generated(all_valid, f"GEODIFF_R{round_id}")
    out_csv = OUT_DIR / f"generated_geodiff_R{round_id}.csv"
    pd.DataFrame({"smiles": all_valid}).to_csv(out_csv, index=False)
    print(f"[GeoDiff] valid={len(all_valid)} -> {out_csv}")
    return out_csv


def brics_augment(df, need=None, seed=SEED, base=None):
    if need is None:
        need = MIN_VALID
    rng = np.random.RandomState(seed)
    base = finalize_generated(base[:] if base else [], "BRICS-base")
    seen = {structure_key(x) for x in base if structure_key(x) is not None}
    mols = [Chem.MolFromSmiles(s) for s in df["SMILES"].tolist()]
    mols = [m for m in mols if m is not None]
    fr = set()
    for m in mols:
        try:
            fr.update(BRICS.BRICSDecompose(m))
        except Exception:
            pass
    fr = list(fr)
    tries = 0
    while len(base) < need and tries < 5000 and len(fr) >= 2:
        a, b = rng.choice(fr, 2, replace=True)
        try:
            for m in BRICS.BRICSBuild([a, b]):
                try:
                    sm = Chem.MolToSmiles(m, canonical=True, isomericSmiles=True)
                except Exception:
                    continue
                _extend_unique_generated(base, seen, [sm], exclude_dm=True)
                if len(base) >= need:
                    break
        except Exception:
            pass
        tries += 1
    return base[:need]


def dataset_random_enum(df, need=None, seed=SEED, base=None):
    if need is None:
        need = MIN_VALID
    base = finalize_generated(base[:] if base else [], "ENUM-base")
    print(
        "[INFO] dataset_random_enum skipped because it only reproduces Dm structures."
    )
    return base[:need]


def extract_top_seeds(pred_csv, tag="MERGED"):
    if pred_csv is None or not os.path.exists(pred_csv):
        return []
    dfp = pd.read_csv(pred_csv)
    scol = (
        "smiles"
        if "smiles" in dfp.columns
        else ("SMILES" if "SMILES" in dfp.columns else None)
    )
    if scol is None:
        return []
    pcol = "D_pred" if "D_pred" in dfp.columns else pick_pred_column(dfp)
    if pcol is None:
        return []
    dfp[pcol] = pd.to_numeric(dfp[pcol], errors="coerce")
    dfp = dfp.dropna(subset=[pcol])
    ascending = TARGET_DIRECTION_RESOLVED == "min"
    dfp = dfp.sort_values(by=pcol, ascending=ascending)
    topk = _auto_seed_count(len(dfp))
    print(f"[AUTO SEEDS:{tag}] candidates={len(dfp)} selected={topk}")
    return dfp[scol].head(topk).astype(str).tolist()


def enrich_by_predictions(ckpt, gen_csv, tag):
    preds, hits = chemprop_predict_paths(ckpt, gen_csv, tag)
    seeds = extract_top_seeds(preds, tag=tag) if preds else []
    return preds, seeds


def merge_and_dedup(csv_list, out_name):
    uniq = []
    seen = set()
    total = 0
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
            if col is None:
                continue
            values = df[col].dropna().astype(str).tolist()
            total += len(values)
            _extend_unique_generated(uniq, seen, values, exclude_dm=True)
        except Exception as e:
            print(f"[WARN] merge failed for {p}: {e}")
    out = OUT_DIR / out_name
    pd.DataFrame({"smiles": uniq}).to_csv(out, index=False)
    print(
        f"[FILTER:MERGE] input={total} kept={len(uniq)} removed={total-len(uniq)} -> {out}"
    )
    return out


def main():
    initialize_excluded_elements()
    df = load_dm_dataframe(DATA_CSV)
    resolve_target_direction(df)
    initialize_chemical_space(df)

    print("[TRAIN] Chemprop fixed 5-fold models")
    ckpt = chemprop_kfold_train(df, 5)

    print("[TRAIN] SELFIES generator once")
    selfies_model = train_selfies_generator(df)
    print("[TRAIN] Graph generator once")
    graph_model = train_graph_generator(df)
    print("[TRAIN] GeoDiff generator once")
    geodiff_model = train_geodiff_generator(df)

    csv_self = selfies_pipeline(selfies_model, df, round_id=0)
    pred_self, seeds_self = enrich_by_predictions(
        ckpt, csv_self, "SELFIES_DIFFUSION_R0"
    )

    csv_graph = graph_pipeline(graph_model, df, seed_smiles=seeds_self, round_id=0)
    pred_graph, seeds_graph = enrich_by_predictions(
        ckpt, csv_graph, "GRAPH_DIFFUSION_R0"
    )

    seeds_all = (seeds_self or []) + (seeds_graph or [])
    csv_geo = geodiff_pipeline(geodiff_model, df, seed_smiles=seeds_all, round_id=0)
    pred_geo, seeds_geo = enrich_by_predictions(ckpt, csv_geo, "GEODIFF_R0")

    merged_prev = merge_and_dedup(
        [csv_self, csv_graph, csv_geo],
        "generated_merged_round0.csv",
    )
    merged_pred, _ = chemprop_predict_paths(ckpt, merged_prev, "MERGED_R0")
    seeds = extract_top_seeds(merged_pred, tag="MERGED_R0")

    for r in range(1, ENRICH_ROUNDS + 1):
        csv_self = selfies_pipeline(selfies_model, df, seed_smiles=seeds, round_id=r)
        pred_self, seeds_self = enrich_by_predictions(
            ckpt, csv_self, f"SELFIES_DIFFUSION_R{r}"
        )

        csv_graph = graph_pipeline(graph_model, df, seed_smiles=seeds, round_id=r)
        pred_graph, seeds_graph = enrich_by_predictions(
            ckpt, csv_graph, f"GRAPH_DIFFUSION_R{r}"
        )

        csv_geo = geodiff_pipeline(geodiff_model, df, seed_smiles=seeds, round_id=r)
        pred_geo, seeds_geo = enrich_by_predictions(ckpt, csv_geo, f"GEODIFF_R{r}")

        merged = merge_and_dedup(
            [merged_prev, csv_self, csv_graph, csv_geo],
            f"generated_merged_round{r}.csv",
        )
        merged_pred, _ = chemprop_predict_paths(ckpt, merged, f"MERGED_R{r}")
        seeds = extract_top_seeds(merged_pred, tag=f"MERGED_R{r}")
        merged_prev = merged

    print("Done.")


if __name__ == "__main__":
    main()

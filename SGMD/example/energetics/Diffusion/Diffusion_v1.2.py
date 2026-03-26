import os, sys, re, random, shutil, subprocess, glob

os.chdir(os.path.split(os.path.realpath(__file__))[0])
from pathlib import Path
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
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
DATA_CSV = "D.csv"
BASE_DIR = Path("runs/mol_chemprop_multi")
CP_DIR = BASE_DIR / "models_chemprop"
OUT_DIR = BASE_DIR / "outputs"
for p in [CP_DIR, OUT_DIR]:
    p.mkdir(parents=True, exist_ok=True)
TARGET_D = 9.5
T_STEPS = 200
BATCH_SIZE = 64
EPOCHS = 30
LR = 1e-3
NUM_GEN = 2048
MIN_VALID = 256
SELFIES_MAX_TOK = 120
SELFIES_PREFIX_MAX = 30
SELFIES_HIDDEN = 384
SELFIES_TOPK = None
SELFIES_TOPP = 0.95
GRAPH_N_MAX = 60
GRAPH_ATOMS = ["C", "H", "N", "O"]
GRAPH_NODE_DIM = len(GRAPH_ATOMS)
GRAPH_HIDDEN = 512
GRAPH_SIGMA_MIN = 0.01
GRAPH_SIGMA_MAX = 1.0
GRAPH_STEPS = 60
GEODIFF_ATOMS = GRAPH_ATOMS
GEODIFF_N_MAX = 60
GEODIFF_HIDDEN = 512
GEODIFF_STEPS = 60
ENRICH_ROUNDS = 5
TOPK_SEED = 256


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


def load_dm_dataframe(path):
    df = pd.read_csv(path)
    smiles_col = None
    d_col = None
    for c in df.columns:
        if re.search("smiles", c, re.I):
            smiles_col = c
        if re.search(r"\bD\b|velocity|target|label", c, re.I):
            d_col = c
    if smiles_col is None or d_col is None:
        raise ValueError("Missing SMILES or target column")
    df = (
        df[[smiles_col, d_col]]
        .rename(columns={smiles_col: "SMILES", d_col: "D"})
        .dropna()
        .reset_index(drop=True)
    )
    return df


def canonicalize_smiles(smiles: str):
    try:
        s = str(smiles).strip()
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            return None
        return Chem.MolToSmiles(mol, canonical=True)
    except Exception:
        return None


def is_neutral_smiles(smiles: str) -> bool:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        return Chem.rdmolops.GetFormalCharge(mol) == 0
    except Exception:
        return False


def canonicalize_and_filter_neutral(smiles_list):
    out = []
    for s in smiles_list:
        cs = canonicalize_smiles(s)
        if cs is None:
            continue
        if not is_neutral_smiles(cs):
            continue
        out.append(cs)
    return list(dict.fromkeys(out))


def simple_kfold_indices(n_splits, n_samples, seed=42):
    rng = np.random.RandomState(seed)
    idx = np.arange(n_samples)
    rng.shuffle(idx)
    folds = []
    sizes = [n_samples // n_splits] * n_splits
    for i in range(n_samples % n_splits):
        sizes[i] += 1
    start = 0
    for s in sizes:
        folds.append(idx[start : start + s])
        start += s
    return folds


def find_fold_checkpoints(save_dir, n_folds=5):
    ckpts = []
    for k in range(n_folds):
        patterns = [
            str((save_dir / f"fold_{k}/model_0/model.pt").as_posix()),
            str((save_dir / f"fold_{k}/model.pt").as_posix()),
        ]
        hit = None
        for p in patterns:
            if os.path.exists(p):
                hit = p
                break
        if hit is None:
            cands = glob.glob(
                str((save_dir / f"fold_{k}/**/*.pt").as_posix()), recursive=True
            )
            hit = cands[0] if cands else None
        ckpts.append(hit)
    return ckpts


def mae(y_true, y_pred):
    y_true = np.array(y_true, dtype=float)
    y_pred = np.array(y_pred, dtype=float)
    return float(np.mean(np.abs(y_true - y_pred)))


def rmse(y_true, y_pred):
    y_true = np.array(y_true, dtype=float)
    y_pred = np.array(y_pred, dtype=float)
    return float(np.sqrt(np.mean((y_true - y_pred) ** 2)))


def r2(y_true, y_pred):
    y_true = np.array(y_true, dtype=float)
    y_pred = np.array(y_pred, dtype=float)
    sse = np.sum((y_true - y_pred) ** 2)
    sst = np.sum((y_true - np.mean(y_true)) ** 2)
    return float(1.0 - sse / sst) if sst > 0 else float("nan")


def pick_pred_column(df: pd.DataFrame):
    lower_cols = {c.lower(): c for c in df.columns}
    for key in ["d_pred", "pred", "prediction", "y_pred", "d"]:
        if key in lower_cols:
            return lower_cols[key]
    num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    if not num_cols:
        return None
    return num_cols[0]


def kfold_eval_and_artifacts(df, save_dir, outs_dir, n_folds=5):
    folds = simple_kfold_indices(n_folds, len(df), seed=42)
    ckpts = find_fold_checkpoints(save_dir, n_folds)
    cols = {}
    per_fold = []
    colors = [plt.get_cmap("tab10")(i) for i in range(n_folds)]
    fig, ax = plt.subplots(figsize=(6, 6), dpi=300)
    xlim = [float("inf"), -float("inf")]
    ylim = [float("inf"), -float("inf")]
    for k in range(n_folds):
        val_idx = sorted(list(folds[k]))
        val_df = df.iloc[val_idx].reset_index(drop=True)
        val_csv = str((outs_dir / f"val_fold{k}.csv").as_posix())
        val_df.rename(columns={"SMILES": "smiles"}).to_csv(val_csv, index=False)
        pred_csv = str((outs_dir / f"preds_val_fold{k}.csv").as_posix())
        ck = ckpts[k]
        if ck and os.path.exists(ck):
            run(
                [
                    CHEMPROP_PREDICT,
                    "--test_path",
                    val_csv,
                    "--checkpoint_path",
                    ck,
                    "--preds_path",
                    pred_csv,
                ]
            )
        else:
            run(
                [
                    CHEMPROP_PREDICT,
                    "--test_path",
                    val_csv,
                    "--checkpoint_dir",
                    str(save_dir),
                    "--preds_path",
                    pred_csv,
                ]
            )
        pr = pd.read_csv(pred_csv)
        y_true = val_df["D"].values.tolist()
        pred_col = pick_pred_column(pr)
        if pred_col is None:
            raise ValueError(f"No numeric prediction column found in {pred_csv}")
        y_pred = pr[pred_col].astype(float).values.tolist()
        m_mae = mae(y_true, y_pred)
        m_rmse = rmse(y_true, y_pred)
        m_r2 = r2(y_true, y_pred)
        per_fold.append((k, m_mae, m_rmse, m_r2, len(y_true)))
        cols[f"fold{k}_true"] = y_true
        cols[f"fold{k}_pred"] = y_pred
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
    ax.set_title("K-fold validation predictions")
    ax.legend(frameon=False, markerscale=1.2, fontsize=8)
    ax.set_xlim(low, high)
    ax.set_ylim(low, high)
    fig.tight_layout()
    fig_path = outs_dir / "kfold_val_pred_300dpi.jpg"
    fig.savefig(fig_path, dpi=300)
    plt.close(fig)
    max_len = max(len(v) for v in cols.values())
    dat = pd.DataFrame(
        {k: (v + ([np.nan] * (max_len - len(v)))) for k, v in cols.items()}
    )
    dat_path = outs_dir / "kfold_val_preds.dat"
    dat.to_csv(dat_path, sep="\t", index=False)
    maes = [p[1] for p in per_fold]
    rmses = [p[2] for p in per_fold]
    r2s = [p[3] for p in per_fold]
    mlog = outs_dir / "verbose.log"
    with open(mlog, "w", encoding="utf-8") as f:
        f.write("5-fold cross validation\n")
        for k, mma, mmr, rrr, n in per_fold:
            f.write(
                f"\tFold {k}: n={n}, MAE = {mma:.6f}, RMSE = {mmr:.6f}, R2 = {rrr:.6f}\n"
            )
        f.write(f"Overall MAE = {np.mean(maes):.6f} +/- {np.std(maes,ddof=1):.6f}\n")
        f.write(f"Overall RMSE = {np.mean(rmses):.6f} +/- {np.std(rmses,ddof=1):.6f}\n")
        f.write(f"Overall R2 = {np.mean(r2s):.6f} +/- {np.std(r2s,ddof=1):.6f}\n")
    print(f"[KFold] metrics log -> {mlog.as_posix()}")
    print(f"[KFold] dat -> {dat_path.as_posix()}")
    print(f"[KFold] figure -> {fig_path.as_posix()}")


def chemprop_kfold_train():
    df = load_dm_dataframe(DATA_CSV)
    work_csv = CP_DIR / "train_all.csv"
    df[["SMILES", "D"]].rename(columns={"SMILES": "smiles"}).to_csv(
        work_csv, index=False
    )
    run(
        [
            CHEMPROP_TRAIN,
            "--data_path",
            str(work_csv),
            "--dataset_type",
            "regression",
            "--target_columns",
            "D",
            "--save_dir",
            str(CP_DIR),
            "--num_folds",
            "5",
            "--epochs",
            "200",
            "--batch_size",
            "64",
            "--metric",
            "rmse",
            "--seed",
            "42",
            "--quiet",
        ]
    )
    try:
        kfold_eval_and_artifacts(df, CP_DIR, OUT_DIR, 5)
    except Exception as e:
        print("[WARN] kfold post-analysis failed:", e)
    return CP_DIR


def chemprop_predict_paths(ckpt_dir, gen_csv, tag):
    out_csv = str((OUT_DIR / f"generated_with_pred_{tag}.csv").as_posix())
    hits_csv = str((OUT_DIR / f"generated_hits_{tag}.csv").as_posix())
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
    if not os.path.exists(out_csv) or os.path.getsize(out_csv) == 0:
        return None, None
    df = pd.read_csv(out_csv)
    if "smiles" not in df.columns:
        if "SMILES" in df.columns:
            df = df.rename(columns={"SMILES": "smiles"})
        elif "smiles_0" in df.columns:
            df = df.rename(columns={"smiles_0": "smiles"})
    pred_col = pick_pred_column(df)
    if pred_col is None:
        return out_csv, None
    df["D_pred"] = df[pred_col].astype(float)
    hits = df[df["D_pred"] >= float(TARGET_D)].copy()
    hits.to_csv(hits_csv, index=False)
    print(f"[OK:{tag}] pred={len(df)} hits={len(hits)}")
    return out_csv, hits_csv


def to_selfies(sm):
    try:
        cs = canonicalize_smiles(sm)
        if cs is None:
            return None
        return sf.encoder(cs)
    except Exception:
        return None


def from_selfies(s):
    try:
        sm = sf.decoder(s)
        cs = canonicalize_smiles(sm)
        return cs
    except Exception:
        return None


def build_selfies_vocab(selfies_list):
    sym = set()
    for s in selfies_list:
        try:
            for tok in sf.split_selfies(s):
                sym.add(tok)
        except:
            pass
    itoks = ["<PAD>", "<BOS>", "<EOS>", "<MASK>"] + sorted(sym)
    stoi = {t: i for i, t in enumerate(itoks)}
    itos = {i: t for t, i in stoi.items()}
    return stoi, itos


def selfies_to_ids(s, stoi, max_len):
    bos, eo, pa, ma = stoi["<BOS>"], stoi["<EOS>"], stoi["<PAD>"], stoi["<MASK>"]
    toks = [bos] + [stoi.get(t, ma) for t in sf.split_selfies(s)] + [eo]
    toks = toks[:max_len]
    x = np.full(max_len, pa, dtype=np.int64)
    x[: len(toks)] = toks
    return x


def selfies_to_prefix_ids(selfies_str, stoi, max_len):
    bos = stoi["<BOS>"]
    pad = stoi["<PAD>"]
    mask = stoi["<MASK>"]
    toks = [bos] + [stoi.get(t, mask) for t in sf.split_selfies(selfies_str)]
    toks = toks[: min(max_len, SELFIES_PREFIX_MAX)]
    if len(toks) < 2:
        return None
    return torch.tensor(toks, dtype=torch.long)


def ids_to_selfies(ids, itos):
    toks = [itos.get(int(i), "") for i in ids if int(i) in itos]
    toks = [t for t in toks if t not in ["<PAD>", "<BOS>", "<EOS>", "<MASK>"]]
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
    def __init__(self, vocab, hidden):
        super().__init__()
        self.emb = nn.Embedding(vocab, hidden, padding_idx=0)
        self.pos = nn.Embedding(SELFIES_MAX_TOK, hidden)
        self.time = TimestepEmbed(hidden)
        self.gru = nn.GRU(hidden, hidden, batch_first=True, bidirectional=True)
        self.fc = nn.Linear(2 * hidden, vocab)

    def forward(self, x, t):
        b, l = x.shape
        h = self.emb(x) + self.pos(
            torch.arange(l, device=x.device)
            .unsqueeze(0)
            .expand(
                b,
                l,
            )
        )
        te = self.time(t).unsqueeze(1).expand_as(h)
        y, _ = self.gru(h + te)
        return self.fc(y)


def make_beta_schedule(T, low=0.02, high=0.3):
    return torch.linspace(low, high, T)


class SelfiesDiffusion:
    def __init__(self, stoi, itos, vocab, hidden, T=200):
        self.stoi, self.itos, self.vocab = stoi, itos, vocab
        self.model = SelfiesDenoiser(vocab, hidden).to(DEVICE)
        self.opt = torch.optim.Adam(self.model.parameters(), lr=LR)
        self.T = T
        self.beta = make_beta_schedule(T).to(DEVICE)
        self.pad = stoi["<PAD>"]
        self.mask = stoi["<MASK>"]
        self.eos = stoi["<EOS>"]
        self.bos = stoi["<BOS>"]

    def corrupt(self, x0, t):
        b, l = x0.shape
        beta_t = self.beta[t - 1].view(-1, 1).clamp(0, 0.99)
        noise = torch.randint(0, self.vocab, (b, l), device=x0.device)
        keep = (torch.rand(b, l, device=x0.device) > beta_t).long()
        xt = keep * x0 + (1 - keep) * noise
        xt = x0 * (x0 == self.pad).long() + xt * (x0 != self.pad).long()
        return xt

    def train_epoch(self, loader):
        self.model.train()
        losses = []
        for x0 in loader:
            x0 = _as_tensor(x0).to(DEVICE)
            t = torch.randint(1, self.T + 1, (x0.size(0),), device=DEVICE)
            xt = self.corrupt(x0, t)
            logits = self.model(xt, t)
            loss = F.cross_entropy(
                logits.reshape(-1, logits.size(-1)),
                x0.reshape(-1),
                ignore_index=self.pad,
            )
            self.opt.zero_grad()
            loss.backward()
            self.opt.step()
            losses.append(loss.item())
        return float(np.mean(losses))

    @torch.no_grad()
    def sample(self, n, max_len, topk=SELFIES_TOPK, topp=SELFIES_TOPP, seed_ids=None):
        self.model.eval()
        x = torch.full((n, max_len), self.pad, device=DEVICE, dtype=torch.long)
        x[:, 0] = self.bos
        if seed_ids is not None:
            k = min(seed_ids.size(0), n)
            for i in range(k):
                pref = seed_ids[i].to(DEVICE)
                pad_pos = (pref == self.pad).nonzero(as_tuple=False)
                L = int(pad_pos[0].item()) if pad_pos.numel() > 0 else int(pref.numel())
                L = max(1, min(L, max_len))
                x[i, :L] = pref[:L]
                x[i, 0] = self.bos
        for t in range(self.T, 0, -1):
            logits = self.model(x, torch.full((n,), t, device=DEVICE))
            probs = F.softmax(logits, dim=-1)
            ban = [self.pad, self.bos, self.mask, self.eos]
            for bid in ban:
                if 0 <= bid < probs.size(-1):
                    probs[..., bid] = 0.0
            probs = probs / (probs.sum(dim=-1, keepdim=True) + 1e-12)
            if topp is not None and topp < 1.0:
                sp, si = torch.sort(probs, dim=-1, descending=True)
                cum = sp.cumsum(-1)
                mask = cum > topp
                mask[:, :, 0] = False
                sp[mask] = 0
                probs = torch.zeros_like(probs).scatter_(-1, si, sp)
                probs = probs / (probs.sum(dim=-1, keepdim=True) + 1e-12)
            if topk is not None:
                k = min(topk, probs.size(-1))
                sp, si = torch.topk(probs, k, dim=-1)
                sp = sp / (sp.sum(-1, keepdim=True) + 1e-12)
                idx = si.gather(
                    -1, torch.multinomial(sp.view(-1, k), 1).view(n, max_len, 1)
                ).squeeze(-1)
            else:
                idx = torch.multinomial(probs.view(-1, probs.size(-1)), 1).view(
                    n, max_len
                )
            fill = x == self.pad
            x = torch.where(fill, idx, x)
            x[:, 0] = self.bos
        out = []
        for row in x.tolist():
            s = ids_to_selfies(row, self.itos)
            sm = from_selfies(s)
            if sm is not None:
                out.append(sm)
        out = canonicalize_and_filter_neutral(out)
        return out


def selfies_pipeline(df, seed_smiles=None, round_id=0):
    base_smiles = df["SMILES"].tolist()
    selfies = [to_selfies(s) for s in base_smiles]
    selfies = [s for s in selfies if s]
    stoi, itos = build_selfies_vocab(selfies)
    vocab = len(stoi)
    X = np.stack([selfies_to_ids(s, stoi, SELFIES_MAX_TOK) for s in selfies], axis=0)
    dataset = torch.utils.data.TensorDataset(torch.tensor(X, dtype=torch.long))
    loader = torch.utils.data.DataLoader(
        dataset, batch_size=BATCH_SIZE, shuffle=True, drop_last=True
    )
    dm = SelfiesDiffusion(stoi, itos, vocab, SELFIES_HIDDEN, T=T_STEPS)
    for ep in range(EPOCHS):
        loss = dm.train_epoch(loader)
        print(f"[SELFIES-Diffusion] epoch {ep+1}/{EPOCHS} loss={loss:.4f}")
    seed_ids = None
    if seed_smiles:
        seed_pref = []
        for sm in seed_smiles:
            ss = to_selfies(sm)
            if not ss:
                continue
            pref = selfies_to_prefix_ids(ss, stoi, SELFIES_MAX_TOK)
            if pref is not None and pref.numel() > 1:
                seed_pref.append(pref)
        if seed_pref:
            seed_ids = torch.nn.utils.rnn.pad_sequence(
                seed_pref, batch_first=True, padding_value=stoi["<PAD>"]
            )
    all_valid = []
    attempts = 0
    while len(all_valid) < MIN_VALID and attempts < 10:
        batch = dm.sample(NUM_GEN, SELFIES_MAX_TOK, seed_ids=seed_ids)
        all_valid.extend(batch)
        all_valid = list(dict.fromkeys(all_valid))
        attempts += 1
    if len(all_valid) < MIN_VALID:
        all_valid = brics_augment(df, need=MIN_VALID, base=all_valid)
    if len(all_valid) < MIN_VALID:
        all_valid = dataset_random_enum(df, need=MIN_VALID, base=all_valid)
    all_valid = canonicalize_and_filter_neutral(all_valid)
    out_csv = OUT_DIR / f"generated_selfies_diffusion_R{round_id}.csv"
    pd.DataFrame({"smiles": all_valid}).to_csv(out_csv, index=False)
    print(f"[SELFIES-Diffusion] valid={len(all_valid)} -> {out_csv}")
    return out_csv


def mol_to_graph(smiles, n_max=GRAPH_N_MAX, atom_vocab=GRAPH_ATOMS):
    m = Chem.MolFromSmiles(smiles)
    if not m:
        return None, None
    m = Chem.AddHs(m)
    n = min(m.GetNumAtoms(), n_max)
    X = np.zeros((n_max, len(atom_vocab)), dtype=np.float32)
    for i, a in enumerate(list(m.GetAtoms())[:n]):
        sym = a.GetSymbol()
        if sym in atom_vocab:
            X[i, atom_vocab.index(sym)] = 1.0
        else:
            X[i, atom_vocab.index("C")] = 1.0
    A = np.zeros((n_max, n_max), dtype=np.float32)
    for b in m.GetBonds():
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()
        if i < n and j < n:
            A[i, j] = 1.0
            A[j, i] = 1.0
    return X, A


def graph_to_mol(X, A, atom_vocab=GRAPH_ATOMS):
    try:
        n = X.shape[0]
        from rdkit.Chem import RWMol

        mol = RWMol()
        idx_map = []
        for i in range(n):
            if X[i].sum() < 0.5:
                continue
            t = int(np.argmax(X[i]))
            idx_map.append(mol.AddAtom(Chem.Atom(atom_vocab[t])))
        if len(idx_map) == 0:
            return None
        for i in range(len(idx_map)):
            for j in range(i + 1, len(idx_map)):
                if A[i, j] > 0.5:
                    try:
                        mol.AddBond(i, j, Chem.BondType.SINGLE)
                    except:
                        pass
        m = mol.GetMol()
        Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        m = Chem.RemoveHs(m)
        sm = Chem.MolToSmiles(m, canonical=True)
        return sm
    except:
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
    def __init__(self, n_max, atom_dim, hidden):
        super().__init__()
        self.time = TimeMLP(hidden)
        in_dim = n_max * atom_dim + n_max * n_max
        self.pre = nn.Linear(in_dim, hidden)
        self.gnn = nn.Sequential(nn.SiLU(), nn.Linear(hidden, hidden), nn.SiLU())
        self.out_x = nn.Linear(hidden, n_max * atom_dim)
        self.out_a = nn.Linear(hidden, n_max * n_max)

    def forward(self, X, A, t):
        b = X.size(0)
        h = torch.cat([X.view(b, -1), A.view(b, -1)], dim=1)
        h = self.pre(h)
        te = self.time(t)
        h = self.gnn(h + te)
        ex = self.out_x(h).view(b, GRAPH_N_MAX, GRAPH_NODE_DIM)
        ea = self.out_a(h).view(b, GRAPH_N_MAX, GRAPH_N_MAX)
        ea = (ea + ea.transpose(1, 2)) / 2.0
        return ex, ea


class GraphDiffusion:
    def __init__(self):
        self.model = GraphDenoiser(GRAPH_N_MAX, GRAPH_NODE_DIM, GRAPH_HIDDEN).to(DEVICE)
        self.opt = torch.optim.Adam(self.model.parameters(), lr=LR)
        self.sig_min = GRAPH_SIGMA_MIN
        self.sig_max = GRAPH_SIGMA_MAX
        self.T = T_STEPS

    def train_epoch(self, loader):
        self.model.train()
        losses = []
        for X0, A0 in loader:
            X0 = _as_tensor(X0).to(DEVICE)
            A0 = _as_tensor(A0).to(DEVICE)
            t = torch.randint(1, self.T + 1, (X0.size(0),), device=DEVICE)
            sigma = (self.sig_min + (self.sig_max - self.sig_min) * (t / self.T)).view(
                -1, 1, 1
            )
            nx = torch.randn_like(X0)
            na = torch.randn_like(A0)
            Xt = X0 + sigma * nx
            At = A0 + sigma * na
            ex, ea = self.model(Xt, At, t)
            loss = ((ex - X0) ** 2).mean() + ((ea - A0) ** 2).mean()
            self.opt.zero_grad()
            loss.backward()
            self.opt.step()
            losses.append(loss.item())
        return float(np.mean(losses))

    @torch.no_grad()
    def sample(self, n, steps=GRAPH_STEPS, seed_graphs=None):
        self.model.eval()
        X = torch.randn(n, GRAPH_N_MAX, GRAPH_NODE_DIM, device=DEVICE)
        A = torch.randn(n, GRAPH_N_MAX, GRAPH_N_MAX, device=DEVICE)
        if seed_graphs is not None:
            k = min(seed_graphs[0].shape[0], n)
            X[:k] = seed_graphs[0][:k].to(DEVICE)
            A[:k] = seed_graphs[1][:k].to(DEVICE)
        for t in range(steps, 0, -1):
            tt = torch.full((n,), t, device=DEVICE)
            ex, ea = self.model(X, A, tt)
            alpha = 0.2
            X = X - alpha * (X - ex)
            A = A - alpha * (A - ea)
            if t > 1:
                X = X + torch.randn_like(X) * 0.01
                A = A + torch.randn_like(A) * 0.01
        X = torch.softmax(X, dim=-1)
        A = torch.sigmoid((A + A.transpose(1, 2)) / 2.0)
        X = (X > 0.5).float()
        A = torch.triu((A > 0.5).float(), diagonal=1)
        A = A + A.transpose(1, 2)
        out = []
        for i in range(n):
            sm = graph_to_mol(X[i].cpu().numpy(), A[i].cpu().numpy())
            if sm:
                out.append(sm)
        return list(dict.fromkeys(out))


class GeoDenoiser(nn.Module):
    def __init__(self, n_max, atom_dim, hidden):
        super().__init__()
        self.time = TimeMLP(hidden)
        in_dim = n_max * 3 + n_max * atom_dim
        self.pre = nn.Linear(in_dim, hidden)
        self.encoder = nn.Sequential(nn.SiLU(), nn.Linear(hidden, hidden), nn.SiLU())
        self.out_x = nn.Linear(hidden, n_max * 3)
        self.out_t = nn.Linear(hidden, n_max * atom_dim)

    def forward(self, coords, types, t):
        b = coords.size(0)
        h = torch.cat([coords.view(b, -1), types.view(b, -1)], dim=1)
        h = self.pre(h)
        te = self.time(t)
        h = self.encoder(h + te)
        dx = self.out_x(h).view(b, GEODIFF_N_MAX, 3)
        dt = self.out_t(h).view(b, GEODIFF_N_MAX, len(GEODIFF_ATOMS))
        return dx, dt


class GeoDiffusion:
    def __init__(self):
        self.model = GeoDenoiser(GEODIFF_N_MAX, len(GEODIFF_ATOMS), GEODIFF_HIDDEN).to(
            DEVICE
        )
        self.opt = torch.optim.Adam(self.model.parameters(), lr=LR)
        self.T = T_STEPS

    def train_epoch(self, loader):
        self.model.train()
        losses = []
        for C0, T0 in loader:
            C0 = _as_tensor(C0).to(DEVICE)
            T0 = _as_tensor(T0).to(DEVICE)
            t = torch.randint(1, self.T + 1, (C0.size(0),), device=DEVICE)
            sigma = (0.01 + 0.99 * (t / self.T)).view(-1, 1, 1)
            nc = torch.randn_like(C0)
            nt = torch.randn_like(T0)
            Ct = C0 + sigma * nc
            Tt = T0 + sigma * nt
            dx, dt = self.model(Ct, Tt, t)
            loss = ((dx - (C0 - Ct)) ** 2).mean() + ((dt - (T0 - Tt)) ** 2).mean()
            self.opt.zero_grad()
            loss.backward()
            self.opt.step()
            losses.append(loss.item())
        return float(np.mean(losses))

    @torch.no_grad()
    def sample(self, n, steps=GEODIFF_STEPS, seed=None):
        self.model.eval()
        C = torch.randn(n, GEODIFF_N_MAX, 3, device=DEVICE) * 0.5
        T = torch.zeros(n, GEODIFF_N_MAX, len(GEODIFF_ATOMS), device=DEVICE)
        T[:, :, GEODIFF_ATOMS.index("C")] = 1.0
        if seed is not None:
            k = min(seed[0].shape[0], n)
            C[:k] = seed[0][:k].to(DEVICE)
            T[:k] = seed[1][:k].to(DEVICE)
        for t in range(steps, 0, -1):
            tt = torch.full((n,), t, device=DEVICE)
            dx, dt = self.model(C, T, tt)
            alpha = 0.2
            C = C + alpha * dx
            T = T + alpha * dt
            if t > 1:
                C = C + torch.randn_like(C) * 0.01
                T = T + torch.randn_like(T) * 0.01
        T = torch.softmax(T, dim=-1)
        return C, T


def coords_to_smiles(C, T):
    try:
        n = C.shape[0]
        from rdkit.Chem import RWMol

        mol = RWMol()
        idx_map = []
        types = np.argmax(T, axis=1)
        for i in range(n):
            if T[i].sum() < 0.5:
                continue
            sym = GEODIFF_ATOMS[int(types[i])]
            idx_map.append(mol.AddAtom(Chem.Atom(sym)))
        if len(idx_map) < 2:
            return None
        D = np.sqrt(((C[None, :] - C[:, None]) ** 2).sum(-1))
        for i in range(len(idx_map)):
            for j in range(i + 1, len(idx_map)):
                d = D[i, j]
                if d < 1.25:
                    bt = Chem.BondType.TRIPLE
                elif d < 1.55:
                    bt = Chem.BondType.DOUBLE
                elif d < 1.9:
                    bt = Chem.BondType.SINGLE
                else:
                    continue
                try:
                    mol.AddBond(i, j, bt)
                except:
                    pass
        m = mol.GetMol()
        try:
            Chem.SanitizeMol(m)
        except:
            Chem.SanitizeMol(m, sanitizeOps=Chem.SanitizeFlags.SANITIZE_NONE)
        try:
            sm = Chem.MolToSmiles(m, canonical=True)
        except:
            return None
        if Chem.MolFromSmiles(sm) is None:
            return None
        return sm
    except:
        return None


def prepare_graph_dataset(df):
    Xs = []
    As = []
    for sm in df["SMILES"].tolist():
        x, a = mol_to_graph(sm, GRAPH_N_MAX, GRAPH_ATOMS)
        if x is None:
            continue
        Xs.append(x)
        As.append(a)
    Xs = np.stack(Xs, axis=0)
    As = np.stack(As, axis=0)
    X = torch.tensor(Xs, dtype=torch.float32)
    A = torch.tensor(As, dtype=torch.float32)
    loader = torch.utils.data.DataLoader(
        torch.utils.data.TensorDataset(X, A),
        batch_size=BATCH_SIZE,
        shuffle=True,
        drop_last=True,
    )
    return loader, (X, A)


def graph_pipeline(df, seed_smiles=None, round_id=0):
    loader, (Xseed, Aseed) = prepare_graph_dataset(df)
    gd = GraphDiffusion()
    for ep in range(EPOCHS):
        loss = gd.train_epoch(loader)
        print(f"[Graph-Diffusion] epoch {ep+1}/{EPOCHS} loss={loss:.6f}")
    seed_graphs = None
    if seed_smiles:
        Xs = []
        As = []
        for s in seed_smiles:
            x, a = mol_to_graph(s, GRAPH_N_MAX, GRAPH_ATOMS)
            if x is not None:
                Xs.append(x)
                As.append(a)
        if Xs:
            Xs = torch.tensor(np.stack(Xs), dtype=torch.float32)
            As = torch.tensor(np.stack(As), dtype=torch.float32)
            seed_graphs = (Xs, As)
    all_valid = []
    attempts = 0
    while len(all_valid) < MIN_VALID and attempts < 20:
        batch = gd.sample(NUM_GEN, seed_graphs=seed_graphs)
        all_valid.extend(batch)
        all_valid = list(dict.fromkeys(all_valid))
        attempts += 1
    if len(all_valid) < MIN_VALID:
        all_valid = brics_augment(df, need=MIN_VALID, base=all_valid)
    if len(all_valid) < MIN_VALID:
        all_valid = dataset_random_enum(df, need=MIN_VALID, base=all_valid)
    all_valid = canonicalize_and_filter_neutral(all_valid)
    out_csv = OUT_DIR / f"generated_graph_diffusion_R{round_id}.csv"
    pd.DataFrame({"smiles": all_valid}).to_csv(out_csv, index=False)
    print(f"[Graph-Diffusion] valid={len(all_valid)} -> {out_csv}")
    return out_csv


def prepare_geodiff_dataset(df):
    Xs = []
    Cs = []
    for sm in df["SMILES"].tolist():
        m = Chem.MolFromSmiles(sm)
        if not m:
            continue
        m = Chem.AddHs(m)
        try:
            AllChem.EmbedMolecule(m, randomSeed=SEED)
            AllChem.UFFOptimizeMolecule(m, maxIters=100)
        except:
            continue
        conf = m.GetConformer()
        n = min(m.GetNumAtoms(), GEODIFF_N_MAX)
        C = np.zeros((GEODIFF_N_MAX, 3), dtype=np.float32)
        T = np.zeros((GEODIFF_N_MAX, len(GEODIFF_ATOMS)), dtype=np.float32)
        for i in range(n):
            p = conf.GetAtomPosition(i)
            C[i] = [p.x, p.y, p.z]
            sym = m.GetAtomWithIdx(i).GetSymbol()
            sym = sym if sym in GEODIFF_ATOMS else "C"
            T[i, GEODIFF_ATOMS.index(sym)] = 1.0
        Cs.append(C)
        Xs.append(T)
    if not Xs:
        raise RuntimeError("Empty GeoDiff dataset")
    C = torch.tensor(np.stack(Cs), dtype=torch.float32)
    T = torch.tensor(np.stack(Xs), dtype=torch.float32)
    loader = torch.utils.data.DataLoader(
        torch.utils.data.TensorDataset(C, T),
        batch_size=BATCH_SIZE,
        shuffle=True,
        drop_last=True,
    )
    return loader, (C, T)


def geodiff_pipeline(df, seed_smiles=None, round_id=0):
    loader, (Cseed, Tseed) = prepare_geodiff_dataset(df)
    gd = GeoDiffusion()
    for ep in range(EPOCHS):
        loss = gd.train_epoch(loader)
        print(f"[GeoDiff] epoch {ep+1}/{EPOCHS} loss={loss:.6f}")
    seed = None
    if seed_smiles:
        Cs = []
        Ts = []
        for s in seed_smiles:
            try:
                m = Chem.MolFromSmiles(s)
                m = Chem.AddHs(m)
                AllChem.EmbedMolecule(m, randomSeed=SEED)
                AllChem.UFFOptimizeMolecule(m, maxIters=50)
                conf = m.GetConformer()
                n = min(m.GetNumAtoms(), GEODIFF_N_MAX)
                C = np.zeros((GEODIFF_N_MAX, 3), dtype=np.float32)
                T = np.zeros((GEODIFF_N_MAX, len(GEODIFF_ATOMS)), dtype=np.float32)
                for i in range(n):
                    p = conf.GetAtomPosition(i)
                    C[i] = [p.x, p.y, p.z]
                    sym = m.GetAtomWithIdx(i).GetSymbol()
                    sym = sym if sym in GEODIFF_ATOMS else "C"
                    T[i, GEODIFF_ATOMS.index(sym)] = 1.0
                Cs.append(C)
                Ts.append(T)
            except:
                pass
        if Cs:
            seed = (
                torch.tensor(np.stack(Cs), dtype=torch.float32),
                torch.tensor(np.stack(Ts), dtype=torch.float32),
            )
    all_valid = []
    tries = 0
    while len(all_valid) < MIN_VALID and tries < 30:
        C, T = gd.sample(NUM_GEN, seed=seed)
        C = C.cpu().numpy()
        T = T.cpu().numpy()
        for i in range(C.shape[0]):
            sm = coords_to_smiles(C[i], T[i])
            if sm:
                all_valid.append(sm)
        all_valid = list(dict.fromkeys(all_valid))
        tries += 1
    if len(all_valid) < MIN_VALID:
        all_valid = brics_augment(df, need=MIN_VALID, base=all_valid)
    if len(all_valid) < MIN_VALID:
        all_valid = dataset_random_enum(df, need=MIN_VALID, base=all_valid)
    all_valid = canonicalize_and_filter_neutral(all_valid)
    out_csv = OUT_DIR / f"generated_geodiff_R{round_id}.csv"
    pd.DataFrame({"smiles": all_valid}).to_csv(out_csv, index=False)
    print(f"[GeoDiff] valid={len(all_valid)} -> {out_csv}")
    return out_csv


def brics_augment(df, need=MIN_VALID, seed=SEED, base=None):
    rng = np.random.RandomState(seed)
    base = base[:] if base else []
    mols = [
        Chem.MolFromSmiles(s) for s in df["SMILES"].tolist() if Chem.MolFromSmiles(s)
    ]
    fr = set()
    for m in mols:
        try:
            fr |= set(BRICS.BRICSDecompose(m))
        except:
            pass
    fr = list(fr)
    tries = 0
    while len(base) < need and tries < 5000 and len(fr) >= 2:
        a, b = rng.choice(fr, 2, replace=True)
        try:
            for m in BRICS.BRICSBuild([a, b]):
                sm = Chem.MolToSmiles(
                    Chem.MolFromSmiles(Chem.MolToSmiles(m, canonical=True)),
                    canonical=True,
                )
                base.append(sm)
        except:
            pass
        base = list(dict.fromkeys(base))
        tries += 1
    base = canonicalize_and_filter_neutral(base)
    return base[:need]


def dataset_random_enum(df, need=MIN_VALID, seed=SEED, base=None):
    rng = np.random.RandomState(seed)
    base = base[:] if base else []
    src = df["SMILES"].tolist()
    rng.shuffle(src)
    for s in src:
        try:
            m = Chem.MolFromSmiles(s)
            if not m:
                continue
            for _ in range(10):
                sm = Chem.MolToSmiles(m, canonical=False, doRandom=True)
                base.append(Chem.MolToSmiles(Chem.MolFromSmiles(sm), canonical=True))
                if len(base) >= need:
                    break
        except:
            pass
        base = list(dict.fromkeys(base))
        if len(base) >= need:
            break
    base = canonicalize_and_filter_neutral(base)
    return base[:need]


def enrich_by_predictions(ckpt, gen_csv, tag, topk=TOPK_SEED):
    preds, hits = chemprop_predict_paths(ckpt, gen_csv, tag)
    if preds is None:
        return None, None
    dfp = pd.read_csv(preds)
    scol = (
        "smiles"
        if "smiles" in dfp.columns
        else ("SMILES" if "SMILES" in dfp.columns else None)
    )
    if scol is None:
        return preds, hits
    pcol = pick_pred_column(dfp)
    if pcol is None:
        return preds, []
    dfp = dfp.sort_values(by=pcol, ascending=False)
    seeds = dfp[scol].head(topk).tolist()
    return preds, seeds


def merge_and_dedup(csv_list, out_name):
    all_sm = []
    for p in csv_list:
        if p and os.path.exists(p):
            try:
                df = pd.read_csv(p)
                col = (
                    "smiles"
                    if "smiles" in df.columns
                    else ("SMILES" if "SMILES" in df.columns else None)
                )
                if col is None:
                    continue
                all_sm += df[col].dropna().astype(str).tolist()
            except:
                pass
    uniq = canonicalize_and_filter_neutral(all_sm)
    out = OUT_DIR / out_name
    pd.DataFrame({"smiles": uniq}).to_csv(out, index=False)
    return out


if __name__ == "__main__":
    df = load_dm_dataframe(DATA_CSV)
    ckpt = chemprop_kfold_train()
    csv_self = selfies_pipeline(df, round_id=0)
    pred_self, seeds_self = enrich_by_predictions(
        ckpt, csv_self, "SELFIES_DIFFUSION_R0"
    )
    csv_graph = graph_pipeline(df, seed_smiles=seeds_self, round_id=0)
    pred_graph, seeds_graph = enrich_by_predictions(
        ckpt, csv_graph, "GRAPH_DIFFUSION_R0"
    )
    seeds_all = (seeds_self or []) + (seeds_graph or [])
    csv_geo = geodiff_pipeline(df, seed_smiles=seeds_all, round_id=0)
    pred_geo, seeds_geo = enrich_by_predictions(ckpt, csv_geo, "GEODIFF_R0")
    merged1 = merge_and_dedup(
        [csv_self, csv_graph, csv_geo], "generated_merged_round0.csv"
    )
    chemprop_predict_paths(ckpt, merged1, "MERGED_R0")
    for r in range(ENRICH_ROUNDS):
        seeds = (seeds_self or []) + (seeds_graph or []) + (seeds_geo or [])
        csv_self = selfies_pipeline(df, seed_smiles=seeds, round_id=r + 1)
        pred_self, seeds_self = enrich_by_predictions(
            ckpt, csv_self, f"SELFIES_DIFFUSION_R{r+1}"
        )
        csv_graph = graph_pipeline(df, seed_smiles=seeds, round_id=r + 1)
        pred_graph, seeds_graph = enrich_by_predictions(
            ckpt, csv_graph, f"GRAPH_DIFFUSION_R{r+1}"
        )
        csv_geo = geodiff_pipeline(df, seed_smiles=seeds, round_id=r + 1)
        pred_geo, seeds_geo = enrich_by_predictions(ckpt, csv_geo, f"GEODIFF_R{r+1}")
        merged = merge_and_dedup(
            [csv_self, csv_graph, csv_geo, merged1], f"generated_merged_round{r+1}.csv"
        )
        chemprop_predict_paths(ckpt, merged, f"MERGED_R{r+1}")
    print("Done.")

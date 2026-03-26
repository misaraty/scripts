import os, sys, re, random, shutil, subprocess, glob

os.chdir(os.path.split(os.path.realpath(__file__))[0])
from pathlib import Path
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
# Basic config
DATA_CSV = "D.csv"
BASE_DIR = Path("runs/mol_chemprop_multi")
CP_DIR = BASE_DIR / "models_chemprop"
OUT_DIR = BASE_DIR / "outputs"
for p in [CP_DIR, OUT_DIR]:
    p.mkdir(parents=True, exist_ok=True)
SEED = 42
torch.manual_seed(SEED)
np.random.seed(SEED)
random.seed(SEED)
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"
EPOCHS_GEN = 30
BATCH_SIZE_GEN = 64
LATENT_DIM = 384
MAX_LEN = 128
TEMPERATURE = 0.9
NUM_GEN = 2048
TARGET = 10.6
TOP_K = 100
TOP_P = 0.9
REPEAT_PENALTY = 1.15
ENRICH_ROUNDS = 5
TOPK_SEED = 256
EPOCHS_VAE = 30
LATENT_DIM_VAE = 128
BETA_VAE = 1.0
KL_ANNEAL = True


# Chemprop
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


def load_dm_dataframe(path):
    df = pd.read_csv(path)
    smiles_col = None
    d_col = None
    for c in df.columns:
        if smiles_col is None and re.search("smiles", c, re.I):
            smiles_col = c
        if d_col is None and re.search(
            r"(^|[^A-Za-z])D([^A-Za-z]|$)|velocity|label|target", c, re.I
        ):
            d_col = c
    if smiles_col is None or d_col is None:
        raise ValueError("Cannot find SMILES/D column.")
    df = (
        df[[smiles_col, d_col]]
        .rename(columns={smiles_col: "SMILES", d_col: "D"})
        .dropna()
        .reset_index(drop=True)
    )
    return df


def simple_kfold_indices(n_splits, n_samples, seed=42):
    rng = np.random.RandomState(seed)
    idx = np.arange(n_samples)
    rng.shuffle(idx)
    folds = []
    sizes = [n_samples // n_splits] * n_splits
    for i in range(n_samples % n_splits):
        sizes[i] += 1
    s = 0
    for size in sizes:
        folds.append(idx[s : s + size])
        s += size
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


def kfold_eval_and_artifacts(df, save_dir, outs_dir, n_folds=5):
    print(
        "[NOTE] kfold_eval_and_artifacts uses external fold split; ensure it matches chemprop internal split for strict OOF evaluation."
    )
    folds = simple_kfold_indices(n_splits=n_folds, n_samples=len(df), seed=42)
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
        val_df[["SMILES"]].rename(columns={"SMILES": "smiles"}).to_csv(
            val_csv, index=False
        )
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
        lower_cols = {c.lower(): c for c in pr.columns}
        preferred = []
        for key in ["d_pred", "pred", "prediction", "y_pred", "d"]:
            if key in lower_cols:
                preferred.append(lower_cols[key])
        if preferred:
            pred_col = preferred[0]
        else:
            num_cols = [c for c in pr.columns if pd.api.types.is_numeric_dtype(pr[c])]
            if not num_cols:
                raise ValueError(f"No numeric prediction column found in {pred_csv}")
            pred_col = num_cols[0]
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
        f.write(f"Overall R2 = {np.nanmean(r2s):.6f} +/- {np.nanstd(r2s,ddof=1):.6f}\n")
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
    if not os.path.exists(gen_csv) or os.path.getsize(gen_csv) == 0:
        return None, None
    try:
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
    hits = df[df["D_pred"] >= float(TARGET)].copy()
    hits.to_csv(hits_csv, index=False)
    print(f"[OK:{tag}] pred={len(df)} hits={len(hits)}")
    return out_csv, hits_csv


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
        return preds, []
    if "D_pred" in dfp.columns:
        sort_col = "D_pred"
    else:
        sort_col = pick_pred_column(dfp)
    if sort_col is None:
        return preds, []
    dfp = dfp.sort_values(by=sort_col, ascending=False)
    seeds = dfp[scol].head(topk).astype(str).tolist()
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
                if col:
                    all_sm += df[col].dropna().astype(str).tolist()
            except:
                pass
    uniq = list(dict.fromkeys(all_sm))
    out = OUT_DIR / out_name
    pd.DataFrame({"smiles": uniq}).to_csv(out, index=False)
    return out


def canonicalize_smiles_list(smiles_list):
    out = []
    bad = 0
    for s in smiles_list:
        try:
            s = str(s).strip()
            mol = Chem.MolFromSmiles(s)
            if mol is None:
                bad += 1
                continue
            out.append(Chem.MolToSmiles(mol, canonical=True))
        except Exception:
            bad += 1
    if bad > 0:
        print(f"[WARN] invalid SMILES filtered: {bad}")
    return out


def pick_pred_column(df: pd.DataFrame):
    lower_cols = {c.lower(): c for c in df.columns}
    for key in ["d_pred", "pred", "prediction", "y_pred", "d"]:
        if key in lower_cols:
            return lower_cols[key]
    num_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    if not num_cols:
        return None
    return num_cols[0]


# Tokenization
PAD_ID = 0


def selfies_tokenize(smiles_list):
    selfies_list = []
    bad = 0
    for s in smiles_list:
        try:
            selfies_list.append(sf.encoder(s))
        except Exception:
            bad += 1
    if bad > 0:
        print(f"[WARN] SELFIES encoding failed: {bad}")
    if not selfies_list:
        raise ValueError(
            "No valid SELFIES after filtering/encoding. Check your SMILES dataset."
        )
    tokens = []
    for s in selfies_list:
        tokens.extend(sf.split_selfies(s))
    vocab = sorted(list(set(tokens)))
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
            toks = sf.split_selfies(s)
            ids = [stoi[t] for t in toks if t in stoi]
        except Exception:
            bad += 1
            continue
        ids_inp = [bos_id] + ids
        ids_tgt = ids + [eos_id]
        if len(ids_inp) < max_len:
            ids_inp = ids_inp + [PAD_ID] * (max_len - len(ids_inp))
        else:
            ids_inp = ids_inp[:max_len]
        if len(ids_tgt) < max_len:
            ids_tgt = ids_tgt + [PAD_ID] * (max_len - len(ids_tgt))
        else:
            ids_tgt = ids_tgt[:max_len]
        seqs.append((ids_inp, ids_tgt))
    if bad > 0:
        print(f"[WARN] selfies_to_tensor skipped {bad} SMILES due to encoder failure.")
    if not seqs:
        raise ValueError(
            "No sequences created in selfies_to_tensor. Check SMILES filtering/tokenization."
        )
    x = torch.tensor([a for a, b in seqs], dtype=torch.long)
    y = torch.tensor([b for a, b in seqs], dtype=torch.long)
    return x, y


# Models
class MolRNN(nn.Module):
    def __init__(self, vocab_size, hidden=LATENT_DIM, extra_tokens=2):
        super().__init__()
        self.embed = nn.Embedding(vocab_size + 1 + extra_tokens, hidden)
        self.gru = nn.GRU(hidden, hidden, batch_first=True)
        self.fc = nn.Linear(hidden, vocab_size + 1 + extra_tokens)

    def forward(self, x, hidden=None):
        emb = self.embed(x)
        out, _ = self.gru(emb, hidden)
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
        self.embed = nn.Embedding(vocab_size + 1 + extra_tokens, d_model)
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
        h = self.embed(x)
        h = self.pos(h)
        if attn_mask is None:
            T = x.size(1)
            attn_mask = torch.triu(torch.ones(T, T, device=x.device), diagonal=1).bool()
        out = self.encoder(h, mask=attn_mask)
        return self.lm_head(out)


class MolVAE(nn.Module):
    def __init__(
        self, vocab_size, hidden=LATENT_DIM, latent_dim=LATENT_DIM_VAE, extra_tokens=2
    ):
        super().__init__()
        self.vocab_size = vocab_size + 1 + extra_tokens
        self.hidden = hidden
        self.latent_dim = latent_dim
        self.embed = nn.Embedding(self.vocab_size, hidden)
        self.enc_rnn = nn.GRU(hidden, hidden, batch_first=True)
        self.fc_mu = nn.Linear(hidden, latent_dim)
        self.fc_logvar = nn.Linear(hidden, latent_dim)
        self.dec_in = nn.Linear(latent_dim, hidden)
        self.dec_rnn = nn.GRU(hidden, hidden, batch_first=True)
        self.fc_out = nn.Linear(hidden, self.vocab_size)

    def encode(self, x):
        emb = self.embed(x)
        _, h = self.enc_rnn(emb)
        h = h[-1]
        mu = self.fc_mu(h)
        logvar = self.fc_logvar(h)
        return mu, logvar

    def reparameterize(self, mu, logvar):
        std = torch.exp(0.5 * logvar)
        eps = torch.randn_like(std)
        return mu + eps * std

    def decode(self, z, x_dec):
        h0 = self.dec_in(z).unsqueeze(0)
        emb = self.embed(x_dec)
        out, _ = self.dec_rnn(emb, h0)
        logits = self.fc_out(out)
        return logits

    def forward(self, x_enc, x_dec):
        mu, logvar = self.encode(x_enc)
        z = self.reparameterize(mu, logvar)
        logits = self.decode(z, x_dec)
        return logits, mu, logvar


# Train
def train_rnn_model(df):
    smiles = canonicalize_smiles_list(df["SMILES"].tolist())
    stoi, itos, vocab, bos_id, eos_id = selfies_tokenize(smiles)
    x, y = selfies_to_tensor(smiles, stoi, MAX_LEN, bos_id, eos_id)
    model = MolRNN(len(vocab)).to(DEVICE)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    for epoch in range(EPOCHS_GEN):
        idx = np.random.permutation(len(x))
        x = x[idx]
        y = y[idx]
        losses = []
        for i in range(0, len(x), BATCH_SIZE_GEN):
            xb = x[i : i + BATCH_SIZE_GEN].to(DEVICE)
            yb = y[i : i + BATCH_SIZE_GEN].to(DEVICE)
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
    torch.save(
        {
            "model": model.state_dict(),
            "stoi": stoi,
            "itos": itos,
            "bos_id": bos_id,
            "eos_id": eos_id,
        },
        OUT_DIR / "molrnn.pt",
    )
    return model, stoi, itos, bos_id, eos_id


def train_transformer_model(df):
    smiles = canonicalize_smiles_list(df["SMILES"].tolist())
    stoi, itos, vocab, bos_id, eos_id = selfies_tokenize(smiles)
    x, y = selfies_to_tensor(smiles, stoi, MAX_LEN, bos_id, eos_id)
    model = MolTransformer(len(vocab)).to(DEVICE)
    opt = torch.optim.AdamW(
        model.parameters(), lr=2e-4, betas=(0.9, 0.98), weight_decay=0.01
    )
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, T_max=EPOCHS_GEN)
    for epoch in range(EPOCHS_GEN):
        idx = np.random.permutation(len(x))
        x = x[idx]
        y = y[idx]
        losses = []
        for i in range(0, len(x), BATCH_SIZE_GEN):
            xb = x[i : i + BATCH_SIZE_GEN].to(DEVICE)
            yb = y[i : i + BATCH_SIZE_GEN].to(DEVICE)
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
    torch.save(
        {
            "model": model.state_dict(),
            "stoi": stoi,
            "itos": itos,
            "bos_id": bos_id,
            "eos_id": eos_id,
        },
        OUT_DIR / "moltransformer.pt",
    )
    return model, stoi, itos, bos_id, eos_id


def kl_divergence(mu, logvar):
    return -0.5 * torch.sum(1 + logvar - mu.pow(2) - logvar.exp(), dim=1).mean()


def train_vae_model(df):
    smiles = canonicalize_smiles_list(df["SMILES"].tolist())
    stoi, itos, vocab, bos_id, eos_id = selfies_tokenize(smiles)
    x_enc, x_tgt = selfies_to_tensor(smiles, stoi, MAX_LEN, bos_id, eos_id)
    model = MolVAE(len(vocab)).to(DEVICE)
    opt = torch.optim.Adam(model.parameters(), lr=1e-3)
    for epoch in range(EPOCHS_VAE):
        idx = np.random.permutation(len(x_enc))
        x_enc = x_enc[idx]
        x_tgt = x_tgt[idx]
        losses = []
        for i in range(0, len(x_enc), BATCH_SIZE_GEN):
            xb_enc = x_enc[i : i + BATCH_SIZE_GEN].to(DEVICE)
            yb = x_tgt[i : i + BATCH_SIZE_GEN].to(DEVICE)
            xb_dec = xb_enc
            logits, mu, logvar = model(xb_enc, xb_dec)
            rec_loss = F.cross_entropy(
                logits.reshape(-1, logits.size(-1)),
                yb.reshape(-1),
                ignore_index=PAD_ID,
                label_smoothing=0.05,
            )
            kl = kl_divergence(mu, logvar)
            if KL_ANNEAL:
                coef = BETA_VAE * min(1.0, float(epoch + 1) / float(EPOCHS_VAE))
            else:
                coef = BETA_VAE
            loss = rec_loss + coef * kl
            opt.zero_grad()
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), 1.0)
            opt.step()
            losses.append(loss.item())
        print(f"[VAE] epoch {epoch+1}/{EPOCHS_VAE} loss={np.mean(losses):.4f}")
    torch.save(
        {
            "model": model.state_dict(),
            "stoi": stoi,
            "itos": itos,
            "bos_id": bos_id,
            "eos_id": eos_id,
        },
        OUT_DIR / "molvae.pt",
    )
    return model, stoi, itos, bos_id, eos_id


# Sampling
def top_k_top_p_filtering(logits, top_k=0, top_p=0.0, min_tokens_to_keep=1):
    B, V = logits.shape
    if top_k and top_k > 0:
        k = int(min(max(top_k, 1), V - min_tokens_to_keep))
        if k > 0:
            kth = torch.topk(logits, k, dim=-1).values[:, -1].unsqueeze(-1)
            logits = torch.where(
                logits < kth, torch.full_like(logits, float("-inf")), logits
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
        logits = torch.full_like(logits, float("-inf"))
        logits.scatter_(1, sorted_idx, sorted_logits)
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
    try:
        sm = sf.decoder(selfies)
        mol = Chem.MolFromSmiles(sm)
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


def smiles_to_prefix_ids(smiles, stoi, bos_id, max_len):
    try:
        s = sf.encoder(smiles)
        toks = [bos_id] + [stoi[t] for t in sf.split_selfies(s) if t in stoi]
        toks = toks[:max_len]
        return torch.tensor(toks, dtype=torch.long)
    except:
        return None


def seed_prefix_pack(
    seed_smiles, stoi, bos_id, max_len, num_gen, device, round_id=0, topk=TOPK_SEED
):
    if not seed_smiles:
        return None
    seed_smiles = sorted(set(seed_smiles))
    prefixes = []
    for s in seed_smiles:
        p = smiles_to_prefix_ids(s, stoi, bos_id, max_len)
        if p is not None and len(p) > 1:
            prefixes.append(p)
    if not prefixes:
        return None
    k = min(len(prefixes), min(num_gen, int(topk)))
    rng = random.Random(SEED + round_id)
    prefixes = rng.sample(prefixes, k)
    pref = torch.nn.utils.rnn.pad_sequence(
        prefixes, batch_first=True, padding_value=PAD_ID
    )
    return pref.to(device)


def sample_rnn(
    model,
    stoi,
    itos,
    bos_id,
    eos_id,
    num_gen=NUM_GEN,
    temp=TEMPERATURE,
    top_k=TOP_K,
    top_p=TOP_P,
    repeat_penalty=REPEAT_PENALTY,
    seed_smiles=None,
    seed_topk=TOPK_SEED,
    round_id=0,
):
    model = model.to(DEVICE).eval()
    seed_pref = seed_prefix_pack(
        seed_smiles,
        stoi,
        bos_id,
        MAX_LEN,
        num_gen,
        DEVICE,
        round_id=round_id,
        topk=seed_topk,
    )
    if seed_pref is not None:
        k = seed_pref.size(0)
        L = seed_pref.size(1)
        rest_n = max(num_gen - k, 0)
        start_rest = torch.full((rest_n, L), PAD_ID, dtype=torch.long, device=DEVICE)
        if rest_n > 0:
            start_rest[:, 0] = bos_id
        start = torch.cat([seed_pref, start_rest], dim=0)
    else:
        start = torch.full((num_gen, 1), bos_id, dtype=torch.long, device=DEVICE)
    finished = torch.zeros(start.size(0), dtype=torch.bool, device=DEVICE)
    with torch.no_grad():
        while start.size(1) < MAX_LEN:
            logits = model(start)[:, -1, :] / max(1e-6, temp)
            if repeat_penalty is not None and repeat_penalty != 1.0:
                penal = torch.zeros_like(logits, dtype=torch.bool)
                tok = start.clamp(min=0)
                penal.scatter_(1, tok, torch.ones_like(tok, dtype=torch.bool))
                penal[:, PAD_ID] = False
                penal[:, bos_id] = False
                logits = logits / (penal.float() * (repeat_penalty - 1.0) + 1.0)
            V = logits.size(-1)
            dyn_k = min(int(top_k), max(1, V - 1))
            logits = top_k_top_p_filtering(logits, top_k=dyn_k, top_p=top_p)
            probs = F.softmax(logits, dim=-1)
            nxt = torch.multinomial(probs, 1)
            start = torch.cat([start, nxt], dim=1)
            finished |= nxt.squeeze(1) == eos_id
            if finished.all():
                break
    smiles = []
    for row in start.detach().cpu().numpy():
        sm = decode_selfies_sequence(row, itos, bos_id, eos_id)
        if sm is not None and is_neutral_smiles(sm):
            smiles.append(sm)
    smiles = list(dict.fromkeys(smiles))
    out_csv = OUT_DIR / f"generated_smiles_RNN_R{round_id}.csv"
    pd.DataFrame({"smiles": smiles}).to_csv(out_csv, index=False)
    return out_csv


def sample_transformer(
    model,
    stoi,
    itos,
    bos_id,
    eos_id,
    num_gen=NUM_GEN,
    temp=TEMPERATURE,
    top_k=TOP_K,
    top_p=TOP_P,
    repeat_penalty=REPEAT_PENALTY,
    seed_smiles=None,
    seed_topk=TOPK_SEED,
    round_id=0,
):
    model = model.to(DEVICE).eval()
    seed_pref = seed_prefix_pack(
        seed_smiles,
        stoi,
        bos_id,
        MAX_LEN,
        num_gen,
        DEVICE,
        round_id=round_id,
        topk=seed_topk,
    )
    if seed_pref is not None:
        k = seed_pref.size(0)
        L = seed_pref.size(1)
        rest_n = max(num_gen - k, 0)
        start_rest = torch.full((rest_n, L), PAD_ID, dtype=torch.long, device=DEVICE)
        if rest_n > 0:
            start_rest[:, 0] = bos_id
        start = torch.cat([seed_pref, start_rest], dim=0)
    else:
        start = torch.full((num_gen, 1), bos_id, dtype=torch.long, device=DEVICE)
    finished = torch.zeros(start.size(0), dtype=torch.bool, device=DEVICE)
    with torch.no_grad():
        while start.size(1) < MAX_LEN:
            Tlen = start.size(1)
            attn_mask = torch.triu(
                torch.ones(Tlen, Tlen, device=DEVICE), diagonal=1
            ).bool()
            logits = model(start, attn_mask=attn_mask)[:, -1, :] / max(1e-6, temp)
            if repeat_penalty is not None and repeat_penalty != 1.0:
                penal = torch.zeros_like(logits, dtype=torch.bool)
                tok = start.clamp(min=0)
                penal.scatter_(1, tok, torch.ones_like(tok, dtype=torch.bool))
                penal[:, PAD_ID] = False
                penal[:, bos_id] = False
                logits = logits / (penal.float() * (repeat_penalty - 1.0) + 1.0)
            V = logits.size(-1)
            dyn_k = min(int(top_k), max(1, V - 1))
            logits = top_k_top_p_filtering(logits, top_k=dyn_k, top_p=top_p)
            probs = F.softmax(logits, dim=-1)
            nxt = torch.multinomial(probs, 1)
            start = torch.cat([start, nxt], dim=1)
            finished |= nxt.squeeze(1) == eos_id
            if finished.all():
                break
    smiles = []
    for row in start.detach().cpu().numpy():
        sm = decode_selfies_sequence(row, itos, bos_id, eos_id)
        if sm is not None and is_neutral_smiles(sm):
            smiles.append(sm)
    smiles = list(dict.fromkeys(smiles))
    out_csv = OUT_DIR / f"generated_smiles_Transformer_R{round_id}.csv"
    pd.DataFrame({"smiles": smiles}).to_csv(out_csv, index=False)
    return out_csv


def sample_vae(
    model,
    stoi,
    itos,
    bos_id,
    eos_id,
    num_gen=NUM_GEN,
    temp=TEMPERATURE,
    top_k=TOP_K,
    top_p=TOP_P,
    repeat_penalty=REPEAT_PENALTY,
    seed_smiles=None,
    seed_topk=TOPK_SEED,
    round_id=0,
):
    model = model.to(DEVICE).eval()
    seed_pref = seed_prefix_pack(
        seed_smiles,
        stoi,
        bos_id,
        MAX_LEN,
        num_gen,
        DEVICE,
        round_id=round_id,
        topk=seed_topk,
    )
    with torch.no_grad():
        if seed_pref is not None:
            k = seed_pref.size(0)
            mu, logvar = model.encode(seed_pref)
            z = model.reparameterize(mu, logvar)
            if k < num_gen:
                extra = num_gen - k
                z_extra = torch.randn(extra, model.latent_dim, device=DEVICE)
                z = torch.cat([z, z_extra], dim=0)
            L = seed_pref.size(1)
            rest_n = max(num_gen - k, 0)
            start_rest = torch.full(
                (rest_n, L), PAD_ID, dtype=torch.long, device=DEVICE
            )
            if rest_n > 0:
                start_rest[:, 0] = bos_id
            start = torch.cat([seed_pref, start_rest], dim=0)
        else:
            z = torch.randn(num_gen, model.latent_dim, device=DEVICE)
            start = torch.full((num_gen, 1), bos_id, dtype=torch.long, device=DEVICE)
        finished = torch.zeros(start.size(0), dtype=torch.bool, device=DEVICE)
        while start.size(1) < MAX_LEN:
            logits = model.decode(z, start)[:, -1, :] / max(1e-6, temp)
            if repeat_penalty is not None and repeat_penalty != 1.0:
                penal = torch.zeros_like(logits, dtype=torch.bool)
                tok = start.clamp(min=0)
                penal.scatter_(1, tok, torch.ones_like(tok, dtype=torch.bool))
                penal[:, PAD_ID] = False
                penal[:, bos_id] = False
                logits = logits / (penal.float() * (repeat_penalty - 1.0) + 1.0)
            V = logits.size(-1)
            dyn_k = min(int(top_k), max(1, V - 1))
            logits = top_k_top_p_filtering(logits, top_k=dyn_k, top_p=top_p)
            probs = F.softmax(logits, dim=-1)
            nxt = torch.multinomial(probs, 1)
            start = torch.cat([start, nxt], dim=1)
            finished |= nxt.squeeze(1) == eos_id
            if finished.all():
                break
    smiles = []
    for row in start.detach().cpu().numpy():
        sm = decode_selfies_sequence(row, itos, bos_id, eos_id)
        if sm is not None and is_neutral_smiles(sm):
            smiles.append(sm)
    smiles = list(dict.fromkeys(smiles))
    out_csv = OUT_DIR / f"generated_smiles_VAE_R{round_id}.csv"
    pd.DataFrame({"smiles": smiles}).to_csv(out_csv, index=False)
    return out_csv


# Main
if __name__ == "__main__":
    df = load_dm_dataframe(DATA_CSV)
    ckpt_root = chemprop_kfold_train()
    rnn, stoi_rnn, itos_rnn, bos_rnn, eos_rnn = train_rnn_model(df)
    tf, stoi_tf, itos_tf, bos_tf, eos_tf = train_transformer_model(df)
    vae, stoi_vae, itos_vae, bos_vae, eos_vae = train_vae_model(df)
    gen_csv_rnn = sample_rnn(
        rnn,
        stoi_rnn,
        itos_rnn,
        bos_rnn,
        eos_rnn,
        num_gen=NUM_GEN,
        temp=TEMPERATURE,
        top_k=TOP_K,
        top_p=TOP_P,
        repeat_penalty=REPEAT_PENALTY,
        seed_smiles=None,
        round_id=0,
    )
    gen_csv_tf = sample_transformer(
        tf,
        stoi_tf,
        itos_tf,
        bos_tf,
        eos_tf,
        num_gen=NUM_GEN,
        temp=TEMPERATURE,
        top_k=TOP_K,
        top_p=TOP_P,
        repeat_penalty=REPEAT_PENALTY,
        seed_smiles=None,
        round_id=0,
    )
    gen_csv_vae = sample_vae(
        vae,
        stoi_vae,
        itos_vae,
        bos_vae,
        eos_vae,
        num_gen=NUM_GEN,
        temp=TEMPERATURE,
        top_k=TOP_K,
        top_p=TOP_P,
        repeat_penalty=REPEAT_PENALTY,
        seed_smiles=None,
        round_id=0,
    )
    pred_rnn, seeds_rnn = enrich_by_predictions(ckpt_root, gen_csv_rnn, "RNN")
    pred_tf, seeds_tf = enrich_by_predictions(ckpt_root, gen_csv_tf, "Transformer")
    pred_vae, seeds_vae = enrich_by_predictions(ckpt_root, gen_csv_vae, "VAE")
    merged0 = merge_and_dedup(
        [gen_csv_rnn, gen_csv_tf, gen_csv_vae], "generated_merged_round0.csv"
    )
    chemprop_predict_paths(ckpt_root, merged0, "MERGED_R0")
    for r in range(ENRICH_ROUNDS):
        seeds = (seeds_rnn or []) + (seeds_tf or []) + (seeds_vae or [])
        gen_csv_rnn = sample_rnn(
            rnn,
            stoi_rnn,
            itos_rnn,
            bos_rnn,
            eos_rnn,
            num_gen=NUM_GEN,
            temp=TEMPERATURE,
            top_k=TOP_K,
            top_p=TOP_P,
            repeat_penalty=REPEAT_PENALTY,
            seed_smiles=seeds,
            round_id=r + 1,
        )
        gen_csv_tf = sample_transformer(
            tf,
            stoi_tf,
            itos_tf,
            bos_tf,
            eos_tf,
            num_gen=NUM_GEN,
            temp=TEMPERATURE,
            top_k=TOP_K,
            top_p=TOP_P,
            repeat_penalty=REPEAT_PENALTY,
            seed_smiles=seeds,
            round_id=r + 1,
        )
        gen_csv_vae = sample_vae(
            vae,
            stoi_vae,
            itos_vae,
            bos_vae,
            eos_vae,
            num_gen=NUM_GEN,
            temp=TEMPERATURE,
            top_k=TOP_K,
            top_p=TOP_P,
            repeat_penalty=REPEAT_PENALTY,
            seed_smiles=seeds,
            round_id=r + 1,
        )
        pred_rnn, seeds_rnn = enrich_by_predictions(
            ckpt_root, gen_csv_rnn, f"RNN_R{r+1}"
        )
        pred_tf, seeds_tf = enrich_by_predictions(
            ckpt_root, gen_csv_tf, f"Transformer_R{r+1}"
        )
        pred_vae, seeds_vae = enrich_by_predictions(
            ckpt_root, gen_csv_vae, f"VAE_R{r+1}"
        )
        merged = merge_and_dedup(
            [gen_csv_rnn, gen_csv_tf, gen_csv_vae, merged0],
            f"generated_merged_round{r+1}.csv",
        )
        chemprop_predict_paths(ckpt_root, merged, f"MERGED_R{r+1}")
    print("Done.")

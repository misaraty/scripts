import os
import datetime
import warnings
import random
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pymatgen.core.structure import Structure
from sklearn.model_selection import train_test_split, KFold
from sklearn.metrics import mean_squared_error, mean_absolute_error, r2_score
from megnet.models import MEGNetModel
from megnet.data.crystal import CrystalGraph
from megnet.data.graph import GaussianDistance
import optuna
import tensorflow as tf

os.environ["TF_CPP_MIN_LOG_LEVEL"] = "3"
warnings.filterwarnings("ignore")
p = plt.rcParams
p["font.sans-serif"] = ["Arial"]
plt.rcParams["font.size"] = 16
SCRIPT_DIR = os.path.split(os.path.realpath(__file__))[0]
os.chdir(SCRIPT_DIR)
STRUCT_DIR = os.path.join(SCRIPT_DIR, "cif")
META_FILE = os.path.join(SCRIPT_DIR, "data.xlsx")
ID_COL = "cif"
TARGET_COL = "tc"
USE_EN_GLOBAL = True
method = "megnet_default"  # 'megnet_default' / 'megnet_tuned'
epochs = 200
n_trials = 20
batch_size = 256
patience = 200
n_folds = 5  # 'none' / 5
train_ratio = 0.7
val_ratio = 0.0
test_ratio = 0.3
lr = 1e-3
random_state = 42
LR_SCHEDULE = "step"  # 'none' / 'step' / 'multistep'
INTERVAL_MIN = 0
INTERVAL_MAX = 200
BINS = np.linspace(INTERVAL_MIN, INTERVAL_MAX, 17)
EN_MAP = {
    "H": "2.20",
    "Li": "0.98",
    "Be": "1.57",
    "B": "2.04",
    "C": "2.55",
    "N": "3.04",
    "O": "3.44",
    "F": "3.98",
    "Ne": "4.48",
    "Na": "0.93",
    "Mg": "1.31",
    "Al": "1.61",
    "Si": "1.90",
    "P": "2.19",
    "S": "2.58",
    "Cl": "3.16",
    "Ar": "3.24",
    "K": "0.82",
    "Ca": "1.00",
    "Sc": "1.36",
    "Ti": "1.54",
    "V": "1.63",
    "Cr": "1.66",
    "Mn": "1.55",
    "Fe": "1.83",
    "Co": "1.88",
    "Ni": "1.91",
    "Cu": "1.90",
    "Zn": "1.65",
    "Ga": "1.81",
    "Ge": "2.01",
    "As": "2.18",
    "Se": "2.55",
    "Br": "2.96",
    "Kr": "3.00",
    "Rb": "0.82",
    "Sr": "0.95",
    "Y": "1.22",
    "Zr": "1.33",
    "Nb": "1.60",
    "Mo": "2.16",
    "Tc": "1.90",
    "Ru": "2.20",
    "Rh": "2.28",
    "Pd": "2.20",
    "Ag": "1.93",
    "Cd": "1.69",
    "In": "1.78",
    "Sn": "1.96",
    "Sb": "2.05",
    "Te": "2.10",
    "I": "2.66",
    "Xe": "2.60",
    "Cs": "0.79",
    "Ba": "0.89",
    "La": "1.10",
    "Ce": "1.12",
    "Pr": "1.13",
    "Nd": "1.14",
    "Pm": "1.13",
    "Sm": "1.17",
    "Eu": "1.20",
    "Gd": "1.20",
    "Tb": "1.22",
    "Dy": "1.23",
    "Ho": "1.24",
    "Er": "1.24",
    "Tm": "1.25",
    "Yb": "1.10",
    "Lu": "1.27",
    "Hf": "1.30",
    "Ta": "1.50",
    "W": "2.36",
    "Re": "1.90",
    "Os": "2.20",
    "Ir": "2.20",
    "Pt": "2.28",
    "Au": "2.54",
    "Hg": "2.00",
    "Tl": "1.62",
    "Pb": "2.33",
    "Bi": "2.02",
    "Po": "2.00",
    "At": "2.20",
    "Rn": "2.60",
    "Fr": "0.70",
    "Ra": "0.89",
    "Ac": "1.10",
    "Th": "1.30",
    "Pa": "1.50",
    "U": "1.38",
    "Np": "1.36",
    "Pu": "1.28",
    "Am": "1.30",
    "Cm": "1.30",
    "Bk": "1.30",
    "Cf": "1.30",
    "Es": "1.30",
    "Fm": "1.30",
    "Md": "1.30",
    "No": "1.30",
}
gpus = tf.config.list_physical_devices("GPU")
if gpus:
    try:
        for _gpu in gpus:
            tf.config.experimental.set_memory_growth(_gpu, True)
        print(f"Using GPU: {len(gpus)} device(s) available.")
    except Exception as e:
        print(f"Warning: failed to set GPU memory growth: {e}")
else:
    print("No GPU detected. Running on CPU.")
random.seed(random_state)
np.random.seed(random_state)
tf.random.set_seed(random_state)
try:
    tf.config.experimental.enable_op_determinism(True)
except Exception:
    pass


def read_structure_by_id(root_dir, _id):
    cand = [
        os.path.join(root_dir, f"{_id}.cif"),
        os.path.join(root_dir, f"{_id}.vasp"),
        os.path.join(root_dir, f"{_id}.POSCAR"),
        os.path.join(root_dir, f"{_id}.CONTCAR"),
        os.path.join(root_dir, f"{_id}.poscar"),
        os.path.join(root_dir, f"{_id}.contcar"),
    ]
    for pth in cand:
        if os.path.isfile(pth):
            return Structure.from_file(pth)
    pats = [
        os.path.join(root_dir, f"{_id}.*"),
        os.path.join(root_dir, f"{_id}*"),
    ]
    for pat in pats:
        for fp in sorted(glob.glob(pat)):
            try:
                return Structure.from_file(fp)
            except Exception:
                continue
    raise FileNotFoundError(f"No structure file found for id={_id}")


def en_state_for_structure(struct: Structure):
    vals = []
    for site in struct.sites:
        try:
            sym = str(site.specie.symbol)
        except Exception:
            sym = str(site.species_string).split()[0]
        vals.append(float(EN_MAP.get(sym, "0.0")))
    if not vals:
        vals = [0.0]
    arr = np.array(vals, dtype=float)
    mean = float(np.mean(arr))
    std = float(np.std(arr))
    return [[mean, std]]


def inject_state(structs, use_en):
    if use_en:
        for s in structs:
            s.state = en_state_for_structure(s)
    else:
        for s in structs:
            if not hasattr(s, "state") or s.state is None:
                s.state = [[0.0, 0.0]]


def build_megnet(n1=64, n2=32, n3=16, lr=1e-3):
    return MEGNetModel(
        nfeat_edge=100,
        nfeat_global=2,
        npass=3,
        nblocks=3,
        n1=n1,
        n2=n2,
        n3=n3,
        graph_converter=CrystalGraph(
            bond_converter=GaussianDistance(
                centers=np.linspace(0, 6.0, 100), width=0.5
            ),
            cutoff=5.0,
        ),
        lr=lr,
        embedding_dim=16,
        ntarget=1,
    )


def predict_batch(model: MEGNetModel, structs):
    return np.array(model.predict_structures(structs)).ravel()


def evaluate_and_dump(split_name, y_true, y_pred, fp, logf):
    y_true = np.asarray(y_true).ravel()
    y_pred = np.asarray(y_pred).ravel()
    rmse = mean_squared_error(y_true, y_pred, squared=False)
    mae = mean_absolute_error(y_true, y_pred)
    r2 = r2_score(y_true, y_pred)
    print(f"{split_name} RMSE: {rmse:.5f}")
    print(f"{split_name} MAE: {mae:.5f}")
    print(f"{split_name} R^2: {r2:.5f}")
    print(f"{split_name} RMSE: {rmse:.5f}", file=logf)
    print(f"{split_name} MAE: {mae:.5f}", file=logf)
    print(f"{split_name} R^2: {r2:.5f}", file=logf)
    np.savetxt(
        fp,
        np.column_stack([y_true, y_pred]),
        header=f"------{split_name}_true------|----{split_name}_predicted-----",
    )
    return rmse, mae, r2


def plot_scatter_with_marginals(split: str, dat_path: str, out_png: str):
    data = np.loadtxt(dat_path, skiprows=1)
    X = data[:, 0]
    Y = data[:, 1]
    color_map = {"train": "tab:blue", "vali": "tab:green", "test": "tab:purple"}
    sca_color = color_map.get(split, "tab:blue")
    fig = plt.figure(constrained_layout=True, figsize=(5, 5))
    gspec = gridspec.GridSpec(
        ncols=2, nrows=2, figure=fig, height_ratios=[1, 7.5], width_ratios=[7.5, 1]
    )
    ax = plt.subplot(gspec[1, 0])
    ax.scatter(X, Y, s=50, c=sca_color, edgecolor="white", alpha=1)
    ax.plot(
        [INTERVAL_MIN, INTERVAL_MAX],
        [INTERVAL_MIN, INTERVAL_MAX],
        c="tab:orange",
        linestyle="dashed",
        lw=1,
        alpha=1,
    )
    ax.set_xlim(INTERVAL_MIN, INTERVAL_MAX)
    ax.set_xticks(np.linspace(INTERVAL_MIN, INTERVAL_MAX, 5))
    ax.set_ylim(INTERVAL_MIN, INTERVAL_MAX)
    ax.set_yticks(np.linspace(INTERVAL_MIN, INTERVAL_MAX, 5))
    ax.set_xlabel("DFT (eV)")
    ax.set_ylabel("Predicted (eV)")
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax = plt.subplot(gspec[0, 0])
    ax.set_xlim(INTERVAL_MIN, INTERVAL_MAX)
    ax.set_xticks(np.linspace(INTERVAL_MIN, INTERVAL_MAX, 5))
    ax.set_xticklabels([])
    ax.set_yticks([])
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.hist(X, bins=BINS, facecolor="tab:orange", edgecolor="white", linewidth=2)
    ax = plt.subplot(gspec[1, 1])
    ax.set_ylim(INTERVAL_MIN, INTERVAL_MAX)
    ax.set_yticks(np.linspace(INTERVAL_MIN, INTERVAL_MAX, 5))
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.spines["right"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.hist(
        Y,
        bins=BINS,
        facecolor=sca_color,
        edgecolor="white",
        orientation="horizontal",
        linewidth=2,
    )
    plt.savefig(out_png, dpi=600)
    plt.close()


def _read_loss_dat(loss_path):
    if not os.path.isfile(loss_path):
        return None
    try:
        dat = np.loadtxt(loss_path, skiprows=1)
        if dat.ndim == 1:
            dat = dat[None, :]
        if dat.shape[1] == 3:
            ep, tr, va = dat[:, 0], dat[:, 1], dat[:, 2]
        elif dat.shape[1] == 2:
            ep, tr = dat[:, 0], dat[:, 1]
            va = np.full_like(tr, np.nan)
        else:
            return None
        return ep, tr, va
    except Exception:
        return None


def plot_cv_mean_std(loss_files, out_png="cv_loss_mean_std.jpg", ylabel="RMSE (eV)"):
    eps, trs, vas = [], [], []
    for p in loss_files:
        r = _read_loss_dat(p)
        if r is None:
            continue
        ep, tr, va = r
        eps.append(ep)
        trs.append(tr)
        vas.append(va)
    if len(trs) == 0:
        print("Warning: no valid loss.dat found, skip CV mean±std plot.")
        return
    min_len = min(len(x) for x in trs)
    ep = eps[0][:min_len]
    tr_arr = np.stack([x[:min_len] for x in trs], axis=0)
    va_arr = np.stack([x[:min_len] for x in vas], axis=0)
    tr_mean, tr_std = np.nanmean(tr_arr, axis=0), np.nanstd(tr_arr, axis=0)
    va_mean, va_std = np.nanmean(va_arr, axis=0), np.nanstd(va_arr, axis=0)
    out_dat = os.path.splitext(out_png)[0] + ".dat"
    with open(out_dat, "w+", encoding="utf-8") as f:
        f.write("epoch train_mean train_std val_mean val_std\n")
        for i in range(len(ep)):
            f.write(
                f"{int(ep[i])} {tr_mean[i]:.6f} {tr_std[i]:.6f} {va_mean[i]:.6f} {va_std[i]:.6f}\n"
            )
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(ep, tr_mean, label="Training subset")
    ax.fill_between(ep, tr_mean - tr_std, tr_mean + tr_std, alpha=0.25)
    ax.plot(ep, va_mean, label="Validation subset")
    ax.fill_between(ep, va_mean - va_std, va_mean + va_std, alpha=0.25)
    ax.set_xlabel("Epochs")
    ax.set_ylabel(ylabel)
    ax.legend()
    ax.tick_params(which="both", top=False, right=False)
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


class RMSECallback(tf.keras.callbacks.Callback):
    def __init__(self, model_ref, X_tr, y_tr, X_va, y_va, save_prefix):
        super().__init__()
        self.model_ref = model_ref
        self.X_tr = X_tr
        self.y_tr = y_tr
        self.X_va = X_va
        self.y_va = y_va
        self.logs = []
        self.save_prefix = save_prefix

    def on_epoch_end(self, epoch, logs=None):
        ypt = predict_batch(self.model_ref, self.X_tr)
        ypv = predict_batch(self.model_ref, self.X_va)
        tr = mean_squared_error(self.y_tr, ypt, squared=False)
        va = mean_squared_error(self.y_va, ypv, squared=False)
        self.logs.append((epoch + 1, tr, va))

    def save(self):
        if not self.logs:
            return
        with open(f"{self.save_prefix}_loss.dat", "w+", encoding="utf-8") as f:
            print("epoch train_rmse val_rmse", file=f)
            for e, tr, va in self.logs:
                print(f"{e} {tr:.6f} {va:.6f}", file=f)
        dat = np.array([[e, tr, va] for e, tr, va in self.logs])
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(dat[:, 0], dat[:, 1], label="Training subset")
        ax.plot(dat[:, 0], dat[:, 2], label="Validation subset")
        ax.set_xlabel("Epochs")
        ax.set_ylabel("RMSE (eV)")
        ax.legend()
        ax.tick_params(which="both", top=False, right=False)
        plt.tight_layout()
        plt.savefig(f"{self.save_prefix}_loss_curve.jpg", dpi=300)
        plt.close()


def train_with_monitor(
    model,
    Xtr,
    ytr,
    Xva,
    yva,
    epochs,
    batch_size,
    patience,
    save_prefix="",
    early_stop=True,
):
    cb = []
    rmse_cb = None
    has_val = isinstance(Xva, list) and len(Xva) > 0
    if has_val and save_prefix:
        rmse_cb = RMSECallback(model, Xtr, ytr, Xva, yva, save_prefix)
        cb.append(rmse_cb)
    if has_val and early_stop:
        cb.append(
            tf.keras.callbacks.EarlyStopping(
                monitor="val_loss", patience=patience, restore_best_weights=True
            )
        )
    if LR_SCHEDULE == "step":
        cb.append(make_step_lr_callback(step_size=50, gamma=0.5))
    elif LR_SCHEDULE == "multistep":
        try:
            base_lr = float(tf.keras.backend.get_value(model.optimizer.learning_rate))
        except Exception:
            base_lr = None
        cb.append(
            make_multistep_lr_callback(
                milestones=[100, 150], gamma=0.3, base_lr=base_lr
            )
        )
    if has_val:
        history = model.train(
            train_structures=Xtr,
            train_targets=ytr,
            validation_structures=Xva,
            validation_targets=yva,
            epochs=epochs,
            verbose=1,
            batch_size=batch_size,
            callbacks=cb,
            save_checkpoint=False,
        )
    else:
        history = model.train(
            train_structures=Xtr,
            train_targets=ytr,
            epochs=epochs,
            verbose=1,
            batch_size=batch_size,
            callbacks=cb,
            save_checkpoint=False,
        )
    if rmse_cb:
        rmse_cb.save()
    return history


def run_and_dump_fold(model, out_dir, prefix, Xtr, ytr, Xva, yva, Xte, yte):
    os.makedirs(out_dir, exist_ok=True)
    y_pred_train = predict_batch(model, Xtr)
    np.savetxt(
        os.path.join(out_dir, "scatter_train.dat"),
        np.column_stack([ytr, y_pred_train]),
        header="------train_true------|----train_predicted-----",
    )
    if isinstance(Xva, list) and len(Xva) > 0:
        y_pred_val = predict_batch(model, Xva)
        np.savetxt(
            os.path.join(out_dir, "scatter_vali.dat"),
            np.column_stack([yva, y_pred_val]),
            header="------vali_true------|----vali_predicted-----",
        )
        plot_scatter_with_marginals(
            "vali",
            os.path.join(out_dir, "scatter_vali.dat"),
            os.path.join(out_dir, "scatter_plot_vali.jpg"),
        )
    else:
        y_pred_val = np.array([])
    y_pred_test = predict_batch(model, Xte)
    np.savetxt(
        os.path.join(out_dir, "scatter_test.dat"),
        np.column_stack([yte, y_pred_test]),
        header="------test_true------|----test_predicted-----",
    )
    plot_scatter_with_marginals(
        "train",
        os.path.join(out_dir, "scatter_train.dat"),
        os.path.join(out_dir, "scatter_plot_train.jpg"),
    )
    plot_scatter_with_marginals(
        "test",
        os.path.join(out_dir, "scatter_test.dat"),
        os.path.join(out_dir, "scatter_plot_test.jpg"),
    )
    with open(os.path.join(out_dir, f"{prefix}.log"), "w+", encoding="utf-8") as logf:
        tr = evaluate_and_dump(
            "train", ytr, y_pred_train, os.path.join(out_dir, "scatter_train.dat"), logf
        )
        if isinstance(Xva, list) and len(Xva) > 0:
            va = evaluate_and_dump(
                "vali", yva, y_pred_val, os.path.join(out_dir, "scatter_vali.dat"), logf
            )
        else:
            va = (np.nan, np.nan, np.nan)
        te = evaluate_and_dump(
            "test", yte, y_pred_test, os.path.join(out_dir, "scatter_test.dat"), logf
        )
    return tr, va, te


def _first_exists(paths):
    for p in paths:
        if os.path.exists(p):
            return p
    return None


def _save_split_csvs_and_log(
    model, Xtr, ytr, Xva, yva, Xte, yte, out_dir, out_prefix, log_file
):
    tr, va, te = run_and_dump_fold(
        model, out_dir, out_prefix, Xtr, ytr, Xva, yva, Xte, yte
    )
    np.savetxt(
        os.path.join(out_dir, "train_results.csv"),
        np.column_stack([ytr, predict_batch(model, Xtr)]),
        delimiter=",",
        header="true,pred",
        comments="",
    )
    np.savetxt(
        os.path.join(out_dir, "vali_results.csv"),
        np.column_stack([yva, predict_batch(model, Xva)]),
        delimiter=",",
        header="true,pred",
        comments="",
    )
    np.savetxt(
        os.path.join(out_dir, "test_results.csv"),
        np.column_stack([yte, predict_batch(model, Xte)]),
        delimiter=",",
        header="true,pred",
        comments="",
    )
    with open(log_file, "a+", encoding="utf-8") as logf:
        print(f"Final evaluation on best fold:", file=logf)
        print(f"  Train RMSE/MAE/R2: {tr}", file=logf)
        print(f"  Val   RMSE/MAE/R2: {va}", file=logf)
        print(f"  Test  RMSE/MAE/R2: {te}", file=logf)


def seconds_str(delta):
    return f"{int(delta.total_seconds())}s"


def make_step_lr_callback(step_size=50, gamma=0.5):
    def schedule(epoch, lr):
        if (epoch + 1) % step_size == 0:
            return lr * gamma
        return lr

    return tf.keras.callbacks.LearningRateScheduler(schedule, verbose=0)


def make_multistep_lr_callback(milestones=(100, 150), gamma=0.3, base_lr=None):
    milestones = sorted(list(milestones))

    def schedule(epoch, lr):
        k = sum(1 for m in milestones if (epoch + 1) >= m)
        if base_lr is None:
            return lr * (gamma if ((epoch + 1) in milestones) else 1.0)
        else:
            return base_lr * (gamma**k)

    return tf.keras.callbacks.LearningRateScheduler(schedule, verbose=0)


# ---- Data loading ----
ext = os.path.splitext(META_FILE)[1].lower()
df = pd.read_excel(META_FILE) if ext in (".xlsx", ".xls") else pd.read_csv(META_FILE)
ids = df[ID_COL].astype(str).tolist()
y_all = df[TARGET_COL].to_numpy().astype(float).ravel()
structures = [read_structure_by_id(STRUCT_DIR, _id) for _id in ids]
if USE_EN_GLOBAL:
    inject_state(structures, use_en=True)
assert (
    abs(train_ratio + val_ratio + test_ratio - 1.0) < 1e-8
), f"train/val/test ratio must sum to 1, got {train_ratio}+{val_ratio}+{test_ratio}"
X_pool, X_test, y_pool, y_test = train_test_split(
    structures, y_all, test_size=test_ratio, shuffle=True, random_state=random_state
)
if n_folds == "none":
    if val_ratio > 0:
        val_frac = val_ratio / (1.0 - test_ratio)
        X_train, X_val, y_train, y_val = train_test_split(
            X_pool, y_pool, test_size=val_frac, shuffle=True, random_state=random_state
        )
    else:
        X_train, X_val, y_train, y_val = X_pool, [], y_pool, []
start_time_all = datetime.datetime.now()
# ---- Main workflows ----
if method == "megnet_tuned":
    if n_folds == "none":
        X_tr, X_te, y_tr, y_te = X_pool, X_test, y_pool, y_test
        if val_ratio > 0:
            X_va, y_va = X_val, y_val
            X_tr, y_tr = X_train, y_train
        else:
            X_va, y_va = [], []

        def objective(trial):
            trial_dir = f"optuna_trial_{trial.number}"
            os.makedirs(trial_dir, exist_ok=True)
            params = {
                "n1": trial.suggest_int("n1", 48, 95, step=8),
                "n2": trial.suggest_int("n2", 24, 47, step=4),
                "n3": trial.suggest_int("n3", 12, 23, step=2),
                "lr": 1e-3,
            }
            model = build_megnet(**params)
            with open(
                os.path.join(trial_dir, f"optuna_trial_{trial.number}.log"),
                "w+",
                encoding="utf-8",
            ) as flog:
                model.model.summary(print_fn=lambda x: print(x, file=flog))
            start = datetime.datetime.now()
            _ = train_with_monitor(
                model,
                X_tr,
                y_tr,
                X_va,
                y_va,
                epochs,
                batch_size,
                patience,
                save_prefix=os.path.join(trial_dir, f"optuna_trial_{trial.number}"),
                early_stop=True,
            )
            end = datetime.datetime.now()
            try:
                model.save_model(
                    os.path.join(trial_dir, f"model_trial_{trial.number}.hdf5")
                )
            except Exception:
                try:
                    model.model.save(
                        os.path.join(trial_dir, f"model_trial_{trial.number}.h5")
                    )
                except Exception:
                    pass
            tr, va, te = run_and_dump_fold(
                model,
                trial_dir,
                f"optuna_trial_{trial.number}",
                X_tr,
                y_tr,
                X_va,
                y_va,
                X_te,
                y_te,
            )
            if isinstance(X_va, list) and len(X_va) > 0:
                y_pred_val = predict_batch(model, X_va)
                val_rmse = mean_squared_error(y_va, y_pred_val, squared=False)
            else:
                y_pred_val = predict_batch(model, X_tr)
                val_rmse = mean_squared_error(y_tr, y_pred_val, squared=False)
            with open(
                os.path.join(trial_dir, f"optuna_trial_{trial.number}.log"),
                "a+",
                encoding="utf-8",
            ) as flog:
                print(f"Best Value (val RMSE): {val_rmse:.12g}", file=flog)
                print(f"Best Params: {params}", file=flog)
                print(f"Runtime: {end - start}", file=flog)
                print(f"RuntimeSeconds: {seconds_str(end - start)}", file=flog)
            return val_rmse

        study = optuna.create_study(
            direction="minimize", sampler=optuna.samplers.TPESampler(seed=random_state)
        )
        study.optimize(objective, n_trials=n_trials)
        with open("megnet_tuned.log", "w+", encoding="utf-8") as logf:
            print(f"Best Number: {study.best_trial.number}", file=logf)
            print(f"Best Value (val RMSE): {study.best_trial.value}", file=logf)
            print(f"Best Params: {study.best_params}", file=logf)
            print(
                f"Best Trial Folder: optuna_trial_{study.best_trial.number}", file=logf
            )
            study.best_trial.user_attrs
            model = build_megnet(**study.best_params)
            model.model.summary(print_fn=lambda x: print(x, file=logf))
            end_time_all = datetime.datetime.now()
            print(f"Runtime: {end_time_all - start_time_all}", file=logf)
            print(
                f"RuntimeSeconds: {seconds_str(end_time_all - start_time_all)}",
                file=logf,
            )
    else:
        k = int(n_folds)
        kf = KFold(n_splits=k, shuffle=True, random_state=random_state)

        def objective(trial):
            trial_dir = f"optuna_trial_{trial.number}"
            os.makedirs(trial_dir, exist_ok=True)
            params = {
                "n1": trial.suggest_int("n1", 16, 64),
                "n2": trial.suggest_int("n2", 8, 32),
                "n3": trial.suggest_int("n3", 4, 16),
                "lr": trial.suggest_float("lr", 1e-4, 5e-3, log=True),
            }
            model_tmp = build_megnet(**params)
            with open(
                os.path.join(trial_dir, f"optuna_trial_{trial.number}.log"),
                "w+",
                encoding="utf-8",
            ) as flog:
                model_tmp.model.summary(print_fn=lambda x: print(x, file=flog))
            val_rmses = []
            start_trial = datetime.datetime.now()
            for fi, (tr_idx, va_idx) in enumerate(kf.split(X_pool)):
                Xtr = [X_pool[i] for i in tr_idx]
                ytr = y_pool[tr_idx].ravel()
                Xva = [X_pool[i] for i in va_idx]
                yva = y_pool[va_idx].ravel()
                fold_dir = os.path.join(trial_dir, f"fold_{fi}")
                os.makedirs(fold_dir, exist_ok=True)
                model = build_megnet(**params)
                start = datetime.datetime.now()
                _ = train_with_monitor(
                    model,
                    Xtr,
                    ytr,
                    Xva,
                    yva,
                    epochs,
                    batch_size,
                    patience,
                    save_prefix=os.path.join(
                        fold_dir, f"optuna_trial_{trial.number}_fold_{fi}"
                    ),
                    early_stop=True,
                )
                end = datetime.datetime.now()
                try:
                    model.save_model(
                        os.path.join(
                            fold_dir, f"model_trial_{trial.number}_fold_{fi}.hdf5"
                        )
                    )
                except Exception:
                    try:
                        model.model.save(
                            os.path.join(
                                fold_dir, f"model_trial_{trial.number}_fold_{fi}.h5"
                            )
                        )
                    except Exception:
                        pass
                run_and_dump_fold(
                    model,
                    fold_dir,
                    f"optuna_trial_{trial.number}_fold_{fi}",
                    Xtr,
                    ytr,
                    Xva,
                    yva,
                    X_test,
                    y_test,
                )
                y_pred_val = predict_batch(model, Xva)
                val_rmse = mean_squared_error(yva, y_pred_val, squared=False)
                val_rmses.append(val_rmse)
                with open(
                    os.path.join(fold_dir, "val_metric.txt"), "w+", encoding="utf-8"
                ) as f:
                    f.write(f"val_rmse={val_rmse:.6f}\n")
                    f.write(f"Runtime: {end - start}\n")
                    f.write(f"RuntimeSeconds: {seconds_str(end - start)}\n")
            avg_rmse = float(np.mean(val_rmses))
            best_idx = int(np.nanargmin(val_rmses))
            best_val_rmse = float(val_rmses[best_idx])
            with open(
                os.path.join(trial_dir, "cv_summary.txt"), "w+", encoding="utf-8"
            ) as f:
                f.write("fold,val_rmse\n")
                for i, v in enumerate(val_rmses):
                    f.write(f"{i},{v:.6f}\n")
                f.write(f"avg_val_rmse,{avg_rmse:.6f}\n")
                f.write(f"best_fold,fold_{best_idx}\n")
                f.write(f"best_val_rmse,{best_val_rmse:.6f}\n")
            loss_files = [
                os.path.join(
                    trial_dir,
                    f"fold_{i}",
                    f"optuna_trial_{trial.number}_fold_{i}_loss.dat",
                )
                for i in range(k)
            ]
            plot_cv_mean_std(
                loss_files,
                out_png=os.path.join(trial_dir, "cv_loss_mean_std.jpg"),
                ylabel="RMSE (eV)",
            )
            tr_idx, va_idx = list(kf.split(X_pool))[best_idx]
            Xtr = [X_pool[i] for i in tr_idx]
            ytr = y_pool[tr_idx].ravel()
            Xva = [X_pool[i] for i in va_idx]
            yva = y_pool[va_idx].ravel()
            model_path = _first_exists(
                [
                    os.path.join(
                        trial_dir,
                        f"fold_{best_idx}",
                        f"optuna_trial_{trial.number}_fold_{best_idx}.hdf5",
                    ),
                    os.path.join(
                        trial_dir,
                        f"fold_{best_idx}",
                        f"optuna_trial_{trial.number}_fold_{best_idx}.h5",
                    ),
                    os.path.join(
                        trial_dir,
                        f"fold_{best_idx}",
                        f"model_trial_{trial.number}_fold_{best_idx}.hdf5",
                    ),
                    os.path.join(
                        trial_dir,
                        f"fold_{best_idx}",
                        f"model_trial_{trial.number}_fold_{best_idx}.h5",
                    ),
                ]
            )
            best_model = build_megnet(**params)
            if model_path:
                best_model = MEGNetModel.from_file(model_path)
            _save_split_csvs_and_log(
                best_model,
                Xtr,
                ytr,
                Xva,
                yva,
                X_test,
                y_test,
                SCRIPT_DIR,
                f"best_fold_trial_{trial.number}",
                os.path.join(trial_dir, f"optuna_trial_{trial.number}.log"),
            )
            end_trial = datetime.datetime.now()
            with open(
                os.path.join(trial_dir, f"optuna_trial_{trial.number}.log"),
                "a+",
                encoding="utf-8",
            ) as flog:
                print(f"Best Value (avg val RMSE): {avg_rmse:.12g}", file=flog)
                print(f"Best Params: {params}", file=flog)
                print(f"Runtime: {end_trial - start_trial}", file=flog)
                print(
                    f"RuntimeSeconds: {seconds_str(end_trial - start_trial)}", file=flog
                )
            return avg_rmse

        study = optuna.create_study(
            direction="minimize", sampler=optuna.samplers.TPESampler(seed=random_state)
        )
        study.optimize(objective, n_trials=n_trials)
        with open("megnet_tuned.log", "w+", encoding="utf-8") as logf:
            print(f"Best Number: {study.best_trial.number}", file=logf)
            print(f"Best Value (avg val RMSE): {study.best_trial.value}", file=logf)
            print(f"Best Params: {study.best_params}", file=logf)
            print(
                f"Best Trial Folder: optuna_trial_{study.best_trial.number}", file=logf
            )
            best_params = study.best_params
            best_model = build_megnet(**best_params)
            best_model.model.summary(print_fn=lambda x: print(x, file=logf))
            end_time_all = datetime.datetime.now()
            print(f"Runtime: {end_time_all - start_time_all}", file=logf)
            print(
                f"RuntimeSeconds: {seconds_str(end_time_all - start_time_all)}",
                file=logf,
            )
elif method == "megnet_default":
    if n_folds == "none":
        X_tr, X_te, y_tr, y_te = X_pool, X_test, y_pool, y_test
        if val_ratio > 0:
            X_va, y_va = X_val, y_val
            X_tr, y_tr = X_train, y_train
        else:
            X_va, y_va = [], []
        with open("megnet_default.log", "w+", encoding="utf-8") as logf:
            model = build_megnet(n1=64, n2=32, n3=16, lr=lr)
            model.model.summary(print_fn=lambda x: print(x, file=logf))
            start = datetime.datetime.now()
            _ = train_with_monitor(
                model,
                X_tr,
                y_tr,
                X_va,
                y_va,
                epochs,
                batch_size,
                patience,
                save_prefix="default",
                early_stop=True,
            )
            end = datetime.datetime.now()
            try:
                model.save_model("megnet_default.hdf5")
            except Exception:
                try:
                    model.model.save("megnet_default.h5")
                except Exception:
                    pass
            y_pred_train = predict_batch(model, X_tr)
            evaluate_and_dump("train", y_tr, y_pred_train, "scatter_train.dat", logf)
            if isinstance(X_va, list) and len(X_va) > 0:
                y_pred_val = predict_batch(model, X_va)
                evaluate_and_dump("vali", y_va, y_pred_val, "scatter_vali.dat", logf)
            y_pred_test = predict_batch(model, X_te)
            evaluate_and_dump("test", y_te, y_pred_test, "scatter_test.dat", logf)
            print(f"Runtime: {end - start}", file=logf)
            print(f"RuntimeSeconds: {seconds_str(end - start)}", file=logf)
        plot_scatter_with_marginals(
            "train", "scatter_train.dat", "scatter_plot_train.jpg"
        )
        if isinstance(X_va, list) and len(X_va) > 0:
            plot_scatter_with_marginals(
                "vali", "scatter_vali.dat", "scatter_plot_vali.jpg"
            )
        plot_scatter_with_marginals("test", "scatter_test.dat", "scatter_plot_test.jpg")
    else:
        k = int(n_folds)
        kf = KFold(n_splits=k, shuffle=True, random_state=random_state)
        fold_metrics = []
        start_all = datetime.datetime.now()
        for fi, (tr_idx, va_idx) in enumerate(kf.split(X_pool)):
            Xtr = [X_pool[i] for i in tr_idx]
            ytr = y_pool[tr_idx].ravel()
            Xva = [X_pool[i] for i in va_idx]
            yva = y_pool[va_idx].ravel()
            fold_dir = os.path.join(SCRIPT_DIR, f"fold_{fi}")
            os.makedirs(fold_dir, exist_ok=True)
            with open(
                os.path.join(fold_dir, f"megnet_default_fold_{fi}.log"),
                "w+",
                encoding="utf-8",
            ) as logf:
                model = build_megnet(n1=64, n2=32, n3=16, lr=lr)
                model.model.summary(print_fn=lambda x: print(x, file=logf))
                start = datetime.datetime.now()
                _ = train_with_monitor(
                    model,
                    Xtr,
                    ytr,
                    Xva,
                    yva,
                    epochs,
                    batch_size,
                    patience,
                    save_prefix=os.path.join(fold_dir, f"default_fold_{fi}"),
                    early_stop=True,
                )
                end = datetime.datetime.now()
                try:
                    model.save_model(
                        os.path.join(fold_dir, f"megnet_default_fold_{fi}.hdf5")
                    )
                except Exception:
                    try:
                        model.model.save(
                            os.path.join(fold_dir, f"megnet_default_fold_{fi}.h5")
                        )
                    except Exception:
                        pass
                tr, va, te = run_and_dump_fold(
                    model, fold_dir, "default", Xtr, ytr, Xva, yva, X_test, y_test
                )
                fold_metrics.append((fi, va[0], va[1], va[2], te[0], te[1], te[2]))
                print(f"Runtime: {end - start}", file=logf)
                print(f"RuntimeSeconds: {seconds_str(end - start)}", file=logf)
        end_all = datetime.datetime.now()
        with open("megnet_default_cv_summary.csv", "w+", encoding="utf-8") as f:
            f.write("fold,val_RMSE,val_MAE,val_R2,test_RMSE,test_MAE,test_R2\n")
            for m in fold_metrics:
                row = [str(m[0])] + [
                    f"{x:.6f}" if isinstance(x, float) else str(x) for x in m[1:]
                ]
                f.write(",".join(row) + "\n")
            val_rmse_avg = np.mean([m[1] for m in fold_metrics])
            f.write(f"avg,{val_rmse_avg:.6f},,,,\n")
        with open("megnet_default.log", "w+", encoding="utf-8") as logf:
            print(f"Runtime: {end_all - start_all}", file=logf)
            print(f"RuntimeSeconds: {seconds_str(end_all - start_all)}", file=logf)
        val_rmse_list = [m[1] for m in fold_metrics]
        best_idx = int(np.nanargmin(val_rmse_list))
        best_val_rmse = float(val_rmse_list[best_idx])
        print(f"Best fold: fold_{best_idx} (val RMSE={best_val_rmse:.6f})")
        with open("megnet_default.log", "a+", encoding="utf-8") as logf:
            print(
                f"Best fold: fold_{best_idx} (val RMSE={best_val_rmse:.6f})", file=logf
            )
        loss_files = [
            os.path.join(SCRIPT_DIR, f"fold_{i}", f"default_fold_{i}_loss.dat")
            for i in range(k)
        ]
        plot_cv_mean_std(loss_files, out_png="cv_loss_mean_std.jpg", ylabel="RMSE (eV)")
        tr_idx, va_idx = list(kf.split(X_pool))[best_idx]
        Xtr = [X_pool[i] for i in tr_idx]
        ytr = y_pool[tr_idx].ravel()
        Xva = [X_pool[i] for i in va_idx]
        yva = y_pool[va_idx].ravel()
        model_path = _first_exists(
            [
                os.path.join(
                    SCRIPT_DIR,
                    f"fold_{best_idx}",
                    f"megnet_default_fold_{best_idx}.hdf5",
                ),
                os.path.join(
                    SCRIPT_DIR, f"fold_{best_idx}", f"megnet_default_fold_{best_idx}.h5"
                ),
            ]
        )
        best_model = build_megnet(n1=64, n2=32, n3=16, lr=lr)
        if model_path:
            best_model = MEGNetModel.from_file(model_path)
        _save_split_csvs_and_log(
            best_model,
            Xtr,
            ytr,
            Xva,
            yva,
            X_test,
            y_test,
            SCRIPT_DIR,
            "best_fold",
            "megnet_default.log",
        )
else:
    print("please select method!")

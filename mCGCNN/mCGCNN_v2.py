from __future__ import annotations

import base64
import copy
import hashlib
import json
import logging
import math
import os
import random
import re
import time
import warnings
import zlib
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from pymatgen.core import Structure
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, Dataset


# ============================== User configuration ==============================

MODEL_NAME = "CGCNN"
RUN_VERSION = "v2"

EXCEL_PATH = "./data.xlsx"
CIF_DIR = "./cif"

TARGET_NAME = "bandgap"
TARGET_UNIT = "eV"

SEED = 42
TRAIN_RATIO = 0.80
VAL_RATIO = 0.10
TEST_RATIO = 0.10

RADIUS = 8.0
MAX_NUM_NEIGHBORS = 12
GAUSSIAN_DMIN = 0.0
GAUSSIAN_STEP = 0.2

ATOM_FEATURE_LENGTH = 64
N_CONV = 3
HIDDEN_FEATURE_LENGTH = 128
N_HIDDEN_LAYERS = 1

BATCH_SIZE = 64
MAX_EPOCHS = 300
LEARNING_RATE = 1.0e-3
WEIGHT_DECAY = 1.0e-5
PATIENCE = 40
MIN_DELTA = 1.0e-5
LR_PATIENCE = 12
LR_FACTOR = 0.5
MIN_LEARNING_RATE = 1.0e-6
NUM_WORKERS = 0
USE_AMP = True

USE_OPTUNA = False
OPTUNA_N_TRIALS = 30
OPTUNA_TIMEOUT = None
OPTUNA_MAX_EPOCHS = 120

FIG_DPI = 600
FONT_FAMILY = "Arial"
TITLE_FONTSIZE = 16
LABEL_FONTSIZE = 14
TICK_FONTSIZE = 12
LEGEND_FONTSIZE = 11
ANNOTATION_FONTSIZE = 10
USE_GRID = True

COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:red"]

OUTPUT_ROOT = Path(f"./{MODEL_NAME}_{RUN_VERSION}")
FIGURE_DIR = OUTPUT_ROOT / "figure"
DAT_DIR = OUTPUT_ROOT / "dat"
TABLE_DIR = OUTPUT_ROOT / "table"
LOG_DIR = OUTPUT_ROOT / "log"
SPLIT_DIR = OUTPUT_ROOT / "split"
CACHE_DIR = OUTPUT_ROOT / "cache"
CHECKPOINT_PATH = OUTPUT_ROOT / f"{MODEL_NAME}_best.pt"
SPLIT_PATH = SPLIT_DIR / f"{MODEL_NAME}_split.csv"
LOG_PATH = LOG_DIR / f"{MODEL_NAME}_training.log"

# ================================================================================


# The original CGCNN element descriptors are embedded to keep this file standalone.
_ATOM_FEATURES_B85 = (
    "c-rlp(T>|76hwa|+UI>>jInWlrRtYd{r4)g>UJBG0CxtQ)Y3jhb`+nz7lvW*_cy$KxqNl+7rb9w{U7Qtfb9!w??Ty?|9<}M{qi>b"
    "M$`6#(1m8VFUKb|-c&WY8_{|%NC>K`X6M5zn)#z@)Wz)DYFpN@_pZF62^9^6Ue6Cw2i^a{kI*c2V+}~{4Oz>S_tYaaYxOV;#^?a$"
    "fYQ+BdoxsQnW12N0AfWnv{f;ht#R^(Vblb5KCUgpVl;3jVdFH~MmnXJvw>CFSOf}d?4o=>ZrDs#X){XraZ7f$ZoF07SOgr*gNL0x"
    "vJ*IFtGH3Hk>Ww=Y5~EjZkitvBj)d_v3Rk{8|jisaZNd!)v9lrpMWl0jH2Uqo{Cm|qw~YMEiW?04K{+lX<kIMVTP>zKdwoBpCjm-"
    "W(lX3f9=LjVuMdTh3?WABZg@djm6WAQH+jk)cOX^m73V>>0<O+HbC>}RX;4Ri>Hh(#sSz+Y|}<dG_lF+VjO^tz3-fX&0SrL1F(^<"
    "W$&|jq>FXaO^3QP_Edo`hUQ<3aWbKhO_^Ve(RIUBiDhhjv0sdF;|)uISI*%c1mBO&rl;?h*M5E?_<o9g-2?c3-V=O3(Vn+x*-*}A"
    "Cis2>ZBdX0n7y^84+P&&5!}*63W|L{zYu&swe`oH?}{2Zu~|93AIl;=g70^|EU}3k-_KY^5EPRsDwNo49N*8P${?G-(LB|rq@UyX"
    "estScx`r?FL}ze(zah3LNNjQm(TU^x4Y5VVp6Kj3zMpO>j&W>ej_)V6DZswHkWJwDerkIQBiJk)-|w*YmZhgwj_=op4f@ti<oJHX"
    "yn1Kf@6H^)aeP0nzLDzUVviy?j_>zGwrDvUaC|?7_7hc4v40aLj_)TLi@W)LyKFqi_iL*FX^^DBBKn+73Ny#|)3Kp%RusWbYy!>q"
    "vz-3;bNi8ezp&7JKZ`AjE;Cf6?-!CTrWB3wP*u+pL65N|Ht}SWWwSjc8-+V?!)9?WUN?<9zljYbY}~=vRNtQ1O#emQwDpvK@w#bT"
    "{UtW@SvT3b39moipq0&{VuPpGH_Q6sY%Dz$Pc~UL+v~Gw*%D7LcyX^h7PnhFi4AD^8`c~?<dDd2-AuIn4L5#ufOB`Nr{!<x(KL)^"
    "Gt=@ntj+4f<ZsqM%ik#2hzhHa%|grHIN6&q>HB+yY*t$S#wp(#iV<ugEq~+0Z&W45lM2~vwERuO&mVOxv5ZZ0wET@|EN0=U*M1SS"
    "{0(u(p#S#WE}My#zv1qt=;BGmIg03M`5XF96GpR{Y55z)s&Q!kCI(vmMzo@)E2N@Z5}SpVzj3m6i%NsIci61V{EcgUv%EXlxw~I%"
    "A~S!}98Wq<{^qaTR;2pHM`L5#gV~&f%}eYvY&hBEbn!N+*f`}&=AHUfa$@7py2<v`e6q>S;ad>))U(&ka<a*?Sx+|Eo{A@%ESv3Q"
    "lVbztPB^=8e*6Udvl?m"
)


def embedded_atom_features() -> Dict[int, np.ndarray]:
    raw = zlib.decompress(base64.b85decode(_ATOM_FEATURES_B85.encode("ascii")))
    parsed = json.loads(raw.decode("utf-8"))
    return {
        int(key): np.asarray(value, dtype=np.float32) for key, value in parsed.items()
    }


ATOM_FEATURES = embedded_atom_features()
ORIGINAL_ATOM_FEATURE_LENGTH = len(next(iter(ATOM_FEATURES.values())))


def create_directories() -> None:
    for path in (
        OUTPUT_ROOT,
        FIGURE_DIR,
        DAT_DIR,
        TABLE_DIR,
        LOG_DIR,
        SPLIT_DIR,
        CACHE_DIR,
    ):
        path.mkdir(parents=True, exist_ok=True)


def setup_logger() -> logging.Logger:
    logger = logging.getLogger(MODEL_NAME)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    file_handler = logging.FileHandler(LOG_PATH, mode="w", encoding="utf-8")
    stream_handler = logging.StreamHandler()
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def set_global_seed(seed: int) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def format_duration(seconds: float) -> str:
    seconds = max(0, int(round(seconds)))
    hours, remainder = divmod(seconds, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def torch_load_compat(path: Path | str, map_location: Any = "cpu") -> Any:
    try:
        return torch.load(path, map_location=map_location, weights_only=False)
    except TypeError:
        return torch.load(path, map_location=map_location)


def normalize_cif_name(value: Any) -> str:
    if pd.isna(value):
        raise ValueError("empty CIF filename")
    if isinstance(value, (int, np.integer)):
        name = str(int(value))
    elif isinstance(value, (float, np.floating)) and float(value).is_integer():
        name = str(int(value))
    else:
        name = str(value).strip()
        if re.fullmatch(r"[+-]?\d+\.0+", name):
            name = str(int(float(name)))
    if not name:
        raise ValueError("empty CIF filename")
    if not name.lower().endswith(".cif"):
        name += ".cif"
    return name


@dataclass(frozen=True)
class Record:
    row_id: int
    cif_name: str
    cif_path: str
    target: float
    cache_path: str


class GaussianDistance:
    def __init__(
        self, dmin: float, dmax: float, step: float, var: Optional[float] = None
    ):
        if not (dmin < dmax and dmax - dmin > step):
            raise ValueError("Invalid Gaussian distance settings.")
        self.filter = np.arange(dmin, dmax + step, step, dtype=np.float32)
        self.var = float(step if var is None else var)

    def expand(self, distances: np.ndarray) -> np.ndarray:
        return np.exp(-((distances[..., np.newaxis] - self.filter) ** 2) / self.var**2)


def neighbor_index_and_distance(neighbor: Any) -> Tuple[int, float]:
    index = getattr(neighbor, "index", None)
    distance = getattr(neighbor, "nn_distance", None)
    if index is None:
        index = neighbor[2]
    if distance is None:
        distance = neighbor[1]
    return int(index), float(distance)


def build_crystal_graph(
    cif_path: Path | str,
    atom_features: Dict[int, np.ndarray] = ATOM_FEATURES,
    radius: float = RADIUS,
    max_num_neighbors: int = MAX_NUM_NEIGHBORS,
    dmin: float = GAUSSIAN_DMIN,
    step: float = GAUSSIAN_STEP,
) -> Dict[str, torch.Tensor]:
    structure = Structure.from_file(str(cif_path))
    if len(structure) == 0:
        raise ValueError("The structure contains no atomic sites.")
    if not structure.is_ordered:
        raise ValueError(
            "Disordered or partially occupied structures are not supported."
        )

    atomic_numbers: List[int] = []
    for site in structure:
        atomic_number = int(site.specie.Z)
        if atomic_number not in atom_features:
            raise ValueError(
                f"Atomic number {atomic_number} is absent from the embedded CGCNN descriptors."
            )
        atomic_numbers.append(atomic_number)
    atom_fea = np.vstack([atom_features[z] for z in atomic_numbers]).astype(np.float32)

    all_neighbors = structure.get_all_neighbors(radius, include_index=True)
    neighbor_indices: List[List[int]] = []
    neighbor_distances: List[List[float]] = []
    for neighbors in all_neighbors:
        parsed = sorted(
            (neighbor_index_and_distance(nbr) for nbr in neighbors),
            key=lambda item: item[1],
        )
        parsed = parsed[:max_num_neighbors]
        indices = [item[0] for item in parsed]
        distances = [item[1] for item in parsed]
        missing = max_num_neighbors - len(parsed)
        if missing:
            indices.extend([0] * missing)
            distances.extend([radius + 1.0] * missing)
        neighbor_indices.append(indices)
        neighbor_distances.append(distances)

    gdf = GaussianDistance(dmin=dmin, dmax=radius, step=step)
    nbr_fea = gdf.expand(np.asarray(neighbor_distances, dtype=np.float32)).astype(
        np.float32
    )
    nbr_fea_idx = np.asarray(neighbor_indices, dtype=np.int64)
    return {
        "atom_fea": torch.from_numpy(atom_fea),
        "nbr_fea": torch.from_numpy(nbr_fea),
        "nbr_fea_idx": torch.from_numpy(nbr_fea_idx),
    }


def cache_file_for(cif_path: Path) -> Path:
    stat = cif_path.stat()
    signature = (
        f"{cif_path.resolve()}|{stat.st_size}|{stat.st_mtime_ns}|{RADIUS}|"
        f"{MAX_NUM_NEIGHBORS}|{GAUSSIAN_DMIN}|{GAUSSIAN_STEP}|cgcnn-original-atom-features"
    )
    digest = hashlib.sha256(signature.encode("utf-8")).hexdigest()[:20]
    safe_stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", cif_path.stem)[:60]
    return CACHE_DIR / f"{safe_stem}_{digest}.pt"


def load_and_validate_records(
    logger: logging.Logger,
) -> Tuple[List[Record], List[Tuple[Any, str]]]:
    excel_path = Path(EXCEL_PATH)
    cif_dir = Path(CIF_DIR)
    if not excel_path.is_file():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    if not cif_dir.is_dir():
        raise NotADirectoryError(f"CIF directory not found: {cif_dir}")

    frame = pd.read_excel(excel_path)
    if frame.shape[1] < 2:
        raise ValueError("The Excel file must contain at least two columns.")

    records: List[Record] = []
    failures: List[Tuple[Any, str]] = []
    for row_id, row in frame.iterrows():
        raw_name = row.iloc[0]
        raw_target = row.iloc[1]
        try:
            cif_name = normalize_cif_name(raw_name)
            target = float(raw_target)
            if not math.isfinite(target):
                raise ValueError("target is NaN or infinite")
            cif_path = cif_dir / cif_name
            if not cif_path.is_file():
                raise FileNotFoundError(f"CIF file not found: {cif_path}")
            cache_path = cache_file_for(cif_path)
            if not cache_path.is_file():
                graph = build_crystal_graph(cif_path)
                torch.save(graph, cache_path)
            else:
                graph = torch_load_compat(cache_path)
                required = {"atom_fea", "nbr_fea", "nbr_fea_idx"}
                if not required.issubset(graph):
                    graph = build_crystal_graph(cif_path)
                    torch.save(graph, cache_path)
            records.append(
                Record(
                    row_id=int(row_id),
                    cif_name=cif_name,
                    cif_path=str(cif_path),
                    target=target,
                    cache_path=str(cache_path),
                )
            )
        except Exception as exc:
            failures.append((raw_name, str(exc)))
            logger.warning("Skipped row %s (%s): %s", row_id, raw_name, exc)

    if len(records) < 10:
        raise ValueError(
            f"Only {len(records)} valid samples remain; at least 10 are required."
        )
    return records, failures


def assign_splits(
    records: Sequence[Record], logger: logging.Logger
) -> Dict[str, List[Record]]:
    ratios = np.asarray([TRAIN_RATIO, VAL_RATIO, TEST_RATIO], dtype=float)
    if not np.isclose(ratios.sum(), 1.0) or np.any(ratios <= 0):
        raise ValueError(
            "TRAIN_RATIO, VAL_RATIO, and TEST_RATIO must be positive and sum to 1."
        )

    if SPLIT_PATH.is_file():
        existing = pd.read_csv(SPLIT_PATH)
        required = {"row_id", "cif_name", "target", "split"}
        if required.issubset(existing.columns) and len(existing) == len(records):
            names_match = existing["cif_name"].astype(str).tolist() == [
                r.cif_name for r in records
            ]
            targets_match = np.allclose(
                existing["target"].to_numpy(float), [r.target for r in records]
            )
            splits_valid = set(existing["split"].astype(str)) == {
                "train",
                "val",
                "test",
            }
            if names_match and targets_match and splits_valid:
                split_labels = existing["split"].astype(str).tolist()
                logger.info("Reusing split file: %s", SPLIT_PATH)
                return {
                    name: [
                        record
                        for record, label in zip(records, split_labels)
                        if label == name
                    ]
                    for name in ("train", "val", "test")
                }
        logger.warning(
            "Existing split file does not match the current valid dataset and will be replaced."
        )

    indices = np.arange(len(records))
    train_idx, holdout_idx = train_test_split(
        indices,
        train_size=TRAIN_RATIO,
        random_state=SEED,
        shuffle=True,
    )
    relative_test_ratio = TEST_RATIO / (VAL_RATIO + TEST_RATIO)
    val_idx, test_idx = train_test_split(
        holdout_idx,
        test_size=relative_test_ratio,
        random_state=SEED,
        shuffle=True,
    )
    label_by_index = {int(i): "train" for i in train_idx}
    label_by_index.update({int(i): "val" for i in val_idx})
    label_by_index.update({int(i): "test" for i in test_idx})
    rows = [
        {
            "row_id": record.row_id,
            "cif_name": record.cif_name,
            "target": record.target,
            "split": label_by_index[i],
        }
        for i, record in enumerate(records)
    ]
    pd.DataFrame(rows).to_csv(SPLIT_PATH, index=False)
    logger.info("Created split file: %s", SPLIT_PATH)
    return {
        name: [record for i, record in enumerate(records) if label_by_index[i] == name]
        for name in ("train", "val", "test")
    }


class CrystalDataset(Dataset):
    def __init__(self, records: Sequence[Record]):
        self.records = list(records)

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(
        self, index: int
    ) -> Tuple[Dict[str, torch.Tensor], torch.Tensor, str]:
        record = self.records[index]
        graph = torch_load_compat(record.cache_path)
        target = torch.tensor([record.target], dtype=torch.float32)
        return graph, target, record.cif_name


def collate_crystals(
    batch: Sequence[Tuple[Dict[str, torch.Tensor], torch.Tensor, str]]
) -> Tuple[
    Tuple[torch.Tensor, torch.Tensor, torch.Tensor, List[torch.Tensor]],
    torch.Tensor,
    List[str],
]:
    atom_features: List[torch.Tensor] = []
    neighbor_features: List[torch.Tensor] = []
    neighbor_indices: List[torch.Tensor] = []
    crystal_atom_indices: List[torch.Tensor] = []
    targets: List[torch.Tensor] = []
    names: List[str] = []
    base_index = 0

    for graph, target, name in batch:
        atom_fea = graph["atom_fea"].float()
        nbr_fea = graph["nbr_fea"].float()
        nbr_fea_idx = graph["nbr_fea_idx"].long()
        atom_count = atom_fea.shape[0]
        atom_features.append(atom_fea)
        neighbor_features.append(nbr_fea)
        neighbor_indices.append(nbr_fea_idx + base_index)
        crystal_atom_indices.append(
            torch.arange(base_index, base_index + atom_count, dtype=torch.long)
        )
        targets.append(target.float())
        names.append(name)
        base_index += atom_count

    return (
        (
            torch.cat(atom_features, dim=0),
            torch.cat(neighbor_features, dim=0),
            torch.cat(neighbor_indices, dim=0),
            crystal_atom_indices,
        ),
        torch.stack(targets, dim=0).view(-1),
        names,
    )


class SafeBatchNorm1d(nn.BatchNorm1d):
    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        if self.training and inputs.shape[0] == 1:
            return F.batch_norm(
                inputs,
                self.running_mean,
                self.running_var,
                self.weight,
                self.bias,
                False,
                self.momentum,
                self.eps,
            )
        return super().forward(inputs)


class ConvLayer(nn.Module):
    def __init__(self, atom_fea_len: int, nbr_fea_len: int):
        super().__init__()
        self.atom_fea_len = atom_fea_len
        self.fc_full = nn.Linear(2 * atom_fea_len + nbr_fea_len, 2 * atom_fea_len)
        self.bn1 = SafeBatchNorm1d(2 * atom_fea_len)
        self.bn2 = SafeBatchNorm1d(atom_fea_len)
        self.softplus = nn.Softplus()

    def forward(
        self,
        atom_in_fea: torch.Tensor,
        nbr_fea: torch.Tensor,
        nbr_fea_idx: torch.Tensor,
    ) -> torch.Tensor:
        atom_count, max_neighbors = nbr_fea_idx.shape
        atom_nbr_fea = atom_in_fea[nbr_fea_idx]
        center_fea = atom_in_fea.unsqueeze(1).expand(
            atom_count, max_neighbors, self.atom_fea_len
        )
        total_nbr_fea = torch.cat([center_fea, atom_nbr_fea, nbr_fea], dim=2)
        total_gated_fea = self.fc_full(total_nbr_fea)
        total_gated_fea = self.bn1(total_gated_fea.reshape(-1, 2 * self.atom_fea_len))
        total_gated_fea = total_gated_fea.reshape(
            atom_count, max_neighbors, 2 * self.atom_fea_len
        )
        nbr_filter, nbr_core = total_gated_fea.chunk(2, dim=2)
        nbr_sum = torch.sum(torch.sigmoid(nbr_filter) * self.softplus(nbr_core), dim=1)
        nbr_sum = self.bn2(nbr_sum)
        return self.softplus(atom_in_fea + nbr_sum)


class CrystalGraphConvNet(nn.Module):
    def __init__(
        self,
        orig_atom_fea_len: int,
        nbr_fea_len: int,
        atom_fea_len: int = 64,
        n_conv: int = 3,
        h_fea_len: int = 128,
        n_h: int = 1,
    ):
        super().__init__()
        self.embedding = nn.Linear(orig_atom_fea_len, atom_fea_len)
        self.convs = nn.ModuleList(
            [ConvLayer(atom_fea_len, nbr_fea_len) for _ in range(n_conv)]
        )
        self.conv_to_fc = nn.Linear(atom_fea_len, h_fea_len)
        self.softplus = nn.Softplus()
        self.fcs = nn.ModuleList(
            [nn.Linear(h_fea_len, h_fea_len) for _ in range(max(0, n_h - 1))]
        )
        self.fc_out = nn.Linear(h_fea_len, 1)

    @staticmethod
    def pooling(
        atom_fea: torch.Tensor, crystal_atom_idx: Sequence[torch.Tensor]
    ) -> torch.Tensor:
        if sum(len(indices) for indices in crystal_atom_idx) != atom_fea.shape[0]:
            raise RuntimeError(
                "Crystal-to-atom indices do not match the batched atom tensor."
            )
        return torch.cat(
            [
                torch.mean(atom_fea[indices], dim=0, keepdim=True)
                for indices in crystal_atom_idx
            ],
            dim=0,
        )

    def forward(
        self,
        atom_fea: torch.Tensor,
        nbr_fea: torch.Tensor,
        nbr_fea_idx: torch.Tensor,
        crystal_atom_idx: Sequence[torch.Tensor],
    ) -> torch.Tensor:
        atom_fea = self.embedding(atom_fea)
        for conv in self.convs:
            atom_fea = conv(atom_fea, nbr_fea, nbr_fea_idx)
        crystal_fea = self.pooling(atom_fea, crystal_atom_idx)
        crystal_fea = self.softplus(self.conv_to_fc(self.softplus(crystal_fea)))
        for layer in self.fcs:
            crystal_fea = self.softplus(layer(crystal_fea))
        return self.fc_out(crystal_fea).view(-1)


@dataclass
class TargetNormalizer:
    mean: float
    std: float

    @classmethod
    def from_targets(cls, targets: Sequence[float]) -> "TargetNormalizer":
        values = np.asarray(targets, dtype=np.float64)
        std = float(values.std())
        if not math.isfinite(std) or std < 1.0e-12:
            std = 1.0
        return cls(mean=float(values.mean()), std=std)

    def normalize(self, values: torch.Tensor) -> torch.Tensor:
        return (values - self.mean) / self.std

    def denormalize(self, values: torch.Tensor) -> torch.Tensor:
        return values * self.std + self.mean

    def state_dict(self) -> Dict[str, float]:
        return {"mean": self.mean, "std": self.std}


def make_model(model_config: Dict[str, Any]) -> CrystalGraphConvNet:
    return CrystalGraphConvNet(**model_config)


def make_loaders(
    splits: Dict[str, List[Record]], batch_size: int
) -> Dict[str, DataLoader]:
    generator = torch.Generator()
    generator.manual_seed(SEED)
    common = {
        "batch_size": batch_size,
        "num_workers": NUM_WORKERS,
        "collate_fn": collate_crystals,
        "pin_memory": torch.cuda.is_available(),
        "persistent_workers": NUM_WORKERS > 0,
    }
    return {
        "train": DataLoader(
            CrystalDataset(splits["train"]), shuffle=True, generator=generator, **common
        ),
        "train_eval": DataLoader(
            CrystalDataset(splits["train"]), shuffle=False, **common
        ),
        "val": DataLoader(CrystalDataset(splits["val"]), shuffle=False, **common),
        "test": DataLoader(CrystalDataset(splits["test"]), shuffle=False, **common),
    }


def move_batch(
    inputs: Tuple[torch.Tensor, torch.Tensor, torch.Tensor, List[torch.Tensor]],
    targets: torch.Tensor,
    device: torch.device,
) -> Tuple[
    Tuple[torch.Tensor, torch.Tensor, torch.Tensor, List[torch.Tensor]], torch.Tensor
]:
    atom_fea, nbr_fea, nbr_fea_idx, crystal_atom_idx = inputs
    moved = (
        atom_fea.to(device, non_blocking=True),
        nbr_fea.to(device, non_blocking=True),
        nbr_fea_idx.to(device, non_blocking=True),
        [indices.to(device, non_blocking=True) for indices in crystal_atom_idx],
    )
    return moved, targets.to(device, non_blocking=True)


def amp_settings(device: torch.device) -> Tuple[bool, torch.dtype, Any]:
    enabled = bool(USE_AMP and device.type == "cuda")
    dtype = (
        torch.bfloat16 if enabled and torch.cuda.is_bf16_supported() else torch.float16
    )
    if enabled and dtype == torch.float16:
        try:
            scaler = torch.amp.GradScaler("cuda", enabled=True)
        except (AttributeError, TypeError):
            scaler = torch.cuda.amp.GradScaler(enabled=True)
    else:
        scaler = None
    return enabled, dtype, scaler


def autocast_context(device: torch.device, enabled: bool, dtype: torch.dtype):
    if enabled:
        return torch.autocast(device_type="cuda", dtype=dtype)
    return nullcontext()


def train_one_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    normalizer: TargetNormalizer,
    device: torch.device,
    amp_enabled: bool,
    amp_dtype: torch.dtype,
    scaler: Any,
) -> float:
    model.train()
    total_loss = 0.0
    total_samples = 0
    for inputs, targets, _ in loader:
        inputs, targets = move_batch(inputs, targets, device)
        normalized_targets = normalizer.normalize(targets).view(-1)
        optimizer.zero_grad(set_to_none=True)
        with autocast_context(device, amp_enabled, amp_dtype):
            predictions = model(*inputs).view(-1)
            loss = F.mse_loss(predictions, normalized_targets)
        if not torch.isfinite(loss):
            raise FloatingPointError("A non-finite training loss was encountered.")
        if scaler is not None:
            scaler.scale(loss).backward()
            scaler.unscale_(optimizer)
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            scaler.step(optimizer)
            scaler.update()
        else:
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=5.0)
            optimizer.step()
        count = targets.numel()
        total_loss += float(loss.detach().cpu()) * count
        total_samples += count
    return total_loss / max(total_samples, 1)


@torch.inference_mode()
def evaluate(
    model: nn.Module,
    loader: DataLoader,
    normalizer: TargetNormalizer,
    device: torch.device,
    amp_enabled: bool,
    amp_dtype: torch.dtype,
) -> Tuple[np.ndarray, np.ndarray, List[str], float]:
    model.eval()
    true_values: List[float] = []
    predictions: List[float] = []
    names: List[str] = []
    normalized_mse_sum = 0.0
    sample_count = 0
    for inputs, targets, batch_names in loader:
        inputs, targets = move_batch(inputs, targets, device)
        with autocast_context(device, amp_enabled, amp_dtype):
            normalized_predictions = model(*inputs).view(-1)
            normalized_targets = normalizer.normalize(targets).view(-1)
            loss = F.mse_loss(
                normalized_predictions, normalized_targets, reduction="sum"
            )
        physical_predictions = normalizer.denormalize(normalized_predictions.float())
        true_values.extend(targets.detach().cpu().numpy().astype(float).tolist())
        predictions.extend(
            physical_predictions.detach().cpu().numpy().astype(float).tolist()
        )
        names.extend(batch_names)
        normalized_mse_sum += float(loss.detach().cpu())
        sample_count += targets.numel()
    return (
        np.asarray(true_values, dtype=float),
        np.asarray(predictions, dtype=float),
        names,
        normalized_mse_sum / max(sample_count, 1),
    )


def calculate_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    return {
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(math.sqrt(mean_squared_error(y_true, y_pred))),
        "r2": float(r2_score(y_true, y_pred)),
    }


def model_and_graph_config(
    hparams: Optional[Dict[str, Any]] = None
) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    hparams = hparams or {}
    nbr_fea_len = len(np.arange(GAUSSIAN_DMIN, RADIUS + GAUSSIAN_STEP, GAUSSIAN_STEP))
    model_config = {
        "orig_atom_fea_len": ORIGINAL_ATOM_FEATURE_LENGTH,
        "nbr_fea_len": nbr_fea_len,
        "atom_fea_len": int(hparams.get("atom_fea_len", ATOM_FEATURE_LENGTH)),
        "n_conv": int(hparams.get("n_conv", N_CONV)),
        "h_fea_len": int(hparams.get("h_fea_len", HIDDEN_FEATURE_LENGTH)),
        "n_h": int(hparams.get("n_h", N_HIDDEN_LAYERS)),
    }
    graph_config = {
        "radius": RADIUS,
        "max_num_neighbors": MAX_NUM_NEIGHBORS,
        "dmin": GAUSSIAN_DMIN,
        "step": GAUSSIAN_STEP,
        "atom_feature_source": "embedded original CGCNN atom_init.json",
        "atom_feature_length": ORIGINAL_ATOM_FEATURE_LENGTH,
    }
    return model_config, graph_config


def save_checkpoint(
    model: nn.Module,
    model_config: Dict[str, Any],
    graph_config: Dict[str, Any],
    normalizer: TargetNormalizer,
    best_epoch: int,
    best_val_rmse: float,
    history: List[Dict[str, float]],
    hparams: Dict[str, Any],
) -> None:
    checkpoint = {
        "model_name": MODEL_NAME,
        "run_version": RUN_VERSION,
        "model_state_dict": copy.deepcopy(model.state_dict()),
        "model_config": model_config,
        "graph_config": graph_config,
        "normalization": normalizer.state_dict(),
        "best_epoch": best_epoch,
        "best_val_rmse": best_val_rmse,
        "history": history,
        "hyperparameters": hparams,
        "seed": SEED,
        "target_name": TARGET_NAME,
        "target_unit": TARGET_UNIT,
        "atom_features": {key: value.tolist() for key, value in ATOM_FEATURES.items()},
    }
    torch.save(checkpoint, CHECKPOINT_PATH)


def train_model(
    loaders: Dict[str, DataLoader],
    normalizer: TargetNormalizer,
    device: torch.device,
    logger: logging.Logger,
    hparams: Dict[str, Any],
    max_epochs: int,
    save_best: bool,
    verbose: bool,
) -> Tuple[nn.Module, List[Dict[str, float]], int, float]:
    model_config, graph_config = model_and_graph_config(hparams)
    model = make_model(model_config).to(device)
    learning_rate = float(hparams.get("learning_rate", LEARNING_RATE))
    weight_decay = float(hparams.get("weight_decay", WEIGHT_DECAY))
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=learning_rate, weight_decay=weight_decay
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="min",
        factor=LR_FACTOR,
        patience=LR_PATIENCE,
        min_lr=MIN_LEARNING_RATE,
    )
    amp_enabled, amp_dtype, scaler = amp_settings(device)
    history: List[Dict[str, float]] = []
    best_val_rmse = float("inf")
    best_epoch = 0
    best_state: Optional[Dict[str, torch.Tensor]] = None
    stale_epochs = 0

    for epoch in range(1, max_epochs + 1):
        train_mse_loss = train_one_epoch(
            model,
            loaders["train"],
            optimizer,
            normalizer,
            device,
            amp_enabled,
            amp_dtype,
            scaler,
        )
        train_true, train_pred, _, _ = evaluate(
            model, loaders["train_eval"], normalizer, device, amp_enabled, amp_dtype
        )
        val_true, val_pred, _, val_mse_loss = evaluate(
            model, loaders["val"], normalizer, device, amp_enabled, amp_dtype
        )
        train_rmse = calculate_metrics(train_true, train_pred)["rmse"]
        val_rmse = calculate_metrics(val_true, val_pred)["rmse"]
        scheduler.step(val_rmse)
        current_lr = float(optimizer.param_groups[0]["lr"])
        history.append(
            {
                "epoch": float(epoch),
                "train_rmse_eV": train_rmse,
                "val_rmse_eV": val_rmse,
                "learning_rate": current_lr,
                "train_mse_loss": train_mse_loss,
                "val_mse_loss": val_mse_loss,
            }
        )
        if verbose:
            logger.info(
                "Epoch %04d | train RMSE %.6f eV | val RMSE %.6f eV | lr %.3e",
                epoch,
                train_rmse,
                val_rmse,
                current_lr,
            )
        if val_rmse < best_val_rmse - MIN_DELTA:
            best_val_rmse = val_rmse
            best_epoch = epoch
            best_state = copy.deepcopy(model.state_dict())
            stale_epochs = 0
            if save_best:
                save_checkpoint(
                    model,
                    model_config,
                    graph_config,
                    normalizer,
                    best_epoch,
                    best_val_rmse,
                    history,
                    hparams,
                )
        else:
            stale_epochs += 1
        if stale_epochs >= PATIENCE:
            if verbose:
                logger.info("Early stopping triggered at epoch %d.", epoch)
            break

    if best_state is None:
        raise RuntimeError("Training ended without a valid checkpoint state.")
    model.load_state_dict(best_state)
    if save_best:
        save_checkpoint(
            model,
            model_config,
            graph_config,
            normalizer,
            best_epoch,
            best_val_rmse,
            history,
            hparams,
        )
    return model, history, best_epoch, best_val_rmse


def run_optuna(
    splits: Dict[str, List[Record]],
    normalizer: TargetNormalizer,
    device: torch.device,
    logger: logging.Logger,
) -> Dict[str, Any]:
    try:
        import optuna
    except ImportError as exc:
        raise ImportError("USE_OPTUNA=True requires the optuna package.") from exc

    optuna.logging.set_verbosity(optuna.logging.WARNING)

    def objective(trial: Any) -> float:
        set_global_seed(SEED)
        params = {
            "atom_fea_len": trial.suggest_categorical("atom_fea_len", [64, 96, 128]),
            "n_conv": trial.suggest_int("n_conv", 3, 5),
            "h_fea_len": trial.suggest_categorical("h_fea_len", [128, 192, 256]),
            "n_h": trial.suggest_int("n_h", 1, 3),
            "learning_rate": trial.suggest_float(
                "learning_rate", 2.0e-4, 3.0e-3, log=True
            ),
            "weight_decay": trial.suggest_float(
                "weight_decay", 1.0e-7, 1.0e-3, log=True
            ),
        }
        trial_loaders = make_loaders(splits, BATCH_SIZE)
        _, _, _, best_rmse = train_model(
            trial_loaders,
            normalizer,
            device,
            logger,
            params,
            OPTUNA_MAX_EPOCHS,
            save_best=False,
            verbose=False,
        )
        return best_rmse

    study = optuna.create_study(
        direction="minimize", sampler=optuna.samplers.TPESampler(seed=SEED)
    )
    study.optimize(objective, n_trials=OPTUNA_N_TRIALS, timeout=OPTUNA_TIMEOUT)
    result = {"best_value": float(study.best_value), "best_params": study.best_params}
    output_path = TABLE_DIR / f"{MODEL_NAME}_optuna_best_params.json"
    output_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    logger.info("Optuna best validation RMSE: %.6f eV", study.best_value)
    logger.info("Optuna best parameters: %s", study.best_params)
    return dict(study.best_params)


def load_trained_model(
    checkpoint_path: Path | str,
    device: Optional[torch.device] = None,
) -> Tuple[
    CrystalGraphConvNet, TargetNormalizer, Dict[str, Any], Dict[int, np.ndarray]
]:
    device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint = torch_load_compat(checkpoint_path, map_location=device)
    if checkpoint.get("model_name") != MODEL_NAME:
        warnings.warn(
            f"Checkpoint model name is {checkpoint.get('model_name')}, expected {MODEL_NAME}."
        )
    model = make_model(checkpoint["model_config"]).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    normalizer = TargetNormalizer(**checkpoint["normalization"])
    atom_features = {
        int(key): np.asarray(value, dtype=np.float32)
        for key, value in checkpoint.get("atom_features", ATOM_FEATURES).items()
    }
    return model, normalizer, checkpoint, atom_features


@torch.inference_mode()
def predict_cifs(
    checkpoint_path: Path | str,
    cif_paths: Sequence[Path | str],
    device: Optional[torch.device] = None,
) -> pd.DataFrame:
    device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, normalizer, checkpoint, atom_features = load_trained_model(
        checkpoint_path, device
    )
    graph_config = checkpoint["graph_config"]
    items = []
    for cif_path in cif_paths:
        path = Path(cif_path)
        graph = build_crystal_graph(
            path,
            atom_features=atom_features,
            radius=float(graph_config["radius"]),
            max_num_neighbors=int(graph_config["max_num_neighbors"]),
            dmin=float(graph_config["dmin"]),
            step=float(graph_config["step"]),
        )
        items.append((graph, torch.tensor([0.0], dtype=torch.float32), path.name))
    inputs, _, names = collate_crystals(items)
    inputs, _ = move_batch(inputs, torch.zeros(len(items)), device)
    amp_enabled, amp_dtype, _ = amp_settings(device)
    with autocast_context(device, amp_enabled, amp_dtype):
        normalized_predictions = model(*inputs).view(-1)
    predictions = normalizer.denormalize(normalized_predictions.float()).cpu().numpy()
    return pd.DataFrame(
        {"cif_name": names, f"predicted_{TARGET_NAME}_{TARGET_UNIT}": predictions}
    )


def save_prediction_dat(
    split: str, names: Sequence[str], y_true: np.ndarray, y_pred: np.ndarray
) -> None:
    frame = pd.DataFrame(
        {
            "cif_name": names,
            "true_bandgap_eV": y_true,
            "predicted_bandgap_eV": y_pred,
            "error_eV": y_pred - y_true,
            "absolute_error_eV": np.abs(y_pred - y_true),
        }
    )
    frame.to_csv(
        DAT_DIR / f"{MODEL_NAME}_parity_{split}.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )


def apply_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [FONT_FAMILY, "DejaVu Sans"],
            "axes.titlesize": TITLE_FONTSIZE,
            "axes.labelsize": LABEL_FONTSIZE,
            "xtick.labelsize": TICK_FONTSIZE,
            "ytick.labelsize": TICK_FONTSIZE,
            "legend.fontsize": LEGEND_FONTSIZE,
        }
    )


def parity_limits(
    all_true: Iterable[np.ndarray], all_pred: Iterable[np.ndarray]
) -> Tuple[float, float]:
    values = np.concatenate(
        [np.asarray(x, dtype=float) for x in [*all_true, *all_pred]]
    )
    low = float(np.min(values))
    high = float(np.max(values))
    margin = max(0.05 * (high - low), 0.05)
    return low - margin, high + margin


def plot_parity(
    split: str,
    y_true: np.ndarray,
    y_pred: np.ndarray,
    metrics: Dict[str, float],
    limits: Tuple[float, float],
    color: str,
) -> None:
    figure, axis = plt.subplots(figsize=(6.2, 6.0))
    axis.scatter(y_true, y_pred, s=24, alpha=0.78, color=color, edgecolors="none")
    axis.plot(limits, limits, linestyle="--", color="black", linewidth=1.2)
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})")
    axis.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})")
    axis.set_title(f"{MODEL_NAME} {split.capitalize()} Parity")
    annotation = (
        f"MAE = {metrics['mae']:.4f} eV\n"
        f"RMSE = {metrics['rmse']:.4f} eV\n"
        f"R² = {metrics['r2']:.4f}"
    )
    axis.text(
        0.04,
        0.96,
        annotation,
        transform=axis.transAxes,
        va="top",
        fontsize=ANNOTATION_FONTSIZE,
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "alpha": 0.85,
            "edgecolor": "0.75",
        },
    )
    if USE_GRID:
        axis.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)
    figure.tight_layout()
    figure.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_parity_{split}.jpg",
        dpi=FIG_DPI,
        bbox_inches="tight",
    )
    plt.close(figure)


def plot_combined_parity(
    predictions: Dict[str, Tuple[np.ndarray, np.ndarray, List[str]]],
    metrics: Dict[str, Dict[str, float]],
    limits: Tuple[float, float],
) -> None:
    figure, axis = plt.subplots(figsize=(6.5, 6.2))
    for index, split in enumerate(("train", "val", "test")):
        y_true, y_pred, _ = predictions[split]
        axis.scatter(
            y_true,
            y_pred,
            s=22,
            alpha=0.72,
            color=COLORS[index],
            edgecolors="none",
            label=f"{split.capitalize()} (RMSE={metrics[split]['rmse']:.4f})",
        )
    axis.plot(
        limits, limits, linestyle="--", color="black", linewidth=1.2, label="y = x"
    )
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})")
    axis.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})")
    axis.set_title(f"{MODEL_NAME} Combined Parity")
    axis.legend()
    if USE_GRID:
        axis.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)
    figure.tight_layout()
    figure.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_parity_all.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(figure)


def save_combined_dat(
    predictions: Dict[str, Tuple[np.ndarray, np.ndarray, List[str]]]
) -> None:
    frames = []
    for split in ("train", "val", "test"):
        y_true, y_pred, names = predictions[split]
        frames.append(
            pd.DataFrame(
                {
                    "split": split,
                    "cif_name": names,
                    "true_bandgap_eV": y_true,
                    "predicted_bandgap_eV": y_pred,
                    "error_eV": y_pred - y_true,
                    "absolute_error_eV": np.abs(y_pred - y_true),
                }
            )
        )
    pd.concat(frames, ignore_index=True).to_csv(
        DAT_DIR / f"{MODEL_NAME}_parity_all.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )


def save_and_plot_history(history: Sequence[Dict[str, float]], best_epoch: int) -> None:
    frame = pd.DataFrame(history)
    frame["epoch"] = frame["epoch"].astype(int)
    frame.to_csv(
        DAT_DIR / f"{MODEL_NAME}_rmse_curve.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )

    figure, axis = plt.subplots(figsize=(7.2, 5.2))
    axis.plot(
        frame["epoch"], frame["train_rmse_eV"], color=COLORS[0], label="Train RMSE"
    )
    axis.plot(
        frame["epoch"], frame["val_rmse_eV"], color=COLORS[1], label="Validation RMSE"
    )
    axis.axvline(
        best_epoch,
        color=COLORS[3],
        linestyle="--",
        linewidth=1.2,
        label=f"Best epoch: {best_epoch}",
    )
    axis.set_xlabel("Epoch")
    axis.set_ylabel(f"RMSE ({TARGET_UNIT})")
    axis.set_title(f"{MODEL_NAME} Training History")
    axis.legend()
    if USE_GRID:
        axis.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)
    figure.tight_layout()
    figure.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_rmse_curve.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(figure)


def save_metrics_table(
    metrics: Dict[str, Dict[str, float]], predictions: Dict[str, Any]
) -> None:
    rows = []
    for split in ("train", "val", "test"):
        rows.append(
            {
                "split": split,
                "n_samples": len(predictions[split][0]),
                "mae_eV": metrics[split]["mae"],
                "rmse_eV": metrics[split]["rmse"],
                "r2": metrics[split]["r2"],
            }
        )
    pd.DataFrame(rows).to_csv(
        TABLE_DIR / f"{MODEL_NAME}_metrics.dat",
        sep="\t",
        index=False,
        float_format="%.6f",
    )


def log_environment(logger: logging.Logger, device: torch.device) -> None:
    logger.info("Model name: %s", MODEL_NAME)
    logger.info("Run version: %s", RUN_VERSION)
    logger.info("Output root: %s", OUTPUT_ROOT)
    logger.info("Excel path: %s", EXCEL_PATH)
    logger.info("CIF directory: %s", CIF_DIR)
    logger.info("Excel columns: first column = CIF filename; second column = target")
    logger.info("Target: %s (%s)", TARGET_NAME, TARGET_UNIT)
    logger.info("Seed: %d", SEED)
    logger.info(
        "Split ratios: train=%.2f, val=%.2f, test=%.2f",
        TRAIN_RATIO,
        VAL_RATIO,
        TEST_RATIO,
    )
    logger.info("Device: %s", device)
    logger.info("PyTorch version: %s", torch.__version__)
    if device.type == "cuda":
        props = torch.cuda.get_device_properties(device)
        logger.info("GPU: %s", torch.cuda.get_device_name(device))
        logger.info("GPU memory: %.2f GiB", props.total_memory / 1024**3)
        logger.info(
            "AMP requested: %s; BF16 supported: %s",
            USE_AMP,
            torch.cuda.is_bf16_supported(),
        )
    logger.info("Loss function: MSELoss")
    logger.info("Optimizer: AdamW")


def main() -> None:
    total_start = time.perf_counter()
    create_directories()
    logger = setup_logger()
    set_global_seed(SEED)
    apply_plot_style()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    log_environment(logger, device)

    data_start = time.perf_counter()
    records, failures = load_and_validate_records(logger)
    splits = assign_splits(records, logger)
    data_seconds = time.perf_counter() - data_start
    logger.info("Valid samples: %d; invalid samples: %d", len(records), len(failures))
    logger.info(
        "Split sizes: train=%d, val=%d, test=%d",
        len(splits["train"]),
        len(splits["val"]),
        len(splits["test"]),
    )
    logger.info(
        "Data preparation time: %.2f s (%s)",
        data_seconds,
        format_duration(data_seconds),
    )

    normalizer = TargetNormalizer.from_targets(
        [record.target for record in splits["train"]]
    )
    logger.info(
        "Training target normalization: mean=%.8f, std=%.8f",
        normalizer.mean,
        normalizer.std,
    )

    hparams: Dict[str, Any] = {
        "atom_fea_len": ATOM_FEATURE_LENGTH,
        "n_conv": N_CONV,
        "h_fea_len": HIDDEN_FEATURE_LENGTH,
        "n_h": N_HIDDEN_LAYERS,
        "learning_rate": LEARNING_RATE,
        "weight_decay": WEIGHT_DECAY,
    }
    if USE_OPTUNA:
        hparams.update(run_optuna(splits, normalizer, device, logger))
    else:
        logger.info("Optuna enabled: False")

    loaders = make_loaders(splits, BATCH_SIZE)
    train_start = time.perf_counter()
    model, history, best_epoch, best_val_rmse = train_model(
        loaders,
        normalizer,
        device,
        logger,
        hparams,
        MAX_EPOCHS,
        save_best=True,
        verbose=True,
    )
    train_seconds = time.perf_counter() - train_start
    logger.info("Best epoch: %d", best_epoch)
    logger.info("Best validation RMSE: %.6f eV", best_val_rmse)
    logger.info(
        "Training time: %.2f s (%s)", train_seconds, format_duration(train_seconds)
    )

    eval_start = time.perf_counter()
    model, normalizer, checkpoint, _ = load_trained_model(CHECKPOINT_PATH, device)
    amp_enabled, amp_dtype, _ = amp_settings(device)
    predictions: Dict[str, Tuple[np.ndarray, np.ndarray, List[str]]] = {}
    metrics: Dict[str, Dict[str, float]] = {}
    loader_names = {"train": "train_eval", "val": "val", "test": "test"}
    for split, loader_name in loader_names.items():
        y_true, y_pred, names, _ = evaluate(
            model, loaders[loader_name], normalizer, device, amp_enabled, amp_dtype
        )
        predictions[split] = (y_true, y_pred, names)
        metrics[split] = calculate_metrics(y_true, y_pred)
        logger.info(
            "%s metrics | MAE %.6f eV | RMSE %.6f eV | R2 %.6f",
            split.capitalize(),
            metrics[split]["mae"],
            metrics[split]["rmse"],
            metrics[split]["r2"],
        )
        save_prediction_dat(split, names, y_true, y_pred)

    limits = parity_limits(
        [predictions[name][0] for name in ("train", "val", "test")],
        [predictions[name][1] for name in ("train", "val", "test")],
    )
    for index, split in enumerate(("train", "val", "test")):
        plot_parity(
            split,
            predictions[split][0],
            predictions[split][1],
            metrics[split],
            limits,
            COLORS[index],
        )
    plot_combined_parity(predictions, metrics, limits)
    save_combined_dat(predictions)
    save_and_plot_history(checkpoint["history"], int(checkpoint["best_epoch"]))
    save_metrics_table(metrics, predictions)
    eval_seconds = time.perf_counter() - eval_start
    logger.info(
        "Evaluation and plotting time: %.2f s (%s)",
        eval_seconds,
        format_duration(eval_seconds),
    )

    model_config = checkpoint["model_config"]
    total_parameters = sum(parameter.numel() for parameter in model.parameters())
    trainable_parameters = sum(
        parameter.numel() for parameter in model.parameters() if parameter.requires_grad
    )
    logger.info("Model configuration: %s", model_config)
    logger.info("Graph configuration: %s", checkpoint["graph_config"])
    logger.info("Model structure:\n%s", model)
    logger.info("Total parameters: %d", total_parameters)
    logger.info("Trainable parameters: %d", trainable_parameters)
    logger.info("Checkpoint: %s", CHECKPOINT_PATH)
    total_seconds = time.perf_counter() - total_start
    logger.info(
        "Total runtime: %.2f s (%s)", total_seconds, format_duration(total_seconds)
    )


if __name__ == "__main__":
    main()

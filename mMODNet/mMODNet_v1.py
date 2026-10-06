from __future__ import annotations

import json
import logging
import math
import os
import random
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from sklearn.metrics import (
    mean_absolute_error,
    mean_squared_error,
    normalized_mutual_info_score,
    r2_score,
)
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, TensorDataset


# ==================== User configuration ====================

MODEL_NAME = "MODNet"
RUN_VERSION = "v1"

EXCEL_PATH = "./data.xlsx"
CIF_DIR = "./cif"

TARGET_NAME = "bandgap"
TARGET_UNIT = "eV"

SEED = 42
TRAIN_RATIO = 0.80
VAL_RATIO = 0.10
TEST_RATIO = 0.10

N_FEATURES = 128
MAX_RR_CANDIDATES = 256
NMI_BINS = 16
CONTINUOUS_ONLY = True
CACHE_DESCRIPTORS = True
REBUILD_DESCRIPTOR_CACHE = False

COMMON_DIMS = (256, 128)
INTERMEDIATE_DIMS = (128, 64)
PROPERTY_DIMS = (64, 32)
DROPOUT = 0.10

BATCH_SIZE = 128
MAX_EPOCHS = 300
LEARNING_RATE = 1.0e-3
WEIGHT_DECAY = 1.0e-5
PATIENCE = 60
MIN_DELTA = 1.0e-5
LR_FACTOR = 0.5
LR_PATIENCE = 15
MIN_LR = 1.0e-6
NUM_WORKERS = 0
USE_AMP = True
AMP_DTYPE = "bfloat16"

USE_OPTUNA = False
OPTUNA_N_TRIALS = 30
OPTUNA_TIMEOUT = None
OPTUNA_MAX_EPOCHS = 250
OPTUNA_PATIENCE = 35

FIG_DPI = 600
TITLE_FONTSIZE = 14
LABEL_FONTSIZE = 12
TICK_FONTSIZE = 10
LEGEND_FONTSIZE = 10
ANNOTATION_FONTSIZE = 10
USE_GRID = True

COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:red"]

# ============================================================


OUTPUT_DIR = Path(f"{MODEL_NAME}_{RUN_VERSION}")
FIGURE_DIR = OUTPUT_DIR / "figure"
DAT_DIR = OUTPUT_DIR / "dat"
TABLE_DIR = OUTPUT_DIR / "table"
LOG_DIR = OUTPUT_DIR / "log"
SPLIT_DIR = OUTPUT_DIR / "split"
CACHE_DIR = OUTPUT_DIR / "cache"
CHECKPOINT_PATH = OUTPUT_DIR / f"{MODEL_NAME}_best.pt"
SPLIT_PATH = SPLIT_DIR / f"{MODEL_NAME}_split.csv"
DESCRIPTOR_CACHE_PATH = CACHE_DIR / f"{MODEL_NAME}_descriptors.pkl"
LOG_PATH = LOG_DIR / f"{MODEL_NAME}_training.log"


@dataclass
class Record:
    cif_name: str
    cif_path: str
    target: float


@dataclass
class TrainConfig:
    n_features: int
    common_dims: Tuple[int, ...]
    intermediate_dims: Tuple[int, ...]
    property_dims: Tuple[int, ...]
    dropout: float
    batch_size: int
    learning_rate: float
    weight_decay: float
    max_epochs: int
    patience: int
    min_delta: float


def ensure_directories() -> None:
    for path in (
        OUTPUT_DIR,
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
    stream_handler = logging.StreamHandler(sys.stdout)
    file_handler.setFormatter(formatter)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def set_global_seed(seed: int) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    try:
        torch.use_deterministic_algorithms(True, warn_only=True)
    except TypeError:
        torch.use_deterministic_algorithms(True)


def format_duration(seconds: float) -> str:
    seconds_int = max(0, int(round(seconds)))
    hours, rem = divmod(seconds_int, 3600)
    minutes, secs = divmod(rem, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def get_device(logger: logging.Logger) -> torch.device:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info("Device: %s", device)
    logger.info("PyTorch version: %s", torch.__version__)
    if device.type == "cuda":
        props = torch.cuda.get_device_properties(0)
        logger.info("GPU: %s", torch.cuda.get_device_name(0))
        logger.info("GPU memory: %.3f GiB", props.total_memory / 1024**3)
        logger.info("CUDA version: %s", torch.version.cuda)
    return device


def load_excel_records(logger: logging.Logger) -> List[Record]:
    excel_path = Path(EXCEL_PATH)
    cif_dir = Path(CIF_DIR)
    if not excel_path.is_file():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    if not cif_dir.is_dir():
        raise FileNotFoundError(f"CIF directory not found: {cif_dir}")

    frame = pd.read_excel(excel_path)
    if frame.shape[1] < 2:
        raise ValueError("The Excel file must contain at least two columns.")

    records: List[Record] = []
    for row_idx, row in frame.iloc[:, :2].iterrows():
        raw_name, raw_target = row.iloc[0], row.iloc[1]
        if pd.isna(raw_name):
            logger.warning("Skipped Excel row %d: empty CIF name.", row_idx + 2)
            continue
        cif_name = f"{int(raw_name)}.cif"
        if not cif_name.lower().endswith(".cif"):
            cif_name = f"{cif_name}.cif"
        try:
            target = float(raw_target)
        except (TypeError, ValueError):
            logger.warning("Skipped %s: target is not numeric.", cif_name)
            continue
        if not np.isfinite(target):
            logger.warning("Skipped %s: target is not finite.", cif_name)
            continue
        cif_path = cif_dir / cif_name
        if not cif_path.is_file():
            logger.warning("Skipped %s: CIF file not found.", cif_name)
            continue
        records.append(Record(cif_name, str(cif_path), target))

    if len(records) < 20:
        raise ValueError(
            f"Only {len(records)} valid rows were found. At least 20 are required."
        )
    if len({record.cif_name for record in records}) != len(records):
        raise ValueError("Duplicate CIF names were found in the Excel file.")
    return records


def _to_numeric(value: Any) -> float:
    if value is None:
        return np.nan
    if isinstance(value, (bool, np.bool_)):
        return float(value)
    if isinstance(value, str):
        crystal_system = {
            "cubic": 1.0,
            "tetragonal": 2.0,
            "orthorhombic": 3.0,
            "orthorombic": 3.0,
            "hexagonal": 4.0,
            "trigonal": 5.0,
            "monoclinic": 6.0,
            "triclinic": 7.0,
        }
        return crystal_system.get(value.strip().lower(), np.nan)
    try:
        out = float(value)
    except (TypeError, ValueError, OverflowError):
        return np.nan
    return out if np.isfinite(out) else np.nan


def build_featurizers() -> List[Tuple[str, str, Any]]:
    try:
        from matminer.featurizers.composition import (
            BandCenter,
            ElementFraction,
            ElementProperty,
            Stoichiometry,
            TMetalFraction,
            ValenceOrbital,
        )
        from matminer.featurizers.structure import (
            DensityFeatures,
            EwaldEnergy,
            GlobalSymmetryFeatures,
            StructuralComplexity,
        )
        from matminer.utils.data import DemlData, PymatgenData
    except ImportError as exc:
        raise ImportError(
            "matminer is required. Install it with: pip install matminer"
        ) from exc

    magpie = ElementProperty.from_preset("magpie")
    magpie.stats = ["mean", "avg_dev"]
    pymatgen_features = [
        "block",
        "mendeleev_no",
        "electrical_resistivity",
        "velocity_of_sound",
        "thermal_conductivity",
        "bulk_modulus",
        "coefficient_of_linear_thermal_expansion",
    ]
    deml_features = [
        "atom_radius",
        "molar_vol",
        "heat_fusion",
        "boiling_point",
        "heat_cap",
        "first_ioniz",
        "electric_pol",
        "GGAU_Etot",
        "mus_fere",
        "FERE correction",
    ]
    return [
        ("composition", "BandCenter", BandCenter()),
        ("composition", "ElementFraction", ElementFraction()),
        ("composition", "Magpie", magpie),
        (
            "composition",
            "PymatgenData",
            ElementProperty(
                data_source=PymatgenData(),
                stats=["mean", "avg_dev"],
                features=pymatgen_features,
            ),
        ),
        (
            "composition",
            "DemlData",
            ElementProperty(
                data_source=DemlData(),
                stats=["mean", "avg_dev"],
                features=deml_features,
            ),
        ),
        ("composition", "Stoichiometry", Stoichiometry(p_list=[2, 3, 5, 7, 10])),
        ("composition", "TMetalFraction", TMetalFraction()),
        ("composition", "ValenceOrbital", ValenceOrbital(props=["frac"])),
        ("structure", "DensityFeatures", DensityFeatures()),
        ("structure", "EwaldEnergy", EwaldEnergy()),
        ("structure", "GlobalSymmetryFeatures", GlobalSymmetryFeatures()),
        ("structure", "StructuralComplexity", StructuralComplexity()),
    ]


def _feature_schema(
    featurizers: Sequence[Tuple[str, str, Any]],
) -> List[Tuple[str, str, Any, List[str]]]:
    schema = []
    for input_kind, prefix, featurizer in featurizers:
        labels = [f"{prefix}|{label}" for label in featurizer.feature_labels()]
        schema.append((input_kind, prefix, featurizer, labels))
    return schema


def featurize_structure(
    cif_path: str,
    schema: Sequence[Tuple[str, str, Any, List[str]]],
) -> Tuple[np.ndarray, List[str]]:
    from pymatgen.core import Structure

    structure = Structure.from_file(cif_path)
    composition = structure.composition
    values: List[float] = []
    names: List[str] = []
    for input_kind, _, featurizer, labels in schema:
        source = composition if input_kind == "composition" else structure
        try:
            raw_values = list(featurizer.featurize(source))
        except Exception:
            raw_values = [np.nan] * len(labels)
        if len(raw_values) != len(labels):
            raw_values = [np.nan] * len(labels)
        values.extend(_to_numeric(value) for value in raw_values)
        names.extend(labels)
    return np.asarray(values, dtype=np.float64), names


def _drop_discontinuous_features(frame: pd.DataFrame) -> pd.DataFrame:
    if not CONTINUOUS_ONLY:
        return frame
    suffixes = (
        "mean electric_pol",
        "mean FERE correction",
        "mean GGAU_Etot",
        "mean heat_fusion",
        "mean mus_fere",
    )
    columns = [
        column
        for column in frame.columns
        if not any(column.endswith(x) for x in suffixes)
    ]
    return frame.loc[:, columns]


def compute_descriptors(
    records: Sequence[Record], logger: logging.Logger
) -> Tuple[pd.DataFrame, np.ndarray, List[Record]]:
    expected_names = [record.cif_name for record in records]
    expected_targets = np.asarray(
        [record.target for record in records], dtype=np.float64
    )
    if (
        CACHE_DESCRIPTORS
        and DESCRIPTOR_CACHE_PATH.is_file()
        and not REBUILD_DESCRIPTOR_CACHE
    ):
        try:
            cached = pd.read_pickle(DESCRIPTOR_CACHE_PATH)
            if cached.get("cif_names") == expected_names and np.allclose(
                cached.get("targets"), expected_targets
            ):
                frame = cached["features"]
                valid_names = cached["valid_names"]
                record_map = {record.cif_name: record for record in records}
                valid_records = [record_map[name] for name in valid_names]
                targets = np.asarray([record.target for record in valid_records])
                logger.info("Loaded descriptor cache: %s", DESCRIPTOR_CACHE_PATH)
                return frame, targets, valid_records
        except Exception as exc:
            logger.warning("Descriptor cache was ignored: %s", exc)

    featurizers = build_featurizers()
    schema = _feature_schema(featurizers)
    rows: List[np.ndarray] = []
    valid_records: List[Record] = []
    feature_names: Optional[List[str]] = None
    for idx, record in enumerate(records, start=1):
        try:
            values, names = featurize_structure(record.cif_path, schema)
            if feature_names is None:
                feature_names = names
            elif names != feature_names:
                raise RuntimeError("Descriptor schema changed during featurization.")
            rows.append(values)
            valid_records.append(record)
        except Exception as exc:
            logger.warning(
                "Skipped %s: descriptor generation failed: %s", record.cif_name, exc
            )
        if idx == 1 or idx % 25 == 0 or idx == len(records):
            logger.info("Featurized %d/%d structures.", idx, len(records))

    if feature_names is None or len(valid_records) < 20:
        raise ValueError("Fewer than 20 structures produced valid descriptors.")
    frame = pd.DataFrame(
        np.vstack(rows),
        index=[r.cif_name for r in valid_records],
        columns=feature_names,
    )
    frame = _drop_discontinuous_features(frame)
    frame = frame.replace([np.inf, -np.inf], np.nan)
    targets = np.asarray([record.target for record in valid_records], dtype=np.float64)
    if CACHE_DESCRIPTORS:
        pd.to_pickle(
            {
                "cif_names": expected_names,
                "targets": expected_targets,
                "valid_names": [r.cif_name for r in valid_records],
                "features": frame,
                "featurizer": "Matminer2024Fast-compatible",
                "continuous_only": CONTINUOUS_ONLY,
            },
            DESCRIPTOR_CACHE_PATH,
        )
        logger.info("Saved descriptor cache: %s", DESCRIPTOR_CACHE_PATH)
    return frame, targets, valid_records


def create_or_load_split(
    records: Sequence[Record], logger: logging.Logger
) -> Dict[str, np.ndarray]:
    names = [record.cif_name for record in records]
    target_map = {record.cif_name: record.target for record in records}
    if SPLIT_PATH.is_file():
        split_frame = pd.read_csv(SPLIT_PATH)
        required = {"cif_name", "target", "split"}
        if required.issubset(split_frame.columns) and set(
            split_frame["cif_name"]
        ) == set(names):
            if set(split_frame["split"]).issubset({"train", "val", "test"}):
                logger.info("Reused split file: %s", SPLIT_PATH)
                return {
                    split: np.asarray(
                        [
                            names.index(name)
                            for name in split_frame.loc[
                                split_frame["split"] == split, "cif_name"
                            ]
                        ],
                        dtype=np.int64,
                    )
                    for split in ("train", "val", "test")
                }
        logger.warning("Existing split file is incompatible and will be replaced.")

    indices = np.arange(len(records))
    train_idx, temp_idx = train_test_split(
        indices,
        test_size=1.0 - TRAIN_RATIO,
        random_state=SEED,
        shuffle=True,
    )
    relative_test = TEST_RATIO / (VAL_RATIO + TEST_RATIO)
    val_idx, test_idx = train_test_split(
        temp_idx,
        test_size=relative_test,
        random_state=SEED,
        shuffle=True,
    )
    split_map = {"train": train_idx, "val": val_idx, "test": test_idx}
    rows = []
    for split in ("train", "val", "test"):
        for idx in split_map[split]:
            name = names[int(idx)]
            rows.append({"cif_name": name, "target": target_map[name], "split": split})
    pd.DataFrame(rows).to_csv(SPLIT_PATH, index=False)
    logger.info("Created split file: %s", SPLIT_PATH)
    return {key: np.asarray(value, dtype=np.int64) for key, value in split_map.items()}


def _quantile_bins(values: np.ndarray, n_bins: int) -> np.ndarray:
    finite = np.isfinite(values)
    if finite.sum() == 0:
        return np.zeros(len(values), dtype=np.int32)
    clean = values.copy()
    median = np.nanmedian(clean[finite])
    clean[~finite] = median
    unique_count = np.unique(clean).size
    if unique_count <= 1:
        return np.zeros(len(clean), dtype=np.int32)
    bins = min(n_bins, unique_count)
    try:
        codes = pd.qcut(clean, q=bins, labels=False, duplicates="drop")
        return np.asarray(codes, dtype=np.int32)
    except ValueError:
        return np.zeros(len(clean), dtype=np.int32)


def rank_relevance_redundancy(
    x_train: np.ndarray,
    y_train: np.ndarray,
    feature_names: Sequence[str],
    n_features: int,
    max_candidates: int,
    n_bins: int,
    logger: logging.Logger,
) -> Tuple[List[int], List[Dict[str, float]]]:
    y_disc = _quantile_bins(y_train, n_bins)
    x_disc = np.column_stack(
        [_quantile_bins(x_train[:, idx], n_bins) for idx in range(x_train.shape[1])]
    )
    relevance = np.asarray(
        [
            normalized_mutual_info_score(
                y_disc, x_disc[:, idx], average_method="arithmetic"
            )
            for idx in range(x_disc.shape[1])
        ],
        dtype=np.float64,
    )
    relevance = np.nan_to_num(relevance, nan=0.0, posinf=0.0, neginf=0.0)
    candidate_count = min(max_candidates, len(feature_names))
    candidates = np.argsort(-relevance, kind="stable")[:candidate_count].tolist()
    n_select = min(n_features, candidate_count)
    if n_select < 1:
        raise ValueError("No descriptor is available for feature selection.")

    selected = [candidates.pop(0)]
    ranking = [
        {
            "rank": 1,
            "feature_index": int(selected[0]),
            "feature": feature_names[selected[0]],
            "relevance_nmi": float(relevance[selected[0]]),
            "rr_score": float("nan"),
        }
    ]
    redundancy_cache: Dict[Tuple[int, int], float] = {}
    while len(selected) < n_select:
        n_chosen = len(selected)
        p_value = max(0.1, 4.5 - n_chosen**0.4)
        c_value = 1.0e-6 * n_chosen**3
        best_idx = None
        best_score = -np.inf
        for candidate in candidates:
            max_redundancy = 0.0
            for chosen in selected:
                key = (min(candidate, chosen), max(candidate, chosen))
                if key not in redundancy_cache:
                    redundancy_cache[key] = normalized_mutual_info_score(
                        x_disc[:, candidate],
                        x_disc[:, chosen],
                        average_method="arithmetic",
                    )
                max_redundancy = max(max_redundancy, redundancy_cache[key])
            score = relevance[candidate] / (max_redundancy**p_value + c_value)
            if score > best_score:
                best_score = score
                best_idx = candidate
        if best_idx is None:
            break
        selected.append(best_idx)
        candidates.remove(best_idx)
        ranking.append(
            {
                "rank": len(selected),
                "feature_index": int(best_idx),
                "feature": feature_names[best_idx],
                "relevance_nmi": float(relevance[best_idx]),
                "rr_score": float(best_score),
            }
        )
        if len(selected) % 25 == 0 or len(selected) == n_select:
            logger.info(
                "Selected %d/%d descriptors by RR ranking.", len(selected), n_select
            )
    return selected, ranking


def fit_preprocessor(
    features: pd.DataFrame,
    targets: np.ndarray,
    train_idx: np.ndarray,
    logger: logging.Logger,
) -> Tuple[np.ndarray, Dict[str, Any]]:
    train_frame = features.iloc[train_idx]
    valid_columns = [
        column for column in features.columns if train_frame[column].notna().any()
    ]
    if not valid_columns:
        raise ValueError("All descriptors are missing in the training set.")
    work = features.loc[:, valid_columns].to_numpy(dtype=np.float64)
    medians = np.nanmedian(work[train_idx], axis=0)
    medians = np.where(np.isfinite(medians), medians, 0.0)
    work = np.where(np.isfinite(work), work, medians[None, :])
    variances = np.var(work[train_idx], axis=0)
    keep = variances > 1.0e-14
    work = work[:, keep]
    medians = medians[keep]
    kept_names = [name for name, flag in zip(valid_columns, keep) if flag]
    if work.shape[1] == 0:
        raise ValueError("All training descriptors have zero variance.")

    selected_local, ranking = rank_relevance_redundancy(
        work[train_idx],
        targets[train_idx],
        kept_names,
        N_FEATURES,
        MAX_RR_CANDIDATES,
        NMI_BINS,
        logger,
    )
    work = work[:, selected_local]
    selected_names = [kept_names[idx] for idx in selected_local]
    means = work[train_idx].mean(axis=0)
    scales = work[train_idx].std(axis=0)
    scales = np.where(scales > 1.0e-12, scales, 1.0)
    transformed = (work - means[None, :]) / scales[None, :]
    y_mean = float(targets[train_idx].mean())
    y_scale = float(targets[train_idx].std())
    if y_scale <= 1.0e-12:
        y_scale = 1.0

    pd.DataFrame(ranking).to_csv(
        TABLE_DIR / f"{MODEL_NAME}_selected_features.dat", sep="\t", index=False
    )
    preprocessing = {
        "all_feature_names": list(features.columns),
        "valid_feature_names": valid_columns,
        "nonconstant_mask": keep.astype(bool).tolist(),
        "imputer_medians": medians.tolist(),
        "selected_local_indices": [int(idx) for idx in selected_local],
        "selected_feature_names": selected_names,
        "feature_means": means.tolist(),
        "feature_scales": scales.tolist(),
        "target_mean": y_mean,
        "target_scale": y_scale,
        "rr_ranking": ranking,
    }
    return transformed.astype(np.float32), preprocessing


def transform_features(
    features: pd.DataFrame,
    preprocessing: Dict[str, Any],
) -> np.ndarray:
    all_names = preprocessing["all_feature_names"]
    aligned = features.reindex(columns=all_names)
    valid_names = preprocessing["valid_feature_names"]
    work = aligned.loc[:, valid_names].to_numpy(dtype=np.float64)
    nonconstant = np.asarray(preprocessing["nonconstant_mask"], dtype=bool)
    medians = np.asarray(preprocessing["imputer_medians"], dtype=np.float64)
    work = work[:, nonconstant]
    work = np.where(np.isfinite(work), work, medians[None, :])
    selected = np.asarray(preprocessing["selected_local_indices"], dtype=np.int64)
    work = work[:, selected]
    means = np.asarray(preprocessing["feature_means"], dtype=np.float64)
    scales = np.asarray(preprocessing["feature_scales"], dtype=np.float64)
    return ((work - means[None, :]) / scales[None, :]).astype(np.float32)


class DenseBlock(nn.Module):
    def __init__(self, input_dim: int, dims: Sequence[int], dropout: float):
        super().__init__()
        layers: List[nn.Module] = []
        current = input_dim
        for dim in dims:
            layers.append(nn.Linear(current, dim))
            layers.append(nn.ReLU())
            if dropout > 0:
                layers.append(nn.Dropout(dropout))
            current = dim
        self.layers = nn.Sequential(*layers)
        self.output_dim = current

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.layers(x)


class MODNetRegressor(nn.Module):
    def __init__(
        self,
        n_features: int,
        common_dims: Sequence[int],
        intermediate_dims: Sequence[int],
        property_dims: Sequence[int],
        dropout: float,
    ):
        super().__init__()
        self.common = DenseBlock(n_features, common_dims, dropout)
        self.intermediate = DenseBlock(
            self.common.output_dim, intermediate_dims, dropout
        )
        self.property_block = DenseBlock(
            self.intermediate.output_dim, property_dims, dropout
        )
        self.output = nn.Linear(self.property_block.output_dim, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.common(x)
        x = self.intermediate(x)
        x = self.property_block(x)
        return self.output(x).view(-1)


def make_loader(
    x: np.ndarray,
    y: np.ndarray,
    indices: np.ndarray,
    batch_size: int,
    shuffle: bool,
) -> DataLoader:
    dataset = TensorDataset(
        torch.from_numpy(x[indices]),
        torch.from_numpy(y[indices].astype(np.float32)),
    )
    generator = torch.Generator()
    generator.manual_seed(SEED)
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=NUM_WORKERS,
        pin_memory=torch.cuda.is_available(),
        generator=generator if shuffle else None,
        persistent_workers=NUM_WORKERS > 0,
    )


def _autocast_dtype() -> torch.dtype:
    return torch.bfloat16 if AMP_DTYPE.lower() == "bfloat16" else torch.float16


def predict_loader(
    model: nn.Module,
    loader: DataLoader,
    device: torch.device,
    target_mean: float,
    target_scale: float,
) -> Tuple[np.ndarray, np.ndarray]:
    model.eval()
    predictions: List[np.ndarray] = []
    targets: List[np.ndarray] = []
    amp_enabled = USE_AMP and device.type == "cuda"
    with torch.inference_mode():
        for batch_x, batch_y in loader:
            batch_x = batch_x.to(device, non_blocking=True)
            with torch.autocast(
                device_type=device.type,
                dtype=_autocast_dtype(),
                enabled=amp_enabled,
            ):
                pred_scaled = model(batch_x)
            pred = pred_scaled.float().cpu().numpy() * target_scale + target_mean
            predictions.append(pred)
            targets.append(batch_y.numpy())
    return np.concatenate(targets), np.concatenate(predictions)


def train_one_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    criterion: nn.Module,
    device: torch.device,
    target_mean: float,
    target_scale: float,
) -> float:
    model.train()
    amp_enabled = USE_AMP and device.type == "cuda"
    scaler_enabled = amp_enabled and _autocast_dtype() == torch.float16
    scaler = torch.amp.GradScaler("cuda", enabled=scaler_enabled)
    total_loss = 0.0
    total_count = 0
    for batch_x, batch_y in loader:
        batch_x = batch_x.to(device, non_blocking=True)
        batch_y = batch_y.to(device, non_blocking=True)
        scaled_y = (batch_y - target_mean) / target_scale
        optimizer.zero_grad(set_to_none=True)
        with torch.autocast(
            device_type=device.type,
            dtype=_autocast_dtype(),
            enabled=amp_enabled,
        ):
            pred_scaled = model(batch_x)
            loss = criterion(pred_scaled.view(-1), scaled_y.view(-1))
        if not torch.isfinite(loss):
            raise FloatingPointError("A non-finite training loss was encountered.")
        scaler.scale(loss).backward()
        scaler.step(optimizer)
        scaler.update()
        total_loss += float(loss.detach().cpu()) * len(batch_y)
        total_count += len(batch_y)
    return total_loss / max(total_count, 1)


def calculate_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    return {
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(math.sqrt(mean_squared_error(y_true, y_pred))),
        "r2": float(r2_score(y_true, y_pred)),
    }


def train_model(
    x: np.ndarray,
    y: np.ndarray,
    splits: Dict[str, np.ndarray],
    config: TrainConfig,
    preprocessing: Dict[str, Any],
    device: torch.device,
    logger: logging.Logger,
    checkpoint_path: Optional[Path],
    verbose: bool = True,
) -> Tuple[nn.Module, List[Dict[str, float]], int, float, bool]:
    x_use = x[:, : config.n_features]
    loaders = {
        split: make_loader(
            x_use,
            y,
            indices,
            config.batch_size,
            shuffle=(split == "train"),
        )
        for split, indices in splits.items()
    }
    model = MODNetRegressor(
        config.n_features,
        config.common_dims,
        config.intermediate_dims,
        config.property_dims,
        config.dropout,
    ).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=config.learning_rate,
        weight_decay=config.weight_decay,
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="min",
        factor=LR_FACTOR,
        patience=LR_PATIENCE,
        min_lr=MIN_LR,
    )
    criterion = nn.MSELoss()
    target_mean = float(preprocessing["target_mean"])
    target_scale = float(preprocessing["target_scale"])
    best_val_rmse = float("inf")
    best_epoch = 0
    epochs_without_improvement = 0
    history: List[Dict[str, float]] = []
    best_state: Optional[Dict[str, torch.Tensor]] = None
    early_stopped = False

    for epoch in range(1, config.max_epochs + 1):
        train_mse_scaled = train_one_epoch(
            model,
            loaders["train"],
            optimizer,
            criterion,
            device,
            target_mean,
            target_scale,
        )
        y_train, p_train = predict_loader(
            model, loaders["train"], device, target_mean, target_scale
        )
        y_val, p_val = predict_loader(
            model, loaders["val"], device, target_mean, target_scale
        )
        train_rmse = math.sqrt(mean_squared_error(y_train, p_train))
        val_rmse = math.sqrt(mean_squared_error(y_val, p_val))
        lr = float(optimizer.param_groups[0]["lr"])
        history.append(
            {
                "epoch": epoch,
                "train_rmse_eV": train_rmse,
                "val_rmse_eV": val_rmse,
                "learning_rate": lr,
                "train_mse_loss_scaled": train_mse_scaled,
            }
        )
        scheduler.step(val_rmse)
        if verbose:
            logger.info(
                "Epoch %04d | train RMSE %.6f eV | val RMSE %.6f eV | lr %.6e",
                epoch,
                train_rmse,
                val_rmse,
                lr,
            )

        if val_rmse < best_val_rmse - config.min_delta:
            best_val_rmse = val_rmse
            best_epoch = epoch
            epochs_without_improvement = 0
            best_state = {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            }
            if checkpoint_path is not None:
                checkpoint = {
                    "model_name": MODEL_NAME,
                    "run_version": RUN_VERSION,
                    "model_state_dict": best_state,
                    "model_config": asdict(config),
                    "best_epoch": best_epoch,
                    "best_val_rmse": best_val_rmse,
                    "seed": SEED,
                    "target_name": TARGET_NAME,
                    "target_unit": TARGET_UNIT,
                    "preprocessing": preprocessing,
                    "featurizer_config": {
                        "name": "Matminer2024Fast-compatible",
                        "continuous_only": CONTINUOUS_ONLY,
                    },
                }
                torch.save(checkpoint, checkpoint_path)
        else:
            epochs_without_improvement += 1
            if epochs_without_improvement >= config.patience:
                early_stopped = True
                break

    if best_state is None:
        raise RuntimeError("Training did not produce a valid checkpoint.")
    model.load_state_dict(best_state)
    return model, history, best_epoch, best_val_rmse, early_stopped


def default_train_config(n_available_features: int) -> TrainConfig:
    return TrainConfig(
        n_features=min(N_FEATURES, n_available_features),
        common_dims=tuple(COMMON_DIMS),
        intermediate_dims=tuple(INTERMEDIATE_DIMS),
        property_dims=tuple(PROPERTY_DIMS),
        dropout=DROPOUT,
        batch_size=BATCH_SIZE,
        learning_rate=LEARNING_RATE,
        weight_decay=WEIGHT_DECAY,
        max_epochs=MAX_EPOCHS,
        patience=PATIENCE,
        min_delta=MIN_DELTA,
    )


def run_optuna(
    x: np.ndarray,
    y: np.ndarray,
    splits: Dict[str, np.ndarray],
    preprocessing: Dict[str, Any],
    device: torch.device,
    logger: logging.Logger,
) -> TrainConfig:
    try:
        import optuna
    except ImportError as exc:
        raise ImportError("Optuna is enabled but is not installed.") from exc

    feature_options = sorted(
        set(
            min(value, x.shape[1])
            for value in (32, 64, 96, 128, 192)
            if min(value, x.shape[1]) >= 8
        )
    )

    def objective(trial: Any) -> float:
        width = trial.suggest_categorical("width", [64, 128, 256, 384])
        depth = trial.suggest_int("depth", 2, 4)
        hidden = [max(16, width // (2**idx)) for idx in range(depth)]
        if depth == 2:
            common_dims, intermediate_dims, property_dims = (
                (hidden[0],),
                (hidden[1],),
                (max(8, hidden[1] // 2),),
            )
        elif depth == 3:
            common_dims, intermediate_dims, property_dims = (
                (hidden[0],),
                (hidden[1],),
                (hidden[2],),
            )
        else:
            common_dims, intermediate_dims, property_dims = (
                (hidden[0], hidden[1]),
                (hidden[2],),
                (hidden[3],),
            )
        config = TrainConfig(
            n_features=trial.suggest_categorical("n_features", feature_options),
            common_dims=tuple(common_dims),
            intermediate_dims=tuple(intermediate_dims),
            property_dims=tuple(property_dims),
            dropout=trial.suggest_float("dropout", 0.0, 0.35),
            batch_size=trial.suggest_categorical("batch_size", [64, 128, 256]),
            learning_rate=trial.suggest_float(
                "learning_rate", 2.0e-4, 3.0e-3, log=True
            ),
            weight_decay=trial.suggest_float("weight_decay", 1.0e-7, 1.0e-3, log=True),
            max_epochs=OPTUNA_MAX_EPOCHS,
            patience=OPTUNA_PATIENCE,
            min_delta=MIN_DELTA,
        )
        set_global_seed(SEED)
        _, _, _, best_rmse, _ = train_model(
            x,
            y,
            splits,
            config,
            preprocessing,
            device,
            logger,
            checkpoint_path=None,
            verbose=False,
        )
        return best_rmse

    study = optuna.create_study(
        direction="minimize", sampler=optuna.samplers.TPESampler(seed=SEED)
    )
    study.optimize(objective, n_trials=OPTUNA_N_TRIALS, timeout=OPTUNA_TIMEOUT)
    result = {"best_value": study.best_value, "best_params": study.best_params}
    with open(
        TABLE_DIR / f"{MODEL_NAME}_optuna_best_params.json", "w", encoding="utf-8"
    ) as handle:
        json.dump(result, handle, indent=2)
    logger.info("Optuna best validation RMSE: %.6f eV", study.best_value)
    logger.info("Optuna best parameters: %s", study.best_params)

    params = study.best_params
    width = int(params["width"])
    depth = int(params["depth"])
    hidden = [max(16, width // (2**idx)) for idx in range(depth)]
    if depth == 2:
        common_dims, intermediate_dims, property_dims = (
            (hidden[0],),
            (hidden[1],),
            (max(8, hidden[1] // 2),),
        )
    elif depth == 3:
        common_dims, intermediate_dims, property_dims = (
            (hidden[0],),
            (hidden[1],),
            (hidden[2],),
        )
    else:
        common_dims, intermediate_dims, property_dims = (
            (hidden[0], hidden[1]),
            (hidden[2],),
            (hidden[3],),
        )
    return TrainConfig(
        n_features=int(params["n_features"]),
        common_dims=tuple(common_dims),
        intermediate_dims=tuple(intermediate_dims),
        property_dims=tuple(property_dims),
        dropout=float(params["dropout"]),
        batch_size=int(params["batch_size"]),
        learning_rate=float(params["learning_rate"]),
        weight_decay=float(params["weight_decay"]),
        max_epochs=MAX_EPOCHS,
        patience=PATIENCE,
        min_delta=MIN_DELTA,
    )


def load_trained_model(
    checkpoint_path: str | Path,
    device: Optional[torch.device] = None,
) -> Tuple[MODNetRegressor, Dict[str, Any]]:
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint = torch.load(checkpoint_path, map_location=device, weights_only=False)
    config = checkpoint["model_config"]
    model = MODNetRegressor(
        n_features=int(config["n_features"]),
        common_dims=tuple(config["common_dims"]),
        intermediate_dims=tuple(config["intermediate_dims"]),
        property_dims=tuple(config["property_dims"]),
        dropout=float(config["dropout"]),
    ).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    return model, checkpoint


def descriptors_for_cifs(cif_paths: Sequence[str | Path]) -> pd.DataFrame:
    schema = _feature_schema(build_featurizers())
    rows = []
    names = None
    index = []
    for path in cif_paths:
        values, current_names = featurize_structure(str(path), schema)
        if names is None:
            names = current_names
        elif current_names != names:
            raise RuntimeError("Descriptor schema changed during prediction.")
        rows.append(values)
        index.append(Path(path).name)
    if names is None:
        return pd.DataFrame()
    return _drop_discontinuous_features(pd.DataFrame(rows, index=index, columns=names))


def predict_cifs(
    checkpoint_path: str | Path,
    cif_paths: Sequence[str | Path],
    device: Optional[torch.device] = None,
) -> pd.DataFrame:
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, checkpoint = load_trained_model(checkpoint_path, device)
    frame = descriptors_for_cifs(cif_paths)
    x = transform_features(frame, checkpoint["preprocessing"])
    n_features = int(checkpoint["model_config"]["n_features"])
    tensor = torch.from_numpy(x[:, :n_features]).to(device)
    with torch.inference_mode():
        with torch.autocast(
            device_type=device.type,
            dtype=_autocast_dtype(),
            enabled=USE_AMP and device.type == "cuda",
        ):
            scaled = model(tensor)
    pred = scaled.float().cpu().numpy()
    pred = pred * float(checkpoint["preprocessing"]["target_scale"]) + float(
        checkpoint["preprocessing"]["target_mean"]
    )
    return pd.DataFrame(
        {"cif_name": frame.index, f"predicted_{TARGET_NAME}_{TARGET_UNIT}": pred}
    )


def save_prediction_dat(
    path: Path,
    split: str,
    names: Sequence[str],
    y_true: np.ndarray,
    y_pred: np.ndarray,
    include_split: bool = False,
) -> None:
    data: Dict[str, Any] = {}
    if include_split:
        data["split"] = [split] * len(y_true)
    data.update(
        {
            "cif_name": list(names),
            "true_bandgap_eV": y_true,
            "predicted_bandgap_eV": y_pred,
            "error_eV": y_pred - y_true,
            "absolute_error_eV": np.abs(y_pred - y_true),
        }
    )
    pd.DataFrame(data).to_csv(path, sep="\t", index=False, float_format="%.8f")


def _plot_limits(all_true: np.ndarray, all_pred: np.ndarray) -> Tuple[float, float]:
    low = float(min(all_true.min(), all_pred.min()))
    high = float(max(all_true.max(), all_pred.max()))
    margin = max(0.05 * (high - low), 0.05)
    return low - margin, high + margin


def plot_parity(
    path: Path,
    split: str,
    y_true: np.ndarray,
    y_pred: np.ndarray,
    color: str,
    limits: Tuple[float, float],
) -> None:
    metrics = calculate_metrics(y_true, y_pred)
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    ax.scatter(y_true, y_pred, s=22, alpha=0.78, color=color, edgecolors="none")
    ax.plot(limits, limits, "--", color="black", linewidth=1.2)
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"{MODEL_NAME}: {split.capitalize()} Set", fontsize=TITLE_FONTSIZE)
    ax.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.tick_params(labelsize=TICK_FONTSIZE)
    text = f"MAE = {metrics['mae']:.4f} eV\nRMSE = {metrics['rmse']:.4f} eV\n$R^2$ = {metrics['r2']:.4f}"
    ax.text(
        0.04,
        0.96,
        text,
        transform=ax.transAxes,
        va="top",
        fontsize=ANNOTATION_FONTSIZE,
        bbox={"facecolor": "white", "alpha": 0.85, "edgecolor": "0.8"},
    )
    if USE_GRID:
        ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)
    fig.tight_layout()
    fig.savefig(path, dpi=FIG_DPI, format="jpg", bbox_inches="tight")
    plt.close(fig)


def plot_combined_parity(
    path: Path,
    results: Dict[str, Dict[str, Any]],
    limits: Tuple[float, float],
) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 6.1))
    for color, split in zip(COLORS, ("train", "val", "test")):
        item = results[split]
        metrics = item["metrics"]
        label = f"{split.capitalize()} ($R^2$={metrics['r2']:.3f})"
        ax.scatter(
            item["true"],
            item["pred"],
            s=21,
            alpha=0.72,
            color=color,
            label=label,
            edgecolors="none",
        )
    ax.plot(limits, limits, "--", color="black", linewidth=1.2)
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"{MODEL_NAME}: All Splits", fontsize=TITLE_FONTSIZE)
    ax.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.tick_params(labelsize=TICK_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE)
    if USE_GRID:
        ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)
    fig.tight_layout()
    fig.savefig(path, dpi=FIG_DPI, format="jpg", bbox_inches="tight")
    plt.close(fig)


def plot_rmse_curve(history: Sequence[Dict[str, float]], best_epoch: int) -> None:
    frame = pd.DataFrame(history)
    frame.to_csv(
        DAT_DIR / f"{MODEL_NAME}_rmse_curve.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(frame["epoch"], frame["train_rmse_eV"], color=COLORS[0], label="Train RMSE")
    ax.plot(
        frame["epoch"], frame["val_rmse_eV"], color=COLORS[1], label="Validation RMSE"
    )
    ax.axvline(
        best_epoch,
        color=COLORS[3],
        linestyle="--",
        linewidth=1.2,
        label=f"Best epoch: {best_epoch}",
    )
    ax.set_title(f"{MODEL_NAME} Learning Curve", fontsize=TITLE_FONTSIZE)
    ax.set_xlabel("Epoch", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"RMSE ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.tick_params(labelsize=TICK_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE)
    if USE_GRID:
        ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_rmse_curve.jpg",
        dpi=FIG_DPI,
        format="jpg",
        bbox_inches="tight",
    )
    plt.close(fig)


def evaluate_and_save(
    x: np.ndarray,
    y: np.ndarray,
    records: Sequence[Record],
    splits: Dict[str, np.ndarray],
    checkpoint_path: Path,
    device: torch.device,
    history: Sequence[Dict[str, float]],
    best_epoch: int,
    logger: logging.Logger,
) -> Dict[str, Dict[str, Any]]:
    model, checkpoint = load_trained_model(checkpoint_path, device)
    n_features = int(checkpoint["model_config"]["n_features"])
    target_mean = float(checkpoint["preprocessing"]["target_mean"])
    target_scale = float(checkpoint["preprocessing"]["target_scale"])
    results: Dict[str, Dict[str, Any]] = {}
    all_names: List[str] = []
    all_true: List[np.ndarray] = []
    all_pred: List[np.ndarray] = []
    all_splits: List[str] = []

    for split in ("train", "val", "test"):
        indices = splits[split]
        loader = make_loader(x[:, :n_features], y, indices, BATCH_SIZE, False)
        y_true, y_pred = predict_loader(
            model, loader, device, target_mean, target_scale
        )
        names = [records[int(idx)].cif_name for idx in indices]
        metrics = calculate_metrics(y_true, y_pred)
        results[split] = {
            "true": y_true,
            "pred": y_pred,
            "metrics": metrics,
            "names": names,
        }
        save_prediction_dat(
            DAT_DIR / f"{MODEL_NAME}_parity_{split}.dat",
            split,
            names,
            y_true,
            y_pred,
        )
        all_names.extend(names)
        all_true.append(y_true)
        all_pred.append(y_pred)
        all_splits.extend([split] * len(names))
        logger.info(
            "%s | n=%d | MAE %.6f eV | RMSE %.6f eV | R2 %.6f",
            split.capitalize(),
            len(y_true),
            metrics["mae"],
            metrics["rmse"],
            metrics["r2"],
        )

    true_all = np.concatenate(all_true)
    pred_all = np.concatenate(all_pred)
    pd.DataFrame(
        {
            "split": all_splits,
            "cif_name": all_names,
            "true_bandgap_eV": true_all,
            "predicted_bandgap_eV": pred_all,
            "error_eV": pred_all - true_all,
            "absolute_error_eV": np.abs(pred_all - true_all),
        }
    ).to_csv(
        DAT_DIR / f"{MODEL_NAME}_parity_all.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )

    limits = _plot_limits(true_all, pred_all)
    for color, split in zip(COLORS, ("train", "val", "test")):
        plot_parity(
            FIGURE_DIR / f"{MODEL_NAME}_parity_{split}.jpg",
            split,
            results[split]["true"],
            results[split]["pred"],
            color,
            limits,
        )
    plot_combined_parity(FIGURE_DIR / f"{MODEL_NAME}_parity_all.jpg", results, limits)
    plot_rmse_curve(history, best_epoch)

    metric_rows = []
    for split in ("train", "val", "test"):
        metric_rows.append(
            {
                "split": split,
                "n_samples": len(results[split]["true"]),
                "mae_eV": results[split]["metrics"]["mae"],
                "rmse_eV": results[split]["metrics"]["rmse"],
                "r2": results[split]["metrics"]["r2"],
            }
        )
    pd.DataFrame(metric_rows).to_csv(
        TABLE_DIR / f"{MODEL_NAME}_metrics.dat",
        sep="\t",
        index=False,
        float_format="%.6f",
    )
    return results


def main() -> None:
    total_start = time.perf_counter()
    ensure_directories()
    logger = setup_logger()
    set_global_seed(SEED)
    device = get_device(logger)
    logger.info("Model: %s", MODEL_NAME)
    logger.info("Run version: %s", RUN_VERSION)
    logger.info("Output directory: %s", OUTPUT_DIR)
    logger.info("Excel path: %s", EXCEL_PATH)
    logger.info("CIF directory: %s", CIF_DIR)
    logger.info("Target: %s (%s)", TARGET_NAME, TARGET_UNIT)
    logger.info("Seed: %d", SEED)
    logger.info("Split ratios: %.2f/%.2f/%.2f", TRAIN_RATIO, VAL_RATIO, TEST_RATIO)
    logger.info("Loss: MSELoss")
    logger.info("Optuna enabled: %s", USE_OPTUNA)

    data_start = time.perf_counter()
    records = load_excel_records(logger)
    features, targets, records = compute_descriptors(records, logger)
    splits = create_or_load_split(records, logger)
    x, preprocessing = fit_preprocessor(features, targets, splits["train"], logger)
    data_seconds = time.perf_counter() - data_start
    logger.info("Raw descriptor count: %d", features.shape[1])
    logger.info("Selected descriptor count: %d", x.shape[1])
    logger.info(
        "Samples: train=%d, val=%d, test=%d",
        len(splits["train"]),
        len(splits["val"]),
        len(splits["test"]),
    )
    logger.info(
        "Data preparation time: %.3f s (%s)",
        data_seconds,
        format_duration(data_seconds),
    )

    if USE_OPTUNA:
        config = run_optuna(x, targets, splits, preprocessing, device, logger)
    else:
        config = default_train_config(x.shape[1])
    logger.info("Training configuration: %s", asdict(config))

    probe_model = MODNetRegressor(
        config.n_features,
        config.common_dims,
        config.intermediate_dims,
        config.property_dims,
        config.dropout,
    )
    total_params = sum(parameter.numel() for parameter in probe_model.parameters())
    trainable_params = sum(
        parameter.numel()
        for parameter in probe_model.parameters()
        if parameter.requires_grad
    )
    logger.info("Model architecture:\n%s", probe_model)
    logger.info("Total parameters: %d", total_params)
    logger.info("Trainable parameters: %d", trainable_params)

    set_global_seed(SEED)
    train_start = time.perf_counter()
    _, history, best_epoch, best_val_rmse, early_stopped = train_model(
        x,
        targets,
        splits,
        config,
        preprocessing,
        device,
        logger,
        checkpoint_path=CHECKPOINT_PATH,
        verbose=True,
    )
    train_seconds = time.perf_counter() - train_start
    logger.info("Best epoch: %d", best_epoch)
    logger.info("Best validation RMSE: %.6f eV", best_val_rmse)
    logger.info("Early stopping triggered: %s", early_stopped)
    logger.info("Best checkpoint: %s", CHECKPOINT_PATH)
    logger.info(
        "Training time: %.3f s (%s)", train_seconds, format_duration(train_seconds)
    )

    evaluation_start = time.perf_counter()
    evaluate_and_save(
        x,
        targets,
        records,
        splits,
        CHECKPOINT_PATH,
        device,
        history,
        best_epoch,
        logger,
    )
    evaluation_seconds = time.perf_counter() - evaluation_start
    total_seconds = time.perf_counter() - total_start
    logger.info(
        "Evaluation and plotting time: %.3f s (%s)",
        evaluation_seconds,
        format_duration(evaluation_seconds),
    )
    logger.info(
        "Total runtime: %.3f s (%s)", total_seconds, format_duration(total_seconds)
    )


if __name__ == "__main__":
    main()

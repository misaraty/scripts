from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import random
import time
import warnings
from concurrent.futures import ProcessPoolExecutor, as_completed
from contextlib import nullcontext
from dataclasses import asdict, dataclass
from datetime import timedelta
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
from scipy.special import gamma as gamma_function
from scipy.special import gammainc, gammaincc
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from torch_geometric.data import Batch, Data
from torch_geometric.loader import DataLoader
from torch_geometric.nn import MessagePassing, global_mean_pool
from tqdm import tqdm


# User configuration
MODEL_NAME = "PotNet"
RUN_VERSION = "v1"

EXCEL_PATH = "./data.xlsx"
CIF_DIR = "./cif"

TARGET_NAME = "bandgap"
TARGET_UNIT = "eV"

SEED = 42
TRAIN_RATIO = 0.8
VAL_RATIO = 0.1
TEST_RATIO = 0.1

LOCAL_CUTOFF = 4.0
MAX_NEIGHBORS = 16
POTENTIAL_R = 3
POTENTIAL_EPS = 1.0e-12
POTENTIAL_NAMES = ("coulomb", "dispersion", "pauli")
POTENTIAL_PARAMS = (0.5, 3.0, 3.0)
POTENTIAL_COEFFICIENTS = (-0.801, -0.074, 0.145)

HIDDEN_DIM = 256  # 256->128
NUM_CONV_LAYERS = 3
LOCAL_RBF_BINS = 256
POTENTIAL_RBF_BINS = 64  # 64->32
RBF_MIN = -4.0
RBF_MAX = 4.0
DROPOUT = 0.0

BATCH_SIZE = 32  # 32->4
MAX_EPOCHS = 300
LEARNING_RATE = 1.0e-3
WEIGHT_DECAY = 0.0
PATIENCE = 60
MIN_DELTA = 1.0e-6
LR_FACTOR = 0.5
LR_PATIENCE = 15
MIN_LR = 1.0e-6
NUM_WORKERS = 4
PREPROCESS_WORKERS = 4
PIN_MEMORY = True
NORMALIZE_TARGET = True
USE_AMP = True
AMP_DTYPE = "bfloat16"

USE_OPTUNA = False
OPTUNA_N_TRIALS = 30
OPTUNA_TIMEOUT = None
OPTUNA_MAX_EPOCHS = 180

FIG_DPI = 600
TITLE_FONTSIZE = 14
LABEL_FONTSIZE = 12
TICK_FONTSIZE = 10
LEGEND_FONTSIZE = 10
ANNOTATION_FONTSIZE = 10
USE_GRID = True
COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:red"]

OUTPUT_ROOT = Path(f"{MODEL_NAME}_{RUN_VERSION}")
FIGURE_DIR = OUTPUT_ROOT / "figure"
DAT_DIR = OUTPUT_ROOT / "dat"
TABLE_DIR = OUTPUT_ROOT / "table"
LOG_DIR = OUTPUT_ROOT / "log"
SPLIT_DIR = OUTPUT_ROOT / "split"
CACHE_DIR = OUTPUT_ROOT / "cache"
CHECKPOINT_PATH = OUTPUT_ROOT / f"{MODEL_NAME}_best.pt"
SPLIT_PATH = SPLIT_DIR / f"{MODEL_NAME}_split.csv"
LOG_PATH = LOG_DIR / f"{MODEL_NAME}_training.log"


@dataclass(frozen=True)
class GraphConfig:
    local_cutoff: float = LOCAL_CUTOFF
    max_neighbors: int = MAX_NEIGHBORS
    potential_r: int = POTENTIAL_R
    potential_eps: float = POTENTIAL_EPS
    potential_names: Tuple[str, ...] = POTENTIAL_NAMES
    potential_params: Tuple[float, ...] = POTENTIAL_PARAMS


@dataclass(frozen=True)
class ModelConfig:
    max_atomic_number: int = 118
    hidden_dim: int = HIDDEN_DIM
    num_conv_layers: int = NUM_CONV_LAYERS
    local_rbf_bins: int = LOCAL_RBF_BINS
    potential_rbf_bins: int = POTENTIAL_RBF_BINS
    rbf_min: float = RBF_MIN
    rbf_max: float = RBF_MAX
    potential_coefficients: Tuple[float, ...] = POTENTIAL_COEFFICIENTS
    dropout: float = DROPOUT


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


def set_global_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


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


def format_duration(seconds: float) -> str:
    return f"{seconds:.3f} s ({str(timedelta(seconds=int(round(seconds))))})"


def select_device(logger: logging.Logger) -> torch.device:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info("Device: %s", device)
    if device.type == "cuda":
        props = torch.cuda.get_device_properties(device)
        logger.info("GPU: %s", props.name)
        logger.info("GPU memory: %.2f GiB", props.total_memory / 1024**3)
        logger.info("CUDA version: %s", torch.version.cuda)
    return device


def normalize_cif_name(raw_name: Any) -> str:
    if pd.isna(raw_name):
        raise ValueError("Empty CIF filename")
    if isinstance(raw_name, (int, np.integer)):
        name = str(int(raw_name))
    elif isinstance(raw_name, (float, np.floating)) and float(raw_name).is_integer():
        name = str(int(raw_name))
    else:
        name = str(raw_name).strip()
    if not name:
        raise ValueError("Empty CIF filename")
    if not name.lower().endswith(".cif"):
        name += ".cif"
    return name


def load_excel_records(logger: logging.Logger) -> List[Dict[str, Any]]:
    excel_path = Path(EXCEL_PATH)
    if not excel_path.is_file():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    cif_dir = Path(CIF_DIR)
    if not cif_dir.is_dir():
        raise FileNotFoundError(f"CIF directory not found: {cif_dir}")
    frame = pd.read_excel(excel_path)
    if frame.shape[1] < 2:
        raise ValueError("The Excel file must contain at least two columns")

    records: List[Dict[str, Any]] = []
    seen: set[str] = set()
    for row_number, (_, row) in enumerate(frame.iterrows(), start=2):
        try:
            cif_name = normalize_cif_name(row.iloc[0])
            if cif_name in seen:
                raise ValueError(f"Duplicate CIF filename: {cif_name}")
            target = float(row.iloc[1])
            if not np.isfinite(target):
                raise ValueError("Target is not finite")
            cif_path = cif_dir / cif_name
            if not cif_path.is_file():
                raise FileNotFoundError(f"CIF file not found: {cif_path}")
            records.append(
                {"cif_name": cif_name, "cif_path": str(cif_path), "target": target}
            )
            seen.add(cif_name)
        except Exception as exc:
            logger.warning("Skipped Excel row %d: %s", row_number, exc)
    if len(records) < 10:
        raise ValueError(
            f"Only {len(records)} valid rows were found; at least 10 are required"
        )
    return records


def split_records(
    records: List[Dict[str, Any]], logger: logging.Logger
) -> Dict[str, List[Dict[str, Any]]]:
    expected = {record["cif_name"]: record for record in records}
    if SPLIT_PATH.is_file():
        try:
            saved = pd.read_csv(SPLIT_PATH)
            required = {"cif_name", "target", "split"}
            if not required.issubset(saved.columns):
                raise ValueError("Missing required split columns")
            if saved["cif_name"].duplicated().any():
                raise ValueError("Duplicate names in the saved split")
            if set(saved["cif_name"]) != set(expected):
                raise ValueError("Saved split does not match the current dataset")
            if not set(saved["split"]).issubset({"train", "val", "test"}):
                raise ValueError("Invalid split label")
            result = {key: [] for key in ("train", "val", "test")}
            for row in saved.itertuples(index=False):
                record = expected[str(row.cif_name)]
                if not math.isclose(
                    float(row.target), record["target"], rel_tol=0.0, abs_tol=1.0e-10
                ):
                    raise ValueError(f"Target changed for {row.cif_name}")
                result[str(row.split)].append(record)
            if any(len(result[key]) == 0 for key in result):
                raise ValueError("Saved split contains an empty subset")
            logger.info("Reused split: %s", SPLIT_PATH)
            return result
        except Exception as exc:
            logger.warning("Saved split was not reused: %s", exc)

    if not math.isclose(TRAIN_RATIO + VAL_RATIO + TEST_RATIO, 1.0, abs_tol=1.0e-12):
        raise ValueError("TRAIN_RATIO + VAL_RATIO + TEST_RATIO must equal 1")
    train_records, remainder = train_test_split(
        records, train_size=TRAIN_RATIO, random_state=SEED, shuffle=True
    )
    relative_val = VAL_RATIO / (VAL_RATIO + TEST_RATIO)
    val_records, test_records = train_test_split(
        remainder, train_size=relative_val, random_state=SEED, shuffle=True
    )
    result = {"train": train_records, "val": val_records, "test": test_records}
    rows = []
    for split_name in ("train", "val", "test"):
        rows.extend(
            {
                "cif_name": item["cif_name"],
                "target": item["target"],
                "split": split_name,
            }
            for item in result[split_name]
        )
    pd.DataFrame(rows).to_csv(SPLIT_PATH, index=False)
    logger.info("Created split: %s", SPLIT_PATH)
    return result


def _upper_incomplete_gamma(a: float, x: np.ndarray) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    if a > 0:
        return gammaincc(a, x) * gamma_function(a)
    shifted = a
    steps = 0
    while shifted <= 0:
        shifted += 1.0
        steps += 1
    value = gammaincc(shifted, x) * gamma_function(shifted)
    for _ in range(steps):
        shifted -= 1.0
        value = (value - np.power(x, shifted) * np.exp(-x)) / shifted
    return value


def _upper_bessel_k(nu: float, x: np.ndarray, y: float, eps: float) -> np.ndarray:
    x = np.asarray(x, dtype=np.float64)
    original_shape = x.shape
    flat_x = x.reshape(-1)
    result = np.zeros_like(flat_x)
    zero_x = np.isclose(flat_x, 0.0, atol=1.0e-15)

    if np.any(zero_x):
        if y <= 0 or nu <= 0:
            result[zero_x] = np.inf
        else:
            result[zero_x] = np.power(y, -nu) * gammainc(nu, y) * gamma_function(nu)

    active_indices = np.flatnonzero(~zero_x)
    if active_indices.size == 0:
        return result.reshape(original_shape)
    active_x = flat_x[active_indices]

    if abs(y) < 1.0e-15:
        result[active_indices] = np.power(active_x, nu) * _upper_incomplete_gamma(
            -nu, active_x
        )
        return result.reshape(original_shape)

    bounded = np.ones(active_x.size, dtype=bool)
    if -21.0 <= nu <= 21.0:
        bounded &= active_x <= 111.0
        bounded &= ~((active_x < y) & (active_x * y > 58.0**2))
    work_indices = active_indices[bounded]
    work_x = active_x[bounded]
    if work_x.size == 0:
        return result.reshape(original_shape)

    with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
        d0 = np.exp(work_x + y)
        n0 = np.zeros_like(work_x)
        n1 = np.ones_like(work_x)
        n2 = 0.5 * (work_x + nu + 3.0 - y) * n1
        n3 = ((work_x + nu + 5.0 - y) * n2 + (2.0 * y - nu - 2.0) * n1) / 3.0
        d1 = (work_x + nu + 1.0 - y) * d0
        d2 = 0.5 * (work_x + nu + 3.0 - y) * d1 + 0.5 * (2.0 * y - nu - 1.0) * d0
        d3 = ((work_x + nu + 5.0 - y) * d2 + (2.0 * y - nu - 2.0) * d1 - y * d0) / 3.0
        old = n2 / d2
        new = n3 / d3
        converged = np.abs(new - old) <= eps
        for n in range(3, 128):
            if np.all(converged):
                break
            old = np.where(converged, new, old)
            n0, n1, n2 = n1, n2, n3
            d0, d1, d2 = d1, d2, d3
            n3_candidate = (
                (work_x + nu + 1.0 + 2.0 * n - y) * n2
                + (2.0 * y - nu - n) * n1
                - y * n0
            ) / (n + 1.0)
            d3_candidate = (
                (work_x + nu + 1.0 + 2.0 * n - y) * d2
                + (2.0 * y - nu - n) * d1
                - y * d0
            ) / (n + 1.0)
            candidate = n3_candidate / d3_candidate
            valid = np.isfinite(candidate) & (np.abs(candidate) > 0.0)
            n3 = np.where(converged, n2, n3_candidate)
            d3 = np.where(converged, d2, d3_candidate)
            candidate = np.where(valid, candidate, new)
            delta = np.abs(candidate - new)
            new = np.where(converged, new, candidate)
            converged |= (delta <= eps) | (np.abs(new) < eps) | ~valid
    result[work_indices] = new
    return result.reshape(original_shape)


def _integer_grid(radius: int, include_zero: bool) -> np.ndarray:
    axes = [range(-radius, radius + 1)] * 3
    grid = np.asarray(
        [(i, j, k) for i in axes[0] for j in axes[1] for k in axes[2]], dtype=np.float64
    )
    if not include_zero:
        grid = grid[np.any(grid != 0.0, axis=1)]
    return grid


def periodic_zeta(
    displacements: np.ndarray,
    lattice: np.ndarray,
    param: float,
    radius: int,
    eps: float,
) -> np.ndarray:
    d = 3
    omega = np.asarray(lattice, dtype=np.float64)
    det = abs(float(np.linalg.det(omega)))
    if det <= 1.0e-14:
        raise ValueError("Singular lattice matrix")
    gamma_norm = det ** (1.0 / d)
    omega = omega / gamma_norm
    omega_inv = np.linalg.inv(omega).T
    v = np.asarray(displacements, dtype=np.float64) / gamma_norm
    products = _integer_grid(radius, include_zero=False)
    real_vectors = products @ omega
    reciprocal_vectors = products @ omega_inv

    real_x = np.sum(real_vectors**2, axis=1)
    real_sum = float(np.sum(_upper_bessel_k(-param, real_x, 0.0, eps)))
    reciprocal_x = math.pi**2 * np.sum(reciprocal_vectors**2, axis=1)
    reciprocal_kernel = _upper_bessel_k(param - d / 2.0, reciprocal_x, 0.0, eps)
    phases = 2.0 * math.pi * (v @ reciprocal_vectors.T)
    reciprocal_sum = math.pi ** (d / 2.0) * (np.cos(phases) @ reciprocal_kernel)

    v_norm2 = np.sum(v**2, axis=1)
    zero_displacement = np.isclose(v_norm2, 0.0, atol=1.0e-15)
    direct_zero = np.empty_like(v_norm2)
    direct_zero[zero_displacement] = -1.0 / param
    if np.any(~zero_displacement):
        direct_zero[~zero_displacement] = _upper_bessel_k(
            -param, v_norm2[~zero_displacement], 0.0, eps
        )
    reciprocal_zero = -math.pi ** (d / 2.0) / (d / 2.0 - param)
    values = real_sum + reciprocal_sum + direct_zero + reciprocal_zero
    return np.asarray(
        values * gamma_norm ** (-2.0 * param) / gamma_function(param), dtype=np.float64
    )


def periodic_exp(
    displacements: np.ndarray,
    lattice: np.ndarray,
    param: float,
    radius: int,
    eps: float,
) -> np.ndarray:
    d = 3
    omega = np.asarray(lattice, dtype=np.float64)
    det = abs(float(np.linalg.det(omega)))
    if det <= 1.0e-14:
        raise ValueError("Singular lattice matrix")
    gamma_norm = det ** (1.0 / d)
    omega = omega / gamma_norm
    omega_inv = np.linalg.inv(omega).T
    v = np.asarray(displacements, dtype=np.float64) / gamma_norm
    scaled_param = param * math.sqrt(gamma_norm)
    b_value = scaled_param**2 / (4.0 * math.pi)
    products = _integer_grid(radius, include_zero=True)
    real_vectors = products @ omega
    reciprocal_vectors = products @ omega_inv

    reciprocal_x = b_value + math.pi * np.sum(reciprocal_vectors**2, axis=1)
    reciprocal_kernel = _upper_bessel_k(-0.5 - d / 2.0, reciprocal_x, 0.0, eps)
    phases = 2.0 * math.pi * (v @ reciprocal_vectors.T)
    reciprocal_sum = np.cos(phases) @ reciprocal_kernel

    shifted = real_vectors[None, :, :] + v[:, None, :]
    real_x = math.pi * np.sum(shifted**2, axis=2)
    real_sum = np.sum(_upper_bessel_k(0.5, real_x, b_value, eps), axis=1)
    return np.asarray(
        (reciprocal_sum + real_sum) * scaled_param / (2.0 * math.pi), dtype=np.float64
    )


def compute_periodic_potentials(
    displacements: np.ndarray, lattice: np.ndarray, config: GraphConfig
) -> np.ndarray:
    values = []
    for name, param in zip(config.potential_names, config.potential_params):
        if name in {"coulomb", "dispersion", "zeta"}:
            values.append(
                periodic_zeta(
                    displacements,
                    lattice,
                    param,
                    config.potential_r,
                    config.potential_eps,
                )
            )
        elif name in {"pauli", "exp"}:
            values.append(
                periodic_exp(
                    displacements,
                    lattice,
                    param,
                    config.potential_r,
                    config.potential_eps,
                )
            )
        else:
            raise ValueError(f"Unsupported potential: {name}")
    result = np.stack(values, axis=1)
    if not np.all(np.isfinite(result)):
        raise FloatingPointError(
            "Periodic potential calculation produced non-finite values"
        )
    return result.astype(np.float32)


def _limit_neighbors(
    centers: np.ndarray,
    neighbors: np.ndarray,
    distances: np.ndarray,
    max_neighbors: int,
) -> np.ndarray:
    keep: List[int] = []
    for center in np.unique(centers):
        indices = np.flatnonzero(centers == center)
        order = indices[np.argsort(distances[indices], kind="stable")]
        keep.extend(order[:max_neighbors].tolist())
    return np.asarray(keep, dtype=np.int64)


def build_potnet_graph(cif_path: str, graph_config: GraphConfig) -> Data:
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"Issues encountered while parsing CIF:.*",
            category=UserWarning,
        )
        structure = Structure.from_file(cif_path)
    if len(structure) == 0:
        raise ValueError("Structure contains no atoms")
    if not structure.is_ordered:
        raise ValueError("Disordered structures are not supported")
    atomic_numbers = np.asarray([site.specie.Z for site in structure], dtype=np.int64)
    if np.any((atomic_numbers < 1) | (atomic_numbers > 118)):
        raise ValueError("Atomic number is outside the supported range 1-118")

    centers, neighbors, _, distances = structure.get_neighbor_list(
        graph_config.local_cutoff
    )
    centers = np.asarray(centers, dtype=np.int64)
    neighbors = np.asarray(neighbors, dtype=np.int64)
    distances = np.asarray(distances, dtype=np.float64)
    nonzero = distances > 1.0e-8
    centers, neighbors, distances = (
        centers[nonzero],
        neighbors[nonzero],
        distances[nonzero],
    )
    if centers.size == 0:
        raise ValueError("No local periodic neighbors were found")
    keep = _limit_neighbors(centers, neighbors, distances, graph_config.max_neighbors)
    local_edge_index = np.stack([neighbors[keep], centers[keep]], axis=0)
    local_edge_distance = distances[keep].astype(np.float32)

    n_atoms = len(structure)
    sources = np.repeat(np.arange(n_atoms, dtype=np.int64), n_atoms)
    targets = np.tile(np.arange(n_atoms, dtype=np.int64), n_atoms)
    cart_coords = np.asarray(structure.cart_coords, dtype=np.float64)
    displacements = cart_coords[sources] - cart_coords[targets]
    potential_features = compute_periodic_potentials(
        displacements,
        np.asarray(structure.lattice.matrix, dtype=np.float64),
        graph_config,
    )
    infinite_edge_index = np.stack([sources, targets], axis=0)

    return Data(
        z=torch.from_numpy(atomic_numbers),
        edge_index=torch.from_numpy(local_edge_index).long(),
        edge_distance=torch.from_numpy(local_edge_distance).float(),
        inf_edge_index=torch.from_numpy(infinite_edge_index).long(),
        inf_edge_attr=torch.from_numpy(potential_features).float(),
        num_nodes=n_atoms,
    )


def graph_to_numpy_payload(graph: Data) -> Dict[str, Any]:
    return {
        "z": graph.z.cpu().numpy(),
        "edge_index": graph.edge_index.cpu().numpy(),
        "edge_distance": graph.edge_distance.cpu().numpy(),
        "inf_edge_index": graph.inf_edge_index.cpu().numpy(),
        "inf_edge_attr": graph.inf_edge_attr.cpu().numpy(),
        "num_nodes": int(graph.num_nodes),
    }


def numpy_payload_to_graph(payload: Dict[str, Any]) -> Data:
    return Data(
        z=torch.from_numpy(np.asarray(payload["z"], dtype=np.int64)).long(),
        edge_index=torch.from_numpy(
            np.asarray(payload["edge_index"], dtype=np.int64)
        ).long(),
        edge_distance=torch.from_numpy(
            np.asarray(payload["edge_distance"], dtype=np.float32)
        ).float(),
        inf_edge_index=torch.from_numpy(
            np.asarray(payload["inf_edge_index"], dtype=np.int64)
        ).long(),
        inf_edge_attr=torch.from_numpy(
            np.asarray(payload["inf_edge_attr"], dtype=np.float32)
        ).float(),
        num_nodes=int(payload["num_nodes"]),
    )


def _build_graph_worker(
    payload: Tuple[int, Dict[str, Any], GraphConfig]
) -> Tuple[int, Optional[Dict[str, Any]], Optional[str]]:
    index, record, graph_config = payload
    try:
        graph = build_potnet_graph(record["cif_path"], graph_config)
        return index, graph_to_numpy_payload(graph), None
    except Exception as exc:
        return index, None, f"{type(exc).__name__}: {exc}"


def graph_cache_path(
    records: Sequence[Dict[str, Any]], graph_config: GraphConfig
) -> Path:
    payload = {
        "records": [
            (item["cif_name"], item["target"], os.path.getmtime(item["cif_path"]))
            for item in records
        ],
        "graph_config": asdict(graph_config),
        "format": 2,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode("utf-8")
    ).hexdigest()[:16]
    return CACHE_DIR / f"{MODEL_NAME}_graphs_{digest}.pt"


def safe_torch_load(path: Path, map_location: Any = "cpu") -> Any:
    try:
        return torch.load(path, map_location=map_location, weights_only=False)
    except TypeError:
        return torch.load(path, map_location=map_location)


def build_or_load_graphs(
    records: List[Dict[str, Any]], graph_config: GraphConfig, logger: logging.Logger
) -> Tuple[List[Dict[str, Any]], List[Data]]:
    cache_path = graph_cache_path(records, graph_config)
    if cache_path.is_file():
        try:
            payload = safe_torch_load(cache_path)
            source_names = [item["cif_name"] for item in records]
            if payload["source_cif_names"] == source_names:
                record_by_name = {item["cif_name"]: item for item in records}
                valid_records = [
                    record_by_name[name] for name in payload["valid_cif_names"]
                ]
                logger.info(
                    "Loaded %d graphs from %s", len(payload["graphs"]), cache_path
                )
                return valid_records, payload["graphs"]
        except Exception as exc:
            logger.warning("Graph cache was not used: %s", exc)

    tasks = [(index, record, graph_config) for index, record in enumerate(records)]
    completed: Dict[int, Data] = {}
    failures: Dict[int, str] = {}
    if PREPROCESS_WORKERS > 1:
        try:
            with ProcessPoolExecutor(max_workers=PREPROCESS_WORKERS) as executor:
                futures = [executor.submit(_build_graph_worker, task) for task in tasks]
                for future in tqdm(
                    as_completed(futures),
                    total=len(futures),
                    desc="Building PotNet graphs",
                ):
                    index, graph_payload, error = future.result()
                    if graph_payload is None:
                        failures[index] = error or "Unknown graph error"
                    else:
                        completed[index] = numpy_payload_to_graph(graph_payload)
        except Exception as exc:
            logger.warning(
                "Parallel graph construction failed (%s: %s). Remaining graphs will be built sequentially.",
                type(exc).__name__,
                exc,
            )
            remaining_tasks = [
                task
                for task in tasks
                if task[0] not in completed and task[0] not in failures
            ]
            for task in tqdm(remaining_tasks, desc="Building remaining graphs"):
                index, graph_payload, error = _build_graph_worker(task)
                if graph_payload is None:
                    failures[index] = error or "Unknown graph error"
                else:
                    completed[index] = numpy_payload_to_graph(graph_payload)
    else:
        for task in tqdm(tasks, desc="Building PotNet graphs"):
            index, graph_payload, error = _build_graph_worker(task)
            if graph_payload is None:
                failures[index] = error or "Unknown graph error"
            else:
                completed[index] = numpy_payload_to_graph(graph_payload)

    valid_records: List[Dict[str, Any]] = []
    graphs: List[Data] = []
    for index, record in enumerate(records):
        if index in completed:
            valid_records.append(record)
            graphs.append(completed[index])
        else:
            logger.warning(
                "Skipped %s during graph construction: %s",
                record["cif_name"],
                failures[index],
            )
    if len(graphs) < 10:
        raise ValueError(f"Only {len(graphs)} CIF files produced valid PotNet graphs")
    torch.save(
        {
            "source_cif_names": [item["cif_name"] for item in records],
            "valid_cif_names": [item["cif_name"] for item in valid_records],
            "graphs": graphs,
        },
        cache_path,
    )
    logger.info("Saved %d graphs to %s", len(graphs), cache_path)
    return valid_records, graphs


class PotNetDataset(torch.utils.data.Dataset):
    def __init__(
        self,
        records: Sequence[Dict[str, Any]],
        graph_by_name: Dict[str, Data],
        sample_id_by_name: Dict[str, int],
        target_mean: float,
        target_std: float,
    ) -> None:
        self.records = list(records)
        self.graph_by_name = graph_by_name
        self.sample_id_by_name = sample_id_by_name
        self.target_mean = float(target_mean)
        self.target_std = float(target_std)

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, index: int) -> Data:
        record = self.records[index]
        graph = self.graph_by_name[record["cif_name"]].clone()
        target = float(record["target"])
        graph.y = torch.tensor(
            [(target - self.target_mean) / self.target_std], dtype=torch.float32
        )
        graph.y_raw = torch.tensor([target], dtype=torch.float32)
        graph.sample_id = torch.tensor(
            [self.sample_id_by_name[record["cif_name"]]], dtype=torch.long
        )
        return graph


class RBFExpansion(nn.Module):
    def __init__(
        self, vmin: float, vmax: float, bins: int, basis_type: str = "gaussian"
    ) -> None:
        super().__init__()
        centers = torch.linspace(vmin, vmax, bins)
        self.register_buffer("centers", centers)
        spacing = float(centers[1] - centers[0]) if bins > 1 else 1.0
        self.gamma = 1.0 / spacing
        self.basis_type = basis_type

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        base = self.gamma * (values.reshape(-1, 1) - self.centers)
        if self.basis_type == "gaussian":
            return torch.exp(-(base**2))
        if self.basis_type == "multiquadric":
            return torch.sqrt(1.0 + base**2)
        raise ValueError(f"Unsupported RBF type: {self.basis_type}")


class SafeBatchNorm1d(nn.BatchNorm1d):
    def forward(self, values: torch.Tensor) -> torch.Tensor:
        if self.training and values.shape[0] <= 1:
            return F.batch_norm(
                values,
                self.running_mean,
                self.running_var,
                self.weight,
                self.bias,
                False,
                self.momentum,
                self.eps,
            )
        return super().forward(values)


class PotNetConv(MessagePassing):
    def __init__(self, hidden_dim: int) -> None:
        super().__init__(aggr="add", node_dim=0)
        self.output_norm = SafeBatchNorm1d(hidden_dim)
        self.gate_norm = SafeBatchNorm1d(hidden_dim)
        self.gate_network = nn.Sequential(
            nn.Linear(3 * hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.message_network = nn.Sequential(
            nn.Linear(3 * hidden_dim, hidden_dim),
            nn.SiLU(),
            nn.Linear(hidden_dim, hidden_dim),
        )

    def forward(
        self, x: torch.Tensor, edge_index: torch.Tensor, edge_attr: torch.Tensor
    ) -> torch.Tensor:
        output = self.propagate(
            edge_index, x=x, edge_attr=edge_attr, size=(x.size(0), x.size(0))
        )
        return F.relu(x + self.output_norm(output))

    def message(
        self, x_i: torch.Tensor, x_j: torch.Tensor, edge_attr: torch.Tensor
    ) -> torch.Tensor:
        merged = torch.cat((x_i, x_j, edge_attr), dim=-1)
        gate = torch.sigmoid(self.gate_norm(self.gate_network(merged)))
        return gate * self.message_network(merged)


class ShiftedSoftplus(nn.Module):
    def forward(self, values: torch.Tensor) -> torch.Tensor:
        return F.softplus(values) - math.log(2.0)


class PotNet(nn.Module):
    def __init__(self, config: ModelConfig) -> None:
        super().__init__()
        self.config = config
        self.atom_embedding = nn.Embedding(
            config.max_atomic_number + 1, config.hidden_dim, padding_idx=0
        )
        self.local_edge_embedding = nn.Sequential(
            RBFExpansion(
                config.rbf_min, config.rbf_max, config.local_rbf_bins, "gaussian"
            ),
            nn.Linear(config.local_rbf_bins, config.hidden_dim),
            nn.SiLU(),
        )
        self.potential_rbf = RBFExpansion(
            config.rbf_min, config.rbf_max, config.potential_rbf_bins, "multiquadric"
        )
        self.potential_projection = nn.Linear(
            config.potential_rbf_bins, config.hidden_dim
        )
        self.potential_norm = SafeBatchNorm1d(config.hidden_dim)
        self.convolutions = nn.ModuleList(
            [PotNetConv(config.hidden_dim) for _ in range(config.num_conv_layers)]
        )
        self.readout = nn.Sequential(
            nn.Linear(config.hidden_dim, config.hidden_dim),
            ShiftedSoftplus(),
            nn.Dropout(config.dropout),
            nn.Linear(config.hidden_dim, 1),
        )
        self.register_buffer(
            "potential_coefficients",
            torch.tensor(config.potential_coefficients, dtype=torch.float32),
        )
        self.reset_parameters()

    def reset_parameters(self) -> None:
        with torch.no_grad():
            nn.init.xavier_uniform_(self.atom_embedding.weight[1:])
            self.atom_embedding.weight[0].zero_()

    def forward(self, data: Batch) -> torch.Tensor:
        node_features = self.atom_embedding(data.z.long())
        local_values = -0.75 / data.edge_distance.clamp_min(1.0e-8)
        local_edge_features = self.local_edge_embedding(local_values)
        combined_potential = torch.sum(
            data.inf_edge_attr * self.potential_coefficients, dim=-1
        )
        infinite_edge_features = self.potential_rbf(combined_potential)
        infinite_edge_features = self.potential_norm(
            F.softplus(self.potential_projection(infinite_edge_features))
        )
        edge_index = torch.cat((data.edge_index, data.inf_edge_index), dim=1)
        edge_features = torch.cat((local_edge_features, infinite_edge_features), dim=0)
        for convolution in self.convolutions:
            node_features = convolution(node_features, edge_index, edge_features)
        graph_features = global_mean_pool(node_features, data.batch)
        return self.readout(graph_features).view(-1)


def count_parameters(model: nn.Module) -> Tuple[int, int]:
    total = sum(parameter.numel() for parameter in model.parameters())
    trainable = sum(
        parameter.numel() for parameter in model.parameters() if parameter.requires_grad
    )
    return total, trainable


def amp_context(device: torch.device):
    if not (USE_AMP and device.type == "cuda"):
        return nullcontext()
    dtype = torch.bfloat16 if AMP_DTYPE.lower() == "bfloat16" else torch.float16
    return torch.autocast(device_type="cuda", dtype=dtype)


def train_one_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    criterion: nn.Module,
    device: torch.device,
    target_mean: float,
    target_std: float,
) -> Tuple[float, float]:
    model.train()
    squared_error_sum = 0.0
    normalized_loss_sum = 0.0
    sample_count = 0
    for batch in loader:
        batch = batch.to(device, non_blocking=True)
        optimizer.zero_grad(set_to_none=True)
        with amp_context(device):
            prediction = model(batch).view(-1)
            target = batch.y.view(-1)
            loss = criterion(prediction, target)
        if not torch.isfinite(loss):
            raise FloatingPointError("Non-finite training loss")
        loss.backward()
        torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=10.0)
        optimizer.step()
        batch_size = target.numel()
        prediction_raw = prediction.detach().float() * target_std + target_mean
        target_raw = target.detach().float() * target_std + target_mean
        squared_error_sum += float(torch.sum((prediction_raw - target_raw) ** 2).cpu())
        normalized_loss_sum += float(loss.detach().cpu()) * batch_size
        sample_count += batch_size
    return (
        math.sqrt(squared_error_sum / sample_count),
        normalized_loss_sum / sample_count,
    )


@torch.inference_mode()
def evaluate_loader(
    model: nn.Module,
    loader: DataLoader,
    criterion: nn.Module,
    device: torch.device,
    target_mean: float,
    target_std: float,
) -> Dict[str, Any]:
    model.eval()
    predictions: List[np.ndarray] = []
    targets: List[np.ndarray] = []
    sample_ids: List[np.ndarray] = []
    normalized_loss_sum = 0.0
    sample_count = 0
    for batch in loader:
        batch = batch.to(device, non_blocking=True)
        with amp_context(device):
            prediction = model(batch).view(-1)
            target = batch.y.view(-1)
            loss = criterion(prediction, target)
        prediction_raw = prediction.float() * target_std + target_mean
        target_raw = target.float() * target_std + target_mean
        predictions.append(prediction_raw.cpu().numpy())
        targets.append(target_raw.cpu().numpy())
        sample_ids.append(batch.sample_id.view(-1).cpu().numpy())
        batch_size = target.numel()
        normalized_loss_sum += float(loss.float().cpu()) * batch_size
        sample_count += batch_size
    y_pred = np.concatenate(predictions)
    y_true = np.concatenate(targets)
    indices = np.concatenate(sample_ids).astype(np.int64)
    return {
        "y_true": y_true,
        "y_pred": y_pred,
        "sample_ids": indices,
        "rmse": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "normalized_mse": normalized_loss_sum / sample_count,
    }


def make_loaders(
    datasets: Dict[str, PotNetDataset], batch_size: int, device: torch.device
) -> Dict[str, DataLoader]:
    loader_kwargs = {
        "batch_size": batch_size,
        "num_workers": NUM_WORKERS,
        "pin_memory": PIN_MEMORY and device.type == "cuda",
        "persistent_workers": NUM_WORKERS > 0,
    }
    generator = torch.Generator().manual_seed(SEED)
    return {
        "train": DataLoader(
            datasets["train"], shuffle=True, generator=generator, **loader_kwargs
        ),
        "val": DataLoader(datasets["val"], shuffle=False, **loader_kwargs),
        "test": DataLoader(datasets["test"], shuffle=False, **loader_kwargs),
    }


def save_checkpoint(
    path: Path,
    model: nn.Module,
    model_config: ModelConfig,
    graph_config: GraphConfig,
    best_epoch: int,
    best_val_rmse: float,
    target_mean: float,
    target_std: float,
    training_config: Dict[str, Any],
) -> None:
    torch.save(
        {
            "model_name": MODEL_NAME,
            "run_version": RUN_VERSION,
            "model_state_dict": model.state_dict(),
            "model_config": asdict(model_config),
            "graph_config": asdict(graph_config),
            "best_epoch": best_epoch,
            "best_val_rmse": best_val_rmse,
            "seed": SEED,
            "target_name": TARGET_NAME,
            "target_unit": TARGET_UNIT,
            "target_mean": target_mean,
            "target_std": target_std,
            "training_config": training_config,
        },
        path,
    )


def train_model(
    model_config: ModelConfig,
    graph_config: GraphConfig,
    datasets: Dict[str, PotNetDataset],
    device: torch.device,
    target_mean: float,
    target_std: float,
    logger: logging.Logger,
    training_config: Dict[str, Any],
    checkpoint_path: Optional[Path],
    trial: Any = None,
) -> Tuple[nn.Module, List[Dict[str, float]], int, float, bool]:
    set_global_seed(SEED)
    model = PotNet(model_config).to(device)
    loaders = make_loaders(datasets, int(training_config["batch_size"]), device)
    criterion = nn.MSELoss()
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=float(training_config["learning_rate"]),
        weight_decay=float(training_config["weight_decay"]),
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="min",
        factor=LR_FACTOR,
        patience=LR_PATIENCE,
        min_lr=MIN_LR,
    )
    history: List[Dict[str, float]] = []
    best_val_rmse = math.inf
    best_epoch = 0
    epochs_without_improvement = 0
    early_stopped = False
    best_state: Optional[Dict[str, torch.Tensor]] = None

    for epoch in range(1, int(training_config["max_epochs"]) + 1):
        train_rmse, train_mse = train_one_epoch(
            model,
            loaders["train"],
            optimizer,
            criterion,
            device,
            target_mean,
            target_std,
        )
        val_output = evaluate_loader(
            model, loaders["val"], criterion, device, target_mean, target_std
        )
        val_rmse = float(val_output["rmse"])
        learning_rate = float(optimizer.param_groups[0]["lr"])
        history.append(
            {
                "epoch": float(epoch),
                "train_rmse_eV": train_rmse,
                "val_rmse_eV": val_rmse,
                "learning_rate": learning_rate,
                "train_mse_loss": train_mse,
                "val_mse_loss": float(val_output["normalized_mse"]),
            }
        )
        logger.info(
            "Epoch %04d | train RMSE %.6f %s | val RMSE %.6f %s | lr %.6e",
            epoch,
            train_rmse,
            TARGET_UNIT,
            val_rmse,
            TARGET_UNIT,
            learning_rate,
        )
        scheduler.step(val_rmse)

        if val_rmse < best_val_rmse - MIN_DELTA:
            best_val_rmse = val_rmse
            best_epoch = epoch
            epochs_without_improvement = 0
            best_state = {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            }
            if checkpoint_path is not None:
                save_checkpoint(
                    checkpoint_path,
                    model,
                    model_config,
                    graph_config,
                    best_epoch,
                    best_val_rmse,
                    target_mean,
                    target_std,
                    training_config,
                )
        else:
            epochs_without_improvement += 1

        if trial is not None:
            trial.report(best_val_rmse, epoch)
            if trial.should_prune():
                import optuna

                raise optuna.TrialPruned()
        if epochs_without_improvement >= PATIENCE:
            early_stopped = True
            logger.info("Early stopping triggered at epoch %d", epoch)
            break

    if best_state is None:
        raise RuntimeError("Training ended without a valid model state")
    model.load_state_dict(best_state)
    return model, history, best_epoch, best_val_rmse, early_stopped


def run_optuna(
    base_model_config: ModelConfig,
    graph_config: GraphConfig,
    datasets: Dict[str, PotNetDataset],
    device: torch.device,
    target_mean: float,
    target_std: float,
    logger: logging.Logger,
) -> Tuple[ModelConfig, Dict[str, Any]]:
    try:
        import optuna
    except ImportError as exc:
        raise ImportError("Optuna is required when USE_OPTUNA=True") from exc

    def objective(trial: Any) -> float:
        model_config = ModelConfig(
            hidden_dim=trial.suggest_categorical("hidden_dim", [128, 192, 256, 320]),
            num_conv_layers=trial.suggest_int("num_conv_layers", 2, 5),
            local_rbf_bins=base_model_config.local_rbf_bins,
            potential_rbf_bins=base_model_config.potential_rbf_bins,
            rbf_min=base_model_config.rbf_min,
            rbf_max=base_model_config.rbf_max,
            potential_coefficients=base_model_config.potential_coefficients,
            dropout=trial.suggest_float("dropout", 0.0, 0.25),
        )
        training_config = {
            "batch_size": trial.suggest_categorical("batch_size", [16, 32, 64]),
            "learning_rate": trial.suggest_float(
                "learning_rate", 2.0e-4, 3.0e-3, log=True
            ),
            "weight_decay": trial.suggest_float(
                "weight_decay", 1.0e-7, 1.0e-3, log=True
            ),
            "max_epochs": OPTUNA_MAX_EPOCHS,
        }
        _, _, _, best_rmse, _ = train_model(
            model_config,
            graph_config,
            datasets,
            device,
            target_mean,
            target_std,
            logger,
            training_config,
            checkpoint_path=None,
            trial=trial,
        )
        if device.type == "cuda":
            torch.cuda.empty_cache()
        return best_rmse

    sampler = optuna.samplers.TPESampler(seed=SEED)
    study = optuna.create_study(direction="minimize", sampler=sampler)
    study.optimize(objective, n_trials=OPTUNA_N_TRIALS, timeout=OPTUNA_TIMEOUT)
    best_params = dict(study.best_params)
    best_params["best_val_rmse"] = float(study.best_value)
    with open(
        TABLE_DIR / f"{MODEL_NAME}_optuna_best_params.json", "w", encoding="utf-8"
    ) as handle:
        json.dump(best_params, handle, indent=2)
    logger.info("Optuna best validation RMSE: %.6f %s", study.best_value, TARGET_UNIT)
    logger.info("Optuna best parameters: %s", study.best_params)
    tuned_model_config = ModelConfig(
        hidden_dim=int(study.best_params["hidden_dim"]),
        num_conv_layers=int(study.best_params["num_conv_layers"]),
        local_rbf_bins=base_model_config.local_rbf_bins,
        potential_rbf_bins=base_model_config.potential_rbf_bins,
        rbf_min=base_model_config.rbf_min,
        rbf_max=base_model_config.rbf_max,
        potential_coefficients=base_model_config.potential_coefficients,
        dropout=float(study.best_params["dropout"]),
    )
    training_config = {
        "batch_size": int(study.best_params["batch_size"]),
        "learning_rate": float(study.best_params["learning_rate"]),
        "weight_decay": float(study.best_params["weight_decay"]),
        "max_epochs": MAX_EPOCHS,
    }
    return tuned_model_config, training_config


def load_trained_model(
    checkpoint_path: str | Path, device: Optional[torch.device] = None
) -> Tuple[PotNet, Dict[str, Any], torch.device]:
    selected_device = device or torch.device(
        "cuda" if torch.cuda.is_available() else "cpu"
    )
    checkpoint = safe_torch_load(Path(checkpoint_path), map_location=selected_device)
    config_data = dict(checkpoint["model_config"])
    config_data["potential_coefficients"] = tuple(config_data["potential_coefficients"])
    model = PotNet(ModelConfig(**config_data)).to(selected_device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    return model, checkpoint, selected_device


def calculate_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    return {
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "r2": float(r2_score(y_true, y_pred)),
    }


def evaluate_splits(
    model: PotNet,
    loaders: Dict[str, DataLoader],
    device: torch.device,
    target_mean: float,
    target_std: float,
    records: List[Dict[str, Any]],
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Dict[str, float]]]:
    criterion = nn.MSELoss()
    outputs: Dict[str, Dict[str, Any]] = {}
    metrics: Dict[str, Dict[str, float]] = {}
    for split_name in ("train", "val", "test"):
        output = evaluate_loader(
            model, loaders[split_name], criterion, device, target_mean, target_std
        )
        order = np.argsort(output["sample_ids"])
        output = {
            key: (
                value[order]
                if isinstance(value, np.ndarray) and value.shape[0] == order.shape[0]
                else value
            )
            for key, value in output.items()
        }
        output["cif_names"] = [
            records[index]["cif_name"] for index in output["sample_ids"]
        ]
        outputs[split_name] = output
        metrics[split_name] = calculate_metrics(output["y_true"], output["y_pred"])
    return outputs, metrics


def save_metrics_table(
    metrics: Dict[str, Dict[str, float]], outputs: Dict[str, Dict[str, Any]]
) -> None:
    rows = []
    for split_name in ("train", "val", "test"):
        rows.append(
            {
                "split": split_name,
                "n_samples": len(outputs[split_name]["y_true"]),
                "mae_eV": metrics[split_name]["mae"],
                "rmse_eV": metrics[split_name]["rmse"],
                "r2": metrics[split_name]["r2"],
            }
        )
    pd.DataFrame(rows).to_csv(
        TABLE_DIR / f"{MODEL_NAME}_metrics.dat",
        sep="\t",
        index=False,
        float_format="%.6f",
    )


def save_prediction_dat(split_name: str, output: Dict[str, Any]) -> None:
    frame = pd.DataFrame(
        {
            "cif_name": output["cif_names"],
            "true_bandgap_eV": output["y_true"],
            "predicted_bandgap_eV": output["y_pred"],
            "error_eV": output["y_pred"] - output["y_true"],
            "absolute_error_eV": np.abs(output["y_pred"] - output["y_true"]),
        }
    )
    frame.to_csv(
        DAT_DIR / f"{MODEL_NAME}_parity_{split_name}.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )


def save_all_prediction_dat(outputs: Dict[str, Dict[str, Any]]) -> None:
    frames = []
    for split_name in ("train", "val", "test"):
        output = outputs[split_name]
        frames.append(
            pd.DataFrame(
                {
                    "split": split_name,
                    "cif_name": output["cif_names"],
                    "true_bandgap_eV": output["y_true"],
                    "predicted_bandgap_eV": output["y_pred"],
                    "error_eV": output["y_pred"] - output["y_true"],
                    "absolute_error_eV": np.abs(output["y_pred"] - output["y_true"]),
                }
            )
        )
    pd.concat(frames, ignore_index=True).to_csv(
        DAT_DIR / f"{MODEL_NAME}_parity_all.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )


def parity_limits(outputs: Iterable[Dict[str, Any]]) -> Tuple[float, float]:
    values = np.concatenate(
        [
            np.concatenate((np.asarray(output["y_true"]), np.asarray(output["y_pred"])))
            for output in outputs
        ]
    )
    low, high = float(np.min(values)), float(np.max(values))
    margin = max(0.05 * (high - low), 0.05)
    return low - margin, high + margin


def style_axes(ax: plt.Axes) -> None:
    ax.tick_params(labelsize=TICK_FONTSIZE)
    if USE_GRID:
        ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)


def plot_parity(
    split_name: str,
    output: Dict[str, Any],
    metric: Dict[str, float],
    limits: Tuple[float, float],
    color: str,
) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    ax.scatter(
        output["y_true"],
        output["y_pred"],
        s=22,
        alpha=0.75,
        color=color,
        edgecolors="none",
    )
    ax.plot(limits, limits, linestyle="--", color="black", linewidth=1.2)
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_title(
        f"{MODEL_NAME}: {split_name.capitalize()} parity", fontsize=TITLE_FONTSIZE
    )
    annotation = (
        f"MAE = {metric['mae']:.4f} {TARGET_UNIT}\n"
        f"RMSE = {metric['rmse']:.4f} {TARGET_UNIT}\n"
        f"R2 = {metric['r2']:.4f}"
    )
    ax.text(
        0.04,
        0.96,
        annotation,
        transform=ax.transAxes,
        va="top",
        fontsize=ANNOTATION_FONTSIZE,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85},
    )
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_parity_{split_name}.jpg",
        dpi=FIG_DPI,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_combined_parity(
    outputs: Dict[str, Dict[str, Any]],
    metrics: Dict[str, Dict[str, float]],
    limits: Tuple[float, float],
) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 6.2))
    for index, split_name in enumerate(("train", "val", "test")):
        output = outputs[split_name]
        ax.scatter(
            output["y_true"],
            output["y_pred"],
            s=20,
            alpha=0.7,
            color=COLORS[index],
            edgecolors="none",
            label=(
                f"{split_name.capitalize()} "
                f"(RMSE={metrics[split_name]['rmse']:.4f}, R2={metrics[split_name]['r2']:.4f})"
            ),
        )
    ax.plot(limits, limits, linestyle="--", color="black", linewidth=1.2, label="y = x")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_title(f"{MODEL_NAME}: Combined parity", fontsize=TITLE_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_parity_all.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(fig)


def save_and_plot_history(history: List[Dict[str, float]], best_epoch: int) -> None:
    frame = pd.DataFrame(history)
    frame["epoch"] = frame["epoch"].astype(int)
    frame.to_csv(
        DAT_DIR / f"{MODEL_NAME}_rmse_curve.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
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
    ax.set_xlabel("Epoch", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"RMSE ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_title(f"{MODEL_NAME}: Training history", fontsize=TITLE_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_rmse_curve.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(fig)


@torch.inference_mode()
def predict_cifs(
    checkpoint_path: str | Path,
    cif_paths: Sequence[str | Path],
    device: Optional[torch.device] = None,
) -> pd.DataFrame:
    model, checkpoint, selected_device = load_trained_model(checkpoint_path, device)
    graph_data = dict(checkpoint["graph_config"])
    graph_data["potential_names"] = tuple(graph_data["potential_names"])
    graph_data["potential_params"] = tuple(graph_data["potential_params"])
    graph_config = GraphConfig(**graph_data)
    graphs = []
    names = []
    for index, path in enumerate(cif_paths):
        graph = build_potnet_graph(str(path), graph_config)
        graph.sample_id = torch.tensor([index], dtype=torch.long)
        graphs.append(graph)
        names.append(Path(path).name)
    loader = DataLoader(graphs, batch_size=BATCH_SIZE, shuffle=False)
    predictions = []
    indices = []
    for batch in loader:
        batch = batch.to(selected_device)
        with amp_context(selected_device):
            normalized = model(batch).view(-1)
        raw = normalized.float() * float(checkpoint["target_std"]) + float(
            checkpoint["target_mean"]
        )
        predictions.extend(raw.cpu().tolist())
        indices.extend(batch.sample_id.view(-1).cpu().tolist())
    order = np.argsort(indices)
    return pd.DataFrame(
        {
            "cif_name": [names[index] for index in order],
            f"predicted_{TARGET_NAME}_{TARGET_UNIT}": np.asarray(predictions)[order],
        }
    )


def log_configuration(
    logger: logging.Logger,
    graph_config: GraphConfig,
    model_config: ModelConfig,
    splits: Dict[str, List[Dict[str, Any]]],
    target_mean: float,
    target_std: float,
) -> None:
    logger.info("Model name: %s", MODEL_NAME)
    logger.info("Run version: %s", RUN_VERSION)
    logger.info("Output root: %s", OUTPUT_ROOT)
    logger.info("Excel path: %s", EXCEL_PATH)
    logger.info("CIF directory: %s", CIF_DIR)
    logger.info("Excel columns: first column=CIF filename, second column=target")
    logger.info("Target: %s (%s)", TARGET_NAME, TARGET_UNIT)
    logger.info("Seed: %d", SEED)
    logger.info("Split ratios: %.3f / %.3f / %.3f", TRAIN_RATIO, VAL_RATIO, TEST_RATIO)
    logger.info(
        "Split sizes: train=%d, val=%d, test=%d",
        len(splits["train"]),
        len(splits["val"]),
        len(splits["test"]),
    )
    logger.info("Graph configuration: %s", asdict(graph_config))
    logger.info("Model configuration: %s", asdict(model_config))
    logger.info("Target normalization: %s", NORMALIZE_TARGET)
    logger.info("Target mean/std: %.8f / %.8f", target_mean, target_std)
    logger.info("Loss: MSELoss")
    logger.info("Optuna enabled: %s", USE_OPTUNA)


def main() -> None:
    total_start = time.perf_counter()
    create_directories()
    logger = setup_logger()
    set_global_seed(SEED)
    device = select_device(logger)
    graph_config = GraphConfig()
    base_model_config = ModelConfig()

    data_start = time.perf_counter()
    raw_records = load_excel_records(logger)
    valid_records, graphs = build_or_load_graphs(raw_records, graph_config, logger)
    splits = split_records(valid_records, logger)
    graph_by_name = {
        record["cif_name"]: graph for record, graph in zip(valid_records, graphs)
    }
    sample_id_by_name = {
        record["cif_name"]: index for index, record in enumerate(valid_records)
    }
    train_targets = np.asarray(
        [record["target"] for record in splits["train"]], dtype=np.float64
    )
    target_mean = float(np.mean(train_targets)) if NORMALIZE_TARGET else 0.0
    target_std = float(np.std(train_targets)) if NORMALIZE_TARGET else 1.0
    if target_std < 1.0e-12:
        raise ValueError("Training targets have near-zero standard deviation")
    datasets = {
        split_name: PotNetDataset(
            split_records_list,
            graph_by_name,
            sample_id_by_name,
            target_mean,
            target_std,
        )
        for split_name, split_records_list in splits.items()
    }
    data_time = time.perf_counter() - data_start
    log_configuration(
        logger, graph_config, base_model_config, splits, target_mean, target_std
    )
    logger.info("Data preparation time: %s", format_duration(data_time))

    if USE_OPTUNA:
        model_config, training_config = run_optuna(
            base_model_config,
            graph_config,
            datasets,
            device,
            target_mean,
            target_std,
            logger,
        )
    else:
        model_config = base_model_config
        training_config = {
            "batch_size": BATCH_SIZE,
            "learning_rate": LEARNING_RATE,
            "weight_decay": WEIGHT_DECAY,
            "max_epochs": MAX_EPOCHS,
        }

    preview_model = PotNet(model_config)
    total_parameters, trainable_parameters = count_parameters(preview_model)
    logger.info("Model architecture:\n%s", repr(preview_model))
    logger.info("Total parameters: %d", total_parameters)
    logger.info("Trainable parameters: %d", trainable_parameters)
    logger.info("Optimizer: AdamW, parameters=%s", training_config)
    logger.info(
        "Scheduler: ReduceLROnPlateau(factor=%s, patience=%s, min_lr=%s)",
        LR_FACTOR,
        LR_PATIENCE,
        MIN_LR,
    )
    del preview_model

    training_start = time.perf_counter()
    _, history, best_epoch, best_val_rmse, early_stopped = train_model(
        model_config,
        graph_config,
        datasets,
        device,
        target_mean,
        target_std,
        logger,
        training_config,
        CHECKPOINT_PATH,
    )
    training_time = time.perf_counter() - training_start
    logger.info("Best epoch: %d", best_epoch)
    logger.info("Best validation RMSE: %.6f %s", best_val_rmse, TARGET_UNIT)
    logger.info("Early stopping triggered: %s", early_stopped)
    logger.info("Training time: %s", format_duration(training_time))
    logger.info("Best checkpoint: %s", CHECKPOINT_PATH)

    evaluation_start = time.perf_counter()
    best_model, checkpoint, device = load_trained_model(CHECKPOINT_PATH, device)
    loaders = make_loaders(datasets, int(training_config["batch_size"]), device)
    outputs, metrics = evaluate_splits(
        best_model, loaders, device, target_mean, target_std, valid_records
    )
    save_metrics_table(metrics, outputs)
    for split_name in ("train", "val", "test"):
        save_prediction_dat(split_name, outputs[split_name])
    save_all_prediction_dat(outputs)
    limits = parity_limits(outputs.values())
    for index, split_name in enumerate(("train", "val", "test")):
        plot_parity(
            split_name, outputs[split_name], metrics[split_name], limits, COLORS[index]
        )
    plot_combined_parity(outputs, metrics, limits)
    save_and_plot_history(history, int(checkpoint["best_epoch"]))
    evaluation_time = time.perf_counter() - evaluation_start

    for split_name in ("train", "val", "test"):
        metric = metrics[split_name]
        logger.info(
            "%s metrics | MAE %.6f %s | RMSE %.6f %s | R2 %.6f",
            split_name.capitalize(),
            metric["mae"],
            TARGET_UNIT,
            metric["rmse"],
            TARGET_UNIT,
            metric["r2"],
        )
    logger.info("Evaluation and plotting time: %s", format_duration(evaluation_time))
    logger.info("Total runtime: %s", format_duration(time.perf_counter() - total_start))


if __name__ == "__main__":
    main()

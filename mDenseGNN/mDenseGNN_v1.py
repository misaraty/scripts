from __future__ import annotations

import contextlib
import hashlib
import json
import logging
import math
import os
import platform
import random
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple, Union

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from pymatgen.core import Element, Structure
from scipy.spatial import ConvexHull, Voronoi
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from torch import Tensor, nn
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch_geometric.utils import scatter
from tqdm.auto import tqdm


# User configuration

MODEL_NAME = "DenseGNN"
RUN_VERSION = "v1"

EXCEL_PATH = "./data.xlsx"
CIF_DIR = "./cif"

TARGET_NAME = "bandgap"
TARGET_UNIT = "eV"

SEED = 42
TRAIN_RATIO = 0.80
VAL_RATIO = 0.10
TEST_RATIO = 0.10

# Graph and model settings
VORONOI_MIN_RIDGE_AREA = 0.10
MAX_ATOMIC_NUMBER = 118
CACHE_GRAPHS = True

HIDDEN_DIM = 128
DEPTH = 5
ATOM_EMBEDDING_DIM = 128
DISTANCE_BINS = 32
DISTANCE_MAX = 8.0
RIDGE_AREA_BINS = 25
RIDGE_AREA_MAX = 32.0
EDGE_POOL = "sum"
GRAPH_POOL = "mean"

# Training settings
BATCH_SIZE = 32
MAX_EPOCHS = 300
LEARNING_RATE = 1.0e-3
WEIGHT_DECAY = 1.0e-5
PATIENCE = 50
MIN_DELTA = 1.0e-5
NUM_WORKERS = 0
STANDARDIZE_TARGET = True
USE_AMP = True
GRADIENT_CLIP_NORM = 5.0
LR_FACTOR = 0.5
LR_PATIENCE = 15
MIN_LEARNING_RATE = 1.0e-6
DETERMINISTIC_ALGORITHMS = True

# Optional Optuna search
USE_OPTUNA = False
OPTUNA_N_TRIALS = 30
OPTUNA_TIMEOUT: Optional[int] = None
OPTUNA_EPOCHS = 100
OPTUNA_PATIENCE = 20

# Plot settings
FIG_DPI = 600
FIG_SIZE = (6.2, 5.4)
TITLE_FONTSIZE = 14
LABEL_FONTSIZE = 12
TICK_FONTSIZE = 10
LEGEND_FONTSIZE = 10
ANNOTATION_FONTSIZE = 9
MARKER_SIZE = 24
USE_GRID = True
FONT_FAMILY = ["Arial", "DejaVu Sans"]
COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:red"]

# Output paths
OUTPUT_ROOT = f"./{MODEL_NAME}_{RUN_VERSION}"
FIGURE_DIR = f"{OUTPUT_ROOT}/figure"
DAT_DIR = f"{OUTPUT_ROOT}/dat"
TABLE_DIR = f"{OUTPUT_ROOT}/table"
LOG_DIR = f"{OUTPUT_ROOT}/log"
SPLIT_DIR = f"{OUTPUT_ROOT}/split"
CACHE_DIR = f"{OUTPUT_ROOT}/cache/DenseGNN_graphs"
CHECKPOINT_PATH = f"{OUTPUT_ROOT}/{MODEL_NAME}_best.pt"
SPLIT_FILE = f"{SPLIT_DIR}/{MODEL_NAME}_split.csv"


@dataclass(frozen=True)
class GraphConfig:

    min_ridge_area: float = VORONOI_MIN_RIDGE_AREA
    max_atomic_number: int = MAX_ATOMIC_NUMBER


@dataclass(frozen=True)
class ModelConfig:

    hidden_dim: int = HIDDEN_DIM
    depth: int = DEPTH
    atom_embedding_dim: int = ATOM_EMBEDDING_DIM
    distance_bins: int = DISTANCE_BINS
    distance_max: float = DISTANCE_MAX
    ridge_area_bins: int = RIDGE_AREA_BINS
    ridge_area_max: float = RIDGE_AREA_MAX
    edge_pool: str = EDGE_POOL
    graph_pool: str = GRAPH_POOL


@dataclass
class TargetNormalizer:

    mean: float = 0.0
    std: float = 1.0
    enabled: bool = True

    @classmethod
    def fit(cls, values: Sequence[float], enabled: bool = True) -> "TargetNormalizer":
        arr = np.asarray(values, dtype=np.float64)
        if not enabled:
            return cls(0.0, 1.0, False)
        std = float(np.std(arr, ddof=0))
        if not np.isfinite(std) or std < 1.0e-12:
            raise ValueError(
                "The training-target standard deviation is too small for normalization."
            )
        return cls(float(np.mean(arr)), std, True)

    def normalize_tensor(self, value: Tensor) -> Tensor:
        if not self.enabled:
            return value
        return (value - self.mean) / self.std

    def denormalize_tensor(self, value: Tensor) -> Tensor:
        if not self.enabled:
            return value
        return value * self.std + self.mean

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)


class CrystalGraphData(Data):
    pass


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
    if DETERMINISTIC_ALGORITHMS:
        try:
            torch.use_deterministic_algorithms(True, warn_only=True)
        except TypeError:
            torch.use_deterministic_algorithms(True)


def seed_worker(worker_id: int) -> None:

    del worker_id
    worker_seed = torch.initial_seed() % (2**32)
    np.random.seed(worker_seed)
    random.seed(worker_seed)


def create_directories() -> None:

    for path in [FIGURE_DIR, DAT_DIR, TABLE_DIR, LOG_DIR, SPLIT_DIR]:
        Path(path).mkdir(parents=True, exist_ok=True)
    if CACHE_GRAPHS:
        Path(CACHE_DIR).mkdir(parents=True, exist_ok=True)


def setup_logger() -> logging.Logger:

    logger = logging.getLogger(MODEL_NAME)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    logger.handlers.clear()
    formatter = logging.Formatter(
        fmt="%(asctime)s | %(levelname)s | %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
    file_handler = logging.FileHandler(
        Path(LOG_DIR) / f"{MODEL_NAME}_training.log", mode="w", encoding="utf-8"
    )
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def format_duration(seconds: float) -> str:

    total = int(round(seconds))
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def load_excel_records(logger: logging.Logger) -> List[Dict[str, Any]]:

    excel_path = Path(EXCEL_PATH)
    if not excel_path.is_file():
        raise FileNotFoundError(f"Excel file not found: {excel_path.resolve()}")
    frame = pd.read_excel(excel_path)
    if frame.shape[1] < 2:
        raise ValueError(
            "The Excel file must contain at least two columns: CIF filename and bandgap."
        )
    cif_col = frame.columns[0]
    target_col = frame.columns[1]

    cif_root = Path(CIF_DIR)
    records: List[Dict[str, Any]] = []
    seen: Dict[str, float] = {}
    invalid = 0
    for row_index, row in frame.iterrows():
        raw_name = row[cif_col]
        raw_target = row[target_col]
        if pd.isna(raw_name) or not str(raw_name).strip():
            logger.warning("Skipping Excel row %s: empty CIF filename.", row_index + 2)
            invalid += 1
            continue
        cif_name = f"{int(raw_name)}.cif"
        try:
            target = float(raw_target)
        except (TypeError, ValueError):
            logger.warning(
                "Skipping %s: invalid bandgap value (%r).", cif_name, raw_target
            )
            invalid += 1
            continue
        if not np.isfinite(target):
            logger.warning("Skipping %s: bandgap is NaN or infinite.", cif_name)
            invalid += 1
            continue

        input_path = Path(cif_name)
        cif_path = input_path if input_path.is_absolute() else cif_root / input_path
        if not cif_path.is_file():
            logger.warning("Skipping %s: CIF file not found at %s.", cif_name, cif_path)
            invalid += 1
            continue

        canonical_name = cif_name.replace("\\", "/")
        if canonical_name in seen:
            if not math.isclose(
                seen[canonical_name], target, rel_tol=0.0, abs_tol=1e-10
            ):
                logger.warning(
                    "Skipping duplicate CIF %s: conflicting targets %.8g and %.8g.",
                    canonical_name,
                    seen[canonical_name],
                    target,
                )
            else:
                logger.warning("Skipping duplicate CIF %s.", canonical_name)
            invalid += 1
            continue
        seen[canonical_name] = target
        records.append(
            {
                "cif_name": canonical_name,
                "cif_path": str(cif_path.resolve()),
                "target": target,
            }
        )

    logger.info(
        "Excel loaded: rows=%d, valid candidates=%d, initially invalid=%d, "
        "CIF column=%r, target column=%r.",
        len(frame),
        len(records),
        invalid,
        cif_col,
        target_col,
    )
    if len(records) < 10:
        raise ValueError(
            "Fewer than 10 valid samples remain; the dataset cannot be split reliably."
        )
    return records


def _ridge_area(points: np.ndarray) -> float:

    ridge_points = np.asarray(points, dtype=np.float64)
    if ridge_points.size == 0:
        return 0.0
    while ridge_points.shape[0] <= 3:
        ridge_points = np.append(ridge_points, ridge_points[:1], axis=0)
    return float(ConvexHull(ridge_points, qhull_options="QJ").area / 2.0)


def _supercell_fractional_coordinates(frac: np.ndarray) -> np.ndarray:

    expanded = np.zeros((3, 3, 3, len(frac), 3), dtype=np.float64)
    for ix, dx in enumerate(range(-1, 2)):
        for iy, dy in enumerate(range(-1, 2)):
            for iz, dz in enumerate(range(-1, 2)):
                expanded[ix, iy, iz] = frac + np.array([dx, dy, dz], dtype=np.float64)
    return expanded


def _voronoi_edges(
    positions: np.ndarray, lattice: np.ndarray, min_ridge_area: float
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:

    fractional = positions @ np.linalg.inv(lattice)
    expanded_cart = _supercell_fractional_coordinates(fractional) @ lattice
    flat_cart = expanded_cart.reshape(-1, 3)
    if flat_cart.shape[0] < 5:
        raise ValueError(
            "Too few periodic points for a three-dimensional Voronoi graph."
        )
    voronoi = Voronoi(flat_cart)
    unraveled = np.array(
        np.unravel_index(voronoi.ridge_points, expanded_cart.shape[:-1])
    )
    unraveled = np.moveaxis(unraveled, np.arange(3), np.roll(np.arange(3), 1))
    source_center = np.flatnonzero(np.all(unraveled[:, 0, :3] == 1, axis=-1))
    target_center = np.flatnonzero(np.all(unraveled[:, 1, :3] == 1, axis=-1))
    if source_center.size + target_center.size == 0:
        raise ValueError("No center-cell Voronoi ridges were generated.")
    edge_info = np.vstack(
        [unraveled[source_center][:, [1, 0]], unraveled[target_center]]
    )
    ridge_ids = np.concatenate([source_center, target_center])
    edge_pairs = edge_info[:, :, -1].astype(np.int64)
    distances = np.asarray(
        [
            np.linalg.norm(
                expanded_cart[tuple(info[0])] - expanded_cart[tuple(info[1])]
            )
            for info in edge_info
        ],
        dtype=np.float64,
    )
    areas: List[float] = []
    for ridge_id in ridge_ids:
        vertex_ids = voronoi.ridge_vertices[int(ridge_id)]
        if not vertex_ids or any(index < 0 for index in vertex_ids):
            areas.append(0.0)
            continue
        try:
            areas.append(_ridge_area(voronoi.vertices[vertex_ids]))
        except Exception:
            areas.append(0.0)
    ridge_areas = np.asarray(areas, dtype=np.float64)
    valid = (
        np.isfinite(distances)
        & np.isfinite(ridge_areas)
        & (distances > 1.0e-8)
        & (ridge_areas > float(min_ridge_area))
    )
    edge_pairs = edge_pairs[valid]
    distances = distances[valid]
    ridge_areas = ridge_areas[valid]
    if edge_pairs.shape[0] == 0:
        raise ValueError("No Voronoi edges remain after ridge-area filtering.")
    return edge_pairs[:, [1, 0]].T, distances, ridge_areas


def build_graph_tensors(
    cif_path: Union[str, Path], config: GraphConfig
) -> Dict[str, Tensor]:

    structure = Structure.from_file(str(cif_path), primitive=False, sort=False)
    if len(structure) == 0:
        raise ValueError("The CIF structure contains no atoms.")
    if not structure.is_ordered:
        raise ValueError(
            "Disordered or partially occupied structures are not supported."
        )
    atomic_numbers = np.asarray(
        [int(site.specie.Z) for site in structure], dtype=np.int64
    )
    if np.any(atomic_numbers < 1) or np.any(atomic_numbers > config.max_atomic_number):
        raise ValueError(
            f"Atomic numbers outside the supported range 1-{config.max_atomic_number}: "
            f"{np.unique(atomic_numbers).tolist()}"
        )
    positions = np.asarray(structure.cart_coords, dtype=np.float64)
    lattice = np.asarray(structure.lattice.matrix, dtype=np.float64)
    edge_index, distances, ridge_areas = _voronoi_edges(
        positions, lattice, config.min_ridge_area
    )
    edge_attr = np.stack([distances, ridge_areas], axis=-1)
    return {
        "z": torch.as_tensor(atomic_numbers, dtype=torch.long),
        "edge_index": torch.as_tensor(edge_index, dtype=torch.long),
        "edge_attr": torch.as_tensor(edge_attr, dtype=torch.float32),
        "num_nodes": torch.tensor(len(structure), dtype=torch.long),
    }


def graph_cache_path(cif_path: Union[str, Path], config: GraphConfig) -> Path:

    path = Path(cif_path).resolve()
    stat = path.stat()
    payload = {
        "path": str(path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "graph_config": asdict(config),
        "cache_version": 2,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, ensure_ascii=False).encode("utf-8")
    ).hexdigest()[:24]
    return Path(CACHE_DIR) / f"{digest}.pt"


def safe_torch_load(path: Union[str, Path], map_location: Any = "cpu") -> Any:

    try:
        return torch.load(path, map_location=map_location, weights_only=True)
    except TypeError:
        return torch.load(path, map_location=map_location)


def load_or_build_graph(
    cif_path: Union[str, Path], config: GraphConfig
) -> Dict[str, Tensor]:

    if not CACHE_GRAPHS:
        return build_graph_tensors(cif_path, config)
    cache_path = graph_cache_path(cif_path, config)
    if cache_path.is_file():
        try:
            cached = safe_torch_load(cache_path, map_location="cpu")
            required = {"z", "edge_index", "edge_attr", "num_nodes"}
            if required.issubset(cached):
                return cached
        except Exception:
            pass
    graph = build_graph_tensors(cif_path, config)
    temporary = cache_path.with_suffix(".tmp")
    torch.save(graph, temporary)
    os.replace(temporary, cache_path)
    return graph


def validate_and_cache_records(
    records: Sequence[Dict[str, Any]], config: GraphConfig, logger: logging.Logger
) -> List[Dict[str, Any]]:

    valid_records: List[Dict[str, Any]] = []
    for record in tqdm(
        records, desc="Validating CIFs / building DenseGNN graphs", unit="cif"
    ):
        try:
            graph = load_or_build_graph(record["cif_path"], config)
            if graph["edge_index"].numel() == 0:
                raise ValueError("The crystal graph is empty.")
            item = dict(record)
            item["cache_path"] = str(graph_cache_path(record["cif_path"], config))
            valid_records.append(item)
        except Exception as exc:
            logger.warning(
                "Skipping %s: CIF parsing or graph construction failed: %s",
                record["cif_name"],
                exc,
            )
    logger.info(
        "CIF graph construction complete: valid=%d, failed=%d.",
        len(valid_records),
        len(records) - len(valid_records),
    )
    if len(valid_records) < 10:
        raise ValueError(
            "Fewer than 10 samples produced valid graphs; training cannot continue."
        )
    return valid_records


def create_or_load_split(
    records: Sequence[Dict[str, Any]], logger: logging.Logger
) -> Dict[str, List[Dict[str, Any]]]:

    ratio_sum = TRAIN_RATIO + VAL_RATIO + TEST_RATIO
    if not math.isclose(ratio_sum, 1.0, rel_tol=0.0, abs_tol=1e-10):
        raise ValueError(
            f"Train/validation/test ratios must sum to 1; received {ratio_sum}."
        )
    split_path = Path(SPLIT_FILE)
    record_by_name = {item["cif_name"]: item for item in records}
    if split_path.is_file():
        saved = pd.read_csv(split_path)
        required = {"cif_name", "target", "split"}
        if not required.issubset(saved.columns):
            raise ValueError(
                f"Existing split file {split_path} is missing columns: "
                f"{required - set(saved.columns)}"
            )
        saved_names = set(saved["cif_name"].astype(str))
        if saved_names != set(record_by_name):
            missing = sorted(set(record_by_name) - saved_names)[:5]
            extra = sorted(saved_names - set(record_by_name))[:5]
            raise ValueError(
                "The existing split file does not match the current valid dataset. "
                "Check the data or move the old split file before rerunning. "
                f"Missing examples={missing}; extra examples={extra}"
            )
        splits: Dict[str, List[Dict[str, Any]]] = {"train": [], "val": [], "test": []}
        for row in saved.itertuples(index=False):
            split_name = str(row.split)
            if split_name not in splits:
                raise ValueError(
                    f"Unknown split value in the split file: {split_name!r}."
                )
            current = record_by_name[str(row.cif_name)]
            if not math.isclose(
                float(row.target), current["target"], rel_tol=0.0, abs_tol=1e-8
            ):
                raise ValueError(
                    f"The current target for {row.cif_name} differs from the split file."
                )
            splits[split_name].append(current)
        if any(len(items) == 0 for items in splits.values()):
            raise ValueError(
                "The existing split file contains an empty train, validation, or test set."
            )
        logger.info("Reusing existing data split: %s", split_path)
        return splits

    indices = np.arange(len(records))
    train_indices, holdout_indices = train_test_split(
        indices,
        test_size=VAL_RATIO + TEST_RATIO,
        random_state=SEED,
        shuffle=True,
    )
    relative_test = TEST_RATIO / (VAL_RATIO + TEST_RATIO)
    val_indices, test_indices = train_test_split(
        holdout_indices,
        test_size=relative_test,
        random_state=SEED,
        shuffle=True,
    )
    index_groups = {"train": train_indices, "val": val_indices, "test": test_indices}
    splits = {
        name: [records[int(index)] for index in values]
        for name, values in index_groups.items()
    }
    rows: List[Dict[str, Any]] = []
    for split_name in ["train", "val", "test"]:
        for item in splits[split_name]:
            rows.append(
                {
                    "cif_name": item["cif_name"],
                    "target": item["target"],
                    "split": split_name,
                }
            )
    pd.DataFrame(rows).to_csv(split_path, index=False, encoding="utf-8-sig")
    logger.info("Created fixed data split: %s", split_path)
    return splits


class BandgapDataset(torch.utils.data.Dataset):

    def __init__(
        self,
        records: Sequence[Dict[str, Any]],
        graph_config: GraphConfig,
        normalizer: TargetNormalizer,
    ) -> None:
        self.records = list(records)
        self.graph_config = graph_config
        self.normalizer = normalizer

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, index: int) -> CrystalGraphData:
        item = self.records[index]
        graph = load_or_build_graph(item["cif_path"], self.graph_config)
        raw_target = torch.tensor([item["target"]], dtype=torch.float32)
        target = self.normalizer.normalize_tensor(raw_target)
        return CrystalGraphData(
            z=graph["z"],
            edge_index=graph["edge_index"],
            edge_attr=graph["edge_attr"],
            y=target,
            cif_name=item["cif_name"],
            num_nodes=int(graph["num_nodes"].item()),
        )


def make_loaders(
    splits: Dict[str, List[Dict[str, Any]]],
    graph_config: GraphConfig,
    normalizer: TargetNormalizer,
    batch_size: int,
) -> Dict[str, DataLoader]:

    generator = torch.Generator()
    generator.manual_seed(SEED)
    common = {
        "batch_size": batch_size,
        "num_workers": NUM_WORKERS,
        "pin_memory": torch.cuda.is_available(),
        "worker_init_fn": seed_worker,
        "persistent_workers": NUM_WORKERS > 0,
    }
    loaders: Dict[str, DataLoader] = {}
    for split_name in ["train", "val", "test"]:
        dataset = BandgapDataset(splits[split_name], graph_config, normalizer)
        loaders[split_name] = DataLoader(
            dataset,
            shuffle=split_name == "train",
            generator=generator if split_name == "train" else None,
            **common,
        )
    return loaders


def _numeric(value: Any) -> float:

    try:
        return float(value)
    except (TypeError, ValueError, OverflowError):
        return float("nan")


def _standardize_element_property(values: np.ndarray) -> Tensor:

    active = np.asarray(values[:MAX_ATOMIC_NUMBER], dtype=np.float64)
    finite = np.isfinite(active)
    fill = float(np.median(active[finite])) if np.any(finite) else 0.0
    active = np.where(finite, active, fill)
    scale = float(np.std(active, ddof=1))
    if not np.isfinite(scale) or scale < 1.0e-12:
        scale = 1.0
    standardized = (active - float(np.mean(active))) / scale
    output = np.zeros(MAX_ATOMIC_NUMBER + 1, dtype=np.float32)
    output[:MAX_ATOMIC_NUMBER] = standardized.astype(np.float32)
    return torch.from_numpy(output)


_PERIODIC_PROPERTIES_CACHE: Optional[Dict[str, Tensor]] = None


def _periodic_properties() -> Dict[str, Tensor]:

    global _PERIODIC_PROPERTIES_CACHE
    if _PERIODIC_PROPERTIES_CACHE is not None:
        return _PERIODIC_PROPERTIES_CACHE
    raw = {
        "mass": np.full(MAX_ATOMIC_NUMBER + 1, np.nan, dtype=np.float64),
        "radius": np.full(MAX_ATOMIC_NUMBER + 1, np.nan, dtype=np.float64),
        "en": np.full(MAX_ATOMIC_NUMBER + 1, np.nan, dtype=np.float64),
        "ie": np.full(MAX_ATOMIC_NUMBER + 1, np.nan, dtype=np.float64),
        "mp": np.full(MAX_ATOMIC_NUMBER + 1, np.nan, dtype=np.float64),
        "density": np.full(MAX_ATOMIC_NUMBER + 1, np.nan, dtype=np.float64),
    }
    oxidation = np.zeros((MAX_ATOMIC_NUMBER + 1, 14), dtype=np.float32)
    for atomic_number in range(1, MAX_ATOMIC_NUMBER + 1):
        element = Element.from_Z(atomic_number)
        index = atomic_number - 1
        raw["mass"][index] = _numeric(element.atomic_mass)
        raw["radius"][index] = _numeric(element.atomic_radius)
        raw["en"][index] = _numeric(element.X)
        ionization = getattr(element, "ionization_energies", ())
        raw["ie"][index] = _numeric(ionization[0] if ionization else None)
        raw["mp"][index] = _numeric(getattr(element, "melting_point", None))
        raw["density"][index] = _numeric(getattr(element, "density_of_solid", None))
        for state in getattr(element, "oxidation_states", ()):
            state_index = int(round(float(state))) + 7
            if 0 <= state_index < oxidation.shape[1]:
                oxidation[index, state_index] = 1.0
    properties = {
        key: _standardize_element_property(values) for key, values in raw.items()
    }
    properties["ox"] = torch.from_numpy(oxidation)
    _PERIODIC_PROPERTIES_CACHE = properties
    return _PERIODIC_PROPERTIES_CACHE


def _initialize_linear_layers(module: nn.Module) -> None:

    for layer in module.modules():
        if isinstance(layer, nn.Linear):
            nn.init.xavier_uniform_(layer.weight)
            if layer.bias is not None:
                nn.init.zeros_(layer.bias)


def _swish_mlp(sizes: Sequence[int]) -> nn.Sequential:

    layers: List[nn.Module] = []
    for input_size, output_size in zip(sizes[:-1], sizes[1:]):
        layers.extend([nn.Linear(input_size, output_size), nn.SiLU()])
    return nn.Sequential(*layers)


class GaussianBasis(nn.Module):

    def __init__(self, bins: int, max_value: float) -> None:
        super().__init__()
        offsets = torch.linspace(0.0, max_value, bins + 1)[1:]
        sigma = math.sqrt(max_value / bins)
        self.register_buffer("offsets", offsets)
        self.register_buffer("sigma_inv2", torch.tensor(1.0 / (2.0 * sigma * sigma)))

    def forward(self, values: Tensor) -> Tensor:
        values = values.reshape(-1)
        return torch.exp(
            -((values.unsqueeze(-1) - self.offsets) ** 2) * self.sigma_inv2
        )


class DenseGNNAtomEmbedding(nn.Module):

    def __init__(self, hidden_dim: int, embedding_dim: int) -> None:
        super().__init__()
        self.embedding = nn.Embedding(MAX_ATOMIC_NUMBER + 1, embedding_dim)
        nn.init.uniform_(self.embedding.weight, -0.05, 0.05)
        for name, values in _periodic_properties().items():
            self.register_buffer(f"pt_{name}", values.float())
        self.projection = nn.Linear(embedding_dim + 6 + 14, hidden_dim)

    def forward(self, atomic_numbers: Tensor) -> Tensor:
        indices = (atomic_numbers.long() - 1).clamp(0, MAX_ATOMIC_NUMBER)
        features = [
            self.embedding(indices),
            self.pt_mass[indices].unsqueeze(-1),
            self.pt_radius[indices].unsqueeze(-1),
            self.pt_en[indices].unsqueeze(-1),
            self.pt_ie[indices].unsqueeze(-1),
            self.pt_ox[indices],
            self.pt_mp[indices].unsqueeze(-1),
            self.pt_density[indices].unsqueeze(-1),
        ]
        return self.projection(torch.cat(features, dim=-1))


class DenseGNNEdgeEmbedding(nn.Module):

    def __init__(self, config: ModelConfig) -> None:
        super().__init__()
        self.distance_basis = GaussianBasis(config.distance_bins, config.distance_max)
        self.area_basis = GaussianBasis(config.ridge_area_bins, config.ridge_area_max)
        self.projection = nn.Linear(
            config.distance_bins + config.ridge_area_bins, config.hidden_dim
        )

    def forward(self, edge_attr: Tensor) -> Tensor:
        distance = self.distance_basis(edge_attr[:, 0])
        ridge_area = self.area_basis(edge_attr[:, 1])
        return self.projection(torch.cat([distance, ridge_area], dim=-1))


class DenseGNNConv(nn.Module):

    def __init__(self, hidden_dim: int, edge_pool: str, graph_pool: str) -> None:
        super().__init__()
        self.edge_pool = edge_pool
        self.graph_pool = graph_pool
        self.edge_mlp = _swish_mlp([4 * hidden_dim, hidden_dim, hidden_dim, hidden_dim])
        self.node_mlp = _swish_mlp([3 * hidden_dim, hidden_dim])
        self.graph_mlp = _swish_mlp([3 * hidden_dim, hidden_dim])

    def forward(
        self,
        node_features: Tensor,
        edge_features: Tensor,
        global_features: Tensor,
        edge_index: Tensor,
        node_batch: Tensor,
        edge_batch: Tensor,
        num_graphs: int,
    ) -> Tuple[Tensor, Tensor, Tensor]:
        target, source = edge_index[0], edge_index[1]
        edge_input = torch.cat(
            [
                node_features[target],
                node_features[source],
                edge_features,
                global_features[edge_batch],
            ],
            dim=-1,
        )
        updated_edges = F.silu(self.edge_mlp(edge_input) + edge_features)
        aggregated_edges = scatter(
            updated_edges,
            target,
            dim=0,
            dim_size=node_features.size(0),
            reduce=self.edge_pool,
        )
        node_input = torch.cat(
            [node_features, aggregated_edges, global_features[node_batch]], dim=-1
        )
        updated_nodes = F.silu(self.node_mlp(node_input) + node_features)
        pooled_edges = scatter(
            updated_edges,
            edge_batch,
            dim=0,
            dim_size=num_graphs,
            reduce=self.graph_pool,
        )
        pooled_nodes = scatter(
            updated_nodes,
            node_batch,
            dim=0,
            dim_size=num_graphs,
            reduce=self.graph_pool,
        )
        global_input = torch.cat([pooled_edges, pooled_nodes, global_features], dim=-1)
        updated_globals = F.silu(self.graph_mlp(global_input) + global_features)
        return updated_nodes, updated_edges, updated_globals


class DenseGNN(nn.Module):
    """DenseGNN crystal regressor with node, edge, and global states."""

    def __init__(self, config: ModelConfig) -> None:
        super().__init__()
        if config.depth < 1:
            raise ValueError("DenseGNN depth must be at least one.")
        if config.edge_pool not in {"sum", "mean", "max"}:
            raise ValueError(f"Unsupported edge pooling: {config.edge_pool}")
        if config.graph_pool not in {"sum", "mean", "max"}:
            raise ValueError(f"Unsupported graph pooling: {config.graph_pool}")
        self.config = config
        self.atom_embedding = DenseGNNAtomEmbedding(
            config.hidden_dim, config.atom_embedding_dim
        )
        self.edge_embedding = DenseGNNEdgeEmbedding(config)
        self.global_embedding = nn.Sequential(
            nn.Linear(1, config.hidden_dim), nn.ReLU()
        )
        self.convolutions = nn.ModuleList(
            [
                DenseGNNConv(config.hidden_dim, config.edge_pool, config.graph_pool)
                for _ in range(config.depth)
            ]
        )
        self.node_transitions = nn.ModuleList(
            [
                nn.Sequential(
                    nn.Linear(config.hidden_dim * (index + 1), config.hidden_dim),
                    nn.SiLU(),
                )
                for index in range(1, config.depth)
            ]
        )
        self.edge_transitions = nn.ModuleList(
            [
                nn.Sequential(
                    nn.Linear(config.hidden_dim * (index + 1), config.hidden_dim),
                    nn.SiLU(),
                )
                for index in range(1, config.depth)
            ]
        )
        self.global_transitions = nn.ModuleList(
            [
                nn.Sequential(
                    nn.Linear(config.hidden_dim * (index + 1), config.hidden_dim),
                    nn.SiLU(),
                )
                for index in range(1, config.depth)
            ]
        )
        self.output = nn.Linear(3 * config.hidden_dim, 1)
        _initialize_linear_layers(self)

    def forward(self, data: CrystalGraphData) -> Tensor:
        node_batch = data.batch
        edge_index = data.edge_index
        edge_batch = node_batch[edge_index[0]]
        num_graphs = int(node_batch.max().item()) + 1
        node_features = self.atom_embedding(data.z)
        edge_features = self.edge_embedding(data.edge_attr)
        charge = torch.zeros(
            (num_graphs, 1), device=node_features.device, dtype=node_features.dtype
        )
        global_features = self.global_embedding(charge)
        node_history = [node_features]
        edge_history = [edge_features]
        global_history = [global_features]
        for index, convolution in enumerate(self.convolutions):
            if index > 0:
                node_features = self.node_transitions[index - 1](
                    torch.cat(node_history, dim=-1)
                )
                edge_features = self.edge_transitions[index - 1](
                    torch.cat(edge_history, dim=-1)
                )
                global_features = self.global_transitions[index - 1](
                    torch.cat(global_history, dim=-1)
                )
            node_features, edge_features, global_features = convolution(
                node_features,
                edge_features,
                global_features,
                edge_index,
                node_batch,
                edge_batch,
                num_graphs,
            )
            node_history.append(node_features)
            edge_history.append(edge_features)
            global_history.append(global_features)
        pooled_edges = scatter(
            edge_features, edge_batch, dim=0, dim_size=num_graphs, reduce="mean"
        )
        pooled_nodes = scatter(
            node_features, node_batch, dim=0, dim_size=num_graphs, reduce="mean"
        )
        graph_features = torch.cat(
            [pooled_edges, pooled_nodes, global_features], dim=-1
        )
        return self.output(graph_features).view(-1)


def calculate_metrics(
    y_true: Sequence[float], y_pred: Sequence[float]
) -> Dict[str, float]:

    true = np.asarray(y_true, dtype=np.float64).reshape(-1)
    pred = np.asarray(y_pred, dtype=np.float64).reshape(-1)
    if true.size != pred.size or true.size == 0:
        raise ValueError("Targets and predictions must have the same nonzero length.")
    r2 = float(r2_score(true, pred)) if true.size >= 2 else float("nan")
    return {
        "mae": float(mean_absolute_error(true, pred)),
        "rmse": float(math.sqrt(mean_squared_error(true, pred))),
        "r2": r2,
    }


def amp_context(device: torch.device):

    if USE_AMP and device.type == "cuda":
        return torch.autocast(device_type="cuda", dtype=torch.bfloat16)
    return contextlib.nullcontext()


def train_one_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    criterion: nn.Module,
    normalizer: TargetNormalizer,
    device: torch.device,
) -> Tuple[float, float]:

    model.train()
    squared_error_sum = 0.0
    sample_count = 0
    raw_true: List[np.ndarray] = []
    raw_pred: List[np.ndarray] = []
    for batch in loader:
        batch = batch.to(device, non_blocking=True)
        optimizer.zero_grad(set_to_none=True)
        with amp_context(device):
            prediction = model(batch).view(-1)
            target = batch.y.view(-1)
            loss = criterion(prediction, target)
        if not torch.isfinite(loss):
            raise FloatingPointError(
                f"Training loss became NaN or infinite: {float(loss.detach().cpu())}"
            )
        loss.backward()
        if GRADIENT_CLIP_NORM is not None and GRADIENT_CLIP_NORM > 0:
            nn.utils.clip_grad_norm_(model.parameters(), GRADIENT_CLIP_NORM)
        optimizer.step()

        count = int(target.numel())
        squared_error_sum += float(loss.detach().cpu()) * count
        sample_count += count
        raw_true.append(
            normalizer.denormalize_tensor(target.detach()).float().cpu().numpy()
        )
        raw_pred.append(
            normalizer.denormalize_tensor(prediction.detach()).float().cpu().numpy()
        )
    train_rmse = calculate_metrics(np.concatenate(raw_true), np.concatenate(raw_pred))[
        "rmse"
    ]
    return squared_error_sum / max(sample_count, 1), train_rmse


@torch.inference_mode()
def evaluate_loader(
    model: nn.Module,
    loader: DataLoader,
    criterion: nn.Module,
    normalizer: TargetNormalizer,
    device: torch.device,
) -> Dict[str, Any]:

    model.eval()
    squared_error_sum = 0.0
    sample_count = 0
    all_true: List[np.ndarray] = []
    all_pred: List[np.ndarray] = []
    all_names: List[str] = []
    for batch in loader:
        names = list(batch.cif_name)
        batch = batch.to(device, non_blocking=True)
        with amp_context(device):
            prediction = model(batch).view(-1)
            target = batch.y.view(-1)
            loss = criterion(prediction, target)
        count = int(target.numel())
        squared_error_sum += float(loss.detach().cpu()) * count
        sample_count += count
        all_true.append(normalizer.denormalize_tensor(target).float().cpu().numpy())
        all_pred.append(normalizer.denormalize_tensor(prediction).float().cpu().numpy())
        all_names.extend(str(name) for name in names)
    true = np.concatenate(all_true)
    pred = np.concatenate(all_pred)
    return {
        "cif_name": all_names,
        "y_true": true,
        "y_pred": pred,
        "mse_loss": squared_error_sum / max(sample_count, 1),
        "metrics": calculate_metrics(true, pred),
    }


def cpu_state_dict(model: nn.Module) -> Dict[str, Tensor]:

    return {key: value.detach().cpu() for key, value in model.state_dict().items()}


def save_checkpoint(
    path: Union[str, Path],
    model: nn.Module,
    model_config: ModelConfig,
    graph_config: GraphConfig,
    normalizer: TargetNormalizer,
    best_epoch: int,
    best_val_rmse: float,
    training_parameters: Dict[str, Any],
) -> None:

    checkpoint = {
        "format_version": 1,
        "model_name": MODEL_NAME,
        "run_version": RUN_VERSION,
        "model_state_dict": cpu_state_dict(model),
        "model_config": asdict(model_config),
        "graph_config": asdict(graph_config),
        "normalization": normalizer.to_dict(),
        "best_epoch": int(best_epoch),
        "best_val_rmse": float(best_val_rmse),
        "seed": SEED,
        "target_name": TARGET_NAME,
        "target_unit": TARGET_UNIT,
        "training_parameters": training_parameters,
        "torch_version": str(torch.__version__),
    }
    torch.save(checkpoint, path)


def load_trained_model(
    checkpoint_path: Union[str, Path] = CHECKPOINT_PATH,
    device: Optional[torch.device] = None,
) -> Tuple[DenseGNN, GraphConfig, TargetNormalizer, Dict[str, Any]]:
    """Load a trained model and its preprocessing metadata."""

    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint = safe_torch_load(checkpoint_path, map_location=device)
    model_config = ModelConfig(**checkpoint["model_config"])
    graph_config = GraphConfig(**checkpoint["graph_config"])
    normalizer = TargetNormalizer(**checkpoint["normalization"])
    model = DenseGNN(model_config).to(device)
    model.load_state_dict(checkpoint["model_state_dict"], strict=True)
    model.eval()
    return model, graph_config, normalizer, checkpoint


def fit_model(
    loaders: Dict[str, DataLoader],
    model_config: ModelConfig,
    graph_config: GraphConfig,
    normalizer: TargetNormalizer,
    device: torch.device,
    learning_rate: float,
    weight_decay: float,
    max_epochs: int,
    patience: int,
    checkpoint_path: Optional[Union[str, Path]],
    logger: Optional[logging.Logger],
    trial: Any = None,
) -> Tuple[float, int, List[Dict[str, float]]]:

    set_global_seed(SEED)
    model = DenseGNN(model_config).to(device)
    criterion = nn.MSELoss()
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
    best_rmse = float("inf")
    best_epoch = 0
    epochs_without_improvement = 0
    history: List[Dict[str, float]] = []
    training_parameters = {
        "loss_function": "MSELoss",
        "optimizer": "AdamW",
        "scheduler": "ReduceLROnPlateau",
        "learning_rate": learning_rate,
        "weight_decay": weight_decay,
        "batch_size": loaders["train"].batch_size,
        "max_epochs": max_epochs,
        "patience": patience,
        "min_delta": MIN_DELTA,
    }

    for epoch in range(1, max_epochs + 1):
        epoch_start = time.perf_counter()
        current_lr = float(optimizer.param_groups[0]["lr"])
        train_mse, train_rmse = train_one_epoch(
            model, loaders["train"], optimizer, criterion, normalizer, device
        )
        val_result = evaluate_loader(
            model, loaders["val"], criterion, normalizer, device
        )
        val_rmse = float(val_result["metrics"]["rmse"])
        val_mse = float(val_result["mse_loss"])
        scheduler.step(val_rmse)
        history.append(
            {
                "epoch": float(epoch),
                "train_rmse_eV": train_rmse,
                "val_rmse_eV": val_rmse,
                "learning_rate": current_lr,
                "train_mse_loss": train_mse,
                "val_mse_loss": val_mse,
            }
        )
        improved = val_rmse < best_rmse - MIN_DELTA
        if improved:
            best_rmse = val_rmse
            best_epoch = epoch
            epochs_without_improvement = 0
            if checkpoint_path is not None:
                save_checkpoint(
                    checkpoint_path,
                    model,
                    model_config,
                    graph_config,
                    normalizer,
                    best_epoch,
                    best_rmse,
                    training_parameters,
                )
        else:
            epochs_without_improvement += 1

        if logger is not None:
            logger.info(
                "Epoch %04d | train_MSE=%.8f | train_RMSE=%.6f eV | "
                "val_MSE=%.8f | val_RMSE=%.6f eV | lr=%.3e | %.2f s%s",
                epoch,
                train_mse,
                train_rmse,
                val_mse,
                val_rmse,
                current_lr,
                time.perf_counter() - epoch_start,
                " | BEST" if improved else "",
            )
        if trial is not None:
            trial.report(val_rmse, step=epoch)
            if trial.should_prune():
                import optuna

                raise optuna.TrialPruned()
        if epochs_without_improvement >= patience:
            if logger is not None:
                logger.info(
                    "Early stopping: validation RMSE did not improve for %d epochs.",
                    patience,
                )
            break

    if best_epoch == 0:
        raise RuntimeError(
            "Training did not produce a valid best model; check the data and numerical stability."
        )
    del model
    if device.type == "cuda":
        torch.cuda.empty_cache()
    return best_rmse, best_epoch, history


def run_optuna(
    splits: Dict[str, List[Dict[str, Any]]],
    graph_config: GraphConfig,
    normalizer: TargetNormalizer,
    device: torch.device,
    logger: logging.Logger,
) -> Dict[str, Any]:

    try:
        import optuna
    except ImportError as exc:
        raise ImportError("Install Optuna before setting USE_OPTUNA=True.") from exc

    def objective(trial: Any) -> float:
        hidden = trial.suggest_categorical("hidden_dim", [96, 128, 192])
        params = {
            "batch_size": trial.suggest_categorical("batch_size", [16, 32, 64]),
            "learning_rate": trial.suggest_float("learning_rate", 1e-4, 5e-3, log=True),
            "weight_decay": trial.suggest_float("weight_decay", 1e-8, 1e-3, log=True),
            "depth": trial.suggest_int("depth", 3, 6),
            "hidden_dim": hidden,
        }
        trial_config = ModelConfig(
            hidden_dim=params["hidden_dim"],
            depth=params["depth"],
        )
        loaders = make_loaders(
            splits, graph_config, normalizer, batch_size=params["batch_size"]
        )
        best_rmse, _, _ = fit_model(
            loaders=loaders,
            model_config=trial_config,
            graph_config=graph_config,
            normalizer=normalizer,
            device=device,
            learning_rate=params["learning_rate"],
            weight_decay=params["weight_decay"],
            max_epochs=OPTUNA_EPOCHS,
            patience=OPTUNA_PATIENCE,
            checkpoint_path=None,
            logger=None,
            trial=trial,
        )
        return float(best_rmse)

    sampler = optuna.samplers.TPESampler(seed=SEED)
    pruner = optuna.pruners.MedianPruner(n_startup_trials=5, n_warmup_steps=10)
    study = optuna.create_study(direction="minimize", sampler=sampler, pruner=pruner)
    logger.info(
        "Starting Optuna: n_trials=%d, timeout=%r.", OPTUNA_N_TRIALS, OPTUNA_TIMEOUT
    )
    study.optimize(objective, n_trials=OPTUNA_N_TRIALS, timeout=OPTUNA_TIMEOUT)
    result = {
        "best_value_val_rmse_eV": float(study.best_value),
        "best_params": study.best_params,
        "n_trials": len(study.trials),
    }
    output_path = Path(TABLE_DIR) / f"{MODEL_NAME}_optuna_best_params.json"
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(result, handle, ensure_ascii=False, indent=2)
    logger.info("Optuna best validation RMSE: %.6f eV", study.best_value)
    logger.info("Optuna best parameters: %s", study.best_params)
    return result


def prediction_frame(result: Dict[str, Any], split_name: str) -> pd.DataFrame:

    true = np.asarray(result["y_true"], dtype=np.float64)
    pred = np.asarray(result["y_pred"], dtype=np.float64)
    return pd.DataFrame(
        {
            "split": split_name,
            "cif_name": result["cif_name"],
            "true_bandgap_eV": true,
            "predicted_bandgap_eV": pred,
            "error_eV": pred - true,
            "absolute_error_eV": np.abs(pred - true),
        }
    )


def save_prediction_dat(results: Dict[str, Dict[str, Any]]) -> Dict[str, pd.DataFrame]:

    frames: Dict[str, pd.DataFrame] = {}
    for split_name in ["train", "val", "test"]:
        frame = prediction_frame(results[split_name], split_name)
        frames[split_name] = frame
        frame.drop(columns="split").to_csv(
            Path(DAT_DIR) / f"{MODEL_NAME}_parity_{split_name}.dat",
            sep="\t",
            index=False,
            float_format="%.8f",
        )
    frames["all"] = pd.concat(
        [frames["train"], frames["val"], frames["test"]], ignore_index=True
    )
    frames["all"].to_csv(
        Path(DAT_DIR) / f"{MODEL_NAME}_parity_all.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )
    return frames


def save_metrics_table(results: Dict[str, Dict[str, Any]]) -> Path:

    rows = []
    for split_name in ["train", "val", "test"]:
        metrics = results[split_name]["metrics"]
        rows.append(
            {
                "split": split_name,
                "n_samples": len(results[split_name]["y_true"]),
                "mae_eV": metrics["mae"],
                "rmse_eV": metrics["rmse"],
                "r2": metrics["r2"],
            }
        )
    output_path = Path(TABLE_DIR) / f"{MODEL_NAME}_metrics.dat"
    pd.DataFrame(rows).to_csv(output_path, sep="\t", index=False, float_format="%.6f")
    return output_path


def _global_axis_limits(frames: Dict[str, pd.DataFrame]) -> Tuple[float, float]:
    values = np.concatenate(
        [
            frames["all"]["true_bandgap_eV"].to_numpy(dtype=float),
            frames["all"]["predicted_bandgap_eV"].to_numpy(dtype=float),
        ]
    )
    minimum, maximum = float(np.min(values)), float(np.max(values))
    span = maximum - minimum
    padding = 0.05 * span if span > 1.0e-12 else 1.0
    return minimum - padding, maximum + padding


def _style_axis(axis: plt.Axes) -> None:
    axis.tick_params(axis="both", labelsize=TICK_FONTSIZE)
    if USE_GRID:
        axis.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)


def plot_parity(
    frame: pd.DataFrame,
    split_name: str,
    metrics: Dict[str, float],
    limits: Tuple[float, float],
) -> None:

    color_map = {"train": COLORS[0], "val": COLORS[1], "test": COLORS[2]}
    label_map = {"train": "Train", "val": "Validation", "test": "Test"}
    fig, axis = plt.subplots(figsize=FIG_SIZE, constrained_layout=True)
    axis.scatter(
        frame["true_bandgap_eV"],
        frame["predicted_bandgap_eV"],
        s=MARKER_SIZE,
        color=color_map[split_name],
        alpha=0.75,
        edgecolors="none",
        label=label_map[split_name],
    )
    axis.plot(
        limits, limits, color="0.25", linestyle="--", linewidth=1.2, label="y = x"
    )
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("True bandgap (eV)", fontsize=LABEL_FONTSIZE)
    axis.set_ylabel("Predicted bandgap (eV)", fontsize=LABEL_FONTSIZE)
    axis.set_title(
        f"{MODEL_NAME}: {label_map[split_name]} parity plot", fontsize=TITLE_FONTSIZE
    )
    annotation = (
        f"MAE = {metrics['mae']:.4f} eV\n"
        f"RMSE = {metrics['rmse']:.4f} eV\n"
        f"$R^2$ = {metrics['r2']:.4f}"
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
            "alpha": 0.82,
            "edgecolor": "0.7",
        },
    )
    axis.legend(fontsize=LEGEND_FONTSIZE, loc="lower right")
    _style_axis(axis)
    fig.savefig(
        Path(FIGURE_DIR) / f"{MODEL_NAME}_parity_{split_name}.jpg",
        dpi=FIG_DPI,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_combined_parity(
    frames: Dict[str, pd.DataFrame],
    results: Dict[str, Dict[str, Any]],
    limits: Tuple[float, float],
) -> None:

    labels = {"train": "Train", "val": "Validation", "test": "Test"}
    fig, axis = plt.subplots(figsize=FIG_SIZE, constrained_layout=True)
    lines: List[str] = []
    for index, split_name in enumerate(["train", "val", "test"]):
        frame = frames[split_name]
        axis.scatter(
            frame["true_bandgap_eV"],
            frame["predicted_bandgap_eV"],
            s=MARKER_SIZE,
            color=COLORS[index],
            alpha=0.72,
            edgecolors="none",
            label=labels[split_name],
        )
        metric = results[split_name]["metrics"]
        lines.append(
            f"{labels[split_name]}: MAE={metric['mae']:.4f}, "
            f"RMSE={metric['rmse']:.4f}, $R^2$={metric['r2']:.4f}"
        )
    axis.plot(
        limits, limits, color="0.25", linestyle="--", linewidth=1.2, label="y = x"
    )
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("True bandgap (eV)", fontsize=LABEL_FONTSIZE)
    axis.set_ylabel("Predicted bandgap (eV)", fontsize=LABEL_FONTSIZE)
    axis.set_title(f"{MODEL_NAME}: combined parity plot", fontsize=TITLE_FONTSIZE)
    axis.text(
        0.03,
        0.97,
        "\n".join(lines),
        transform=axis.transAxes,
        va="top",
        fontsize=ANNOTATION_FONTSIZE,
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "alpha": 0.82,
            "edgecolor": "0.7",
        },
    )
    axis.legend(fontsize=LEGEND_FONTSIZE, loc="lower right")
    _style_axis(axis)
    fig.savefig(
        Path(FIGURE_DIR) / f"{MODEL_NAME}_parity_all.jpg",
        dpi=FIG_DPI,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_rmse_curve(history: Sequence[Dict[str, float]], best_epoch: int) -> None:

    frame = pd.DataFrame(history)
    frame["epoch"] = frame["epoch"].astype(int)
    frame.to_csv(
        Path(DAT_DIR) / f"{MODEL_NAME}_rmse_curve.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )
    fig, axis = plt.subplots(figsize=FIG_SIZE, constrained_layout=True)
    axis.plot(
        frame["epoch"],
        frame["train_rmse_eV"],
        color=COLORS[0],
        linewidth=1.6,
        label="Train RMSE",
    )
    axis.plot(
        frame["epoch"],
        frame["val_rmse_eV"],
        color=COLORS[1],
        linewidth=1.6,
        label="Validation RMSE",
    )
    best_row = frame.loc[frame["epoch"] == best_epoch].iloc[0]
    axis.axvline(
        best_epoch,
        color=COLORS[3],
        linestyle="--",
        linewidth=1.2,
        label=f"Best epoch = {best_epoch}",
    )
    axis.scatter(
        [best_epoch], [best_row["val_rmse_eV"]], color=COLORS[4], s=42, zorder=5
    )
    axis.set_xlabel("Epoch", fontsize=LABEL_FONTSIZE)
    axis.set_ylabel("RMSE (eV)", fontsize=LABEL_FONTSIZE)
    axis.set_title(f"{MODEL_NAME}: training history", fontsize=TITLE_FONTSIZE)
    axis.legend(fontsize=LEGEND_FONTSIZE)
    _style_axis(axis)
    fig.savefig(
        Path(FIGURE_DIR) / f"{MODEL_NAME}_rmse_curve.jpg",
        dpi=FIG_DPI,
        bbox_inches="tight",
    )
    plt.close(fig)


@torch.inference_mode()
def predict_cifs(
    checkpoint_path: Union[str, Path],
    cif_paths: Sequence[Union[str, Path]],
    device: Optional[torch.device] = None,
    batch_size: int = BATCH_SIZE,
) -> pd.DataFrame:
    """Predict bandgaps for new CIF files from a saved checkpoint."""

    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if len(cif_paths) == 0:
        return pd.DataFrame(columns=["cif_name", "predicted_bandgap_eV"])
    model, graph_config, normalizer, _ = load_trained_model(checkpoint_path, device)
    graphs: List[CrystalGraphData] = []
    for raw_path in cif_paths:
        path = Path(raw_path)
        graph = load_or_build_graph(path, graph_config)
        graphs.append(
            CrystalGraphData(
                z=graph["z"],
                edge_index=graph["edge_index"],
                edge_attr=graph["edge_attr"],
                cif_name=str(path),
                num_nodes=int(graph["num_nodes"].item()),
            )
        )
    loader = DataLoader(graphs, batch_size=batch_size, shuffle=False)
    names: List[str] = []
    predictions: List[np.ndarray] = []
    model.eval()
    for batch in loader:
        names.extend(str(name) for name in batch.cif_name)
        batch = batch.to(device)
        with amp_context(device):
            output = model(batch).view(-1)
        predictions.append(normalizer.denormalize_tensor(output).float().cpu().numpy())
    return pd.DataFrame(
        {"cif_name": names, "predicted_bandgap_eV": np.concatenate(predictions)}
    )


def log_environment(logger: logging.Logger, device: torch.device) -> None:

    logger.info("Model name: %s", MODEL_NAME)
    logger.info("Run version: %s", RUN_VERSION)
    logger.info("Output directory: %s", Path(OUTPUT_ROOT).resolve())
    logger.info("Python: %s", sys.version.replace("\n", " "))
    logger.info("Platform: %s", platform.platform())
    logger.info("PyTorch: %s", torch.__version__)
    try:
        import torch_geometric

        logger.info("PyTorch Geometric: %s", torch_geometric.__version__)
    except Exception:
        pass
    logger.info("Compute device: %s", device)
    logger.info("CUDA available: %s", torch.cuda.is_available())
    if device.type == "cuda":
        properties = torch.cuda.get_device_properties(device)
        logger.info("GPU: %s", properties.name)
        logger.info("GPU total memory: %.3f GiB", properties.total_memory / 1024**3)
        logger.info("CUDA Toolkit: %s", torch.version.cuda)
    logger.info("Excel: %s", Path(EXCEL_PATH).resolve())
    logger.info("CIF directory: %s", Path(CIF_DIR).resolve())
    logger.info(
        "Excel columns: first=CIF filename; second=%s (%s)", TARGET_NAME, TARGET_UNIT
    )
    logger.info("Target: %s (%s)", TARGET_NAME, TARGET_UNIT)
    logger.info("Random seed: %d", SEED)
    logger.info("Loss function: nn.MSELoss(); best model selected by validation RMSE")


def main() -> None:

    script_dir = Path(__file__).resolve().parent
    os.chdir(script_dir)
    if not RUN_VERSION.strip() or "/" in RUN_VERSION or "\\" in RUN_VERSION:
        raise ValueError("RUN_VERSION must be a non-empty directory-safe name.")
    create_directories()
    logger = setup_logger()
    total_start = time.perf_counter()
    set_global_seed(SEED)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    plt.rcParams["font.family"] = FONT_FAMILY
    plt.rcParams["axes.unicode_minus"] = False
    log_environment(logger, device)

    graph_config = GraphConfig()
    logger.info("Graph configuration: %s", asdict(graph_config))
    data_start = time.perf_counter()
    raw_records = load_excel_records(logger)
    valid_records = validate_and_cache_records(raw_records, graph_config, logger)
    splits = create_or_load_split(valid_records, logger)
    logger.info(
        "Data split: train=%d, validation=%d, test=%d (ratios %.2f/%.2f/%.2f).",
        len(splits["train"]),
        len(splits["val"]),
        len(splits["test"]),
        TRAIN_RATIO,
        VAL_RATIO,
        TEST_RATIO,
    )
    normalizer = TargetNormalizer.fit(
        [item["target"] for item in splits["train"]], enabled=STANDARDIZE_TARGET
    )
    logger.info("Target normalization: %s", normalizer.to_dict())
    data_seconds = time.perf_counter() - data_start
    logger.info(
        "Data loading and graph construction time: %.3f s (%s)",
        data_seconds,
        format_duration(data_seconds),
    )

    selected = {
        "batch_size": BATCH_SIZE,
        "learning_rate": LEARNING_RATE,
        "weight_decay": WEIGHT_DECAY,
        "hidden_dim": HIDDEN_DIM,
        "depth": DEPTH,
    }
    logger.info("Optuna enabled: %s", USE_OPTUNA)
    if USE_OPTUNA:
        optuna_result = run_optuna(splits, graph_config, normalizer, device, logger)
        selected.update(optuna_result["best_params"])

    model_config = ModelConfig(
        hidden_dim=int(selected["hidden_dim"]),
        depth=int(selected["depth"]),
    )
    loaders = make_loaders(
        splits, graph_config, normalizer, batch_size=int(selected["batch_size"])
    )
    model_for_log = DenseGNN(model_config)
    total_parameters = sum(
        parameter.numel() for parameter in model_for_log.parameters()
    )
    trainable_parameters = sum(
        parameter.numel()
        for parameter in model_for_log.parameters()
        if parameter.requires_grad
    )
    logger.info("Model configuration: %s", asdict(model_config))
    logger.info("Model architecture:\n%s", model_for_log)
    logger.info(
        "Total parameters: %d; trainable parameters: %d",
        total_parameters,
        trainable_parameters,
    )
    logger.info(
        "Optimizer: AdamW; initial learning rate=%.6g; weight_decay=%.6g; "
        "batch_size=%d; scheduler=ReduceLROnPlateau.",
        float(selected["learning_rate"]),
        float(selected["weight_decay"]),
        int(selected["batch_size"]),
    )
    del model_for_log

    train_start = time.perf_counter()
    best_rmse, best_epoch, history = fit_model(
        loaders=loaders,
        model_config=model_config,
        graph_config=graph_config,
        normalizer=normalizer,
        device=device,
        learning_rate=float(selected["learning_rate"]),
        weight_decay=float(selected["weight_decay"]),
        max_epochs=MAX_EPOCHS,
        patience=PATIENCE,
        checkpoint_path=CHECKPOINT_PATH,
        logger=logger,
    )
    train_seconds = time.perf_counter() - train_start
    logger.info("Best epoch: %d; best validation RMSE: %.6f eV", best_epoch, best_rmse)
    logger.info(
        "Training time: %.3f s (%s)", train_seconds, format_duration(train_seconds)
    )

    evaluation_start = time.perf_counter()
    best_model, saved_graph_config, saved_normalizer, checkpoint = load_trained_model(
        CHECKPOINT_PATH, device
    )
    if saved_graph_config != graph_config:
        raise RuntimeError(
            "Graph parameters in the checkpoint do not match the current run."
        )
    criterion = nn.MSELoss()
    results = {
        split_name: evaluate_loader(
            best_model, loaders[split_name], criterion, saved_normalizer, device
        )
        for split_name in ["train", "val", "test"]
    }
    if int(checkpoint["best_epoch"]) != best_epoch:
        raise RuntimeError(
            "The reloaded best epoch does not match the training record."
        )

    frames = save_prediction_dat(results)
    metrics_path = save_metrics_table(results)
    limits = _global_axis_limits(frames)
    for split_name in ["train", "val", "test"]:
        plot_parity(
            frames[split_name], split_name, results[split_name]["metrics"], limits
        )
    plot_combined_parity(frames, results, limits)
    plot_rmse_curve(history, best_epoch)

    for split_name in ["train", "val", "test"]:
        metric = results[split_name]["metrics"]
        logger.info(
            "%s | n=%d | MAE=%.6f eV | RMSE=%.6f eV | R2=%.6f",
            split_name,
            len(results[split_name]["y_true"]),
            metric["mae"],
            metric["rmse"],
            metric["r2"],
        )
    evaluation_seconds = time.perf_counter() - evaluation_start
    total_seconds = time.perf_counter() - total_start
    logger.info("Best model: %s", Path(CHECKPOINT_PATH).resolve())
    logger.info("Metrics table: %s", metrics_path.resolve())
    logger.info(
        "Evaluation, saving, and plotting time: %.3f s (%s)",
        evaluation_seconds,
        format_duration(evaluation_seconds),
    )
    logger.info(
        "Total runtime: %.3f s (%s)", total_seconds, format_duration(total_seconds)
    )


if __name__ == "__main__":
    main()

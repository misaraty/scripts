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
from pymatgen.core import Structure
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from torch import Tensor, nn
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch_geometric.nn import global_mean_pool
from torch_geometric.utils import scatter
from tqdm.auto import tqdm


# User configuration

MODEL_NAME = "FlyCrysNet"
RUN_VERSION = "v3"

EXCEL_PATH = "./data.xlsx"
CIF_DIR = "./cif"

TARGET_NAME = "bandgap"
TARGET_UNIT = "eV"

SEED = 42
TRAIN_RATIO = 0.80
VAL_RATIO = 0.10
TEST_RATIO = 0.10

# Crystal graph settings
CUTOFF = 8.0
MAX_CUTOFF = 20.0
CUTOFF_STEP = 2.0
MAX_NEIGHBORS = 12
MAX_ATOMIC_NUMBER = 118
CACHE_GRAPHS = True

CRYSTAL_LAYERS = 4
ATOM_INPUT_FEATURES = MAX_ATOMIC_NUMBER
EDGE_INPUT_FEATURES = 80
EMBEDDING_FEATURES = 64
HIDDEN_FEATURES = 192
DROPOUT = 0.10

# MaleCNS-derived sparse topology
CONNECTOME_FEATHER = (
    "./fly_connectome/connectome-weights-male-cns-v1.0-minconf-0.5.feather"
)
CONNECTOME_NODES = 256
CONNECTOME_CHANNELS = 8
CONNECTOME_STEPS = 1
TOPOLOGY_GATE_INIT = 0.10
EDGE_WEIGHT_POWER_INIT = 1.00

# Necessary topology ablations
RUN_ABLATIONS = True
ABLATION_RANDOM_SEED = 2026
ABLATION_SWAP_FACTOR = 10

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
FONT_FAMILY = ["Arial", "Microsoft YaHei", "DejaVu Sans"]
COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:red"]

# Output paths
OUTPUT_ROOT = f"./{MODEL_NAME}_{RUN_VERSION}"
FIGURE_DIR = f"{OUTPUT_ROOT}/figure"
DAT_DIR = f"{OUTPUT_ROOT}/dat"
TABLE_DIR = f"{OUTPUT_ROOT}/table"
LOG_DIR = f"{OUTPUT_ROOT}/log"
SPLIT_DIR = f"{OUTPUT_ROOT}/split"
CACHE_DIR = f"{OUTPUT_ROOT}/cache"
ABLATION_CHECKPOINT_DIR = f"{OUTPUT_ROOT}/ablation_checkpoints"
GRAPH_CACHE_DIR = f"{CACHE_DIR}/crystal_graphs"
CONNECTOME_CACHE_PATH = f"{CACHE_DIR}/MaleCNS_top{CONNECTOME_NODES}.npz"
CHECKPOINT_PATH = f"{OUTPUT_ROOT}/{MODEL_NAME}_best.pt"
SPLIT_FILE = f"{SPLIT_DIR}/{MODEL_NAME}_split.csv"


@dataclass(frozen=True)
class GraphConfig:

    cutoff: float = CUTOFF
    max_cutoff: float = MAX_CUTOFF
    cutoff_step: float = CUTOFF_STEP
    max_neighbors: int = MAX_NEIGHBORS
    max_atomic_number: int = MAX_ATOMIC_NUMBER


@dataclass(frozen=True)
class ModelConfig:

    atom_input_features: int = ATOM_INPUT_FEATURES
    edge_cutoff: float = CUTOFF
    edge_input_features: int = EDGE_INPUT_FEATURES
    embedding_features: int = EMBEDDING_FEATURES
    hidden_features: int = HIDDEN_FEATURES
    crystal_layers: int = CRYSTAL_LAYERS
    connectome_nodes: int = CONNECTOME_NODES
    connectome_channels: int = CONNECTOME_CHANNELS
    connectome_steps: int = CONNECTOME_STEPS
    topology_gate_init: float = TOPOLOGY_GATE_INIT
    edge_weight_power_init: float = EDGE_WEIGHT_POWER_INIT
    dropout: float = DROPOUT


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

    for path in [
        FIGURE_DIR,
        DAT_DIR,
        TABLE_DIR,
        LOG_DIR,
        SPLIT_DIR,
        ABLATION_CHECKPOINT_DIR,
    ]:
        Path(path).mkdir(parents=True, exist_ok=True)
    if CACHE_GRAPHS:
        Path(GRAPH_CACHE_DIR).mkdir(parents=True, exist_ok=True)
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


def _adaptive_neighbor_list(
    structure: Structure, config: GraphConfig
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:

    radius = float(config.cutoff)
    n_atoms = len(structure)
    while radius <= config.max_cutoff + 1e-12:
        center, neighbor, images, distances = structure.get_neighbor_list(
            r=radius, numerical_tol=1.0e-8, exclude_self=True
        )
        center = np.asarray(center, dtype=np.int64)
        neighbor = np.asarray(neighbor, dtype=np.int64)
        images = np.asarray(images, dtype=np.float64)
        distances = np.asarray(distances, dtype=np.float64)
        if center.size:
            counts = np.bincount(center, minlength=n_atoms)
            if int(counts.min()) >= config.max_neighbors:
                break
        radius += config.cutoff_step
    else:
        raise ValueError(
            f"At least one atom has fewer than {config.max_neighbors} neighbors "
            f"within the maximum cutoff of {config.max_cutoff:.3f} Å."
        )

    # Keep the full neighbor shell at the k-th-neighbor distance.
    selected: List[np.ndarray] = []
    for atom_index in range(n_atoms):
        candidates = np.flatnonzero(center == atom_index)
        order = candidates[np.argsort(distances[candidates], kind="stable")]
        kth_distance = distances[order[config.max_neighbors - 1]]
        keep = order[distances[order] <= kth_distance + 1.0e-8]
        selected.append(keep)
    indices = np.concatenate(selected)

    # Add the reverse direction of every periodic edge.
    directed: Dict[Tuple[int, int, Tuple[int, int, int]], float] = {}
    for index in indices:
        source = int(center[index])
        target = int(neighbor[index])
        image = tuple(int(value) for value in np.rint(images[index]))
        distance = float(distances[index])
        directed[(source, target, image)] = distance
        reverse_image = tuple(-value for value in image)
        directed[(target, source, reverse_image)] = distance
    ordered = sorted(
        directed.items(), key=lambda item: (item[0][0], item[1], item[0][1], item[0][2])
    )
    selected_center = np.asarray([item[0][0] for item in ordered], dtype=np.int64)
    selected_neighbor = np.asarray([item[0][1] for item in ordered], dtype=np.int64)
    selected_images = np.asarray([item[0][2] for item in ordered], dtype=np.float64)
    selected_distances = np.asarray([item[1] for item in ordered], dtype=np.float64)
    return selected_center, selected_neighbor, selected_images, selected_distances


def build_graph_tensors(
    cif_path: Union[str, Path], config: GraphConfig
) -> Dict[str, Tensor]:

    structure = Structure.from_file(str(cif_path), primitive=False, sort=False)
    if len(structure) == 0:
        raise ValueError("The CIF structure contains no atoms.")
    atomic_numbers = np.asarray(
        [int(site.specie.Z) for site in structure], dtype=np.int64
    )
    if np.any(atomic_numbers < 1) or np.any(atomic_numbers > config.max_atomic_number):
        raise ValueError(
            f"Atomic numbers outside the supported range 1-{config.max_atomic_number}: "
            f"{np.unique(atomic_numbers).tolist()}"
        )

    src, dst, images, _ = _adaptive_neighbor_list(structure, config)
    frac = np.asarray(structure.frac_coords, dtype=np.float64)
    lattice = np.asarray(structure.lattice.matrix, dtype=np.float64)
    edge_vec = (frac[dst] + images - frac[src]) @ lattice
    lengths = np.linalg.norm(edge_vec, axis=1)
    valid = lengths > 1.0e-8
    src, dst, edge_vec, lengths = (
        src[valid],
        dst[valid],
        edge_vec[valid],
        lengths[valid],
    )
    if src.size == 0:
        raise ValueError("No valid crystal edges were generated.")

    return {
        "atomic_numbers": torch.as_tensor(atomic_numbers, dtype=torch.long),
        "edge_index": torch.as_tensor(np.vstack([src, dst]), dtype=torch.long),
        "edge_vec": torch.as_tensor(edge_vec, dtype=torch.float32),
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
        "cache_version": 1,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, ensure_ascii=False).encode("utf-8")
    ).hexdigest()[:24]
    return Path(GRAPH_CACHE_DIR) / f"{digest}.pt"


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
            required = {
                "atomic_numbers",
                "edge_index",
                "edge_vec",
                "num_nodes",
            }
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
        records, desc="Validating CIFs / building crystal graphs", unit="cif"
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
            atomic_numbers=graph["atomic_numbers"],
            edge_index=graph["edge_index"],
            edge_vec=graph["edge_vec"],
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


class RBFExpansion(nn.Module):

    def __init__(self, vmin: float, vmax: float, bins: int) -> None:
        super().__init__()
        centers = torch.linspace(vmin, vmax, bins)
        self.register_buffer("centers", centers)
        spacing = float((centers[1] - centers[0]).item()) if bins > 1 else 1.0
        self.gamma = 1.0 / spacing

    def forward(self, values: Tensor) -> Tensor:
        return torch.exp(-self.gamma * (values.unsqueeze(-1) - self.centers) ** 2)


class SafeBatchNorm1d(nn.BatchNorm1d):

    def forward(self, values: Tensor) -> Tensor:
        if self.training and values.size(0) < 2:
            return F.batch_norm(
                values,
                self.running_mean,
                self.running_var,
                self.weight,
                self.bias,
                training=False,
                momentum=0.0,
                eps=self.eps,
            )
        return super().forward(values)


class MLPLayer(nn.Module):
    """Linear, batch normalization, and SiLU."""

    def __init__(self, in_features: int, out_features: int) -> None:
        super().__init__()
        self.layer = nn.Sequential(
            nn.Linear(in_features, out_features),
            SafeBatchNorm1d(out_features),
            nn.SiLU(),
        )

    def forward(self, values: Tensor) -> Tensor:
        return self.layer(values)


class EdgeGatedGraphConv(nn.Module):

    def __init__(
        self, input_features: int, output_features: int, residual: bool = True
    ) -> None:
        super().__init__()
        self.residual = residual and input_features == output_features
        self.src_gate = nn.Linear(input_features, output_features)
        self.dst_gate = nn.Linear(input_features, output_features)
        self.edge_gate = nn.Linear(input_features, output_features)
        self.bn_edges = SafeBatchNorm1d(output_features)
        self.src_update = nn.Linear(input_features, output_features)
        self.dst_update = nn.Linear(input_features, output_features)
        self.bn_nodes = SafeBatchNorm1d(output_features)

    def forward(
        self, edge_index: Tensor, node_features: Tensor, edge_features: Tensor
    ) -> Tuple[Tensor, Tensor]:
        src, dst = edge_index[0], edge_index[1]
        gate_logits = (
            self.src_gate(node_features[src])
            + self.dst_gate(node_features[dst])
            + self.edge_gate(edge_features)
        )
        gates = torch.sigmoid(gate_logits)
        messages = gates * self.dst_update(node_features[src])
        numerator = scatter(
            messages, dst, dim=0, dim_size=node_features.size(0), reduce="sum"
        )
        denominator = scatter(
            gates, dst, dim=0, dim_size=node_features.size(0), reduce="sum"
        )
        updated_nodes = self.src_update(node_features) + numerator / (
            denominator + 1.0e-6
        )
        updated_nodes = F.silu(self.bn_nodes(updated_nodes))
        updated_edges = F.silu(self.bn_edges(gate_logits))
        if self.residual:
            updated_nodes = node_features + updated_nodes
            updated_edges = edge_features + updated_edges
        return updated_nodes, updated_edges


def _arrow_reader(feather_path: Path):
    try:
        import pyarrow as pa
        import pyarrow.ipc as ipc
    except ImportError as exc:
        raise ImportError(
            "pyarrow is required to read the MaleCNS Feather file."
        ) from exc
    memory_map = pa.memory_map(str(feather_path), "r")
    return memory_map, ipc.open_file(memory_map)


def _membership_positions(
    values: np.ndarray, sorted_ids: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    positions = np.searchsorted(sorted_ids, values)
    valid = positions < len(sorted_ids)
    valid_positions = positions[valid]
    matches = np.zeros(len(values), dtype=bool)
    matches[np.flatnonzero(valid)] = sorted_ids[valid_positions] == values[valid]
    return matches, positions


def build_connectome_topology(
    feather_path: Union[str, Path],
    n_nodes: int,
    cache_path: Union[str, Path],
    logger: logging.Logger,
) -> Dict[str, Tensor]:
    feather_path = Path(feather_path)
    cache_path = Path(cache_path)
    if not feather_path.is_file():
        raise FileNotFoundError(
            f"MaleCNS Feather file not found: {feather_path.resolve()}"
        )
    if n_nodes < 2:
        raise ValueError("CONNECTOME_NODES must be at least 2.")
    feather_stat = feather_path.stat()
    if cache_path.is_file():
        cached = np.load(cache_path)
        cache_matches = (
            int(cached["n_nodes"].item()) == n_nodes
            and "source_size" in cached
            and "source_mtime_ns" in cached
            and int(cached["source_size"].item()) == feather_stat.st_size
            and int(cached["source_mtime_ns"].item()) == feather_stat.st_mtime_ns
        )
        if cache_matches:
            logger.info("Loaded cached MaleCNS topology: %s", cache_path)
            return {
                "edge_index": torch.as_tensor(cached["edge_index"], dtype=torch.long),
                "edge_weight": torch.as_tensor(
                    cached["edge_weight"], dtype=torch.float32
                ),
                "body_ids": torch.as_tensor(cached["body_ids"], dtype=torch.long),
            }

    logger.info("MaleCNS pass 1/2: computing weighted neuron degrees.")
    memory_map, reader = _arrow_reader(feather_path)
    names = reader.schema.names
    required = {"body_pre", "body_post", "weight"}
    if not required.issubset(names):
        memory_map.close()
        raise ValueError(f"Unexpected MaleCNS schema: {names}")
    pre_col = names.index("body_pre")
    post_col = names.index("body_post")
    weight_col = names.index("weight")
    weighted_degree: Dict[int, float] = {}
    for batch_number in range(reader.num_record_batches):
        batch = reader.get_batch(batch_number)
        pre = np.asarray(batch.column(pre_col), dtype=np.int64)
        post = np.asarray(batch.column(post_col), dtype=np.int64)
        weight = np.asarray(batch.column(weight_col), dtype=np.float64)
        valid = np.isfinite(weight) & (weight > 0)
        pre, post, weight = pre[valid], post[valid], weight[valid]
        for node_ids in (pre, post):
            unique_ids, inverse = np.unique(node_ids, return_inverse=True)
            sums = np.bincount(inverse, weights=weight)
            for node_id, value in zip(unique_ids.tolist(), sums.tolist()):
                weighted_degree[int(node_id)] = weighted_degree.get(
                    int(node_id), 0.0
                ) + float(value)
        if (
            batch_number + 1
        ) % 25 == 0 or batch_number + 1 == reader.num_record_batches:
            logger.info(
                "MaleCNS degree batches: %d/%d",
                batch_number + 1,
                reader.num_record_batches,
            )
    memory_map.close()
    if len(weighted_degree) < n_nodes:
        raise ValueError(
            f"MaleCNS contains only {len(weighted_degree)} connected neurons, "
            f"fewer than CONNECTOME_NODES={n_nodes}."
        )

    ranked = sorted(weighted_degree.items(), key=lambda item: (-item[1], item[0]))
    body_ids = np.asarray([item[0] for item in ranked[:n_nodes]], dtype=np.int64)
    sorted_ids = np.sort(body_ids)
    adjacency_sorted = np.zeros((n_nodes, n_nodes), dtype=np.float64)

    logger.info(
        "MaleCNS pass 2/2: extracting the %d-neuron directed induced subgraph.",
        n_nodes,
    )
    memory_map, reader = _arrow_reader(feather_path)
    for batch_number in range(reader.num_record_batches):
        batch = reader.get_batch(batch_number)
        pre = np.asarray(batch.column(pre_col), dtype=np.int64)
        post = np.asarray(batch.column(post_col), dtype=np.int64)
        weight = np.asarray(batch.column(weight_col), dtype=np.float64)
        pre_mask, pre_pos = _membership_positions(pre, sorted_ids)
        post_mask, post_pos = _membership_positions(post, sorted_ids)
        valid = (
            pre_mask & post_mask & (pre != post) & np.isfinite(weight) & (weight > 0)
        )
        if np.any(valid):
            np.add.at(
                adjacency_sorted, (post_pos[valid], pre_pos[valid]), weight[valid]
            )
        if (
            batch_number + 1
        ) % 25 == 0 or batch_number + 1 == reader.num_record_batches:
            logger.info(
                "MaleCNS subgraph batches: %d/%d",
                batch_number + 1,
                reader.num_record_batches,
            )
    memory_map.close()

    lookup = {int(node_id): index for index, node_id in enumerate(sorted_ids.tolist())}
    order = np.asarray([lookup[int(node_id)] for node_id in body_ids], dtype=np.int64)
    adjacency = adjacency_sorted[np.ix_(order, order)]
    adjacency = np.log1p(np.maximum(adjacency, 0.0))
    adjacency += np.eye(n_nodes, dtype=np.float64)
    row_sum = adjacency.sum(axis=1, keepdims=True)
    adjacency = adjacency / np.maximum(row_sum, 1.0e-12)
    dst, src = np.nonzero(adjacency)
    edge_index = np.vstack([src, dst]).astype(np.int64)
    edge_weight = adjacency[dst, src].astype(np.float32)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        cache_path,
        n_nodes=np.asarray(n_nodes, dtype=np.int64),
        source_size=np.asarray(feather_stat.st_size, dtype=np.int64),
        source_mtime_ns=np.asarray(feather_stat.st_mtime_ns, dtype=np.int64),
        edge_index=edge_index,
        edge_weight=edge_weight,
        body_ids=body_ids,
    )
    logger.info(
        "Cached MaleCNS-derived topology: nodes=%d, directed edges including self-loops=%d.",
        n_nodes,
        edge_index.shape[1],
    )
    return {
        "edge_index": torch.as_tensor(edge_index, dtype=torch.long),
        "edge_weight": torch.as_tensor(edge_weight, dtype=torch.float32),
        "body_ids": torch.as_tensor(body_ids, dtype=torch.long),
    }


def _normalize_incoming_weights(
    edge_index: Tensor, edge_weight: Tensor, n_nodes: int
) -> Tensor:
    weights = edge_weight.detach().cpu().float().clamp_min(1.0e-12)
    destination = edge_index[1].detach().cpu().long()
    totals = torch.zeros(n_nodes, dtype=torch.float32)
    totals.scatter_add_(0, destination, weights)
    return weights / totals[destination].clamp_min(1.0e-12)


def _random_topology(topology: Dict[str, Tensor], seed: int) -> Dict[str, Tensor]:
    n_nodes = int(topology["body_ids"].numel())
    edge_index = topology["edge_index"].detach().cpu().long()
    edge_weight = topology["edge_weight"].detach().cpu().float()
    off_diagonal = edge_index[0] != edge_index[1]
    n_off_diagonal = int(off_diagonal.sum().item())
    maximum = n_nodes * (n_nodes - 1)
    if n_off_diagonal > maximum:
        raise ValueError("The requested random topology contains too many edges.")

    rng = np.random.default_rng(seed)
    sampled = rng.choice(maximum, size=n_off_diagonal, replace=False)
    source = sampled // (n_nodes - 1)
    remainder = sampled % (n_nodes - 1)
    destination = remainder + (remainder >= source)
    self_nodes = np.arange(n_nodes, dtype=np.int64)
    random_edges = torch.as_tensor(
        np.vstack(
            [
                np.concatenate([source, self_nodes]),
                np.concatenate([destination, self_nodes]),
            ]
        ),
        dtype=torch.long,
    )

    off_weights = edge_weight[off_diagonal].numpy().copy()
    rng.shuffle(off_weights)
    original_self = edge_index[0] == edge_index[1]
    self_weights = edge_weight[original_self].numpy().copy()
    if self_weights.size != n_nodes:
        self_weights = np.ones(n_nodes, dtype=np.float32)
    random_weights = torch.as_tensor(
        np.concatenate([off_weights, self_weights]), dtype=torch.float32
    )
    random_weights = _normalize_incoming_weights(random_edges, random_weights, n_nodes)
    return {
        "edge_index": random_edges,
        "edge_weight": random_weights,
        "body_ids": topology["body_ids"].detach().cpu().clone(),
    }


def _degree_preserving_topology(
    topology: Dict[str, Tensor], seed: int, swap_factor: int
) -> Dict[str, Tensor]:
    n_nodes = int(topology["body_ids"].numel())
    edge_index = topology["edge_index"].detach().cpu().long()
    edge_weight = topology["edge_weight"].detach().cpu().float()
    off_mask = edge_index[0] != edge_index[1]
    off_edges = [
        (int(source), int(destination))
        for source, destination in edge_index[:, off_mask].t().tolist()
    ]
    edge_set = set(off_edges)
    rng = random.Random(seed)
    target_swaps = max(1, swap_factor * len(off_edges))
    accepted = 0
    attempts = 0
    maximum_attempts = max(1000, 30 * target_swaps)
    while accepted < target_swaps and attempts < maximum_attempts:
        attempts += 1
        first, second = rng.sample(range(len(off_edges)), 2)
        source_a, destination_a = off_edges[first]
        source_b, destination_b = off_edges[second]
        new_a = (source_a, destination_b)
        new_b = (source_b, destination_a)
        if (
            source_a == source_b
            or destination_a == destination_b
            or source_a == destination_b
            or source_b == destination_a
            or new_a in edge_set
            or new_b in edge_set
        ):
            continue
        edge_set.remove(off_edges[first])
        edge_set.remove(off_edges[second])
        edge_set.add(new_a)
        edge_set.add(new_b)
        off_edges[first] = new_a
        off_edges[second] = new_b
        accepted += 1
    if accepted < max(1, len(off_edges)):
        raise RuntimeError(
            "Degree-preserving randomization produced too few successful swaps."
        )

    self_nodes = torch.arange(n_nodes, dtype=torch.long)
    rewired_edges = torch.cat(
        [
            torch.as_tensor(off_edges, dtype=torch.long).t().contiguous(),
            torch.stack([self_nodes, self_nodes]),
        ],
        dim=1,
    )
    rng_np = np.random.default_rng(seed)
    off_weights = edge_weight[off_mask].numpy().copy()
    rng_np.shuffle(off_weights)
    self_weights = edge_weight[~off_mask].numpy().copy()
    if self_weights.size != n_nodes:
        self_weights = np.ones(n_nodes, dtype=np.float32)
    rewired_weights = torch.as_tensor(
        np.concatenate([off_weights, self_weights]), dtype=torch.float32
    )
    rewired_weights = _normalize_incoming_weights(
        rewired_edges, rewired_weights, n_nodes
    )
    return {
        "edge_index": rewired_edges,
        "edge_weight": rewired_weights,
        "body_ids": topology["body_ids"].detach().cpu().clone(),
    }


def build_ablation_topologies(
    topology: Dict[str, Tensor], logger: logging.Logger
) -> Dict[str, Dict[str, Any]]:
    n_nodes = int(topology["body_ids"].numel())
    self_nodes = torch.arange(n_nodes, dtype=torch.long)
    self_edges = torch.stack([self_nodes, self_nodes])
    self_only = {
        "edge_index": self_edges,
        "edge_weight": torch.ones(n_nodes, dtype=torch.float32),
        "body_ids": topology["body_ids"].detach().cpu().clone(),
    }
    unweighted = {
        "edge_index": topology["edge_index"].detach().cpu().clone(),
        "edge_weight": _normalize_incoming_weights(
            topology["edge_index"],
            torch.ones_like(topology["edge_weight"], dtype=torch.float32),
            n_nodes,
        ),
        "body_ids": topology["body_ids"].detach().cpu().clone(),
    }
    variants: Dict[str, Dict[str, Any]] = {
        "full_malecns": {
            "label": "Full MaleCNS",
            "description": "Authentic MaleCNS-derived weighted topology",
            "topology": topology,
        },
        "self_only": {
            "label": "Self-only",
            "description": "No inter-node connectome edges",
            "topology": self_only,
        },
        "random_topology": {
            "label": "Random topology",
            "description": "Random directed topology with matched node and edge counts",
            "topology": _random_topology(topology, ABLATION_RANDOM_SEED),
        },
        "degree_preserving": {
            "label": "Degree-preserving",
            "description": "Degree-preserving randomized MaleCNS topology",
            "topology": _degree_preserving_topology(
                topology, ABLATION_RANDOM_SEED, ABLATION_SWAP_FACTOR
            ),
        },
        "unweighted_malecns": {
            "label": "Unweighted MaleCNS",
            "description": "Authentic MaleCNS edges with uniform incoming weights",
            "topology": unweighted,
        },
    }
    for variant in variants.values():
        variant_topology = variant["topology"]
        logger.info(
            "Ablation topology prepared | %s | nodes=%d | edges=%d | %s",
            variant["label"],
            int(variant_topology["body_ids"].numel()),
            int(variant_topology["edge_weight"].numel()),
            variant["description"],
        )
    return variants


class CrystalEncoder(nn.Module):

    def __init__(self, config: ModelConfig) -> None:
        super().__init__()
        self.config = config
        self.atom_embedding = nn.Embedding(
            config.atom_input_features + 1, config.hidden_features, padding_idx=0
        )
        self.edge_embedding = nn.Sequential(
            RBFExpansion(0.0, config.edge_cutoff, config.edge_input_features),
            MLPLayer(config.edge_input_features, config.embedding_features),
            MLPLayer(config.embedding_features, config.hidden_features),
        )
        self.layers = nn.ModuleList(
            [
                EdgeGatedGraphConv(config.hidden_features, config.hidden_features)
                for _ in range(config.crystal_layers)
            ]
        )
        self.norm = nn.LayerNorm(config.hidden_features)

    def forward(self, data: CrystalGraphData) -> Tensor:
        atom_features = self.atom_embedding(data.atomic_numbers.long())
        bond_lengths = torch.linalg.vector_norm(data.edge_vec, dim=1)
        bond_features = self.edge_embedding(bond_lengths)
        for layer in self.layers:
            atom_features, bond_features = layer(
                data.edge_index, atom_features, bond_features
            )
        return self.norm(global_mean_pool(atom_features, data.batch))


class MaleCNSBlock(nn.Module):

    def __init__(
        self,
        channels: int,
        dropout: float,
        topology_gate_init: float,
        edge_weight_power_init: float,
    ) -> None:
        super().__init__()
        if not 0.0 < topology_gate_init < 1.0:
            raise ValueError("topology_gate_init must be strictly between 0 and 1.")
        if not 0.0 < edge_weight_power_init < 2.0:
            raise ValueError("edge_weight_power_init must be strictly between 0 and 2.")
        self.self_layer = nn.Linear(channels, channels)
        self.message_layer = nn.Linear(channels, channels, bias=False)
        self.pre_norm = nn.LayerNorm(channels)
        self.output_norm = nn.LayerNorm(channels)
        self.dropout = nn.Dropout(dropout)
        gate_logit = math.log(topology_gate_init / (1.0 - topology_gate_init))
        power_ratio = edge_weight_power_init / 2.0
        power_logit = math.log(power_ratio / (1.0 - power_ratio))
        self.topology_gate_logit = nn.Parameter(
            torch.tensor(gate_logit, dtype=torch.float32)
        )
        self.edge_weight_power_logit = nn.Parameter(
            torch.tensor(power_logit, dtype=torch.float32)
        )

    def topology_gate(self) -> Tensor:
        return torch.sigmoid(self.topology_gate_logit)

    def edge_weight_power(self) -> Tensor:
        return 2.0 * torch.sigmoid(self.edge_weight_power_logit)

    def effective_edge_weight(
        self, edge_index: Tensor, edge_weight: Tensor, n_nodes: int
    ) -> Tensor:
        destination = edge_index[1]
        powered = edge_weight.clamp_min(1.0e-12).pow(self.edge_weight_power())
        incoming_sum = scatter(
            powered, destination, dim=0, dim_size=n_nodes, reduce="sum"
        )
        return powered / incoming_sum[destination].clamp_min(1.0e-12)

    def forward(
        self, hidden: Tensor, edge_index: Tensor, edge_weight: Tensor
    ) -> Tensor:
        src, dst = edge_index
        normalized = self.pre_norm(hidden)
        local_update = F.gelu(self.self_layer(normalized))
        local_hidden = hidden + self.dropout(local_update)
        effective_weight = self.effective_edge_weight(
            edge_index, edge_weight, hidden.size(1)
        )
        messages = normalized[:, src, :] * effective_weight.view(1, -1, 1)
        aggregated = scatter(
            messages, dst, dim=1, dim_size=hidden.size(1), reduce="sum"
        )
        topology_delta = aggregated - normalized
        topology_update = F.gelu(self.message_layer(topology_delta))
        gated_topology_update = self.topology_gate() * topology_update
        return self.output_norm(local_hidden + self.dropout(gated_topology_update))


class FlyConnectome(nn.Module):
    """Crystal encoder followed by MaleCNS-derived sparse message passing."""

    def __init__(
        self,
        config: ModelConfig,
        connectome_edge_index: Tensor,
        connectome_edge_weight: Tensor,
    ) -> None:
        super().__init__()
        self.config = config
        if connectome_edge_index.ndim != 2 or connectome_edge_index.size(0) != 2:
            raise ValueError("Connectome edge_index must have shape [2, n_edges].")
        if connectome_edge_weight.numel() != connectome_edge_index.size(1):
            raise ValueError("Connectome edge weights do not match edge_index.")
        self.register_buffer(
            "connectome_edge_index", connectome_edge_index.long().contiguous()
        )
        self.register_buffer(
            "connectome_edge_weight", connectome_edge_weight.float().contiguous()
        )
        self.crystal_encoder = CrystalEncoder(config)
        projected = config.connectome_nodes * config.connectome_channels
        self.input_projection = nn.Sequential(
            nn.LayerNorm(config.hidden_features),
            nn.Linear(config.hidden_features, projected),
        )
        self.connectome_blocks = nn.ModuleList(
            [
                MaleCNSBlock(
                    config.connectome_channels,
                    config.dropout,
                    config.topology_gate_init,
                    config.edge_weight_power_init,
                )
                for _ in range(config.connectome_steps)
            ]
        )
        self.attention = nn.Linear(config.connectome_channels, 1)
        self.global_skip = nn.Linear(config.hidden_features, config.connectome_channels)
        output_hidden = max(32, 4 * config.connectome_channels)
        self.output = nn.Sequential(
            nn.Linear(2 * config.connectome_channels, output_hidden),
            nn.SiLU(),
            nn.Dropout(config.dropout),
            nn.Linear(output_hidden, 1),
        )

    def forward(self, data: CrystalGraphData) -> Tensor:
        crystal_features = self.crystal_encoder(data)
        hidden = self.input_projection(crystal_features).reshape(
            -1, self.config.connectome_nodes, self.config.connectome_channels
        )
        for block in self.connectome_blocks:
            hidden = block(
                hidden, self.connectome_edge_index, self.connectome_edge_weight
            )
        attention = torch.softmax(self.attention(hidden).squeeze(-1), dim=1)
        pooled = torch.sum(hidden * attention.unsqueeze(-1), dim=1)
        skip = F.silu(self.global_skip(crystal_features))
        return self.output(torch.cat([pooled, skip], dim=-1)).view(-1)

    def topology_diagnostics(self) -> List[Dict[str, float]]:
        return [
            {
                "block": float(index),
                "topology_gate": float(block.topology_gate().detach().cpu()),
                "edge_weight_power": float(block.edge_weight_power().detach().cpu()),
            }
            for index, block in enumerate(self.connectome_blocks, start=1)
        ]


def summarize_topology_diagnostics(model: FlyConnectome) -> Dict[str, Any]:
    diagnostics = model.topology_diagnostics()
    gates = [item["topology_gate"] for item in diagnostics]
    powers = [item["edge_weight_power"] for item in diagnostics]
    return {
        "topology_gate_mean": float(np.mean(gates)),
        "edge_weight_power_mean": float(np.mean(powers)),
        "topology_gate_by_block": json.dumps(gates),
        "edge_weight_power_by_block": json.dumps(powers),
    }


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
    connectome_topology: Dict[str, Tensor],
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
        "connectome_topology": {
            key: value.detach().cpu() for key, value in connectome_topology.items()
        },
        "connectome_source": str(Path(CONNECTOME_FEATHER)),
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
) -> Tuple[FlyConnectome, GraphConfig, TargetNormalizer, Dict[str, Any]]:
    """Load a trained model and its preprocessing metadata."""

    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint = safe_torch_load(checkpoint_path, map_location=device)
    model_config = ModelConfig(**checkpoint["model_config"])
    graph_config = GraphConfig(**checkpoint["graph_config"])
    normalizer = TargetNormalizer(**checkpoint["normalization"])
    topology = checkpoint["connectome_topology"]
    model = FlyConnectome(
        model_config,
        topology["edge_index"],
        topology["edge_weight"],
    ).to(device)
    model.load_state_dict(checkpoint["model_state_dict"], strict=True)
    model.eval()
    return model, graph_config, normalizer, checkpoint


def fit_model(
    loaders: Dict[str, DataLoader],
    model_config: ModelConfig,
    graph_config: GraphConfig,
    connectome_topology: Dict[str, Tensor],
    normalizer: TargetNormalizer,
    device: torch.device,
    learning_rate: float,
    weight_decay: float,
    max_epochs: int,
    patience: int,
    checkpoint_path: Optional[Union[str, Path]],
    logger: Optional[logging.Logger],
    trial: Any = None,
    run_label: str = "Main",
) -> Tuple[float, int, List[Dict[str, float]]]:

    set_global_seed(SEED)
    model = FlyConnectome(
        model_config,
        connectome_topology["edge_index"],
        connectome_topology["edge_weight"],
    ).to(device)
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
        "run_label": run_label,
        "connectome_steps": model_config.connectome_steps,
        "topology_gate_init": model_config.topology_gate_init,
        "edge_weight_power_init": model_config.edge_weight_power_init,
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
                    connectome_topology,
                    normalizer,
                    best_epoch,
                    best_rmse,
                    training_parameters,
                )
        else:
            epochs_without_improvement += 1

        if logger is not None:
            logger.info(
                "[%s] Epoch %04d | train_MSE=%.8f | train_RMSE=%.6f eV | "
                "val_MSE=%.8f | val_RMSE=%.6f eV | lr=%.3e | %.2f s%s",
                run_label,
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
            trial.report(best_rmse, step=epoch)
            if trial.should_prune():
                import optuna

                raise optuna.TrialPruned()
        if epochs_without_improvement >= patience:
            if logger is not None:
                logger.info(
                    "[%s] Early stopping: validation RMSE did not improve for %d epochs.",
                    run_label,
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
    connectome_topology: Dict[str, Tensor],
    normalizer: TargetNormalizer,
    device: torch.device,
    logger: logging.Logger,
) -> Dict[str, Any]:

    try:
        import optuna
    except ImportError as exc:
        raise ImportError("Install Optuna before setting USE_OPTUNA=True.") from exc

    def objective(trial: Any) -> float:
        hidden = trial.suggest_categorical("hidden_features", [128, 192, 256])
        params = {
            "batch_size": trial.suggest_categorical("batch_size", [16, 32, 64]),
            "learning_rate": trial.suggest_float("learning_rate", 1e-4, 5e-3, log=True),
            "weight_decay": trial.suggest_float("weight_decay", 1e-8, 1e-3, log=True),
            "crystal_layers": trial.suggest_int("crystal_layers", 2, 5),
            "hidden_features": hidden,
            "connectome_channels": trial.suggest_categorical(
                "connectome_channels", [4, 8, 12, 16]
            ),
            "connectome_steps": trial.suggest_int("connectome_steps", 1, 3),
            "dropout": trial.suggest_categorical("dropout", [0.0, 0.1, 0.2]),
        }
        trial_config = ModelConfig(
            atom_input_features=graph_config.max_atomic_number,
            edge_cutoff=graph_config.cutoff,
            hidden_features=params["hidden_features"],
            crystal_layers=params["crystal_layers"],
            connectome_nodes=int(connectome_topology["body_ids"].numel()),
            connectome_channels=params["connectome_channels"],
            connectome_steps=params["connectome_steps"],
            dropout=params["dropout"],
        )
        loaders = make_loaders(
            splits, graph_config, normalizer, batch_size=params["batch_size"]
        )
        best_rmse, _, _ = fit_model(
            loaders=loaders,
            model_config=trial_config,
            graph_config=graph_config,
            connectome_topology=connectome_topology,
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


def _save_ablation_outputs(
    summary_rows: Sequence[Dict[str, Any]],
    prediction_frames: Sequence[pd.DataFrame],
    history_frames: Sequence[pd.DataFrame],
) -> Dict[str, Path]:
    summary = pd.DataFrame(summary_rows)
    predictions = pd.concat(prediction_frames, ignore_index=True)
    histories = pd.concat(history_frames, ignore_index=True)

    metrics_path = Path(TABLE_DIR) / f"{MODEL_NAME}_ablation_metrics.dat"
    predictions_path = Path(DAT_DIR) / f"{MODEL_NAME}_ablation_predictions.dat"
    curves_path = Path(DAT_DIR) / f"{MODEL_NAME}_ablation_rmse_curves.dat"
    summary.to_csv(metrics_path, sep="\t", index=False, float_format="%.8f")
    predictions.to_csv(predictions_path, sep="\t", index=False, float_format="%.8f")
    histories.to_csv(curves_path, sep="\t", index=False, float_format="%.8f")

    variant_order = list(dict.fromkeys(summary["variant_key"].tolist()))
    test = (
        summary.loc[summary["split"] == "test"]
        .set_index("variant_key")
        .loc[variant_order]
        .reset_index()
    )
    plot_labels = [
        label.replace(" ", "\n", 1) if len(label) > 14 else label
        for label in test["variant"].tolist()
    ]
    x_positions = np.arange(len(test))
    fig, axes = plt.subplots(1, 3, figsize=(15.8, 4.8), constrained_layout=True)
    metric_specs = [
        ("mae_eV", "MAE (eV)"),
        ("rmse_eV", "RMSE (eV)"),
        ("r2", "$R^2$"),
    ]
    for axis, (column, ylabel) in zip(axes, metric_specs):
        axis.bar(
            x_positions,
            test[column].to_numpy(dtype=float),
            color=COLORS[: len(test)],
            alpha=0.86,
        )
        axis.set_xticks(x_positions, plot_labels, fontsize=TICK_FONTSIZE - 1)
        axis.set_ylabel(ylabel, fontsize=LABEL_FONTSIZE)
        axis.set_title(f"Test {ylabel}", fontsize=TITLE_FONTSIZE)
        _style_axis(axis)
        axis.grid(False, axis="x")
    figure_metrics_path = Path(FIGURE_DIR) / f"{MODEL_NAME}_ablation_metrics.jpg"
    fig.savefig(figure_metrics_path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(7.6, 5.6), constrained_layout=True)
    for index, variant_key in enumerate(variant_order):
        frame = histories.loc[histories["variant_key"] == variant_key]
        axis.plot(
            frame["epoch"],
            frame["val_rmse_eV"],
            linewidth=1.5,
            color=COLORS[index % len(COLORS)],
            label=str(frame["variant"].iloc[0]),
        )
    axis.set_xlabel("Epoch", fontsize=LABEL_FONTSIZE)
    axis.set_ylabel("Validation RMSE (eV)", fontsize=LABEL_FONTSIZE)
    axis.set_title(f"{MODEL_NAME}: topology ablations", fontsize=TITLE_FONTSIZE)
    axis.legend(fontsize=LEGEND_FONTSIZE)
    _style_axis(axis)
    figure_curves_path = Path(FIGURE_DIR) / f"{MODEL_NAME}_ablation_val_rmse.jpg"
    fig.savefig(figure_curves_path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(fig)
    return {
        "metrics_table": metrics_path,
        "predictions_dat": predictions_path,
        "curves_dat": curves_path,
        "metrics_figure": figure_metrics_path,
        "curves_figure": figure_curves_path,
    }


def run_ablation_experiments(
    splits: Dict[str, List[Dict[str, Any]]],
    model_config: ModelConfig,
    graph_config: GraphConfig,
    connectome_topology: Dict[str, Tensor],
    normalizer: TargetNormalizer,
    device: torch.device,
    selected: Dict[str, Any],
    full_results: Dict[str, Dict[str, Any]],
    full_history: Sequence[Dict[str, float]],
    full_best_epoch: int,
    full_best_rmse: float,
    full_train_seconds: float,
    full_topology_diagnostics: Dict[str, Any],
    logger: logging.Logger,
) -> Dict[str, Path]:
    variants = build_ablation_topologies(connectome_topology, logger)
    summary_rows: List[Dict[str, Any]] = []
    prediction_frames: List[pd.DataFrame] = []
    history_frames: List[pd.DataFrame] = []

    parameter_model = FlyConnectome(
        model_config,
        connectome_topology["edge_index"],
        connectome_topology["edge_weight"],
    )
    total_parameters = sum(p.numel() for p in parameter_model.parameters())
    trainable_parameters = sum(
        p.numel() for p in parameter_model.parameters() if p.requires_grad
    )
    del parameter_model

    def record_variant(
        variant_key: str,
        results: Dict[str, Dict[str, Any]],
        history: Sequence[Dict[str, float]],
        best_epoch: int,
        best_rmse: float,
        train_seconds: float,
        topology_diagnostics: Dict[str, Any],
    ) -> None:
        variant = variants[variant_key]
        topology = variant["topology"]
        for split_name in ["train", "val", "test"]:
            metrics = results[split_name]["metrics"]
            summary_rows.append(
                {
                    "variant_key": variant_key,
                    "variant": variant["label"],
                    "description": variant["description"],
                    "split": split_name,
                    "n_samples": len(results[split_name]["y_true"]),
                    "topology_nodes": int(topology["body_ids"].numel()),
                    "topology_edges": int(topology["edge_weight"].numel()),
                    "total_parameters": total_parameters,
                    "trainable_parameters": trainable_parameters,
                    "best_epoch": int(best_epoch),
                    "best_val_rmse_eV": float(best_rmse),
                    "training_seconds": float(train_seconds),
                    **topology_diagnostics,
                    "mae_eV": float(metrics["mae"]),
                    "rmse_eV": float(metrics["rmse"]),
                    "r2": float(metrics["r2"]),
                }
            )
            prediction = prediction_frame(results[split_name], split_name)
            prediction.insert(0, "variant", variant["label"])
            prediction.insert(0, "variant_key", variant_key)
            prediction_frames.append(prediction)
        history_frame = pd.DataFrame(history).copy()
        history_frame["epoch"] = history_frame["epoch"].astype(int)
        history_frame.insert(0, "variant", variant["label"])
        history_frame.insert(0, "variant_key", variant_key)
        history_frames.append(history_frame)

    record_variant(
        "full_malecns",
        full_results,
        full_history,
        full_best_epoch,
        full_best_rmse,
        full_train_seconds,
        full_topology_diagnostics,
    )

    for variant_key in [
        "self_only",
        "random_topology",
        "degree_preserving",
        "unweighted_malecns",
    ]:
        variant = variants[variant_key]
        topology = variant["topology"]
        checkpoint_path = (
            Path(ABLATION_CHECKPOINT_DIR) / f"{MODEL_NAME}_{variant_key}_best.pt"
        )
        logger.info(
            "Ablation start | %s | %s | checkpoint=%s",
            variant["label"],
            variant["description"],
            checkpoint_path.resolve(),
        )
        variant_loaders = make_loaders(
            splits,
            graph_config,
            normalizer,
            batch_size=int(selected["batch_size"]),
        )
        start = time.perf_counter()
        best_rmse, best_epoch, history = fit_model(
            loaders=variant_loaders,
            model_config=model_config,
            graph_config=graph_config,
            connectome_topology=topology,
            normalizer=normalizer,
            device=device,
            learning_rate=float(selected["learning_rate"]),
            weight_decay=float(selected["weight_decay"]),
            max_epochs=MAX_EPOCHS,
            patience=PATIENCE,
            checkpoint_path=checkpoint_path,
            logger=logger,
            run_label=str(variant["label"]),
        )
        train_seconds = time.perf_counter() - start
        model, saved_graph_config, saved_normalizer, checkpoint = load_trained_model(
            checkpoint_path, device
        )
        if saved_graph_config != graph_config:
            raise RuntimeError(
                f"Graph parameters do not match for ablation {variant_key}."
            )
        criterion = nn.MSELoss()
        results = {
            split_name: evaluate_loader(
                model,
                variant_loaders[split_name],
                criterion,
                saved_normalizer,
                device,
            )
            for split_name in ["train", "val", "test"]
        }
        if int(checkpoint["best_epoch"]) != best_epoch:
            raise RuntimeError(f"Best-epoch mismatch for ablation {variant_key}.")
        topology_diagnostics = summarize_topology_diagnostics(model)
        record_variant(
            variant_key,
            results,
            history,
            best_epoch,
            best_rmse,
            train_seconds,
            topology_diagnostics,
        )
        logger.info(
            "Ablation complete | %s | best_epoch=%d | best_val_RMSE=%.6f eV | "
            "test_MAE=%.6f eV | test_RMSE=%.6f eV | test_R2=%.6f | "
            "topology_gate=%s | edge_weight_power=%s | %.3f s (%s)",
            variant["label"],
            best_epoch,
            best_rmse,
            results["test"]["metrics"]["mae"],
            results["test"]["metrics"]["rmse"],
            results["test"]["metrics"]["r2"],
            topology_diagnostics["topology_gate_by_block"],
            topology_diagnostics["edge_weight_power_by_block"],
            train_seconds,
            format_duration(train_seconds),
        )
        del model
        if device.type == "cuda":
            torch.cuda.empty_cache()

    paths = _save_ablation_outputs(summary_rows, prediction_frames, history_frames)
    for name, path in paths.items():
        logger.info("Ablation output | %s=%s", name, path.resolve())
    return paths


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
                atomic_numbers=graph["atomic_numbers"],
                edge_index=graph["edge_index"],
                edge_vec=graph["edge_vec"],
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
    logger.info("MaleCNS Feather: %s", Path(CONNECTOME_FEATHER).resolve())
    logger.info(
        "Excel columns: first=CIF filename; second=%s (%s)", TARGET_NAME, TARGET_UNIT
    )
    logger.info("Target: %s (%s)", TARGET_NAME, TARGET_UNIT)
    logger.info("Random seed: %d", SEED)
    logger.info(
        "Topology ablations enabled: %s; ablation random seed: %d; swap factor: %d",
        RUN_ABLATIONS,
        ABLATION_RANDOM_SEED,
        ABLATION_SWAP_FACTOR,
    )
    logger.info(
        "Adaptive topology propagation: steps=%d; initial gate=%.4f; "
        "initial edge-weight power=%.4f",
        CONNECTOME_STEPS,
        TOPOLOGY_GATE_INIT,
        EDGE_WEIGHT_POWER_INIT,
    )
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
    connectome_topology = build_connectome_topology(
        CONNECTOME_FEATHER,
        CONNECTOME_NODES,
        CONNECTOME_CACHE_PATH,
        logger,
    )
    logger.info(
        "MaleCNS-derived topology: selected neurons=%d, directed edges=%d.",
        int(connectome_topology["body_ids"].numel()),
        int(connectome_topology["edge_weight"].numel()),
    )
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
        "crystal_layers": CRYSTAL_LAYERS,
        "hidden_features": HIDDEN_FEATURES,
        "connectome_channels": CONNECTOME_CHANNELS,
        "connectome_steps": CONNECTOME_STEPS,
        "dropout": DROPOUT,
    }
    logger.info("Optuna enabled: %s", USE_OPTUNA)
    if USE_OPTUNA:
        optuna_result = run_optuna(
            splits,
            graph_config,
            connectome_topology,
            normalizer,
            device,
            logger,
        )
        selected.update(optuna_result["best_params"])

    model_config = ModelConfig(
        atom_input_features=graph_config.max_atomic_number,
        edge_cutoff=graph_config.cutoff,
        hidden_features=int(selected["hidden_features"]),
        crystal_layers=int(selected["crystal_layers"]),
        connectome_nodes=int(connectome_topology["body_ids"].numel()),
        connectome_channels=int(selected["connectome_channels"]),
        connectome_steps=int(selected["connectome_steps"]),
        dropout=float(selected["dropout"]),
    )
    loaders = make_loaders(
        splits, graph_config, normalizer, batch_size=int(selected["batch_size"])
    )
    model_for_log = FlyConnectome(
        model_config,
        connectome_topology["edge_index"],
        connectome_topology["edge_weight"],
    )
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
        connectome_topology=connectome_topology,
        normalizer=normalizer,
        device=device,
        learning_rate=float(selected["learning_rate"]),
        weight_decay=float(selected["weight_decay"]),
        max_epochs=MAX_EPOCHS,
        patience=PATIENCE,
        checkpoint_path=CHECKPOINT_PATH,
        logger=logger,
        run_label="Full MaleCNS",
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
    full_topology_diagnostics = summarize_topology_diagnostics(best_model)
    logger.info(
        "Best-model topology diagnostics | gate_by_block=%s | "
        "edge_weight_power_by_block=%s",
        full_topology_diagnostics["topology_gate_by_block"],
        full_topology_diagnostics["edge_weight_power_by_block"],
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
    logger.info("Best model: %s", Path(CHECKPOINT_PATH).resolve())
    logger.info("Metrics table: %s", metrics_path.resolve())
    logger.info(
        "Evaluation, saving, and plotting time: %.3f s (%s)",
        evaluation_seconds,
        format_duration(evaluation_seconds),
    )

    del best_model
    if device.type == "cuda":
        torch.cuda.empty_cache()

    if RUN_ABLATIONS:
        ablation_start = time.perf_counter()
        logger.info(
            "Starting necessary topology ablations. All variants use the same "
            "data split, model configuration, optimizer settings, and validation-RMSE selection."
        )
        run_ablation_experiments(
            splits=splits,
            model_config=model_config,
            graph_config=graph_config,
            connectome_topology=connectome_topology,
            normalizer=normalizer,
            device=device,
            selected=selected,
            full_results=results,
            full_history=history,
            full_best_epoch=best_epoch,
            full_best_rmse=best_rmse,
            full_train_seconds=train_seconds,
            full_topology_diagnostics=full_topology_diagnostics,
            logger=logger,
        )
        ablation_seconds = time.perf_counter() - ablation_start
        logger.info(
            "Ablation runtime: %.3f s (%s)",
            ablation_seconds,
            format_duration(ablation_seconds),
        )

    total_seconds = time.perf_counter() - total_start
    logger.info(
        "Total runtime: %.3f s (%s)", total_seconds, format_duration(total_seconds)
    )


if __name__ == "__main__":
    main()

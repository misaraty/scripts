from __future__ import annotations

import json
import logging
import math
import os
import random
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from pymatgen.core import Structure
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from torch import Tensor, nn
from torch.utils.data import DataLoader, Dataset


# User configuration
MODEL_NAME = "SchNet"
RUN_VERSION = "v1_USE_COSINE_CUTOFF"  # v1_COSINE_CUTOFF, USE_COSINE_CUTOFF = True

EXCEL_PATH = "./data.xlsx"
CIF_DIR = "./cif"

TARGET_NAME = "bandgap"
TARGET_UNIT = "eV"

SEED = 42
TRAIN_RATIO = 0.80
VAL_RATIO = 0.10
TEST_RATIO = 0.10

CUTOFF = 6.0
MAX_NEIGHBORS: Optional[int] = None
NUM_ATOM_TYPES = 119
HIDDEN_DIM = 128
NUM_FILTERS = 128
NUM_GAUSSIANS = 64
NUM_INTERACTIONS = 6
READOUT = "mean"
USE_COSINE_CUTOFF = True

BATCH_SIZE = 32
MAX_EPOCHS = 300
LEARNING_RATE = 1.0e-3
WEIGHT_DECAY = 1.0e-6
PATIENCE = 60
MIN_DELTA = 1.0e-5
LR_FACTOR = 0.5
LR_PATIENCE = 15
MIN_LR = 1.0e-6
NUM_WORKERS = 0
PIN_MEMORY = True
USE_AMP = True
AMP_DTYPE = "bfloat16"
GRAD_CLIP_NORM = 5.0

USE_OPTUNA = False
OPTUNA_N_TRIALS = 30
OPTUNA_TIMEOUT: Optional[int] = None
OPTUNA_MAX_EPOCHS = 180
OPTUNA_PATIENCE = 30

FIG_DPI = 600
TITLE_FONTSIZE = 14
LABEL_FONTSIZE = 12
TICK_FONTSIZE = 10
LEGEND_FONTSIZE = 10
ANNOTATION_FONTSIZE = 9
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
GRAPH_CACHE_PATH = CACHE_DIR / f"{MODEL_NAME}_graphs.pt"


@dataclass(frozen=True)
class GraphConfig:
    cutoff: float = CUTOFF
    max_neighbors: Optional[int] = MAX_NEIGHBORS


@dataclass(frozen=True)
class ModelConfig:
    num_atom_types: int = NUM_ATOM_TYPES
    hidden_dim: int = HIDDEN_DIM
    num_filters: int = NUM_FILTERS
    num_gaussians: int = NUM_GAUSSIANS
    num_interactions: int = NUM_INTERACTIONS
    cutoff: float = CUTOFF
    readout: str = READOUT
    use_cosine_cutoff: bool = USE_COSINE_CUTOFF


@dataclass
class CrystalGraph:
    cif_name: str
    atomic_numbers: Tensor
    edge_index: Tensor
    edge_distances: Tensor
    target: float


@dataclass
class CrystalBatch:
    cif_names: List[str]
    atomic_numbers: Tensor
    edge_index: Tensor
    edge_distances: Tensor
    batch_index: Tensor
    targets: Tensor

    def to(self, device: torch.device) -> "CrystalBatch":
        return CrystalBatch(
            cif_names=self.cif_names,
            atomic_numbers=self.atomic_numbers.to(device, non_blocking=True),
            edge_index=self.edge_index.to(device, non_blocking=True),
            edge_distances=self.edge_distances.to(device, non_blocking=True),
            batch_index=self.batch_index.to(device, non_blocking=True),
            targets=self.targets.to(device, non_blocking=True),
        )


@dataclass(frozen=True)
class TargetScaler:
    mean: float
    std: float

    def transform_tensor(self, values: Tensor) -> Tensor:
        return (values - self.mean) / self.std

    def inverse_numpy(self, values: np.ndarray) -> np.ndarray:
        return values * self.std + self.mean


class CrystalGraphDataset(Dataset[CrystalGraph]):
    def __init__(self, graphs: Sequence[CrystalGraph]) -> None:
        self.graphs = list(graphs)

    def __len__(self) -> int:
        return len(self.graphs)

    def __getitem__(self, index: int) -> CrystalGraph:
        return self.graphs[index]


def create_output_directories() -> None:
    for directory in (
        OUTPUT_ROOT,
        FIGURE_DIR,
        DAT_DIR,
        TABLE_DIR,
        LOG_DIR,
        SPLIT_DIR,
        CACHE_DIR,
    ):
        directory.mkdir(parents=True, exist_ok=True)


def setup_logger() -> logging.Logger:
    logger = logging.getLogger(f"{MODEL_NAME}_{RUN_VERSION}")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    file_handler = logging.FileHandler(
        LOG_DIR / f"{MODEL_NAME}_training.log", mode="w", encoding="utf-8"
    )
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
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
    total = max(0, int(round(seconds)))
    hours, remainder = divmod(total, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def get_device(logger: logging.Logger) -> torch.device:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info("Device: %s", device)
    logger.info("PyTorch version: %s", torch.__version__)
    if device.type == "cuda":
        props = torch.cuda.get_device_properties(device)
        logger.info("GPU: %s", props.name)
        logger.info("GPU memory: %.2f GiB", props.total_memory / 1024**3)
        logger.info("CUDA runtime: %s", torch.version.cuda)
    return device


def normalize_cif_name(raw_name: Any) -> str:
    if pd.isna(raw_name):
        raise ValueError("empty CIF name")
    if isinstance(raw_name, (int, np.integer)):
        name = str(int(raw_name))
    elif isinstance(raw_name, (float, np.floating)) and float(raw_name).is_integer():
        name = str(int(raw_name))
    else:
        name = str(raw_name).strip()
    if not name:
        raise ValueError("empty CIF name")
    if not name.lower().endswith(".cif"):
        name = f"{name}.cif"
    return name


def load_excel_records(logger: logging.Logger) -> List[Tuple[str, float]]:
    excel_path = Path(EXCEL_PATH)
    if not excel_path.is_file():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    frame = pd.read_excel(excel_path)
    if frame.shape[1] < 2:
        raise ValueError("The Excel file must contain at least two columns.")

    records: List[Tuple[str, float]] = []
    seen: set[str] = set()
    for row_number, (raw_name, raw_target) in enumerate(
        frame.iloc[:, :2].itertuples(index=False, name=None), start=2
    ):
        try:
            cif_name = normalize_cif_name(raw_name)
            target = float(raw_target)
            if not np.isfinite(target):
                raise ValueError("target is NaN or infinite")
            if cif_name in seen:
                raise ValueError("duplicate CIF name")
            cif_path = Path(CIF_DIR) / cif_name
            if not cif_path.is_file():
                raise FileNotFoundError(f"CIF file not found: {cif_path}")
            records.append((cif_name, target))
            seen.add(cif_name)
        except Exception as exc:
            logger.warning("Skipped Excel row %d: %s", row_number, exc)

    if len(records) < 10:
        raise ValueError(
            f"Only {len(records)} valid rows were found; at least 10 are required."
        )
    logger.info("Valid Excel rows: %d", len(records))
    return records


def _limit_neighbors(
    center_indices: np.ndarray,
    neighbor_indices: np.ndarray,
    distances: np.ndarray,
    max_neighbors: Optional[int],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    if max_neighbors is None:
        return center_indices, neighbor_indices, distances
    keep: List[int] = []
    for center in np.unique(center_indices):
        candidates = np.flatnonzero(center_indices == center)
        order = candidates[np.argsort(distances[candidates], kind="stable")]
        keep.extend(order[:max_neighbors].tolist())
    keep_array = np.asarray(sorted(keep), dtype=np.int64)
    return (
        center_indices[keep_array],
        neighbor_indices[keep_array],
        distances[keep_array],
    )


def structure_to_graph(
    structure: Structure,
    cif_name: str,
    target: float,
    graph_config: GraphConfig,
) -> CrystalGraph:
    atomic_numbers = np.asarray(
        [int(site.specie.Z) for site in structure], dtype=np.int64
    )
    if atomic_numbers.size == 0:
        raise ValueError("structure contains no atoms")
    if np.any(atomic_numbers <= 0) or np.any(atomic_numbers >= NUM_ATOM_TYPES):
        raise ValueError("atomic number is outside the embedding range")

    centers, neighbors, _, distances = structure.get_neighbor_list(
        r=graph_config.cutoff,
        numerical_tol=1.0e-8,
        exclude_self=True,
    )
    centers = np.asarray(centers, dtype=np.int64)
    neighbors = np.asarray(neighbors, dtype=np.int64)
    distances = np.asarray(distances, dtype=np.float32)
    valid = np.isfinite(distances) & (distances > 1.0e-8)
    centers, neighbors, distances = (
        centers[valid],
        neighbors[valid],
        distances[valid],
    )
    centers, neighbors, distances = _limit_neighbors(
        centers, neighbors, distances, graph_config.max_neighbors
    )
    if distances.size == 0:
        raise ValueError("no periodic neighbors were found inside the cutoff")

    edge_index = np.stack([centers, neighbors], axis=0)
    return CrystalGraph(
        cif_name=cif_name,
        atomic_numbers=torch.from_numpy(atomic_numbers),
        edge_index=torch.from_numpy(edge_index),
        edge_distances=torch.from_numpy(distances),
        target=float(target),
    )


def graph_cache_signature(
    records: Sequence[Tuple[str, float]], config: GraphConfig
) -> Dict[str, Any]:
    source_items = []
    for cif_name, target in records:
        path = Path(CIF_DIR) / cif_name
        stat = path.stat()
        source_items.append((cif_name, float(target), stat.st_size, stat.st_mtime_ns))
    return {
        "excel_path": str(Path(EXCEL_PATH).resolve()),
        "cif_dir": str(Path(CIF_DIR).resolve()),
        "records": source_items,
        "graph_config": asdict(config),
    }


def load_or_build_graphs(
    records: Sequence[Tuple[str, float]],
    graph_config: GraphConfig,
    logger: logging.Logger,
) -> List[CrystalGraph]:
    signature = graph_cache_signature(records, graph_config)
    if GRAPH_CACHE_PATH.is_file():
        try:
            cached = torch.load(
                GRAPH_CACHE_PATH, map_location="cpu", weights_only=False
            )
            if cached.get("signature") == signature:
                graphs = cached["graphs"]
                logger.info("Loaded %d graphs from cache.", len(graphs))
                return graphs
            logger.info("Graph cache signature changed; rebuilding graphs.")
        except Exception as exc:
            logger.warning("Graph cache could not be loaded: %s", exc)

    graphs: List[CrystalGraph] = []
    start = time.perf_counter()
    for index, (cif_name, target) in enumerate(records, start=1):
        try:
            structure = Structure.from_file(Path(CIF_DIR) / cif_name)
            graphs.append(structure_to_graph(structure, cif_name, target, graph_config))
        except Exception as exc:
            logger.warning("Skipped %s during graph construction: %s", cif_name, exc)
        if index % 250 == 0 or index == len(records):
            logger.info("Graph construction: %d/%d", index, len(records))

    if len(graphs) < 10:
        raise ValueError(
            f"Only {len(graphs)} valid crystal graphs were built; at least 10 are required."
        )
    torch.save({"signature": signature, "graphs": graphs}, GRAPH_CACHE_PATH)
    elapsed = time.perf_counter() - start
    logger.info(
        "Built %d graphs in %.3f s (%s).",
        len(graphs),
        elapsed,
        format_duration(elapsed),
    )
    return graphs


def create_or_load_split(
    graphs: Sequence[CrystalGraph], logger: logging.Logger
) -> Dict[str, str]:
    names = [graph.cif_name for graph in graphs]
    targets = {graph.cif_name: graph.target for graph in graphs}
    name_set = set(names)

    if SPLIT_PATH.is_file():
        try:
            frame = pd.read_csv(SPLIT_PATH)
            required = {"cif_name", "target", "split"}
            if not required.issubset(frame.columns):
                raise ValueError("split file has missing columns")
            if set(frame["cif_name"].astype(str)) != name_set:
                raise ValueError("split file does not match the valid graph set")
            if not set(frame["split"]).issubset({"train", "val", "test"}):
                raise ValueError("split file contains invalid split labels")
            for row in frame.itertuples(index=False):
                if not math.isclose(
                    float(row.target),
                    targets[str(row.cif_name)],
                    rel_tol=0.0,
                    abs_tol=1.0e-10,
                ):
                    raise ValueError(
                        "split targets do not match the current Excel file"
                    )
            split_map = dict(zip(frame["cif_name"].astype(str), frame["split"]))
            if all(value in split_map.values() for value in ("train", "val", "test")):
                logger.info("Reused split file: %s", SPLIT_PATH)
                return split_map
            raise ValueError("one or more split subsets are empty")
        except Exception as exc:
            logger.warning("Existing split file was not reused: %s", exc)

    if not math.isclose(TRAIN_RATIO + VAL_RATIO + TEST_RATIO, 1.0, abs_tol=1.0e-8):
        raise ValueError("TRAIN_RATIO + VAL_RATIO + TEST_RATIO must equal 1.0.")
    train_names, remainder_names = train_test_split(
        names,
        test_size=VAL_RATIO + TEST_RATIO,
        random_state=SEED,
        shuffle=True,
    )
    test_fraction = TEST_RATIO / (VAL_RATIO + TEST_RATIO)
    val_names, test_names = train_test_split(
        remainder_names,
        test_size=test_fraction,
        random_state=SEED,
        shuffle=True,
    )
    split_map = {name: "train" for name in train_names}
    split_map.update({name: "val" for name in val_names})
    split_map.update({name: "test" for name in test_names})
    rows = [
        {"cif_name": name, "target": targets[name], "split": split_map[name]}
        for name in names
    ]
    pd.DataFrame(rows).to_csv(SPLIT_PATH, index=False)
    logger.info("Created split file: %s", SPLIT_PATH)
    return split_map


def partition_graphs(
    graphs: Sequence[CrystalGraph], split_map: Dict[str, str]
) -> Tuple[List[CrystalGraph], List[CrystalGraph], List[CrystalGraph]]:
    subsets: Dict[str, List[CrystalGraph]] = {"train": [], "val": [], "test": []}
    for graph in graphs:
        subsets[split_map[graph.cif_name]].append(graph)
    if any(len(subsets[key]) == 0 for key in subsets):
        raise ValueError(
            "The train, validation, and test subsets must all be non-empty."
        )
    return subsets["train"], subsets["val"], subsets["test"]


def collate_crystal_graphs(graphs: Sequence[CrystalGraph]) -> CrystalBatch:
    atomic_numbers: List[Tensor] = []
    edge_indices: List[Tensor] = []
    edge_distances: List[Tensor] = []
    batch_indices: List[Tensor] = []
    targets: List[float] = []
    cif_names: List[str] = []
    node_offset = 0
    for graph_index, graph in enumerate(graphs):
        num_atoms = graph.atomic_numbers.numel()
        atomic_numbers.append(graph.atomic_numbers)
        edge_indices.append(graph.edge_index + node_offset)
        edge_distances.append(graph.edge_distances)
        batch_indices.append(torch.full((num_atoms,), graph_index, dtype=torch.long))
        targets.append(graph.target)
        cif_names.append(graph.cif_name)
        node_offset += num_atoms
    return CrystalBatch(
        cif_names=cif_names,
        atomic_numbers=torch.cat(atomic_numbers, dim=0),
        edge_index=torch.cat(edge_indices, dim=1),
        edge_distances=torch.cat(edge_distances, dim=0),
        batch_index=torch.cat(batch_indices, dim=0),
        targets=torch.tensor(targets, dtype=torch.float32),
    )


def make_loader(
    graphs: Sequence[CrystalGraph],
    batch_size: int,
    shuffle: bool,
    device: torch.device,
) -> DataLoader:
    generator = torch.Generator()
    generator.manual_seed(SEED)
    return DataLoader(
        CrystalGraphDataset(graphs),
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=NUM_WORKERS,
        pin_memory=PIN_MEMORY and device.type == "cuda",
        collate_fn=collate_crystal_graphs,
        generator=generator,
        persistent_workers=NUM_WORKERS > 0,
    )


class ShiftedSoftplus(nn.Module):
    def forward(self, inputs: Tensor) -> Tensor:
        return torch.nn.functional.softplus(inputs) - math.log(2.0)


class GaussianRBF(nn.Module):
    def __init__(self, cutoff: float, num_gaussians: int) -> None:
        super().__init__()
        centers = torch.linspace(0.0, cutoff, num_gaussians)
        spacing = float(centers[1] - centers[0]) if num_gaussians > 1 else cutoff
        self.register_buffer("centers", centers)
        self.gamma = 1.0 / max(spacing, 1.0e-12)

    def forward(self, distances: Tensor) -> Tensor:
        difference = distances.unsqueeze(-1) - self.centers
        return torch.exp(-self.gamma * difference.square())


class CosineCutoff(nn.Module):
    def __init__(self, cutoff: float) -> None:
        super().__init__()
        self.cutoff = float(cutoff)

    def forward(self, distances: Tensor) -> Tensor:
        values = 0.5 * (torch.cos(math.pi * distances / self.cutoff) + 1.0)
        return values * (distances < self.cutoff).to(values.dtype)


class ContinuousFilterConvolution(nn.Module):
    def __init__(
        self,
        hidden_dim: int,
        num_filters: int,
        num_gaussians: int,
        cutoff: float,
        use_cosine_cutoff: bool,
    ) -> None:
        super().__init__()
        activation = ShiftedSoftplus()
        self.atom_projection = nn.Linear(hidden_dim, num_filters, bias=False)
        self.filter_network = nn.Sequential(
            nn.Linear(num_gaussians, num_filters),
            activation,
            nn.Linear(num_filters, num_filters),
            ShiftedSoftplus(),
        )
        self.output_network = nn.Sequential(
            nn.Linear(num_filters, hidden_dim),
            ShiftedSoftplus(),
            nn.Linear(hidden_dim, hidden_dim),
        )
        self.cutoff_network = CosineCutoff(cutoff)
        self.use_cosine_cutoff = use_cosine_cutoff

    def forward(
        self,
        atom_features: Tensor,
        edge_index: Tensor,
        edge_distances: Tensor,
        rbf: Tensor,
    ) -> Tensor:
        center, neighbor = edge_index
        projected = self.atom_projection(atom_features)
        filters = self.filter_network(rbf)
        if self.use_cosine_cutoff:
            filters = filters * self.cutoff_network(edge_distances).unsqueeze(-1)
        messages = projected[neighbor] * filters
        aggregated = atom_features.new_zeros((atom_features.size(0), messages.size(1)))
        aggregated.index_add_(0, center, messages)
        return self.output_network(aggregated)


class SchNetInteraction(nn.Module):
    def __init__(
        self,
        hidden_dim: int,
        num_filters: int,
        num_gaussians: int,
        cutoff: float,
        use_cosine_cutoff: bool,
    ) -> None:
        super().__init__()
        self.convolution = ContinuousFilterConvolution(
            hidden_dim=hidden_dim,
            num_filters=num_filters,
            num_gaussians=num_gaussians,
            cutoff=cutoff,
            use_cosine_cutoff=use_cosine_cutoff,
        )

    def forward(
        self,
        atom_features: Tensor,
        edge_index: Tensor,
        edge_distances: Tensor,
        rbf: Tensor,
    ) -> Tensor:
        return atom_features + self.convolution(
            atom_features, edge_index, edge_distances, rbf
        )


def segment_pool(values: Tensor, batch_index: Tensor, mode: str) -> Tensor:
    num_graphs = int(batch_index.max().item()) + 1
    pooled = values.new_zeros((num_graphs, values.size(-1)))
    pooled.index_add_(0, batch_index, values)
    if mode == "mean":
        counts = torch.bincount(batch_index, minlength=num_graphs).to(values.dtype)
        pooled = pooled / counts.clamp_min(1.0).unsqueeze(-1)
    elif mode != "sum":
        raise ValueError(f"Unsupported readout mode: {mode}")
    return pooled


class SchNetBandgap(nn.Module):
    def __init__(self, config: ModelConfig) -> None:
        super().__init__()
        if config.readout not in {"mean", "sum"}:
            raise ValueError("readout must be 'mean' or 'sum'")
        self.config = config
        self.atom_embedding = nn.Embedding(config.num_atom_types, config.hidden_dim)
        self.rbf = GaussianRBF(config.cutoff, config.num_gaussians)
        self.interactions = nn.ModuleList(
            [
                SchNetInteraction(
                    hidden_dim=config.hidden_dim,
                    num_filters=config.num_filters,
                    num_gaussians=config.num_gaussians,
                    cutoff=config.cutoff,
                    use_cosine_cutoff=config.use_cosine_cutoff,
                )
                for _ in range(config.num_interactions)
            ]
        )
        self.atomwise_output = nn.Sequential(
            nn.Linear(config.hidden_dim, config.hidden_dim // 2),
            ShiftedSoftplus(),
            nn.Linear(config.hidden_dim // 2, 1),
        )
        self.reset_parameters()

    def reset_parameters(self) -> None:
        nn.init.normal_(self.atom_embedding.weight, mean=0.0, std=1.0)
        for module in self.modules():
            if isinstance(module, nn.Linear):
                nn.init.xavier_uniform_(module.weight)
                if module.bias is not None:
                    nn.init.zeros_(module.bias)
        final = self.atomwise_output[-1]
        nn.init.zeros_(final.weight)
        nn.init.zeros_(final.bias)

    def forward(self, batch: CrystalBatch) -> Tensor:
        atom_features = self.atom_embedding(batch.atomic_numbers)
        rbf = self.rbf(batch.edge_distances)
        for interaction in self.interactions:
            atom_features = interaction(
                atom_features, batch.edge_index, batch.edge_distances, rbf
            )
        atom_outputs = self.atomwise_output(atom_features)
        graph_outputs = segment_pool(
            atom_outputs, batch.batch_index, self.config.readout
        )
        return graph_outputs.view(-1)


def count_parameters(model: nn.Module) -> Tuple[int, int]:
    total = sum(parameter.numel() for parameter in model.parameters())
    trainable = sum(
        parameter.numel() for parameter in model.parameters() if parameter.requires_grad
    )
    return total, trainable


def make_target_scaler(train_graphs: Sequence[CrystalGraph]) -> TargetScaler:
    targets = np.asarray([graph.target for graph in train_graphs], dtype=np.float64)
    mean = float(targets.mean())
    std = float(targets.std(ddof=0))
    if not np.isfinite(std) or std < 1.0e-12:
        std = 1.0
    return TargetScaler(mean=mean, std=std)


def autocast_context(device: torch.device):
    enabled = USE_AMP and device.type == "cuda"
    dtype = torch.bfloat16 if AMP_DTYPE.lower() == "bfloat16" else torch.float16
    return torch.autocast(device_type=device.type, dtype=dtype, enabled=enabled)


def train_one_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    criterion: nn.Module,
    scaler: TargetScaler,
    device: torch.device,
    grad_scaler: torch.amp.GradScaler,
) -> float:
    model.train()
    loss_sum = 0.0
    sample_count = 0
    for batch in loader:
        batch = batch.to(device)
        optimizer.zero_grad(set_to_none=True)
        with autocast_context(device):
            predictions = model(batch).view(-1)
            targets = scaler.transform_tensor(batch.targets.view(-1))
            loss = criterion(predictions, targets)
        if not torch.isfinite(loss):
            raise FloatingPointError("A non-finite training loss was encountered.")
        grad_scaler.scale(loss).backward()
        if GRAD_CLIP_NORM is not None and GRAD_CLIP_NORM > 0:
            grad_scaler.unscale_(optimizer)
            torch.nn.utils.clip_grad_norm_(model.parameters(), GRAD_CLIP_NORM)
        grad_scaler.step(optimizer)
        grad_scaler.update()
        batch_size = batch.targets.numel()
        loss_sum += float(loss.detach()) * batch_size
        sample_count += batch_size
    return loss_sum / max(sample_count, 1)


@torch.inference_mode()
def predict_loader(
    model: nn.Module,
    loader: DataLoader,
    scaler: TargetScaler,
    device: torch.device,
) -> Tuple[List[str], np.ndarray, np.ndarray]:
    model.eval()
    all_names: List[str] = []
    all_targets: List[np.ndarray] = []
    all_predictions: List[np.ndarray] = []
    for batch in loader:
        batch = batch.to(device)
        with autocast_context(device):
            scaled_predictions = model(batch).view(-1)
        predictions = scaler.inverse_numpy(
            scaled_predictions.float().cpu().numpy().astype(np.float64)
        )
        all_names.extend(batch.cif_names)
        all_targets.append(batch.targets.float().cpu().numpy().astype(np.float64))
        all_predictions.append(predictions)
    return (
        all_names,
        np.concatenate(all_targets),
        np.concatenate(all_predictions),
    )


def calculate_metrics(targets: np.ndarray, predictions: np.ndarray) -> Dict[str, float]:
    if targets.shape != predictions.shape:
        raise ValueError("Target and prediction shapes do not match.")
    metrics = {
        "mae": float(mean_absolute_error(targets, predictions)),
        "rmse": float(np.sqrt(mean_squared_error(targets, predictions))),
        "r2": (
            float(r2_score(targets, predictions)) if targets.size >= 2 else float("nan")
        ),
    }
    return metrics


def checkpoint_payload(
    model: nn.Module,
    optimizer: torch.optim.Optimizer,
    scheduler: Any,
    model_config: ModelConfig,
    graph_config: GraphConfig,
    scaler: TargetScaler,
    epoch: int,
    best_val_rmse: float,
    hyperparameters: Dict[str, Any],
) -> Dict[str, Any]:
    return {
        "model_name": MODEL_NAME,
        "run_version": RUN_VERSION,
        "model_state_dict": model.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "scheduler_state_dict": scheduler.state_dict(),
        "model_config": asdict(model_config),
        "graph_config": asdict(graph_config),
        "target_scaler": asdict(scaler),
        "best_epoch": epoch,
        "best_val_rmse": best_val_rmse,
        "seed": SEED,
        "target_name": TARGET_NAME,
        "target_unit": TARGET_UNIT,
        "hyperparameters": hyperparameters,
    }


def train_model(
    train_graphs: Sequence[CrystalGraph],
    val_graphs: Sequence[CrystalGraph],
    model_config: ModelConfig,
    graph_config: GraphConfig,
    hyperparameters: Dict[str, Any],
    device: torch.device,
    logger: logging.Logger,
    save_checkpoint: bool,
    max_epochs: int,
    patience: int,
) -> Tuple[nn.Module, TargetScaler, List[Dict[str, float]], int, float, bool]:
    batch_size = int(hyperparameters["batch_size"])
    train_loader = make_loader(train_graphs, batch_size, True, device)
    train_eval_loader = make_loader(train_graphs, batch_size, False, device)
    val_loader = make_loader(val_graphs, batch_size, False, device)
    target_scaler = make_target_scaler(train_graphs)
    model = SchNetBandgap(model_config).to(device)
    criterion = nn.MSELoss()
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=float(hyperparameters["learning_rate"]),
        weight_decay=float(hyperparameters["weight_decay"]),
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="min",
        factor=LR_FACTOR,
        patience=LR_PATIENCE,
        min_lr=MIN_LR,
    )
    use_grad_scaler = (
        USE_AMP and device.type == "cuda" and AMP_DTYPE.lower() == "float16"
    )
    grad_scaler = torch.amp.GradScaler("cuda", enabled=use_grad_scaler)

    history: List[Dict[str, float]] = []
    best_val_rmse = float("inf")
    best_epoch = 0
    epochs_without_improvement = 0
    early_stopped = False
    in_memory_best: Optional[Dict[str, Tensor]] = None

    for epoch in range(1, max_epochs + 1):
        train_mse = train_one_epoch(
            model,
            train_loader,
            optimizer,
            criterion,
            target_scaler,
            device,
            grad_scaler,
        )
        _, train_true, train_pred = predict_loader(
            model, train_eval_loader, target_scaler, device
        )
        _, val_true, val_pred = predict_loader(model, val_loader, target_scaler, device)
        train_rmse = calculate_metrics(train_true, train_pred)["rmse"]
        val_rmse = calculate_metrics(val_true, val_pred)["rmse"]
        current_lr = float(optimizer.param_groups[0]["lr"])
        history.append(
            {
                "epoch": float(epoch),
                "train_mse_loss": train_mse,
                "train_rmse_eV": train_rmse,
                "val_rmse_eV": val_rmse,
                "learning_rate": current_lr,
            }
        )
        logger.info(
            "Epoch %04d | train MSE %.8f | train RMSE %.6f %s | val RMSE %.6f %s | lr %.6e",
            epoch,
            train_mse,
            train_rmse,
            TARGET_UNIT,
            val_rmse,
            TARGET_UNIT,
            current_lr,
        )
        scheduler.step(val_rmse)

        if val_rmse < best_val_rmse - MIN_DELTA:
            best_val_rmse = val_rmse
            best_epoch = epoch
            epochs_without_improvement = 0
            payload = checkpoint_payload(
                model,
                optimizer,
                scheduler,
                model_config,
                graph_config,
                target_scaler,
                epoch,
                best_val_rmse,
                hyperparameters,
            )
            if save_checkpoint:
                torch.save(payload, CHECKPOINT_PATH)
            else:
                in_memory_best = {
                    key: value.detach().cpu().clone()
                    for key, value in model.state_dict().items()
                }
        else:
            epochs_without_improvement += 1

        if epochs_without_improvement >= patience:
            early_stopped = True
            logger.info("Early stopping triggered at epoch %d.", epoch)
            break

    if save_checkpoint:
        checkpoint = torch.load(
            CHECKPOINT_PATH, map_location=device, weights_only=False
        )
        model.load_state_dict(checkpoint["model_state_dict"])
    elif in_memory_best is not None:
        model.load_state_dict(in_memory_best)
    return (
        model,
        target_scaler,
        history,
        best_epoch,
        best_val_rmse,
        early_stopped,
    )


def run_optuna(
    train_graphs: Sequence[CrystalGraph],
    val_graphs: Sequence[CrystalGraph],
    graph_config: GraphConfig,
    device: torch.device,
    logger: logging.Logger,
) -> Tuple[Dict[str, Any], ModelConfig]:
    try:
        import optuna
    except ImportError as exc:
        raise ImportError("Optuna is required when USE_OPTUNA=True.") from exc

    optuna.logging.set_verbosity(optuna.logging.WARNING)

    def objective(trial: Any) -> float:
        set_global_seed(SEED)
        model_config = ModelConfig(
            hidden_dim=trial.suggest_categorical("hidden_dim", [64, 96, 128, 192]),
            num_filters=trial.suggest_categorical("num_filters", [64, 96, 128, 192]),
            num_gaussians=trial.suggest_categorical("num_gaussians", [32, 48, 64, 96]),
            num_interactions=trial.suggest_int("num_interactions", 3, 7),
            cutoff=graph_config.cutoff,
            readout=READOUT,
            use_cosine_cutoff=USE_COSINE_CUTOFF,
        )
        hyperparameters = {
            "batch_size": trial.suggest_categorical("batch_size", [16, 24, 32, 48]),
            "learning_rate": trial.suggest_float(
                "learning_rate", 2.0e-4, 3.0e-3, log=True
            ),
            "weight_decay": trial.suggest_float(
                "weight_decay", 1.0e-8, 1.0e-4, log=True
            ),
        }
        _, _, _, _, best_rmse, _ = train_model(
            train_graphs=train_graphs,
            val_graphs=val_graphs,
            model_config=model_config,
            graph_config=graph_config,
            hyperparameters=hyperparameters,
            device=device,
            logger=logger,
            save_checkpoint=False,
            max_epochs=OPTUNA_MAX_EPOCHS,
            patience=OPTUNA_PATIENCE,
        )
        if device.type == "cuda":
            torch.cuda.empty_cache()
        return best_rmse

    study = optuna.create_study(direction="minimize")
    study.optimize(objective, n_trials=OPTUNA_N_TRIALS, timeout=OPTUNA_TIMEOUT)
    best = dict(study.best_params)
    best["best_val_rmse"] = float(study.best_value)
    with open(
        TABLE_DIR / f"{MODEL_NAME}_optuna_best_params.json", "w", encoding="utf-8"
    ) as handle:
        json.dump(best, handle, indent=2)
    logger.info("Optuna best validation RMSE: %.6f %s", study.best_value, TARGET_UNIT)
    logger.info("Optuna best parameters: %s", study.best_params)
    model_config = ModelConfig(
        hidden_dim=int(best["hidden_dim"]),
        num_filters=int(best["num_filters"]),
        num_gaussians=int(best["num_gaussians"]),
        num_interactions=int(best["num_interactions"]),
        cutoff=graph_config.cutoff,
        readout=READOUT,
        use_cosine_cutoff=USE_COSINE_CUTOFF,
    )
    hyperparameters = {
        "batch_size": int(best["batch_size"]),
        "learning_rate": float(best["learning_rate"]),
        "weight_decay": float(best["weight_decay"]),
    }
    return hyperparameters, model_config


def load_trained_model(
    checkpoint_path: str | Path,
    device: Optional[torch.device] = None,
) -> Tuple[SchNetBandgap, TargetScaler, GraphConfig, Dict[str, Any]]:
    target_device = device or torch.device(
        "cuda" if torch.cuda.is_available() else "cpu"
    )
    checkpoint = torch.load(
        checkpoint_path, map_location=target_device, weights_only=False
    )
    model_config = ModelConfig(**checkpoint["model_config"])
    graph_config = GraphConfig(**checkpoint["graph_config"])
    target_scaler = TargetScaler(**checkpoint["target_scaler"])
    model = SchNetBandgap(model_config).to(target_device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    return model, target_scaler, graph_config, checkpoint


@torch.inference_mode()
def predict_cifs(
    model: SchNetBandgap,
    cif_paths: Sequence[str | Path],
    target_scaler: TargetScaler,
    graph_config: GraphConfig,
    device: Optional[torch.device] = None,
    batch_size: int = BATCH_SIZE,
) -> pd.DataFrame:
    target_device = device or next(model.parameters()).device
    graphs: List[CrystalGraph] = []
    errors: Dict[str, str] = {}
    for raw_path in cif_paths:
        path = Path(raw_path)
        try:
            structure = Structure.from_file(path)
            graphs.append(
                structure_to_graph(structure, path.name, float("nan"), graph_config)
            )
        except Exception as exc:
            errors[str(path)] = str(exc)
    rows: List[Dict[str, Any]] = []
    if graphs:
        loader = make_loader(graphs, batch_size, False, target_device)
        model.eval()
        for batch in loader:
            batch = batch.to(target_device)
            with autocast_context(target_device):
                scaled = model(batch).view(-1)
            predictions = target_scaler.inverse_numpy(
                scaled.float().cpu().numpy().astype(np.float64)
            )
            rows.extend(
                {
                    "cif_name": name,
                    f"predicted_{TARGET_NAME}_{TARGET_UNIT}": float(prediction),
                    "status": "ok",
                }
                for name, prediction in zip(batch.cif_names, predictions)
            )
    rows.extend(
        {
            "cif_name": path,
            f"predicted_{TARGET_NAME}_{TARGET_UNIT}": float("nan"),
            "status": error,
        }
        for path, error in errors.items()
    )
    return pd.DataFrame(rows)


def save_prediction_dat(
    split_name: str,
    names: Sequence[str],
    targets: np.ndarray,
    predictions: np.ndarray,
) -> pd.DataFrame:
    frame = pd.DataFrame(
        {
            "cif_name": names,
            "true_bandgap_eV": targets,
            "predicted_bandgap_eV": predictions,
            "error_eV": predictions - targets,
            "absolute_error_eV": np.abs(predictions - targets),
        }
    )
    frame.to_csv(
        DAT_DIR / f"{MODEL_NAME}_parity_{split_name}.dat",
        sep="\t",
        index=False,
        float_format="%.10f",
    )
    return frame


def metric_text(metrics: Dict[str, float]) -> str:
    return (
        f"MAE = {metrics['mae']:.4f} {TARGET_UNIT}\n"
        f"RMSE = {metrics['rmse']:.4f} {TARGET_UNIT}\n"
        f"R² = {metrics['r2']:.4f}"
    )


def parity_limits(
    all_targets: np.ndarray, all_predictions: np.ndarray
) -> Tuple[float, float]:
    lower = float(min(all_targets.min(), all_predictions.min()))
    upper = float(max(all_targets.max(), all_predictions.max()))
    span = max(upper - lower, 1.0e-6)
    margin = 0.05 * span
    return lower - margin, upper + margin


def style_axes(ax: plt.Axes) -> None:
    ax.tick_params(labelsize=TICK_FONTSIZE)
    if USE_GRID:
        ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)


def plot_parity(
    split_name: str,
    targets: np.ndarray,
    predictions: np.ndarray,
    metrics: Dict[str, float],
    color: str,
    limits: Tuple[float, float],
) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    ax.scatter(targets, predictions, s=20, alpha=0.78, color=color, edgecolors="none")
    ax.plot(limits, limits, color="black", linestyle="--", linewidth=1.2, label="y = x")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(
        f"{MODEL_NAME}: {split_name.capitalize()} Set", fontsize=TITLE_FONTSIZE
    )
    ax.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.text(
        0.04,
        0.96,
        metric_text(metrics),
        transform=ax.transAxes,
        va="top",
        fontsize=ANNOTATION_FONTSIZE,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85},
    )
    ax.legend(fontsize=LEGEND_FONTSIZE)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_parity_{split_name}.jpg",
        dpi=FIG_DPI,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_combined_parity(
    results: Dict[str, Dict[str, Any]], limits: Tuple[float, float]
) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 6.2))
    for index, split_name in enumerate(("train", "val", "test")):
        result = results[split_name]
        ax.scatter(
            result["targets"],
            result["predictions"],
            s=20,
            alpha=0.72,
            color=COLORS[index],
            edgecolors="none",
            label=f"{split_name.capitalize()} (n={len(result['targets'])})",
        )
    ax.plot(limits, limits, color="black", linestyle="--", linewidth=1.2, label="y = x")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"{MODEL_NAME}: All Splits", fontsize=TITLE_FONTSIZE)
    ax.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_parity_all.jpg",
        dpi=FIG_DPI,
        bbox_inches="tight",
    )
    plt.close(fig)


def save_all_prediction_dat(results: Dict[str, Dict[str, Any]]) -> None:
    frames = []
    for split_name in ("train", "val", "test"):
        frame = results[split_name]["frame"].copy()
        frame.insert(0, "split", split_name)
        frames.append(frame)
    pd.concat(frames, ignore_index=True).to_csv(
        DAT_DIR / f"{MODEL_NAME}_parity_all.dat",
        sep="\t",
        index=False,
        float_format="%.10f",
    )


def save_and_plot_history(history: Sequence[Dict[str, float]], best_epoch: int) -> None:
    frame = pd.DataFrame(history)
    frame["epoch"] = frame["epoch"].astype(int)
    frame.to_csv(
        DAT_DIR / f"{MODEL_NAME}_rmse_curve.dat",
        sep="\t",
        index=False,
        float_format="%.10f",
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
    ax.set_title(f"{MODEL_NAME}: RMSE Learning Curve", fontsize=TITLE_FONTSIZE)
    ax.set_xlabel("Epoch", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"RMSE ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_rmse_curve.jpg",
        dpi=FIG_DPI,
        bbox_inches="tight",
    )
    plt.close(fig)


def save_metrics_table(results: Dict[str, Dict[str, Any]]) -> None:
    rows = []
    for split_name in ("train", "val", "test"):
        metrics = results[split_name]["metrics"]
        rows.append(
            {
                "split": split_name,
                "n_samples": len(results[split_name]["targets"]),
                "mae_eV": metrics["mae"],
                "rmse_eV": metrics["rmse"],
                "r2": metrics["r2"],
            }
        )
    pd.DataFrame(rows).to_csv(
        TABLE_DIR / f"{MODEL_NAME}_metrics.dat",
        sep="\t",
        index=False,
        float_format="%.6f",
    )


def log_configuration(
    logger: logging.Logger,
    model: nn.Module,
    model_config: ModelConfig,
    graph_config: GraphConfig,
    hyperparameters: Dict[str, Any],
    subset_sizes: Dict[str, int],
) -> None:
    total, trainable = count_parameters(model)
    logger.info("Model name: %s", MODEL_NAME)
    logger.info("Run version: %s", RUN_VERSION)
    logger.info("Output root: %s", OUTPUT_ROOT)
    logger.info("Excel path: %s", EXCEL_PATH)
    logger.info("CIF directory: %s", CIF_DIR)
    logger.info("Excel columns: first column=CIF name, second column=target")
    logger.info("Target: %s (%s)", TARGET_NAME, TARGET_UNIT)
    logger.info("Seed: %d", SEED)
    logger.info("Split ratios: %.2f/%.2f/%.2f", TRAIN_RATIO, VAL_RATIO, TEST_RATIO)
    logger.info("Subset sizes: %s", subset_sizes)
    logger.info("Graph configuration: %s", asdict(graph_config))
    logger.info("Model configuration: %s", asdict(model_config))
    logger.info("Training hyperparameters: %s", hyperparameters)
    logger.info("Loss function: MSELoss")
    logger.info("Optimizer: AdamW")
    logger.info("Scheduler: ReduceLROnPlateau")
    logger.info("AMP enabled: %s; dtype: %s", USE_AMP, AMP_DTYPE)
    logger.info("Total parameters: %d", total)
    logger.info("Trainable parameters: %d", trainable)
    logger.info("Model structure:\n%s", model)


def main() -> None:
    total_start = time.perf_counter()
    create_output_directories()
    logger = setup_logger()
    set_global_seed(SEED)
    device = get_device(logger)
    logger.info("SchNet uses pure PyTorch tensor operations on the selected device.")

    data_start = time.perf_counter()
    records = load_excel_records(logger)
    graph_config = GraphConfig()
    graphs = load_or_build_graphs(records, graph_config, logger)
    split_map = create_or_load_split(graphs, logger)
    train_graphs, val_graphs, test_graphs = partition_graphs(graphs, split_map)
    data_elapsed = time.perf_counter() - data_start
    logger.info(
        "Data preparation time: %.3f s (%s)",
        data_elapsed,
        format_duration(data_elapsed),
    )

    model_config = ModelConfig()
    hyperparameters: Dict[str, Any] = {
        "batch_size": BATCH_SIZE,
        "learning_rate": LEARNING_RATE,
        "weight_decay": WEIGHT_DECAY,
    }
    logger.info("Optuna enabled: %s", USE_OPTUNA)
    if USE_OPTUNA:
        hyperparameters, model_config = run_optuna(
            train_graphs, val_graphs, graph_config, device, logger
        )

    preview_model = SchNetBandgap(model_config).to(device)
    log_configuration(
        logger,
        preview_model,
        model_config,
        graph_config,
        hyperparameters,
        {
            "train": len(train_graphs),
            "val": len(val_graphs),
            "test": len(test_graphs),
        },
    )
    del preview_model

    train_start = time.perf_counter()
    model, target_scaler, history, best_epoch, best_val_rmse, early_stopped = (
        train_model(
            train_graphs=train_graphs,
            val_graphs=val_graphs,
            model_config=model_config,
            graph_config=graph_config,
            hyperparameters=hyperparameters,
            device=device,
            logger=logger,
            save_checkpoint=True,
            max_epochs=MAX_EPOCHS,
            patience=PATIENCE,
        )
    )
    train_elapsed = time.perf_counter() - train_start
    logger.info("Best epoch: %d", best_epoch)
    logger.info("Best validation RMSE: %.6f %s", best_val_rmse, TARGET_UNIT)
    logger.info("Early stopping triggered: %s", early_stopped)
    logger.info(
        "Training time: %.3f s (%s)", train_elapsed, format_duration(train_elapsed)
    )
    logger.info("Best checkpoint: %s", CHECKPOINT_PATH)

    evaluation_start = time.perf_counter()
    model, target_scaler, graph_config, checkpoint = load_trained_model(
        CHECKPOINT_PATH, device
    )
    eval_batch_size = int(checkpoint["hyperparameters"]["batch_size"])
    subsets = {"train": train_graphs, "val": val_graphs, "test": test_graphs}
    results: Dict[str, Dict[str, Any]] = {}
    for split_name, subset in subsets.items():
        loader = make_loader(subset, eval_batch_size, False, device)
        names, targets, predictions = predict_loader(
            model, loader, target_scaler, device
        )
        metrics = calculate_metrics(targets, predictions)
        frame = save_prediction_dat(split_name, names, targets, predictions)
        results[split_name] = {
            "names": names,
            "targets": targets,
            "predictions": predictions,
            "metrics": metrics,
            "frame": frame,
        }
        logger.info(
            "%s metrics | MAE %.6f %s | RMSE %.6f %s | R2 %.6f",
            split_name.capitalize(),
            metrics["mae"],
            TARGET_UNIT,
            metrics["rmse"],
            TARGET_UNIT,
            metrics["r2"],
        )

    all_targets = np.concatenate(
        [results[key]["targets"] for key in ("train", "val", "test")]
    )
    all_predictions = np.concatenate(
        [results[key]["predictions"] for key in ("train", "val", "test")]
    )
    limits = parity_limits(all_targets, all_predictions)
    for index, split_name in enumerate(("train", "val", "test")):
        result = results[split_name]
        plot_parity(
            split_name,
            result["targets"],
            result["predictions"],
            result["metrics"],
            COLORS[index],
            limits,
        )
    plot_combined_parity(results, limits)
    save_all_prediction_dat(results)
    save_and_plot_history(history, best_epoch)
    save_metrics_table(results)
    evaluation_elapsed = time.perf_counter() - evaluation_start
    logger.info(
        "Evaluation and plotting time: %.3f s (%s)",
        evaluation_elapsed,
        format_duration(evaluation_elapsed),
    )
    total_elapsed = time.perf_counter() - total_start
    logger.info(
        "Total runtime: %.3f s (%s)", total_elapsed, format_duration(total_elapsed)
    )


if __name__ == "__main__":
    main()

from __future__ import annotations

import copy
import hashlib
import json
import logging
import math
import os

os.environ.setdefault("CUBLAS_WORKSPACE_CONFIG", ":4096:8")
import random
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple, Union

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
from torch.utils.data import Dataset
from torch_geometric.data import Batch, Data
from torch_geometric.loader import DataLoader
from torch_geometric.nn import Set2Set
from torch_geometric.utils import scatter
from tqdm.auto import tqdm


# ============================== User configuration ==============================

MODEL_NAME = "MEGNet"
RUN_VERSION = "v2"

EXCEL_PATH = "./data.xlsx"
CIF_DIR = "./cif"

TARGET_NAME = "bandgap"
TARGET_UNIT = "eV"

SEED = 42
TRAIN_RATIO = 0.80
VAL_RATIO = 0.10
TEST_RATIO = 0.10

GRAPH_CUTOFF = 4.0
GAUSSIAN_MIN = 0.0
GAUSSIAN_MAX = 5.0
GAUSSIAN_BASIS_SIZE = 100
GAUSSIAN_WIDTH = 0.5
MAX_ATOMIC_NUMBER = 118

EMBEDDING_DIM = 16
N1 = 64
N2 = 32
N3 = 16
NUM_MEGNET_BLOCKS = 3
SET2SET_STEPS = 3
DROPOUT = 0.0

BATCH_SIZE = 64
MAX_EPOCHS = 300
LEARNING_RATE = 1.0e-3
WEIGHT_DECAY = 0.0
PATIENCE = 80
MIN_DELTA = 1.0e-5
GRAD_CLIP_NORM = 3.0
LR_FACTOR = 0.5
LR_PATIENCE = 20
MIN_LR = 1.0e-6
NUM_WORKERS = 0
PIN_MEMORY = True
USE_AMP = True

USE_OPTUNA = False
OPTUNA_N_TRIALS = 30
OPTUNA_TIMEOUT = None
OPTUNA_MAX_EPOCHS = 250
OPTUNA_PATIENCE = 40

FIG_DPI = 600
TITLE_FONTSIZE = 15
LABEL_FONTSIZE = 13
TICK_FONTSIZE = 11
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


def softplus2(x: torch.Tensor) -> torch.Tensor:
    return F.softplus(x) - math.log(2.0)


def format_duration(seconds: float) -> str:
    seconds_int = max(0, int(round(seconds)))
    hours, remainder = divmod(seconds_int, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


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
    if hasattr(torch, "use_deterministic_algorithms"):
        torch.use_deterministic_algorithms(True, warn_only=True)
    if hasattr(torch, "set_float32_matmul_precision"):
        torch.set_float32_matmul_precision("high")


def create_output_directories() -> None:
    for path in [
        OUTPUT_ROOT,
        FIGURE_DIR,
        DAT_DIR,
        TABLE_DIR,
        LOG_DIR,
        SPLIT_DIR,
        CACHE_DIR,
    ]:
        path.mkdir(parents=True, exist_ok=True)


def setup_logger(log_path: Path) -> logging.Logger:
    logger = logging.getLogger(f"{MODEL_NAME}_{RUN_VERSION}")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    for handler in list(logger.handlers):
        logger.removeHandler(handler)
        handler.close()
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    file_handler = logging.FileHandler(log_path, mode="w", encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def get_device(logger: Optional[logging.Logger] = None) -> torch.device:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if logger is not None:
        logger.info("Device: %s", device)
        logger.info("PyTorch version: %s", torch.__version__)
        if device.type == "cuda":
            properties = torch.cuda.get_device_properties(device)
            logger.info("CUDA version: %s", torch.version.cuda)
            logger.info("GPU: %s", properties.name)
            logger.info("GPU memory: %.2f GiB", properties.total_memory / (1024**3))
        else:
            logger.info("CUDA is unavailable; CPU execution will be used.")
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
        if re.fullmatch(r"[+-]?\d+\.0+", name):
            name = str(int(float(name)))
    if not name:
        raise ValueError("Empty CIF filename")
    if not name.lower().endswith(".cif"):
        name += ".cif"
    return name


@dataclass(frozen=True)
class SampleRecord:
    cif_name: str
    cif_path: str
    target: float


def load_excel_records(
    excel_path: Union[str, Path], cif_dir: Union[str, Path], logger: logging.Logger
) -> Tuple[List[SampleRecord], List[Tuple[str, str]]]:
    excel_path = Path(excel_path)
    cif_dir = Path(cif_dir)
    if not excel_path.is_file():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    if not cif_dir.is_dir():
        raise FileNotFoundError(f"CIF directory not found: {cif_dir}")

    frame = pd.read_excel(excel_path)
    if frame.shape[1] < 2:
        raise ValueError("The Excel file must contain at least two columns.")

    records: List[SampleRecord] = []
    failures: List[Tuple[str, str]] = []
    seen: set[str] = set()
    for row_index, row in frame.iloc[:, :2].iterrows():
        raw_name = row.iloc[0]
        raw_target = row.iloc[1]
        display_name = str(raw_name)
        try:
            cif_name = normalize_cif_name(raw_name)
            if cif_name in seen:
                raise ValueError("Duplicate CIF filename")
            target = float(raw_target)
            if not np.isfinite(target):
                raise ValueError("Target is NaN or infinite")
            cif_path = cif_dir / cif_name
            if not cif_path.is_file():
                raise FileNotFoundError(f"CIF file not found: {cif_path}")
            records.append(SampleRecord(cif_name, str(cif_path), target))
            seen.add(cif_name)
        except Exception as exc:
            reason = f"Excel row {row_index + 2}: {exc}"
            failures.append((display_name, reason))
            logger.warning("Skipped sample %s | %s", display_name, reason)
    if len(records) < 10:
        raise ValueError(
            f"Only {len(records)} valid rows remain; at least 10 are required."
        )
    logger.info("Excel rows: %d", len(frame))
    logger.info("Valid rows before graph validation: %d", len(records))
    logger.info("Invalid Excel rows: %d", len(failures))
    return records, failures


def graph_config() -> Dict[str, Any]:
    return {
        "cutoff": GRAPH_CUTOFF,
        "gaussian_min": GAUSSIAN_MIN,
        "gaussian_max": GAUSSIAN_MAX,
        "gaussian_basis_size": GAUSSIAN_BASIS_SIZE,
        "gaussian_width": GAUSSIAN_WIDTH,
        "max_atomic_number": MAX_ATOMIC_NUMBER,
        "state_dimension": 2,
    }


def model_config() -> Dict[str, Any]:
    return {
        "max_atomic_number": MAX_ATOMIC_NUMBER,
        "gaussian_min": GAUSSIAN_MIN,
        "gaussian_max": GAUSSIAN_MAX,
        "gaussian_basis_size": GAUSSIAN_BASIS_SIZE,
        "gaussian_width": GAUSSIAN_WIDTH,
        "embedding_dim": EMBEDDING_DIM,
        "n1": N1,
        "n2": N2,
        "n3": N3,
        "num_blocks": NUM_MEGNET_BLOCKS,
        "set2set_steps": SET2SET_STEPS,
        "dropout": DROPOUT,
        "state_dim": 2,
    }


def cache_key(cif_path: Union[str, Path], config: Dict[str, Any]) -> str:
    path = Path(cif_path).resolve()
    stat = path.stat()
    payload = {
        "path": str(path),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "graph_config": config,
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode("utf-8")
    ).hexdigest()


def build_crystal_graph(
    cif_path: Union[str, Path],
    target: Optional[float] = None,
    cif_name: Optional[str] = None,
    config: Optional[Dict[str, Any]] = None,
) -> Data:
    config = graph_config() if config is None else config
    structure = Structure.from_file(str(cif_path))
    if len(structure) == 0:
        raise ValueError("Structure contains no atomic sites")
    if not structure.is_ordered:
        raise ValueError(
            "Disordered structures are not supported by integer atom embeddings"
        )

    atomic_numbers = []
    for site in structure:
        z = int(site.specie.Z)
        if z < 1 or z > int(config["max_atomic_number"]):
            raise ValueError(f"Atomic number {z} is outside the configured range")
        atomic_numbers.append(z)

    center, neighbor, offsets, distances = structure.get_neighbor_list(
        r=float(config["cutoff"]), exclude_self=True
    )
    if len(center) == 0:
        raise ValueError("No periodic neighbors were found within the cutoff")

    edge_index = torch.tensor(np.vstack([center, neighbor]), dtype=torch.long)
    data = Data(
        z=torch.tensor(atomic_numbers, dtype=torch.long),
        edge_index=edge_index,
        edge_distance=torch.tensor(distances, dtype=torch.float32),
        edge_shift=torch.tensor(offsets, dtype=torch.float32),
        u=torch.zeros((1, int(config["state_dimension"])), dtype=torch.float32),
    )
    data.num_nodes = len(atomic_numbers)
    data.cif_name = cif_name if cif_name is not None else Path(cif_path).name
    if target is not None:
        data.y = torch.tensor([float(target)], dtype=torch.float32)
    return data


def load_or_build_graph(record: SampleRecord, config: Dict[str, Any]) -> Data:
    key = cache_key(record.cif_path, config)
    cache_path = CACHE_DIR / f"{key}.pt"
    if cache_path.is_file():
        try:
            data = torch.load(cache_path, map_location="cpu", weights_only=False)
            if isinstance(data, Data):
                data.y = torch.tensor([record.target], dtype=torch.float32)
                data.cif_name = record.cif_name
                return data
        except Exception:
            pass
    data = build_crystal_graph(record.cif_path, record.target, record.cif_name, config)
    torch.save(data, cache_path)
    return data


def validate_and_cache_records(
    records: Sequence[SampleRecord], config: Dict[str, Any], logger: logging.Logger
) -> Tuple[List[SampleRecord], Dict[str, Data], List[Tuple[str, str]]]:
    valid_records: List[SampleRecord] = []
    graph_map: Dict[str, Data] = {}
    failures: List[Tuple[str, str]] = []
    for record in tqdm(records, desc="Validating CIF graphs"):
        try:
            data = load_or_build_graph(record, config)
            valid_records.append(record)
            graph_map[record.cif_name] = data
        except Exception as exc:
            reason = f"{type(exc).__name__}: {exc}"
            failures.append((record.cif_name, reason))
            logger.warning("Skipped sample %s | %s", record.cif_name, reason)
    if len(valid_records) < 10:
        raise ValueError(
            f"Only {len(valid_records)} valid graphs remain; at least 10 are required."
        )
    logger.info("Valid crystal graphs: %d", len(valid_records))
    logger.info("Invalid crystal graphs: %d", len(failures))
    return valid_records, graph_map, failures


def _split_counts_are_valid(frame: pd.DataFrame) -> bool:
    counts = frame["split"].value_counts().to_dict()
    return all(counts.get(name, 0) > 0 for name in ["train", "val", "test"])


def create_or_load_split(
    records: Sequence[SampleRecord], split_path: Path, logger: logging.Logger
) -> Dict[str, str]:
    expected = {record.cif_name: record.target for record in records}
    if split_path.is_file():
        try:
            frame = pd.read_csv(split_path)
            required = {"cif_name", "target", "split"}
            if not required.issubset(frame.columns):
                raise ValueError("Required split columns are missing")
            if frame["cif_name"].duplicated().any():
                raise ValueError("Duplicate CIF names exist in the split file")
            split_names = set(frame["cif_name"].astype(str))
            if split_names != set(expected):
                raise ValueError(
                    "Split file samples do not match the current valid dataset"
                )
            target_map = dict(
                zip(frame["cif_name"].astype(str), frame["target"].astype(float))
            )
            for name, target in expected.items():
                if not np.isclose(target_map[name], target, rtol=0.0, atol=1.0e-10):
                    raise ValueError(f"Target mismatch for {name}")
            if not set(frame["split"]).issubset(
                {"train", "val", "test"}
            ) or not _split_counts_are_valid(frame):
                raise ValueError("Invalid split labels or empty split")
            logger.info("Reused split file: %s", split_path)
            return dict(zip(frame["cif_name"].astype(str), frame["split"].astype(str)))
        except Exception as exc:
            logger.warning("Existing split file was not reused: %s", exc)

    if not math.isclose(TRAIN_RATIO + VAL_RATIO + TEST_RATIO, 1.0, abs_tol=1.0e-10):
        raise ValueError("TRAIN_RATIO + VAL_RATIO + TEST_RATIO must equal 1.0")
    names = [record.cif_name for record in records]
    train_names, temp_names = train_test_split(
        names, test_size=VAL_RATIO + TEST_RATIO, random_state=SEED, shuffle=True
    )
    relative_test = TEST_RATIO / (VAL_RATIO + TEST_RATIO)
    val_names, test_names = train_test_split(
        temp_names, test_size=relative_test, random_state=SEED, shuffle=True
    )
    split_map = {name: "train" for name in train_names}
    split_map.update({name: "val" for name in val_names})
    split_map.update({name: "test" for name in test_names})
    frame = pd.DataFrame(
        [
            {"cif_name": r.cif_name, "target": r.target, "split": split_map[r.cif_name]}
            for r in records
        ]
    )
    frame.to_csv(split_path, index=False)
    logger.info("Created split file: %s", split_path)
    return split_map


class CrystalGraphDataset(Dataset):
    def __init__(
        self,
        records: Sequence[SampleRecord],
        graph_map: Dict[str, Data],
        target_mean: float,
        target_std: float,
    ) -> None:
        self.records = list(records)
        self.graph_map = graph_map
        self.target_mean = float(target_mean)
        self.target_std = float(target_std)

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, index: int) -> Data:
        record = self.records[index]
        data = self.graph_map[record.cif_name].clone()
        data.y = torch.tensor(
            [(record.target - self.target_mean) / self.target_std], dtype=torch.float32
        )
        data.y_raw = torch.tensor([record.target], dtype=torch.float32)
        data.cif_name = record.cif_name
        return data


def make_loaders(
    split_records: Dict[str, List[SampleRecord]],
    graph_map: Dict[str, Data],
    target_mean: float,
    target_std: float,
    batch_size: int,
) -> Dict[str, DataLoader]:
    loaders: Dict[str, DataLoader] = {}
    generator = torch.Generator()
    generator.manual_seed(SEED)
    for split_name in ["train", "val", "test"]:
        dataset = CrystalGraphDataset(
            split_records[split_name], graph_map, target_mean, target_std
        )
        loaders[split_name] = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=split_name == "train",
            num_workers=NUM_WORKERS,
            pin_memory=PIN_MEMORY and torch.cuda.is_available(),
            persistent_workers=NUM_WORKERS > 0,
            generator=generator if split_name == "train" else None,
        )
    return loaders


class MLP(nn.Module):
    def __init__(self, dimensions: Sequence[int]) -> None:
        super().__init__()
        layers: List[nn.Module] = []
        for in_dim, out_dim in zip(dimensions[:-1], dimensions[1:]):
            layers.append(nn.Linear(in_dim, out_dim))
            layers.append(Softplus2())
        self.network = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.network(x)


class Softplus2(nn.Module):
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return softplus2(x)


class GaussianExpansion(nn.Module):
    def __init__(
        self, minimum: float, maximum: float, basis_size: int, width: float
    ) -> None:
        super().__init__()
        centers = torch.linspace(minimum, maximum, basis_size)
        self.register_buffer("centers", centers)
        self.width = float(width)

    def forward(self, distances: torch.Tensor) -> torch.Tensor:
        return torch.exp(
            -((distances.view(-1, 1) - self.centers.view(1, -1)) ** 2) / (self.width**2)
        )


class MEGNetBlock(nn.Module):
    def __init__(self, feature_dim: int, n1: int, dropout: float) -> None:
        super().__init__()
        self.edge_mlp = MLP([4 * feature_dim, n1, n1, feature_dim])
        self.node_mlp = MLP([3 * feature_dim, n1, n1, feature_dim])
        self.global_mlp = MLP([3 * feature_dim, n1, n1, feature_dim])
        self.output_dropout = nn.Dropout(dropout) if dropout > 0.0 else nn.Identity()

    def forward(
        self,
        node: torch.Tensor,
        edge: torch.Tensor,
        state: torch.Tensor,
        edge_index: torch.Tensor,
        node_batch: torch.Tensor,
        edge_batch: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor]:
        receiver, sender = edge_index[0], edge_index[1]
        edge_input = torch.cat(
            [node[receiver], node[sender], edge, state[edge_batch]], dim=-1
        )
        edge_updated = self.edge_mlp(edge_input)

        edge_to_node = scatter(
            edge_updated, receiver, dim=0, dim_size=node.size(0), reduce="mean"
        )
        node_input = torch.cat([edge_to_node, node, state[node_batch]], dim=-1)
        node_updated = self.node_mlp(node_input)

        edge_to_global = scatter(
            edge_updated, edge_batch, dim=0, dim_size=state.size(0), reduce="mean"
        )
        node_to_global = scatter(
            node_updated, node_batch, dim=0, dim_size=state.size(0), reduce="mean"
        )
        global_input = torch.cat([edge_to_global, node_to_global, state], dim=-1)
        state_updated = self.global_mlp(global_input)
        return (
            self.output_dropout(node_updated),
            self.output_dropout(edge_updated),
            self.output_dropout(state_updated),
        )


class MEGNetRegressor(nn.Module):
    def __init__(
        self,
        max_atomic_number: int = 118,
        gaussian_min: float = 0.0,
        gaussian_max: float = 5.0,
        gaussian_basis_size: int = 100,
        gaussian_width: float = 0.5,
        embedding_dim: int = 16,
        n1: int = 64,
        n2: int = 32,
        n3: int = 16,
        num_blocks: int = 3,
        set2set_steps: int = 3,
        dropout: float = 0.0,
        state_dim: int = 2,
    ) -> None:
        super().__init__()
        self.atom_embedding = nn.Embedding(
            max_atomic_number + 1, embedding_dim, padding_idx=0
        )
        self.distance_expansion = GaussianExpansion(
            gaussian_min, gaussian_max, gaussian_basis_size, gaussian_width
        )
        self.pre_node = MLP([embedding_dim, n1, n2])
        self.pre_edge = MLP([gaussian_basis_size, n1, n2])
        self.pre_state = MLP([state_dim, n1, n2])
        self.block_pre_node = nn.ModuleList(
            [MLP([n2, n1, n2]) for _ in range(num_blocks - 1)]
        )
        self.block_pre_edge = nn.ModuleList(
            [MLP([n2, n1, n2]) for _ in range(num_blocks - 1)]
        )
        self.block_pre_state = nn.ModuleList(
            [MLP([n2, n1, n2]) for _ in range(num_blocks - 1)]
        )
        self.blocks = nn.ModuleList(
            [MEGNetBlock(n2, n1, dropout) for _ in range(num_blocks)]
        )
        self.node_projection = nn.Linear(n2, n3)
        self.edge_projection = nn.Linear(n2, n3)
        self.node_set2set = Set2Set(n3, processing_steps=set2set_steps)
        self.edge_set2set = Set2Set(n3, processing_steps=set2set_steps)
        self.final_dropout = nn.Dropout(dropout) if dropout > 0.0 else nn.Identity()
        self.readout_0 = nn.Linear(4 * n3 + n2, n2)
        self.readout_1 = nn.Linear(n2, n3)
        self.readout_2 = nn.Linear(n3, 1)

    def forward(self, batch: Batch) -> torch.Tensor:
        node = self.pre_node(self.atom_embedding(batch.z))
        edge = self.pre_edge(self.distance_expansion(batch.edge_distance))
        state = self.pre_state(batch.u.view(batch.num_graphs, -1))
        edge_batch = batch.batch[batch.edge_index[0]]

        for index, block in enumerate(self.blocks):
            block_node, block_edge, block_state = node, edge, state
            if index > 0:
                block_node = self.block_pre_node[index - 1](block_node)
                block_edge = self.block_pre_edge[index - 1](block_edge)
                block_state = self.block_pre_state[index - 1](block_state)
            delta_node, delta_edge, delta_state = block(
                block_node,
                block_edge,
                block_state,
                batch.edge_index,
                batch.batch,
                edge_batch,
            )
            node = node + delta_node
            edge = edge + delta_edge
            state = state + delta_state

        node_vector = self.node_set2set(self.node_projection(node), batch.batch)
        edge_vector = self.edge_set2set(self.edge_projection(edge), edge_batch)
        graph_vector = torch.cat([node_vector, edge_vector, state], dim=-1)
        graph_vector = self.final_dropout(graph_vector)
        graph_vector = softplus2(self.readout_0(graph_vector))
        graph_vector = softplus2(self.readout_1(graph_vector))
        return self.readout_2(graph_vector).view(-1)


@dataclass
class TrainingConfig:
    batch_size: int
    learning_rate: float
    weight_decay: float
    max_epochs: int
    patience: int
    min_delta: float
    grad_clip_norm: float
    lr_factor: float
    lr_patience: int
    min_lr: float


@dataclass
class FitResult:
    best_state: Dict[str, torch.Tensor]
    best_epoch: int
    best_val_rmse: float
    history: List[Dict[str, float]]
    early_stopped: bool
    elapsed_seconds: float


def autocast_context(device: torch.device):
    enabled = bool(USE_AMP and device.type == "cuda" and torch.cuda.is_bf16_supported())
    return torch.autocast(
        device_type=device.type, dtype=torch.bfloat16, enabled=enabled
    )


def normalized_to_raw(
    values: np.ndarray, target_mean: float, target_std: float
) -> np.ndarray:
    return values * target_std + target_mean


def train_one_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    criterion: nn.Module,
    device: torch.device,
    target_mean: float,
    target_std: float,
    grad_clip_norm: float,
) -> Tuple[float, float]:
    model.train()
    squared_error_sum = 0.0
    sample_count = 0
    loss_sum = 0.0
    for batch in loader:
        batch = batch.to(device, non_blocking=True)
        target = batch.y.view(-1)
        optimizer.zero_grad(set_to_none=True)
        with autocast_context(device):
            prediction = model(batch).view(-1)
            loss = criterion(prediction.float(), target.float())
        if not torch.isfinite(loss):
            raise FloatingPointError("Non-finite training loss encountered")
        loss.backward()
        if grad_clip_norm > 0.0:
            torch.nn.utils.clip_grad_norm_(model.parameters(), grad_clip_norm)
        optimizer.step()

        raw_prediction = (
            prediction.detach().float().cpu().numpy() * target_std + target_mean
        )
        raw_target = batch.y_raw.view(-1).detach().float().cpu().numpy()
        squared_error_sum += float(np.square(raw_prediction - raw_target).sum())
        loss_sum += float(loss.item()) * raw_target.size
        sample_count += raw_target.size
    return loss_sum / sample_count, math.sqrt(squared_error_sum / sample_count)


@torch.inference_mode()
def evaluate_loader(
    model: nn.Module,
    loader: DataLoader,
    criterion: nn.Module,
    device: torch.device,
    target_mean: float,
    target_std: float,
) -> Tuple[float, np.ndarray, np.ndarray, List[str]]:
    model.eval()
    predictions: List[np.ndarray] = []
    targets: List[np.ndarray] = []
    names: List[str] = []
    total_loss = 0.0
    sample_count = 0
    for batch in loader:
        batch = batch.to(device, non_blocking=True)
        with autocast_context(device):
            prediction = model(batch).view(-1)
            loss = criterion(prediction.float(), batch.y.view(-1).float())
        raw_prediction = prediction.float().cpu().numpy() * target_std + target_mean
        raw_target = batch.y_raw.view(-1).float().cpu().numpy()
        predictions.append(raw_prediction)
        targets.append(raw_target)
        batch_names = batch.cif_name
        names.extend(
            list(batch_names)
            if isinstance(batch_names, (list, tuple))
            else [str(batch_names)]
        )
        total_loss += float(loss.item()) * raw_target.size
        sample_count += raw_target.size
    return (
        total_loss / sample_count,
        np.concatenate(targets),
        np.concatenate(predictions),
        names,
    )


def fit_model(
    model: nn.Module,
    loaders: Dict[str, DataLoader],
    config: TrainingConfig,
    device: torch.device,
    target_mean: float,
    target_std: float,
    logger: Optional[logging.Logger] = None,
    verbose: bool = True,
) -> FitResult:
    criterion = nn.MSELoss()
    optimizer = torch.optim.Adam(
        model.parameters(), lr=config.learning_rate, weight_decay=config.weight_decay
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="min",
        factor=config.lr_factor,
        patience=config.lr_patience,
        min_lr=config.min_lr,
    )
    best_val_rmse = float("inf")
    best_epoch = 0
    best_state: Dict[str, torch.Tensor] = {}
    history: List[Dict[str, float]] = []
    stale_epochs = 0
    early_stopped = False
    start_time = time.perf_counter()

    epoch_iterator: Iterable[int] = range(1, config.max_epochs + 1)
    if not verbose:
        epoch_iterator = tqdm(epoch_iterator, desc="Optuna epochs", leave=False)
    for epoch in epoch_iterator:
        train_mse, train_rmse = train_one_epoch(
            model,
            loaders["train"],
            optimizer,
            criterion,
            device,
            target_mean,
            target_std,
            config.grad_clip_norm,
        )
        val_mse, val_true, val_pred, _ = evaluate_loader(
            model, loaders["val"], criterion, device, target_mean, target_std
        )
        val_rmse = float(np.sqrt(mean_squared_error(val_true, val_pred)))
        learning_rate = float(optimizer.param_groups[0]["lr"])
        history.append(
            {
                "epoch": float(epoch),
                "train_mse_loss": train_mse,
                "val_mse_loss": val_mse,
                "train_rmse_eV": train_rmse,
                "val_rmse_eV": val_rmse,
                "learning_rate": learning_rate,
            }
        )
        if logger is not None and verbose:
            logger.info(
                "Epoch %04d | train RMSE %.6f eV | val RMSE %.6f eV | lr %.6e",
                epoch,
                train_rmse,
                val_rmse,
                learning_rate,
            )

        improved = val_rmse < best_val_rmse - config.min_delta
        if improved:
            best_val_rmse = val_rmse
            best_epoch = epoch
            best_state = {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            }
            stale_epochs = 0
        else:
            stale_epochs += 1
        scheduler.step(val_rmse)
        if stale_epochs >= config.patience:
            early_stopped = True
            if logger is not None and verbose:
                logger.info("Early stopping triggered at epoch %d.", epoch)
            break

    if not best_state:
        raise RuntimeError("Training did not produce a valid best model state")
    elapsed = time.perf_counter() - start_time
    return FitResult(
        best_state, best_epoch, best_val_rmse, history, early_stopped, elapsed
    )


def default_training_config() -> TrainingConfig:
    return TrainingConfig(
        batch_size=BATCH_SIZE,
        learning_rate=LEARNING_RATE,
        weight_decay=WEIGHT_DECAY,
        max_epochs=MAX_EPOCHS,
        patience=PATIENCE,
        min_delta=MIN_DELTA,
        grad_clip_norm=GRAD_CLIP_NORM,
        lr_factor=LR_FACTOR,
        lr_patience=LR_PATIENCE,
        min_lr=MIN_LR,
    )


def run_optuna(
    split_records: Dict[str, List[SampleRecord]],
    graph_map: Dict[str, Data],
    target_mean: float,
    target_std: float,
    device: torch.device,
    logger: logging.Logger,
) -> Dict[str, Any]:
    try:
        import optuna
    except ImportError as exc:
        raise ImportError(
            "Optuna is required when USE_OPTUNA=True. Install it with: pip install optuna"
        ) from exc

    optuna.logging.set_verbosity(optuna.logging.WARNING)

    def objective(trial: Any) -> float:
        trial_seed = SEED
        set_global_seed(trial_seed)
        trial_model_config = model_config()
        trial_model_config.update(
            {
                "embedding_dim": trial.suggest_categorical("embedding_dim", [16, 32]),
                "n1": trial.suggest_categorical("n1", [64, 96, 128]),
                "n2": trial.suggest_categorical("n2", [32, 64]),
                "n3": trial.suggest_categorical("n3", [16, 32]),
                "num_blocks": trial.suggest_int("num_blocks", 2, 4),
                "dropout": trial.suggest_float("dropout", 0.0, 0.25),
            }
        )
        batch_size = trial.suggest_categorical("batch_size", [32, 64, 96])
        loaders = make_loaders(
            split_records, graph_map, target_mean, target_std, batch_size
        )
        model = MEGNetRegressor(**trial_model_config).to(device)
        config = TrainingConfig(
            batch_size=batch_size,
            learning_rate=trial.suggest_float(
                "learning_rate", 1.0e-4, 3.0e-3, log=True
            ),
            weight_decay=trial.suggest_float("weight_decay", 1.0e-8, 1.0e-3, log=True),
            max_epochs=OPTUNA_MAX_EPOCHS,
            patience=OPTUNA_PATIENCE,
            min_delta=MIN_DELTA,
            grad_clip_norm=GRAD_CLIP_NORM,
            lr_factor=LR_FACTOR,
            lr_patience=LR_PATIENCE,
            min_lr=MIN_LR,
        )
        try:
            result = fit_model(
                model,
                loaders,
                config,
                device,
                target_mean,
                target_std,
                logger=None,
                verbose=False,
            )
            return result.best_val_rmse
        except torch.cuda.OutOfMemoryError:
            if device.type == "cuda":
                torch.cuda.empty_cache()
            raise optuna.TrialPruned("CUDA out of memory")
        finally:
            del model
            if device.type == "cuda":
                torch.cuda.empty_cache()

    study = optuna.create_study(
        direction="minimize", sampler=optuna.samplers.TPESampler(seed=SEED)
    )
    study.optimize(objective, n_trials=OPTUNA_N_TRIALS, timeout=OPTUNA_TIMEOUT)
    best = {
        "best_val_rmse_eV": float(study.best_value),
        "best_params": study.best_params,
    }
    output_path = TABLE_DIR / f"{MODEL_NAME}_optuna_best_params.json"
    with output_path.open("w", encoding="utf-8") as handle:
        json.dump(best, handle, indent=2)
    logger.info("Optuna best validation RMSE: %.6f eV", study.best_value)
    logger.info("Optuna best parameters: %s", study.best_params)
    return best


def apply_optuna_parameters(
    base_model_config: Dict[str, Any],
    base_training_config: TrainingConfig,
    best_params: Dict[str, Any],
) -> Tuple[Dict[str, Any], TrainingConfig]:
    updated_model_config = copy.deepcopy(base_model_config)
    updated_training_config = copy.deepcopy(base_training_config)
    for key in ["embedding_dim", "n1", "n2", "n3", "num_blocks", "dropout"]:
        if key in best_params:
            updated_model_config[key] = best_params[key]
    for key in ["batch_size", "learning_rate", "weight_decay"]:
        if key in best_params:
            setattr(updated_training_config, key, best_params[key])
    return updated_model_config, updated_training_config


def calculate_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    return {
        "mae_eV": float(mean_absolute_error(y_true, y_pred)),
        "rmse_eV": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "r2": float(r2_score(y_true, y_pred)),
    }


def save_checkpoint(
    path: Path,
    model_state: Dict[str, torch.Tensor],
    used_model_config: Dict[str, Any],
    used_training_config: TrainingConfig,
    target_mean: float,
    target_std: float,
    fit_result: FitResult,
) -> None:
    checkpoint = {
        "model_name": MODEL_NAME,
        "run_version": RUN_VERSION,
        "model_state_dict": model_state,
        "model_config": used_model_config,
        "training_config": vars(used_training_config),
        "graph_config": graph_config(),
        "best_epoch": fit_result.best_epoch,
        "best_val_rmse": fit_result.best_val_rmse,
        "seed": SEED,
        "target_name": TARGET_NAME,
        "target_unit": TARGET_UNIT,
        "normalization": {"mean": target_mean, "std": target_std},
        "history": fit_result.history,
    }
    torch.save(checkpoint, path)


def load_trained_model(
    checkpoint_path: Union[str, Path], device: Optional[torch.device] = None
) -> Tuple[MEGNetRegressor, Dict[str, Any]]:
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    try:
        checkpoint = torch.load(
            checkpoint_path, map_location=device, weights_only=False
        )
    except TypeError:
        checkpoint = torch.load(checkpoint_path, map_location=device)
    required = {"model_state_dict", "model_config", "graph_config", "normalization"}
    if not required.issubset(checkpoint):
        raise ValueError(
            f"Checkpoint is missing required keys: {sorted(required - set(checkpoint))}"
        )
    model = MEGNetRegressor(**checkpoint["model_config"]).to(device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    return model, checkpoint


@torch.inference_mode()
def predict_cifs(
    model: MEGNetRegressor,
    cif_paths: Sequence[Union[str, Path]],
    device: torch.device,
    normalization: Dict[str, float],
    graph_settings: Optional[Dict[str, Any]] = None,
    batch_size: int = 64,
) -> pd.DataFrame:
    graphs: List[Data] = []
    valid_names: List[str] = []
    errors: List[Dict[str, Any]] = []
    settings = graph_config() if graph_settings is None else graph_settings
    for cif_path in cif_paths:
        try:
            graph = build_crystal_graph(
                cif_path, cif_name=Path(cif_path).name, config=settings
            )
            graphs.append(graph)
            valid_names.append(Path(cif_path).name)
        except Exception as exc:
            errors.append(
                {
                    "cif_name": Path(cif_path).name,
                    f"predicted_{TARGET_NAME}_{TARGET_UNIT}": np.nan,
                    "error": f"{type(exc).__name__}: {exc}",
                }
            )
    rows: List[Dict[str, Any]] = []
    model.eval()
    for start in range(0, len(graphs), batch_size):
        batch_graphs = graphs[start : start + batch_size]
        batch = Batch.from_data_list(batch_graphs).to(device)
        with autocast_context(device):
            prediction = model(batch).float().cpu().numpy()
        prediction = normalized_to_raw(
            prediction, normalization["mean"], normalization["std"]
        )
        for name, value in zip(valid_names[start : start + batch_size], prediction):
            rows.append(
                {
                    "cif_name": name,
                    f"predicted_{TARGET_NAME}_{TARGET_UNIT}": float(value),
                    "error": "",
                }
            )
    rows.extend(errors)
    return pd.DataFrame(rows)


def save_prediction_dat(
    split_name: str,
    names: Sequence[str],
    y_true: np.ndarray,
    y_pred: np.ndarray,
) -> pd.DataFrame:
    frame = pd.DataFrame(
        {
            "cif_name": list(names),
            "true_bandgap_eV": y_true,
            "predicted_bandgap_eV": y_pred,
            "error_eV": y_pred - y_true,
            "absolute_error_eV": np.abs(y_pred - y_true),
        }
    )
    frame.to_csv(
        DAT_DIR / f"{MODEL_NAME}_parity_{split_name}.dat",
        sep="\t",
        index=False,
        float_format="%.10f",
    )
    return frame


def global_axis_limits(results: Dict[str, Dict[str, Any]]) -> Tuple[float, float]:
    all_values = np.concatenate(
        [
            np.concatenate([entry["y_true"], entry["y_pred"]])
            for entry in results.values()
        ]
    )
    minimum = float(np.min(all_values))
    maximum = float(np.max(all_values))
    span = maximum - minimum
    padding = 0.05 * span if span > 0 else 0.5
    return minimum - padding, maximum + padding


def style_axes(ax: plt.Axes) -> None:
    ax.tick_params(axis="both", labelsize=TICK_FONTSIZE)
    if USE_GRID:
        ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)


def plot_parity(
    split_name: str,
    y_true: np.ndarray,
    y_pred: np.ndarray,
    metrics: Dict[str, float],
    limits: Tuple[float, float],
    color: str,
) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    ax.scatter(y_true, y_pred, s=24, alpha=0.78, color=color, edgecolors="none")
    ax.plot(limits, limits, linestyle="--", color="black", linewidth=1.2, label="y = x")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(
        f"{MODEL_NAME}: {split_name.capitalize()} parity", fontsize=TITLE_FONTSIZE
    )
    ax.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    annotation = (
        f"MAE = {metrics['mae_eV']:.4f} eV\n"
        f"RMSE = {metrics['rmse_eV']:.4f} eV\n"
        f"R² = {metrics['r2']:.4f}"
    )
    ax.text(
        0.04,
        0.96,
        annotation,
        transform=ax.transAxes,
        va="top",
        fontsize=ANNOTATION_FONTSIZE,
        bbox={
            "boxstyle": "round",
            "facecolor": "white",
            "alpha": 0.8,
            "edgecolor": "0.7",
        },
    )
    style_axes(ax)
    ax.legend(fontsize=LEGEND_FONTSIZE)
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
    fig, ax = plt.subplots(figsize=(6.8, 6.2))
    for color, split_name in zip(COLORS[:3], ["train", "val", "test"]):
        entry = results[split_name]
        metrics = entry["metrics"]
        label = (
            f"{split_name.capitalize()} "
            f"(RMSE={metrics['rmse_eV']:.4f}, R²={metrics['r2']:.4f})"
        )
        ax.scatter(
            entry["y_true"],
            entry["y_pred"],
            s=22,
            alpha=0.72,
            color=color,
            label=label,
            edgecolors="none",
        )
    ax.plot(limits, limits, linestyle="--", color="black", linewidth=1.2, label="y = x")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"{MODEL_NAME}: Combined parity", fontsize=TITLE_FONTSIZE)
    ax.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    style_axes(ax)
    ax.legend(fontsize=LEGEND_FONTSIZE)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_parity_all.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(fig)


def plot_rmse_curve(history: Sequence[Dict[str, float]], best_epoch: int) -> None:
    frame = pd.DataFrame(history)
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
    ax.set_title(f"{MODEL_NAME}: RMSE learning curve", fontsize=TITLE_FONTSIZE)
    ax.set_xlabel("Epoch", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"RMSE ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    style_axes(ax)
    ax.legend(fontsize=LEGEND_FONTSIZE)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_rmse_curve.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(fig)


def save_all_predictions_and_plots(
    results: Dict[str, Dict[str, Any]], best_epoch: int, history: List[Dict[str, float]]
) -> None:
    limits = global_axis_limits(results)
    combined_frames: List[pd.DataFrame] = []
    for color, split_name in zip(COLORS[:3], ["train", "val", "test"]):
        entry = results[split_name]
        frame = save_prediction_dat(
            split_name, entry["names"], entry["y_true"], entry["y_pred"]
        )
        frame.insert(0, "split", split_name)
        combined_frames.append(frame)
        plot_parity(
            split_name,
            entry["y_true"],
            entry["y_pred"],
            entry["metrics"],
            limits,
            color,
        )
    pd.concat(combined_frames, ignore_index=True).to_csv(
        DAT_DIR / f"{MODEL_NAME}_parity_all.dat",
        sep="\t",
        index=False,
        float_format="%.10f",
    )
    plot_combined_parity(results, limits)
    plot_rmse_curve(history, best_epoch)


def save_metrics_table(results: Dict[str, Dict[str, Any]]) -> None:
    rows = []
    for split_name in ["train", "val", "test"]:
        entry = results[split_name]
        rows.append(
            {
                "split": split_name,
                "n_samples": len(entry["y_true"]),
                **entry["metrics"],
            }
        )
    pd.DataFrame(rows).to_csv(
        TABLE_DIR / f"{MODEL_NAME}_metrics.dat",
        sep="\t",
        index=False,
        float_format="%.6f",
    )


def log_model_summary(model: nn.Module, logger: logging.Logger) -> None:
    total = sum(parameter.numel() for parameter in model.parameters())
    trainable = sum(
        parameter.numel() for parameter in model.parameters() if parameter.requires_grad
    )
    logger.info("Model structure:\n%s", repr(model))
    logger.info("Total parameters: %d", total)
    logger.info("Trainable parameters: %d", trainable)


def main() -> None:
    total_start = time.perf_counter()
    create_output_directories()
    logger = setup_logger(LOG_PATH)
    set_global_seed(SEED)
    device = get_device(logger)

    logger.info("Model name: %s", MODEL_NAME)
    logger.info("Run version: %s", RUN_VERSION)
    logger.info("Output root: %s", OUTPUT_ROOT.resolve())
    logger.info("Excel path: %s", Path(EXCEL_PATH).resolve())
    logger.info("CIF directory: %s", Path(CIF_DIR).resolve())
    logger.info("Excel mapping: first column=CIF filename, second column=target")
    logger.info("Target: %s (%s)", TARGET_NAME, TARGET_UNIT)
    logger.info("Seed: %d", SEED)
    logger.info(
        "Split ratios: train=%.2f, val=%.2f, test=%.2f",
        TRAIN_RATIO,
        VAL_RATIO,
        TEST_RATIO,
    )
    logger.info("Graph configuration: %s", graph_config())
    logger.info("Loss function: MSELoss")
    logger.info("Optuna enabled: %s", USE_OPTUNA)
    logger.info("AMP requested: %s", USE_AMP)
    logger.info(
        "BF16 autocast active: %s",
        bool(USE_AMP and device.type == "cuda" and torch.cuda.is_bf16_supported()),
    )

    data_start = time.perf_counter()
    records, excel_failures = load_excel_records(EXCEL_PATH, CIF_DIR, logger)
    records, graph_map, graph_failures = validate_and_cache_records(
        records, graph_config(), logger
    )
    split_map = create_or_load_split(records, SPLIT_PATH, logger)
    split_records = {
        name: [record for record in records if split_map[record.cif_name] == name]
        for name in ["train", "val", "test"]
    }
    for split_name, items in split_records.items():
        logger.info("%s samples: %d", split_name.capitalize(), len(items))
    train_targets = np.asarray(
        [record.target for record in split_records["train"]], dtype=np.float64
    )
    target_mean = float(train_targets.mean())
    target_std = float(train_targets.std(ddof=0))
    if not np.isfinite(target_std) or target_std < 1.0e-12:
        raise ValueError("Training target standard deviation is zero or invalid")
    logger.info("Training target mean: %.10f %s", target_mean, TARGET_UNIT)
    logger.info("Training target standard deviation: %.10f %s", target_std, TARGET_UNIT)
    data_elapsed = time.perf_counter() - data_start
    logger.info(
        "Data processing time: %.3f s (%s)", data_elapsed, format_duration(data_elapsed)
    )

    used_model_config = model_config()
    used_training_config = default_training_config()
    if USE_OPTUNA:
        optuna_result = run_optuna(
            split_records, graph_map, target_mean, target_std, device, logger
        )
        used_model_config, used_training_config = apply_optuna_parameters(
            used_model_config, used_training_config, optuna_result["best_params"]
        )

    set_global_seed(SEED)
    loaders = make_loaders(
        split_records,
        graph_map,
        target_mean,
        target_std,
        used_training_config.batch_size,
    )
    model = MEGNetRegressor(**used_model_config).to(device)
    log_model_summary(model, logger)
    logger.info("Model configuration: %s", used_model_config)
    logger.info("Training configuration: %s", vars(used_training_config))
    logger.info("Optimizer: Adam")
    logger.info(
        "Scheduler: ReduceLROnPlateau(factor=%s, patience=%s, min_lr=%s)",
        used_training_config.lr_factor,
        used_training_config.lr_patience,
        used_training_config.min_lr,
    )

    try:
        fit_result = fit_model(
            model,
            loaders,
            used_training_config,
            device,
            target_mean,
            target_std,
            logger=logger,
            verbose=True,
        )
    except torch.cuda.OutOfMemoryError as exc:
        raise RuntimeError(
            "CUDA ran out of memory. Reduce BATCH_SIZE or model dimensions in the user configuration."
        ) from exc

    save_checkpoint(
        CHECKPOINT_PATH,
        fit_result.best_state,
        used_model_config,
        used_training_config,
        target_mean,
        target_std,
        fit_result,
    )
    logger.info("Best epoch: %d", fit_result.best_epoch)
    logger.info("Best validation RMSE: %.6f %s", fit_result.best_val_rmse, TARGET_UNIT)
    logger.info("Early stopping triggered: %s", fit_result.early_stopped)
    logger.info(
        "Training time: %.3f s (%s)",
        fit_result.elapsed_seconds,
        format_duration(fit_result.elapsed_seconds),
    )
    logger.info("Best checkpoint: %s", CHECKPOINT_PATH)

    evaluation_start = time.perf_counter()
    best_model, checkpoint = load_trained_model(CHECKPOINT_PATH, device)
    criterion = nn.MSELoss()
    results: Dict[str, Dict[str, Any]] = {}
    for split_name in ["train", "val", "test"]:
        _, y_true, y_pred, names = evaluate_loader(
            best_model,
            loaders[split_name],
            criterion,
            device,
            target_mean,
            target_std,
        )
        metrics = calculate_metrics(y_true, y_pred)
        results[split_name] = {
            "y_true": y_true,
            "y_pred": y_pred,
            "names": names,
            "metrics": metrics,
        }
        logger.info(
            "%s | MAE %.6f eV | RMSE %.6f eV | R2 %.6f",
            split_name.capitalize(),
            metrics["mae_eV"],
            metrics["rmse_eV"],
            metrics["r2"],
        )

    save_metrics_table(results)
    save_all_predictions_and_plots(
        results, checkpoint["best_epoch"], checkpoint["history"]
    )
    evaluation_elapsed = time.perf_counter() - evaluation_start
    logger.info(
        "Evaluation and plotting time: %.3f s (%s)",
        evaluation_elapsed,
        format_duration(evaluation_elapsed),
    )
    logger.info("Invalid samples: %d", len(excel_failures) + len(graph_failures))
    for sample_name, reason in excel_failures + graph_failures:
        logger.info("Invalid sample | %s | %s", sample_name, reason)
    total_elapsed = time.perf_counter() - total_start
    logger.info(
        "Total runtime: %.3f s (%s)", total_elapsed, format_duration(total_elapsed)
    )


if __name__ == "__main__":
    main()

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import random
import re
import sys
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
import torch.nn.functional as F
from pymatgen.core import Structure
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from torch import Tensor, nn
from torch_geometric.data import Batch, Data
from torch_geometric.loader import DataLoader
from torch_geometric.nn import MessagePassing, global_mean_pool
from torch_geometric.utils import scatter
from tqdm.auto import tqdm


# ============================== User configuration ==============================

MODEL_NAME = "eComFormer"  # "iComFormer" or "eComFormer"
RUN_VERSION = "v1"

EXCEL_PATH = "./data.xlsx"
CIF_DIR = "./cif"

TARGET_NAME = "bandgap"
TARGET_UNIT = "eV"

SEED = 42
TRAIN_RATIO = 0.8
VAL_RATIO = 0.1
TEST_RATIO = 0.1

BATCH_SIZE = 32
MAX_EPOCHS = 300
LEARNING_RATE = 1.0e-3
WEIGHT_DECAY = 1.0e-5
PATIENCE = 50
MIN_DELTA = 1.0e-5
NUM_WORKERS = 0
PIN_MEMORY = True
USE_AMP = False

CUTOFF = 4.0
MAX_NEIGHBORS = 25
MAX_CUTOFF_EXPANSIONS = 8
CACHE_GRAPHS = True

ATOM_FEATURES = 92
HIDDEN_DIM = 256
FC_DIM = 256
DROPOUT = 0.0
ATTENTION_HEADS = 1
EQUIVARIANT_SCALARS = 64
EQUIVARIANT_VECTORS = 8

USE_OPTUNA = False
OPTUNA_N_TRIALS = 30
OPTUNA_TIMEOUT = None
OPTUNA_MAX_EPOCHS = 120

FIG_DPI = 600
TITLE_FONTSIZE = 16
LABEL_FONTSIZE = 14
TICK_FONTSIZE = 12
LEGEND_FONTSIZE = 11
ANNOTATION_FONTSIZE = 11
USE_GRID = True
COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:red"]

# ================================================================================


SUPPORTED_MODELS = ("iComFormer", "eComFormer")
SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_ROOT = SCRIPT_DIR / f"{MODEL_NAME}_{RUN_VERSION}"
FIGURE_DIR = OUTPUT_ROOT / "figure"
DAT_DIR = OUTPUT_ROOT / "dat"
TABLE_DIR = OUTPUT_ROOT / "table"
LOG_DIR = OUTPUT_ROOT / "log"
SPLIT_DIR = OUTPUT_ROOT / "split"
CACHE_DIR = OUTPUT_ROOT / "cache"
CHECKPOINT_PATH = OUTPUT_ROOT / f"{MODEL_NAME}_best.pt"
SPLIT_PATH = SPLIT_DIR / f"{MODEL_NAME}_split.csv"
LOG_PATH = LOG_DIR / f"{MODEL_NAME}_training.log"


@dataclass
class ModelConfig:
    name: str
    atom_features: int = ATOM_FEATURES
    hidden_dim: int = HIDDEN_DIM
    fc_dim: int = FC_DIM
    heads: int = ATTENTION_HEADS
    dropout: float = DROPOUT
    equi_scalars: int = EQUIVARIANT_SCALARS
    equi_vectors: int = EQUIVARIANT_VECTORS


@dataclass
class GraphConfig:
    cutoff: float = CUTOFF
    max_neighbors: int = MAX_NEIGHBORS
    max_cutoff_expansions: int = MAX_CUTOFF_EXPANSIONS


def validate_configuration() -> None:
    if MODEL_NAME not in SUPPORTED_MODELS:
        raise ValueError(
            f"MODEL_NAME must be one of {SUPPORTED_MODELS}, got {MODEL_NAME!r}."
        )
    if not math.isclose(TRAIN_RATIO + VAL_RATIO + TEST_RATIO, 1.0, abs_tol=1e-9):
        raise ValueError("TRAIN_RATIO + VAL_RATIO + TEST_RATIO must equal 1.0.")
    if min(TRAIN_RATIO, VAL_RATIO, TEST_RATIO) <= 0:
        raise ValueError("All split ratios must be positive.")
    if BATCH_SIZE < 1 or MAX_EPOCHS < 1 or MAX_NEIGHBORS < 1 or CUTOFF <= 0:
        raise ValueError("Batch size, epochs, neighbors, and cutoff must be positive.")
    if HIDDEN_DIM % ATTENTION_HEADS != 0:
        raise ValueError("HIDDEN_DIM must be divisible by ATTENTION_HEADS.")


def create_directories() -> None:
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


def setup_logger() -> logging.Logger:
    logger = logging.getLogger("ComFormerBandgap")
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


def elapsed_string(seconds: float) -> str:
    seconds = max(0.0, float(seconds))
    hours, remainder = divmod(int(round(seconds)), 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def get_device(logger: Optional[logging.Logger] = None) -> torch.device:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    if logger is not None:
        logger.info("Device: %s", device)
        logger.info("PyTorch version: %s", torch.__version__)
        if device.type == "cuda":
            props = torch.cuda.get_device_properties(0)
            logger.info("GPU: %s", torch.cuda.get_device_name(0))
            logger.info("GPU memory: %.2f GiB", props.total_memory / 1024**3)
            logger.info("CUDA runtime: %s", torch.version.cuda)
    return device


def load_excel_records(logger: logging.Logger) -> List[Dict[str, Any]]:
    excel_path = (
        (SCRIPT_DIR / EXCEL_PATH).resolve()
        if not Path(EXCEL_PATH).is_absolute()
        else Path(EXCEL_PATH)
    )
    cif_root = (
        (SCRIPT_DIR / CIF_DIR).resolve()
        if not Path(CIF_DIR).is_absolute()
        else Path(CIF_DIR)
    )
    if not excel_path.is_file():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    if not cif_root.is_dir():
        raise FileNotFoundError(f"CIF directory not found: {cif_root}")
    frame = pd.read_excel(excel_path)
    if frame.shape[1] < 2:
        raise ValueError("The Excel file must contain at least two columns.")
    records: List[Dict[str, Any]] = []
    for row_index, row in frame.iloc[:, :2].iterrows():
        raw_name, raw_target = row.iloc[0], row.iloc[1]
        if pd.isna(raw_name) or not str(raw_name).strip():
            logger.warning("Skipping Excel row %d: empty CIF filename.", row_index + 2)
            continue
        cif_name = f"{int(raw_name)}.cif"
        cif_path = Path(cif_name)
        if not cif_path.is_absolute():
            cif_path = cif_root / cif_path
        try:
            target = float(raw_target)
        except (TypeError, ValueError):
            logger.warning(
                "Skipping %s: target is not numeric (%r).", cif_name, raw_target
            )
            continue
        if not np.isfinite(target):
            logger.warning("Skipping %s: target is not finite.", cif_name)
            continue
        if not cif_path.is_file():
            logger.warning("Skipping %s: CIF file not found at %s.", cif_name, cif_path)
            continue
        records.append(
            {
                "cif_name": cif_name,
                "cif_path": str(cif_path.resolve()),
                "target": target,
            }
        )
    if len(records) < 10:
        raise ValueError(
            f"Only {len(records)} valid Excel/file records were found; at least 10 are required."
        )
    duplicate_names = pd.Series([r["cif_name"] for r in records]).duplicated(keep=False)
    if duplicate_names.any():
        duplicates = sorted(
            set(np.asarray([r["cif_name"] for r in records])[duplicate_names])
        )
        raise ValueError(f"Duplicate CIF names are not allowed: {duplicates[:10]}")
    logger.info("Excel path: %s", excel_path)
    logger.info("CIF directory: %s", cif_root)
    logger.info("Excel columns: first column=CIF filename, second column=target")
    logger.info("Valid tabular records: %d", len(records))
    return records


def _canonical_lattice_vectors(lattice_matrix: np.ndarray) -> np.ndarray:
    candidates: List[Tuple[float, Tuple[int, int, int], np.ndarray]] = []
    for i in range(-2, 3):
        for j in range(-2, 3):
            for k in range(-2, 3):
                if i == 0 and j == 0 and k == 0:
                    continue
                coeff = np.asarray([i, j, k], dtype=float)
                vector = coeff @ lattice_matrix
                candidates.append((float(np.linalg.norm(vector)), (i, j, k), vector))
    candidates.sort(key=lambda item: (round(item[0], 10), item[1]))
    chosen: List[np.ndarray] = []
    for _, _, vector in candidates:
        if not chosen:
            chosen.append(vector)
        elif len(chosen) == 1 and np.linalg.norm(np.cross(chosen[0], vector)) > 1e-7:
            chosen.append(vector)
        elif (
            len(chosen) == 2
            and abs(np.linalg.det(np.stack([chosen[0], chosen[1], vector]))) > 1e-7
        ):
            chosen.append(vector)
            break
    if len(chosen) != 3:
        raise ValueError(
            "Could not determine three independent lattice reference vectors."
        )
    refs = np.stack(chosen).astype(np.float64)
    if np.dot(refs[0], refs[1]) < 0:
        refs[1] *= -1
    if np.dot(refs[0], refs[2]) < 0:
        refs[2] *= -1
    if np.linalg.det(refs) < 0:
        refs *= -1
    return refs.astype(np.float32)


def build_graph_from_structure(structure: Structure, graph_config: GraphConfig) -> Data:
    if len(structure) == 0:
        raise ValueError("The structure has no atomic sites.")
    cutoff = float(graph_config.cutoff)
    selected: Optional[List[np.ndarray]] = None
    distances: Optional[np.ndarray] = None
    centers = points = offsets = None
    for _ in range(graph_config.max_cutoff_expansions + 1):
        centers, points, offsets, distances = structure.get_neighbor_list(cutoff)
        groups = [np.where(centers == index)[0] for index in range(len(structure))]
        if groups and min(len(group) for group in groups) >= graph_config.max_neighbors:
            selected = groups
            break
        cutoff *= 1.5
    if (
        selected is None
        or centers is None
        or points is None
        or offsets is None
        or distances is None
    ):
        raise ValueError(
            f"Unable to find {graph_config.max_neighbors} neighbors per atom after "
            f"{graph_config.max_cutoff_expansions + 1} cutoff attempts."
        )
    chosen_indices: List[int] = []
    for group in selected:
        order = group[np.argsort(distances[group], kind="stable")]
        kth_distance = distances[order[graph_config.max_neighbors - 1]]
        keep = order[distances[order] <= kth_distance + 1e-8]
        chosen_indices.extend(keep.tolist())
    chosen = np.asarray(chosen_indices, dtype=np.int64)
    src = centers[chosen].astype(np.int64)
    dst = points[chosen].astype(np.int64)
    image = offsets[chosen].astype(np.float64)
    lattice_matrix = np.asarray(structure.lattice.matrix, dtype=np.float64)
    cart_coords = np.asarray(structure.cart_coords, dtype=np.float64)
    edge_vec = cart_coords[dst] + image @ lattice_matrix - cart_coords[src]
    edge_length = np.linalg.norm(edge_vec, axis=1)
    valid = edge_length > 1e-8
    src, dst, edge_vec = src[valid], dst[valid], edge_vec[valid]
    if len(src) == 0:
        raise ValueError("Graph construction produced no valid edges.")
    atomic_numbers = []
    for site in structure:
        if not site.is_ordered:
            raise ValueError("Disordered/partially occupied sites are not supported.")
        atomic_numbers.append(int(site.specie.Z))
    if max(atomic_numbers) > 118 or min(atomic_numbers) < 1:
        raise ValueError("Atomic numbers must be in the range 1..118.")
    references = _canonical_lattice_vectors(lattice_matrix)
    edge_ref = np.repeat(references[None, :, :], len(src), axis=0)
    return Data(
        z=torch.tensor(atomic_numbers, dtype=torch.long),
        edge_index=torch.tensor(np.stack([src, dst]), dtype=torch.long),
        edge_vec=torch.tensor(edge_vec, dtype=torch.float32),
        edge_ref=torch.tensor(edge_ref, dtype=torch.float32),
        num_nodes=len(structure),
    )


def _cache_path(cif_path: str, graph_config: GraphConfig) -> Path:
    path = Path(cif_path)
    stat = path.stat()
    signature = json.dumps(
        {
            "path": str(path.resolve()),
            "size": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
            "graph": asdict(graph_config),
            "cache_format": 2,
        },
        sort_keys=True,
    )
    digest = hashlib.sha256(signature.encode("utf-8")).hexdigest()[:20]
    safe_stem = re.sub(r"[^A-Za-z0-9_.-]+", "_", path.stem)[:80]
    return CACHE_DIR / f"{safe_stem}_{digest}.pt"


def torch_load_full(path: os.PathLike[str] | str, map_location: Any) -> Any:
    try:
        return torch.load(path, map_location=map_location, weights_only=False)
    except TypeError:
        return torch.load(path, map_location=map_location)


def load_or_build_graph(cif_path: str, graph_config: GraphConfig) -> Data:
    cache_path = _cache_path(cif_path, graph_config)
    if CACHE_GRAPHS and cache_path.is_file():
        try:
            return torch_load_full(cache_path, map_location="cpu")
        except Exception:
            cache_path.unlink(missing_ok=True)
    structure = Structure.from_file(cif_path)
    graph = build_graph_from_structure(structure, graph_config)
    if CACHE_GRAPHS:
        torch.save(graph, cache_path)
    return graph


def validate_and_build_records(
    records: Sequence[Dict[str, Any]], graph_config: GraphConfig, logger: logging.Logger
) -> List[Dict[str, Any]]:
    valid: List[Dict[str, Any]] = []
    for record in tqdm(records, desc="Building periodic graphs", unit="cif"):
        try:
            graph = load_or_build_graph(record["cif_path"], graph_config)
            item = dict(record)
            item["graph"] = graph
            valid.append(item)
        except Exception as exc:
            logger.warning(
                "Skipping %s: graph construction failed: %s", record["cif_name"], exc
            )
    if len(valid) < 10:
        raise ValueError(
            f"Only {len(valid)} CIF files produced valid graphs; at least 10 are required."
        )
    logger.info("Valid periodic graphs: %d", len(valid))
    return valid


def split_records(
    records: Sequence[Dict[str, Any]], logger: logging.Logger
) -> Dict[str, List[Dict[str, Any]]]:
    lookup = {record["cif_name"]: record for record in records}
    if SPLIT_PATH.is_file():
        split_frame = pd.read_csv(SPLIT_PATH)
        required = {"cif_name", "target", "split"}
        if required.issubset(split_frame.columns) and set(
            split_frame["cif_name"]
        ) == set(lookup):
            if set(split_frame["split"]) <= {"train", "val", "test"}:
                output = {
                    name: [
                        lookup[cif]
                        for cif in split_frame.loc[
                            split_frame["split"] == name, "cif_name"
                        ]
                    ]
                    for name in ("train", "val", "test")
                }
                if all(output.values()):
                    logger.info("Reused split file: %s", SPLIT_PATH)
                    return output
        logger.warning("Existing split file is incompatible and will be replaced.")
    indices = np.arange(len(records))
    train_indices, remainder = train_test_split(
        indices, train_size=TRAIN_RATIO, random_state=SEED, shuffle=True
    )
    relative_val = VAL_RATIO / (VAL_RATIO + TEST_RATIO)
    val_indices, test_indices = train_test_split(
        remainder, train_size=relative_val, random_state=SEED, shuffle=True
    )
    mapping = {
        "train": [records[int(index)] for index in train_indices],
        "val": [records[int(index)] for index in val_indices],
        "test": [records[int(index)] for index in test_indices],
    }
    rows = []
    for split_name, subset in mapping.items():
        rows.extend(
            {
                "cif_name": item["cif_name"],
                "target": item["target"],
                "split": split_name,
            }
            for item in subset
        )
    pd.DataFrame(rows).to_csv(SPLIT_PATH, index=False)
    logger.info("Created split file: %s", SPLIT_PATH)
    return mapping


class BandgapDataset(torch.utils.data.Dataset):
    def __init__(
        self, records: Sequence[Dict[str, Any]], target_mean: float, target_std: float
    ):
        self.records = list(records)
        self.target_mean = float(target_mean)
        self.target_std = float(target_std)

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, index: int) -> Data:
        record = self.records[index]
        graph = record["graph"].clone()
        normalized = (float(record["target"]) - self.target_mean) / self.target_std
        graph.y = torch.tensor([normalized], dtype=torch.float32)
        graph.target_raw = torch.tensor([float(record["target"])], dtype=torch.float32)
        graph.cif_name = record["cif_name"]
        return graph


class PredictionDataset(torch.utils.data.Dataset):
    def __init__(self, graphs: Sequence[Data], names: Sequence[str]):
        self.graphs = list(graphs)
        self.names = list(names)

    def __len__(self) -> int:
        return len(self.graphs)

    def __getitem__(self, index: int) -> Data:
        graph = self.graphs[index].clone()
        graph.cif_name = self.names[index]
        return graph


class SafeBatchNorm1d(nn.BatchNorm1d):
    def forward(self, x: Tensor) -> Tensor:
        if self.training and x.shape[0] <= 1:
            return F.batch_norm(
                x,
                self.running_mean,
                self.running_var,
                self.weight,
                self.bias,
                False,
                0.0,
                self.eps,
            )
        return super().forward(x)


class RBFExpansion(nn.Module):
    def __init__(self, vmin: float, vmax: float, bins: int):
        super().__init__()
        self.register_buffer("centers", torch.linspace(vmin, vmax, bins))
        spacing = float((vmax - vmin) / max(1, bins - 1))
        self.gamma = 1.0 / max(spacing, 1e-8)

    def forward(self, values: Tensor) -> Tensor:
        return torch.exp(-self.gamma * (values.unsqueeze(-1) - self.centers) ** 2)


class ComformerConv(MessagePassing):
    def __init__(self, channels: int, heads: int = 1, dropout: float = 0.0):
        super().__init__(aggr="add", node_dim=0)
        if channels % heads != 0:
            raise ValueError("channels must be divisible by heads")
        self.channels = channels
        self.heads = heads
        self.head_dim = channels // heads
        self.dropout = dropout
        self.lin_query = nn.Linear(channels, channels)
        self.lin_key = nn.Linear(channels, channels)
        self.lin_value = nn.Linear(channels, channels)
        self.lin_edge = nn.Linear(channels, channels)
        self.key_update = nn.Sequential(
            nn.Linear(3 * self.head_dim, self.head_dim),
            nn.SiLU(),
            nn.Linear(self.head_dim, self.head_dim),
        )
        self.message_update = nn.Sequential(
            nn.Linear(3 * self.head_dim, self.head_dim),
            nn.SiLU(),
            nn.Linear(self.head_dim, self.head_dim),
        )
        self.output = nn.Linear(channels, channels)
        self.att_norm = SafeBatchNorm1d(channels)
        self.out_norm = SafeBatchNorm1d(channels)
        self.activation = nn.Softplus()

    def forward(self, x: Tensor, edge_index: Tensor, edge_attr: Tensor) -> Tensor:
        query = self.lin_query(x).view(-1, self.heads, self.head_dim)
        key = self.lin_key(x).view(-1, self.heads, self.head_dim)
        value = self.lin_value(x).view(-1, self.heads, self.head_dim)
        output = self.propagate(
            edge_index, query=query, key=key, value=value, edge_attr=edge_attr
        )
        output = self.output(output.reshape(-1, self.channels))
        output = F.dropout(output, p=self.dropout, training=self.training)
        return self.activation(x + self.out_norm(output))

    def message(
        self,
        query_i: Tensor,
        key_i: Tensor,
        key_j: Tensor,
        value_i: Tensor,
        value_j: Tensor,
        edge_attr: Tensor,
    ) -> Tensor:
        edge = self.lin_edge(edge_attr).view(-1, self.heads, self.head_dim)
        updated_key = self.key_update(torch.cat([key_i, key_j, edge], dim=-1))
        attention = query_i * updated_key / math.sqrt(self.head_dim)
        attention = torch.sigmoid(
            self.att_norm(attention.reshape(-1, self.channels)).view_as(attention)
        )
        message = self.message_update(torch.cat([value_i, value_j, edge], dim=-1))
        return message * attention


class ComformerEdgeUpdate(nn.Module):
    def __init__(self, channels: int, heads: int = 1, dropout: float = 0.0):
        super().__init__()
        if channels % heads != 0:
            raise ValueError("channels must be divisible by heads")
        self.channels = channels
        self.heads = heads
        self.head_dim = channels // heads
        self.dropout = dropout
        self.query = nn.Linear(channels, channels)
        self.key_center = nn.Linear(channels, channels)
        self.value_center = nn.Linear(channels, channels)
        self.key_refs = nn.ModuleList([nn.Linear(channels, channels) for _ in range(3)])
        self.value_refs = nn.ModuleList(
            [nn.Linear(channels, channels) for _ in range(3)]
        )
        self.angle_projection = nn.Linear(channels, channels, bias=False)
        self.key_update = nn.Sequential(
            nn.Linear(3 * self.head_dim, self.head_dim),
            nn.SiLU(),
            nn.Linear(self.head_dim, self.head_dim),
        )
        self.message_update = nn.Sequential(
            nn.Linear(3 * self.head_dim, self.head_dim),
            nn.SiLU(),
            nn.Linear(self.head_dim, self.head_dim),
        )
        self.output = nn.Linear(channels, channels)
        self.att_norm = SafeBatchNorm1d(channels)
        self.out_norm = SafeBatchNorm1d(channels)
        self.activation = nn.Softplus()

    def forward(self, edge: Tensor, ref_length: Tensor, ref_angle: Tensor) -> Tensor:
        edge_count = edge.shape[0]
        query = (
            self.query(edge)
            .view(edge_count, 1, self.heads, self.head_dim)
            .expand(-1, 3, -1, -1)
        )
        key_x = (
            self.key_center(edge)
            .view(edge_count, 1, self.heads, self.head_dim)
            .expand(-1, 3, -1, -1)
        )
        value_x = (
            self.value_center(edge)
            .view(edge_count, 1, self.heads, self.head_dim)
            .expand(-1, 3, -1, -1)
        )
        key_y = torch.stack(
            [
                layer(ref_length[:, index]).view(edge_count, self.heads, self.head_dim)
                for index, layer in enumerate(self.key_refs)
            ],
            dim=1,
        )
        value_y = torch.stack(
            [
                layer(ref_length[:, index]).view(edge_count, self.heads, self.head_dim)
                for index, layer in enumerate(self.value_refs)
            ],
            dim=1,
        )
        angle = self.angle_projection(ref_angle).view(
            edge_count, 3, self.heads, self.head_dim
        )
        key = self.key_update(torch.cat([key_x, key_y, angle], dim=-1))
        attention = query * key / math.sqrt(self.head_dim)
        attention = torch.sigmoid(
            self.att_norm(attention.reshape(-1, self.channels)).view_as(attention)
        )
        message = (
            self.message_update(torch.cat([value_x, value_y, angle], dim=-1))
            * attention
        )
        update = self.output(message.reshape(edge_count, 3, self.channels).sum(dim=1))
        update = F.dropout(update, p=self.dropout, training=self.training)
        return self.activation(edge + self.out_norm(update))


_E3NN_O3 = None


def require_e3nn():
    global _E3NN_O3
    if _E3NN_O3 is None:
        try:
            from e3nn import o3
        except ImportError as exc:
            raise ImportError(
                "eComFormer requires e3nn. Install it with: pip install e3nn"
            ) from exc
        _E3NN_O3 = o3
    return _E3NN_O3


class TensorProductConvLayer(nn.Module):
    def __init__(
        self, in_irreps: str, out_irreps: str, edge_channels: int, residual: bool
    ):
        super().__init__()
        o3 = require_e3nn()
        self.in_irreps = o3.Irreps(in_irreps)
        self.sh_irreps = o3.Irreps("1x0e + 1x1o + 1x2e")
        self.out_irreps = o3.Irreps(out_irreps)
        self.residual = residual
        self.tensor_product = o3.FullyConnectedTensorProduct(
            self.in_irreps, self.sh_irreps, self.out_irreps, shared_weights=False
        )
        self.weight_network = nn.Sequential(
            nn.Linear(edge_channels, edge_channels),
            nn.Softplus(),
            nn.Linear(edge_channels, self.tensor_product.weight_numel),
        )

    def forward(
        self, node_attr: Tensor, edge_index: Tensor, edge_attr: Tensor, edge_sh: Tensor
    ) -> Tensor:
        src, dst = edge_index
        messages = self.tensor_product(
            node_attr[src], edge_sh, self.weight_network(edge_attr)
        )
        output = scatter(
            messages, dst, dim=0, dim_size=node_attr.shape[0], reduce="mean"
        )
        if self.residual:
            if output.shape[-1] >= node_attr.shape[-1]:
                output = output + F.pad(
                    node_attr, (0, output.shape[-1] - node_attr.shape[-1])
                )
            else:
                output = output + node_attr[:, : output.shape[-1]]
        return output


class ComformerEquivariantUpdate(nn.Module):
    def __init__(self, channels: int, edge_channels: int, ns: int, nv: int):
        super().__init__()
        require_e3nn()
        hidden_irreps = f"{ns}x0e + {nv}x1o + {nv}x2e"
        self.ns = ns
        self.node_in = nn.Linear(channels, ns)
        self.skip = nn.Linear(channels, channels)
        self.layer1 = TensorProductConvLayer(
            f"{ns}x0e", hidden_irreps, edge_channels, residual=True
        )
        self.layer2 = TensorProductConvLayer(
            hidden_irreps, f"{ns}x0e", edge_channels, residual=False
        )
        self.norm = SafeBatchNorm1d(ns)
        self.node_out = nn.Linear(ns, channels)
        self.activation = nn.Softplus()

    def forward(
        self, data: Batch, node_features: Tensor, edge_features: Tensor
    ) -> Tensor:
        o3 = require_e3nn()
        edge_sh = o3.spherical_harmonics(
            self.layer1.sh_irreps,
            data.edge_vec,
            normalize=True,
            normalization="component",
        )
        scalar = self.node_in(node_features)
        scalar = self.layer1(scalar, data.edge_index, edge_features, edge_sh)
        scalar = self.layer2(scalar, data.edge_index, edge_features, edge_sh)
        scalar = self.node_out(self.activation(self.norm(scalar)))
        return self.activation(scalar + self.skip(node_features))


class BaseComFormer(nn.Module):
    def __init__(self, config: ModelConfig):
        super().__init__()
        self.config = config
        self.atomic_embedding = nn.Embedding(119, config.atom_features, padding_idx=0)
        self.atom_projection = nn.Linear(config.atom_features, config.hidden_dim)
        self.distance_rbf = nn.Sequential(
            RBFExpansion(-4.0, 0.0, config.hidden_dim),
            nn.Linear(config.hidden_dim, config.hidden_dim),
            nn.Softplus(),
        )
        self.readout = nn.Sequential(
            nn.Linear(config.hidden_dim, config.fc_dim),
            nn.SiLU(),
            nn.Dropout(config.dropout),
            nn.Linear(config.fc_dim, 1),
        )

    def initial_features(self, data: Batch) -> Tuple[Tensor, Tensor]:
        node = self.atom_projection(self.atomic_embedding(data.z))
        distance = torch.linalg.vector_norm(data.edge_vec, dim=-1).clamp_min(1e-8)
        edge = self.distance_rbf(-0.75 / distance)
        return node, edge

    def crystal_output(self, node_features: Tensor, batch_index: Tensor) -> Tensor:
        return self.readout(global_mean_pool(node_features, batch_index)).view(-1)


class iComFormer(BaseComFormer):
    def __init__(self, config: ModelConfig):
        super().__init__(config)
        self.angle_rbf = nn.Sequential(
            RBFExpansion(-1.0, 1.0, config.hidden_dim),
            nn.Linear(config.hidden_dim, config.hidden_dim),
            nn.Softplus(),
        )
        self.layers = nn.ModuleList(
            [
                ComformerConv(config.hidden_dim, config.heads, config.dropout)
                for _ in range(4)
            ]
        )
        self.edge_update = ComformerEdgeUpdate(
            config.hidden_dim, config.heads, config.dropout
        )

    def forward(self, data: Batch) -> Tensor:
        node, edge = self.initial_features(data)
        node = self.layers[0](node, data.edge_index, edge)
        ref_distance = torch.linalg.vector_norm(data.edge_ref, dim=-1).clamp_min(1e-8)
        ref_length = self.distance_rbf((-0.75 / ref_distance).reshape(-1)).reshape(
            data.edge_ref.shape[0], 3, -1
        )
        edge_unit = F.normalize(data.edge_vec, dim=-1)
        ref_unit = F.normalize(data.edge_ref, dim=-1)
        cosine = torch.sum(ref_unit * edge_unit.unsqueeze(1), dim=-1).clamp(-1.0, 1.0)
        ref_angle = self.angle_rbf(cosine.reshape(-1)).reshape(
            data.edge_ref.shape[0], 3, -1
        )
        edge = self.edge_update(edge, ref_length, ref_angle)
        for layer in self.layers[1:]:
            node = layer(node, data.edge_index, edge)
        return self.crystal_output(node, data.batch)


class eComFormer(BaseComFormer):
    def __init__(self, config: ModelConfig):
        super().__init__(config)
        self.layers = nn.ModuleList(
            [
                ComformerConv(config.hidden_dim, config.heads, config.dropout)
                for _ in range(3)
            ]
        )
        self.equivariant_update = ComformerEquivariantUpdate(
            config.hidden_dim,
            config.hidden_dim,
            config.equi_scalars,
            config.equi_vectors,
        )

    def forward(self, data: Batch) -> Tensor:
        node, edge = self.initial_features(data)
        node = self.layers[0](node, data.edge_index, edge)
        node = self.equivariant_update(data, node, edge)
        node = self.layers[1](node, data.edge_index, edge)
        node = self.layers[2](node, data.edge_index, edge)
        return self.crystal_output(node, data.batch)


def build_model(config: ModelConfig) -> nn.Module:
    if config.name == "iComFormer":
        return iComFormer(config)
    if config.name == "eComFormer":
        return eComFormer(config)
    raise ValueError(f"Unsupported model name: {config.name}")


def make_loaders(
    splits: Dict[str, List[Dict[str, Any]]],
    target_mean: float,
    target_std: float,
    batch_size: int,
) -> Dict[str, DataLoader]:
    generator = torch.Generator()
    generator.manual_seed(SEED)
    common = dict(
        num_workers=NUM_WORKERS, pin_memory=PIN_MEMORY and torch.cuda.is_available()
    )
    return {
        "train": DataLoader(
            BandgapDataset(splits["train"], target_mean, target_std),
            batch_size=batch_size,
            shuffle=True,
            generator=generator,
            **common,
        ),
        "val": DataLoader(
            BandgapDataset(splits["val"], target_mean, target_std),
            batch_size=batch_size,
            shuffle=False,
            **common,
        ),
        "test": DataLoader(
            BandgapDataset(splits["test"], target_mean, target_std),
            batch_size=batch_size,
            shuffle=False,
            **common,
        ),
    }


def _autocast_context(device: torch.device):
    enabled = bool(USE_AMP and device.type == "cuda")
    return torch.autocast(
        device_type=device.type, dtype=torch.bfloat16, enabled=enabled
    )


def train_one_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    criterion: nn.Module,
    device: torch.device,
    target_std: float,
) -> Tuple[float, float]:
    model.train()
    squared_error = 0.0
    normalized_loss = 0.0
    sample_count = 0
    for batch in loader:
        batch = batch.to(device, non_blocking=True)
        optimizer.zero_grad(set_to_none=True)
        try:
            with _autocast_context(device):
                prediction = model(batch).view(-1)
                target = batch.y.view(-1)
                loss = criterion(prediction, target)
            if not torch.isfinite(loss):
                raise FloatingPointError("Non-finite training loss detected.")
            loss.backward()
            torch.nn.utils.clip_grad_norm_(model.parameters(), max_norm=10.0)
            optimizer.step()
        except torch.cuda.OutOfMemoryError as exc:
            raise RuntimeError(
                "CUDA out of memory. Reduce BATCH_SIZE or HIDDEN_DIM."
            ) from exc
        count = target.numel()
        normalized_loss += float(loss.detach()) * count
        squared_error += (
            float(torch.sum((prediction.detach() - target) ** 2)) * target_std**2
        )
        sample_count += count
    return math.sqrt(squared_error / sample_count), normalized_loss / sample_count


@torch.inference_mode()
def evaluate_loader(
    model: nn.Module,
    loader: DataLoader,
    device: torch.device,
    target_mean: float,
    target_std: float,
) -> Dict[str, Any]:
    model.eval()
    names: List[str] = []
    true_values: List[float] = []
    predictions: List[float] = []
    for batch in loader:
        batch = batch.to(device, non_blocking=True)
        with _autocast_context(device):
            normalized_prediction = model(batch).view(-1)
        prediction = (
            normalized_prediction.float().cpu().numpy() * target_std + target_mean
        )
        target = batch.target_raw.view(-1).float().cpu().numpy()
        predictions.extend(prediction.tolist())
        true_values.extend(target.tolist())
        batch_names = (
            batch.cif_name if isinstance(batch.cif_name, list) else [batch.cif_name]
        )
        names.extend([str(name) for name in batch_names])
    y_true = np.asarray(true_values, dtype=float)
    y_pred = np.asarray(predictions, dtype=float)
    return {
        "names": names,
        "y_true": y_true,
        "y_pred": y_pred,
        "metrics": calculate_metrics(y_true, y_pred),
    }


def calculate_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    if len(y_true) == 0:
        raise ValueError("Cannot calculate metrics for an empty dataset.")
    return {
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(math.sqrt(mean_squared_error(y_true, y_pred))),
        "r2": float(r2_score(y_true, y_pred)) if len(y_true) >= 2 else float("nan"),
    }


def save_checkpoint(
    model: nn.Module,
    model_config: ModelConfig,
    graph_config: GraphConfig,
    target_mean: float,
    target_std: float,
    best_epoch: int,
    best_val_rmse: float,
) -> None:
    torch.save(
        {
            "format_version": 1,
            "model_name": MODEL_NAME,
            "run_version": RUN_VERSION,
            "model_state_dict": model.state_dict(),
            "model_config": asdict(model_config),
            "graph_config": asdict(graph_config),
            "best_epoch": int(best_epoch),
            "best_val_rmse": float(best_val_rmse),
            "seed": SEED,
            "target_name": TARGET_NAME,
            "target_unit": TARGET_UNIT,
            "normalization": {"mean": float(target_mean), "std": float(target_std)},
        },
        CHECKPOINT_PATH,
    )


def load_trained_model(
    checkpoint_path: os.PathLike[str] | str,
    device: Optional[torch.device] = None,
) -> Tuple[nn.Module, Dict[str, Any]]:
    selected_device = device or torch.device(
        "cuda" if torch.cuda.is_available() else "cpu"
    )
    try:
        checkpoint = torch_load_full(checkpoint_path, map_location=selected_device)
    except Exception as exc:
        raise RuntimeError(
            f"Could not load checkpoint {checkpoint_path}: {exc}"
        ) from exc
    model_config = ModelConfig(**checkpoint["model_config"])
    model = build_model(model_config).to(selected_device)
    model.load_state_dict(checkpoint["model_state_dict"])
    model.eval()
    return model, checkpoint


def train_model(
    loaders: Dict[str, DataLoader],
    model_config: ModelConfig,
    graph_config: GraphConfig,
    device: torch.device,
    target_mean: float,
    target_std: float,
    logger: logging.Logger,
    learning_rate: float,
    weight_decay: float,
    max_epochs: int,
    save_best: bool = True,
) -> Tuple[nn.Module, List[Dict[str, float]], int, float, bool]:
    model = build_model(model_config).to(device)
    criterion = nn.MSELoss()
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=learning_rate, weight_decay=weight_decay
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=max(5, PATIENCE // 5), min_lr=1e-7
    )
    logger.info("Model architecture:\n%s", model)
    logger.info(
        "Total parameters: %d",
        sum(parameter.numel() for parameter in model.parameters()),
    )
    logger.info(
        "Trainable parameters: %d",
        sum(
            parameter.numel()
            for parameter in model.parameters()
            if parameter.requires_grad
        ),
    )
    logger.info("Loss: MSELoss")
    logger.info("Optimizer: AdamW(lr=%g, weight_decay=%g)", learning_rate, weight_decay)
    logger.info("Scheduler: ReduceLROnPlateau(factor=0.5)")
    history: List[Dict[str, float]] = []
    best_val_rmse = float("inf")
    best_epoch = 0
    epochs_without_improvement = 0
    stopped_early = False
    best_state: Optional[Dict[str, Tensor]] = None
    for epoch in range(1, max_epochs + 1):
        train_rmse, train_mse_loss = train_one_epoch(
            model, loaders["train"], optimizer, criterion, device, target_std
        )
        val_result = evaluate_loader(
            model, loaders["val"], device, target_mean, target_std
        )
        val_rmse = val_result["metrics"]["rmse"]
        val_mse_loss = (val_rmse / target_std) ** 2
        scheduler.step(val_rmse)
        learning_rate_now = optimizer.param_groups[0]["lr"]
        history.append(
            {
                "epoch": epoch,
                "train_rmse_eV": train_rmse,
                "val_rmse_eV": val_rmse,
                "learning_rate": learning_rate_now,
                "train_mse_loss": train_mse_loss,
                "val_mse_loss": val_mse_loss,
            }
        )
        logger.info(
            "Epoch %04d | train RMSE %.6f eV | val RMSE %.6f eV | lr %.6e",
            epoch,
            train_rmse,
            val_rmse,
            learning_rate_now,
        )
        if val_rmse < best_val_rmse - MIN_DELTA:
            best_val_rmse = val_rmse
            best_epoch = epoch
            epochs_without_improvement = 0
            best_state = {
                key: value.detach().cpu().clone()
                for key, value in model.state_dict().items()
            }
            if save_best:
                save_checkpoint(
                    model,
                    model_config,
                    graph_config,
                    target_mean,
                    target_std,
                    best_epoch,
                    best_val_rmse,
                )
        else:
            epochs_without_improvement += 1
        if epochs_without_improvement >= PATIENCE:
            stopped_early = True
            logger.info("Early stopping triggered at epoch %d.", epoch)
            break
    if best_state is None:
        raise RuntimeError("Training did not produce a valid best model.")
    model.load_state_dict(best_state)
    if save_best and not CHECKPOINT_PATH.is_file():
        raise RuntimeError("Best checkpoint was not created.")
    return model, history, best_epoch, best_val_rmse, stopped_early


def run_optuna(
    splits: Dict[str, List[Dict[str, Any]]],
    graph_config: GraphConfig,
    device: torch.device,
    target_mean: float,
    target_std: float,
    logger: logging.Logger,
) -> Dict[str, Any]:
    try:
        import optuna
    except ImportError as exc:
        raise ImportError(
            "USE_OPTUNA=True requires Optuna: pip install optuna"
        ) from exc

    trial_logger = logging.getLogger("ComFormerBandgap.optuna")
    trial_logger.handlers = [logging.NullHandler()]
    trial_logger.propagate = False

    def objective(trial) -> float:
        trial_seed = SEED
        set_global_seed(trial_seed)
        hidden_dim = trial.suggest_categorical("hidden_dim", [128, 192, 256])
        batch_size = trial.suggest_categorical("batch_size", [16, 32, 64])
        model_config = ModelConfig(
            name=MODEL_NAME,
            hidden_dim=hidden_dim,
            fc_dim=trial.suggest_categorical("fc_dim", [128, 256]),
            dropout=trial.suggest_float("dropout", 0.0, 0.2),
        )
        loaders = make_loaders(splits, target_mean, target_std, batch_size)
        model, _, _, best_rmse, _ = train_model(
            loaders,
            model_config,
            graph_config,
            device,
            target_mean,
            target_std,
            trial_logger,
            learning_rate=trial.suggest_float("learning_rate", 1e-4, 3e-3, log=True),
            weight_decay=trial.suggest_float("weight_decay", 1e-7, 1e-3, log=True),
            max_epochs=min(MAX_EPOCHS, OPTUNA_MAX_EPOCHS),
            save_best=False,
        )
        del model, loaders
        if device.type == "cuda":
            torch.cuda.empty_cache()
        return best_rmse

    study = optuna.create_study(
        direction="minimize", sampler=optuna.samplers.TPESampler(seed=SEED)
    )
    study.optimize(objective, n_trials=OPTUNA_N_TRIALS, timeout=OPTUNA_TIMEOUT)
    result = {"best_value": float(study.best_value), "best_params": study.best_params}
    with open(
        TABLE_DIR / f"{MODEL_NAME}_optuna_best_params.json", "w", encoding="utf-8"
    ) as handle:
        json.dump(result, handle, indent=2)
    logger.info("Optuna best validation RMSE: %.6f eV", study.best_value)
    logger.info("Optuna best parameters: %s", study.best_params)
    return result


def save_prediction_dat(split_name: str, result: Dict[str, Any]) -> None:
    frame = pd.DataFrame(
        {
            "cif_name": result["names"],
            "true_bandgap_eV": result["y_true"],
            "predicted_bandgap_eV": result["y_pred"],
        }
    )
    frame["error_eV"] = frame["predicted_bandgap_eV"] - frame["true_bandgap_eV"]
    frame["absolute_error_eV"] = frame["error_eV"].abs()
    frame.to_csv(
        DAT_DIR / f"{MODEL_NAME}_parity_{split_name}.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )


def _plot_limits(results: Sequence[Dict[str, Any]]) -> Tuple[float, float]:
    values = np.concatenate(
        [np.concatenate([result["y_true"], result["y_pred"]]) for result in results]
    )
    minimum, maximum = float(np.min(values)), float(np.max(values))
    span = max(maximum - minimum, 1e-3)
    margin = 0.05 * span
    return minimum - margin, maximum + margin


def plot_parity(
    split_name: str, result: Dict[str, Any], limits: Tuple[float, float]
) -> None:
    metrics = result["metrics"]
    fig, axis = plt.subplots(figsize=(6.5, 6.0))
    color = {"train": COLORS[0], "val": COLORS[1], "test": COLORS[2]}[split_name]
    axis.scatter(
        result["y_true"],
        result["y_pred"],
        s=28,
        alpha=0.78,
        color=color,
        edgecolors="none",
    )
    axis.plot(limits, limits, linestyle="--", color="black", linewidth=1.2)
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_aspect("equal", adjustable="box")
    axis.set_title(
        f"{MODEL_NAME}: {split_name.capitalize()} parity", fontsize=TITLE_FONTSIZE
    )
    axis.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    axis.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    axis.tick_params(labelsize=TICK_FONTSIZE)
    axis.text(
        0.04,
        0.96,
        f"MAE = {metrics['mae']:.4f} eV\nRMSE = {metrics['rmse']:.4f} eV\n$R^2$ = {metrics['r2']:.4f}",
        transform=axis.transAxes,
        va="top",
        fontsize=ANNOTATION_FONTSIZE,
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.82),
    )
    if USE_GRID:
        axis.grid(True, linestyle="--", alpha=0.25)
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
    fig, axis = plt.subplots(figsize=(6.8, 6.2))
    for index, split_name in enumerate(("train", "val", "test")):
        result = results[split_name]
        axis.scatter(
            result["y_true"],
            result["y_pred"],
            s=27,
            alpha=0.72,
            color=COLORS[index],
            edgecolors="none",
            label=f"{split_name.capitalize()} (RMSE={result['metrics']['rmse']:.4f} eV)",
        )
    axis.plot(
        limits, limits, linestyle="--", color="black", linewidth=1.2, label="y = x"
    )
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_aspect("equal", adjustable="box")
    axis.set_title(f"{MODEL_NAME}: Combined parity", fontsize=TITLE_FONTSIZE)
    axis.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    axis.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    axis.tick_params(labelsize=TICK_FONTSIZE)
    axis.legend(fontsize=LEGEND_FONTSIZE)
    if USE_GRID:
        axis.grid(True, linestyle="--", alpha=0.25)
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
        float_format="%.8f",
    )
    fig, axis = plt.subplots(figsize=(7.4, 5.4))
    axis.plot(
        frame["epoch"], frame["train_rmse_eV"], color=COLORS[0], label="Train RMSE"
    )
    axis.plot(
        frame["epoch"], frame["val_rmse_eV"], color=COLORS[1], label="Validation RMSE"
    )
    best_row = frame.loc[frame["epoch"] == best_epoch].iloc[0]
    axis.scatter(
        [best_epoch],
        [best_row["val_rmse_eV"]],
        color=COLORS[4],
        marker="*",
        s=130,
        zorder=4,
        label=f"Best epoch: {best_epoch}",
    )
    axis.set_title(f"{MODEL_NAME}: RMSE learning curve", fontsize=TITLE_FONTSIZE)
    axis.set_xlabel("Epoch", fontsize=LABEL_FONTSIZE)
    axis.set_ylabel(f"RMSE ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    axis.tick_params(labelsize=TICK_FONTSIZE)
    axis.legend(fontsize=LEGEND_FONTSIZE)
    if USE_GRID:
        axis.grid(True, linestyle="--", alpha=0.25)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_rmse_curve.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(fig)


def save_results(
    results: Dict[str, Dict[str, Any]],
    history: Sequence[Dict[str, float]],
    best_epoch: int,
) -> None:
    metric_rows = []
    combined_rows = []
    limits = _plot_limits(list(results.values()))
    for split_name in ("train", "val", "test"):
        result = results[split_name]
        save_prediction_dat(split_name, result)
        plot_parity(split_name, result, limits)
        metrics = result["metrics"]
        metric_rows.append(
            {
                "split": split_name,
                "n_samples": len(result["y_true"]),
                "mae_eV": metrics["mae"],
                "rmse_eV": metrics["rmse"],
                "r2": metrics["r2"],
            }
        )
        for name, true_value, prediction in zip(
            result["names"], result["y_true"], result["y_pred"]
        ):
            combined_rows.append(
                {
                    "split": split_name,
                    "cif_name": name,
                    "true_bandgap_eV": true_value,
                    "predicted_bandgap_eV": prediction,
                    "error_eV": prediction - true_value,
                    "absolute_error_eV": abs(prediction - true_value),
                }
            )
    pd.DataFrame(metric_rows).to_csv(
        TABLE_DIR / f"{MODEL_NAME}_metrics.dat",
        sep="\t",
        index=False,
        float_format="%.6f",
    )
    pd.DataFrame(combined_rows).to_csv(
        DAT_DIR / f"{MODEL_NAME}_parity_all.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )
    plot_combined_parity(results, limits)
    plot_rmse_curve(history, best_epoch)


@torch.inference_mode()
def predict_cifs(
    checkpoint_path: os.PathLike[str] | str,
    cif_paths: Sequence[os.PathLike[str] | str],
    device: Optional[torch.device] = None,
    batch_size: int = BATCH_SIZE,
) -> pd.DataFrame:
    selected_device = device or torch.device(
        "cuda" if torch.cuda.is_available() else "cpu"
    )
    model, checkpoint = load_trained_model(checkpoint_path, selected_device)
    graph_config = GraphConfig(**checkpoint["graph_config"])
    graphs: List[Data] = []
    names: List[str] = []
    for cif_path in cif_paths:
        path = Path(cif_path).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"CIF file not found: {path}")
        graphs.append(
            build_graph_from_structure(Structure.from_file(path), graph_config)
        )
        names.append(path.name)
    loader = DataLoader(
        PredictionDataset(graphs, names), batch_size=batch_size, shuffle=False
    )
    mean = float(checkpoint["normalization"]["mean"])
    std = float(checkpoint["normalization"]["std"])
    output_names: List[str] = []
    predictions: List[float] = []
    for batch in loader:
        batch = batch.to(selected_device)
        normalized = model(batch).view(-1).float().cpu().numpy()
        predictions.extend((normalized * std + mean).tolist())
        batch_names = (
            batch.cif_name if isinstance(batch.cif_name, list) else [batch.cif_name]
        )
        output_names.extend([str(name) for name in batch_names])
    return pd.DataFrame(
        {
            "cif_name": output_names,
            f"predicted_{TARGET_NAME}_{TARGET_UNIT}": predictions,
        }
    )


def main() -> None:
    total_start = time.perf_counter()
    validate_configuration()
    create_directories()
    set_global_seed(SEED)
    logger = setup_logger()
    device = get_device(logger)
    logger.info("Model: %s", MODEL_NAME)
    logger.info("Run version: %s", RUN_VERSION)
    logger.info("Output root: %s", OUTPUT_ROOT)
    logger.info("Target: %s (%s)", TARGET_NAME, TARGET_UNIT)
    logger.info(
        "Atomic representation: learned %d-dimensional embedding", ATOM_FEATURES
    )
    logger.info("Seed: %d", SEED)
    logger.info(
        "Split ratios: train=%.3f, val=%.3f, test=%.3f",
        TRAIN_RATIO,
        VAL_RATIO,
        TEST_RATIO,
    )
    logger.info("Optuna enabled: %s", USE_OPTUNA)
    graph_config = GraphConfig()
    logger.info("Graph configuration: %s", asdict(graph_config))

    data_start = time.perf_counter()
    raw_records = load_excel_records(logger)
    valid_records = validate_and_build_records(raw_records, graph_config, logger)
    splits = split_records(valid_records, logger)
    data_seconds = time.perf_counter() - data_start
    logger.info(
        "Split sizes: train=%d, val=%d, test=%d",
        len(splits["train"]),
        len(splits["val"]),
        len(splits["test"]),
    )
    train_targets = np.asarray(
        [item["target"] for item in splits["train"]], dtype=float
    )
    target_mean = float(np.mean(train_targets))
    target_std = float(np.std(train_targets))
    if not np.isfinite(target_std) or target_std < 1e-12:
        raise ValueError("Training targets have zero or invalid standard deviation.")
    logger.info(
        "Training-target normalization: mean=%.8f, std=%.8f", target_mean, target_std
    )

    selected = {
        "hidden_dim": HIDDEN_DIM,
        "fc_dim": FC_DIM,
        "dropout": DROPOUT,
        "batch_size": BATCH_SIZE,
        "learning_rate": LEARNING_RATE,
        "weight_decay": WEIGHT_DECAY,
    }
    if USE_OPTUNA:
        optuna_result = run_optuna(
            splits, graph_config, device, target_mean, target_std, logger
        )
        selected.update(optuna_result["best_params"])
        set_global_seed(SEED)
    model_config = ModelConfig(
        name=MODEL_NAME,
        hidden_dim=int(selected["hidden_dim"]),
        fc_dim=int(selected["fc_dim"]),
        dropout=float(selected["dropout"]),
    )
    loaders = make_loaders(splits, target_mean, target_std, int(selected["batch_size"]))
    logger.info("Model configuration: %s", asdict(model_config))

    training_start = time.perf_counter()
    _, history, best_epoch, best_val_rmse, stopped_early = train_model(
        loaders,
        model_config,
        graph_config,
        device,
        target_mean,
        target_std,
        logger,
        learning_rate=float(selected["learning_rate"]),
        weight_decay=float(selected["weight_decay"]),
        max_epochs=MAX_EPOCHS,
        save_best=True,
    )
    training_seconds = time.perf_counter() - training_start
    logger.info("Best epoch: %d", best_epoch)
    logger.info("Best validation RMSE: %.6f eV", best_val_rmse)
    logger.info("Early stopping triggered: %s", stopped_early)

    evaluation_start = time.perf_counter()
    best_model, checkpoint = load_trained_model(CHECKPOINT_PATH, device)
    results = {
        split_name: evaluate_loader(
            best_model, loaders[split_name], device, target_mean, target_std
        )
        for split_name in ("train", "val", "test")
    }
    save_results(results, history, int(checkpoint["best_epoch"]))
    evaluation_seconds = time.perf_counter() - evaluation_start
    for split_name in ("train", "val", "test"):
        metrics = results[split_name]["metrics"]
        logger.info(
            "%s metrics | MAE %.6f eV | RMSE %.6f eV | R2 %.6f",
            split_name.capitalize(),
            metrics["mae"],
            metrics["rmse"],
            metrics["r2"],
        )
    total_seconds = time.perf_counter() - total_start
    logger.info("Data time: %.3f s (%s)", data_seconds, elapsed_string(data_seconds))
    logger.info(
        "Training time: %.3f s (%s)", training_seconds, elapsed_string(training_seconds)
    )
    logger.info(
        "Evaluation/plotting time: %.3f s (%s)",
        evaluation_seconds,
        elapsed_string(evaluation_seconds),
    )
    logger.info(
        "Total runtime: %.3f s (%s)", total_seconds, elapsed_string(total_seconds)
    )
    logger.info("Best checkpoint: %s", CHECKPOINT_PATH)


if __name__ == "__main__":
    main()

import hashlib
import json
import logging
import math
import os
import random
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple, Union

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
from torch_geometric.data import Batch, Data
from torch_geometric.loader import DataLoader
from torch_geometric.nn import MessagePassing, global_mean_pool
from tqdm import tqdm

try:
    from torch_scatter import scatter_add

    SCATTER_BACKEND = "torch_scatter"
except ImportError:
    SCATTER_BACKEND = "torch.index_add_"

    def scatter_add(
        source: torch.Tensor,
        index: torch.Tensor,
        dim: int = 0,
        dim_size: Optional[int] = None,
    ) -> torch.Tensor:
        if index.ndim != 1:
            raise ValueError(
                "The PyTorch scatter fallback requires a one-dimensional index tensor."
            )
        if source.shape[dim] != index.numel():
            raise ValueError(
                "The scatter index length must match the source size along the scatter dimension."
            )
        if dim_size is None:
            dim_size = int(index.max().item()) + 1 if index.numel() else 0
        output_shape = list(source.shape)
        output_shape[dim] = int(dim_size)
        output = source.new_zeros(output_shape)
        return output.index_add_(dim, index, source)


# -------------------- User configuration --------------------

MODEL_NAME = "ReciNet"
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
NUM_WORKERS = 4
PIN_MEMORY = True
GRADIENT_CLIP_NORM = 5.0

CUTOFF = 4.0
MAX_NEIGHBORS = 16
MAX_ATOMIC_NUMBER = 92
NUM_K_VECTORS = 16
K_SEARCH_RADIUS = 2

CONV_LAYERS = 4
HIDDEN_DIM = 256
DOWNPROJECTION_DIM = 64
RECIPROCAL_HIDDEN_LAYERS = 2
DROPOUT = 0.0
RBF_MIN = -4.0
RBF_MAX = 4.0

USE_OPTUNA = False
OPTUNA_N_TRIALS = 30
OPTUNA_TIMEOUT = None
OPTUNA_MAX_EPOCHS = 150

FIG_DPI = 600
TITLE_FONTSIZE = 14
LABEL_FONTSIZE = 12
TICK_FONTSIZE = 10
LEGEND_FONTSIZE = 10
ANNOTATION_FONTSIZE = 10
USE_GRID = True
COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:red"]


SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = SCRIPT_DIR / f"{MODEL_NAME}_{RUN_VERSION}"
FIGURE_DIR = OUTPUT_DIR / "figure"
DAT_DIR = OUTPUT_DIR / "dat"
TABLE_DIR = OUTPUT_DIR / "table"
LOG_DIR = OUTPUT_DIR / "log"
SPLIT_DIR = OUTPUT_DIR / "split"
CACHE_DIR = OUTPUT_DIR / "cache"
CHECKPOINT_PATH = OUTPUT_DIR / f"{MODEL_NAME}_best.pt"
SPLIT_PATH = SPLIT_DIR / f"{MODEL_NAME}_split.csv"
LOG_PATH = LOG_DIR / f"{MODEL_NAME}_training.log"


@dataclass
class ModelConfig:
    atom_input_features: int = MAX_ATOMIC_NUMBER
    hidden_dim: int = HIDDEN_DIM
    conv_layers: int = CONV_LAYERS
    downprojection_dim: int = DOWNPROJECTION_DIM
    reciprocal_hidden_layers: int = RECIPROCAL_HIDDEN_LAYERS
    dropout: float = DROPOUT
    rbf_min: float = RBF_MIN
    rbf_max: float = RBF_MAX


@dataclass
class GraphConfig:
    cutoff: float = CUTOFF
    max_neighbors: int = MAX_NEIGHBORS
    max_atomic_number: int = MAX_ATOMIC_NUMBER
    num_k_vectors: int = NUM_K_VECTORS
    k_search_radius: int = K_SEARCH_RADIUS


def create_output_directories() -> None:
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
    logger = logging.getLogger(f"{MODEL_NAME}_{RUN_VERSION}")
    logger.setLevel(logging.INFO)
    logger.propagate = False
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
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    torch.use_deterministic_algorithms(True, warn_only=True)


def seed_worker(worker_id: int) -> None:
    worker_seed = (torch.initial_seed() + worker_id) % (2**32)
    np.random.seed(worker_seed)
    random.seed(worker_seed)


def format_duration(seconds: float) -> str:
    seconds_int = max(0, int(round(seconds)))
    hours, remainder = divmod(seconds_int, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{seconds:.2f} s ({hours:02d}:{minutes:02d}:{secs:02d})"


def resolve_path(path_value: Union[str, Path]) -> Path:
    path = Path(path_value).expanduser()
    return path if path.is_absolute() else SCRIPT_DIR / path


def safe_torch_load(path: Union[str, Path], map_location=None):
    try:
        return torch.load(path, map_location=map_location, weights_only=False)
    except TypeError:
        return torch.load(path, map_location=map_location)


def device_information(device: torch.device) -> Dict[str, object]:
    info: Dict[str, object] = {
        "device": str(device),
        "cuda_available": torch.cuda.is_available(),
    }
    if device.type == "cuda":
        index = (
            device.index if device.index is not None else torch.cuda.current_device()
        )
        properties = torch.cuda.get_device_properties(index)
        info.update(
            {
                "gpu_name": properties.name,
                "gpu_memory_gb": properties.total_memory / (1024**3),
                "cuda_version": torch.version.cuda,
            }
        )
    return info


def read_excel_records(logger: logging.Logger) -> List[Dict[str, object]]:
    excel_path = resolve_path(EXCEL_PATH)
    cif_dir = resolve_path(CIF_DIR)
    if not excel_path.is_file():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    if not cif_dir.is_dir():
        raise FileNotFoundError(f"CIF directory not found: {cif_dir}")

    frame = pd.read_excel(excel_path)
    if frame.shape[1] < 2:
        raise ValueError("The Excel file must contain at least two columns.")

    records: List[Dict[str, object]] = []
    for row_number, (cif_value, target_value) in enumerate(
        zip(frame.iloc[:, 0], frame.iloc[:, 1]), start=2
    ):
        if pd.isna(cif_value) or str(cif_value).strip() == "":
            logger.warning("Skipped Excel row %d: empty CIF filename.", row_number)
            continue
        cif_name = f"{int(cif_value)}.cif"
        try:
            target = float(target_value)
        except (TypeError, ValueError):
            logger.warning("Skipped %s: target is not numeric.", cif_name)
            continue
        if not np.isfinite(target):
            logger.warning("Skipped %s: target is not finite.", cif_name)
            continue
        cif_path = cif_dir / cif_name
        if not cif_path.is_file():
            logger.warning("Skipped %s: CIF file not found at %s.", cif_name, cif_path)
            continue
        records.append({"cif_name": cif_name, "cif_path": cif_path, "target": target})

    if len(records) < 10:
        raise ValueError(
            f"Only {len(records)} valid Excel rows remain; at least 10 are required."
        )
    return records


def atomic_features(structure: Structure, max_atomic_number: int) -> torch.Tensor:
    features = np.zeros((len(structure), max_atomic_number), dtype=np.float32)
    for site_index, site in enumerate(structure):
        total_occupancy = 0.0
        for specie, occupancy in site.species.items():
            atomic_number = int(specie.Z)
            if atomic_number < 1 or atomic_number > max_atomic_number:
                raise ValueError(
                    f"Atomic number {atomic_number} is outside the supported range 1-{max_atomic_number}."
                )
            features[site_index, atomic_number - 1] += float(occupancy)
            total_occupancy += float(occupancy)
        if total_occupancy <= 0.0:
            raise ValueError(f"Site {site_index} has no positive occupancy.")
        features[site_index] /= total_occupancy
    return torch.from_numpy(features)


def build_periodic_edges(
    structure: Structure, cutoff: float, max_neighbors: int
) -> Tuple[torch.Tensor, torch.Tensor]:
    center, neighbor, offsets, distances = structure.get_neighbor_list(
        r=cutoff, numerical_tol=1.0e-8
    )
    candidates: List[List[Tuple[float, int, np.ndarray]]] = [
        [] for _ in range(len(structure))
    ]
    for center_index, neighbor_index, offset, distance in zip(
        center, neighbor, offsets, distances
    ):
        if float(distance) <= 1.0e-8:
            continue
        candidates[int(center_index)].append(
            (float(distance), int(neighbor_index), np.asarray(offset))
        )

    sources: List[int] = []
    targets: List[int] = []
    selected_distances: List[float] = []
    for target_index, atom_neighbors in enumerate(candidates):
        atom_neighbors.sort(key=lambda item: item[0])
        for distance, source_index, _ in atom_neighbors[:max_neighbors]:
            sources.append(source_index)
            targets.append(target_index)
            selected_distances.append(distance)

    if not sources:
        raise ValueError(
            f"No periodic neighbors were found within cutoff={cutoff:.3f} angstrom."
        )
    edge_index = torch.tensor([sources, targets], dtype=torch.long)
    edge_distance = torch.tensor(selected_distances, dtype=torch.float32)
    return edge_index, edge_distance


def generate_reciprocal_vectors(
    structure: Structure, num_k_vectors: int, search_radius: int
) -> torch.Tensor:
    if num_k_vectors < 1:
        raise ValueError("num_k_vectors must be positive.")
    if search_radius < 1:
        raise ValueError("k_search_radius must be positive.")

    reciprocal_basis = np.asarray(
        structure.lattice.reciprocal_lattice.matrix, dtype=np.float64
    )
    indices = []
    for h in range(-search_radius, search_radius + 1):
        for k in range(-search_radius, search_radius + 1):
            for l in range(-search_radius, search_radius + 1):
                if h == 0 and k == 0 and l == 0:
                    continue
                hkl = np.asarray([h, k, l], dtype=np.float64)
                vector = hkl @ reciprocal_basis
                indices.append((float(np.linalg.norm(vector)), h, k, l, vector))
    indices.sort(key=lambda item: (item[0], item[1], item[2], item[3]))
    if len(indices) < num_k_vectors:
        raise ValueError(
            f"The reciprocal search grid contains only {len(indices)} nonzero vectors, "
            f"but {num_k_vectors} were requested."
        )
    vectors = np.stack([item[4] for item in indices[:num_k_vectors]], axis=0)
    if not np.isfinite(vectors).all():
        raise ValueError("The reciprocal lattice contains non-finite values.")
    return torch.tensor(vectors, dtype=torch.float32).unsqueeze(0)


def graph_cache_path(cif_path: Path, graph_config: GraphConfig) -> Path:
    stat = cif_path.stat()
    token = {
        "path": str(cif_path.resolve()),
        "size": stat.st_size,
        "mtime_ns": stat.st_mtime_ns,
        "graph_config": asdict(graph_config),
        "cache_format": 2,
    }
    digest = hashlib.sha256(
        json.dumps(token, sort_keys=True).encode("utf-8")
    ).hexdigest()[:24]
    return CACHE_DIR / f"{digest}.pt"


def build_graph(
    cif_path: Union[str, Path], cif_name: str, graph_config: GraphConfig
) -> Data:
    structure = Structure.from_file(str(cif_path))
    if len(structure) == 0:
        raise ValueError("The structure contains no atomic sites.")
    if abs(float(structure.lattice.volume)) <= 1.0e-10:
        raise ValueError("The lattice volume is zero or numerically singular.")

    x = atomic_features(structure, graph_config.max_atomic_number)
    edge_index, edge_distance = build_periodic_edges(
        structure, graph_config.cutoff, graph_config.max_neighbors
    )
    pos = torch.tensor(np.asarray(structure.cart_coords), dtype=torch.float32)
    k_vectors = generate_reciprocal_vectors(
        structure, graph_config.num_k_vectors, graph_config.k_search_radius
    )
    lattice = torch.tensor(
        np.asarray(structure.lattice.matrix), dtype=torch.float32
    ).unsqueeze(0)
    return Data(
        x=x,
        edge_index=edge_index,
        edge_attr=edge_distance,
        pos=pos,
        k_vectors=k_vectors,
        lattice=lattice,
        cif_name=cif_name,
    )


def build_or_load_graph(record: Dict[str, object], graph_config: GraphConfig) -> Data:
    cif_path = Path(record["cif_path"])
    cache_path = graph_cache_path(cif_path, graph_config)
    if cache_path.is_file():
        graph = safe_torch_load(cache_path, map_location="cpu")
    else:
        graph = build_graph(cif_path, str(record["cif_name"]), graph_config)
        torch.save(graph, cache_path)
    graph.raw_y = torch.tensor([float(record["target"])], dtype=torch.float32)
    graph.cif_name = str(record["cif_name"])
    return graph


def validate_and_cache_records(
    raw_records: Sequence[Dict[str, object]],
    graph_config: GraphConfig,
    logger: logging.Logger,
) -> List[Dict[str, object]]:
    valid_records: List[Dict[str, object]] = []
    for record in tqdm(raw_records, desc="Building crystal graphs"):
        try:
            graph = build_or_load_graph(record, graph_config)
            valid_records.append({**record, "graph": graph})
        except Exception as exc:
            logger.warning("Skipped %s: %s", record["cif_name"], exc)
    if len(valid_records) < 10:
        raise ValueError(
            f"Only {len(valid_records)} samples produced valid graphs; at least 10 are required."
        )
    return valid_records


def create_or_load_split(
    records: Sequence[Dict[str, object]], logger: logging.Logger
) -> Dict[str, List[int]]:
    if not math.isclose(
        TRAIN_RATIO + VAL_RATIO + TEST_RATIO, 1.0, rel_tol=0.0, abs_tol=1.0e-9
    ):
        raise ValueError("TRAIN_RATIO + VAL_RATIO + TEST_RATIO must equal 1.0.")

    name_to_index = {
        str(record["cif_name"]): index for index, record in enumerate(records)
    }
    if len(name_to_index) != len(records):
        raise ValueError("CIF filenames must be unique after validation.")

    if SPLIT_PATH.is_file():
        split_frame = pd.read_csv(SPLIT_PATH)
        required = {"cif_name", "target", "split"}
        if required.issubset(split_frame.columns):
            saved_names = split_frame["cif_name"].astype(str).tolist()
            current_names = list(name_to_index.keys())
            valid_labels = set(split_frame["split"].astype(str)) == {
                "train",
                "val",
                "test",
            }
            if (
                len(saved_names) == len(set(saved_names))
                and set(saved_names) == set(current_names)
                and valid_labels
            ):
                result = {
                    split: [
                        name_to_index[name]
                        for name in split_frame.loc[
                            split_frame["split"] == split, "cif_name"
                        ].astype(str)
                    ]
                    for split in ("train", "val", "test")
                }
                if all(result.values()):
                    logger.info("Reused split file: %s", SPLIT_PATH)
                    return result
        logger.warning("Existing split file is incompatible and will be replaced.")

    all_indices = np.arange(len(records))
    train_indices, remainder = train_test_split(
        all_indices,
        test_size=VAL_RATIO + TEST_RATIO,
        random_state=SEED,
        shuffle=True,
    )
    relative_test_ratio = TEST_RATIO / (VAL_RATIO + TEST_RATIO)
    val_indices, test_indices = train_test_split(
        remainder,
        test_size=relative_test_ratio,
        random_state=SEED,
        shuffle=True,
    )
    split_indices = {
        "train": train_indices.tolist(),
        "val": val_indices.tolist(),
        "test": test_indices.tolist(),
    }
    rows = []
    for split in ("train", "val", "test"):
        for index in split_indices[split]:
            rows.append(
                {
                    "cif_name": records[index]["cif_name"],
                    "target": float(records[index]["target"]),
                    "split": split,
                }
            )
    pd.DataFrame(rows).to_csv(SPLIT_PATH, index=False)
    logger.info("Created split file: %s", SPLIT_PATH)
    return split_indices


class CrystalDataset(torch.utils.data.Dataset):
    def __init__(
        self,
        records: Sequence[Dict[str, object]],
        indices: Sequence[int],
        target_mean: float,
        target_std: float,
    ):
        self.records = records
        self.indices = list(indices)
        self.target_mean = float(target_mean)
        self.target_std = float(target_std)

    def __len__(self) -> int:
        return len(self.indices)

    def __getitem__(self, item: int) -> Data:
        record = self.records[self.indices[item]]
        graph = record["graph"].clone()
        raw_target = float(record["target"])
        graph.y = torch.tensor(
            [(raw_target - self.target_mean) / self.target_std], dtype=torch.float32
        )
        graph.raw_y = torch.tensor([raw_target], dtype=torch.float32)
        graph.cif_name = str(record["cif_name"])
        return graph


def make_data_loaders(
    records: Sequence[Dict[str, object]],
    split_indices: Dict[str, List[int]],
    batch_size: int,
) -> Tuple[Dict[str, DataLoader], Dict[str, float]]:
    train_targets = np.asarray(
        [float(records[index]["target"]) for index in split_indices["train"]]
    )
    target_mean = float(train_targets.mean())
    target_std = float(train_targets.std(ddof=0))
    if not np.isfinite(target_std) or target_std < 1.0e-12:
        target_std = 1.0
    normalization = {"mean": target_mean, "std": target_std}

    datasets = {
        split: CrystalDataset(records, indices, target_mean, target_std)
        for split, indices in split_indices.items()
    }
    generator = torch.Generator()
    generator.manual_seed(SEED)
    common = {
        "batch_size": batch_size,
        "num_workers": NUM_WORKERS,
        "pin_memory": PIN_MEMORY and torch.cuda.is_available(),
        "worker_init_fn": seed_worker,
        "persistent_workers": NUM_WORKERS > 0,
    }
    loaders = {
        "train": DataLoader(
            datasets["train"], shuffle=True, generator=generator, **common
        ),
        "val": DataLoader(datasets["val"], shuffle=False, **common),
        "test": DataLoader(datasets["test"], shuffle=False, **common),
    }
    return loaders, normalization


class RBFExpansion(nn.Module):
    def __init__(self, vmin: float, vmax: float, bins: int):
        super().__init__()
        centers = torch.linspace(vmin, vmax, bins)
        self.register_buffer("centers", centers)
        spacing = float(centers[1] - centers[0]) if bins > 1 else 1.0
        self.gamma = 1.0 / max(abs(spacing), 1.0e-12)

    def forward(self, values: torch.Tensor) -> torch.Tensor:
        scaled = self.gamma * (values.unsqueeze(-1) - self.centers)
        return torch.exp(-(scaled**2))


class ScaledSiLU(nn.Module):
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return F.silu(x) / 0.6


class Dense(nn.Module):
    def __init__(
        self,
        in_features: int,
        out_features: int,
        bias: bool = False,
        activation: Optional[str] = None,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.linear = nn.Linear(in_features, out_features, bias=bias)
        self.activation = (
            ScaledSiLU() if activation in {"silu", "swish"} else nn.Identity()
        )
        self.dropout = float(dropout)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x = self.activation(self.linear(x))
        return (
            F.dropout(x, p=self.dropout, training=self.training)
            if self.dropout > 0.0
            else x
        )


class ResidualLayer(nn.Module):
    def __init__(self, units: int, activation: str = "silu", dropout: float = 0.0):
        super().__init__()
        self.layers = nn.Sequential(
            Dense(units, units, activation=activation, dropout=dropout),
            Dense(units, units, activation=activation, dropout=dropout),
        )
        self.scale = 1.0 / math.sqrt(2.0)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return (x + self.layers(x)) * self.scale


class SafeBatchNorm1d(nn.BatchNorm1d):
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        if self.training and x.shape[0] <= 1:
            return F.batch_norm(
                x,
                self.running_mean,
                self.running_var,
                self.weight,
                self.bias,
                False,
                self.momentum,
                self.eps,
            )
        return super().forward(x)


class LocalConv(MessagePassing):
    def __init__(self, hidden_dim: int):
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
        update = self.propagate(
            edge_index, x=x, edge_attr=edge_attr, size=(x.size(0), x.size(0))
        )
        return F.relu(x + self.output_norm(update))

    def message(
        self, x_i: torch.Tensor, x_j: torch.Tensor, edge_attr: torch.Tensor
    ) -> torch.Tensor:
        inputs = torch.cat((x_i, x_j, edge_attr), dim=-1)
        gate = torch.sigmoid(self.gate_norm(self.gate_network(inputs)))
        return gate * self.message_network(inputs)


class ReciprocalBlock(nn.Module):
    def __init__(
        self,
        hidden_dim: int,
        downprojection_dim: int,
        num_hidden: int,
        dropout: float,
    ):
        super().__init__()
        self.pre_residual = ResidualLayer(
            hidden_dim, activation="silu", dropout=dropout
        )
        self.down = Dense(
            hidden_dim,
            downprojection_dim,
            bias=False,
            activation="silu",
            dropout=dropout,
        )
        self.up = Dense(downprojection_dim, hidden_dim, bias=False)
        layers: List[nn.Module] = [
            Dense(
                hidden_dim, hidden_dim, bias=False, activation="silu", dropout=dropout
            )
        ]
        layers.extend(
            ResidualLayer(hidden_dim, activation="silu", dropout=dropout)
            for _ in range(num_hidden)
        )
        self.reciprocal_layers = nn.Sequential(*layers)
        self.h_update_scale = nn.Parameter(torch.tensor(0.01, dtype=torch.float32))

    def forward(
        self,
        h: torch.Tensor,
        positions: torch.Tensor,
        k_vectors: torch.Tensor,
        batch: torch.Tensor,
        num_graphs: int,
    ) -> torch.Tensor:
        if k_vectors.ndim != 3 or k_vectors.shape[-1] != 3:
            raise ValueError(
                f"Expected k_vectors with shape [num_graphs, num_k, 3], got {tuple(k_vectors.shape)}."
            )
        if k_vectors.shape[0] != num_graphs:
            raise ValueError(
                f"k-vector batch dimension {k_vectors.shape[0]} does not match {num_graphs} graphs."
            )

        h_residual = self.pre_residual(h)
        projected = self.down(h_residual)
        atom_k = k_vectors.index_select(0, batch)
        phase = torch.einsum("nc,nkc->nk", positions, atom_k)
        cosine = torch.cos(phase)
        sine = torch.sin(phase)

        real_terms = projected.unsqueeze(1) * cosine.unsqueeze(-1)
        imag_terms = projected.unsqueeze(1) * sine.unsqueeze(-1)
        structure_real = scatter_add(real_terms, batch, dim=0, dim_size=num_graphs)
        structure_imag = scatter_add(imag_terms, batch, dim=0, dim_size=num_graphs)

        counts = scatter_add(
            projected.new_ones((projected.shape[0], 1)),
            batch,
            dim=0,
            dim_size=num_graphs,
        ).clamp_min_(1.0)
        structure_real = structure_real / counts.unsqueeze(1)
        structure_imag = structure_imag / counts.unsqueeze(1)

        atom_real = structure_real.index_select(0, batch)
        atom_imag = structure_imag.index_select(0, batch)
        reciprocal_update = (
            atom_real * cosine.unsqueeze(-1) + atom_imag * sine.unsqueeze(-1)
        ).mean(dim=1)
        reciprocal_update = self.up(reciprocal_update)
        reciprocal_update = self.reciprocal_layers(reciprocal_update)
        return self.h_update_scale * reciprocal_update


class ShiftedSoftplus(nn.Module):
    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return F.softplus(x) - math.log(2.0)


class ReciNet(nn.Module):
    def __init__(self, config: ModelConfig):
        super().__init__()
        self.config = config
        hidden_dim = config.hidden_dim
        self.atom_embedding = nn.Linear(config.atom_input_features, hidden_dim)
        self.edge_embedding = nn.Sequential(
            RBFExpansion(config.rbf_min, config.rbf_max, hidden_dim),
            nn.Linear(hidden_dim, hidden_dim),
            nn.SiLU(),
        )
        self.local_modules = nn.ModuleList(
            LocalConv(hidden_dim) for _ in range(config.conv_layers)
        )
        self.reciprocal_blocks = nn.ModuleList(
            ReciprocalBlock(
                hidden_dim,
                config.downprojection_dim,
                config.reciprocal_hidden_layers,
                config.dropout,
            )
            for _ in range(config.conv_layers)
        )
        self.initial_reciprocal_projection = Dense(
            hidden_dim, hidden_dim, bias=True, activation="swish"
        )
        self.readout = nn.Sequential(
            nn.Linear(hidden_dim, hidden_dim), ShiftedSoftplus()
        )
        self.output = nn.Linear(hidden_dim, 1)

    def forward(self, data: Union[Data, Batch]) -> torch.Tensor:
        edge_distance = data.edge_attr.view(-1).clamp_min(1.0e-8)
        edge_features = self.edge_embedding(-0.75 / edge_distance)
        node_features = self.atom_embedding(data.x)
        reciprocal_features = self.initial_reciprocal_projection(node_features)
        batch = getattr(
            data,
            "batch",
            torch.zeros(
                node_features.shape[0], dtype=torch.long, device=node_features.device
            ),
        )
        num_graphs = int(getattr(data, "num_graphs", 1))

        for local_module, reciprocal_block in zip(
            self.local_modules, self.reciprocal_blocks
        ):
            reciprocal_update = reciprocal_block(
                reciprocal_features,
                data.pos,
                data.k_vectors,
                batch,
                num_graphs,
            )
            reciprocal_features = reciprocal_features + reciprocal_update
            local_features = local_module(node_features, data.edge_index, edge_features)
            node_features = reciprocal_features + local_features
            reciprocal_features = node_features

        crystal_features = global_mean_pool(node_features, batch, size=num_graphs)
        return self.output(self.readout(crystal_features)).view(-1)


def count_parameters(model: nn.Module) -> Tuple[int, int]:
    total = sum(parameter.numel() for parameter in model.parameters())
    trainable = sum(
        parameter.numel() for parameter in model.parameters() if parameter.requires_grad
    )
    return total, trainable


def optimizer_parameter_groups(model: nn.Module, weight_decay: float):
    decay, no_decay = [], []
    for name, parameter in model.named_parameters():
        if not parameter.requires_grad:
            continue
        if parameter.ndim == 1 or name.endswith("bias") or "norm" in name.lower():
            no_decay.append(parameter)
        else:
            decay.append(parameter)
    return [
        {"params": decay, "weight_decay": weight_decay},
        {"params": no_decay, "weight_decay": 0.0},
    ]


def calculate_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    if len(y_true) == 0:
        raise ValueError("Cannot calculate metrics for an empty array.")
    mae = float(mean_absolute_error(y_true, y_pred))
    rmse = float(np.sqrt(mean_squared_error(y_true, y_pred)))
    r2 = float(r2_score(y_true, y_pred)) if len(y_true) >= 2 else float("nan")
    return {"mae": mae, "rmse": rmse, "r2": r2}


def train_one_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    criterion: nn.Module,
    device: torch.device,
    normalization: Dict[str, float],
) -> Dict[str, float]:
    model.train()
    total_loss = 0.0
    total_samples = 0
    all_true: List[np.ndarray] = []
    all_pred: List[np.ndarray] = []

    for batch in loader:
        batch = batch.to(device, non_blocking=True)
        optimizer.zero_grad(set_to_none=True)
        with torch.autocast(
            device_type=device.type,
            dtype=torch.bfloat16 if device.type == "cuda" else torch.float32,
            enabled=device.type == "cuda",
        ):
            prediction = model(batch).view(-1)
            target = batch.y.view(-1)
            loss = criterion(prediction, target)
        if not torch.isfinite(loss):
            raise FloatingPointError(
                f"Non-finite training loss detected: {loss.item()}"
            )
        loss.backward()
        if GRADIENT_CLIP_NORM is not None and GRADIENT_CLIP_NORM > 0:
            torch.nn.utils.clip_grad_norm_(model.parameters(), GRADIENT_CLIP_NORM)
        optimizer.step()

        count = int(target.numel())
        total_loss += float(loss.detach()) * count
        total_samples += count
        predicted_raw = (
            prediction.detach().float().cpu().numpy() * normalization["std"]
            + normalization["mean"]
        )
        all_pred.append(predicted_raw)
        all_true.append(batch.raw_y.view(-1).detach().float().cpu().numpy())

    y_true = np.concatenate(all_true)
    y_pred = np.concatenate(all_pred)
    return {
        "mse_loss": total_loss / max(total_samples, 1),
        "rmse": float(np.sqrt(mean_squared_error(y_true, y_pred))),
    }


@torch.inference_mode()
def evaluate(
    model: nn.Module,
    loader: DataLoader,
    criterion: nn.Module,
    device: torch.device,
    normalization: Dict[str, float],
    include_names: bool = False,
) -> Dict[str, object]:
    model.eval()
    total_loss = 0.0
    total_samples = 0
    all_true: List[np.ndarray] = []
    all_pred: List[np.ndarray] = []
    names: List[str] = []

    for batch in loader:
        batch = batch.to(device, non_blocking=True)
        prediction = model(batch).view(-1)
        target = batch.y.view(-1)
        loss = criterion(prediction, target)
        if not torch.isfinite(loss):
            raise FloatingPointError(
                f"Non-finite evaluation loss detected: {loss.item()}"
            )
        count = int(target.numel())
        total_loss += float(loss) * count
        total_samples += count
        all_pred.append(
            prediction.float().cpu().numpy() * normalization["std"]
            + normalization["mean"]
        )
        all_true.append(batch.raw_y.view(-1).float().cpu().numpy())
        if include_names:
            batch_names = batch.cif_name
            names.extend(
                [str(value) for value in batch_names]
                if isinstance(batch_names, (list, tuple))
                else [str(batch_names)]
            )

    y_true = np.concatenate(all_true)
    y_pred = np.concatenate(all_pred)
    metrics = calculate_metrics(y_true, y_pred)
    return {
        "mse_loss": total_loss / max(total_samples, 1),
        "y_true": y_true,
        "y_pred": y_pred,
        "cif_name": names,
        **metrics,
    }


def save_checkpoint(
    path: Path,
    model: nn.Module,
    model_config: ModelConfig,
    graph_config: GraphConfig,
    normalization: Dict[str, float],
    epoch: int,
    val_rmse: float,
    optimizer: torch.optim.Optimizer,
    hyperparameters: Dict[str, object],
) -> None:
    checkpoint = {
        "model_name": MODEL_NAME,
        "run_version": RUN_VERSION,
        "model_state_dict": model.state_dict(),
        "model_config": asdict(model_config),
        "graph_config": asdict(graph_config),
        "normalization": normalization,
        "best_epoch": int(epoch),
        "best_val_rmse": float(val_rmse),
        "seed": SEED,
        "target_name": TARGET_NAME,
        "target_unit": TARGET_UNIT,
        "hyperparameters": hyperparameters,
        "optimizer_state_dict": optimizer.state_dict(),
        "format_version": 2,
    }
    torch.save(checkpoint, path)


def load_trained_model(
    checkpoint_path: Union[str, Path], device: Optional[torch.device] = None
) -> Tuple[ReciNet, Dict[str, object]]:
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint = safe_torch_load(checkpoint_path, map_location=device)
    required = {"model_state_dict", "model_config", "graph_config", "normalization"}
    missing = required.difference(checkpoint)
    if missing:
        raise KeyError(f"Checkpoint is missing required fields: {sorted(missing)}")
    model = ReciNet(ModelConfig(**checkpoint["model_config"]))
    model.load_state_dict(checkpoint["model_state_dict"], strict=True)
    model.to(device)
    model.eval()
    return model, checkpoint


def fit_model(
    loaders: Dict[str, DataLoader],
    normalization: Dict[str, float],
    model_config: ModelConfig,
    graph_config: GraphConfig,
    device: torch.device,
    learning_rate: float,
    weight_decay: float,
    max_epochs: int,
    patience: int,
    logger: logging.Logger,
    write_checkpoint: bool,
    trial=None,
) -> Tuple[float, int, List[Dict[str, float]], bool]:
    model = ReciNet(model_config).to(device)
    criterion = nn.MSELoss()
    optimizer = torch.optim.AdamW(
        optimizer_parameter_groups(model, weight_decay), lr=learning_rate
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=max(5, patience // 5), min_lr=1.0e-6
    )
    hyperparameters = {
        "learning_rate": learning_rate,
        "weight_decay": weight_decay,
        "batch_size": loaders["train"].batch_size,
        "max_epochs": max_epochs,
        "patience": patience,
        "min_delta": MIN_DELTA,
    }

    best_rmse = float("inf")
    best_epoch = 0
    epochs_without_improvement = 0
    history: List[Dict[str, float]] = []
    stopped_early = False

    for epoch in range(1, max_epochs + 1):
        train_result = train_one_epoch(
            model, loaders["train"], optimizer, criterion, device, normalization
        )
        val_result = evaluate(
            model, loaders["val"], criterion, device, normalization, include_names=False
        )
        val_rmse = float(val_result["rmse"])
        learning_rate_now = float(optimizer.param_groups[0]["lr"])
        history.append(
            {
                "epoch": epoch,
                "train_rmse_eV": float(train_result["rmse"]),
                "val_rmse_eV": val_rmse,
                "learning_rate": learning_rate_now,
                "train_mse_loss": float(train_result["mse_loss"]),
                "val_mse_loss": float(val_result["mse_loss"]),
            }
        )
        logger.info(
            "Epoch %04d | train RMSE %.6f eV | val RMSE %.6f eV | lr %.6e",
            epoch,
            train_result["rmse"],
            val_rmse,
            learning_rate_now,
        )
        scheduler.step(val_rmse)

        if val_rmse < best_rmse - MIN_DELTA:
            best_rmse = val_rmse
            best_epoch = epoch
            epochs_without_improvement = 0
            if write_checkpoint:
                save_checkpoint(
                    CHECKPOINT_PATH,
                    model,
                    model_config,
                    graph_config,
                    normalization,
                    epoch,
                    val_rmse,
                    optimizer,
                    hyperparameters,
                )
        else:
            epochs_without_improvement += 1

        if trial is not None:
            trial.report(best_rmse, step=epoch)
            if trial.should_prune():
                import optuna

                raise optuna.TrialPruned()
        if epochs_without_improvement >= patience:
            stopped_early = True
            logger.info("Early stopping triggered at epoch %d.", epoch)
            break

    if best_epoch == 0:
        raise RuntimeError("Training ended without a finite validation RMSE.")
    return best_rmse, best_epoch, history, stopped_early


def run_optuna(
    records: Sequence[Dict[str, object]],
    split_indices: Dict[str, List[int]],
    graph_config: GraphConfig,
    device: torch.device,
    logger: logging.Logger,
) -> Dict[str, object]:
    try:
        import optuna
    except ImportError as exc:
        raise ImportError(
            "Optuna is required when USE_OPTUNA=True. Install it with: pip install optuna"
        ) from exc

    optuna.logging.set_verbosity(optuna.logging.WARNING)

    def objective(trial):
        set_global_seed(SEED)
        batch_size = trial.suggest_categorical("batch_size", [16, 32, 64])
        learning_rate = trial.suggest_float("learning_rate", 1.0e-4, 3.0e-3, log=True)
        weight_decay = trial.suggest_float("weight_decay", 1.0e-7, 1.0e-3, log=True)
        hidden_dim = trial.suggest_categorical("hidden_dim", [128, 192, 256])
        conv_layers = trial.suggest_int("conv_layers", 3, 5)
        downprojection_dim = trial.suggest_categorical(
            "downprojection_dim", [32, 64, 96]
        )
        dropout = trial.suggest_float("dropout", 0.0, 0.2)
        loaders, normalization = make_data_loaders(records, split_indices, batch_size)
        config = ModelConfig(
            atom_input_features=graph_config.max_atomic_number,
            hidden_dim=hidden_dim,
            conv_layers=conv_layers,
            downprojection_dim=downprojection_dim,
            reciprocal_hidden_layers=RECIPROCAL_HIDDEN_LAYERS,
            dropout=dropout,
            rbf_min=RBF_MIN,
            rbf_max=RBF_MAX,
        )
        best_rmse, _, _, _ = fit_model(
            loaders,
            normalization,
            config,
            graph_config,
            device,
            learning_rate,
            weight_decay,
            OPTUNA_MAX_EPOCHS,
            min(PATIENCE, 30),
            logger,
            write_checkpoint=False,
            trial=trial,
        )
        del loaders
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


def save_prediction_dat(split: str, result: Dict[str, object]) -> None:
    y_true = np.asarray(result["y_true"])
    y_pred = np.asarray(result["y_pred"])
    frame = pd.DataFrame(
        {
            "cif_name": result["cif_name"],
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


def plot_parity(
    split: str,
    result: Dict[str, object],
    axis_limits: Tuple[float, float],
) -> None:
    color_index = {"train": 0, "val": 1, "test": 2}[split]
    y_true = np.asarray(result["y_true"])
    y_pred = np.asarray(result["y_pred"])
    fig, axis = plt.subplots(figsize=(6.2, 6.0))
    axis.scatter(
        y_true, y_pred, s=24, alpha=0.78, color=COLORS[color_index], edgecolors="none"
    )
    axis.plot(axis_limits, axis_limits, color="black", linestyle="--", linewidth=1.2)
    axis.set_xlim(axis_limits)
    axis.set_ylim(axis_limits)
    axis.set_aspect("equal", adjustable="box")
    axis.set_title(f"{MODEL_NAME}: {split.capitalize()} Set", fontsize=TITLE_FONTSIZE)
    axis.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    axis.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    axis.tick_params(labelsize=TICK_FONTSIZE)
    if USE_GRID:
        axis.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)
    annotation = (
        f"MAE = {result['mae']:.4f} {TARGET_UNIT}\n"
        f"RMSE = {result['rmse']:.4f} {TARGET_UNIT}\n"
        f"$R^2$ = {result['r2']:.4f}"
    )
    axis.text(
        0.04,
        0.96,
        annotation,
        transform=axis.transAxes,
        va="top",
        fontsize=ANNOTATION_FONTSIZE,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.85},
    )
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_parity_{split}.jpg",
        dpi=FIG_DPI,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_combined_parity(
    results: Dict[str, Dict[str, object]], axis_limits: Tuple[float, float]
) -> None:
    fig, axis = plt.subplots(figsize=(6.4, 6.1))
    for color_index, split in enumerate(("train", "val", "test")):
        result = results[split]
        axis.scatter(
            result["y_true"],
            result["y_pred"],
            s=22,
            alpha=0.72,
            color=COLORS[color_index],
            edgecolors="none",
            label=f"{split.capitalize()} (RMSE={result['rmse']:.4f} {TARGET_UNIT})",
        )
    axis.plot(
        axis_limits,
        axis_limits,
        color="black",
        linestyle="--",
        linewidth=1.2,
        label="y = x",
    )
    axis.set_xlim(axis_limits)
    axis.set_ylim(axis_limits)
    axis.set_aspect("equal", adjustable="box")
    axis.set_title(f"{MODEL_NAME}: All Splits", fontsize=TITLE_FONTSIZE)
    axis.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    axis.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    axis.tick_params(labelsize=TICK_FONTSIZE)
    axis.legend(fontsize=LEGEND_FONTSIZE)
    if USE_GRID:
        axis.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_parity_all.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(fig)

    frames = []
    for split in ("train", "val", "test"):
        result = results[split]
        y_true = np.asarray(result["y_true"])
        y_pred = np.asarray(result["y_pred"])
        frames.append(
            pd.DataFrame(
                {
                    "split": split,
                    "cif_name": result["cif_name"],
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


def plot_rmse_curve(history: Sequence[Dict[str, float]], best_epoch: int) -> None:
    frame = pd.DataFrame(history)
    frame.to_csv(
        DAT_DIR / f"{MODEL_NAME}_rmse_curve.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )
    fig, axis = plt.subplots(figsize=(7.2, 5.2))
    axis.plot(
        frame["epoch"], frame["train_rmse_eV"], color=COLORS[0], label="Train RMSE"
    )
    axis.plot(
        frame["epoch"], frame["val_rmse_eV"], color=COLORS[1], label="Validation RMSE"
    )
    best_row = frame.loc[frame["epoch"] == best_epoch].iloc[0]
    axis.axvline(
        best_epoch,
        color=COLORS[3],
        linestyle="--",
        linewidth=1.2,
        label=f"Best epoch: {best_epoch}",
    )
    axis.scatter(
        [best_epoch], [best_row["val_rmse_eV"]], color=COLORS[3], s=40, zorder=5
    )
    axis.set_title(f"{MODEL_NAME} Training History", fontsize=TITLE_FONTSIZE)
    axis.set_xlabel("Epoch", fontsize=LABEL_FONTSIZE)
    axis.set_ylabel(f"RMSE ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    axis.tick_params(labelsize=TICK_FONTSIZE)
    axis.legend(fontsize=LEGEND_FONTSIZE)
    if USE_GRID:
        axis.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_rmse_curve.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(fig)


def save_metrics_table(results: Dict[str, Dict[str, object]]) -> None:
    rows = []
    for split in ("train", "val", "test"):
        result = results[split]
        rows.append(
            {
                "split": split,
                "n_samples": len(result["y_true"]),
                "mae_eV": result["mae"],
                "rmse_eV": result["rmse"],
                "r2": result["r2"],
            }
        )
    pd.DataFrame(rows).to_csv(
        TABLE_DIR / f"{MODEL_NAME}_metrics.dat",
        sep="\t",
        index=False,
        float_format="%.6f",
    )


def generate_all_outputs(
    results: Dict[str, Dict[str, object]],
    history: Sequence[Dict[str, float]],
    best_epoch: int,
) -> None:
    combined_values = np.concatenate(
        [
            np.asarray(results[split][key])
            for split in ("train", "val", "test")
            for key in ("y_true", "y_pred")
        ]
    )
    minimum = float(np.min(combined_values))
    maximum = float(np.max(combined_values))
    padding = max(0.05 * (maximum - minimum), 0.05)
    axis_limits = (minimum - padding, maximum + padding)
    for split in ("train", "val", "test"):
        save_prediction_dat(split, results[split])
        plot_parity(split, results[split], axis_limits)
    plot_combined_parity(results, axis_limits)
    plot_rmse_curve(history, best_epoch)
    save_metrics_table(results)


@torch.inference_mode()
def predict_cifs(
    checkpoint_path: Union[str, Path],
    cif_paths: Sequence[Union[str, Path]],
    device: Optional[torch.device] = None,
    batch_size: int = BATCH_SIZE,
) -> pd.DataFrame:
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model, checkpoint = load_trained_model(checkpoint_path, device)
    graph_config = GraphConfig(**checkpoint["graph_config"])
    graphs: List[Data] = []
    errors: List[str] = []
    for value in cif_paths:
        path = Path(value).expanduser().resolve()
        try:
            graph = build_graph(path, path.name, graph_config)
            graphs.append(graph)
        except Exception as exc:
            errors.append(f"{path}: {exc}")
    if not graphs:
        raise ValueError("No valid CIF files were supplied. " + " | ".join(errors))

    loader = DataLoader(graphs, batch_size=batch_size, shuffle=False)
    predictions: List[float] = []
    names: List[str] = []
    normalization = checkpoint["normalization"]
    model.eval()
    for batch in loader:
        batch = batch.to(device)
        values = model(batch).float().cpu().numpy()
        values = values * float(normalization["std"]) + float(normalization["mean"])
        predictions.extend(values.tolist())
        batch_names = batch.cif_name
        names.extend(
            [str(name) for name in batch_names]
            if isinstance(batch_names, (list, tuple))
            else [str(batch_names)]
        )
    frame = pd.DataFrame(
        {"cif_name": names, f"predicted_{TARGET_NAME}_{TARGET_UNIT}": predictions}
    )
    frame.attrs["errors"] = errors
    return frame


def log_run_configuration(
    logger: logging.Logger,
    model: nn.Module,
    model_config: ModelConfig,
    graph_config: GraphConfig,
    split_indices: Dict[str, List[int]],
    normalization: Dict[str, float],
    device: torch.device,
) -> None:
    total_parameters, trainable_parameters = count_parameters(model)
    logger.info("Model name: %s", MODEL_NAME)
    logger.info("Run version: %s", RUN_VERSION)
    logger.info("Output directory: %s", OUTPUT_DIR)
    logger.info("Excel path: %s", resolve_path(EXCEL_PATH))
    logger.info("CIF directory: %s", resolve_path(CIF_DIR))
    logger.info("Excel columns: first column = CIF filename; second column = target")
    logger.info("Target: %s (%s)", TARGET_NAME, TARGET_UNIT)
    logger.info("Seed: %d", SEED)
    logger.info(
        "Split ratios: train=%.3f, val=%.3f, test=%.3f",
        TRAIN_RATIO,
        VAL_RATIO,
        TEST_RATIO,
    )
    logger.info(
        "Split sizes: train=%d, val=%d, test=%d",
        len(split_indices["train"]),
        len(split_indices["val"]),
        len(split_indices["test"]),
    )
    logger.info("Target normalization fitted on train only: %s", normalization)
    logger.info("Model configuration: %s", asdict(model_config))
    logger.info("Graph configuration: %s", asdict(graph_config))
    logger.info("Loss: MSELoss")
    logger.info("Optimizer: AdamW")
    logger.info("Scheduler: ReduceLROnPlateau")
    logger.info("Scatter backend: %s", SCATTER_BACKEND)
    logger.info("Device information: %s", device_information(device))
    logger.info("Total parameters: %d", total_parameters)
    logger.info("Trainable parameters: %d", trainable_parameters)
    logger.info("Model architecture:\n%s", model)


def main() -> None:
    os.chdir(SCRIPT_DIR)
    create_output_directories()
    logger = setup_logger()
    total_start = time.perf_counter()
    set_global_seed(SEED)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    graph_config = GraphConfig()

    data_start = time.perf_counter()
    raw_records = read_excel_records(logger)
    records = validate_and_cache_records(raw_records, graph_config, logger)
    split_indices = create_or_load_split(records, logger)
    data_elapsed = time.perf_counter() - data_start

    optuna_result = None
    if USE_OPTUNA:
        optuna_result = run_optuna(records, split_indices, graph_config, device, logger)
        best_params = optuna_result["best_params"]
        batch_size = int(best_params["batch_size"])
        learning_rate = float(best_params["learning_rate"])
        weight_decay = float(best_params["weight_decay"])
        model_config = ModelConfig(
            atom_input_features=graph_config.max_atomic_number,
            hidden_dim=int(best_params["hidden_dim"]),
            conv_layers=int(best_params["conv_layers"]),
            downprojection_dim=int(best_params["downprojection_dim"]),
            reciprocal_hidden_layers=RECIPROCAL_HIDDEN_LAYERS,
            dropout=float(best_params["dropout"]),
            rbf_min=RBF_MIN,
            rbf_max=RBF_MAX,
        )
    else:
        batch_size = BATCH_SIZE
        learning_rate = LEARNING_RATE
        weight_decay = WEIGHT_DECAY
        model_config = ModelConfig(atom_input_features=graph_config.max_atomic_number)

    loaders, normalization = make_data_loaders(records, split_indices, batch_size)
    preview_model = ReciNet(model_config).to(device)
    log_run_configuration(
        logger,
        preview_model,
        model_config,
        graph_config,
        split_indices,
        normalization,
        device,
    )
    del preview_model

    training_start = time.perf_counter()
    best_val_rmse, best_epoch, history, stopped_early = fit_model(
        loaders,
        normalization,
        model_config,
        graph_config,
        device,
        learning_rate,
        weight_decay,
        MAX_EPOCHS,
        PATIENCE,
        logger,
        write_checkpoint=True,
    )
    training_elapsed = time.perf_counter() - training_start
    logger.info("Best epoch: %d", best_epoch)
    logger.info("Best validation RMSE: %.6f eV", best_val_rmse)
    logger.info("Early stopping: %s", stopped_early)

    evaluation_start = time.perf_counter()
    best_model, checkpoint = load_trained_model(CHECKPOINT_PATH, device)
    criterion = nn.MSELoss()
    results = {
        split: evaluate(
            best_model,
            loaders[split],
            criterion,
            device,
            checkpoint["normalization"],
            include_names=True,
        )
        for split in ("train", "val", "test")
    }
    generate_all_outputs(results, history, int(checkpoint["best_epoch"]))
    evaluation_elapsed = time.perf_counter() - evaluation_start

    for split in ("train", "val", "test"):
        result = results[split]
        logger.info(
            "%s metrics | n=%d | MAE=%.6f eV | RMSE=%.6f eV | R2=%.6f",
            split.capitalize(),
            len(result["y_true"]),
            result["mae"],
            result["rmse"],
            result["r2"],
        )
    logger.info("Optuna enabled: %s", USE_OPTUNA)
    if optuna_result is not None:
        logger.info("Optuna result: %s", optuna_result)
    logger.info(
        "Data loading and graph construction time: %s", format_duration(data_elapsed)
    )
    logger.info("Training time: %s", format_duration(training_elapsed))
    logger.info("Evaluation and plotting time: %s", format_duration(evaluation_elapsed))
    logger.info("Total runtime: %s", format_duration(time.perf_counter() - total_start))
    logger.info("Best checkpoint: %s", CHECKPOINT_PATH)


if __name__ == "__main__":
    main()

from __future__ import annotations

import hashlib
import json
import logging
import math
import os
import random
import sys
import time
import warnings
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence, Tuple

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
from torch.optim import AdamW
from torch.optim.lr_scheduler import ReduceLROnPlateau
from torch.utils.data import Dataset
from torch_geometric.data import Data
from torch_geometric.loader import DataLoader
from torch_geometric.utils import scatter
from tqdm.auto import tqdm


# ============================== User configuration ==============================

MODEL_NAME = "M3GNet"
RUN_VERSION = "v1"

EXCEL_PATH = "./data.xlsx"
CIF_DIR = "./cif"

TARGET_NAME = "bandgap"
TARGET_UNIT = "eV"

SEED = 42
TRAIN_RATIO = 0.80
VAL_RATIO = 0.10
TEST_RATIO = 0.10

CUTOFF = 5.0
THREEBODY_CUTOFF = 4.0
MAX_NEIGHBORS: Optional[int] = None
MAX_ATOMIC_NUMBER = 118
MAX_N = 3
MAX_L = 3
N_BLOCKS = 3
UNITS = 64

BATCH_SIZE = 16
MAX_EPOCHS = 300
LEARNING_RATE = 1.0e-3
WEIGHT_DECAY = 1.0e-5
PATIENCE = 60
MIN_DELTA = 1.0e-5
LR_PATIENCE = 15
LR_FACTOR = 0.5
MIN_LR = 1.0e-6
GRAD_CLIP_NORM = 5.0
NUM_WORKERS = 0
PIN_MEMORY = True
USE_BF16 = True

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

# ================================================================================


class M3GNetData(Data):
    """PyG data object with edge-index-aware triplet batching."""

    def __inc__(self, key: str, value: Any, *args: Any, **kwargs: Any) -> Any:
        if key == "triplet_edge_index":
            return int(self.edge_index.size(1))
        return super().__inc__(key, value, *args, **kwargs)


class CrystalGraphDataset(Dataset):
    def __init__(self, records: Sequence[Dict[str, Any]]) -> None:
        self.records = list(records)

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, index: int) -> M3GNetData:
        return self.records[index]["graph"]


class MLP(nn.Module):
    def __init__(
        self,
        in_dim: int,
        hidden_dims: Sequence[int],
        out_dim: int,
        final_activation: bool = False,
        bias: bool = True,
    ) -> None:
        super().__init__()
        dims = [in_dim, *hidden_dims, out_dim]
        layers: List[nn.Module] = []
        for i in range(len(dims) - 1):
            layers.append(nn.Linear(dims[i], dims[i + 1], bias=bias))
            if i < len(dims) - 2 or final_activation:
                layers.append(nn.SiLU())
        self.net = nn.Sequential(*layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class GatedMLP(nn.Module):
    def __init__(
        self, in_dim: int, hidden_dims: Sequence[int], out_dim: int, bias: bool = True
    ) -> None:
        super().__init__()
        self.value = MLP(in_dim, hidden_dims, out_dim, final_activation=True, bias=bias)
        gate_dims = [in_dim, *hidden_dims, out_dim]
        gate_layers: List[nn.Module] = []
        for i in range(len(gate_dims) - 1):
            gate_layers.append(nn.Linear(gate_dims[i], gate_dims[i + 1], bias=bias))
            gate_layers.append(nn.SiLU() if i < len(gate_dims) - 2 else nn.Sigmoid())
        self.gate = nn.Sequential(*gate_layers)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.value(x) * self.gate(x)


def polynomial_cutoff(r: torch.Tensor, cutoff: float) -> torch.Tensor:
    ratio = r / cutoff
    value = 1.0 - 6.0 * ratio.pow(5) + 15.0 * ratio.pow(4) - 10.0 * ratio.pow(3)
    return torch.where(r <= cutoff, value, torch.zeros_like(value))


def spherical_bessel_j(l_value: int, x: torch.Tensor) -> torch.Tensor:
    eps = torch.finfo(x.dtype).eps
    safe_x = x.clamp_min(eps)
    j0 = torch.sin(safe_x) / safe_x
    if l_value == 0:
        return torch.where(x.abs() < 1.0e-4, torch.ones_like(x), j0)
    j1 = torch.sin(safe_x) / safe_x.pow(2) - torch.cos(safe_x) / safe_x
    if l_value == 1:
        near_zero = x / 3.0
        return torch.where(x.abs() < 1.0e-4, near_zero, j1)
    previous, current = j0, j1
    for order in range(1, l_value):
        next_value = (2 * order + 1) * current / safe_x - previous
        previous, current = current, next_value
    if l_value == 2:
        current = torch.where(x.abs() < 1.0e-3, x.pow(2) / 15.0, current)
    return current


def default_spherical_bessel_roots(max_l: int, max_n: int) -> np.ndarray:
    known = np.array(
        [
            [math.pi, 2.0 * math.pi, 3.0 * math.pi, 4.0 * math.pi, 5.0 * math.pi],
            [4.493409458, 7.725251837, 10.904121659, 14.066193913, 17.220755272],
            [5.763459197, 9.095011331, 12.322940971, 15.514603011, 18.689036355],
            [6.987932001, 10.417118547, 13.698023153, 16.923621285, 20.121806174],
            [8.182561453, 11.704907155, 15.039664708, 18.301255960, 21.525417734],
        ],
        dtype=np.float32,
    )
    if max_l > known.shape[0] or max_n > known.shape[1]:
        raise ValueError(
            "MAX_L and MAX_N must not exceed 5 in this standalone implementation."
        )
    return known[:max_l, :max_n]


class SphericalBesselBasis(nn.Module):
    def __init__(self, max_l: int, max_n: int, cutoff: float) -> None:
        super().__init__()
        roots = torch.tensor(default_spherical_bessel_roots(max_l, max_n))
        self.register_buffer("roots", roots)
        self.max_l = max_l
        self.max_n = max_n
        self.cutoff = float(cutoff)

        normalizers = []
        for l_value in range(max_l):
            root = roots[l_value]
            denom = spherical_bessel_j(l_value + 1, root).abs().clamp_min(1.0e-8)
            normalizers.append(math.sqrt(2.0 / cutoff**3) / denom)
        self.register_buffer("normalizers", torch.stack(normalizers, dim=0))

    @property
    def output_dim(self) -> int:
        return self.max_l * self.max_n

    def forward(self, distances: torch.Tensor) -> torch.Tensor:
        outputs = []
        for l_value in range(self.max_l):
            x = distances[:, None] * self.roots[l_value][None, :] / self.cutoff
            outputs.append(
                spherical_bessel_j(l_value, x) * self.normalizers[l_value][None, :]
            )
        return torch.cat(outputs, dim=-1)


def legendre_polynomial(l_value: int, x: torch.Tensor) -> torch.Tensor:
    if l_value == 0:
        return torch.ones_like(x)
    if l_value == 1:
        return x
    p_previous = torch.ones_like(x)
    p_current = x
    for order in range(1, l_value):
        p_next = ((2 * order + 1) * x * p_current - order * p_previous) / (order + 1)
        p_previous, p_current = p_current, p_next
    return p_current


class ThreeBodyBasis(nn.Module):
    def __init__(self, max_l: int, max_n: int, cutoff: float) -> None:
        super().__init__()
        self.radial = SphericalBesselBasis(max_l=max_l, max_n=max_n, cutoff=cutoff)
        self.max_l = max_l
        self.max_n = max_n

    @property
    def output_dim(self) -> int:
        return self.max_l * self.max_n

    def forward(
        self,
        edge_dist: torch.Tensor,
        edge_vec: torch.Tensor,
        triplet_edge_index: torch.Tensor,
    ) -> torch.Tensor:
        if triplet_edge_index.numel() == 0:
            return edge_dist.new_zeros((0, self.output_dim))
        edge_ij, edge_ik = triplet_edge_index[0], triplet_edge_index[1]
        v_ij = edge_vec[edge_ij]
        v_ik = edge_vec[edge_ik]
        cosine = F.cosine_similarity(v_ij, v_ik, dim=-1, eps=1.0e-8).clamp(-1.0, 1.0)
        radial = self.radial(edge_dist[edge_ik]).view(-1, self.max_l, self.max_n)
        angular = []
        for l_value in range(self.max_l):
            norm = math.sqrt((2 * l_value + 1) / (4.0 * math.pi))
            angular.append(norm * legendre_polynomial(l_value, cosine))
        angular_tensor = torch.stack(angular, dim=1).unsqueeze(-1)
        return (radial * angular_tensor).reshape(-1, self.output_dim)


class ThreeDInteraction(nn.Module):
    def __init__(self, units: int, basis_dim: int) -> None:
        super().__init__()
        self.atom_gate = nn.Sequential(nn.Linear(units, basis_dim), nn.Sigmoid())
        self.update = GatedMLP(basis_dim, [], units, bias=False)

    def forward(
        self,
        atoms: torch.Tensor,
        bonds: torch.Tensor,
        edge_index: torch.Tensor,
        triplet_edge_index: torch.Tensor,
        three_basis: torch.Tensor,
        three_cutoff: torch.Tensor,
    ) -> torch.Tensor:
        if triplet_edge_index.numel() == 0:
            return bonds
        edge_ij, edge_ik = triplet_edge_index[0], triplet_edge_index[1]
        end_atoms = edge_index[1, edge_ik]
        basis = three_basis * self.atom_gate(atoms[end_atoms])
        weights = three_cutoff[edge_ij] * three_cutoff[edge_ik]
        messages = basis * weights[:, None]
        aggregated = scatter(
            messages, edge_ij, dim=0, dim_size=bonds.size(0), reduce="sum"
        )
        return bonds + self.update(aggregated)


class M3GNetBlock(nn.Module):
    def __init__(self, units: int, radial_dim: int) -> None:
        super().__init__()
        self.bond_update = GatedMLP(3 * units, [units], units)
        self.bond_weight = nn.Linear(radial_dim, units, bias=False)
        self.atom_update = GatedMLP(3 * units, [units], units)
        self.atom_weight = nn.Sequential(
            nn.Linear(radial_dim, units, bias=False), nn.SiLU()
        )

    def forward(
        self,
        atoms: torch.Tensor,
        bonds: torch.Tensor,
        edge_index: torch.Tensor,
        radial_basis: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        center, neighbor = edge_index[0], edge_index[1]
        bond_input = torch.cat([atoms[center], atoms[neighbor], bonds], dim=-1)
        bonds = bonds + self.bond_update(bond_input) * self.bond_weight(radial_basis)

        atom_input = torch.cat([atoms[center], atoms[neighbor], bonds], dim=-1)
        messages = self.atom_update(atom_input) * self.atom_weight(radial_basis)
        atoms = atoms + scatter(
            messages, center, dim=0, dim_size=atoms.size(0), reduce="sum"
        )
        return atoms, bonds


class WeightedAtomReadout(nn.Module):
    def __init__(self, units: int) -> None:
        super().__init__()
        self.value = MLP(units, [units], units, final_activation=True)
        self.weight = nn.Sequential(
            nn.Linear(units, units), nn.SiLU(), nn.Linear(units, 1), nn.Sigmoid()
        )

    def forward(
        self, atoms: torch.Tensor, batch: torch.Tensor, num_graphs: int
    ) -> torch.Tensor:
        values = self.value(atoms)
        weights = self.weight(atoms)
        denominator = scatter(
            weights, batch, dim=0, dim_size=num_graphs, reduce="sum"
        ).clamp_min(1.0e-8)
        normalized = weights / denominator[batch]
        return scatter(
            normalized * values, batch, dim=0, dim_size=num_graphs, reduce="sum"
        )


class M3GNetRegressor(nn.Module):
    """M3GNet-style three-body graph network for intensive regression."""

    def __init__(
        self,
        max_atomic_number: int = MAX_ATOMIC_NUMBER,
        max_n: int = MAX_N,
        max_l: int = MAX_L,
        n_blocks: int = N_BLOCKS,
        units: int = UNITS,
        cutoff: float = CUTOFF,
        threebody_cutoff: float = THREEBODY_CUTOFF,
    ) -> None:
        super().__init__()
        self.model_config = {
            "max_atomic_number": int(max_atomic_number),
            "max_n": int(max_n),
            "max_l": int(max_l),
            "n_blocks": int(n_blocks),
            "units": int(units),
            "cutoff": float(cutoff),
            "threebody_cutoff": float(threebody_cutoff),
        }
        self.atom_embedding = nn.Embedding(max_atomic_number + 1, units, padding_idx=0)
        self.radial_basis = SphericalBesselBasis(
            max_l=max_l, max_n=max_n, cutoff=cutoff
        )
        self.three_basis = ThreeBodyBasis(max_l=max_l, max_n=max_n, cutoff=cutoff)
        radial_dim = self.radial_basis.output_dim
        basis_dim = self.three_basis.output_dim
        self.bond_featurizer = nn.Sequential(
            nn.Linear(radial_dim, units, bias=False), nn.SiLU()
        )
        self.three_interactions = nn.ModuleList(
            [ThreeDInteraction(units, basis_dim) for _ in range(n_blocks)]
        )
        self.graph_blocks = nn.ModuleList(
            [M3GNetBlock(units, radial_dim) for _ in range(n_blocks)]
        )
        self.readout = WeightedAtomReadout(units)
        self.output_mlp = MLP(units, [units, units], 1)

    def forward(self, data: Data) -> torch.Tensor:
        atoms = self.atom_embedding(data.z)
        radial = self.radial_basis(data.edge_dist)
        bonds = self.bond_featurizer(radial)
        three_basis = self.three_basis(
            data.edge_dist, data.edge_vec, data.triplet_edge_index
        )
        three_cutoff = polynomial_cutoff(
            data.edge_dist, self.model_config["threebody_cutoff"]
        )

        for interaction, graph_block in zip(self.three_interactions, self.graph_blocks):
            bonds = interaction(
                atoms,
                bonds,
                data.edge_index,
                data.triplet_edge_index,
                three_basis,
                three_cutoff,
            )
            atoms, bonds = graph_block(atoms, bonds, data.edge_index, radial)

        batch = getattr(data, "batch", None)
        if batch is None:
            batch = torch.zeros(atoms.size(0), dtype=torch.long, device=atoms.device)
        num_graphs = int(data.num_graphs) if hasattr(data, "num_graphs") else 1
        graph_features = self.readout(atoms, batch, num_graphs)
        return self.output_mlp(graph_features).view(-1)


def ensure_directories() -> None:
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


def setup_logger() -> logging.Logger:
    logger = logging.getLogger(f"{MODEL_NAME}_{RUN_VERSION}")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    file_handler = logging.FileHandler(LOG_PATH, mode="w", encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(sys.stdout)
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
    try:
        torch.use_deterministic_algorithms(True, warn_only=True)
    except TypeError:
        torch.use_deterministic_algorithms(True)


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
        props = torch.cuda.get_device_properties(0)
        logger.info("GPU: %s", torch.cuda.get_device_name(0))
        logger.info("GPU memory: %.3f GiB", props.total_memory / 1024**3)
        logger.info("CUDA runtime: %s", torch.version.cuda)
        logger.info(
            "BF16 enabled: %s", bool(USE_BF16 and torch.cuda.is_bf16_supported())
        )
    return device


def read_excel_rows(logger: logging.Logger) -> List[Dict[str, Any]]:
    path = Path(EXCEL_PATH)
    if not path.is_file():
        raise FileNotFoundError(f"Excel file not found: {path}")
    frame = pd.read_excel(path)
    if frame.shape[1] < 2:
        raise ValueError("The Excel file must contain at least two columns.")

    rows: List[Dict[str, Any]] = []
    seen: set[str] = set()
    for row_number, (_, row) in enumerate(frame.iterrows(), start=2):
        cif_value = row.iloc[0]
        target_value = row.iloc[1]
        if pd.isna(cif_value):
            logger.warning("Skipped Excel row %d: empty CIF filename.", row_number)
            continue
        cif_name = f"{int(float(cif_value))}.cif"
        if not cif_name:
            logger.warning("Skipped Excel row %d: empty CIF filename.", row_number)
            continue
        if cif_name in seen:
            logger.warning(
                "Skipped duplicate CIF entry at row %d: %s", row_number, cif_name
            )
            continue
        try:
            target = float(target_value)
        except (TypeError, ValueError):
            logger.warning("Skipped %s: target is not numeric.", cif_name)
            continue
        if not np.isfinite(target):
            logger.warning("Skipped %s: target is NaN or infinite.", cif_name)
            continue
        cif_path = Path(CIF_DIR) / cif_name
        if not cif_path.is_file():
            logger.warning("Skipped %s: CIF file not found at %s", cif_name, cif_path)
            continue
        seen.add(cif_name)
        rows.append({"cif_name": cif_name, "cif_path": str(cif_path), "target": target})
    if len(rows) < 10:
        raise ValueError(
            f"Only {len(rows)} valid rows remain; at least 10 are required."
        )
    return rows


def limit_neighbors(
    centers: np.ndarray,
    neighbors: np.ndarray,
    offsets: np.ndarray,
    distances: np.ndarray,
    max_neighbors: Optional[int],
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    if max_neighbors is None:
        return centers, neighbors, offsets, distances
    keep: List[int] = []
    for center in np.unique(centers):
        local = np.where(centers == center)[0]
        ordered = local[np.argsort(distances[local], kind="stable")]
        keep.extend(ordered[:max_neighbors].tolist())
    keep_array = np.asarray(keep, dtype=np.int64)
    return (
        centers[keep_array],
        neighbors[keep_array],
        offsets[keep_array],
        distances[keep_array],
    )


def build_triplets(
    edge_index: np.ndarray, edge_dist: np.ndarray, threebody_cutoff: float
) -> np.ndarray:
    valid_edges = np.where(edge_dist <= threebody_cutoff + 1.0e-8)[0]
    by_center: Dict[int, List[int]] = {}
    for edge_id in valid_edges.tolist():
        by_center.setdefault(int(edge_index[0, edge_id]), []).append(edge_id)
    first: List[int] = []
    second: List[int] = []
    for edges in by_center.values():
        for edge_ij in edges:
            for edge_ik in edges:
                if edge_ij != edge_ik:
                    first.append(edge_ij)
                    second.append(edge_ik)
    if not first:
        return np.empty((2, 0), dtype=np.int64)
    return np.asarray([first, second], dtype=np.int64)


def structure_to_graph(
    structure: Structure,
    target: float,
    cif_name: str,
    cutoff: float,
    threebody_cutoff: float,
    max_neighbors: Optional[int],
) -> M3GNetData:
    atomic_numbers = np.asarray([site.specie.Z for site in structure], dtype=np.int64)
    if atomic_numbers.size == 0:
        raise ValueError("The structure contains no atoms.")
    if int(atomic_numbers.max()) > MAX_ATOMIC_NUMBER:
        raise ValueError(
            f"Atomic number exceeds MAX_ATOMIC_NUMBER={MAX_ATOMIC_NUMBER}."
        )

    centers, neighbors, offsets, distances = structure.get_neighbor_list(
        r=cutoff,
        numerical_tol=1.0e-8,
        exclude_self=True,
    )
    centers = np.asarray(centers, dtype=np.int64)
    neighbors = np.asarray(neighbors, dtype=np.int64)
    offsets = np.asarray(offsets, dtype=np.int64)
    distances = np.asarray(distances, dtype=np.float32)
    if centers.size == 0:
        raise ValueError(
            f"No periodic neighbors were found within {cutoff:.3f} angstrom."
        )

    centers, neighbors, offsets, distances = limit_neighbors(
        centers, neighbors, offsets, distances, max_neighbors
    )
    order = np.lexsort(
        (offsets[:, 2], offsets[:, 1], offsets[:, 0], distances, neighbors, centers)
    )
    centers, neighbors = centers[order], neighbors[order]
    offsets, distances = offsets[order], distances[order]
    edge_index = np.stack([centers, neighbors], axis=0)

    cart_coords = np.asarray(structure.cart_coords, dtype=np.float32)
    lattice = np.asarray(structure.lattice.matrix, dtype=np.float32)
    translated_neighbors = cart_coords[neighbors] + offsets.astype(np.float32) @ lattice
    edge_vec = translated_neighbors - cart_coords[centers]
    vector_distances = np.linalg.norm(edge_vec, axis=1)
    if not np.allclose(vector_distances, distances, rtol=1.0e-4, atol=1.0e-5):
        distances = vector_distances.astype(np.float32)

    triplets = build_triplets(edge_index, distances, threebody_cutoff)
    graph = M3GNetData(
        z=torch.from_numpy(atomic_numbers),
        edge_index=torch.from_numpy(edge_index).long(),
        edge_dist=torch.from_numpy(distances).float(),
        edge_vec=torch.from_numpy(edge_vec.astype(np.float32)).float(),
        triplet_edge_index=torch.from_numpy(triplets).long(),
        y=torch.tensor([target], dtype=torch.float32),
    )
    graph.cif_name = cif_name
    graph.num_nodes = int(atomic_numbers.size)
    return graph


def graph_cache_path() -> Path:
    excel = Path(EXCEL_PATH).resolve()
    stat = excel.stat()
    signature = {
        "excel": str(excel),
        "excel_size": stat.st_size,
        "excel_mtime_ns": stat.st_mtime_ns,
        "cif_dir": str(Path(CIF_DIR).resolve()),
        "cutoff": CUTOFF,
        "threebody_cutoff": THREEBODY_CUTOFF,
        "max_neighbors": MAX_NEIGHBORS,
    }
    digest = hashlib.sha256(
        json.dumps(signature, sort_keys=True).encode("utf-8")
    ).hexdigest()[:16]
    return CACHE_DIR / f"{MODEL_NAME}_graphs_{digest}.pt"


def load_or_build_graphs(
    rows: Sequence[Dict[str, Any]], logger: logging.Logger
) -> List[Dict[str, Any]]:
    cache_path = graph_cache_path()
    if cache_path.is_file():
        try:
            payload = torch.load(cache_path, map_location="cpu", weights_only=False)
            cached = payload["records"]
            cached_keys = [(item["cif_name"], float(item["target"])) for item in cached]
            row_keys = [(item["cif_name"], float(item["target"])) for item in rows]
            if cached_keys == row_keys:
                logger.info("Loaded %d graphs from cache: %s", len(cached), cache_path)
                return cached
            logger.warning(
                "Graph cache does not match the current valid rows; rebuilding it."
            )
        except Exception as exc:
            logger.warning("Could not load graph cache %s: %s", cache_path, exc)

    records: List[Dict[str, Any]] = []
    for row in tqdm(rows, desc="Building crystal graphs", leave=False):
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                structure = Structure.from_file(row["cif_path"])
            graph = structure_to_graph(
                structure=structure,
                target=float(row["target"]),
                cif_name=row["cif_name"],
                cutoff=CUTOFF,
                threebody_cutoff=THREEBODY_CUTOFF,
                max_neighbors=MAX_NEIGHBORS,
            )
            records.append(
                {
                    "cif_name": row["cif_name"],
                    "target": float(row["target"]),
                    "graph": graph,
                }
            )
        except Exception as exc:
            logger.warning(
                "Skipped %s during CIF parsing or graph construction: %s",
                row["cif_name"],
                exc,
            )
    if len(records) < 10:
        raise ValueError(
            f"Only {len(records)} valid crystal graphs remain; at least 10 are required."
        )
    torch.save({"records": records}, cache_path)
    logger.info("Built and cached %d crystal graphs: %s", len(records), cache_path)
    return records


def split_records(
    records: Sequence[Dict[str, Any]], logger: logging.Logger
) -> Tuple[List[Dict[str, Any]], List[Dict[str, Any]], List[Dict[str, Any]]]:
    if not math.isclose(
        TRAIN_RATIO + VAL_RATIO + TEST_RATIO, 1.0, rel_tol=0.0, abs_tol=1.0e-9
    ):
        raise ValueError("TRAIN_RATIO + VAL_RATIO + TEST_RATIO must equal 1.0.")
    record_map = {item["cif_name"]: item for item in records}

    if SPLIT_PATH.is_file():
        try:
            split_frame = pd.read_csv(SPLIT_PATH)
            required = {"cif_name", "target", "split"}
            if not required.issubset(split_frame.columns):
                raise ValueError("missing required columns")
            if set(split_frame["cif_name"].astype(str)) != set(record_map):
                raise ValueError("CIF names differ from the current valid data")
            target_map = dict(
                zip(
                    split_frame["cif_name"].astype(str),
                    split_frame["target"].astype(float),
                )
            )
            for name, item in record_map.items():
                if not np.isclose(
                    target_map[name], float(item["target"]), rtol=0.0, atol=1.0e-10
                ):
                    raise ValueError(f"target changed for {name}")
            output: Dict[str, List[Dict[str, Any]]] = {
                "train": [],
                "val": [],
                "test": [],
            }
            for _, row in split_frame.iterrows():
                split_name = str(row["split"])
                if split_name not in output:
                    raise ValueError(f"unknown split label: {split_name}")
                output[split_name].append(record_map[str(row["cif_name"])])
            if any(len(output[key]) == 0 for key in output):
                raise ValueError("one or more splits are empty")
            logger.info("Reused data split: %s", SPLIT_PATH)
            return output["train"], output["val"], output["test"]
        except Exception as exc:
            logger.warning(
                "Existing split could not be reused (%s); creating a new split.", exc
            )

    indices = np.arange(len(records))
    train_indices, remainder = train_test_split(
        indices,
        test_size=1.0 - TRAIN_RATIO,
        random_state=SEED,
        shuffle=True,
    )
    relative_test = TEST_RATIO / (VAL_RATIO + TEST_RATIO)
    val_indices, test_indices = train_test_split(
        remainder,
        test_size=relative_test,
        random_state=SEED,
        shuffle=True,
    )
    split_lookup: Dict[int, str] = {}
    split_lookup.update({int(index): "train" for index in train_indices})
    split_lookup.update({int(index): "val" for index in val_indices})
    split_lookup.update({int(index): "test" for index in test_indices})
    split_frame = pd.DataFrame(
        [
            {
                "cif_name": item["cif_name"],
                "target": float(item["target"]),
                "split": split_lookup[index],
            }
            for index, item in enumerate(records)
        ]
    )
    split_frame.to_csv(SPLIT_PATH, index=False)
    train = [records[int(index)] for index in train_indices]
    val = [records[int(index)] for index in val_indices]
    test = [records[int(index)] for index in test_indices]
    logger.info("Created data split: %s", SPLIT_PATH)
    return train, val, test


def make_loader(
    records: Sequence[Dict[str, Any]], batch_size: int, shuffle: bool
) -> DataLoader:
    generator = torch.Generator()
    generator.manual_seed(SEED)
    return DataLoader(
        CrystalGraphDataset(records),
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=NUM_WORKERS,
        pin_memory=bool(PIN_MEMORY and torch.cuda.is_available()),
        persistent_workers=bool(NUM_WORKERS > 0),
        generator=generator,
    )


def model_config_from_parameters(
    parameters: Optional[Dict[str, Any]] = None
) -> Dict[str, Any]:
    parameters = parameters or {}
    return {
        "max_atomic_number": MAX_ATOMIC_NUMBER,
        "max_n": MAX_N,
        "max_l": MAX_L,
        "n_blocks": int(parameters.get("n_blocks", N_BLOCKS)),
        "units": int(parameters.get("units", UNITS)),
        "cutoff": CUTOFF,
        "threebody_cutoff": THREEBODY_CUTOFF,
    }


def autocast_context(device: torch.device):
    enabled = bool(
        device.type == "cuda" and USE_BF16 and torch.cuda.is_bf16_supported()
    )
    return torch.autocast(
        device_type=device.type, dtype=torch.bfloat16, enabled=enabled
    )


def train_one_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    device: torch.device,
    target_mean: float,
    target_std: float,
) -> Tuple[float, float]:
    model.train()
    criterion = nn.MSELoss()
    squared_error_sum = 0.0
    scaled_loss_sum = 0.0
    sample_count = 0
    for batch in loader:
        batch = batch.to(device, non_blocking=True)
        raw_target = batch.y.view(-1)
        scaled_target = (raw_target - target_mean) / target_std
        optimizer.zero_grad(set_to_none=True)
        with autocast_context(device):
            prediction_scaled = model(batch).view(-1)
            loss = criterion(prediction_scaled, scaled_target)
        if not torch.isfinite(loss):
            raise FloatingPointError("A non-finite training loss was encountered.")
        loss.backward()
        if GRAD_CLIP_NORM > 0:
            torch.nn.utils.clip_grad_norm_(model.parameters(), GRAD_CLIP_NORM)
        optimizer.step()
        prediction = prediction_scaled.detach().float() * target_std + target_mean
        squared_error_sum += float(
            torch.sum((prediction - raw_target.float()) ** 2).item()
        )
        scaled_loss_sum += float(loss.detach().item()) * raw_target.numel()
        sample_count += raw_target.numel()
    return math.sqrt(squared_error_sum / sample_count), scaled_loss_sum / sample_count


@torch.inference_mode()
def evaluate_model(
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
    scaled_loss_sum = 0.0
    sample_count = 0
    for batch in loader:
        batch_names = batch.cif_name
        if isinstance(batch_names, str):
            batch_names = [batch_names]
        batch = batch.to(device, non_blocking=True)
        target = batch.y.view(-1)
        scaled_target = (target - target_mean) / target_std
        with autocast_context(device):
            prediction_scaled = model(batch).view(-1)
            loss = F.mse_loss(prediction_scaled, scaled_target)
        prediction = prediction_scaled.float() * target_std + target_mean
        names.extend([str(name) for name in batch_names])
        true_values.extend(target.detach().cpu().float().numpy().tolist())
        predictions.extend(prediction.detach().cpu().numpy().tolist())
        scaled_loss_sum += float(loss.item()) * target.numel()
        sample_count += target.numel()
    y_true = np.asarray(true_values, dtype=np.float64)
    y_pred = np.asarray(predictions, dtype=np.float64)
    return {
        "names": names,
        "y_true": y_true,
        "y_pred": y_pred,
        "rmse": float(np.sqrt(np.mean((y_true - y_pred) ** 2))),
        "scaled_mse": scaled_loss_sum / max(1, sample_count),
    }


def save_checkpoint(
    model: M3GNetRegressor,
    optimizer: torch.optim.Optimizer,
    epoch: int,
    best_val_rmse: float,
    target_mean: float,
    target_std: float,
    hyperparameters: Dict[str, Any],
) -> None:
    checkpoint = {
        "model_name": MODEL_NAME,
        "run_version": RUN_VERSION,
        "model_state_dict": model.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "model_config": model.model_config,
        "graph_config": {
            "cutoff": CUTOFF,
            "threebody_cutoff": THREEBODY_CUTOFF,
            "max_neighbors": MAX_NEIGHBORS,
        },
        "best_epoch": int(epoch),
        "best_val_rmse": float(best_val_rmse),
        "seed": SEED,
        "target_name": TARGET_NAME,
        "target_unit": TARGET_UNIT,
        "normalization": {"mean": float(target_mean), "std": float(target_std)},
        "hyperparameters": hyperparameters,
    }
    torch.save(checkpoint, CHECKPOINT_PATH)


def torch_load_checkpoint(
    path: Path, map_location: torch.device | str
) -> Dict[str, Any]:
    try:
        return torch.load(path, map_location=map_location, weights_only=False)
    except TypeError:
        return torch.load(path, map_location=map_location)


def load_trained_model(
    checkpoint_path: str | Path,
    device: Optional[torch.device] = None,
) -> Tuple[M3GNetRegressor, Dict[str, Any]]:
    device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")
    checkpoint = torch_load_checkpoint(Path(checkpoint_path), device)
    model = M3GNetRegressor(**checkpoint["model_config"])
    model.load_state_dict(checkpoint["model_state_dict"], strict=True)
    model.to(device)
    model.eval()
    return model, checkpoint


def fit_model(
    train_records: Sequence[Dict[str, Any]],
    val_records: Sequence[Dict[str, Any]],
    device: torch.device,
    logger: logging.Logger,
    parameters: Optional[Dict[str, Any]] = None,
    save_best: bool = True,
    max_epochs: int = MAX_EPOCHS,
    patience: int = PATIENCE,
    verbose: bool = True,
) -> Tuple[M3GNetRegressor, List[Dict[str, float]], float, int, bool, float, float]:
    set_global_seed(SEED)
    parameters = dict(parameters or {})
    batch_size = int(parameters.get("batch_size", BATCH_SIZE))
    learning_rate = float(parameters.get("learning_rate", LEARNING_RATE))
    weight_decay = float(parameters.get("weight_decay", WEIGHT_DECAY))
    model = M3GNetRegressor(**model_config_from_parameters(parameters)).to(device)
    train_loader = make_loader(train_records, batch_size=batch_size, shuffle=True)
    val_loader = make_loader(val_records, batch_size=batch_size, shuffle=False)
    optimizer = AdamW(model.parameters(), lr=learning_rate, weight_decay=weight_decay)
    scheduler = ReduceLROnPlateau(
        optimizer,
        mode="min",
        factor=LR_FACTOR,
        patience=LR_PATIENCE,
        min_lr=MIN_LR,
    )
    targets = np.asarray(
        [float(item["target"]) for item in train_records], dtype=np.float64
    )
    target_mean = float(targets.mean())
    target_std = float(targets.std(ddof=0))
    if not np.isfinite(target_std) or target_std < 1.0e-12:
        target_std = 1.0

    best_val_rmse = float("inf")
    best_epoch = 0
    epochs_without_improvement = 0
    history: List[Dict[str, float]] = []
    best_state: Optional[Dict[str, torch.Tensor]] = None
    early_stopped = False

    for epoch in range(1, max_epochs + 1):
        train_rmse, train_mse = train_one_epoch(
            model, train_loader, optimizer, device, target_mean, target_std
        )
        val_result = evaluate_model(model, val_loader, device, target_mean, target_std)
        val_rmse = float(val_result["rmse"])
        current_lr = float(optimizer.param_groups[0]["lr"])
        history.append(
            {
                "epoch": float(epoch),
                "train_rmse_eV": train_rmse,
                "val_rmse_eV": val_rmse,
                "learning_rate": current_lr,
                "train_mse_loss": train_mse,
                "val_mse_loss": float(val_result["scaled_mse"]),
            }
        )
        if verbose:
            logger.info(
                "Epoch %04d | train RMSE %.6f %s | val RMSE %.6f %s | lr %.6e",
                epoch,
                train_rmse,
                TARGET_UNIT,
                val_rmse,
                TARGET_UNIT,
                current_lr,
            )
        if not np.isfinite(val_rmse):
            raise FloatingPointError("A non-finite validation RMSE was encountered.")
        improved = val_rmse < best_val_rmse - MIN_DELTA
        if improved:
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
                    optimizer,
                    epoch,
                    best_val_rmse,
                    target_mean,
                    target_std,
                    {
                        **parameters,
                        "batch_size": batch_size,
                        "learning_rate": learning_rate,
                        "weight_decay": weight_decay,
                        "max_epochs": max_epochs,
                        "patience": patience,
                    },
                )
        else:
            epochs_without_improvement += 1
        scheduler.step(val_rmse)
        if epochs_without_improvement >= patience:
            early_stopped = True
            if verbose:
                logger.info("Early stopping triggered at epoch %d.", epoch)
            break

    if best_state is None:
        raise RuntimeError("Training did not produce a valid model state.")
    model.load_state_dict(best_state)
    model.to(device)
    return (
        model,
        history,
        best_val_rmse,
        best_epoch,
        early_stopped,
        target_mean,
        target_std,
    )


def run_optuna(
    train_records: Sequence[Dict[str, Any]],
    val_records: Sequence[Dict[str, Any]],
    device: torch.device,
    logger: logging.Logger,
) -> Dict[str, Any]:
    try:
        import optuna
    except ImportError as exc:
        raise ImportError("USE_OPTUNA=True requires the optuna package.") from exc

    optuna.logging.set_verbosity(optuna.logging.WARNING)

    def objective(trial: Any) -> float:
        parameters = {
            "units": trial.suggest_categorical("units", [48, 64, 96, 128]),
            "n_blocks": trial.suggest_int("n_blocks", 2, 4),
            "batch_size": trial.suggest_categorical("batch_size", [8, 16, 32]),
            "learning_rate": trial.suggest_float(
                "learning_rate", 2.0e-4, 3.0e-3, log=True
            ),
            "weight_decay": trial.suggest_float(
                "weight_decay", 1.0e-7, 1.0e-3, log=True
            ),
        }
        try:
            _, _, best_rmse, _, _, _, _ = fit_model(
                train_records,
                val_records,
                device,
                logger,
                parameters=parameters,
                save_best=False,
                max_epochs=OPTUNA_MAX_EPOCHS,
                patience=OPTUNA_PATIENCE,
                verbose=False,
            )
            return best_rmse
        except torch.cuda.OutOfMemoryError:
            if torch.cuda.is_available():
                torch.cuda.empty_cache()
            raise optuna.TrialPruned("CUDA out of memory")

    study = optuna.create_study(
        direction="minimize", sampler=optuna.samplers.TPESampler(seed=SEED)
    )
    study.optimize(
        objective, n_trials=OPTUNA_N_TRIALS, timeout=OPTUNA_TIMEOUT, gc_after_trial=True
    )
    result = {
        "best_value": float(study.best_value),
        "best_params": dict(study.best_params),
        "n_trials": len(study.trials),
    }
    output_path = TABLE_DIR / f"{MODEL_NAME}_optuna_best_params.json"
    output_path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    logger.info(
        "Optuna best validation RMSE: %.6f %s", result["best_value"], TARGET_UNIT
    )
    logger.info("Optuna best parameters: %s", result["best_params"])
    return dict(study.best_params)


def calculate_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    r2 = float(r2_score(y_true, y_pred)) if len(y_true) >= 2 else float("nan")
    return {
        "mae_eV": float(mean_absolute_error(y_true, y_pred)),
        "rmse_eV": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "r2": r2,
    }


def save_prediction_dat(split_name: str, result: Dict[str, Any]) -> None:
    y_true = result["y_true"]
    y_pred = result["y_pred"]
    frame = pd.DataFrame(
        {
            "cif_name": result["names"],
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


def save_all_prediction_dat(results: Dict[str, Dict[str, Any]]) -> None:
    frames = []
    for split_name in ["train", "val", "test"]:
        result = results[split_name]
        y_true = result["y_true"]
        y_pred = result["y_pred"]
        frames.append(
            pd.DataFrame(
                {
                    "split": split_name,
                    "cif_name": result["names"],
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
        float_format="%.10f",
    )


def parity_limits(results: Dict[str, Dict[str, Any]]) -> Tuple[float, float]:
    values = []
    for result in results.values():
        values.extend(result["y_true"].tolist())
        values.extend(result["y_pred"].tolist())
    minimum, maximum = float(np.min(values)), float(np.max(values))
    span = max(maximum - minimum, 1.0e-6)
    margin = 0.05 * span
    return minimum - margin, maximum + margin


def style_axes(ax: plt.Axes) -> None:
    ax.tick_params(axis="both", labelsize=TICK_FONTSIZE)
    if USE_GRID:
        ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)


def plot_parity(
    split_name: str,
    result: Dict[str, Any],
    metrics: Dict[str, float],
    limits: Tuple[float, float],
    color: str,
) -> None:
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    ax.scatter(
        result["y_true"],
        result["y_pred"],
        s=22,
        alpha=0.75,
        color=color,
        edgecolors="none",
    )
    ax.plot(limits, limits, linestyle="--", linewidth=1.2, color="black", label="y = x")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_title(
        f"{MODEL_NAME} {split_name.capitalize()} Parity", fontsize=TITLE_FONTSIZE
    )
    annotation = (
        f"MAE = {metrics['mae_eV']:.4f} {TARGET_UNIT}\n"
        f"RMSE = {metrics['rmse_eV']:.4f} {TARGET_UNIT}\n"
        f"R² = {metrics['r2']:.4f}"
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
    results: Dict[str, Dict[str, Any]],
    metrics: Dict[str, Dict[str, float]],
    limits: Tuple[float, float],
) -> None:
    fig, ax = plt.subplots(figsize=(6.6, 6.2))
    for index, split_name in enumerate(["train", "val", "test"]):
        result = results[split_name]
        label = f"{split_name.capitalize()} (RMSE={metrics[split_name]['rmse_eV']:.4f})"
        ax.scatter(
            result["y_true"],
            result["y_pred"],
            s=22,
            alpha=0.72,
            color=COLORS[index],
            edgecolors="none",
            label=label,
        )
    ax.plot(limits, limits, linestyle="--", linewidth=1.2, color="black", label="y = x")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_title(f"{MODEL_NAME} Combined Parity", fontsize=TITLE_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_parity_all.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(fig)


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
        linewidth=1.1,
        label=f"Best epoch: {best_epoch}",
    )
    ax.set_xlabel("Epoch", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"RMSE ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_title(f"{MODEL_NAME} Training History", fontsize=TITLE_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_rmse_curve.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(fig)


def save_metrics_table(
    metrics: Dict[str, Dict[str, float]], results: Dict[str, Dict[str, Any]]
) -> None:
    rows = []
    for split_name in ["train", "val", "test"]:
        rows.append(
            {
                "split": split_name,
                "n_samples": len(results[split_name]["y_true"]),
                **metrics[split_name],
            }
        )
    pd.DataFrame(rows).to_csv(
        TABLE_DIR / f"{MODEL_NAME}_metrics.dat",
        sep="\t",
        index=False,
        float_format="%.6f",
    )


@torch.inference_mode()
def predict_cifs(
    model: M3GNetRegressor,
    cif_paths: Sequence[str | Path],
    device: torch.device,
    normalization: Dict[str, float],
    graph_config: Optional[Dict[str, Any]] = None,
) -> pd.DataFrame:
    graph_config = graph_config or {
        "cutoff": CUTOFF,
        "threebody_cutoff": THREEBODY_CUTOFF,
        "max_neighbors": MAX_NEIGHBORS,
    }
    graphs: List[M3GNetData] = []
    valid_names: List[str] = []
    failures: Dict[str, str] = {}
    for path_value in cif_paths:
        path = Path(path_value)
        try:
            structure = Structure.from_file(path)
            graph = structure_to_graph(
                structure,
                target=0.0,
                cif_name=path.name,
                cutoff=float(graph_config["cutoff"]),
                threebody_cutoff=float(graph_config["threebody_cutoff"]),
                max_neighbors=graph_config.get("max_neighbors"),
            )
            graphs.append(graph)
            valid_names.append(path.name)
        except Exception as exc:
            failures[path.name] = str(exc)
    predictions: Dict[str, float] = {}
    if graphs:
        loader = DataLoader(graphs, batch_size=BATCH_SIZE, shuffle=False, num_workers=0)
        model.eval()
        cursor = 0
        for batch in loader:
            batch = batch.to(device)
            with autocast_context(device):
                scaled = model(batch).view(-1)
            values = scaled.float() * float(normalization["std"]) + float(
                normalization["mean"]
            )
            for value in values.cpu().numpy().tolist():
                predictions[valid_names[cursor]] = float(value)
                cursor += 1
    rows = []
    for path_value in cif_paths:
        name = Path(path_value).name
        rows.append(
            {
                "cif_name": name,
                f"predicted_{TARGET_NAME}_{TARGET_UNIT}": predictions.get(name, np.nan),
                "status": "ok" if name in predictions else failures.get(name, "failed"),
            }
        )
    return pd.DataFrame(rows)


def log_configuration(logger: logging.Logger) -> None:
    logger.info("Model name: %s", MODEL_NAME)
    logger.info("Run version: %s", RUN_VERSION)
    logger.info("Output root: %s", OUTPUT_ROOT)
    logger.info("Excel path: %s", EXCEL_PATH)
    logger.info("CIF directory: %s", CIF_DIR)
    logger.info("Excel columns: first=CIF filename, second=target")
    logger.info("Target: %s (%s)", TARGET_NAME, TARGET_UNIT)
    logger.info("Seed: %d", SEED)
    logger.info(
        "Split ratios: train=%.3f, val=%.3f, test=%.3f",
        TRAIN_RATIO,
        VAL_RATIO,
        TEST_RATIO,
    )
    logger.info(
        "Graph configuration: cutoff=%.3f, threebody_cutoff=%.3f, max_neighbors=%s",
        CUTOFF,
        THREEBODY_CUTOFF,
        MAX_NEIGHBORS,
    )
    logger.info(
        "Base model configuration: max_n=%d, max_l=%d, n_blocks=%d, units=%d",
        MAX_N,
        MAX_L,
        N_BLOCKS,
        UNITS,
    )
    logger.info("Loss: MSELoss on standardized training targets")
    logger.info(
        "Optimizer: AdamW(lr=%.6e, weight_decay=%.6e)", LEARNING_RATE, WEIGHT_DECAY
    )
    logger.info(
        "Scheduler: ReduceLROnPlateau(factor=%.3f, patience=%d, min_lr=%.6e)",
        LR_FACTOR,
        LR_PATIENCE,
        MIN_LR,
    )
    logger.info("Early stopping: patience=%d, min_delta=%.6e", PATIENCE, MIN_DELTA)
    logger.info("Optuna enabled: %s", USE_OPTUNA)


def main() -> None:
    total_start = time.perf_counter()
    ensure_directories()
    logger = setup_logger()
    set_global_seed(SEED)
    log_configuration(logger)
    device = get_device(logger)

    data_start = time.perf_counter()
    rows = read_excel_rows(logger)
    logger.info("Valid Excel rows before CIF parsing: %d", len(rows))
    records = load_or_build_graphs(rows, logger)
    train_records, val_records, test_records = split_records(records, logger)
    data_time = time.perf_counter() - data_start
    logger.info(
        "Dataset sizes: train=%d, val=%d, test=%d",
        len(train_records),
        len(val_records),
        len(test_records),
    )
    logger.info(
        "Data preparation time: %.3f s (%s)", data_time, format_duration(data_time)
    )

    selected_parameters: Dict[str, Any] = {}
    if USE_OPTUNA:
        selected_parameters = run_optuna(train_records, val_records, device, logger)

    training_start = time.perf_counter()
    (
        model,
        history,
        best_val_rmse,
        best_epoch,
        early_stopped,
        target_mean,
        target_std,
    ) = fit_model(
        train_records,
        val_records,
        device,
        logger,
        parameters=selected_parameters,
        save_best=True,
        max_epochs=MAX_EPOCHS,
        patience=PATIENCE,
        verbose=True,
    )
    training_time = time.perf_counter() - training_start
    logger.info("Best epoch: %d", best_epoch)
    logger.info("Best validation RMSE: %.6f %s", best_val_rmse, TARGET_UNIT)
    logger.info("Early stopping triggered: %s", early_stopped)
    logger.info(
        "Training time: %.3f s (%s)", training_time, format_duration(training_time)
    )
    logger.info("Best checkpoint: %s", CHECKPOINT_PATH)

    evaluation_start = time.perf_counter()
    model, checkpoint = load_trained_model(CHECKPOINT_PATH, device)
    normalization = checkpoint["normalization"]
    target_mean = float(normalization["mean"])
    target_std = float(normalization["std"])
    evaluation_batch_size = int(
        checkpoint.get("hyperparameters", {}).get("batch_size", BATCH_SIZE)
    )
    loaders = {
        "train": make_loader(train_records, evaluation_batch_size, shuffle=False),
        "val": make_loader(val_records, evaluation_batch_size, shuffle=False),
        "test": make_loader(test_records, evaluation_batch_size, shuffle=False),
    }
    results = {
        split_name: evaluate_model(model, loader, device, target_mean, target_std)
        for split_name, loader in loaders.items()
    }
    metrics = {
        split_name: calculate_metrics(result["y_true"], result["y_pred"])
        for split_name, result in results.items()
    }

    for split_name in ["train", "val", "test"]:
        metric = metrics[split_name]
        logger.info(
            "%s metrics | MAE %.6f %s | RMSE %.6f %s | R2 %.6f",
            split_name.capitalize(),
            metric["mae_eV"],
            TARGET_UNIT,
            metric["rmse_eV"],
            TARGET_UNIT,
            metric["r2"],
        )
        save_prediction_dat(split_name, results[split_name])

    save_all_prediction_dat(results)
    save_metrics_table(metrics, results)
    save_and_plot_history(history, best_epoch)
    limits = parity_limits(results)
    for index, split_name in enumerate(["train", "val", "test"]):
        plot_parity(
            split_name, results[split_name], metrics[split_name], limits, COLORS[index]
        )
    plot_combined_parity(results, metrics, limits)
    evaluation_time = time.perf_counter() - evaluation_start
    logger.info(
        "Evaluation and plotting time: %.3f s (%s)",
        evaluation_time,
        format_duration(evaluation_time),
    )

    parameter_count = sum(parameter.numel() for parameter in model.parameters())
    trainable_count = sum(
        parameter.numel() for parameter in model.parameters() if parameter.requires_grad
    )
    logger.info("Model structure:\n%s", repr(model))
    logger.info("Total parameters: %d", parameter_count)
    logger.info("Trainable parameters: %d", trainable_count)
    total_time = time.perf_counter() - total_start
    logger.info("Total runtime: %.3f s (%s)", total_time, format_duration(total_time))


if __name__ == "__main__":
    try:
        main()
    except torch.cuda.OutOfMemoryError as exc:
        raise RuntimeError(
            "CUDA ran out of memory. Reduce BATCH_SIZE, UNITS, N_BLOCKS, or MAX_NEIGHBORS."
        ) from exc

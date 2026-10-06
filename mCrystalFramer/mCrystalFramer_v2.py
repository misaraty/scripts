from __future__ import annotations

import json
import logging
import math
import os
import random
import re
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
import torch.nn.functional as F
from pymatgen.core import Element, Structure
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, Dataset
from tqdm.auto import tqdm


# User configuration
MODEL_NAME = "CrystalFramer"
RUN_VERSION = "v2"

EXCEL_PATH = "./data.xlsx"
CIF_DIR = "./cif"

TARGET_NAME = "bandgap"
TARGET_UNIT = "eV"

SEED = 42
TRAIN_RATIO = 0.8
VAL_RATIO = 0.1
TEST_RATIO = 0.1

BATCH_SIZE = 8
MAX_EPOCHS = 300
LEARNING_RATE = 5.0e-4
WEIGHT_DECAY = 1.0e-5
PATIENCE = 60
MIN_DELTA = 1.0e-6
NUM_WORKERS = 0
GRAD_CLIP_NORM = 1.0
USE_AMP = True

MODEL_DIM = 128
NUM_HEADS = 8
NUM_LAYERS = 4
FF_DIM = 512
DROPOUT = 0.0
ATOM_FEATURE_DIM = 98
LATTICE_RANGE = 1
SCALE_REAL = 1.4
GAUSS_LOWER_BOUND = 0.5
DISTANCE_BASIS_DIM = 64
DISTANCE_BASIS_MAX = 14.0
DISTANCE_WIDTH_SCALE = 1.0
ANGLE_BASIS_DIM = 48
ANGLE_WIDTH_SCALE = 4.0
DISTANCE_VALUE_COEF = 1.0
ANGLE_VALUE_COEF = 1.0
FRAME_METHOD = "max"
MAX_ATOMS = 320

USE_OPTUNA = False
OPTUNA_N_TRIALS = 20
OPTUNA_TIMEOUT = None
OPTUNA_EPOCHS = 120

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


@dataclass
class ModelConfig:
    atom_feature_dim: int = ATOM_FEATURE_DIM
    model_dim: int = MODEL_DIM
    num_heads: int = NUM_HEADS
    num_layers: int = NUM_LAYERS
    ff_dim: int = FF_DIM
    dropout: float = DROPOUT
    lattice_range: int = LATTICE_RANGE
    scale_real: float = SCALE_REAL
    gauss_lower_bound: float = GAUSS_LOWER_BOUND
    distance_basis_dim: int = DISTANCE_BASIS_DIM
    distance_basis_max: float = DISTANCE_BASIS_MAX
    distance_width_scale: float = DISTANCE_WIDTH_SCALE
    angle_basis_dim: int = ANGLE_BASIS_DIM
    angle_width_scale: float = ANGLE_WIDTH_SCALE
    distance_value_coef: float = DISTANCE_VALUE_COEF
    angle_value_coef: float = ANGLE_VALUE_COEF
    frame_method: str = FRAME_METHOD


@dataclass
class CrystalRecord:
    cif_name: str
    target: float
    atom_features: torch.Tensor
    positions: torch.Tensor
    lattice: torch.Tensor


def ensure_directories() -> None:
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
    logger = logging.getLogger(f"{MODEL_NAME}_{RUN_VERSION}")
    logger.setLevel(logging.INFO)
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
    if hasattr(torch.backends.cuda.matmul, "fp32_precision"):
        torch.backends.cuda.matmul.fp32_precision = "ieee"
        torch.backends.cudnn.conv.fp32_precision = "ieee"
    else:
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False


def seconds_to_hms(seconds: float) -> str:
    seconds = int(round(seconds))
    return f"{seconds // 3600:02d}:{(seconds % 3600) // 60:02d}:{seconds % 60:02d}"


def select_device(logger: logging.Logger) -> torch.device:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info("Device: %s", device)
    if device.type == "cuda":
        props = torch.cuda.get_device_properties(device)
        logger.info("GPU: %s", props.name)
        logger.info("GPU memory: %.2f GiB", props.total_memory / 1024**3)
        logger.info("CUDA runtime: %s", torch.version.cuda)
    return device


def normalize_cif_name(value: Any) -> str:
    if pd.isna(value):
        raise ValueError("empty CIF name")
    if isinstance(value, (int, np.integer)):
        return f"{int(value)}.cif"
    if isinstance(value, (float, np.floating)) and float(value).is_integer():
        return f"{int(value)}.cif"
    name = str(value).strip()
    if not name:
        raise ValueError("empty CIF name")
    if re.fullmatch(r"[+-]?\d+\.0+", name):
        name = str(int(float(name)))
    if not name.lower().endswith(".cif"):
        name += ".cif"
    return name


def site_feature_vector(structure: Structure, feature_dim: int) -> torch.Tensor:
    features = torch.zeros((len(structure), feature_dim), dtype=torch.float32)
    for site_index, site in enumerate(structure):
        total_occupancy = 0.0
        for species, occupancy in site.species.items():
            element = species.element if hasattr(species, "element") else species
            if not isinstance(element, Element):
                element = Element(str(element))
            atomic_number = int(element.Z)
            if atomic_number < 1 or atomic_number > feature_dim:
                raise ValueError(
                    f"atomic number {atomic_number} exceeds feature dimension {feature_dim}"
                )
            occupancy_value = float(occupancy)
            features[site_index, atomic_number - 1] += occupancy_value
            total_occupancy += occupancy_value
        if total_occupancy <= 0.0:
            raise ValueError(f"site {site_index} has no valid occupancy")
    return features


def parse_structure(cif_path: Path, target: float, cif_name: str) -> CrystalRecord:
    structure = Structure.from_file(
        str(cif_path), primitive=False, sort=False, merge_tol=0.0
    )
    if len(structure) < 1:
        raise ValueError("structure contains no sites")
    if len(structure) > MAX_ATOMS:
        raise ValueError(
            f"structure has {len(structure)} atoms, exceeding MAX_ATOMS={MAX_ATOMS}"
        )
    lattice = torch.tensor(np.asarray(structure.lattice.matrix), dtype=torch.float32)
    positions = torch.tensor(np.asarray(structure.cart_coords), dtype=torch.float32)
    if not torch.isfinite(lattice).all() or not torch.isfinite(positions).all():
        raise ValueError("structure contains non-finite coordinates")
    if abs(float(torch.det(lattice))) < 1.0e-8:
        raise ValueError("singular lattice matrix")
    return CrystalRecord(
        cif_name=cif_name,
        target=float(target),
        atom_features=site_feature_vector(structure, ATOM_FEATURE_DIM),
        positions=positions,
        lattice=lattice,
    )


def load_records(
    logger: logging.Logger,
) -> Tuple[List[CrystalRecord], List[Tuple[str, str]]]:
    excel_path = Path(EXCEL_PATH)
    cif_dir = Path(CIF_DIR)
    if not excel_path.is_file():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    if not cif_dir.is_dir():
        raise NotADirectoryError(f"CIF directory not found: {cif_dir}")
    frame = pd.read_excel(excel_path)
    if frame.shape[1] < 2:
        raise ValueError("Excel file must contain at least two columns")

    records: List[CrystalRecord] = []
    failures: List[Tuple[str, str]] = []
    for row_index, row in tqdm(
        frame.iloc[:, :2].iterrows(), total=len(frame), desc="Reading CIF files"
    ):
        raw_name, raw_target = row.iloc[0], row.iloc[1]
        display_name = str(raw_name)
        try:
            cif_name = normalize_cif_name(raw_name)
            target = float(raw_target)
            if not math.isfinite(target):
                raise ValueError("target is NaN or infinite")
            cif_path = cif_dir / cif_name
            if not cif_path.is_file():
                raise FileNotFoundError(f"file not found: {cif_path}")
            records.append(parse_structure(cif_path, target, cif_name))
        except Exception as exc:
            reason = f"row {row_index + 2}: {type(exc).__name__}: {exc}"
            failures.append((display_name, reason))
            logger.warning("Skipped %s | %s", display_name, reason)

    if len(records) < 10:
        raise ValueError(
            f"Only {len(records)} valid samples were found; at least 10 are required"
        )
    logger.info("Excel rows: %d", len(frame))
    logger.info("Valid structures: %d", len(records))
    logger.info("Invalid structures: %d", len(failures))
    return records, failures


def split_records(
    records: Sequence[CrystalRecord], logger: logging.Logger
) -> Dict[str, List[CrystalRecord]]:
    split_path = SPLIT_DIR / f"{MODEL_NAME}_split.csv"
    by_name = {record.cif_name: record for record in records}
    if split_path.is_file():
        saved = pd.read_csv(split_path)
        required = {"cif_name", "target", "split"}
        if required.issubset(saved.columns) and set(saved["cif_name"]) == set(by_name):
            split_map: Dict[str, List[CrystalRecord]] = {
                "train": [],
                "val": [],
                "test": [],
            }
            valid = True
            for row in saved.itertuples(index=False):
                if row.split not in split_map or row.cif_name not in by_name:
                    valid = False
                    break
                record = by_name[row.cif_name]
                if not math.isclose(
                    float(row.target), record.target, rel_tol=0.0, abs_tol=1.0e-10
                ):
                    valid = False
                    break
                split_map[row.split].append(record)
            if valid and all(split_map.values()):
                logger.info("Reused split file: %s", split_path)
                return split_map
        logger.warning("Existing split file is incompatible and will be replaced")

    if not math.isclose(TRAIN_RATIO + VAL_RATIO + TEST_RATIO, 1.0, abs_tol=1.0e-9):
        raise ValueError("TRAIN_RATIO + VAL_RATIO + TEST_RATIO must equal 1")
    indices = np.arange(len(records))
    train_idx, holdout_idx = train_test_split(
        indices, train_size=TRAIN_RATIO, random_state=SEED, shuffle=True
    )
    relative_test = TEST_RATIO / (VAL_RATIO + TEST_RATIO)
    val_idx, test_idx = train_test_split(
        holdout_idx, test_size=relative_test, random_state=SEED, shuffle=True
    )
    split_indices = {"train": train_idx, "val": val_idx, "test": test_idx}
    split_map = {
        key: [records[int(i)] for i in value] for key, value in split_indices.items()
    }
    rows = []
    for split_name in ("train", "val", "test"):
        rows.extend(
            {"cif_name": record.cif_name, "target": record.target, "split": split_name}
            for record in split_map[split_name]
        )
    pd.DataFrame(rows).to_csv(split_path, index=False)
    logger.info("Saved split file: %s", split_path)
    return split_map


class CrystalDataset(Dataset):
    def __init__(
        self, records: Sequence[CrystalRecord], target_mean: float, target_std: float
    ):
        self.records = list(records)
        self.target_mean = float(target_mean)
        self.target_std = float(target_std)

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, index: int) -> Dict[str, Any]:
        record = self.records[index]
        return {
            "cif_name": record.cif_name,
            "atom_features": record.atom_features,
            "positions": record.positions,
            "lattice": record.lattice,
            "target": torch.tensor(
                (record.target - self.target_mean) / self.target_std,
                dtype=torch.float32,
            ),
            "target_raw": torch.tensor(record.target, dtype=torch.float32),
        }


def collate_crystals(samples: Sequence[Dict[str, Any]]) -> Dict[str, Any]:
    return {
        "cif_names": [sample["cif_name"] for sample in samples],
        "graphs": [
            {
                "atom_features": sample["atom_features"],
                "positions": sample["positions"],
                "lattice": sample["lattice"],
            }
            for sample in samples
        ],
        "targets": torch.stack([sample["target"] for sample in samples]),
        "targets_raw": torch.stack([sample["target_raw"] for sample in samples]),
    }


def make_loaders(
    splits: Dict[str, List[CrystalRecord]],
    batch_size: int,
    target_mean: float,
    target_std: float,
) -> Dict[str, DataLoader]:
    loaders: Dict[str, DataLoader] = {}
    for split_name in ("train", "val", "test"):
        dataset = CrystalDataset(splits[split_name], target_mean, target_std)
        generator = torch.Generator().manual_seed(SEED)
        loaders[split_name] = DataLoader(
            dataset,
            batch_size=batch_size,
            shuffle=split_name == "train",
            num_workers=NUM_WORKERS,
            pin_memory=torch.cuda.is_available(),
            persistent_workers=NUM_WORKERS > 0,
            collate_fn=collate_crystals,
            generator=generator,
        )
    return loaders


def move_graphs(
    graphs: Sequence[Dict[str, torch.Tensor]], device: torch.device
) -> List[Dict[str, torch.Tensor]]:
    return [
        {key: value.to(device, non_blocking=True) for key, value in graph.items()}
        for graph in graphs
    ]


def segment_softmax(logits: torch.Tensor, dim: int) -> torch.Tensor:
    return torch.softmax(logits, dim=dim)


class PeriodicGeometry:
    def __init__(
        self, positions: torch.Tensor, lattice: torch.Tensor, lattice_range: int
    ):
        values = torch.arange(
            -lattice_range,
            lattice_range + 1,
            device=positions.device,
            dtype=positions.dtype,
        )
        grid = torch.stack(
            torch.meshgrid(values, values, values, indexing="ij"), dim=-1
        ).reshape(-1, 3)
        translations = grid @ lattice
        vectors = (
            positions[None, :, None, :]
            + translations[None, None, :, :]
            - positions[:, None, None, :]
        )
        self.vectors = vectors
        self.dist2 = vectors.square().sum(dim=-1).clamp_min(1.0e-12)
        self.distances = self.dist2.sqrt()
        self.central_image = int(torch.argmin(grid.square().sum(dim=-1)).item())


def gaussian_basis(
    values: torch.Tensor, count: int, lower: float, upper: float, width_scale: float
) -> torch.Tensor:
    centers = torch.linspace(
        lower, upper, count, device=values.device, dtype=values.dtype
    )
    spacing = (upper - lower) / max(count - 1, 1)
    width = max(spacing * width_scale, 1.0e-6)
    return torch.exp(-0.5 * ((values[..., None] - centers) / width).square())


def build_max_frames(vectors: torch.Tensor, scores: torch.Tensor) -> torch.Tensor:
    local_vectors = vectors.permute(0, 2, 1, 3)
    local_scores = scores.detach().permute(0, 2, 1)
    norms = local_vectors.norm(dim=-1)
    valid = norms > 1.0e-7
    candidate_scores = local_scores.masked_fill(~valid, -torch.inf)
    all_invalid = ~valid.any(dim=-1)
    candidate_scores = torch.where(
        all_invalid[..., None], torch.zeros_like(candidate_scores), candidate_scores
    )

    first_index = candidate_scores.argmax(dim=-1)
    gather_index = first_index[..., None, None].expand(-1, -1, 1, 3)
    primary_raw = torch.gather(local_vectors, 2, gather_index).squeeze(2)
    primary = F.normalize(primary_raw, dim=-1, eps=1.0e-8)
    primary = torch.where(
        all_invalid[..., None],
        torch.tensor([1.0, 0.0, 0.0], device=vectors.device, dtype=vectors.dtype),
        primary,
    )

    unit_vectors = local_vectors / norms[..., None].clamp_min(1.0e-8)
    projection = torch.einsum("ihjc,ihc->ihj", unit_vectors, primary)
    perpendicular = 1.0 - projection.square()
    rank_score = candidate_scores - candidate_scores.amax(dim=-1, keepdim=True)
    rank_score = rank_score + torch.log(perpendicular.clamp_min(1.0e-8))
    rank_score = rank_score.masked_fill(~valid, -torch.inf)
    rank_score = torch.where(
        all_invalid[..., None], torch.zeros_like(rank_score), rank_score
    )
    second_index = rank_score.argmax(dim=-1)
    second_gather = second_index[..., None, None].expand(-1, -1, 1, 3)
    second_vector = torch.gather(unit_vectors, 2, second_gather).squeeze(2)
    secondary_raw = (
        second_vector - (second_vector * primary).sum(dim=-1, keepdim=True) * primary
    )

    axis_index = primary.abs().argmin(dim=-1)
    fallback = F.one_hot(axis_index, num_classes=3).to(primary.dtype)
    fallback = fallback - (fallback * primary).sum(dim=-1, keepdim=True) * primary
    use_fallback = (secondary_raw.norm(dim=-1) < 1.0e-6) | all_invalid
    secondary_raw = torch.where(use_fallback[..., None], fallback, secondary_raw)
    secondary = F.normalize(secondary_raw, dim=-1, eps=1.0e-8)
    tertiary = F.normalize(torch.cross(primary, secondary, dim=-1), dim=-1, eps=1.0e-8)
    secondary = F.normalize(torch.cross(tertiary, primary, dim=-1), dim=-1, eps=1.0e-8)
    return torch.stack([primary, secondary, tertiary], dim=2)


class CrystalFramerAttention(nn.Module):
    def __init__(self, config: ModelConfig):
        super().__init__()
        if config.model_dim % config.num_heads != 0:
            raise ValueError("model_dim must be divisible by num_heads")
        if config.angle_basis_dim % 3 != 0:
            raise ValueError("angle_basis_dim must be divisible by 3")
        self.config = config
        self.num_heads = config.num_heads
        self.head_dim = config.model_dim // config.num_heads
        self.q_proj = nn.Linear(config.model_dim, config.model_dim)
        self.k_proj = nn.Linear(config.model_dim, config.model_dim)
        self.v_proj = nn.Linear(config.model_dim, config.model_dim)
        self.out_proj = nn.Linear(config.model_dim, config.model_dim)
        self.gaussian_selector = nn.Parameter(
            torch.empty(config.num_heads, self.head_dim)
        )
        self.distance_projection = nn.Parameter(
            torch.empty(config.num_heads, config.distance_basis_dim, self.head_dim)
        )
        self.angle_projection = nn.Parameter(
            torch.empty(config.num_heads, config.angle_basis_dim, self.head_dim)
        )
        self.dropout = nn.Dropout(config.dropout)
        self.reset_parameters()

    def reset_parameters(self) -> None:
        for layer in (self.q_proj, self.k_proj, self.v_proj, self.out_proj):
            nn.init.xavier_uniform_(layer.weight)
            nn.init.zeros_(layer.bias)
        nn.init.normal_(self.gaussian_selector, std=self.head_dim**-0.5)
        nn.init.xavier_uniform_(self.distance_projection)
        nn.init.xavier_uniform_(self.angle_projection)

    def forward(self, x: torch.Tensor, geometry: PeriodicGeometry) -> torch.Tensor:
        atom_count = x.shape[0]
        heads = self.num_heads
        head_dim = self.head_dim
        q = self.q_proj(x).view(atom_count, heads, head_dim)
        k = self.k_proj(x).view(atom_count, heads, head_dim)
        v = self.v_proj(x).view(atom_count, heads, head_dim)

        alpha_raw = torch.einsum("ihd,hd->ih", q, self.gaussian_selector)
        alpha = F.softplus(alpha_raw) + self.config.gauss_lower_bound
        image_logits = (
            -0.5
            * alpha[:, None, None, :]
            * geometry.dist2[..., None]
            / self.config.scale_real**2
        )
        positional_bias = torch.logsumexp(image_logits, dim=2)
        pair_logits = torch.einsum("ihd,jhd->ijh", q, k) / math.sqrt(head_dim)
        pair_logits = pair_logits + positional_bias
        attention = segment_softmax(pair_logits, dim=1)
        attention = self.dropout(attention)

        image_weights = segment_softmax(image_logits, dim=2)
        radial = gaussian_basis(
            geometry.distances,
            self.config.distance_basis_dim,
            0.0,
            self.config.distance_basis_max,
            self.config.distance_width_scale,
        )
        radial_mean = torch.einsum("ijrh,ijrk->ijhk", image_weights, radial)
        distance_values = torch.einsum(
            "ijhk,hkd->ijhd", radial_mean, self.distance_projection
        )

        effective_vectors = torch.einsum(
            "ijrh,ijrc->ijhc", image_weights, geometry.vectors
        )
        frames = build_max_frames(effective_vectors, pair_logits).detach()
        periodic_unit_vectors = geometry.vectors / geometry.distances[
            ..., None
        ].clamp_min(1.0e-8)
        basis_per_axis = self.config.angle_basis_dim // 3
        angular_parts = []
        for axis_index in range(3):
            directional_cosines = torch.einsum(
                "ijrc,ihc->ijrh", periodic_unit_vectors, frames[:, :, axis_index, :]
            ).clamp(-1.0, 1.0)
            axis_basis = gaussian_basis(
                directional_cosines,
                basis_per_axis,
                -1.0,
                1.0,
                self.config.angle_width_scale,
            )
            angular_parts.append(
                torch.einsum("ijrh,ijrhk->ijhk", image_weights, axis_basis)
            )
        angular = torch.cat(angular_parts, dim=-1)
        angle_values = torch.einsum("ijhk,hkd->ijhd", angular, self.angle_projection)

        edge_values = (
            v[None, :, :, :]
            + self.config.distance_value_coef * distance_values
            + self.config.angle_value_coef * angle_values
        )
        output = torch.einsum("ijh,ijhd->ihd", attention, edge_values)
        return self.out_proj(output.reshape(atom_count, heads * head_dim))


class CrystalFramerLayer(nn.Module):
    def __init__(self, config: ModelConfig):
        super().__init__()
        self.norm1 = nn.LayerNorm(config.model_dim)
        self.attention = CrystalFramerAttention(config)
        self.norm2 = nn.LayerNorm(config.model_dim)
        self.feed_forward = nn.Sequential(
            nn.Linear(config.model_dim, config.ff_dim),
            nn.GELU(),
            nn.Dropout(config.dropout),
            nn.Linear(config.ff_dim, config.model_dim),
            nn.Dropout(config.dropout),
        )

    def forward(self, x: torch.Tensor, geometry: PeriodicGeometry) -> torch.Tensor:
        x = x + self.attention(self.norm1(x), geometry)
        x = x + self.feed_forward(self.norm2(x))
        return x


class CrystalFramerRegressor(nn.Module):
    def __init__(self, config: ModelConfig):
        super().__init__()
        if config.frame_method != "max":
            raise ValueError(
                "This standalone implementation supports FRAME_METHOD='max' only"
            )
        self.config = config
        self.atom_embedding = nn.Linear(
            config.atom_feature_dim, config.model_dim, bias=False
        )
        self.layers = nn.ModuleList(
            [CrystalFramerLayer(config) for _ in range(config.num_layers)]
        )
        self.final_norm = nn.LayerNorm(config.model_dim)
        self.regression_head = nn.Sequential(
            nn.Linear(config.model_dim, config.model_dim),
            nn.ReLU(),
            nn.Linear(config.model_dim, 1),
        )
        nn.init.normal_(self.atom_embedding.weight, std=config.model_dim**-0.5)

    def forward_graph(self, graph: Dict[str, torch.Tensor]) -> torch.Tensor:
        x = self.atom_embedding(graph["atom_features"])
        geometry = PeriodicGeometry(
            graph["positions"], graph["lattice"], self.config.lattice_range
        )
        for layer in self.layers:
            x = layer(x, geometry)
        pooled = self.final_norm(x).mean(dim=0)
        return self.regression_head(pooled).squeeze(-1)

    def forward(self, graphs: Sequence[Dict[str, torch.Tensor]]) -> torch.Tensor:
        return torch.stack([self.forward_graph(graph) for graph in graphs], dim=0)


def parameter_counts(model: nn.Module) -> Tuple[int, int]:
    total = sum(parameter.numel() for parameter in model.parameters())
    trainable = sum(
        parameter.numel() for parameter in model.parameters() if parameter.requires_grad
    )
    return total, trainable


def autocast_context(device: torch.device):
    enabled = USE_AMP and device.type == "cuda"
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
    squared_error_sum = 0.0
    normalized_loss_sum = 0.0
    sample_count = 0
    for batch in loader:
        targets = batch["targets"].to(device, non_blocking=True)
        targets_raw = batch["targets_raw"].to(device, non_blocking=True)
        graphs = move_graphs(batch["graphs"], device)
        optimizer.zero_grad(set_to_none=True)
        with autocast_context(device):
            predictions = model(graphs).view(-1)
            loss = criterion(predictions, targets.view(-1))
        if not torch.isfinite(loss):
            raise FloatingPointError("non-finite training loss detected")
        loss.backward()
        if GRAD_CLIP_NORM > 0:
            nn.utils.clip_grad_norm_(model.parameters(), GRAD_CLIP_NORM)
        optimizer.step()
        predictions_raw = predictions.detach().float() * target_std + (
            targets_raw - targets * target_std
        )
        squared_error_sum += F.mse_loss(
            predictions_raw, targets_raw, reduction="sum"
        ).item()
        normalized_loss_sum += loss.detach().item() * targets.numel()
        sample_count += targets.numel()
    return (
        math.sqrt(squared_error_sum / sample_count),
        normalized_loss_sum / sample_count,
    )


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
    predicted_values: List[float] = []
    for batch in loader:
        graphs = move_graphs(batch["graphs"], device)
        with autocast_context(device):
            normalized_predictions = model(graphs).view(-1)
        predictions = (
            normalized_predictions.float().cpu().numpy() * target_std + target_mean
        )
        names.extend(batch["cif_names"])
        true_values.extend(batch["targets_raw"].numpy().astype(float).tolist())
        predicted_values.extend(predictions.astype(float).tolist())
    y_true = np.asarray(true_values, dtype=np.float64)
    y_pred = np.asarray(predicted_values, dtype=np.float64)
    metrics = calculate_metrics(y_true, y_pred)
    return {"names": names, "true": y_true, "pred": y_pred, "metrics": metrics}


def calculate_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    return {
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(math.sqrt(mean_squared_error(y_true, y_pred))),
        "r2": float(r2_score(y_true, y_pred)) if len(y_true) > 1 else float("nan"),
    }


def save_checkpoint(
    model: nn.Module,
    optimizer: torch.optim.Optimizer,
    model_config: ModelConfig,
    epoch: int,
    val_rmse: float,
    target_mean: float,
    target_std: float,
) -> None:
    checkpoint = {
        "model_name": MODEL_NAME,
        "run_version": RUN_VERSION,
        "model_state_dict": model.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "model_config": asdict(model_config),
        "best_epoch": int(epoch),
        "best_val_rmse": float(val_rmse),
        "seed": SEED,
        "target_name": TARGET_NAME,
        "target_unit": TARGET_UNIT,
        "target_mean": float(target_mean),
        "target_std": float(target_std),
        "graph_config": {
            "lattice_range": model_config.lattice_range,
            "frame_method": model_config.frame_method,
            "max_atoms": MAX_ATOMS,
        },
        "implementation_note": (
            "Pure PyTorch CrystalFramer-style implementation with finite periodic-image summation; "
            "it does not use the original fused CuPy/CUDA kernels."
        ),
    }
    torch.save(checkpoint, CHECKPOINT_PATH)


def load_trained_model(
    checkpoint_path: str | Path,
    device: Optional[torch.device] = None,
) -> Tuple[CrystalFramerRegressor, Dict[str, Any]]:
    if device is None:
        device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    try:
        checkpoint = torch.load(
            checkpoint_path, map_location=device, weights_only=False
        )
    except TypeError:
        checkpoint = torch.load(checkpoint_path, map_location=device)
    config = ModelConfig(**checkpoint["model_config"])
    model = CrystalFramerRegressor(config).to(device)
    model.load_state_dict(checkpoint["model_state_dict"], strict=True)
    model.eval()
    return model, checkpoint


def train_model(
    splits: Dict[str, List[CrystalRecord]],
    model_config: ModelConfig,
    device: torch.device,
    logger: logging.Logger,
    learning_rate: float = LEARNING_RATE,
    weight_decay: float = WEIGHT_DECAY,
    batch_size: int = BATCH_SIZE,
    max_epochs: int = MAX_EPOCHS,
    save_best: bool = True,
    verbose: bool = True,
) -> Tuple[
    CrystalFramerRegressor, List[Dict[str, float]], int, float, bool, float, float
]:
    train_targets = np.asarray(
        [record.target for record in splits["train"]], dtype=np.float64
    )
    target_mean = float(train_targets.mean())
    target_std = float(train_targets.std())
    if target_std < 1.0e-12:
        target_std = 1.0
    loaders = make_loaders(splits, batch_size, target_mean, target_std)
    model = CrystalFramerRegressor(model_config).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=learning_rate,
        weight_decay=weight_decay,
        betas=(0.9, 0.98),
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer, mode="min", factor=0.5, patience=max(5, PATIENCE // 4), min_lr=1.0e-7
    )
    criterion = nn.MSELoss()
    history: List[Dict[str, float]] = []
    best_val_rmse = float("inf")
    best_epoch = 0
    epochs_without_improvement = 0
    early_stopped = False
    best_state: Optional[Dict[str, torch.Tensor]] = None

    for epoch in range(1, max_epochs + 1):
        train_rmse, train_mse = train_one_epoch(
            model, loaders["train"], optimizer, criterion, device, target_std
        )
        val_result = evaluate_loader(
            model, loaders["val"], device, target_mean, target_std
        )
        val_rmse = val_result["metrics"]["rmse"]
        learning_rate_now = optimizer.param_groups[0]["lr"]
        history.append(
            {
                "epoch": epoch,
                "train_rmse_eV": train_rmse,
                "val_rmse_eV": val_rmse,
                "train_mse_loss": train_mse,
                "learning_rate": learning_rate_now,
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
                    optimizer,
                    model_config,
                    epoch,
                    val_rmse,
                    target_mean,
                    target_std,
                )
        else:
            epochs_without_improvement += 1
        scheduler.step(val_rmse)
        if epochs_without_improvement >= PATIENCE:
            early_stopped = True
            if verbose:
                logger.info("Early stopping triggered at epoch %d", epoch)
            break

    if best_state is None:
        raise RuntimeError("training ended without a valid checkpoint")
    model.load_state_dict(best_state)
    return (
        model,
        history,
        best_epoch,
        best_val_rmse,
        early_stopped,
        target_mean,
        target_std,
    )


def run_optuna(
    splits: Dict[str, List[CrystalRecord]],
    base_config: ModelConfig,
    device: torch.device,
    logger: logging.Logger,
) -> Dict[str, Any]:
    try:
        import optuna
    except ImportError as exc:
        raise ImportError("USE_OPTUNA=True requires optuna") from exc

    def objective(trial: Any) -> float:
        set_global_seed(SEED)
        config = ModelConfig(**asdict(base_config))
        config.dropout = trial.suggest_float("dropout", 0.0, 0.2)
        learning_rate = trial.suggest_float("learning_rate", 1.0e-5, 2.0e-3, log=True)
        weight_decay = trial.suggest_float("weight_decay", 1.0e-7, 1.0e-3, log=True)
        batch_size = trial.suggest_categorical("batch_size", [2, 4, 8, 16])
        _, _, _, best_rmse, _, _, _ = train_model(
            splits,
            config,
            device,
            logger,
            learning_rate=learning_rate,
            weight_decay=weight_decay,
            batch_size=batch_size,
            max_epochs=OPTUNA_EPOCHS,
            save_best=False,
            verbose=False,
        )
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
    logger.info("Optuna best validation RMSE: %.6f %s", study.best_value, TARGET_UNIT)
    logger.info("Optuna best parameters: %s", study.best_params)
    return result


def save_prediction_dat(split_name: str, result: Dict[str, Any]) -> None:
    frame = pd.DataFrame(
        {
            "cif_name": result["names"],
            "true_bandgap_eV": result["true"],
            "predicted_bandgap_eV": result["pred"],
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


def shared_plot_limits(results: Dict[str, Dict[str, Any]]) -> Tuple[float, float]:
    values = np.concatenate(
        [
            np.concatenate([result["true"], result["pred"]])
            for result in results.values()
        ]
    )
    low, high = float(values.min()), float(values.max())
    margin = max((high - low) * 0.06, 0.05)
    return low - margin, high + margin


def style_axes(ax: plt.Axes) -> None:
    ax.tick_params(labelsize=TICK_FONTSIZE)
    if USE_GRID:
        ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)


def plot_parity(
    split_name: str, result: Dict[str, Any], limits: Tuple[float, float], color: str
) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 6.0))
    ax.scatter(
        result["true"], result["pred"], s=34, alpha=0.82, color=color, edgecolors="none"
    )
    ax.plot(limits, limits, "--", color="black", linewidth=1.2, label="y = x")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_title(
        f"{MODEL_NAME}: {split_name.capitalize()} Set", fontsize=TITLE_FONTSIZE
    )
    metrics = result["metrics"]
    annotation = f"MAE = {metrics['mae']:.4f} {TARGET_UNIT}\nRMSE = {metrics['rmse']:.4f} {TARGET_UNIT}\n$R^2$ = {metrics['r2']:.4f}"
    ax.text(
        0.04,
        0.96,
        annotation,
        transform=ax.transAxes,
        va="top",
        fontsize=ANNOTATION_FONTSIZE,
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
    fig, ax = plt.subplots(figsize=(6.4, 6.1))
    for split_name, color in zip(("train", "val", "test"), COLORS[:3]):
        result = results[split_name]
        ax.scatter(
            result["true"],
            result["pred"],
            s=30,
            alpha=0.76,
            color=color,
            edgecolors="none",
            label=f"{split_name.capitalize()} ($R^2$={result['metrics']['r2']:.3f})",
        )
    ax.plot(limits, limits, "--", color="black", linewidth=1.2, label="y = x")
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_title(f"{MODEL_NAME}: All Splits", fontsize=TITLE_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_parity_all.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(fig)

    frames = []
    for split_name in ("train", "val", "test"):
        result = results[split_name]
        frame = pd.DataFrame(
            {
                "split": split_name,
                "cif_name": result["names"],
                "true_bandgap_eV": result["true"],
                "predicted_bandgap_eV": result["pred"],
            }
        )
        frame["error_eV"] = frame["predicted_bandgap_eV"] - frame["true_bandgap_eV"]
        frame["absolute_error_eV"] = frame["error_eV"].abs()
        frames.append(frame)
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
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    ax.plot(frame["epoch"], frame["train_rmse_eV"], color=COLORS[0], label="Train RMSE")
    ax.plot(
        frame["epoch"], frame["val_rmse_eV"], color=COLORS[1], label="Validation RMSE"
    )
    best_row = frame.loc[frame["epoch"] == best_epoch].iloc[0]
    ax.scatter(
        [best_epoch],
        [best_row["val_rmse_eV"]],
        color=COLORS[3],
        s=55,
        zorder=4,
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


def save_metrics(results: Dict[str, Dict[str, Any]]) -> None:
    rows = []
    for split_name in ("train", "val", "test"):
        metrics = results[split_name]["metrics"]
        rows.append(
            {
                "split": split_name,
                "n_samples": len(results[split_name]["true"]),
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


@torch.inference_mode()
def predict_cifs(
    model: CrystalFramerRegressor,
    cif_paths: Sequence[str | Path],
    device: torch.device,
    target_mean: float,
    target_std: float,
) -> pd.DataFrame:
    rows = []
    model.eval()
    for value in cif_paths:
        path = Path(value)
        record = parse_structure(path, 0.0, path.name)
        graph = {
            "atom_features": record.atom_features.to(device),
            "positions": record.positions.to(device),
            "lattice": record.lattice.to(device),
        }
        with autocast_context(device):
            normalized = model([graph]).view(-1)[0]
        prediction = float(normalized.float().cpu()) * target_std + target_mean
        rows.append(
            {
                "cif_name": path.name,
                f"predicted_{TARGET_NAME}_{TARGET_UNIT}": prediction,
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    total_start = time.perf_counter()
    ensure_directories()
    logger = setup_logger()
    set_global_seed(SEED)
    device = select_device(logger)
    logger.info("Model: %s", MODEL_NAME)
    logger.info("Run version: %s", RUN_VERSION)
    logger.info("Excel path: %s", EXCEL_PATH)
    logger.info("CIF directory: %s", CIF_DIR)
    logger.info("Excel columns: first column=CIF name, second column=target")
    logger.info("Target: %s (%s)", TARGET_NAME, TARGET_UNIT)
    logger.info("Seed: %d", SEED)
    logger.info("Split ratios: %.3f / %.3f / %.3f", TRAIN_RATIO, VAL_RATIO, TEST_RATIO)
    logger.info("Loss: MSELoss")
    logger.info("Pure PyTorch path: enabled")
    logger.info("Fused CuPy/CUDA kernels: disabled")
    logger.info(
        "The pure PyTorch path can be slower than the original fused-kernel implementation"
    )

    data_start = time.perf_counter()
    records, failures = load_records(logger)
    splits = split_records(records, logger)
    data_time = time.perf_counter() - data_start
    for split_name in ("train", "val", "test"):
        logger.info("%s samples: %d", split_name.capitalize(), len(splits[split_name]))

    model_config = ModelConfig()
    selected_lr = LEARNING_RATE
    selected_weight_decay = WEIGHT_DECAY
    selected_batch_size = BATCH_SIZE
    logger.info("Optuna enabled: %s", USE_OPTUNA)
    if USE_OPTUNA:
        optuna_result = run_optuna(splits, model_config, device, logger)
        best_params = optuna_result["best_params"]
        model_config.dropout = float(best_params.get("dropout", model_config.dropout))
        selected_lr = float(best_params.get("learning_rate", selected_lr))
        selected_weight_decay = float(
            best_params.get("weight_decay", selected_weight_decay)
        )
        selected_batch_size = int(best_params.get("batch_size", selected_batch_size))

    probe_model = CrystalFramerRegressor(model_config)
    total_parameters, trainable_parameters = parameter_counts(probe_model)
    logger.info(
        "Model configuration: %s", json.dumps(asdict(model_config), sort_keys=True)
    )
    logger.info("Model architecture:\n%s", probe_model)
    logger.info("Total parameters: %d", total_parameters)
    logger.info("Trainable parameters: %d", trainable_parameters)
    logger.info(
        "Optimizer: AdamW | lr %.6e | weight_decay %.6e | betas (0.9, 0.98)",
        selected_lr,
        selected_weight_decay,
    )
    logger.info(
        "Scheduler: ReduceLROnPlateau | factor 0.5 | patience %d | min_lr 1e-7",
        max(5, PATIENCE // 4),
    )
    logger.info(
        "Batch size: %d | maximum epochs: %d | early-stop patience: %d",
        selected_batch_size,
        MAX_EPOCHS,
        PATIENCE,
    )
    del probe_model

    training_start = time.perf_counter()
    (
        model,
        history,
        best_epoch,
        best_val_rmse,
        early_stopped,
        target_mean,
        target_std,
    ) = train_model(
        splits,
        model_config,
        device,
        logger,
        learning_rate=selected_lr,
        weight_decay=selected_weight_decay,
        batch_size=selected_batch_size,
        max_epochs=MAX_EPOCHS,
        save_best=True,
        verbose=True,
    )
    training_time = time.perf_counter() - training_start
    logger.info("Best epoch: %d", best_epoch)
    logger.info("Best validation RMSE: %.6f %s", best_val_rmse, TARGET_UNIT)
    logger.info("Early stopping triggered: %s", early_stopped)
    logger.info(
        "Training time: %.3f s (%s)", training_time, seconds_to_hms(training_time)
    )
    logger.info("Best checkpoint: %s", CHECKPOINT_PATH)

    evaluation_start = time.perf_counter()
    model, checkpoint = load_trained_model(CHECKPOINT_PATH, device)
    loaders = make_loaders(
        splits, selected_batch_size, checkpoint["target_mean"], checkpoint["target_std"]
    )
    results = {
        split_name: evaluate_loader(
            model,
            loaders[split_name],
            device,
            checkpoint["target_mean"],
            checkpoint["target_std"],
        )
        for split_name in ("train", "val", "test")
    }
    for split_name in ("train", "val", "test"):
        metrics = results[split_name]["metrics"]
        logger.info(
            "%s | MAE %.6f %s | RMSE %.6f %s | R2 %.6f",
            split_name.capitalize(),
            metrics["mae"],
            TARGET_UNIT,
            metrics["rmse"],
            TARGET_UNIT,
            metrics["r2"],
        )
        save_prediction_dat(split_name, results[split_name])
    limits = shared_plot_limits(results)
    for split_name, color in zip(("train", "val", "test"), COLORS[:3]):
        plot_parity(split_name, results[split_name], limits, color)
    plot_combined_parity(results, limits)
    plot_rmse_curve(history, best_epoch)
    save_metrics(results)
    evaluation_time = time.perf_counter() - evaluation_start

    logger.info("Invalid sample count: %d", len(failures))
    logger.info("Data reading time: %.3f s (%s)", data_time, seconds_to_hms(data_time))
    logger.info(
        "Evaluation and plotting time: %.3f s (%s)",
        evaluation_time,
        seconds_to_hms(evaluation_time),
    )
    total_time = time.perf_counter() - total_start
    logger.info("Total runtime: %.3f s (%s)", total_time, seconds_to_hms(total_time))


if __name__ == "__main__":
    main()

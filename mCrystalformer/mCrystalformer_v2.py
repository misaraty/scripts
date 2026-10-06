import copy
import gc
import json
import logging
import math
import os
import random
import sys
import time
from contextlib import nullcontext
from collections import Counter
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
from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from torch.utils.data import DataLoader, Dataset
from tqdm.auto import tqdm


# ==================== USER CONFIGURATION ====================

MODEL_NAME = "Crystalformer"
RUN_VERSION = "v2"

EXCEL_PATH = "./data.xlsx"
CIF_DIR = "./cif"

TARGET_NAME = "bandgap"
TARGET_UNIT = "eV"

SEED = 42
TRAIN_RATIO = 0.80
VAL_RATIO = 0.10
TEST_RATIO = 0.10

BATCH_SIZE = 16
MAX_EPOCHS = 300
LEARNING_RATE = 5.0e-4
WEIGHT_DECAY = 1.0e-5
ADAM_BETAS = (0.9, 0.98)
LR_DECAY_STEPS = 4000.0
PATIENCE = 60
MIN_DELTA = 1.0e-5
GRAD_CLIP_NORM = 1.0
NUM_WORKERS = 0
USE_AMP = False

USE_OPTUNA = False
OPTUNA_N_TRIALS = 30
OPTUNA_TIMEOUT = None
OPTUNA_MAX_EPOCHS = 120

MODEL_DIM = 128
NUM_LAYERS = 4
HEAD_NUM = 8
FF_DIM = 512
DROPOUT = 0.0
DOMAIN = "real-reci"
LATTICE_RANGE = 2
SCALE_REAL = 1.4
SCALE_RECI = 2.2
GAUSS_LB_REAL = 0.5
GAUSS_LB_RECI = 0.5
VALUE_PE_DIST_REAL = 64
VALUE_PE_DIST_MAX = 14.0
VALUE_PE_WIDTH_SCALE = 1.0
EXCLUDE_SELF = False
USE_T_FIXUP = True
POOLING = "average"
MAX_ATOMIC_NUMBER = 118
CELL_FORMAT = "primitive"

FIG_DPI = 600
TITLE_FONTSIZE = 15
LABEL_FONTSIZE = 13
TICK_FONTSIZE = 11
LEGEND_FONTSIZE = 11
ANNOTATION_FONTSIZE = 10
USE_GRID = True
COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:red"]

# ============================================================


SCRIPT_DIR = Path(os.path.realpath(__file__)).parent
os.chdir(SCRIPT_DIR)
OUTPUT_ROOT = SCRIPT_DIR / f"{MODEL_NAME}_{RUN_VERSION}"
FIGURE_DIR = OUTPUT_ROOT / "figure"
DAT_DIR = OUTPUT_ROOT / "dat"
TABLE_DIR = OUTPUT_ROOT / "table"
LOG_DIR = OUTPUT_ROOT / "log"
SPLIT_DIR = OUTPUT_ROOT / "split"
CACHE_DIR = OUTPUT_ROOT / "cache"
CHECKPOINT_PATH = OUTPUT_ROOT / f"{MODEL_NAME}_best.pt"
SPLIT_PATH = SPLIT_DIR / f"{MODEL_NAME}_split.csv"
CACHE_PATH = CACHE_DIR / f"{MODEL_NAME}_structures.pt"


@dataclass
class ModelConfig:
    model_dim: int = MODEL_DIM
    num_layers: int = NUM_LAYERS
    head_num: int = HEAD_NUM
    ff_dim: int = FF_DIM
    dropout: float = DROPOUT
    domain: str = DOMAIN
    lattice_range: int = LATTICE_RANGE
    scale_real: float = SCALE_REAL
    scale_reci: float = SCALE_RECI
    gauss_lb_real: float = GAUSS_LB_REAL
    gauss_lb_reci: float = GAUSS_LB_RECI
    value_pe_dist_real: int = VALUE_PE_DIST_REAL
    value_pe_dist_max: float = VALUE_PE_DIST_MAX
    value_pe_width_scale: float = VALUE_PE_WIDTH_SCALE
    exclude_self: bool = EXCLUDE_SELF
    use_t_fixup: bool = USE_T_FIXUP
    pooling: str = POOLING
    max_atomic_number: int = MAX_ATOMIC_NUMBER

    def validate(self) -> None:
        if self.model_dim % self.head_num != 0:
            raise ValueError("MODEL_DIM must be divisible by HEAD_NUM.")
        if self.domain not in {"real", "reci", "real-reci", "reci-real", "multihead"}:
            raise ValueError(f"Unsupported DOMAIN: {self.domain}")
        if self.domain == "multihead" and self.head_num % 2 != 0:
            raise ValueError("HEAD_NUM must be even when DOMAIN='multihead'.")
        if self.lattice_range < 0:
            raise ValueError("LATTICE_RANGE must be non-negative.")
        if self.value_pe_dist_real < 0:
            raise ValueError("VALUE_PE_DIST_REAL must be non-negative.")
        if self.pooling not in {"average", "max"}:
            raise ValueError("POOLING must be 'average' or 'max'.")


@dataclass
class StructureRecord:
    cif_name: str
    target: float
    x: torch.Tensor
    pos: torch.Tensor
    lattice: torch.Tensor


class BandgapDataset(Dataset):
    def __init__(self, records: Sequence[StructureRecord]):
        self.records = list(records)

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, index: int) -> StructureRecord:
        return self.records[index]


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


def set_global_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.benchmark = False
    torch.backends.cudnn.deterministic = True
    if hasattr(torch.backends, "cuda") and hasattr(torch.backends.cuda, "matmul"):
        torch.backends.cuda.matmul.allow_tf32 = False
    torch.backends.cudnn.allow_tf32 = False


def seed_worker(worker_id: int) -> None:
    worker_seed = (SEED + worker_id) % (2**32)
    np.random.seed(worker_seed)
    random.seed(worker_seed)


def setup_logger() -> logging.Logger:
    logger = logging.getLogger(MODEL_NAME)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    file_handler = logging.FileHandler(
        LOG_DIR / f"{MODEL_NAME}_training.log", mode="w", encoding="utf-8"
    )
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def format_duration(seconds: float) -> str:
    seconds_int = max(0, int(round(seconds)))
    hours, remainder = divmod(seconds_int, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{seconds:.3f} s ({hours:02d}:{minutes:02d}:{secs:02d})"


def resolve_path(path_string: str) -> Path:
    path = Path(path_string).expanduser()
    return path if path.is_absolute() else SCRIPT_DIR / path


def read_excel_rows(logger: logging.Logger) -> List[Tuple[str, float]]:
    excel_path = resolve_path(EXCEL_PATH)
    if not excel_path.is_file():
        raise FileNotFoundError(f"Excel file not found: {excel_path}")
    frame = pd.read_excel(excel_path)
    if frame.shape[1] < 2:
        raise ValueError("The Excel file must contain at least two columns.")
    rows: List[Tuple[str, float]] = []
    rejected = 0
    for row_index, row in frame.iloc[:, :2].iterrows():
        raw_name, raw_target = row.iloc[0], row.iloc[1]
        if pd.isna(raw_name) or str(raw_name).strip() == "":
            logger.warning("Skipped Excel row %d: empty CIF name.", row_index + 2)
            rejected += 1
            continue
        try:
            target = float(raw_target)
        except (TypeError, ValueError):
            logger.warning(
                "Skipped %s: target is not numeric (%r).", raw_name, raw_target
            )
            rejected += 1
            continue
        if not math.isfinite(target):
            logger.warning(
                "Skipped %s: target is not finite (%r).", raw_name, raw_target
            )
            rejected += 1
            continue
        rows.append((f"{int(raw_name)}.cif", target))
    if not rows:
        raise ValueError(
            "No valid CIF names and targets were found in the first two Excel columns."
        )
    names = [name for name, _ in rows]
    duplicates = sorted(name for name, count in Counter(names).items() if count > 1)
    if duplicates:
        raise ValueError(f"Duplicate CIF names are not allowed: {duplicates[:10]}")
    logger.info(
        "Excel rows accepted: %d; rejected before CIF parsing: %d", len(rows), rejected
    )
    return rows


def site_species_vector(structure: Structure, max_atomic_number: int) -> torch.Tensor:
    features = torch.zeros((len(structure), max_atomic_number), dtype=torch.float32)
    for site_index, site in enumerate(structure):
        for specie, occupancy in site.species.items():
            atomic_number = int(specie.Z)
            if atomic_number < 1 or atomic_number > max_atomic_number:
                raise ValueError(
                    f"Atomic number {atomic_number} exceeds MAX_ATOMIC_NUMBER={max_atomic_number}."
                )
            features[site_index, atomic_number - 1] += float(occupancy)
    if not torch.isfinite(features).all() or torch.any(features.sum(dim=1) <= 0):
        raise ValueError("Invalid site occupancy vectors were generated.")
    return features


def standardize_structure(structure: Structure) -> Structure:
    mode = CELL_FORMAT.lower().strip()
    if mode == "raw":
        return structure
    analyzer = SpacegroupAnalyzer(structure)
    if mode == "primitive":
        return analyzer.get_primitive_standard_structure()
    if mode == "conventional":
        return analyzer.get_conventional_standard_structure()
    raise ValueError("CELL_FORMAT must be 'raw', 'primitive', or 'conventional'.")


def parse_cif(
    cif_path: Path,
    target: float,
    model_config: ModelConfig,
    record_name: Optional[str] = None,
) -> StructureRecord:
    structure = Structure.from_file(str(cif_path))
    structure = standardize_structure(structure)
    if len(structure) == 0:
        raise ValueError("The structure contains no atomic sites.")
    x = site_species_vector(structure, model_config.max_atomic_number)
    pos = torch.as_tensor(np.asarray(structure.cart_coords), dtype=torch.float32)
    lattice = torch.as_tensor(np.asarray(structure.lattice.matrix), dtype=torch.float32)
    if not torch.isfinite(pos).all() or not torch.isfinite(lattice).all():
        raise ValueError(
            "Non-finite Cartesian coordinates or lattice values were found."
        )
    volume = abs(float(torch.det(lattice).item()))
    if volume <= 1.0e-8:
        raise ValueError("The lattice volume is zero or numerically singular.")
    return StructureRecord(record_name or cif_path.name, float(target), x, pos, lattice)


def _torch_load(path: Path, map_location: Any = "cpu") -> Any:
    try:
        return torch.load(path, map_location=map_location, weights_only=False)
    except TypeError:
        return torch.load(path, map_location=map_location)


def cache_signature(
    rows: Sequence[Tuple[str, float]], model_config: ModelConfig
) -> Dict[str, Any]:
    cif_dir = resolve_path(CIF_DIR)
    items = []
    for name, target in rows:
        path = cif_dir / name
        stat = path.stat() if path.is_file() else None
        items.append(
            (
                name,
                target,
                None if stat is None else stat.st_size,
                None if stat is None else stat.st_mtime_ns,
            )
        )
    return {
        "items": items,
        "cell_format": CELL_FORMAT,
        "max_atomic_number": model_config.max_atomic_number,
    }


def load_or_build_records(
    rows: Sequence[Tuple[str, float]], model_config: ModelConfig, logger: logging.Logger
) -> List[StructureRecord]:
    signature = cache_signature(rows, model_config)
    if CACHE_PATH.is_file():
        try:
            cached = _torch_load(CACHE_PATH)
            if cached.get("signature") == signature:
                records = [StructureRecord(**item) for item in cached["records"]]
                logger.info(
                    "Loaded %d parsed structures from %s", len(records), CACHE_PATH
                )
                return records
            logger.info("Structure cache is stale and will be rebuilt.")
        except Exception as exc:
            logger.warning("Could not load structure cache: %s", exc)

    cif_dir = resolve_path(CIF_DIR)
    if not cif_dir.is_dir():
        raise FileNotFoundError(f"CIF directory not found: {cif_dir}")
    records: List[StructureRecord] = []
    rejected: List[Tuple[str, str]] = []
    for cif_name, target in tqdm(rows, desc="Parsing CIF files"):
        cif_path = cif_dir / cif_name
        if not cif_path.is_file():
            reason = f"file not found: {cif_path}"
            logger.warning("Skipped %s: %s", cif_name, reason)
            rejected.append((cif_name, reason))
            continue
        try:
            records.append(
                parse_cif(cif_path, target, model_config, record_name=cif_name)
            )
        except Exception as exc:
            reason = f"{type(exc).__name__}: {exc}"
            logger.warning("Skipped %s: %s", cif_name, reason)
            rejected.append((cif_name, reason))
    if len(records) < 20:
        raise ValueError(
            f"Only {len(records)} valid structures remain. At least 20 are required for an 80/10/10 split."
        )
    serializable = [
        {
            "cif_name": r.cif_name,
            "target": r.target,
            "x": r.x,
            "pos": r.pos,
            "lattice": r.lattice,
        }
        for r in records
    ]
    torch.save(
        {"signature": signature, "records": serializable, "rejected": rejected},
        CACHE_PATH,
    )
    logger.info(
        "Parsed structures accepted: %d; rejected: %d", len(records), len(rejected)
    )
    return records


def create_or_load_split(
    records: Sequence[StructureRecord], logger: logging.Logger
) -> Dict[str, str]:
    record_map = {record.cif_name: record for record in records}
    if SPLIT_PATH.is_file():
        split_frame = pd.read_csv(SPLIT_PATH)
        required = {"cif_name", "target", "split"}
        if required.issubset(split_frame.columns):
            names = split_frame["cif_name"].astype(str).tolist()
            split_values = set(split_frame["split"].astype(str))
            if (
                set(names) == set(record_map)
                and len(names) == len(record_map)
                and split_values <= {"train", "val", "test"}
            ):
                mapping = dict(zip(names, split_frame["split"].astype(str)))
                if all(
                    any(value == split for value in mapping.values())
                    for split in ("train", "val", "test")
                ):
                    logger.info("Reused fixed split from %s", SPLIT_PATH)
                    return mapping
        logger.warning(
            "Existing split is incompatible with the current valid dataset and will be replaced."
        )

    names = np.asarray([record.cif_name for record in records])
    train_names, remainder = train_test_split(
        names, train_size=TRAIN_RATIO, random_state=SEED, shuffle=True
    )
    relative_val = VAL_RATIO / (VAL_RATIO + TEST_RATIO)
    val_names, test_names = train_test_split(
        remainder, train_size=relative_val, random_state=SEED, shuffle=True
    )
    mapping = {name: "train" for name in train_names}
    mapping.update({name: "val" for name in val_names})
    mapping.update({name: "test" for name in test_names})
    frame = pd.DataFrame(
        {
            "cif_name": sorted(mapping),
            "target": [record_map[name].target for name in sorted(mapping)],
            "split": [mapping[name] for name in sorted(mapping)],
        }
    )
    frame.to_csv(SPLIT_PATH, index=False)
    logger.info("Created fixed split at %s", SPLIT_PATH)
    return mapping


def collate_records(records: Sequence[StructureRecord]) -> Dict[str, Any]:
    sizes = torch.tensor([record.x.shape[0] for record in records], dtype=torch.long)
    batch = torch.repeat_interleave(torch.arange(len(records), dtype=torch.long), sizes)
    return {
        "x": torch.cat([record.x for record in records], dim=0),
        "pos": torch.cat([record.pos for record in records], dim=0),
        "lattice": torch.stack([record.lattice for record in records], dim=0),
        "batch": batch,
        "sizes": sizes,
        "target": torch.tensor(
            [record.target for record in records], dtype=torch.float32
        ),
        "cif_name": [record.cif_name for record in records],
    }


def make_loader(dataset: Dataset, batch_size: int, shuffle: bool) -> DataLoader:
    generator = torch.Generator()
    generator.manual_seed(SEED)
    return DataLoader(
        dataset,
        batch_size=batch_size,
        shuffle=shuffle,
        num_workers=NUM_WORKERS,
        collate_fn=collate_records,
        worker_init_fn=seed_worker,
        generator=generator,
        pin_memory=torch.cuda.is_available(),
        persistent_workers=NUM_WORKERS > 0,
    )


def move_batch(batch: Dict[str, Any], device: torch.device) -> Dict[str, Any]:
    return {
        key: value.to(device, non_blocking=True) if torch.is_tensor(value) else value
        for key, value in batch.items()
    }


def lattice_grid(radius: int, device: torch.device, dtype: torch.dtype) -> torch.Tensor:
    values = torch.arange(-radius, radius + 1, device=device, dtype=dtype)
    return torch.cartesian_prod(values, values, values)


def reciprocal_half_grid(
    radius: int, device: torch.device, dtype: torch.dtype
) -> Tuple[torch.Tensor, torch.Tensor]:
    full = lattice_grid(radius, device, dtype)
    integer_full = full.to(torch.long)
    mask = (
        (integer_full[:, 0] > 0)
        | ((integer_full[:, 0] == 0) & (integer_full[:, 1] > 0))
        | (
            (integer_full[:, 0] == 0)
            & (integer_full[:, 1] == 0)
            & (integer_full[:, 2] >= 0)
        )
    )
    half = full[mask]
    multiplicity = torch.where(
        torch.all(half == 0, dim=1),
        torch.ones(half.shape[0], device=device, dtype=dtype),
        torch.full((half.shape[0],), 2.0, device=device, dtype=dtype),
    )
    return half, multiplicity


def positive_elu(
    value: torch.Tensor, lower_bound: float, beta: float = 0.1
) -> torch.Tensor:
    return (1.0 - lower_bound) * F.elu(value * (beta / (1.0 - lower_bound))) + 1.0


class PeriodicMultiheadAttention(nn.Module):
    def __init__(self, config: ModelConfig, domain: str):
        super().__init__()
        self.config = copy.deepcopy(config)
        self.domain = domain
        self.num_heads = config.head_num
        self.head_dim = config.model_dim // config.head_num
        self.q_proj = nn.Linear(config.model_dim, config.model_dim)
        self.k_proj = nn.Linear(config.model_dim, config.model_dim)
        self.v_proj = nn.Linear(config.model_dim, config.model_dim)
        self.out_proj = nn.Linear(config.model_dim, config.model_dim)
        self.gaussian_weight = nn.Parameter(torch.empty(config.head_num, self.head_dim))
        self.dropout = nn.Dropout(config.dropout)
        self.register_buffer("gaussian_scale", torch.tensor(0.0))
        self.register_buffer("gaussian_shift", torch.tensor(0.0))

        if config.value_pe_dist_real > 0:
            centers = torch.linspace(
                0.0, config.value_pe_dist_max, config.value_pe_dist_real
            )
            self.register_buffer("rbf_centers", centers)
            spacing = config.value_pe_dist_max / max(1, config.value_pe_dist_real - 1)
            self.rbf_gamma = 1.0 / max(
                (spacing * config.value_pe_width_scale) ** 2, 1.0e-8
            )
            self.value_pe_proj = nn.Parameter(
                torch.empty(config.head_num, config.value_pe_dist_real, self.head_dim)
            )
        else:
            self.register_buffer("rbf_centers", torch.empty(0))
            self.rbf_gamma = 1.0
            self.value_pe_proj = None
        self.reset_parameters()

    def reset_parameters(self) -> None:
        for layer in (self.q_proj, self.k_proj, self.v_proj, self.out_proj):
            nn.init.xavier_uniform_(layer.weight)
            if layer.bias is not None:
                nn.init.zeros_(layer.bias)
        nn.init.normal_(self.gaussian_weight, std=self.head_dim**-0.5)
        if self.value_pe_proj is not None:
            nn.init.xavier_uniform_(
                self.value_pe_proj.view(self.num_heads * self.head_dim, -1)
            )

    def normalized_alpha(self, q: torch.Tensor) -> torch.Tensor:
        raw = torch.einsum("ihd,hd->ih", q, self.gaussian_weight)
        if self.training and float(self.gaussian_scale.item()) <= 0.0:
            with torch.no_grad():
                scale = raw.detach().std(unbiased=False).clamp_min(1.0e-6).reciprocal()
                shift = -raw.detach().mean()
                self.gaussian_scale.copy_(scale)
                self.gaussian_shift.copy_(shift)
        scale = self.gaussian_scale.clamp_min(1.0e-6)
        return (raw + self.gaussian_shift) * scale

    def real_periodic_terms(
        self, pos: torch.Tensor, lattice: torch.Tensor, alpha_raw: torch.Tensor
    ) -> Tuple[torch.Tensor, Optional[torch.Tensor]]:
        grid = lattice_grid(self.config.lattice_range, pos.device, pos.dtype)
        translations = grid @ lattice
        delta = pos.unsqueeze(0) - pos.unsqueeze(1)
        images = delta.unsqueeze(2) + translations.view(1, 1, -1, 3)
        dist2 = torch.sum(images * images, dim=-1)
        if self.config.exclude_self:
            center_index = grid.shape[0] // 2
            atom_index = torch.arange(pos.shape[0], device=pos.device)
            dist2[atom_index, atom_index, center_index] = torch.inf
        rho = positive_elu(alpha_raw, self.config.gauss_lb_real)
        alpha = -0.5 * rho * (self.config.scale_real**-2)
        exponent = alpha[:, :, None, None] * dist2[:, None, :, :]
        log_kernel = torch.logsumexp(exponent, dim=-1)

        value_encoding = None
        if self.value_pe_proj is not None:
            image_weights = torch.softmax(exponent, dim=-1)
            distances = torch.sqrt(dist2.clamp_min(0.0))
            rbf = torch.exp(
                -self.rbf_gamma
                * (distances.unsqueeze(-1) - self.rbf_centers.view(1, 1, 1, -1)) ** 2
            )
            expected_rbf = torch.einsum("ihjm,ijmk->ihjk", image_weights, rbf)
            value_projection = self.value_pe_proj[: alpha_raw.shape[1]]
            value_encoding = torch.einsum(
                "ihjk,hkd->ihjd", expected_rbf, value_projection
            )
        return log_kernel, value_encoding

    def reciprocal_periodic_terms(
        self, pos: torch.Tensor, lattice: torch.Tensor, alpha_raw: torch.Tensor
    ) -> torch.Tensor:
        volume = torch.abs(torch.det(lattice)).clamp_min(1.0e-8)
        reciprocal = 2.0 * math.pi * torch.linalg.inv(lattice).transpose(0, 1)
        half_grid, multiplicity = reciprocal_half_grid(
            self.config.lattice_range, pos.device, pos.dtype
        )
        k_vectors = half_grid @ reciprocal
        k2 = torch.sum(k_vectors * k_vectors, dim=-1)
        delta = pos.unsqueeze(0) - pos.unsqueeze(1)
        cos_kr = torch.cos(torch.einsum("ijc,kc->ijk", delta, k_vectors))
        cos_kr = cos_kr * multiplicity.view(1, 1, -1)
        rho = positive_elu(alpha_raw, self.config.gauss_lb_reci)
        alpha = -0.5 * rho * (self.config.scale_reci**2)
        spectral = torch.exp(alpha[:, :, None] * k2.view(1, 1, -1))
        kernel = torch.einsum("ihk,ijk->ihj", spectral, cos_kr).clamp_min(1.0e-8)
        correction = 1.5 * torch.log(
            (-4.0 * math.pi * alpha).clamp_min(1.0e-8)
        ) - torch.log(volume)
        return torch.log(kernel) + correction[:, :, None]

    def forward_one(
        self,
        q: torch.Tensor,
        k: torch.Tensor,
        v: torch.Tensor,
        alpha_raw: torch.Tensor,
        pos: torch.Tensor,
        lattice: torch.Tensor,
    ) -> torch.Tensor:
        atom_count = q.shape[0]
        dot_logits = torch.einsum("ihd,jhd->ihj", q, k) / math.sqrt(self.head_dim)

        if self.domain == "real":
            periodic_logits, value_encoding = self.real_periodic_terms(
                pos, lattice, alpha_raw
            )
        elif self.domain == "reci":
            periodic_logits = self.reciprocal_periodic_terms(pos, lattice, alpha_raw)
            value_encoding = None
        elif self.domain == "multihead":
            if self.num_heads % 2 != 0:
                raise RuntimeError(
                    "Multihead real/reciprocal mode requires an even number of heads."
                )
            split = self.num_heads // 2
            real_logits, real_values = self.real_periodic_terms(
                pos, lattice, alpha_raw[:, :split]
            )
            reci_logits = self.reciprocal_periodic_terms(
                pos, lattice, alpha_raw[:, split:]
            )
            periodic_logits = torch.cat((real_logits, reci_logits), dim=1)
            if real_values is None:
                value_encoding = None
            else:
                zeros = torch.zeros(
                    atom_count,
                    self.num_heads - split,
                    atom_count,
                    self.head_dim,
                    device=q.device,
                    dtype=q.dtype,
                )
                value_encoding = torch.cat((real_values, zeros), dim=1)
        else:
            raise RuntimeError(f"Unexpected attention domain: {self.domain}")

        attention = torch.softmax(dot_logits + periodic_logits, dim=-1)
        attention = self.dropout(attention)
        output = torch.einsum("ihj,jhd->ihd", attention, v)
        if value_encoding is not None:
            output = output + torch.einsum("ihj,ihjd->ihd", attention, value_encoding)
        return output.reshape(atom_count, -1)

    def forward(
        self,
        x: torch.Tensor,
        pos: torch.Tensor,
        lattice: torch.Tensor,
        sizes: torch.Tensor,
    ) -> torch.Tensor:
        atom_count = x.shape[0]
        q = self.q_proj(x).view(atom_count, self.num_heads, self.head_dim)
        k = self.k_proj(x).view(atom_count, self.num_heads, self.head_dim)
        v = self.v_proj(x).view(atom_count, self.num_heads, self.head_dim)
        alpha_raw = self.normalized_alpha(q)
        q_chunks = torch.split(q, sizes.tolist())
        k_chunks = torch.split(k, sizes.tolist())
        v_chunks = torch.split(v, sizes.tolist())
        alpha_chunks = torch.split(alpha_raw, sizes.tolist())
        pos_chunks = torch.split(pos, sizes.tolist())
        outputs = [
            self.forward_one(
                q_part, k_part, v_part, alpha_part, pos_part, lattice[index]
            )
            for index, (q_part, k_part, v_part, alpha_part, pos_part) in enumerate(
                zip(q_chunks, k_chunks, v_chunks, alpha_chunks, pos_chunks)
            )
        ]
        return self.out_proj(torch.cat(outputs, dim=0))


class LatticeformerEncoderLayer(nn.Module):
    def __init__(self, config: ModelConfig, domain: str):
        super().__init__()
        self.self_attn = PeriodicMultiheadAttention(config, domain)
        self.linear1 = nn.Linear(config.model_dim, config.ff_dim)
        self.linear2 = nn.Linear(config.ff_dim, config.model_dim)
        self.dropout = nn.Dropout(config.dropout)
        self.dropout1 = nn.Dropout(config.dropout)
        self.dropout2 = nn.Dropout(config.dropout)
        self.norm1 = (
            nn.Identity() if config.use_t_fixup else nn.LayerNorm(config.model_dim)
        )
        self.norm2 = (
            nn.Identity() if config.use_t_fixup else nn.LayerNorm(config.model_dim)
        )
        self.reset_parameters()

    def reset_parameters(self) -> None:
        nn.init.xavier_uniform_(self.linear1.weight)
        nn.init.xavier_uniform_(self.linear2.weight)
        nn.init.zeros_(self.linear1.bias)
        nn.init.zeros_(self.linear2.bias)

    def apply_t_fixup(self, num_layers: int) -> None:
        scale = 0.67 * (num_layers**-0.25)
        with torch.no_grad():
            self.linear1.weight.mul_(scale)
            self.linear2.weight.mul_(scale)
            self.self_attn.v_proj.weight.mul_(scale)
            self.self_attn.out_proj.weight.mul_(scale)

    def forward(
        self,
        x: torch.Tensor,
        pos: torch.Tensor,
        lattice: torch.Tensor,
        sizes: torch.Tensor,
    ) -> torch.Tensor:
        x = self.norm1(x + self.dropout1(self.self_attn(x, pos, lattice, sizes)))
        feed_forward = self.linear2(self.dropout(F.relu(self.linear1(x))))
        return self.norm2(x + self.dropout2(feed_forward))


class LatticeformerEncoder(nn.Module):
    def __init__(self, config: ModelConfig):
        super().__init__()
        domains = self.layer_domains(config)
        self.layers = nn.ModuleList(
            [LatticeformerEncoderLayer(config, domain) for domain in domains]
        )
        if config.use_t_fixup:
            for layer in self.layers:
                layer.apply_t_fixup(config.num_layers)

    @staticmethod
    def layer_domains(config: ModelConfig) -> List[str]:
        if config.domain in {"real", "reci", "multihead"}:
            return [config.domain] * config.num_layers
        sequence = config.domain.split("-")
        return [sequence[index % len(sequence)] for index in range(config.num_layers)]

    def forward(
        self,
        x: torch.Tensor,
        pos: torch.Tensor,
        lattice: torch.Tensor,
        sizes: torch.Tensor,
    ) -> torch.Tensor:
        for layer in self.layers:
            x = layer(x, pos, lattice, sizes)
        return x


class CrystalformerRegressor(nn.Module):
    def __init__(self, config: ModelConfig):
        super().__init__()
        config.validate()
        self.config = copy.deepcopy(config)
        self.input_embedding = nn.Linear(
            config.max_atomic_number, config.model_dim, bias=False
        )
        embedding_std = config.model_dim**-0.5
        if config.use_t_fixup:
            embedding_std *= (9.0 * config.num_layers) ** -0.25
        nn.init.normal_(self.input_embedding.weight, std=embedding_std)
        self.encoder = LatticeformerEncoder(config)
        self.regression_head = nn.Sequential(
            nn.Linear(config.model_dim, config.model_dim),
            nn.ReLU(),
            nn.Linear(config.model_dim, 1),
        )

    def pool(self, x: torch.Tensor, sizes: torch.Tensor) -> torch.Tensor:
        chunks = torch.split(x, sizes.tolist())
        if self.config.pooling == "average":
            return torch.stack([chunk.mean(dim=0) for chunk in chunks], dim=0)
        return torch.stack([chunk.max(dim=0).values for chunk in chunks], dim=0)

    def forward(self, batch: Dict[str, Any]) -> torch.Tensor:
        x = self.input_embedding(batch["x"])
        x = self.encoder(x, batch["pos"], batch["lattice"], batch["sizes"])
        pooled = self.pool(x, batch["sizes"])
        return self.regression_head(pooled).view(-1)


def calculate_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> Dict[str, float]:
    y_true = np.asarray(y_true, dtype=np.float64).reshape(-1)
    y_pred = np.asarray(y_pred, dtype=np.float64).reshape(-1)
    if y_true.size != y_pred.size or y_true.size == 0:
        raise ValueError("Metric arrays must be non-empty and have identical lengths.")
    return {
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "r2": float(r2_score(y_true, y_pred)) if y_true.size >= 2 else float("nan"),
    }


def autocast_context(device: torch.device):
    enabled = USE_AMP and device.type == "cuda"
    if not enabled:
        return nullcontext()
    return torch.autocast(device_type="cuda", dtype=torch.bfloat16)


def train_one_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    scheduler: torch.optim.lr_scheduler.LambdaLR,
    criterion: nn.Module,
    target_mean: float,
    target_std: float,
    device: torch.device,
) -> float:
    model.train()
    total_loss = 0.0
    total_samples = 0
    for batch in tqdm(loader, desc="Training", leave=False):
        batch = move_batch(batch, device)
        normalized_target = (batch["target"] - target_mean) / target_std
        optimizer.zero_grad(set_to_none=True)
        with autocast_context(device):
            prediction = model(batch).view(-1)
            loss = criterion(prediction, normalized_target.view(-1))
        if not torch.isfinite(loss):
            raise FloatingPointError("A non-finite training loss was encountered.")
        loss.backward()
        if GRAD_CLIP_NORM > 0:
            nn.utils.clip_grad_norm_(model.parameters(), GRAD_CLIP_NORM)
        optimizer.step()
        scheduler.step()
        count = int(batch["target"].numel())
        total_loss += float(loss.item()) * count
        total_samples += count
    return total_loss / max(total_samples, 1)


@torch.inference_mode()
def evaluate(
    model: nn.Module,
    loader: DataLoader,
    target_mean: float,
    target_std: float,
    device: torch.device,
) -> Dict[str, Any]:
    model.eval()
    true_values: List[float] = []
    predictions: List[float] = []
    names: List[str] = []
    for batch in loader:
        batch = move_batch(batch, device)
        with autocast_context(device):
            normalized_prediction = model(batch).view(-1)
        prediction = normalized_prediction.float() * target_std + target_mean
        true_values.extend(
            batch["target"].detach().cpu().numpy().astype(float).tolist()
        )
        predictions.extend(prediction.detach().cpu().numpy().astype(float).tolist())
        names.extend(batch["cif_name"])
    y_true = np.asarray(true_values, dtype=np.float64)
    y_pred = np.asarray(predictions, dtype=np.float64)
    return {
        "names": names,
        "y_true": y_true,
        "y_pred": y_pred,
        "metrics": calculate_metrics(y_true, y_pred),
    }


def save_checkpoint(
    path: Path,
    model: nn.Module,
    model_config: ModelConfig,
    best_epoch: int,
    best_val_rmse: float,
    target_mean: float,
    target_std: float,
    optimizer: torch.optim.Optimizer,
) -> None:
    torch.save(
        {
            "model_name": MODEL_NAME,
            "run_version": RUN_VERSION,
            "model_state_dict": model.state_dict(),
            "model_config": asdict(model_config),
            "best_epoch": best_epoch,
            "best_val_rmse": best_val_rmse,
            "seed": SEED,
            "target_name": TARGET_NAME,
            "target_unit": TARGET_UNIT,
            "target_mean": target_mean,
            "target_std": target_std,
            "optimizer_name": optimizer.__class__.__name__,
            "pure_pytorch_periodic_attention": True,
            "fused_cuda_kernels": False,
            "cell_format": CELL_FORMAT,
        },
        path,
    )


def load_trained_model(
    checkpoint_path: str = str(CHECKPOINT_PATH), device: Optional[torch.device] = None
) -> Tuple[CrystalformerRegressor, Dict[str, Any], torch.device]:
    selected_device = device or torch.device(
        "cuda" if torch.cuda.is_available() else "cpu"
    )
    checkpoint = _torch_load(Path(checkpoint_path), map_location=selected_device)
    if checkpoint.get("model_name") != MODEL_NAME:
        raise ValueError(f"Checkpoint model mismatch: {checkpoint.get('model_name')!r}")
    config = ModelConfig(**checkpoint["model_config"])
    model = CrystalformerRegressor(config).to(selected_device)
    model.load_state_dict(checkpoint["model_state_dict"], strict=True)
    model.eval()
    return model, checkpoint, selected_device


def train_model(
    model_config: ModelConfig,
    loaders: Dict[str, DataLoader],
    target_mean: float,
    target_std: float,
    device: torch.device,
    logger: logging.Logger,
    max_epochs: int,
    learning_rate: float,
    weight_decay: float,
    patience: int,
    save_best: bool,
) -> Tuple[CrystalformerRegressor, List[Dict[str, float]], int, float, bool]:
    model = CrystalformerRegressor(model_config).to(device)
    criterion = nn.MSELoss()
    optimizer = torch.optim.AdamW(
        model.parameters(),
        lr=learning_rate,
        betas=ADAM_BETAS,
        weight_decay=weight_decay,
    )
    scheduler = torch.optim.lr_scheduler.LambdaLR(
        optimizer,
        lr_lambda=lambda step: math.sqrt(LR_DECAY_STEPS / (LR_DECAY_STEPS + step)),
    )
    history: List[Dict[str, float]] = []
    best_val_rmse = float("inf")
    best_epoch = 0
    wait_count = 0
    stopped_early = False
    best_state: Optional[Dict[str, torch.Tensor]] = None

    for epoch in range(1, max_epochs + 1):
        train_mse = train_one_epoch(
            model,
            loaders["train"],
            optimizer,
            scheduler,
            criterion,
            target_mean,
            target_std,
            device,
        )
        train_result = evaluate(
            model, loaders["train_eval"], target_mean, target_std, device
        )
        val_result = evaluate(model, loaders["val"], target_mean, target_std, device)
        train_rmse = train_result["metrics"]["rmse"]
        val_rmse = val_result["metrics"]["rmse"]
        learning_rate_now = float(optimizer.param_groups[0]["lr"])
        history.append(
            {
                "epoch": epoch,
                "train_rmse_eV": train_rmse,
                "val_rmse_eV": val_rmse,
                "train_mse_loss": train_mse,
                "learning_rate": learning_rate_now,
            }
        )
        logger.info(
            "Epoch %04d | train RMSE %.6f eV | val RMSE %.6f eV | normalized train MSE %.6f | lr %.8g",
            epoch,
            train_rmse,
            val_rmse,
            train_mse,
            learning_rate_now,
        )
        if val_rmse < best_val_rmse - MIN_DELTA:
            best_val_rmse = val_rmse
            best_epoch = epoch
            wait_count = 0
            best_state = copy.deepcopy(model.state_dict())
            if save_best:
                save_checkpoint(
                    CHECKPOINT_PATH,
                    model,
                    model_config,
                    best_epoch,
                    best_val_rmse,
                    target_mean,
                    target_std,
                    optimizer,
                )
        else:
            wait_count += 1
            if wait_count >= patience:
                stopped_early = True
                logger.info("Early stopping triggered at epoch %d.", epoch)
                break
    if best_state is None:
        raise RuntimeError("Training ended without a valid checkpoint state.")
    model.load_state_dict(best_state)
    return model, history, best_epoch, best_val_rmse, stopped_early


def run_optuna(
    loaders_by_batch_size: Dict[int, Dict[str, DataLoader]],
    target_mean: float,
    target_std: float,
    device: torch.device,
    logger: logging.Logger,
) -> Dict[str, Any]:
    try:
        import optuna
    except ImportError as exc:
        raise ImportError("USE_OPTUNA=True requires the optuna package.") from exc

    def objective(trial: Any) -> float:
        set_global_seed(SEED)
        model_dim = trial.suggest_categorical("model_dim", [96, 128, 192])
        head_num = trial.suggest_categorical("head_num", [4, 8])
        if model_dim % head_num != 0:
            raise optuna.TrialPruned()
        batch_size = trial.suggest_categorical(
            "batch_size", sorted(loaders_by_batch_size)
        )
        config = ModelConfig(
            model_dim=model_dim,
            num_layers=trial.suggest_int("num_layers", 3, 6),
            head_num=head_num,
            ff_dim=trial.suggest_categorical("ff_dim", [256, 512, 768]),
            dropout=trial.suggest_float("dropout", 0.0, 0.2),
            scale_real=trial.suggest_float("scale_real", 1.0, 2.2),
            scale_reci=trial.suggest_float("scale_reci", 1.5, 3.0),
            value_pe_dist_real=trial.suggest_categorical(
                "value_pe_dist_real", [32, 64]
            ),
        )
        learning_rate = trial.suggest_float("learning_rate", 1.0e-4, 1.0e-3, log=True)
        weight_decay = trial.suggest_float("weight_decay", 1.0e-7, 1.0e-3, log=True)
        _, _, _, best_rmse, _ = train_model(
            config,
            loaders_by_batch_size[batch_size],
            target_mean,
            target_std,
            device,
            logger,
            max_epochs=OPTUNA_MAX_EPOCHS,
            learning_rate=learning_rate,
            weight_decay=weight_decay,
            patience=max(10, min(PATIENCE, OPTUNA_MAX_EPOCHS // 3)),
            save_best=False,
        )
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
        gc.collect()
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
    logger.info("Optuna best validation RMSE: %.6f eV", result["best_value"])
    logger.info("Optuna best parameters: %s", result["best_params"])
    return result


def save_prediction_dat(split: str, result: Dict[str, Any]) -> None:
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
        DAT_DIR / f"{MODEL_NAME}_parity_{split}.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )


def parity_limits(results: Iterable[Dict[str, Any]]) -> Tuple[float, float]:
    values = np.concatenate(
        [np.concatenate((result["y_true"], result["y_pred"])) for result in results]
    )
    lower, upper = float(np.min(values)), float(np.max(values))
    margin = max(0.05 * (upper - lower), 0.05)
    return lower - margin, upper + margin


def style_axes(ax: plt.Axes) -> None:
    ax.tick_params(labelsize=TICK_FONTSIZE)
    if USE_GRID:
        ax.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)


def plot_parity(
    split: str, result: Dict[str, Any], limits: Tuple[float, float], color: str
) -> None:
    metrics = result["metrics"]
    fig, ax = plt.subplots(figsize=(6.4, 6.0))
    ax.scatter(
        result["y_true"],
        result["y_pred"],
        s=28,
        alpha=0.78,
        color=color,
        edgecolors="none",
    )
    ax.plot(limits, limits, linestyle="--", color="black", linewidth=1.2)
    ax.set_xlim(limits)
    ax.set_ylim(limits)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(f"{MODEL_NAME}: {split.capitalize()} Set", fontsize=TITLE_FONTSIZE)
    ax.set_xlabel(f"True {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel(f"Predicted {TARGET_NAME} ({TARGET_UNIT})", fontsize=LABEL_FONTSIZE)
    ax.text(
        0.04,
        0.96,
        f"MAE = {metrics['mae']:.4f} eV\nRMSE = {metrics['rmse']:.4f} eV\nR² = {metrics['r2']:.4f}",
        transform=ax.transAxes,
        va="top",
        fontsize=ANNOTATION_FONTSIZE,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.82},
    )
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_parity_{split}.jpg",
        dpi=FIG_DPI,
        bbox_inches="tight",
    )
    plt.close(fig)


def plot_combined_parity(
    results: Dict[str, Dict[str, Any]], limits: Tuple[float, float]
) -> None:
    fig, ax = plt.subplots(figsize=(6.8, 6.2))
    for index, split in enumerate(("train", "val", "test")):
        result = results[split]
        ax.scatter(
            result["y_true"],
            result["y_pred"],
            s=26,
            alpha=0.72,
            color=COLORS[index],
            label=f"{split.capitalize()} (RMSE={result['metrics']['rmse']:.4f} eV)",
            edgecolors="none",
        )
    ax.plot(limits, limits, linestyle="--", color="black", linewidth=1.2, label="y = x")
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
        FIGURE_DIR / f"{MODEL_NAME}_parity_all.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(fig)


def save_combined_dat(results: Dict[str, Dict[str, Any]]) -> None:
    frames = []
    for split in ("train", "val", "test"):
        result = results[split]
        frame = pd.DataFrame(
            {
                "split": split,
                "cif_name": result["names"],
                "true_bandgap_eV": result["y_true"],
                "predicted_bandgap_eV": result["y_pred"],
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


def save_metrics_table(results: Dict[str, Dict[str, Any]]) -> None:
    rows = []
    for split in ("train", "val", "test"):
        metrics = results[split]["metrics"]
        rows.append(
            {
                "split": split,
                "n_samples": len(results[split]["y_true"]),
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


def save_and_plot_history(history: Sequence[Dict[str, float]], best_epoch: int) -> None:
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
        color=COLORS[4],
        zorder=3,
        label=f"Best epoch: {best_epoch}",
    )
    ax.set_title(f"{MODEL_NAME} Training Curve", fontsize=TITLE_FONTSIZE)
    ax.set_xlabel("Epoch", fontsize=LABEL_FONTSIZE)
    ax.set_ylabel("RMSE (eV)", fontsize=LABEL_FONTSIZE)
    ax.legend(fontsize=LEGEND_FONTSIZE)
    style_axes(ax)
    fig.tight_layout()
    fig.savefig(
        FIGURE_DIR / f"{MODEL_NAME}_rmse_curve.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(fig)


@torch.inference_mode()
def predict_cifs(
    model: CrystalformerRegressor,
    cif_paths: Sequence[str],
    checkpoint: Dict[str, Any],
    device: torch.device,
    batch_size: int = BATCH_SIZE,
) -> pd.DataFrame:
    records = []
    for path_string in cif_paths:
        path = Path(path_string).expanduser().resolve()
        if not path.is_file():
            raise FileNotFoundError(f"CIF file not found: {path}")
        records.append(parse_cif(path, 0.0, model.config))
    loader = make_loader(BandgapDataset(records), batch_size, shuffle=False)
    names: List[str] = []
    predictions: List[float] = []
    for batch in loader:
        batch = move_batch(batch, device)
        normalized = model(batch).view(-1)
        values = normalized * float(checkpoint["target_std"]) + float(
            checkpoint["target_mean"]
        )
        names.extend(batch["cif_name"])
        predictions.extend(values.cpu().numpy().astype(float).tolist())
    return pd.DataFrame(
        {"cif_name": names, f"predicted_{TARGET_NAME}_{TARGET_UNIT}": predictions}
    )


def log_environment(
    logger: logging.Logger, device: torch.device, config: ModelConfig
) -> None:
    logger.info("Model name: %s", MODEL_NAME)
    logger.info("Run version: %s", RUN_VERSION)
    logger.info("Output root: %s", OUTPUT_ROOT)
    logger.info("Excel path: %s", resolve_path(EXCEL_PATH))
    logger.info("CIF directory: %s", resolve_path(CIF_DIR))
    logger.info("Excel columns: first column=CIF name, second column=target")
    logger.info("Target: %s (%s)", TARGET_NAME, TARGET_UNIT)
    logger.info("Seed: %d", SEED)
    logger.info(
        "Split ratios: train=%.3f, val=%.3f, test=%.3f",
        TRAIN_RATIO,
        VAL_RATIO,
        TEST_RATIO,
    )
    logger.info("Device: %s", device)
    logger.info("PyTorch version: %s", torch.__version__)
    logger.info("Pure PyTorch periodic attention: enabled")
    logger.info("Fused custom CUDA kernels: disabled")
    logger.info(
        "Expected speed: slower than the upstream fused CUDA-kernel implementation"
    )
    if device.type == "cuda":
        properties = torch.cuda.get_device_properties(device)
        logger.info("GPU: %s", properties.name)
        logger.info("GPU memory: %.2f GiB", properties.total_memory / 1024**3)
        logger.info("CUDA runtime: %s", torch.version.cuda)
    logger.info("Model configuration: %s", asdict(config))
    logger.info("Training loss: MSELoss")
    logger.info(
        "Optimizer: AdamW(lr=%g, betas=%s, weight_decay=%g)",
        LEARNING_RATE,
        ADAM_BETAS,
        WEIGHT_DECAY,
    )
    logger.info(
        "Scheduler: inverse-square-root decay with decay_steps=%g", LR_DECAY_STEPS
    )
    logger.info("Optuna enabled: %s", USE_OPTUNA)


def main() -> None:
    total_start = time.perf_counter()
    ensure_directories()
    set_global_seed(SEED)
    logger = setup_logger()
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    model_config = ModelConfig()
    model_config.validate()
    if not math.isclose(TRAIN_RATIO + VAL_RATIO + TEST_RATIO, 1.0, abs_tol=1.0e-8):
        raise ValueError("TRAIN_RATIO + VAL_RATIO + TEST_RATIO must equal 1.0.")
    log_environment(logger, device, model_config)

    data_start = time.perf_counter()
    rows = read_excel_rows(logger)
    records = load_or_build_records(rows, model_config, logger)
    split_mapping = create_or_load_split(records, logger)
    split_records = {
        split: [record for record in records if split_mapping[record.cif_name] == split]
        for split in ("train", "val", "test")
    }
    logger.info(
        "Samples: train=%d, val=%d, test=%d",
        len(split_records["train"]),
        len(split_records["val"]),
        len(split_records["test"]),
    )
    train_targets = np.asarray(
        [record.target for record in split_records["train"]], dtype=np.float64
    )
    target_mean = float(np.mean(train_targets))
    target_std = float(np.sqrt(np.mean((train_targets - target_mean) ** 2)))
    if not math.isfinite(target_std) or target_std <= 1.0e-12:
        raise ValueError(
            "The training targets have zero or invalid standard deviation."
        )
    logger.info(
        "Training-target normalization: mean=%.8f eV, std=%.8f eV",
        target_mean,
        target_std,
    )
    logger.info(
        "Data preparation time: %s", format_duration(time.perf_counter() - data_start)
    )

    datasets = {
        split: BandgapDataset(values) for split, values in split_records.items()
    }

    def build_loaders(batch_size: int) -> Dict[str, DataLoader]:
        return {
            "train": make_loader(datasets["train"], batch_size, shuffle=True),
            "train_eval": make_loader(datasets["train"], batch_size, shuffle=False),
            "val": make_loader(datasets["val"], batch_size, shuffle=False),
            "test": make_loader(datasets["test"], batch_size, shuffle=False),
        }

    learning_rate = LEARNING_RATE
    weight_decay = WEIGHT_DECAY
    selected_batch_size = BATCH_SIZE
    if USE_OPTUNA:
        candidate_batch_sizes = sorted(
            {max(1, BATCH_SIZE // 2), BATCH_SIZE, BATCH_SIZE * 2}
        )
        optuna_loaders = {size: build_loaders(size) for size in candidate_batch_sizes}
        best = run_optuna(optuna_loaders, target_mean, target_std, device, logger)
        params = best["best_params"]
        selected_batch_size = int(params.pop("batch_size"))
        learning_rate = float(params.pop("learning_rate"))
        weight_decay = float(params.pop("weight_decay"))
        model_config = ModelConfig(**{**asdict(model_config), **params})
        model_config.validate()
        set_global_seed(SEED)

    loaders = build_loaders(selected_batch_size)
    preview_model = CrystalformerRegressor(model_config)
    total_parameters = sum(
        parameter.numel() for parameter in preview_model.parameters()
    )
    trainable_parameters = sum(
        parameter.numel()
        for parameter in preview_model.parameters()
        if parameter.requires_grad
    )
    logger.info("Model structure:\n%s", repr(preview_model))
    logger.info("Total parameters: %d", total_parameters)
    logger.info("Trainable parameters: %d", trainable_parameters)
    del preview_model

    training_start = time.perf_counter()
    try:
        _, history, best_epoch, best_val_rmse, stopped_early = train_model(
            model_config,
            loaders,
            target_mean,
            target_std,
            device,
            logger,
            max_epochs=MAX_EPOCHS,
            learning_rate=learning_rate,
            weight_decay=weight_decay,
            patience=PATIENCE,
            save_best=True,
        )
    except torch.cuda.OutOfMemoryError as exc:
        raise RuntimeError(
            "CUDA ran out of memory. Reduce BATCH_SIZE, MODEL_DIM, LATTICE_RANGE, or VALUE_PE_DIST_REAL."
        ) from exc
    training_time = time.perf_counter() - training_start
    logger.info("Best epoch: %d", best_epoch)
    logger.info("Best validation RMSE: %.6f eV", best_val_rmse)
    logger.info("Early stopping triggered: %s", stopped_early)
    logger.info("Training time: %s", format_duration(training_time))

    evaluation_start = time.perf_counter()
    best_model, checkpoint, device = load_trained_model(str(CHECKPOINT_PATH), device)
    results = {
        split: evaluate(
            best_model,
            loaders["train_eval" if split == "train" else split],
            target_mean,
            target_std,
            device,
        )
        for split in ("train", "val", "test")
    }
    for split, result in results.items():
        metrics = result["metrics"]
        logger.info(
            "%s metrics | MAE %.6f eV | RMSE %.6f eV | R2 %.6f",
            split,
            metrics["mae"],
            metrics["rmse"],
            metrics["r2"],
        )
        save_prediction_dat(split, result)
    save_combined_dat(results)
    save_metrics_table(results)
    limits = parity_limits(results.values())
    for index, split in enumerate(("train", "val", "test")):
        plot_parity(split, results[split], limits, COLORS[index])
    plot_combined_parity(results, limits)
    save_and_plot_history(history, int(checkpoint["best_epoch"]))
    evaluation_time = time.perf_counter() - evaluation_start
    logger.info("Evaluation and plotting time: %s", format_duration(evaluation_time))
    logger.info("Checkpoint: %s", CHECKPOINT_PATH)
    logger.info("Total runtime: %s", format_duration(time.perf_counter() - total_start))


if __name__ == "__main__":
    main()

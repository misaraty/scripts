from __future__ import annotations

import json
import hashlib
import logging
import os
import random
import re
import sys
import time
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import torch
from rdkit import Chem, RDLogger
from rdkit.Chem.rdchem import BondType, HybridizationType
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from torch import Tensor, nn
from torch.utils.data import DataLoader, Dataset


# User configuration (edit these values)
MODEL_NAME = "Chemprop"
RUN_VERSION = "v3"

DATA_PATH = "./data.xlsx"
# Optional units keyed by exact Excel target-column names.
TARGET_UNITS: dict[str, str] = {}
# TARGET_UNITS = {
    # "MIC_Ecoli": "ug/mL",
    # "MIC_Saureus": "ug/mL",
    # "Cytotoxicity": "uM",
# }

SEED = 42
TRAIN_RATIO = 0.8
VAL_RATIO = 0.1
TEST_RATIO = 0.1

BATCH_SIZE = 50
MAX_EPOCHS = 100
LEARNING_RATE = 1.0e-3
WEIGHT_DECAY = 0.0
PATIENCE = 30
MIN_DELTA = 1.0e-6
NUM_WORKERS = 0

MESSAGE_HIDDEN_DIM = 300
MESSAGE_DEPTH = 3
MESSAGE_DROPOUT = 0.0
AGGREGATION = "norm"
AGGREGATION_NORM = 100.0
FFN_HIDDEN_DIM = 300
FFN_NUM_LAYERS = 1
FFN_DROPOUT = 0.0

USE_AMP = True
GRADIENT_CLIP_NORM = 5.0
LR_REDUCTION_FACTOR = 0.5
LR_PATIENCE = 10
MIN_LEARNING_RATE = 1.0e-6

USE_OPTUNA = False
OPTUNA_N_TRIALS = 30
OPTUNA_TIMEOUT = None
OPTUNA_MAX_EPOCHS = 80

FIG_DPI = 600
TITLE_FONTSIZE = 15
LABEL_FONTSIZE = 13
TICK_FONTSIZE = 11
LEGEND_FONTSIZE = 11
ANNOTATION_FONTSIZE = 10
USE_GRID = True
COLORS = ["tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:red"]


SCRIPT_DIR = Path(__file__).resolve().parent
os.chdir(SCRIPT_DIR)
OUTPUT_DIR = SCRIPT_DIR / f"{MODEL_NAME}_{RUN_VERSION}"
FIGURE_DIR = OUTPUT_DIR / "figure"
DAT_DIR = OUTPUT_DIR / "dat"
TABLE_DIR = OUTPUT_DIR / "table"
LOG_DIR = OUTPUT_DIR / "log"
SPLIT_DIR = OUTPUT_DIR / "split"
CHECKPOINT_PATH = OUTPUT_DIR / f"{MODEL_NAME}_best.pt"
SPLIT_PATH = SPLIT_DIR / f"{MODEL_NAME}_split.csv"
LOG_PATH = LOG_DIR / f"{MODEL_NAME}_training.log"

ATOM_FDIM = 72
BOND_FDIM = 14


@dataclass(frozen=True)
class MoleculeRecord:
    sample_id: int
    smiles: str
    canonical_smiles: str
    targets: np.ndarray
    mol: Chem.Mol


@dataclass(frozen=True)
class MolGraph:
    atom_features: np.ndarray
    bond_features: np.ndarray
    edge_index: np.ndarray
    reverse_edge_index: np.ndarray


@dataclass
class BatchMolGraph:
    atom_features: Tensor
    bond_features: Tensor
    edge_index: Tensor
    reverse_edge_index: Tensor
    atom_batch: Tensor
    num_graphs: int

    def to(self, device: torch.device) -> "BatchMolGraph":
        self.atom_features = self.atom_features.to(device, non_blocking=True)
        self.bond_features = self.bond_features.to(device, non_blocking=True)
        self.edge_index = self.edge_index.to(device, non_blocking=True)
        self.reverse_edge_index = self.reverse_edge_index.to(device, non_blocking=True)
        self.atom_batch = self.atom_batch.to(device, non_blocking=True)
        return self


class MoleculeDataset(Dataset):
    def __init__(self, records: Sequence[MoleculeRecord]):
        self.records = list(records)
        self.graphs = [mol_to_graph(record.mol) for record in self.records]

    def __len__(self) -> int:
        return len(self.records)

    def __getitem__(self, index: int) -> tuple[MolGraph, np.ndarray, int]:
        record = self.records[index]
        return self.graphs[index], record.targets, record.sample_id


class TargetScaler:
    def __init__(self, mean: Sequence[float], std: Sequence[float]):
        self.mean = np.asarray(mean, dtype=np.float32)
        self.std = np.asarray(std, dtype=np.float32)
        if self.mean.ndim != 1 or self.std.shape != self.mean.shape:
            raise ValueError(
                "Target scaler arrays must be one-dimensional and aligned."
            )

    @classmethod
    def fit(cls, values: np.ndarray) -> "TargetScaler":
        array = np.asarray(values, dtype=np.float64)
        if array.ndim != 2:
            raise ValueError("Training targets must have shape [samples, tasks].")
        if np.any(np.sum(np.isfinite(array), axis=0) == 0):
            raise ValueError("Every task needs at least one labeled training sample.")
        mean = np.nanmean(array, axis=0)
        std = np.nanstd(array, axis=0, ddof=0)
        std[~np.isfinite(std) | (std < 1.0e-12)] = 1.0
        return cls(mean, std)

    def transform_tensor(self, values: Tensor) -> Tensor:
        mean = torch.as_tensor(self.mean, dtype=values.dtype, device=values.device)
        std = torch.as_tensor(self.std, dtype=values.dtype, device=values.device)
        return (values - mean) / std

    def inverse_tensor(self, values: Tensor) -> Tensor:
        mean = torch.as_tensor(self.mean, dtype=values.dtype, device=values.device)
        std = torch.as_tensor(self.std, dtype=values.dtype, device=values.device)
        return values * std + mean

    def as_dict(self) -> dict[str, list[float]]:
        return {"mean": self.mean.tolist(), "std": self.std.tolist()}


def set_global_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def create_directories() -> None:
    for path in [OUTPUT_DIR, FIGURE_DIR, DAT_DIR, TABLE_DIR, LOG_DIR, SPLIT_DIR]:
        path.mkdir(parents=True, exist_ok=True)


def setup_logger() -> logging.Logger:
    logger = logging.getLogger(MODEL_NAME)
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


def format_duration(seconds: float) -> str:
    seconds_int = max(0, int(round(seconds)))
    hours, remainder = divmod(seconds_int, 3600)
    minutes, secs = divmod(remainder, 60)
    return f"{hours:02d}:{minutes:02d}:{secs:02d}"


def get_device(logger: logging.Logger) -> torch.device:
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    logger.info("Device: %s", device)
    if device.type == "cuda":
        properties = torch.cuda.get_device_properties(0)
        logger.info("GPU: %s", properties.name)
        logger.info("GPU memory: %.2f GiB", properties.total_memory / 1024**3)
        logger.info("CUDA version: %s", torch.version.cuda)
    return device


def resolve_target_units(task_names: Sequence[str]) -> list[str]:
    unknown = sorted(set(TARGET_UNITS) - set(task_names))
    if unknown:
        raise ValueError(f"TARGET_UNITS contains unknown tasks: {unknown}")
    return [str(TARGET_UNITS.get(name, "")) for name in task_names]


def load_records(
    path: str | Path, logger: logging.Logger
) -> tuple[list[MoleculeRecord], list[str], list[str]]:
    path = Path(path)
    if not path.is_absolute():
        path = SCRIPT_DIR / path
    if not path.exists():
        raise FileNotFoundError(f"Data file not found: {path}")

    frame = pd.read_excel(path)
    if frame.shape[1] < 2:
        raise ValueError(
            "data.xlsx must contain SMILES and at least one target column."
        )
    task_names = [str(name).strip() for name in frame.columns[1:]]
    if any(not name for name in task_names) or len(set(task_names)) != len(task_names):
        raise ValueError("Target column names must be non-empty and unique.")
    if len({safe_name(name) for name in task_names}) != len(task_names):
        raise ValueError("Target column names produce duplicate output-safe names.")
    task_units = resolve_target_units(task_names)
    smiles_series = frame.iloc[:, 0]
    target_frame = frame.iloc[:, 1:].apply(pd.to_numeric, errors="coerce")
    target_values = target_frame.to_numpy(dtype=np.float64)
    target_values[~np.isfinite(target_values)] = np.nan
    records: list[MoleculeRecord] = []
    invalid: list[tuple[int, str]] = []
    RDLogger.DisableLog("rdApp.error")

    for row_offset, raw_smiles in enumerate(smiles_series.tolist()):
        row_number = row_offset + 2
        row_targets = target_values[row_offset]
        if pd.isna(raw_smiles) or not str(raw_smiles).strip():
            invalid.append((row_number, "empty SMILES"))
            continue
        smiles = str(raw_smiles).strip()
        if not np.any(np.isfinite(row_targets)):
            invalid.append((row_number, f"no valid target for SMILES {smiles!r}"))
            continue
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or mol.GetNumAtoms() == 0:
            invalid.append((row_number, f"RDKit could not parse SMILES {smiles!r}"))
            continue
        canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        records.append(
            MoleculeRecord(
                sample_id=row_offset,
                smiles=smiles,
                canonical_smiles=canonical,
                targets=row_targets.astype(np.float32),
                mol=mol,
            )
        )

    RDLogger.EnableLog("rdApp.error")
    for row_number, reason in invalid:
        logger.warning("Skipped Excel row %d: %s", row_number, reason)
    logger.info("Tasks (%d): %s", len(task_names), task_names)
    logger.info(
        "Loaded %d valid molecules; skipped %d invalid rows.",
        len(records),
        len(invalid),
    )
    if len(records) < 10:
        raise ValueError(
            "At least 10 valid molecules are required for an 80/10/10 split."
        )
    label_counts = np.sum(
        np.isfinite(np.stack([record.targets for record in records])), axis=0
    )
    for name, count in zip(task_names, label_counts):
        logger.info("Valid labels for %s: %d", name, int(count))
    return records, task_names, task_units


def _one_hot(value: Any, choices: Sequence[Any]) -> list[float]:
    output = [0.0] * (len(choices) + 1)
    try:
        output[choices.index(value)] = 1.0
    except ValueError:
        output[-1] = 1.0
    return output


def atom_features(atom: Chem.Atom) -> np.ndarray:
    atomic_numbers = list(range(1, 37)) + [53]
    degrees = list(range(6))
    formal_charges = [-1, -2, 1, 2, 0]
    chiral_tags = list(range(4))
    hydrogen_counts = list(range(5))
    hybridizations = [
        HybridizationType.S,
        HybridizationType.SP,
        HybridizationType.SP2,
        HybridizationType.SP2D,
        HybridizationType.SP3,
        HybridizationType.SP3D,
        HybridizationType.SP3D2,
    ]
    features = (
        _one_hot(atom.GetAtomicNum(), atomic_numbers)
        + _one_hot(atom.GetTotalDegree(), degrees)
        + _one_hot(atom.GetFormalCharge(), formal_charges)
        + _one_hot(int(atom.GetChiralTag()), chiral_tags)
        + _one_hot(int(atom.GetTotalNumHs()), hydrogen_counts)
        + _one_hot(atom.GetHybridization(), hybridizations)
        + [float(atom.GetIsAromatic()), 0.01 * atom.GetMass()]
    )
    array = np.asarray(features, dtype=np.float32)
    if array.shape[0] != ATOM_FDIM:
        raise RuntimeError(f"Unexpected atom feature dimension: {array.shape[0]}")
    return array


def bond_features(bond: Chem.Bond) -> np.ndarray:
    bond_types = [BondType.SINGLE, BondType.DOUBLE, BondType.TRIPLE, BondType.AROMATIC]
    stereos = list(range(6))
    features = [0.0]
    bond_type_bits = [0.0] * len(bond_types)
    if bond.GetBondType() in bond_types:
        bond_type_bits[bond_types.index(bond.GetBondType())] = 1.0
    features.extend(bond_type_bits)
    features.extend([float(bond.GetIsConjugated()), float(bond.IsInRing())])
    features.extend(_one_hot(int(bond.GetStereo()), stereos))
    array = np.asarray(features, dtype=np.float32)
    if array.shape[0] != BOND_FDIM:
        raise RuntimeError(f"Unexpected bond feature dimension: {array.shape[0]}")
    return array


def mol_to_graph(mol: Chem.Mol) -> MolGraph:
    atoms = np.stack([atom_features(atom) for atom in mol.GetAtoms()]).astype(
        np.float32
    )
    bonds: list[np.ndarray] = []
    sources: list[int] = []
    destinations: list[int] = []

    for bond in mol.GetBonds():
        features = bond_features(bond)
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        sources.extend([begin, end])
        destinations.extend([end, begin])
        bonds.extend([features, features])

    if bonds:
        bond_array = np.stack(bonds).astype(np.float32)
        edge_index = np.asarray([sources, destinations], dtype=np.int64)
        reverse = np.arange(len(bonds), dtype=np.int64).reshape(-1, 2)[:, ::-1].ravel()
    else:
        bond_array = np.empty((0, BOND_FDIM), dtype=np.float32)
        edge_index = np.empty((2, 0), dtype=np.int64)
        reverse = np.empty((0,), dtype=np.int64)

    return MolGraph(atoms, bond_array, edge_index, reverse)


def collate_graphs(
    items: Iterable[tuple[MolGraph, np.ndarray, int]]
) -> tuple[BatchMolGraph, Tensor, Tensor]:
    graphs, targets, sample_ids = zip(*items)
    atom_arrays: list[np.ndarray] = []
    bond_arrays: list[np.ndarray] = []
    edge_arrays: list[np.ndarray] = []
    reverse_arrays: list[np.ndarray] = []
    atom_batches: list[np.ndarray] = []
    atom_offset = 0
    bond_offset = 0

    for graph_index, graph in enumerate(graphs):
        atom_arrays.append(graph.atom_features)
        bond_arrays.append(graph.bond_features)
        edge_arrays.append(graph.edge_index + atom_offset)
        reverse_arrays.append(graph.reverse_edge_index + bond_offset)
        atom_batches.append(
            np.full(graph.atom_features.shape[0], graph_index, dtype=np.int64)
        )
        atom_offset += graph.atom_features.shape[0]
        bond_offset += graph.bond_features.shape[0]

    batch_graph = BatchMolGraph(
        atom_features=torch.from_numpy(np.concatenate(atom_arrays, axis=0)).float(),
        bond_features=torch.from_numpy(np.concatenate(bond_arrays, axis=0)).float(),
        edge_index=torch.from_numpy(np.concatenate(edge_arrays, axis=1)).long(),
        reverse_edge_index=torch.from_numpy(
            np.concatenate(reverse_arrays, axis=0)
        ).long(),
        atom_batch=torch.from_numpy(np.concatenate(atom_batches, axis=0)).long(),
        num_graphs=len(graphs),
    )
    return (
        batch_graph,
        torch.from_numpy(np.stack(targets)).float(),
        torch.tensor(sample_ids, dtype=torch.long),
    )


class ChempropDMPNN(nn.Module):
    def __init__(
        self,
        hidden_dim: int = MESSAGE_HIDDEN_DIM,
        depth: int = MESSAGE_DEPTH,
        message_dropout: float = MESSAGE_DROPOUT,
        aggregation: str = AGGREGATION,
        aggregation_norm: float = AGGREGATION_NORM,
        ffn_hidden_dim: int = FFN_HIDDEN_DIM,
        ffn_num_layers: int = FFN_NUM_LAYERS,
        ffn_dropout: float = FFN_DROPOUT,
        num_tasks: int = 1,
    ):
        super().__init__()
        if depth < 1:
            raise ValueError("MESSAGE_DEPTH must be at least 1.")
        if ffn_num_layers < 0:
            raise ValueError("FFN_NUM_LAYERS cannot be negative.")
        if aggregation not in {"mean", "sum", "norm"}:
            raise ValueError("AGGREGATION must be 'mean', 'sum', or 'norm'.")
        if aggregation == "norm" and aggregation_norm <= 0:
            raise ValueError("AGGREGATION_NORM must be positive.")
        if num_tasks < 1:
            raise ValueError("num_tasks must be at least 1.")
        self.hidden_dim = hidden_dim
        self.depth = depth
        self.aggregation = aggregation
        self.aggregation_norm = aggregation_norm
        self.num_tasks = num_tasks
        self.input_layer = nn.Linear(ATOM_FDIM + BOND_FDIM, hidden_dim, bias=False)
        self.message_layer = nn.Linear(hidden_dim, hidden_dim, bias=False)
        self.atom_output_layer = nn.Linear(ATOM_FDIM + hidden_dim, hidden_dim)
        self.message_dropout = nn.Dropout(message_dropout)
        self.activation = nn.ReLU()
        self.ffn = self._build_ffn(
            hidden_dim, ffn_hidden_dim, ffn_num_layers, ffn_dropout, num_tasks
        )

    @staticmethod
    def _build_ffn(
        input_dim: int,
        hidden_dim: int,
        num_hidden_layers: int,
        dropout: float,
        output_dim: int,
    ) -> nn.Sequential:
        if num_hidden_layers == 0:
            return nn.Sequential(nn.Linear(input_dim, output_dim))
        layers: list[nn.Module] = [nn.Linear(input_dim, hidden_dim)]
        for _ in range(num_hidden_layers - 1):
            layers.extend(
                [nn.ReLU(), nn.Dropout(dropout), nn.Linear(hidden_dim, hidden_dim)]
            )
        layers.extend(
            [nn.ReLU(), nn.Dropout(dropout), nn.Linear(hidden_dim, output_dim)]
        )
        return nn.Sequential(*layers)

    def forward(self, graph: BatchMolGraph) -> Tensor:
        atoms = graph.atom_features
        bonds = graph.bond_features
        sources = graph.edge_index[0]
        destinations = graph.edge_index[1]

        if bonds.shape[0] > 0:
            initial = self.input_layer(torch.cat([atoms[sources], bonds], dim=1))
            hidden = self.activation(initial)
            for _ in range(1, self.depth):
                incoming = torch.zeros(
                    atoms.shape[0],
                    self.hidden_dim,
                    dtype=hidden.dtype,
                    device=hidden.device,
                )
                incoming.index_add_(0, destinations, hidden)
                messages = incoming[sources] - hidden[graph.reverse_edge_index]
                hidden = self.activation(initial + self.message_layer(messages))
                hidden = self.message_dropout(hidden)
            atom_messages = torch.zeros(
                atoms.shape[0],
                self.hidden_dim,
                dtype=hidden.dtype,
                device=hidden.device,
            )
            atom_messages.index_add_(0, destinations, hidden)
        else:
            atom_messages = torch.zeros(
                atoms.shape[0], self.hidden_dim, dtype=atoms.dtype, device=atoms.device
            )

        atom_hidden = self.activation(
            self.atom_output_layer(torch.cat([atoms, atom_messages], dim=1))
        )
        atom_hidden = self.message_dropout(atom_hidden)
        molecular = torch.zeros(
            graph.num_graphs,
            self.hidden_dim,
            dtype=atom_hidden.dtype,
            device=atom_hidden.device,
        )
        molecular.index_add_(0, graph.atom_batch, atom_hidden)
        if self.aggregation == "mean":
            counts = torch.bincount(
                graph.atom_batch, minlength=graph.num_graphs
            ).clamp_min(1)
            molecular = molecular / counts.to(atom_hidden.dtype).unsqueeze(1)
        elif self.aggregation == "norm":
            molecular = molecular / self.aggregation_norm
        return self.ffn(molecular).view(graph.num_graphs, self.num_tasks)


def model_config(
    num_tasks: int, overrides: dict[str, Any] | None = None
) -> dict[str, Any]:
    config: dict[str, Any] = {
        "hidden_dim": MESSAGE_HIDDEN_DIM,
        "depth": MESSAGE_DEPTH,
        "message_dropout": MESSAGE_DROPOUT,
        "aggregation": AGGREGATION,
        "aggregation_norm": AGGREGATION_NORM,
        "ffn_hidden_dim": FFN_HIDDEN_DIM,
        "ffn_num_layers": FFN_NUM_LAYERS,
        "ffn_dropout": FFN_DROPOUT,
        "num_tasks": num_tasks,
    }
    if overrides:
        config.update(overrides)
    return config


def validate_ratios() -> None:
    ratios = np.asarray([TRAIN_RATIO, VAL_RATIO, TEST_RATIO], dtype=float)
    if np.any(ratios <= 0) or not np.isclose(ratios.sum(), 1.0):
        raise ValueError(
            "TRAIN_RATIO, VAL_RATIO, and TEST_RATIO must be positive and sum to 1."
        )


def validate_task_coverage(
    train: Sequence[MoleculeRecord],
    val: Sequence[MoleculeRecord],
    task_names: Sequence[str],
) -> None:
    for split_name, subset in [("train", train), ("validation", val)]:
        counts = np.sum(
            np.isfinite(np.stack([record.targets for record in subset])), axis=0
        )
        missing = [name for name, count in zip(task_names, counts) if count == 0]
        if missing:
            raise ValueError(
                f"The {split_name} split has no labels for tasks {missing}. "
                "Use more labeled data or remove these target columns."
            )


def split_records(
    records: Sequence[MoleculeRecord],
    task_names: Sequence[str],
    logger: logging.Logger,
) -> tuple[list[MoleculeRecord], list[MoleculeRecord], list[MoleculeRecord]]:
    validate_ratios()
    by_id = {record.sample_id: record for record in records}
    target_columns = [f"target::{name}" for name in task_names]

    if SPLIT_PATH.exists():
        try:
            saved = pd.read_csv(SPLIT_PATH)
            required = {
                "sample_id",
                "smiles",
                "canonical_smiles",
                "split",
                *target_columns,
            }
            if not required.issubset(saved.columns):
                raise ValueError("missing required columns")
            if len(saved) != len(records) or set(saved["sample_id"].astype(int)) != set(
                by_id
            ):
                raise ValueError("sample identities differ from the current data")
            for row in saved.itertuples(index=False):
                record = by_id[int(row.sample_id)]
                stored = (
                    saved.loc[saved["sample_id"] == record.sample_id, target_columns]
                    .iloc[0]
                    .to_numpy(dtype=float)
                )
                if record.canonical_smiles != str(
                    row.canonical_smiles
                ) or not np.allclose(
                    record.targets, stored, rtol=1.0e-10, atol=1.0e-12, equal_nan=True
                ):
                    raise ValueError("SMILES or targets differ from the current data")
            split_map = dict(
                zip(saved["sample_id"].astype(int), saved["split"].astype(str))
            )
            train = [
                record
                for record in records
                if split_map.get(record.sample_id) == "train"
            ]
            val = [
                record for record in records if split_map.get(record.sample_id) == "val"
            ]
            test = [
                record
                for record in records
                if split_map.get(record.sample_id) == "test"
            ]
            if (
                not train
                or not val
                or not test
                or len(train) + len(val) + len(test) != len(records)
            ):
                raise ValueError("invalid split labels")
            validate_task_coverage(train, val, task_names)
            logger.info("Reused split: %s", SPLIT_PATH)
            return train, val, test
        except Exception as error:
            logger.warning("Existing split was not reusable: %s", error)

    indices = np.arange(len(records))
    train_indices, temporary_indices = train_test_split(
        indices, test_size=VAL_RATIO + TEST_RATIO, random_state=SEED, shuffle=True
    )
    relative_test_ratio = TEST_RATIO / (VAL_RATIO + TEST_RATIO)
    val_indices, test_indices = train_test_split(
        temporary_indices,
        test_size=relative_test_ratio,
        random_state=SEED,
        shuffle=True,
    )
    train = [records[index] for index in train_indices]
    val = [records[index] for index in val_indices]
    test = [records[index] for index in test_indices]
    validate_task_coverage(train, val, task_names)
    split_labels = np.empty(len(records), dtype=object)
    split_labels[train_indices] = "train"
    split_labels[val_indices] = "val"
    split_labels[test_indices] = "test"
    split_data: dict[str, Any] = {
        "sample_id": [record.sample_id for record in records],
        "smiles": [record.smiles for record in records],
        "canonical_smiles": [record.canonical_smiles for record in records],
        "split": split_labels,
    }
    for task_index, column in enumerate(target_columns):
        split_data[column] = [record.targets[task_index] for record in records]
    pd.DataFrame(split_data).to_csv(SPLIT_PATH, index=False)
    logger.info("Created split: %s", SPLIT_PATH)
    return train, val, test


def make_loader(
    records: Sequence[MoleculeRecord], shuffle: bool, seed: int
) -> DataLoader:
    generator = torch.Generator()
    generator.manual_seed(seed)
    return DataLoader(
        MoleculeDataset(records),
        batch_size=BATCH_SIZE,
        shuffle=shuffle,
        num_workers=NUM_WORKERS,
        pin_memory=torch.cuda.is_available(),
        collate_fn=collate_graphs,
        generator=generator,
        persistent_workers=NUM_WORKERS > 0,
    )


def amp_context(device: torch.device):
    if USE_AMP and device.type == "cuda":
        return torch.autocast(device_type="cuda", dtype=torch.float16)
    return nullcontext()


def make_grad_scaler(device: torch.device):
    enabled = USE_AMP and device.type == "cuda"
    try:
        return torch.amp.GradScaler("cuda", enabled=enabled)
    except (AttributeError, TypeError):
        return torch.cuda.amp.GradScaler(enabled=enabled)


def calculate_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    if y_true.shape != y_pred.shape:
        raise ValueError(f"Metric shape mismatch: {y_true.shape} versus {y_pred.shape}")
    mask = np.isfinite(y_true) & np.isfinite(y_pred)
    y_true = y_true[mask]
    y_pred = y_pred[mask]
    if len(y_true) == 0:
        return {
            "n_samples": 0,
            "mae": float("nan"),
            "rmse": float("nan"),
            "r2": float("nan"),
        }
    return {
        "n_samples": int(len(y_true)),
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "r2": float(r2_score(y_true, y_pred)) if len(y_true) > 1 else float("nan"),
    }


def task_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> list[dict[str, float]]:
    if y_true.ndim != 2 or y_true.shape != y_pred.shape:
        raise ValueError(
            "Predictions and targets must have aligned [samples, tasks] shapes."
        )
    return [
        calculate_metrics(y_true[:, index], y_pred[:, index])
        for index in range(y_true.shape[1])
    ]


def macro_rmse(y_true: np.ndarray, y_pred: np.ndarray) -> float:
    values = np.asarray(
        [item["rmse"] for item in task_metrics(y_true, y_pred)], dtype=float
    )
    if not np.all(np.isfinite(values)):
        raise ValueError("Every task needs at least one valid label for macro RMSE.")
    return float(values.mean())


def finite_mean(values: Sequence[float]) -> float:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    return float(array.mean()) if len(array) else float("nan")


def train_one_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    scaler: Any,
    target_scaler: TargetScaler,
    device: torch.device,
) -> tuple[float, float]:
    model.train()
    criterion = nn.MSELoss(reduction="sum")
    normalized_loss_sum = 0.0
    labeled_count = 0
    y_true: list[np.ndarray] = []
    y_pred: list[np.ndarray] = []

    for graph, targets, _ in loader:
        graph = graph.to(device)
        targets = targets.to(device, non_blocking=True)
        mask = torch.isfinite(targets)
        if not bool(mask.any()):
            continue
        normalized_targets = target_scaler.transform_tensor(targets)
        optimizer.zero_grad(set_to_none=True)
        with amp_context(device):
            normalized_predictions = model(graph)
            loss = (
                criterion(normalized_predictions[mask], normalized_targets[mask])
                / mask.sum()
            )
        if not torch.isfinite(loss):
            raise FloatingPointError("A non-finite training loss was encountered.")
        scaler.scale(loss).backward()
        scaler.unscale_(optimizer)
        if GRADIENT_CLIP_NORM > 0:
            nn.utils.clip_grad_norm_(model.parameters(), GRADIENT_CLIP_NORM)
        scaler.step(optimizer)
        scaler.update()

        batch_labeled = int(mask.sum().item())
        normalized_loss_sum += float(loss.detach()) * batch_labeled
        labeled_count += batch_labeled
        predictions = target_scaler.inverse_tensor(normalized_predictions.detach())
        y_true.append(targets.detach().cpu().numpy())
        y_pred.append(predictions.float().cpu().numpy())

    truth = np.concatenate(y_true, axis=0)
    prediction = np.concatenate(y_pred, axis=0)
    return normalized_loss_sum / max(labeled_count, 1), macro_rmse(truth, prediction)


@torch.inference_mode()
def evaluate_loader(
    model: nn.Module,
    loader: DataLoader,
    target_scaler: TargetScaler,
    device: torch.device,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    model.eval()
    criterion = nn.MSELoss(reduction="sum")
    normalized_loss_sum = 0.0
    labeled_count = 0
    truths: list[np.ndarray] = []
    predictions: list[np.ndarray] = []
    sample_ids: list[np.ndarray] = []

    for graph, targets, ids in loader:
        graph = graph.to(device)
        targets = targets.to(device, non_blocking=True)
        mask = torch.isfinite(targets)
        with amp_context(device):
            normalized_predictions = model(graph)
            normalized_targets = target_scaler.transform_tensor(targets)
            loss = criterion(normalized_predictions[mask], normalized_targets[mask])
        raw_predictions = target_scaler.inverse_tensor(normalized_predictions)
        normalized_loss_sum += float(loss.detach())
        labeled_count += int(mask.sum().item())
        truths.append(targets.float().cpu().numpy())
        predictions.append(raw_predictions.float().cpu().numpy())
        sample_ids.append(ids.numpy())

    return (
        np.concatenate(truths),
        np.concatenate(predictions),
        np.concatenate(sample_ids),
        normalized_loss_sum / max(labeled_count, 1),
    )


def save_checkpoint(
    model: nn.Module,
    config: dict[str, Any],
    target_scaler: TargetScaler,
    best_epoch: int,
    best_val_rmse: float,
    path: Path,
    task_names: Sequence[str],
    task_units: Sequence[str],
) -> None:
    checkpoint = {
        "model_name": MODEL_NAME,
        "run_version": RUN_VERSION,
        "model_state_dict": model.state_dict(),
        "model_config": config,
        "best_epoch": int(best_epoch),
        "best_val_rmse": float(best_val_rmse),
        "seed": SEED,
        "task_names": list(task_names),
        "task_units": list(task_units),
        "num_tasks": len(task_names),
        "target_scaler": target_scaler.as_dict(),
        "atom_feature_mode": "Chemprop v2 multi-hot",
        "atom_fdim": ATOM_FDIM,
        "bond_fdim": BOND_FDIM,
        "chemprop_source_version": "2.3.1",
    }
    torch.save(checkpoint, path)


def torch_load_checkpoint(path: str | Path, device: torch.device) -> dict[str, Any]:
    try:
        return torch.load(path, map_location=device, weights_only=False)
    except TypeError:
        return torch.load(path, map_location=device)


def load_trained_model(
    checkpoint_path: str | Path = CHECKPOINT_PATH,
    device: torch.device | str | None = None,
) -> tuple[ChempropDMPNN, TargetScaler, dict[str, Any]]:
    selected_device = (
        torch.device(device)
        if device is not None
        else torch.device("cuda" if torch.cuda.is_available() else "cpu")
    )
    checkpoint = torch_load_checkpoint(checkpoint_path, selected_device)
    if checkpoint.get("model_name") != MODEL_NAME:
        raise ValueError(
            f"Unexpected model name in checkpoint: {checkpoint.get('model_name')}"
        )
    model = ChempropDMPNN(**checkpoint["model_config"]).to(selected_device)
    model.load_state_dict(checkpoint["model_state_dict"], strict=True)
    model.eval()
    scaler_info = checkpoint["target_scaler"]
    target_scaler = TargetScaler(scaler_info["mean"], scaler_info["std"])
    return model, target_scaler, checkpoint


def train_model(
    train_loader: DataLoader,
    val_loader: DataLoader,
    target_scaler: TargetScaler,
    device: torch.device,
    config: dict[str, Any],
    logger: logging.Logger,
    task_names: Sequence[str],
    task_units: Sequence[str],
    max_epochs: int = MAX_EPOCHS,
    checkpoint_path: Path | None = CHECKPOINT_PATH,
    trial: Any = None,
) -> tuple[ChempropDMPNN, list[dict[str, float]], int, float, bool]:
    model = ChempropDMPNN(**config).to(device)
    optimizer = torch.optim.AdamW(
        model.parameters(), lr=LEARNING_RATE, weight_decay=WEIGHT_DECAY
    )
    scheduler = torch.optim.lr_scheduler.ReduceLROnPlateau(
        optimizer,
        mode="min",
        factor=LR_REDUCTION_FACTOR,
        patience=LR_PATIENCE,
        min_lr=MIN_LEARNING_RATE,
    )
    scaler = make_grad_scaler(device)
    history: list[dict[str, float]] = []
    best_val_rmse = float("inf")
    best_epoch = 0
    epochs_without_improvement = 0
    stopped_early = False
    best_state: dict[str, Tensor] | None = None

    for epoch in range(1, max_epochs + 1):
        train_loss, train_rmse = train_one_epoch(
            model, train_loader, optimizer, scaler, target_scaler, device
        )
        val_true, val_pred, _, val_loss = evaluate_loader(
            model, val_loader, target_scaler, device
        )
        val_rmse = macro_rmse(val_true, val_pred)
        current_lr = float(optimizer.param_groups[0]["lr"])
        scheduler.step(val_rmse)
        train_true, train_pred, _, _ = evaluate_loader(
            model, train_loader, target_scaler, device
        )
        train_rmse = macro_rmse(train_true, train_pred)
        train_task_values = task_metrics(train_true, train_pred)
        val_task_values = task_metrics(val_true, val_pred)
        epoch_row = {
            "epoch": float(epoch),
            "train_rmse_macro": train_rmse,
            "val_rmse_macro": val_rmse,
            "train_mse_loss": train_loss,
            "val_mse_loss": val_loss,
            "learning_rate": current_lr,
        }
        for task_index, task_name in enumerate(task_names):
            key = safe_name(task_name)
            epoch_row[f"train_rmse::{key}"] = train_task_values[task_index]["rmse"]
            epoch_row[f"val_rmse::{key}"] = val_task_values[task_index]["rmse"]
        history.append(epoch_row)
        logger.info(
            "Epoch %04d | train macro RMSE %.6f | val macro RMSE %.6f | lr %.3e",
            epoch,
            train_rmse,
            val_rmse,
            current_lr,
        )

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
                    model,
                    config,
                    target_scaler,
                    best_epoch,
                    best_val_rmse,
                    checkpoint_path,
                    task_names,
                    task_units,
                )
        else:
            epochs_without_improvement += 1

        if trial is not None:
            trial.report(val_rmse, epoch)
            if trial.should_prune():
                import optuna

                raise optuna.TrialPruned()

        if epochs_without_improvement >= PATIENCE:
            stopped_early = True
            logger.info("Early stopping at epoch %d.", epoch)
            break

    if best_state is None:
        raise RuntimeError("Training did not produce a finite validation score.")
    model.load_state_dict(best_state)
    return model, history, best_epoch, best_val_rmse, stopped_early


def run_optuna(
    train_loader: DataLoader,
    val_loader: DataLoader,
    target_scaler: TargetScaler,
    device: torch.device,
    logger: logging.Logger,
    task_names: Sequence[str],
    task_units: Sequence[str],
) -> dict[str, Any]:
    try:
        import optuna
    except ImportError as error:
        raise ImportError(
            "Optuna is enabled but not installed. Run: pip install optuna"
        ) from error

    optuna.logging.set_verbosity(optuna.logging.WARNING)

    def objective(trial: Any) -> float:
        set_global_seed(SEED)
        overrides = {
            "hidden_dim": trial.suggest_categorical("hidden_dim", [200, 300, 400, 500]),
            "depth": trial.suggest_int("depth", 2, 5),
            "message_dropout": trial.suggest_float(
                "message_dropout", 0.0, 0.3, step=0.05
            ),
            "ffn_hidden_dim": trial.suggest_categorical(
                "ffn_hidden_dim", [200, 300, 400, 500]
            ),
            "ffn_num_layers": trial.suggest_int("ffn_num_layers", 1, 3),
            "ffn_dropout": trial.suggest_float("ffn_dropout", 0.0, 0.4, step=0.05),
        }
        _, _, _, score, _ = train_model(
            train_loader,
            val_loader,
            target_scaler,
            device,
            model_config(len(task_names), overrides),
            logger,
            task_names,
            task_units,
            max_epochs=OPTUNA_MAX_EPOCHS,
            checkpoint_path=None,
            trial=trial,
        )
        if device.type == "cuda":
            torch.cuda.empty_cache()
        return score

    study = optuna.create_study(
        direction="minimize", sampler=optuna.samplers.TPESampler(seed=SEED)
    )
    study.optimize(objective, n_trials=OPTUNA_N_TRIALS, timeout=OPTUNA_TIMEOUT)
    result = {"best_value": float(study.best_value), "best_params": study.best_params}
    path = TABLE_DIR / f"{MODEL_NAME}_optuna_best_params.json"
    path.write_text(json.dumps(result, indent=2), encoding="utf-8")
    logger.info("Optuna best macro validation RMSE: %.6f", study.best_value)
    logger.info("Optuna best parameters: %s", study.best_params)
    return study.best_params


def safe_name(value: str) -> str:
    original = str(value)
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", original).strip("._-")
    if not cleaned:
        cleaned = "target"
    if cleaned != original:
        digest = hashlib.sha1(original.encode("utf-8")).hexdigest()[:8]
        cleaned = f"{cleaned}_{digest}"
    return cleaned


def output_stem(task_name: str, num_tasks: int) -> str:
    return MODEL_NAME if num_tasks == 1 else f"{MODEL_NAME}_{safe_name(task_name)}"


def prediction_frame(
    records_by_id: dict[int, MoleculeRecord],
    sample_ids: np.ndarray,
    truth: np.ndarray,
    prediction: np.ndarray,
    split: str,
    task_names: Sequence[str],
    task_units: Sequence[str],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for row_index, sample_id in enumerate(sample_ids):
        record = records_by_id[int(sample_id)]
        for task_index, (task_name, task_unit) in enumerate(
            zip(task_names, task_units)
        ):
            true_value = truth[row_index, task_index]
            if not np.isfinite(true_value):
                continue
            predicted_value = prediction[row_index, task_index]
            rows.append(
                {
                    "split": split,
                    "sample_id": int(sample_id),
                    "smiles": record.smiles,
                    "canonical_smiles": record.canonical_smiles,
                    "task": task_name,
                    "unit": task_unit,
                    "true_target": float(true_value),
                    "predicted_target": float(predicted_value),
                    "error": float(predicted_value - true_value),
                    "absolute_error": float(abs(predicted_value - true_value)),
                }
            )
    return pd.DataFrame(rows)


def save_prediction_data(
    frames: dict[str, pd.DataFrame], task_names: Sequence[str]
) -> None:
    combined = pd.concat(list(frames.values()), ignore_index=True)
    combined.to_csv(
        DAT_DIR / f"{MODEL_NAME}_parity_all_tasks.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )
    for task_name in task_names:
        stem = output_stem(task_name, len(task_names))
        task_frames = []
        for split, frame in frames.items():
            selected = frame.loc[frame["task"] == task_name].copy()
            selected.drop(columns=["split"]).to_csv(
                DAT_DIR / f"{stem}_parity_{split}.dat",
                sep="\t",
                index=False,
                float_format="%.8f",
            )
            task_frames.append(selected)
        pd.concat(task_frames, ignore_index=True).to_csv(
            DAT_DIR / f"{stem}_parity_all.dat",
            sep="\t",
            index=False,
            float_format="%.8f",
        )


def apply_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "axes.titlesize": TITLE_FONTSIZE,
            "axes.labelsize": LABEL_FONTSIZE,
            "xtick.labelsize": TICK_FONTSIZE,
            "ytick.labelsize": TICK_FONTSIZE,
            "legend.fontsize": LEGEND_FONTSIZE,
        }
    )


def target_axis_label(prefix: str, task_name: str, task_unit: str) -> str:
    suffix = f" ({task_unit})" if task_unit else ""
    return f"{prefix} {task_name}{suffix}"


def parity_limits(
    frames: dict[str, pd.DataFrame], task_name: str
) -> tuple[float, float]:
    values = np.concatenate(
        [
            frame.loc[frame["task"] == task_name, ["true_target", "predicted_target"]]
            .to_numpy()
            .ravel()
            for frame in frames.values()
        ]
    )
    minimum = float(np.nanmin(values))
    maximum = float(np.nanmax(values))
    span = maximum - minimum
    margin = 0.05 * span if span > 0 else max(abs(maximum) * 0.05, 0.1)
    return minimum - margin, maximum + margin


def plot_single_parity(
    frame: pd.DataFrame,
    split: str,
    color: str,
    limits: tuple[float, float],
    task_name: str,
    task_unit: str,
    stem: str,
) -> None:
    metrics = calculate_metrics(
        frame["true_target"].to_numpy(), frame["predicted_target"].to_numpy()
    )
    figure, axis = plt.subplots(figsize=(6.5, 6.0))
    axis.scatter(
        frame["true_target"],
        frame["predicted_target"],
        s=28,
        alpha=0.78,
        color=color,
        edgecolors="none",
    )
    axis.plot(
        limits, limits, linestyle="--", color="black", linewidth=1.2, label="y = x"
    )
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(target_axis_label("True", task_name, task_unit))
    axis.set_ylabel(target_axis_label("Predicted", task_name, task_unit))
    axis.set_title(f"{MODEL_NAME} {split.capitalize()} Parity Plot: {task_name}")
    unit_text = f" {task_unit}" if task_unit else ""
    axis.text(
        0.04,
        0.96,
        f"MAE = {metrics['mae']:.4f}{unit_text}\nRMSE = {metrics['rmse']:.4f}{unit_text}\n$R^2$ = {metrics['r2']:.4f}",
        transform=axis.transAxes,
        va="top",
        fontsize=ANNOTATION_FONTSIZE,
        bbox={"boxstyle": "round", "facecolor": "white", "alpha": 0.8},
    )
    if USE_GRID:
        axis.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)
    axis.legend()
    figure.tight_layout()
    figure.savefig(
        FIGURE_DIR / f"{stem}_parity_{split}.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(figure)


def plot_combined_parity(
    frames: dict[str, pd.DataFrame],
    limits: tuple[float, float],
    task_name: str,
    task_unit: str,
    stem: str,
) -> None:
    figure, axis = plt.subplots(figsize=(6.8, 6.2))
    for index, split in enumerate(["train", "val", "test"]):
        frame = frames[split].loc[frames[split]["task"] == task_name]
        metrics = calculate_metrics(
            frame["true_target"].to_numpy(), frame["predicted_target"].to_numpy()
        )
        axis.scatter(
            frame["true_target"],
            frame["predicted_target"],
            s=27,
            alpha=0.72,
            color=COLORS[index],
            edgecolors="none",
            label=f"{split.capitalize()} (RMSE={metrics['rmse']:.4f})",
        )
    axis.plot(
        limits, limits, linestyle="--", color="black", linewidth=1.2, label="y = x"
    )
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel(target_axis_label("True", task_name, task_unit))
    axis.set_ylabel(target_axis_label("Predicted", task_name, task_unit))
    axis.set_title(f"{MODEL_NAME} Combined Parity Plot: {task_name}")
    if USE_GRID:
        axis.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)
    axis.legend()
    figure.tight_layout()
    figure.savefig(
        FIGURE_DIR / f"{stem}_parity_all.jpg", dpi=FIG_DPI, bbox_inches="tight"
    )
    plt.close(figure)


def plot_rmse_curve(
    frame: pd.DataFrame,
    best_epoch: int,
    train_column: str,
    val_column: str,
    title: str,
    ylabel: str,
    path: Path,
) -> None:
    figure, axis = plt.subplots(figsize=(7.2, 5.2))
    axis.plot(frame["epoch"], frame[train_column], color=COLORS[0], label="Train RMSE")
    axis.plot(
        frame["epoch"], frame[val_column], color=COLORS[1], label="Validation RMSE"
    )
    best_value = float(frame.loc[frame["epoch"] == best_epoch, val_column].iloc[0])
    axis.scatter(
        [best_epoch], [best_value], color=COLORS[4], s=48, zorder=4, label="Best epoch"
    )
    axis.axvline(best_epoch, color=COLORS[4], linestyle="--", linewidth=1.0, alpha=0.75)
    axis.set_xlabel("Epoch")
    axis.set_ylabel(ylabel)
    axis.set_title(title)
    if USE_GRID:
        axis.grid(True, linestyle="--", linewidth=0.6, alpha=0.35)
    axis.legend()
    figure.tight_layout()
    figure.savefig(path, dpi=FIG_DPI, bbox_inches="tight")
    plt.close(figure)


def save_and_plot_history(
    history: Sequence[dict[str, float]],
    best_epoch: int,
    task_names: Sequence[str],
    task_units: Sequence[str],
) -> None:
    frame = pd.DataFrame(history)
    frame["epoch"] = frame["epoch"].astype(int)
    frame.to_csv(
        DAT_DIR / f"{MODEL_NAME}_rmse_curve.dat",
        sep="\t",
        index=False,
        float_format="%.8f",
    )
    plot_rmse_curve(
        frame,
        best_epoch,
        "train_rmse_macro",
        "val_rmse_macro",
        f"{MODEL_NAME} Macro-Averaged Training Curve",
        "Macro-averaged RMSE",
        FIGURE_DIR / f"{MODEL_NAME}_rmse_curve.jpg",
    )
    if len(task_names) > 1:
        for task_name, task_unit in zip(task_names, task_units):
            key = safe_name(task_name)
            ylabel = f"RMSE ({task_unit})" if task_unit else "RMSE"
            plot_rmse_curve(
                frame,
                best_epoch,
                f"train_rmse::{key}",
                f"val_rmse::{key}",
                f"{MODEL_NAME} Training Curve: {task_name}",
                ylabel,
                FIGURE_DIR
                / f"{output_stem(task_name, len(task_names))}_rmse_curve.jpg",
            )


def save_metrics_table(
    frames: dict[str, pd.DataFrame],
    task_names: Sequence[str],
    task_units: Sequence[str],
) -> dict[str, dict[str, dict[str, float]]]:
    rows: list[dict[str, Any]] = []
    result: dict[str, dict[str, dict[str, float]]] = {}
    for split in ["train", "val", "test"]:
        result[split] = {}
        rmse_values = []
        for task_name, task_unit in zip(task_names, task_units):
            selected = frames[split].loc[frames[split]["task"] == task_name]
            metrics = calculate_metrics(
                selected["true_target"].to_numpy(),
                selected["predicted_target"].to_numpy(),
            )
            result[split][task_name] = metrics
            if np.isfinite(metrics["rmse"]):
                rmse_values.append(metrics["rmse"])
            rows.append(
                {"split": split, "task": task_name, "unit": task_unit, **metrics}
            )
        rows.append(
            {
                "split": split,
                "task": "__macro__",
                "unit": "",
                "n_samples": int(
                    sum(result[split][name]["n_samples"] for name in task_names)
                ),
                "mae": finite_mean([result[split][name]["mae"] for name in task_names]),
                "rmse": float(np.mean(rmse_values)) if rmse_values else float("nan"),
                "r2": finite_mean([result[split][name]["r2"] for name in task_names]),
            }
        )
    pd.DataFrame(rows).to_csv(
        TABLE_DIR / f"{MODEL_NAME}_metrics.dat",
        sep="\t",
        index=False,
        float_format="%.6f",
    )
    return result


@torch.inference_mode()
def predict_smiles(
    smiles_list: Sequence[str],
    checkpoint_path: str | Path = CHECKPOINT_PATH,
    device: torch.device | str | None = None,
    batch_size: int = BATCH_SIZE,
) -> pd.DataFrame:
    selected_device = (
        torch.device(device)
        if device is not None
        else torch.device("cuda" if torch.cuda.is_available() else "cpu")
    )
    model, target_scaler, checkpoint = load_trained_model(
        checkpoint_path, selected_device
    )
    task_names = list(checkpoint["task_names"])
    task_units = list(checkpoint.get("task_units", [""] * len(task_names)))
    valid_records: list[MoleculeRecord] = []
    results: list[dict[str, Any]] = []
    for index, raw_smiles in enumerate(smiles_list):
        smiles = str(raw_smiles).strip()
        mol = Chem.MolFromSmiles(smiles)
        if mol is None or mol.GetNumAtoms() == 0:
            row = {"input_index": index, "smiles": smiles, "canonical_smiles": ""}
            row.update({f"prediction::{name}": np.nan for name in task_names})
            results.append(row)
            continue
        canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        valid_records.append(
            MoleculeRecord(
                index,
                smiles,
                canonical,
                np.zeros(len(task_names), dtype=np.float32),
                mol,
            )
        )

    if valid_records:
        loader = DataLoader(
            MoleculeDataset(valid_records),
            batch_size=batch_size,
            shuffle=False,
            num_workers=0,
            collate_fn=collate_graphs,
        )
        predictions_by_id: dict[int, np.ndarray] = {}
        model.eval()
        for graph, _, ids in loader:
            graph = graph.to(selected_device)
            normalized = model(graph)
            predictions = target_scaler.inverse_tensor(normalized).float().cpu().numpy()
            predictions_by_id.update(
                {int(i): p for i, p in zip(ids.numpy(), predictions)}
            )
        for record in valid_records:
            row = {
                "input_index": record.sample_id,
                "smiles": record.smiles,
                "canonical_smiles": record.canonical_smiles,
            }
            row.update(
                {
                    f"prediction::{name}": float(
                        predictions_by_id[record.sample_id][index]
                    )
                    for index, name in enumerate(task_names)
                }
            )
            results.append(row)

    frame = pd.DataFrame(results).sort_values("input_index").reset_index(drop=True)
    for name, unit in zip(task_names, task_units):
        frame[f"unit::{name}"] = unit
    return frame


def log_configuration(
    logger: logging.Logger,
    model: nn.Module,
    config: dict[str, Any],
    train_size: int,
    val_size: int,
    test_size: int,
    target_scaler: TargetScaler,
    task_names: Sequence[str],
    task_units: Sequence[str],
) -> None:
    total_parameters = sum(parameter.numel() for parameter in model.parameters())
    trainable_parameters = sum(
        parameter.numel() for parameter in model.parameters() if parameter.requires_grad
    )
    logger.info("Model: %s", MODEL_NAME)
    logger.info("Run version: %s", RUN_VERSION)
    logger.info("Output directory: %s", OUTPUT_DIR)
    logger.info("Data path: %s", DATA_PATH)
    logger.info("Excel columns: first column = SMILES; remaining columns = targets")
    logger.info("Tasks: %s", list(task_names))
    logger.info("Task units: %s", dict(zip(task_names, task_units)))
    logger.info("Seed: %d", SEED)
    logger.info("Split ratios: %.2f / %.2f / %.2f", TRAIN_RATIO, VAL_RATIO, TEST_RATIO)
    logger.info(
        "Split sizes: train=%d, val=%d, test=%d", train_size, val_size, test_size
    )
    logger.info("Target scaler fitted on training data: %s", target_scaler.as_dict())
    logger.info("Model configuration: %s", config)
    logger.info("Batch size: %d", BATCH_SIZE)
    logger.info("Maximum epochs: %d", MAX_EPOCHS)
    logger.info("Optimizer: AdamW(lr=%g, weight_decay=%g)", LEARNING_RATE, WEIGHT_DECAY)
    logger.info("Loss: masked MSELoss in per-task standardized target space")
    logger.info(
        "Scheduler: ReduceLROnPlateau(factor=%g, patience=%d)",
        LR_REDUCTION_FACTOR,
        LR_PATIENCE,
    )
    logger.info("Early stopping: patience=%d, min_delta=%g", PATIENCE, MIN_DELTA)
    logger.info("AMP enabled: %s", USE_AMP)
    logger.info("Optuna enabled: %s", USE_OPTUNA)
    logger.info("Total parameters: %d", total_parameters)
    logger.info("Trainable parameters: %d", trainable_parameters)
    logger.info("Model architecture:\n%s", model)


def main() -> None:
    total_start = time.perf_counter()
    create_directories()
    logger = setup_logger()
    set_global_seed(SEED)
    apply_plot_style()

    try:
        device = get_device(logger)
        data_start = time.perf_counter()
        records, task_names, task_units = load_records(DATA_PATH, logger)
        train_records, val_records, test_records = split_records(
            records, task_names, logger
        )
        records_by_id = {record.sample_id: record for record in records}
        target_scaler = TargetScaler.fit(
            np.stack([record.targets for record in train_records])
        )
        train_loader = make_loader(train_records, shuffle=True, seed=SEED)
        train_eval_loader = make_loader(train_records, shuffle=False, seed=SEED)
        val_loader = make_loader(val_records, shuffle=False, seed=SEED)
        test_loader = make_loader(test_records, shuffle=False, seed=SEED)
        data_elapsed = time.perf_counter() - data_start

        config = model_config(len(task_names))
        if USE_OPTUNA:
            best_params = run_optuna(
                train_loader,
                val_loader,
                target_scaler,
                device,
                logger,
                task_names,
                task_units,
            )
            config = model_config(len(task_names), best_params)
            set_global_seed(SEED)

        preview_model = ChempropDMPNN(**config).to(device)
        log_configuration(
            logger,
            preview_model,
            config,
            len(train_records),
            len(val_records),
            len(test_records),
            target_scaler,
            task_names,
            task_units,
        )
        del preview_model

        train_start = time.perf_counter()
        _, history, best_epoch, best_val_rmse, stopped_early = train_model(
            train_loader,
            val_loader,
            target_scaler,
            device,
            config,
            logger,
            task_names,
            task_units,
        )
        train_elapsed = time.perf_counter() - train_start

        evaluation_start = time.perf_counter()
        best_model, loaded_scaler, checkpoint = load_trained_model(
            CHECKPOINT_PATH, device
        )
        datasets = {
            "train": (train_eval_loader, train_records),
            "val": (val_loader, val_records),
            "test": (test_loader, test_records),
        }
        frames: dict[str, pd.DataFrame] = {}
        for split, (loader, _) in datasets.items():
            truth, prediction, ids, _ = evaluate_loader(
                best_model, loader, loaded_scaler, device
            )
            frames[split] = prediction_frame(
                records_by_id, ids, truth, prediction, split, task_names, task_units
            )

        save_prediction_data(frames, task_names)
        for task_name, task_unit in zip(task_names, task_units):
            limits = parity_limits(frames, task_name)
            stem = output_stem(task_name, len(task_names))
            for index, split in enumerate(["train", "val", "test"]):
                selected = frames[split].loc[frames[split]["task"] == task_name]
                if not selected.empty:
                    plot_single_parity(
                        selected,
                        split,
                        COLORS[index],
                        limits,
                        task_name,
                        task_unit,
                        stem,
                    )
            plot_combined_parity(frames, limits, task_name, task_unit, stem)
        save_and_plot_history(
            history, int(checkpoint["best_epoch"]), task_names, task_units
        )
        metrics_by_split = save_metrics_table(frames, task_names, task_units)
        evaluation_elapsed = time.perf_counter() - evaluation_start

        logger.info("Best epoch: %d", best_epoch)
        logger.info("Best macro validation RMSE: %.6f", best_val_rmse)
        logger.info("Early stopping triggered: %s", stopped_early)
        for split in ["train", "val", "test"]:
            for task_name, task_unit in zip(task_names, task_units):
                metrics = metrics_by_split[split][task_name]
                logger.info(
                    "%s | %s | n %d | MAE %.6f %s | RMSE %.6f %s | R2 %.6f",
                    split.capitalize(),
                    task_name,
                    metrics["n_samples"],
                    metrics["mae"],
                    task_unit,
                    metrics["rmse"],
                    task_unit,
                    metrics["r2"],
                )

        total_elapsed = time.perf_counter() - total_start
        logger.info(
            "Data processing time: %.2f s (%s)",
            data_elapsed,
            format_duration(data_elapsed),
        )
        logger.info(
            "Training time: %.2f s (%s)", train_elapsed, format_duration(train_elapsed)
        )
        logger.info(
            "Evaluation and plotting time: %.2f s (%s)",
            evaluation_elapsed,
            format_duration(evaluation_elapsed),
        )
        logger.info(
            "Total runtime: %.2f s (%s)", total_elapsed, format_duration(total_elapsed)
        )
        logger.info("Best checkpoint: %s", CHECKPOINT_PATH)
    except torch.cuda.OutOfMemoryError as error:
        logger.exception(
            "CUDA ran out of memory. Reduce BATCH_SIZE or hidden dimensions."
        )
        raise RuntimeError("CUDA out of memory") from error
    except Exception:
        logger.exception("Execution failed.")
        raise


if __name__ == "__main__":
    main()

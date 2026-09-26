from __future__ import annotations

import argparse
import hashlib
import json
import logging
import math
import os
import platform
import random
import shutil
import subprocess
import sys
import time
import traceback
import warnings
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from pathlib import Path
from typing import Any, Sequence
from xml.etree import ElementTree as ET

import joblib
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter, map_coordinates
from scipy.optimize import minimize
from scipy.special import sph_harm_y
from scipy.stats import wilcoxon
from sklearn.calibration import calibration_curve
from sklearn.metrics import (
    average_precision_score,
    balanced_accuracy_score,
    brier_score_loss,
    f1_score,
    log_loss,
    matthews_corrcoef,
    mean_absolute_error,
    mean_squared_error,
    precision_recall_curve,
    r2_score,
    roc_auc_score,
    roc_curve,
)
from sklearn.model_selection import GroupShuffleSplit


os.chdir(os.path.split(os.path.realpath(__file__))[0])


MP20_ROOT = Path("./MP-20-Charge")
OUTPUT_ROOT = Path("./results_mp20_periodic_soed_v13")
CACHE_ROOT = Path("./cache_mp20_periodic_soed")
TARGETS = ("band_gap",)
SCRIPT_BUILD = "v13"

RANDOM_SEED = 142
TRAIN_RATIO = 0.8
VALID_RATIO = 0.1
TEST_RATIO = 0.1
MAX_SAMPLES = 0
SPLIT_GROUP_BY_REDUCED_FORMULA = True
SPLIT_SEARCH_CANDIDATES = 256
DESCRIPTOR_WORKERS = max(1, min(8, (os.cpu_count() or 2) - 1))
DESCRIPTOR_CHUNK_SIZE = 64
FORCE_RECOMPUTE_STRUCTURES = False
REUSE_PERSISTENT_CACHE = True
CACHE_SCHEMA_VERSION = 5
REUSE_DESCRIPTOR_CACHE = True
FORCE_RECOMPUTE_FEATURES = False
DESCRIPTOR_CACHE_SCHEMA_VERSION = 4
DESCRIPTOR_CACHE_MMAP = True
EXCLUDE_UNDEFINED_PAULING_ELEMENTS = True
UNDEFINED_PAULING_ATOMIC_NUMBERS = (2, 10, 18, 36, 54, 86, 118)
UNDEFINED_PAULING_SYMBOLS = ("He", "Ne", "Ar", "Kr", "Xe", "Rn", "Og")

BAND_GAP_ZERO_THRESHOLD_EV = 0.01
BAND_GAP_SENSITIVITY_THRESHOLDS_EV = (0.0, 0.01, 0.05, 0.10)
BAND_GAP_STRATIFICATION_EDGES_EV = (-np.inf, 0.01, 0.5, 1.0, 2.0, 3.0, 6.0, np.inf)
BAND_GAP_STRATIFICATION_LABELS = (
    "metal_or_zero",
    "0.01-0.5",
    "0.5-1",
    "1-2",
    "2-3",
    "3-6",
    ">6",
)

RUN_DIRECT_REGRESSION = True
RUN_HURDLE_AUXILIARY = True
RUN_SOFT_GATE = True
RUN_HARD_GATE = True
CLIP_NEGATIVE_GAP_PREDICTIONS = True
USE_BALANCED_CLASS_WEIGHTS = False
CLASSIFIER_THRESHOLD_METRIC = "mcc"
CLASSIFIER_THRESHOLD_GRID = tuple(np.linspace(0.05, 0.95, 181))

USE_TAIL_SAMPLE_WEIGHTS = True
TAIL_THRESHOLD_EV = 3.0
TAIL_FULL_WEIGHT_EV = 6.0
TAIL_WEIGHT_MULTIPLIER = 2.5
TAIL_METRIC_BINS_EV = (-np.inf, 0.01, 0.5, 1.0, 2.0, 3.0, 6.0, np.inf)
TAIL_METRIC_LABELS = ("<=0.01", "0.01-0.5", "0.5-1", "1-2", "2-3", "3-6", ">6")
RUN_TAIL_WEIGHT_ABLATION = True
USE_MIDGAP_EXPERT = True
MIDGAP_LOWER_EV = 1.0
MIDGAP_UPPER_EV = 2.0
MIDGAP_CENTER_EV = 1.5
MIDGAP_WIDTH_EV = 0.45
MIDGAP_WEIGHT_MULTIPLIER = 2.0
MIDGAP_OBJECTIVE_WEIGHT = 0.10
CATASTROPHIC_ERROR_THRESHOLD_EV = 1.25
CATASTROPHIC_OBJECTIVE_WEIGHT = 0.05
EXPERT_BLEND_L2 = 2.0e-3

SOED_R_CUT = 5.0
SOED_N_MAX = 6
SOED_L_MAX = 4
SOED_N_RADIAL = 20
SOED_N_ANGULAR = 74
SOAP_ALPHA_GRID = (0.1, 0.2, 0.4, 0.8, 1.6, 3.2, 6.4)
SOAP_INCLUDE_LOCAL_STD = True
SOAP_ATOM_DEPOSITION = "cloud_in_cell"
SOAP_REPRESENTATIONS = (
    "psoap_scaffold_alpha_0p1",
    "psoap_scaffold_alpha_0p2",
    "psoap_scaffold_alpha_0p4",
    "psoap_scaffold_alpha_0p8",
    "psoap_scaffold_alpha_1p6",
    "psoap_scaffold_alpha_3p2",
    "psoap_scaffold_alpha_6p4",
)
SOAP_ALPHA_BY_NAME = dict(zip(SOAP_REPRESENTATIONS, SOAP_ALPHA_GRID))
DEFAULT_SOAP_NAME = "psoap_scaffold_alpha_0p4"
PRIMARY_SOAP_NAME = "psoap_validation_selected"
COMPOSITION_ONLY_NAME = "composition_only"
SOAP_COMPOSITION_REPRESENTATIONS = tuple(
    f"{name}_plus_composition" for name in SOAP_REPRESENTATIONS
)
SOAP_COMPOSITION_SOURCE = dict(
    zip(SOAP_COMPOSITION_REPRESENTATIONS, SOAP_REPRESENTATIONS)
)
PRIMARY_CHEMICAL_SOAP_NAME = "psoap_validation_selected_plus_composition"
DENSITY_COMPOSITION_NAME = "psoed_density_plus_composition"
SOAP_LOCAL_CHEMICAL_REPRESENTATIONS = tuple(
    f"{name}_local_chemical" for name in SOAP_REPRESENTATIONS
)
SOAP_LOCAL_CHEMICAL_SOURCE = dict(
    zip(SOAP_LOCAL_CHEMICAL_REPRESENTATIONS, SOAP_REPRESENTATIONS)
)
SOAP_LOCAL_CHEMICAL_BY_SOURCE = {
    source: name for name, source in SOAP_LOCAL_CHEMICAL_SOURCE.items()
}
PRIMARY_LOCAL_CHEMICAL_SOAP_NAME = "psoap_validation_selected_local_chemical"
SOED_INTERPOLATION_ORDER = 1
SOED_DENSITY_NORMALIZATION = "mean"
SOED_CONTRAST_SMOOTHING_FRACTION = 1.0 / 32.0
SOED_RECIPROCAL_GRID = 16
SOED_RECIPROCAL_BINS = 24
SOED_CHEMICAL_PROJECTION_SIZE = 128
SOED_BLOCK_PROJECTION_SIZE = 128
SOED_ELECTRON_NUCLEAR_ALPHA = 1.6
SOED_CHANNEL_SETS = {
    "psoed_density": ("density",),
    "psoed_log_density": ("log_density",),
    "psoed_density_contrast": ("density_contrast",),
    "psoed_density_gradient": ("density", "gradient"),
    "psoed_physics_multichannel": (
        "density",
        "log_density",
        "density_contrast",
        "gradient",
        "laplacian",
    ),
    "psoed_chemical_multichannel": (
        "density",
        "log_density",
        "density_contrast",
        "gradient",
        "laplacian",
    ),
    "psoed_density_local_chemical": ("density",),
    "psoed_density_local_chemical_global": ("density",),
    "psoed_density_local_chemical_radial": ("density",),
    "psoed_density_local_chemical_radial_global": ("density",),
    "psoed_compact_chemical_multichannel": (
        "density",
        "log_density",
        "density_contrast",
        "gradient",
        "laplacian",
        "electron_nuclear_contrast",
    ),
}
SOED_CHEMICAL_DESCRIPTORS = (
    "psoed_chemical_multichannel",
    "psoed_compact_chemical_multichannel",
)
SOED_LOCAL_CHEMICAL_DESCRIPTORS = (
    "psoed_density_local_chemical",
    "psoed_density_local_chemical_global",
    "psoed_density_local_chemical_radial",
    "psoed_density_local_chemical_radial_global",
)
SOED_RADIAL_ELECTRONIC_DESCRIPTORS = (
    "psoed_density_local_chemical_radial",
    "psoed_density_local_chemical_radial_global",
)
SOED_COMPACT_GLOBAL_DESCRIPTORS = (
    "psoed_density_local_chemical_global",
    "psoed_density_local_chemical_radial_global",
)
SOED_SEPARABLE_DESCRIPTORS = ("psoed_compact_chemical_multichannel",)
SOED_CHEMICAL_PROPERTIES = (
    "atomic_number",
    "period",
    "group",
    "electronegativity",
    "covalent_radius",
    "atomic_mass",
    "first_ionization_energy",
    "electron_affinity",
    "valence_s",
    "valence_p",
    "valence_d",
    "valence_f",
)
MATCHED_SOED_NAME = "psoed_density"
MATCHED_CHEMICAL_SOED_NAME = "psoed_density_local_chemical"
LEGACY_CHEMICAL_SOED_NAME = "psoed_chemical_multichannel"
COMPACT_MULTICHANNEL_SOED_NAME = "psoed_compact_chemical_multichannel"
SOED_ENHANCED_CANDIDATES = SOED_LOCAL_CHEMICAL_DESCRIPTORS
PRIMARY_SOED_NAME = "psoed_validation_selected_electronic"
SOED_CANDIDATE_LABELS = {
    "psoed_density_local_chemical": r"Electronic SOED ($\rho$)",
    "psoed_density_local_chemical_global": r"Electronic SOED ($\rho$ + global)",
    "psoed_density_local_chemical_radial": r"Electronic SOED ($\rho$ + radial charge)",
    "psoed_density_local_chemical_radial_global": r"Electronic SOED ($\rho$ + radial + global)",
}
TAIL_ABLATION_REPRESENTATIONS = (
    PRIMARY_LOCAL_CHEMICAL_SOAP_NAME,
    *SOED_ENHANCED_CANDIDATES,
)
MIDGAP_EXPERT_REPRESENTATIONS = (
    *SOAP_LOCAL_CHEMICAL_REPRESENTATIONS,
    *SOED_ENHANCED_CANDIDATES,
)
HURDLE_REPRESENTATIONS = SOED_ENHANCED_CANDIDATES
MAIN_REPRESENTATIONS = (
    COMPOSITION_ONLY_NAME,
    PRIMARY_SOAP_NAME,
    MATCHED_SOED_NAME,
    PRIMARY_LOCAL_CHEMICAL_SOAP_NAME,
    PRIMARY_SOED_NAME,
)

MODEL_NAME = "xgboost"
PREFER_GPU = True
GPU_ID = 0
EARLY_STOPPING_ROUNDS = 200
USE_VALIDATION_WEIGHT_BLEND = True
WEIGHT_BLEND_GRID = tuple(np.linspace(0.0, 1.0, 41))
USE_SOED_CANDIDATE_ENSEMBLE = True
SOED_ENSEMBLE_L2 = 1.0e-2
SOED_ENSEMBLE_MIN_WEIGHT = 0.0
XGBOOST_PARAMS = {
    "n_estimators": 12000,
    "max_depth": 3,
    "learning_rate": 0.02,
    "min_child_weight": 25.0,
    "subsample": 0.75,
    "colsample_bytree": 0.55,
    "reg_alpha": 0.80,
    "reg_lambda": 40.0,
    "gamma": 0.01,
    "max_bin": 256,
}

ENABLE_OPTUNA = False
OPTUNA_TRIALS = 30
OPTUNA_TIMEOUT_SECONDS = None
OPTUNA_REPRESENTATIONS = (
    DEFAULT_SOAP_NAME,
    *SOED_ENHANCED_CANDIDATES,
)
OPTUNA_MAX_TRAIN_SAMPLES = 20000

BOOTSTRAP_REPEATS = 10000
BOOTSTRAP_CONFIDENCE = 0.95
RUN_REPEATED_SPLIT_ROBUSTNESS = True
CONFIRMATORY_PROTOCOL_FROZEN = True
EXPLORATORY_LEGACY_SEED = 42
ROBUSTNESS_SPLIT_SEEDS = (142, 242, 342, 442, 542)

PLOT_DPI = 600
FONT_FAMILY_PRIORITY = (
    "Arial",
    "Microsoft YaHei",
    "Noto Sans CJK SC",
    "SimHei",
    "DejaVu Sans",
)
TITLE_FONT_SIZE = 14
AXIS_LABEL_FONT_SIZE = 12
TICK_LABEL_FONT_SIZE = 10
LEGEND_FONT_SIZE = 10
ANNOTATION_FONT_SIZE = 9
LINE_WIDTH = 1.8
MARKER_SIZE = 18
TAB_COLORS = {
    "train": "tab:blue",
    "validation": "tab:orange",
    "test": "tab:green",
    "soap": "tab:purple",
    "soed": "tab:brown",
    "direct": "tab:blue",
    "soft_gate": "tab:cyan",
    "hard_gate": "tab:red",
    "reference": "tab:gray",
}
TAB_PALETTE = (
    "tab:blue",
    "tab:orange",
    "tab:green",
    "tab:red",
    "tab:purple",
    "tab:brown",
    "tab:pink",
    "tab:gray",
    "tab:olive",
    "tab:cyan",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=MP20_ROOT)
    parser.add_argument("--output", type=Path, default=OUTPUT_ROOT)
    parser.add_argument("--cache", type=Path, default=CACHE_ROOT)
    parser.add_argument("--limit", type=int, default=MAX_SAMPLES)
    parser.add_argument("--workflow-only", action="store_true")
    parser.add_argument("--force-structures", action="store_true")
    parser.add_argument("--force-features", action="store_true")
    return parser.parse_args()


def ensure_dir(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    return path


def save_dat(frame: pd.DataFrame, path: Path) -> None:
    ensure_dir(path.parent)
    frame.to_csv(path, sep="\t", index=False, float_format="%.10g")


def json_default(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, (np.integer, np.floating)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, tuple):
        return list(value)
    return str(value)


def setup_logger(output: Path) -> logging.Logger:
    ensure_dir(output / "logs")
    logger = logging.getLogger("mp20-periodic-soed-xgb")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    file_handler = logging.FileHandler(
        output / "logs" / "run.log", mode="w", encoding="utf-8"
    )
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def seed_everything(seed: int) -> None:
    os.environ["PYTHONHASHSEED"] = str(seed)
    random.seed(seed)
    np.random.seed(seed)


def detect_device(logger: logging.Logger) -> dict[str, Any]:
    info: dict[str, Any] = {
        "platform": platform.platform(),
        "python": sys.version.replace("\n", " "),
        "cpu_count": os.cpu_count(),
        "gpu_requested": PREFER_GPU,
        "gpu_available": False,
        "gpu_name": None,
    }
    if PREFER_GPU and shutil.which("nvidia-smi"):
        try:
            command = [
                "nvidia-smi",
                f"--id={GPU_ID}",
                "--query-gpu=name,memory.total,memory.free",
                "--format=csv,noheader,nounits",
            ]
            result = subprocess.run(
                command, capture_output=True, text=True, check=True, timeout=10
            )
            fields = [value.strip() for value in result.stdout.strip().split(",")]
            info["gpu_available"] = bool(fields and fields[0])
            info["gpu_name"] = fields[0] if fields else None
            if len(fields) >= 3:
                info["gpu_memory_total_mb"] = float(fields[1])
                info["gpu_memory_free_mb"] = float(fields[2])
        except Exception as exc:
            info["gpu_detection_error"] = repr(exc)
    try:
        import xgboost as xgb

        info["xgboost_version"] = xgb.__version__
        info["xgboost_build_info"] = xgb.build_info()
    except Exception as exc:
        info["xgboost_detection_error"] = repr(exc)
    logger.info(
        "Hardware: %s", json.dumps(info, ensure_ascii=False, default=json_default)
    )
    return info


def config_snapshot() -> dict[str, Any]:
    return {
        "script_build": SCRIPT_BUILD,
        "dataset": {
            "targets": TARGETS,
            "seed": RANDOM_SEED,
            "ratios": [TRAIN_RATIO, VALID_RATIO, TEST_RATIO],
            "max_samples": MAX_SAMPLES,
            "group_by_reduced_formula": SPLIT_GROUP_BY_REDUCED_FORMULA,
            "split_search_candidates": SPLIT_SEARCH_CANDIDATES,
            "zero_threshold_ev": BAND_GAP_ZERO_THRESHOLD_EV,
            "tail_threshold_ev": TAIL_THRESHOLD_EV,
            "tail_weight_multiplier": TAIL_WEIGHT_MULTIPLIER,
            "confirmatory_protocol_frozen": CONFIRMATORY_PROTOCOL_FROZEN,
            "exploratory_legacy_seed_not_reused": EXPLORATORY_LEGACY_SEED,
            "exclude_undefined_pauling_elements": EXCLUDE_UNDEFINED_PAULING_ELEMENTS,
            "undefined_pauling_atomic_numbers": UNDEFINED_PAULING_ATOMIC_NUMBERS,
            "undefined_pauling_symbols": UNDEFINED_PAULING_SYMBOLS,
        },
        "soap": {
            "implementation": "periodic Gaussian scaffold density with the shared SOED invariant pipeline",
            "reference_equations": "Gugler-Reiher JCTC 2022, Eqs. 13-14, 26, 38, and 40",
            "alpha_grid_angstrom_minus_2": SOAP_ALPHA_GRID,
            "alpha_selection": "minimum validation RMSE after validation-only weighting selection",
            "periodization": "reciprocal_space_Gaussian_convolution",
            "atom_deposition": SOAP_ATOM_DEPOSITION,
            "r_cut": SOED_R_CUT,
            "n_max": SOED_N_MAX,
            "l_max": SOED_L_MAX,
            "n_radial": SOED_N_RADIAL,
            "n_angular": SOED_N_ANGULAR,
            "local_pooling": "mean_and_std" if SOAP_INCLUDE_LOCAL_STD else "mean",
            "primary": PRIMARY_SOAP_NAME,
            "local_chemical_primary": PRIMARY_LOCAL_CHEMICAL_SOAP_NAME,
        },
        "chemistry_fair_controls": {
            "composition_dimension": PeriodicSOEDFamily.composition_feature_count,
            "composition_only": COMPOSITION_ONLY_NAME,
            "soap_plus_composition": PRIMARY_CHEMICAL_SOAP_NAME,
            "density_soed_plus_composition": DENSITY_COMPOSITION_NAME,
            "local_chemical_soap": PRIMARY_LOCAL_CHEMICAL_SOAP_NAME,
            "local_chemical_density_soed": MATCHED_CHEMICAL_SOED_NAME,
            "shared_local_chemical_projection": {
                "properties": SOED_CHEMICAL_PROPERTIES,
                "property_projection_size": SOED_CHEMICAL_PROJECTION_SIZE,
                "orbital_block_projection_size": SOED_BLOCK_PROJECTION_SIZE,
                "composition_dimension": PeriodicSOEDFamily.composition_feature_count,
            },
            "selection": "SOAP alpha and tail weighting selected by validation RMSE",
            "soap_soed_fusion": False,
        },
        "soed": {
            "r_cut": SOED_R_CUT,
            "n_max": SOED_N_MAX,
            "l_max": SOED_L_MAX,
            "n_radial": SOED_N_RADIAL,
            "n_angular": SOED_N_ANGULAR,
            "interpolation_order": SOED_INTERPOLATION_ORDER,
            "density_normalization": SOED_DENSITY_NORMALIZATION,
            "contrast_smoothing_fraction": SOED_CONTRAST_SMOOTHING_FRACTION,
            "reciprocal_grid": SOED_RECIPROCAL_GRID,
            "reciprocal_bins": SOED_RECIPROCAL_BINS,
            "chemical_projection_size": SOED_CHEMICAL_PROJECTION_SIZE,
            "block_projection_size": SOED_BLOCK_PROJECTION_SIZE,
            "channel_sets": SOED_CHANNEL_SETS,
            "chemical_descriptors": SOED_CHEMICAL_DESCRIPTORS,
            "local_chemical_descriptors": SOED_LOCAL_CHEMICAL_DESCRIPTORS,
            "radial_electronic_descriptors": SOED_RADIAL_ELECTRONIC_DESCRIPTORS,
            "compact_global_descriptors": SOED_COMPACT_GLOBAL_DESCRIPTORS,
            "enhanced_candidates": SOED_ENHANCED_CANDIDATES,
            "candidate_prediction_ensemble": USE_SOED_CANDIDATE_ENSEMBLE,
            "candidate_ensemble_l2": SOED_ENSEMBLE_L2,
            "electronic_radial_feature_count": PeriodicSOEDFamily.electronic_radial_feature_count,
            "compact_global_feature_count": PeriodicSOEDFamily.compact_global_feature_count,
            "separable_descriptors": SOED_SEPARABLE_DESCRIPTORS,
            "electron_nuclear_alpha": SOED_ELECTRON_NUCLEAR_ALPHA,
            "chemical_properties": SOED_CHEMICAL_PROPERTIES,
            "primary": PRIMARY_SOED_NAME,
        },
        "model": {
            "name": MODEL_NAME,
            "weighting_selection": "validation-only convex blend of unweighted, high-gap, and predefined 1-2 eV experts for fair SOAP and SOED candidates",
            "weight_blend_enabled": USE_VALIDATION_WEIGHT_BLEND,
            "weight_blend_optimizer": "SLSQP nonnegative simplex with L2 regularization",
            "midgap_expert": {
                "enabled": USE_MIDGAP_EXPERT,
                "representations": MIDGAP_EXPERT_REPRESENTATIONS,
                "interval_ev": [MIDGAP_LOWER_EV, MIDGAP_UPPER_EV],
                "center_ev": MIDGAP_CENTER_EV,
                "width_ev": MIDGAP_WIDTH_EV,
                "sample_weight_multiplier": MIDGAP_WEIGHT_MULTIPLIER,
                "validation_objective_weight": MIDGAP_OBJECTIVE_WEIGHT,
                "catastrophic_error_threshold_ev": CATASTROPHIC_ERROR_THRESHOLD_EV,
                "catastrophic_objective_weight": CATASTROPHIC_OBJECTIVE_WEIGHT,
                "expert_blend_l2": EXPERT_BLEND_L2,
            },
            "prefer_gpu": PREFER_GPU,
            "gpu_id": GPU_ID,
            "early_stopping_rounds": EARLY_STOPPING_ROUNDS,
            "parameters": XGBOOST_PARAMS,
        },
        "robustness": {
            "group_bootstrap_repeats": BOOTSTRAP_REPEATS,
            "bootstrap_unit": "reduced_formula",
            "repeat_group_split": RUN_REPEATED_SPLIT_ROBUSTNESS,
            "split_seeds": ROBUSTNESS_SPLIT_SEEDS,
            "protocol": "fresh frozen formula-group outer splits; all choices use each outer training/validation partition only",
        },
        "optuna": {
            "enabled": ENABLE_OPTUNA,
            "trials": OPTUNA_TRIALS,
            "representations": OPTUNA_REPRESENTATIONS,
        },
        "statistics": {
            "bootstrap_repeats": BOOTSTRAP_REPEATS,
            "confidence": BOOTSTRAP_CONFIDENCE,
        },
        "cache": {
            "root": str(CACHE_ROOT),
            "structure_schema_version": CACHE_SCHEMA_VERSION,
            "descriptor_schema_version": DESCRIPTOR_CACHE_SCHEMA_VERSION,
            "reuse": REUSE_PERSISTENT_CACHE,
            "reuse_descriptor_cache": REUSE_DESCRIPTOR_CACHE,
            "force_recompute_structures": FORCE_RECOMPUTE_STRUCTURES,
            "force_recompute_features": FORCE_RECOMPUTE_FEATURES,
            "descriptor_mmap": DESCRIPTOR_CACHE_MMAP,
        },
        "plot": {
            "dpi": PLOT_DPI,
            "font_priority": FONT_FAMILY_PRIORITY,
            "title": TITLE_FONT_SIZE,
            "axis": AXIS_LABEL_FONT_SIZE,
            "tick": TICK_LABEL_FONT_SIZE,
            "legend": LEGEND_FONT_SIZE,
            "annotation": ANNOTATION_FONT_SIZE,
            "colors": TAB_COLORS,
        },
    }


def dataset_digest(frame: pd.DataFrame) -> str:
    digest = hashlib.sha256(
        f"schema={CACHE_SCHEMA_VERSION};n={len(frame)}".encode("utf-8")
    )
    for material_id, cif in frame[["material_id", "cif"]].itertuples(
        index=False, name=None
    ):
        digest.update(str(material_id).encode("utf-8"))
        digest.update(hashlib.sha256(str(cif).encode("utf-8")).digest())
    return digest.hexdigest()[:16]


def canonical_composition_key(structure: Any) -> str:
    from pymatgen.core import Element

    pairs: list[tuple[int, Fraction]] = []
    for symbol, amount in structure.composition.get_el_amt_dict().items():
        fraction = Fraction(str(round(float(amount), 8))).limit_denominator(10000)
        pairs.append((int(Element(symbol).Z), fraction))
    denominator = 1
    for _, fraction in pairs:
        denominator = math.lcm(denominator, fraction.denominator)
    integers = [int(fraction * denominator) for _, fraction in pairs]
    divisor = max(
        1, reduce(math.gcd, (abs(value) for value in integers if value != 0), 0)
    )
    return "|".join(
        f"{atomic_number}:{integer // divisor}"
        for (atomic_number, _), integer in sorted(
            zip(pairs, integers), key=lambda item: item[0][0]
        )
    )


def resolve_dataset_root(root: Path) -> Path:
    root = root.expanduser().resolve()
    for candidate in (root, root / "MP-20-Charge"):
        if (candidate / "structure").is_dir() and (
            candidate / "charge_density"
        ).is_dir():
            return candidate
    raise FileNotFoundError(f"Cannot find MP-20-Charge under {root}")


def load_metadata(root: Path, limit: int, logger: logging.Logger) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for source_split in ("train", "val", "test"):
        path = root / "structure" / f"{source_split}.csv"
        if not path.exists():
            raise FileNotFoundError(path)
        frame = pd.read_csv(path, keep_default_na=False)
        missing = {"material_id", "cif"} - set(frame.columns)
        if missing:
            raise ValueError(f"{path} missing columns: {sorted(missing)}")
        frame["source_split"] = source_split
        frames.append(frame)
    source_ids = {
        str(x["source_split"].iloc[0]): set(x["material_id"].astype(str))
        for x in frames
    }
    overlaps = {
        f"{a}_{b}": len(source_ids[a] & source_ids[b])
        for a, b in (("train", "val"), ("train", "test"), ("val", "test"))
    }
    if any(overlaps.values()):
        logger.info("Official source-split ID overlaps detected: %s", overlaps)
    frame = pd.concat(frames, ignore_index=True)
    frame["material_id"] = frame["material_id"].astype(str)
    before = len(frame)
    duplicate_rows = frame.loc[frame.duplicated("material_id", keep=False)]
    inconsistent = []
    for material_id, group in duplicate_rows.groupby("material_id", sort=False):
        target_values = pd.to_numeric(group["band_gap"], errors="coerce").to_numpy(
            dtype=float
        )
        same_target = bool(np.allclose(target_values, target_values[0], equal_nan=True))
        same_cif = (
            group["cif"]
            .astype(str)
            .map(lambda value: hashlib.sha256(value.encode("utf-8")).hexdigest())
            .nunique()
            == 1
        )
        if not (same_target and same_cif):
            inconsistent.append(str(material_id))
    if inconsistent:
        raise ValueError(f"Conflicting duplicated material IDs: {inconsistent[:10]}")
    frame = frame.drop_duplicates("material_id", keep="first").reset_index(drop=True)
    if len(frame) != before:
        logger.info(
            "Collapsed %d exact duplicated material_id rows.", before - len(frame)
        )
    for target in TARGETS:
        if target not in frame:
            raise ValueError(
                f"Missing target {target}; available columns: {list(frame.columns)}"
            )
        frame[target] = pd.to_numeric(frame[target], errors="coerce")
    exists = frame["material_id"].map(
        lambda value: (root / "charge_density" / f"{value}.npy").exists()
    )
    if not bool(exists.all()):
        logger.warning(
            "Dropping %d rows with missing density files.", int((~exists).sum())
        )
        frame = frame.loc[exists].reset_index(drop=True)
    if limit > 0:
        frame = frame.sample(
            n=min(limit, len(frame)), random_state=RANDOM_SEED
        ).reset_index(drop=True)
    logger.info("Loaded %d unique MP-20-Charge records.", len(frame))
    return frame


def parse_structure_from_cif(cif: str):
    from pymatgen.core import Structure

    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"Issues encountered while parsing CIF:.*",
            category=UserWarning,
        )
        warnings.filterwarnings(
            "ignore",
            message=r"No Pauling electronegativity for .*",
            category=UserWarning,
        )
        return Structure.from_str(str(cif), fmt="cif")


@dataclass
class StructureStore:
    cache_dir: Path
    digest: str
    material_ids: np.ndarray
    formula_groups: np.ndarray
    cells: np.ndarray
    atom_offsets: np.ndarray
    atomic_numbers: np.ndarray
    frac_positions: np.ndarray

    def atom_slice(self, index: int) -> slice:
        return slice(int(self.atom_offsets[index]), int(self.atom_offsets[index + 1]))


def prepare_structure_cache(
    frame: pd.DataFrame,
    cache_root: Path,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, StructureStore]:
    digest = dataset_digest(frame)
    cache_dir = ensure_dir(
        cache_root.expanduser().resolve() / f"dataset_{digest}" / "structures"
    )
    manifest_path = cache_dir / "manifest.json"
    index_path = cache_dir / "structure_index.csv"
    array_paths = {
        "cells": cache_dir / "cells.npy",
        "atom_offsets": cache_dir / "atom_offsets.npy",
        "atomic_numbers": cache_dir / "atomic_numbers.npy",
        "frac_positions": cache_dir / "frac_positions.npy",
    }
    cache_ready = (
        manifest_path.exists()
        and index_path.exists()
        and all(path.exists() for path in array_paths.values())
    )
    if cache_ready and not FORCE_RECOMPUTE_STRUCTURES:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        index_frame = pd.read_csv(
            index_path, dtype={"material_id": str, "formula_group": str}
        )
        expected_ids = frame["material_id"].astype(str).to_numpy()
        cached_ids = index_frame["material_id"].astype(str).to_numpy()
        cache_ready = (
            manifest.get("dataset_digest") == digest
            and int(manifest.get("n_structures", -1)) == len(frame)
            and np.array_equal(expected_ids, cached_ids)
        )
    if not cache_ready or FORCE_RECOMPUTE_STRUCTURES:
        cells: list[np.ndarray] = []
        numbers: list[np.ndarray] = []
        positions: list[np.ndarray] = []
        offsets = [0]
        index_rows = []
        failures = []
        for index, row in frame.iterrows():
            try:
                structure = parse_structure_from_cif(str(row["cif"]))
                atomic_numbers = np.asarray(
                    [site.specie.Z for site in structure], dtype=np.int16
                )
                fractional = np.mod(
                    np.asarray(structure.frac_coords, dtype=np.float32), 1.0
                )
                cells.append(np.asarray(structure.lattice.matrix, dtype=np.float64))
                numbers.append(atomic_numbers)
                positions.append(fractional)
                offsets.append(offsets[-1] + len(atomic_numbers))
                index_rows.append(
                    {
                        "material_id": str(row["material_id"]),
                        "formula_group": canonical_composition_key(structure),
                        "n_atoms": len(atomic_numbers),
                    }
                )
            except Exception as exc:
                failures.append((str(row["material_id"]), repr(exc)))
                break
            if (index + 1) % 5000 == 0:
                logger.info("Structure cache progress: %d/%d", index + 1, len(frame))
        if failures:
            raise RuntimeError(f"Structure parsing failed: {failures[:3]}")
        np.save(array_paths["cells"], np.asarray(cells, dtype=np.float64))
        np.save(array_paths["atom_offsets"], np.asarray(offsets, dtype=np.int64))
        np.save(
            array_paths["atomic_numbers"],
            np.concatenate(numbers).astype(np.int16, copy=False),
        )
        np.save(
            array_paths["frac_positions"],
            np.concatenate(positions).astype(np.float32, copy=False),
        )
        index_frame = pd.DataFrame(index_rows)
        index_frame.to_csv(index_path, index=False)
        manifest = {
            "schema_version": CACHE_SCHEMA_VERSION,
            "dataset_digest": digest,
            "n_structures": len(frame),
            "n_atoms": int(offsets[-1]),
            "material_id_sha256": hashlib.sha256(
                "\n".join(frame["material_id"].astype(str)).encode("utf-8")
            ).hexdigest(),
        }
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
        logger.info("Structure cache created: %s", cache_dir)
    else:
        logger.info("Structure cache hit: %s (%d structures)", cache_dir, len(frame))
    index_frame = pd.read_csv(
        index_path, dtype={"material_id": str, "formula_group": str}
    )
    result = frame.copy()
    result["structure_index"] = np.arange(len(result), dtype=np.int64)
    result["reduced_formula"] = index_frame["formula_group"].astype(str).to_numpy()
    store = StructureStore(
        cache_dir=cache_dir,
        digest=digest,
        material_ids=index_frame["material_id"].astype(str).to_numpy(),
        formula_groups=index_frame["formula_group"].astype(str).to_numpy(),
        cells=np.load(array_paths["cells"], mmap_mode="r"),
        atom_offsets=np.load(array_paths["atom_offsets"], mmap_mode="r"),
        atomic_numbers=np.load(array_paths["atomic_numbers"], mmap_mode="r"),
        frac_positions=np.load(array_paths["frac_positions"], mmap_mode="r"),
    )
    return result, store


def exclude_undefined_pauling_structures(
    frame: pd.DataFrame,
    store: StructureStore,
    output: Path,
    logger: logging.Logger,
) -> pd.DataFrame:
    from pymatgen.core import Element

    excluded_numbers = set(map(int, UNDEFINED_PAULING_ATOMIC_NUMBERS))
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message=r"No Pauling electronegativity for .*",
            category=UserWarning,
        )
        for atomic_number in np.unique(np.asarray(store.atomic_numbers, dtype=int)):
            try:
                value = float(Element.from_Z(int(atomic_number)).X)
            except (AttributeError, KeyError, TypeError, ValueError):
                value = math.nan
            if not np.isfinite(value):
                excluded_numbers.add(int(atomic_number))
    checked_symbols = tuple(
        Element.from_Z(value).symbol for value in sorted(excluded_numbers)
    )
    excluded_rows: list[dict[str, Any]] = []
    keep = np.ones(len(frame), dtype=bool)
    for row_position, row in enumerate(frame.itertuples(index=False)):
        structure_index = int(getattr(row, "structure_index"))
        atom_slice = store.atom_slice(structure_index)
        present = sorted(
            excluded_numbers.intersection(map(int, store.atomic_numbers[atom_slice]))
        )
        if not present:
            continue
        keep[row_position] = False
        excluded_rows.append(
            {
                "material_id": str(getattr(row, "material_id")),
                "excluded_elements": ",".join(
                    Element.from_Z(value).symbol for value in present
                ),
                "excluded_atomic_numbers": ",".join(map(str, present)),
                "band_gap": float(getattr(row, "band_gap")),
                "source_split": str(getattr(row, "source_split")),
                "reason": "undefined_Pauling_electronegativity",
            }
        )
    quality_dir = ensure_dir(output / "data_quality")
    excluded = pd.DataFrame(
        excluded_rows,
        columns=(
            "material_id",
            "excluded_elements",
            "excluded_atomic_numbers",
            "band_gap",
            "source_split",
            "reason",
        ),
    )
    excluded.to_csv(
        quality_dir / "excluded_undefined_pauling_elements.csv", index=False
    )
    save_dat(excluded, quality_dir / "excluded_undefined_pauling_elements.dat")
    element_counts: dict[str, int] = {}
    for values in excluded.get("excluded_elements", pd.Series(dtype=str)).astype(str):
        for symbol in filter(None, values.split(",")):
            element_counts[symbol] = element_counts.get(symbol, 0) + 1
    summary = pd.DataFrame(
        [
            {
                "n_before": len(frame),
                "n_excluded": int((~keep).sum()),
                "n_retained": int(keep.sum()),
                "excluded_fraction": float((~keep).mean()),
                "elements_checked": ",".join(checked_symbols),
                "element_structure_counts": json.dumps(element_counts, sort_keys=True),
                "raw_files_deleted": False,
            }
        ]
    )
    summary.to_csv(quality_dir / "undefined_pauling_filter_summary.csv", index=False)
    save_dat(summary, quality_dir / "undefined_pauling_filter_summary.dat")
    logger.info(
        "Undefined-Pauling filter: checked=%s before=%d excluded=%d retained=%d counts=%s; raw structure and density files were not deleted.",
        checked_symbols,
        len(frame),
        int((~keep).sum()),
        int(keep.sum()),
        element_counts,
    )
    if not EXCLUDE_UNDEFINED_PAULING_ELEMENTS:
        logger.warning(
            "Undefined-Pauling exclusion is disabled; %d affected structures are retained.",
            int((~keep).sum()),
        )
        return frame.reset_index(drop=True)
    return frame.loc[keep].reset_index(drop=True)


def add_split_metadata(frame: pd.DataFrame, logger: logging.Logger) -> pd.DataFrame:
    result = frame.copy()
    if "reduced_formula" not in result:
        raise ValueError("Structure cache metadata is missing reduced_formula")
    gap = pd.to_numeric(result["band_gap"], errors="coerce")
    result["band_gap_stratum"] = pd.cut(
        gap,
        bins=BAND_GAP_STRATIFICATION_EDGES_EV,
        labels=BAND_GAP_STRATIFICATION_LABELS,
        include_lowest=True,
        right=True,
    ).astype("object")
    result.loc[gap.isna(), "band_gap_stratum"] = "missing"
    result["band_gap_stratum"] = result["band_gap_stratum"].astype(str)
    result["is_nonmetal"] = (gap > BAND_GAP_ZERO_THRESHOLD_EV).astype(np.int8)
    logger.info(
        "Split metadata: groups=%d failures=%d strata=%s",
        result["reduced_formula"].nunique(),
        0,
        result["band_gap_stratum"].value_counts(dropna=False).to_dict(),
    )
    return result


def distribution_distance(
    all_strata: np.ndarray,
    held_strata: np.ndarray,
    actual_fraction: float,
    requested_fraction: float,
) -> float:
    labels = sorted(set(all_strata.tolist()))
    overall = (
        pd.Series(all_strata)
        .value_counts(normalize=True)
        .reindex(labels, fill_value=0.0)
    )
    held = (
        pd.Series(held_strata)
        .value_counts(normalize=True)
        .reindex(labels, fill_value=0.0)
    )
    missing = float(((overall > 0.0) & (held == 0.0)).sum())
    return (
        8.0 * abs(actual_fraction - requested_fraction)
        + float(np.mean(np.abs(overall - held)))
        + 0.25 * missing
    )


def best_group_holdout(
    indices: np.ndarray,
    groups: np.ndarray,
    strata: np.ndarray,
    holdout_fraction: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    splitter = GroupShuffleSplit(
        n_splits=SPLIT_SEARCH_CANDIDATES,
        test_size=holdout_fraction,
        random_state=seed,
    )
    best_score = math.inf
    best_pair: tuple[np.ndarray, np.ndarray] | None = None
    for retained, held in splitter.split(np.zeros(len(indices)), strata, groups):
        score = distribution_distance(
            strata, strata[held], len(held) / len(indices), holdout_fraction
        )
        if score < best_score:
            best_score = score
            best_pair = (indices[retained], indices[held])
    if best_pair is None:
        raise RuntimeError("Could not construct group-disjoint split")
    return best_pair


def save_split_diagnostics(result: pd.DataFrame, output: Path) -> None:
    split_dir = ensure_dir(output / "splits")
    rows = []
    for split in ("train", "validation", "test"):
        subset = result.loc[result["split"] == split]
        gap = subset["band_gap"].to_numpy(dtype=float)
        rows.append(
            {
                "split": split,
                "n_samples": len(subset),
                "fraction": len(subset) / len(result),
                "n_reduced_formulas": subset["reduced_formula"].nunique(),
                "n_nonmetal": int((gap > BAND_GAP_ZERO_THRESHOLD_EV).sum()),
                "nonmetal_fraction": float(np.mean(gap > BAND_GAP_ZERO_THRESHOLD_EV)),
                "gap_mean_ev": float(np.nanmean(gap)),
                "gap_median_ev": float(np.nanmedian(gap)),
                "gap_max_ev": float(np.nanmax(gap)),
            }
        )
    summary = pd.DataFrame(rows)
    summary.to_csv(split_dir / "split_summary.csv", index=False)
    save_dat(summary, split_dir / "split_summary.dat")
    strata = pd.crosstab(
        result["band_gap_stratum"], result["split"], margins=True
    ).reset_index()
    strata.to_csv(split_dir / "band_gap_strata_counts.csv", index=False)
    save_dat(strata, split_dir / "band_gap_strata_counts.dat")
    sensitivity = []
    for threshold in BAND_GAP_SENSITIVITY_THRESHOLDS_EV:
        for split in ("train", "validation", "test", "all"):
            subset = result if split == "all" else result.loc[result["split"] == split]
            gap = subset["band_gap"].to_numpy(dtype=float)
            sensitivity.append(
                {
                    "threshold_ev": threshold,
                    "split": split,
                    "n_samples": len(subset),
                    "n_zero_or_metal": int((gap <= threshold).sum()),
                    "n_positive_or_nonmetal": int((gap > threshold).sum()),
                    "positive_fraction": float(np.mean(gap > threshold)),
                }
            )
    sensitivity_frame = pd.DataFrame(sensitivity)
    sensitivity_frame.to_csv(
        split_dir / "band_gap_threshold_sensitivity.csv", index=False
    )
    save_dat(sensitivity_frame, split_dir / "band_gap_threshold_sensitivity.dat")


def assign_group_stratified_splits(
    frame: pd.DataFrame, output: Path, logger: logging.Logger
) -> pd.DataFrame:
    if not np.isclose(TRAIN_RATIO + VALID_RATIO + TEST_RATIO, 1.0):
        raise ValueError("Split ratios must sum to one")
    result = add_split_metadata(frame, logger)
    indices = np.arange(len(result))
    groups = (
        result["reduced_formula"].to_numpy(dtype=str)
        if SPLIT_GROUP_BY_REDUCED_FORMULA
        else result["material_id"].to_numpy(dtype=str)
    )
    strata = result["band_gap_stratum"].to_numpy(dtype=str)
    train_valid, test = best_group_holdout(
        indices, groups, strata, TEST_RATIO, RANDOM_SEED
    )
    train, valid = best_group_holdout(
        train_valid,
        groups[train_valid],
        strata[train_valid],
        VALID_RATIO / (TRAIN_RATIO + VALID_RATIO),
        RANDOM_SEED + 1,
    )
    result["split"] = ""
    result.loc[train, "split"] = "train"
    result.loc[valid, "split"] = "validation"
    result.loc[test, "split"] = "test"
    keep = [
        "material_id",
        "source_split",
        "split",
        "reduced_formula",
        "band_gap_stratum",
        "is_nonmetal",
        *TARGETS,
    ]
    split_dir = ensure_dir(output / "splits")
    result[keep].to_csv(split_dir / "all_splits.csv", index=False)
    save_dat(result[keep], split_dir / "all_splits.dat")
    for split in ("train", "validation", "test"):
        subset = result.loc[result["split"] == split, keep]
        subset.to_csv(split_dir / f"{split}.csv", index=False)
        save_dat(subset, split_dir / f"{split}.dat")
        positive = subset.loc[subset["band_gap"] > BAND_GAP_ZERO_THRESHOLD_EV]
        positive.to_csv(split_dir / f"{split}_positive_gap.csv", index=False)
        save_dat(positive, split_dir / f"{split}_positive_gap.dat")
    group_sets = {
        split: set(result.loc[result["split"] == split, "reduced_formula"])
        for split in ("train", "validation", "test")
    }
    overlaps = {
        "train_validation": len(group_sets["train"] & group_sets["validation"]),
        "train_test": len(group_sets["train"] & group_sets["test"]),
        "validation_test": len(group_sets["validation"] & group_sets["test"]),
    }
    if SPLIT_GROUP_BY_REDUCED_FORMULA and any(overlaps.values()):
        raise RuntimeError(f"Reduced-formula leakage: {overlaps}")
    save_split_diagnostics(result, output)
    logger.info(
        "Grouped split seed=%d counts=%s ratios=%.2f:%.2f:%.2f overlaps=%s",
        RANDOM_SEED,
        result["split"].value_counts().to_dict(),
        TRAIN_RATIO,
        VALID_RATIO,
        TEST_RATIO,
        overlaps,
    )
    return result


@dataclass
class PeriodicDensity:
    cell: np.ndarray
    frac_positions: np.ndarray
    atomic_numbers: np.ndarray
    density: np.ndarray
    identifier: str

    @property
    def volume(self) -> float:
        return float(abs(np.linalg.det(self.cell)))

    @property
    def n_atoms(self) -> int:
        return int(len(self.atomic_numbers))


def row_to_sample(
    index: int, row: pd.Series, root: Path, store: StructureStore
) -> PeriodicDensity:
    structure_index = int(row["structure_index"])
    atom_slice = store.atom_slice(structure_index)
    path = root / "charge_density" / f"{row['material_id']}.npy"
    density = np.asarray(np.load(path, mmap_mode="r"), dtype=np.float32)
    if density.ndim != 3 or min(density.shape) < 2 or not np.all(np.isfinite(density)):
        raise ValueError(f"Invalid density grid {path}: {density.shape}")
    return PeriodicDensity(
        cell=np.asarray(store.cells[structure_index], dtype=np.float64),
        frac_positions=np.asarray(store.frac_positions[atom_slice], dtype=np.float64),
        atomic_numbers=np.asarray(store.atomic_numbers[atom_slice], dtype=np.int32),
        density=density,
        identifier=str(row["material_id"]),
    )


def fibonacci_sphere(n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    index = np.arange(n, dtype=np.float64)
    z_value = 1.0 - 2.0 * (index + 0.5) / n
    theta = np.arccos(np.clip(z_value, -1.0, 1.0))
    phi = np.mod(index * np.pi * (3.0 - np.sqrt(5.0)), 2.0 * np.pi)
    sine = np.sin(theta)
    directions = np.column_stack((sine * np.cos(phi), sine * np.sin(phi), z_value))
    return directions, theta, phi


def real_spherical_harmonics(
    l_value: int, theta: np.ndarray, phi: np.ndarray
) -> np.ndarray:
    rows = []
    for m_value in range(-l_value, l_value + 1):
        if m_value < 0:
            harmonic = sph_harm_y(l_value, -m_value, theta, phi)
            rows.append(np.sqrt(2.0) * ((-1) ** m_value) * harmonic.imag)
        elif m_value == 0:
            rows.append(sph_harm_y(l_value, 0, theta, phi).real)
        else:
            harmonic = sph_harm_y(l_value, m_value, theta, phi)
            rows.append(np.sqrt(2.0) * ((-1) ** m_value) * harmonic.real)
    return np.asarray(rows, dtype=np.float64)


def inverse_sqrt(matrix: np.ndarray, eps: float = 1e-12) -> np.ndarray:
    values, vectors = np.linalg.eigh(matrix)
    if float(np.min(values)) <= eps:
        raise ValueError("Singular SOED radial basis")
    return (vectors * values[None, :] ** -0.5) @ vectors.T


def periodic_gradient(field: np.ndarray, cell: np.ndarray) -> np.ndarray:
    derivatives = [
        0.5
        * field.shape[axis]
        * (np.roll(field, -1, axis=axis) - np.roll(field, 1, axis=axis))
        for axis in range(3)
    ]
    fractional = np.stack(derivatives, axis=-1)
    cartesian = np.einsum(
        "...i,ij->...j", fractional, np.linalg.inv(cell).T, optimize=True
    )
    return np.linalg.norm(cartesian, axis=-1)


def periodic_laplacian(field: np.ndarray, cell: np.ndarray) -> np.ndarray:
    metric_inverse = np.linalg.inv(cell @ cell.T)
    result = np.zeros_like(field, dtype=np.float64)
    shape = field.shape
    for i_value in range(3):
        second = shape[i_value] ** 2 * (
            np.roll(field, -1, axis=i_value)
            - 2.0 * field
            + np.roll(field, 1, axis=i_value)
        )
        result += metric_inverse[i_value, i_value] * second
        for j_value in range(i_value + 1, 3):
            pp = np.roll(np.roll(field, -1, axis=i_value), -1, axis=j_value)
            pm = np.roll(np.roll(field, -1, axis=i_value), 1, axis=j_value)
            mp = np.roll(np.roll(field, 1, axis=i_value), -1, axis=j_value)
            mm = np.roll(np.roll(field, 1, axis=i_value), 1, axis=j_value)
            mixed = 0.25 * shape[i_value] * shape[j_value] * (pp - pm - mp + mm)
            result += 2.0 * metric_inverse[i_value, j_value] * mixed
    return result


def normalize_block(vector: np.ndarray) -> np.ndarray:
    vector = np.nan_to_num(np.asarray(vector, dtype=np.float64), copy=False)
    norm = float(np.linalg.norm(vector))
    return vector if norm <= 1e-14 else vector / norm


class PeriodicSOEDFamily:
    allowed_channels = {
        "density",
        "log_density",
        "density_contrast",
        "gradient",
        "laplacian",
        "electron_nuclear_contrast",
    }
    composition_feature_count = 5 * len(SOED_CHEMICAL_PROPERTIES) + 5
    global_feature_count = (
        5 * len(SOED_CHEMICAL_PROPERTIES)
        + 5
        + 10
        + 16
        + 4
        + 2 * (SOED_RECIPROCAL_BINS + 4)
    )
    compact_global_feature_count = global_feature_count - composition_feature_count
    electronic_radial_feature_count = (
        6 * SOED_N_RADIAL
        + len(SOED_CHEMICAL_PROPERTIES) * SOED_N_RADIAL
        + 4 * SOED_N_RADIAL
        + 16
    )

    def __init__(
        self,
        channel_sets: dict[str, Sequence[str]],
        active_atomic_numbers: Sequence[int] | None = None,
    ) -> None:
        self.channel_sets = {name: tuple(value) for name, value in channel_sets.items()}
        unknown = (
            set().union(*map(set, self.channel_sets.values())) - self.allowed_channels
        )
        if unknown:
            raise ValueError(f"Unknown SOED channels: {sorted(unknown)}")
        self.channels = tuple(
            dict.fromkeys(
                channel for values in self.channel_sets.values() for channel in values
            )
        )
        nodes, weights = np.polynomial.legendre.leggauss(SOED_N_RADIAL)
        self.radii = 0.5 * SOED_R_CUT * (nodes + 1.0)
        self.radial_weights = 0.5 * SOED_R_CUT * weights
        self.directions, theta, phi = fibonacci_sphere(SOED_N_ANGULAR)
        self.angular_weights = np.full(SOED_N_ANGULAR, 4.0 * np.pi / SOED_N_ANGULAR)
        self.harmonics = [
            real_spherical_harmonics(l_value, theta, phi)
            for l_value in range(SOED_L_MAX + 1)
        ]
        self.radial_basis = [
            self._build_radial_basis(l_value) for l_value in range(SOED_L_MAX + 1)
        ]
        self.displacements = self.radii[:, None, None] * self.directions[None, :, :]
        self.element_properties = self._build_element_property_table(
            active_atomic_numbers
        )
        self._hash_maps: dict[tuple[int, int, int], tuple[np.ndarray, np.ndarray]] = {}

    def _build_radial_basis(self, l_value: int) -> np.ndarray:
        centers = np.linspace(SOED_R_CUT / (SOED_N_MAX + 1), SOED_R_CUT, SOED_N_MAX)
        alphas = -np.log(1e-3) / centers**2
        raw = self.radii[None, :] ** l_value * np.exp(
            -alphas[:, None] * self.radii[None, :] ** 2
        )
        raw *= 0.5 * (np.cos(np.pi * self.radii / SOED_R_CUT) + 1.0)[None, :]
        measure = self.radial_weights * self.radii**2
        return inverse_sqrt((raw * measure[None, :]) @ raw.T) @ raw

    @staticmethod
    def _build_element_property_table(
        active_atomic_numbers: Sequence[int] | None = None,
    ) -> np.ndarray:
        from pymatgen.core import Element

        def safe_attribute(element: Any, name: str, default: Any = 0.0) -> Any:
            try:
                value = getattr(element, name)
            except (AttributeError, KeyError, TypeError, ValueError):
                return default
            return default if value is None else value

        def finite_value(value: Any, scale: float) -> float:
            try:
                number = float(value)
            except (AttributeError, KeyError, TypeError, ValueError):
                return 0.0
            return number / scale if np.isfinite(number) else 0.0

        def valence_channels(element: Any) -> tuple[float, float, float, float]:
            raw_configuration = safe_attribute(element, "full_electronic_structure", ())
            try:
                configuration = [
                    item for item in list(raw_configuration or ()) if len(item) >= 3
                ]
            except (TypeError, ValueError):
                configuration = []
            if not configuration:
                return 0.0, 0.0, 0.0, 0.0
            try:
                outer_n = max(int(item[0]) for item in configuration)
            except (TypeError, ValueError):
                return 0.0, 0.0, 0.0, 0.0
            values = {"s": 0.0, "p": 0.0, "d": 0.0, "f": 0.0}
            capacities = {"s": 2.0, "p": 6.0, "d": 10.0, "f": 14.0}
            for n_value, orbital, electrons in configuration:
                try:
                    n_integer = int(n_value)
                    orbital_name = str(orbital).lower()
                    electron_count = float(electrons)
                    include = (
                        (orbital_name in ("s", "p") and n_integer == outer_n)
                        or (orbital_name == "d" and n_integer >= outer_n - 1)
                        or (orbital_name == "f" and n_integer >= outer_n - 2)
                    )
                    if (
                        include
                        and orbital_name in values
                        and np.isfinite(electron_count)
                    ):
                        values[orbital_name] += (
                            electron_count / capacities[orbital_name]
                        )
                except (TypeError, ValueError):
                    continue
            return values["s"], values["p"], values["d"], values["f"]

        table = np.zeros((119, len(SOED_CHEMICAL_PROPERTIES)), dtype=np.float64)
        missing_x_elements = {2, 10, 18, 36, 54, 86, 118}
        requested_numbers = (
            range(1, 119)
            if active_atomic_numbers is None
            else sorted(set(map(int, active_atomic_numbers)))
        )
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            for atomic_number in requested_numbers:
                if atomic_number < 1 or atomic_number > 118:
                    raise ValueError(
                        f"Atomic number outside supported range: {atomic_number}"
                    )
                element = Element.from_Z(atomic_number)
                radius = safe_attribute(element, "atomic_radius_calculated", 0.0)
                if not finite_value(radius, 1.0):
                    radius = safe_attribute(element, "atomic_radius", 0.0)
                electronegativity = (
                    0.0
                    if atomic_number in missing_x_elements
                    else finite_value(safe_attribute(element, "X", 0.0), 4.0)
                )
                valence = valence_channels(element)
                table[atomic_number] = (
                    atomic_number / 100.0,
                    finite_value(safe_attribute(element, "row", 0.0), 7.0),
                    finite_value(safe_attribute(element, "group", 0.0), 18.0),
                    electronegativity,
                    finite_value(radius, 3.0),
                    finite_value(safe_attribute(element, "atomic_mass", 0.0), 250.0),
                    finite_value(
                        safe_attribute(element, "ionization_energy", 0.0), 25.0
                    ),
                    finite_value(
                        safe_attribute(element, "electron_affinity", 0.0), 4.0
                    ),
                    *valence,
                )
        return table

    def _normalized_density(self, sample: PeriodicDensity) -> np.ndarray:
        density = np.asarray(sample.density, dtype=np.float64)
        if SOED_DENSITY_NORMALIZATION == "none":
            return density
        scale = float(np.mean(np.abs(density)))
        if scale <= 1e-12:
            raise ValueError(f"All-zero density for {sample.identifier}")
        return density / scale

    @staticmethod
    def _deposit_weighted_atoms(
        shape: tuple[int, int, int],
        fractional_positions: np.ndarray,
        weights: np.ndarray,
    ) -> np.ndarray:
        grid = np.zeros(shape, dtype=np.float64)
        scaled = np.mod(fractional_positions, 1.0) * np.asarray(shape, dtype=np.float64)
        lower = np.floor(scaled).astype(int)
        fraction = scaled - lower
        for dx in (0, 1):
            for dy in (0, 1):
                for dz in (0, 1):
                    corner = np.asarray((dx, dy, dz), dtype=int)
                    indices = np.mod(
                        lower + corner[None, :], np.asarray(shape, dtype=int)
                    )
                    interpolation = np.prod(
                        np.where(corner[None, :] == 1, fraction, 1.0 - fraction), axis=1
                    )
                    np.add.at(
                        grid,
                        tuple(indices[:, axis] for axis in range(3)),
                        interpolation * weights,
                    )
        return grid

    def _electron_nuclear_contrast(
        self,
        sample: PeriodicDensity,
        density: np.ndarray,
    ) -> np.ndarray:
        shape = tuple(map(int, density.shape))
        nuclear_grid = self._deposit_weighted_atoms(
            shape,
            sample.frac_positions,
            sample.atomic_numbers.astype(np.float64),
        )
        frequencies = [np.fft.fftfreq(size) * size for size in shape]
        hkl = np.stack(np.meshgrid(*frequencies, indexing="ij"), axis=-1)
        reciprocal = (
            2.0
            * np.pi
            * np.einsum(
                "...i,ij->...j", hkl, np.linalg.inv(sample.cell).T, optimize=True
            )
        )
        reciprocal_squared = np.einsum(
            "...i,...i->...", reciprocal, reciprocal, optimize=True
        )
        kernel = np.exp(-reciprocal_squared / (4.0 * SOED_ELECTRON_NUCLEAR_ALPHA))
        reference = np.fft.ifftn(np.fft.fftn(nuclear_grid) * kernel).real
        reference /= max(float(np.mean(np.abs(reference))), 1e-12)
        contrast = density - reference
        return contrast / max(float(np.std(contrast)), 1e-12)

    def _channel_grids(self, sample: PeriodicDensity) -> dict[str, np.ndarray]:
        density = self._normalized_density(sample)
        grids: dict[str, np.ndarray] = {}
        if "density" in self.channels:
            grids["density"] = density
        if "log_density" in self.channels:
            grids["log_density"] = np.sign(density) * np.log1p(np.abs(density))
        if "density_contrast" in self.channels:
            sigma = tuple(
                max(0.75, SOED_CONTRAST_SMOOTHING_FRACTION * size)
                for size in density.shape
            )
            smooth = gaussian_filter(density, sigma=sigma, mode="wrap")
            contrast = density - smooth
            grids["density_contrast"] = contrast / max(float(np.std(contrast)), 1e-12)
        if "gradient" in self.channels:
            grids["gradient"] = SOED_R_CUT * periodic_gradient(density, sample.cell)
        if "laplacian" in self.channels:
            grids["laplacian"] = SOED_R_CUT**2 * periodic_laplacian(
                density, sample.cell
            )
        if "electron_nuclear_contrast" in self.channels:
            grids["electron_nuclear_contrast"] = self._electron_nuclear_contrast(
                sample, density
            )
        return grids

    def _sample_grid(
        self,
        grid: np.ndarray,
        center_fractional: np.ndarray,
        inverse_cell: np.ndarray,
    ) -> np.ndarray:
        fractional = (
            center_fractional[None, None, :] + self.displacements @ inverse_cell
        )
        shape = np.asarray(grid.shape, dtype=np.float64)
        coordinates = np.moveaxis(np.mod(fractional, 1.0) * shape[None, None, :], -1, 0)
        sampled = map_coordinates(
            grid,
            coordinates.reshape(3, -1),
            order=SOED_INTERPOLATION_ORDER,
            mode="grid-wrap",
            prefilter=SOED_INTERPOLATION_ORDER > 1,
        )
        return sampled.reshape(SOED_N_RADIAL, SOED_N_ANGULAR)

    def _power_spectrum(self, sampled: np.ndarray) -> np.ndarray:
        features: list[np.ndarray] = []
        radial_measure = self.radial_weights * self.radii**2
        for l_value in range(SOED_L_MAX + 1):
            coefficients = np.einsum(
                "cra,nr,ma,r,a->cnm",
                sampled,
                self.radial_basis[l_value],
                self.harmonics[l_value],
                radial_measure,
                self.angular_weights,
                optimize=True,
            )
            combined = coefficients.reshape(
                sampled.shape[0] * SOED_N_MAX, 2 * l_value + 1
            )
            power = combined @ combined.T
            features.append(power[np.triu_indices_from(power)])
        return np.concatenate(features)

    @staticmethod
    def _cell_angles(cell: np.ndarray) -> np.ndarray:
        lengths = np.linalg.norm(cell, axis=1)
        pairs = ((1, 2), (0, 2), (0, 1))
        return np.asarray(
            [
                np.degrees(
                    np.arccos(
                        np.clip(
                            np.dot(cell[i_value], cell[j_value])
                            / max(lengths[i_value] * lengths[j_value], 1e-14),
                            -1.0,
                            1.0,
                        )
                    )
                )
                for i_value, j_value in pairs
            ]
        )

    def _global_features(
        self, sample: PeriodicDensity, grids: dict[str, np.ndarray]
    ) -> np.ndarray:
        blocks = self._composition_features_raw(sample).tolist()
        lengths = np.linalg.norm(sample.cell, axis=1)
        angles = self._cell_angles(sample.cell)
        volume_per_atom = sample.volume / max(sample.n_atoms, 1)
        blocks.extend(
            (
                sample.n_atoms / 100.0,
                math.log1p(sample.n_atoms) / 5.0,
                volume_per_atom / 100.0,
                sample.n_atoms / max(sample.volume, 1e-12) * 10.0,
                *(lengths / 20.0),
                *(angles / 180.0),
            )
        )
        raw_density = np.asarray(sample.density, dtype=np.float64)
        density = grids["density"]
        raw_quantiles = np.quantile(raw_density, (0.05, 0.25, 0.50, 0.75, 0.95))
        signed_log = lambda value: math.copysign(
            math.log1p(abs(float(value))), float(value)
        )
        integral = float(np.mean(raw_density) * sample.volume)
        total_nuclear_charge = max(float(np.sum(sample.atomic_numbers)), 1.0)
        blocks.extend(
            (
                signed_log(np.mean(raw_density)),
                math.log1p(float(np.std(raw_density))),
                signed_log(np.min(raw_density)),
                signed_log(np.max(raw_density)),
                *(signed_log(value) for value in raw_quantiles),
                float(np.mean(density)),
                float(np.std(density)),
                float(np.min(density)),
                float(np.max(density)),
                signed_log(integral / total_nuclear_charge),
                signed_log(integral / max(sample.n_atoms, 1)),
                float(np.std(raw_density))
                / max(float(np.mean(np.abs(raw_density))), 1e-12),
            )
        )
        gradient = grids["gradient"]
        laplacian = grids["laplacian"]
        blocks.extend(
            (
                float(np.mean(gradient)),
                float(np.std(gradient)),
                float(np.mean(laplacian)),
                float(np.std(laplacian)),
            )
        )
        blocks.extend(self._reciprocal_features(sample, grids))
        vector = np.nan_to_num(np.asarray(blocks, dtype=np.float64))
        if len(vector) != self.global_feature_count:
            raise RuntimeError(
                f"Expected {self.global_feature_count} global SOED features, got {len(vector)}"
            )
        return vector

    def _composition_features_raw(self, sample: PeriodicDensity) -> np.ndarray:
        properties = self.element_properties[sample.atomic_numbers]
        blocks: list[float] = []
        for column in properties.T:
            blocks.extend(
                (
                    float(np.mean(column)),
                    float(np.std(column)),
                    float(np.min(column)),
                    float(np.max(column)),
                    float(np.ptp(column)),
                )
            )
        unique_numbers, counts = np.unique(sample.atomic_numbers, return_counts=True)
        fractions = counts / counts.sum()
        entropy = -float(np.sum(fractions * np.log(np.clip(fractions, 1e-14, None))))
        entropy /= max(math.log(len(fractions)), 1.0)
        blocks.extend(
            (
                len(fractions) / 10.0,
                entropy,
                float(np.max(fractions)),
                float(np.sum(fractions**2)),
                float(np.mean(sample.atomic_numbers)) / 100.0,
            )
        )
        vector = np.nan_to_num(np.asarray(blocks, dtype=np.float64))
        if len(vector) != self.composition_feature_count:
            raise RuntimeError(
                f"Expected {self.composition_feature_count} composition features, got {len(vector)}"
            )
        return vector

    def composition_features(self, sample: PeriodicDensity) -> np.ndarray:
        return normalize_block(self._composition_features_raw(sample)).astype(
            np.float32
        )

    def _reciprocal_features(
        self, sample: PeriodicDensity, grids: dict[str, np.ndarray]
    ) -> np.ndarray:
        grid_size = SOED_RECIPROCAL_GRID
        axes = [
            np.linspace(0.0, size, grid_size, endpoint=False)
            for size in grids["density"].shape
        ]
        coordinates = np.asarray(
            np.meshgrid(*axes, indexing="ij"), dtype=np.float64
        ).reshape(3, -1)
        frequencies = np.fft.fftfreq(grid_size) * grid_size
        hkl = np.stack(
            np.meshgrid(frequencies, frequencies, frequencies, indexing="ij"), axis=-1
        )
        reciprocal = (
            2.0
            * np.pi
            * np.einsum(
                "...i,ij->...j", hkl, np.linalg.inv(sample.cell).T, optimize=True
            )
        )
        magnitude = np.linalg.norm(reciprocal, axis=-1).reshape(-1)
        edges = np.linspace(
            0.0, max(float(np.max(magnitude)), 1e-12), SOED_RECIPROCAL_BINS + 1
        )
        bin_index = np.clip(
            np.digitize(magnitude, edges[1:-1], right=False),
            0,
            SOED_RECIPROCAL_BINS - 1,
        )
        features: list[float] = []
        for channel in ("density", "density_contrast"):
            sampled = map_coordinates(
                grids[channel], coordinates, order=1, mode="grid-wrap"
            ).reshape((grid_size,) * 3)
            sampled = sampled - float(np.mean(sampled))
            power = np.abs(np.fft.fftn(sampled)) ** 2 / sampled.size
            log_power = np.log1p(power.reshape(-1))
            features.extend(
                (
                    float(np.mean(log_power[bin_index == index]))
                    if np.any(bin_index == index)
                    else 0.0
                )
                for index in range(SOED_RECIPROCAL_BINS)
            )
            weights = power.reshape(-1)
            weights = weights / max(float(np.sum(weights)), 1e-12)
            spectral_entropy = -float(
                np.sum(weights * np.log(np.clip(weights, 1e-15, None)))
            ) / math.log(weights.size)
            features.extend(
                (
                    float(np.mean(log_power)),
                    float(np.std(log_power)),
                    float(np.max(log_power)),
                    spectral_entropy,
                )
            )
        return np.asarray(features, dtype=np.float64)

    def _hash_project(self, vector: np.ndarray, size: int, seed: int) -> np.ndarray:
        key = (len(vector), size, seed)
        if key not in self._hash_maps:
            indices = np.arange(len(vector), dtype=np.uint64)
            buckets = (
                (indices * np.uint64(2654435761) + np.uint64(seed * 2246822519))
                % np.uint64(size)
            ).astype(np.int64)
            signs = np.where(
                ((indices * np.uint64(3266489917) + np.uint64(seed)) & np.uint64(1))
                == 0,
                1.0,
                -1.0,
            )
            self._hash_maps[key] = buckets, signs
        buckets, signs = self._hash_maps[key]
        projected = np.bincount(
            buckets,
            weights=np.asarray(vector, dtype=np.float64) * signs,
            minlength=size,
        )
        return projected / math.sqrt(max(len(vector) / size, 1.0))

    def _atomic_block_masks(self, atomic_numbers: np.ndarray) -> tuple[np.ndarray, ...]:
        groups = np.rint(self.element_properties[atomic_numbers, 2] * 18.0).astype(int)
        f_block = ((atomic_numbers >= 57) & (atomic_numbers <= 71)) | (
            (atomic_numbers >= 89) & (atomic_numbers <= 103)
        )
        d_block = (groups >= 3) & (groups <= 12) & ~f_block
        s_block = ((groups <= 2) | (atomic_numbers == 2)) & ~f_block
        p_block = ~(s_block | d_block | f_block)
        return s_block, p_block, d_block, f_block

    def _local_chemical_blocks(
        self,
        sample: PeriodicDensity,
        array: np.ndarray,
        mean: np.ndarray,
    ) -> list[np.ndarray]:
        blocks: list[np.ndarray] = []
        properties = self.element_properties[sample.atomic_numbers]
        for property_index, property_values in enumerate(properties.T):
            centered = property_values - float(np.mean(property_values))
            scale = float(np.sum(np.abs(centered)))
            contrast = np.einsum("a,af->f", centered, array, optimize=True) / max(
                scale, 1e-12
            )
            projected = self._hash_project(
                contrast,
                SOED_CHEMICAL_PROJECTION_SIZE,
                101 + property_index,
            )
            blocks.append(normalize_block(projected))
        for block_index, mask in enumerate(
            self._atomic_block_masks(sample.atomic_numbers)
        ):
            block_mean = (
                np.mean(array[mask], axis=0) if np.any(mask) else np.zeros_like(mean)
            )
            projected = self._hash_project(
                block_mean - mean,
                SOED_BLOCK_PROJECTION_SIZE,
                401 + block_index,
            )
            blocks.append(normalize_block(projected))
        return blocks

    def _electronic_radial_features(
        self,
        sample: PeriodicDensity,
        radial_mean: np.ndarray,
        radial_std: np.ndarray,
    ) -> np.ndarray:
        radial_mean = np.asarray(radial_mean, dtype=np.float64)
        radial_std = np.asarray(radial_std, dtype=np.float64)
        radial_measure = 4.0 * np.pi * self.radial_weights * self.radii**2
        cumulative = np.cumsum(radial_mean * radial_measure[None, :], axis=1)
        blocks: list[np.ndarray] = []
        for array in (radial_mean, radial_std, cumulative):
            blocks.append(normalize_block(np.mean(array, axis=0)))
            blocks.append(normalize_block(np.std(array, axis=0)))
        properties = self.element_properties[sample.atomic_numbers]
        for property_values in properties.T:
            centered = property_values - float(np.mean(property_values))
            scale = float(np.sum(np.abs(centered)))
            contrast = np.einsum("a,ar->r", centered, radial_mean, optimize=True) / max(
                scale, 1e-12
            )
            blocks.append(normalize_block(contrast))
        overall = np.mean(radial_mean, axis=0)
        for mask in self._atomic_block_masks(sample.atomic_numbers):
            block_mean = (
                np.mean(radial_mean[mask], axis=0)
                if np.any(mask)
                else np.zeros_like(overall)
            )
            blocks.append(normalize_block(block_mean - overall))
        charge = cumulative[:, -1]
        charge_scaled = charge / max(float(np.mean(np.abs(charge))), 1e-12)
        nuclear_scaled = sample.atomic_numbers.astype(np.float64)
        nuclear_scaled /= max(float(np.mean(nuclear_scaled)), 1e-12)
        imbalance = charge_scaled - nuclear_scaled

        def seven_statistics(values: np.ndarray) -> list[float]:
            return [
                float(np.mean(values)),
                float(np.std(values)),
                float(np.min(values)),
                float(np.max(values)),
                *map(float, np.quantile(values, (0.25, 0.50, 0.75))),
            ]

        def safe_correlation(left: np.ndarray, right: np.ndarray) -> float:
            left_centered = left - float(np.mean(left))
            right_centered = right - float(np.mean(right))
            denominator = float(
                np.linalg.norm(left_centered) * np.linalg.norm(right_centered)
            )
            return (
                0.0
                if denominator <= 1e-12
                else float(np.dot(left_centered, right_centered) / denominator)
            )

        electronegativity = properties[
            :, SOED_CHEMICAL_PROPERTIES.index("electronegativity")
        ]
        summary = np.asarray(
            [
                *seven_statistics(charge_scaled),
                *seven_statistics(imbalance),
                safe_correlation(charge_scaled, nuclear_scaled),
                safe_correlation(charge_scaled, electronegativity),
            ],
            dtype=np.float64,
        )
        blocks.append(normalize_block(summary))
        vector = np.concatenate(blocks)
        if len(vector) != self.electronic_radial_feature_count:
            raise RuntimeError(
                f"Expected {self.electronic_radial_feature_count} radial electronic features, got {len(vector)}"
            )
        return vector

    def _compact_global_density_features(
        self,
        sample: PeriodicDensity,
        grids: dict[str, np.ndarray],
    ) -> np.ndarray:
        vector = self._global_features(sample, grids)[self.composition_feature_count :]
        if len(vector) != self.compact_global_feature_count:
            raise RuntimeError(
                f"Expected {self.compact_global_feature_count} compact global features, got {len(vector)}"
            )
        return vector

    def local_chemical_descriptor(
        self,
        sample: PeriodicDensity,
        local_array: np.ndarray,
    ) -> np.ndarray:
        array = np.asarray(local_array, dtype=np.float64)
        mean = np.mean(array, axis=0)
        blocks = [normalize_block(mean), normalize_block(np.std(array, axis=0))]
        blocks.extend(self._local_chemical_blocks(sample, array, mean))
        blocks.append(normalize_block(self._composition_features_raw(sample)))
        return normalize_block(np.concatenate(blocks)).astype(np.float32)

    def create_all(
        self, sample: PeriodicDensity, return_timings: bool = False
    ) -> dict[str, np.ndarray] | tuple[dict[str, np.ndarray], dict[str, float]]:
        grids = self._channel_grids(sample)
        inverse_cell = np.linalg.inv(sample.cell)
        local: dict[str, list[np.ndarray]] = {name: [] for name in self.channel_sets}
        descriptor_seconds = {name: 0.0 for name in self.channel_sets}
        radial_mean_rows: list[np.ndarray] = []
        radial_std_rows: list[np.ndarray] = []
        for center in sample.frac_positions:
            sampled_all: dict[str, np.ndarray] = {}
            sampled_seconds: dict[str, float] = {}
            for channel, grid in grids.items():
                started = time.perf_counter()
                sampled_all[channel] = self._sample_grid(grid, center, inverse_cell)
                sampled_seconds[channel] = time.perf_counter() - started
            radial_mean_rows.append(np.mean(sampled_all["density"], axis=1))
            radial_std_rows.append(np.std(sampled_all["density"], axis=1))
            spectrum_cache: dict[
                tuple[tuple[str, ...], bool], tuple[np.ndarray, float]
            ] = {}
            for descriptor_name, channels in self.channel_sets.items():
                started = time.perf_counter()
                separable = descriptor_name in SOED_SEPARABLE_DESCRIPTORS
                cache_key = (tuple(channels), separable)
                if cache_key not in spectrum_cache:
                    if separable:
                        spectrum = np.concatenate(
                            [
                                self._power_spectrum(sampled_all[channel][None, :, :])
                                for channel in channels
                            ]
                        )
                    else:
                        sampled = np.stack(
                            [sampled_all[channel] for channel in channels], axis=0
                        )
                        spectrum = self._power_spectrum(sampled)
                    spectrum_cache[cache_key] = (
                        spectrum,
                        time.perf_counter() - started,
                    )
                spectrum, spectrum_seconds = spectrum_cache[cache_key]
                local[descriptor_name].append(spectrum)
                descriptor_seconds[descriptor_name] += spectrum_seconds + sum(
                    sampled_seconds[channel] for channel in channels
                )
        radial_mean = np.asarray(radial_mean_rows, dtype=np.float64)
        radial_std = np.asarray(radial_std_rows, dtype=np.float64)
        started = time.perf_counter()
        radial_features = self._electronic_radial_features(
            sample, radial_mean, radial_std
        )
        radial_seconds = time.perf_counter() - started
        started = time.perf_counter()
        compact_global = self._compact_global_density_features(sample, grids)
        compact_global_seconds = time.perf_counter() - started
        for name in SOED_RADIAL_ELECTRONIC_DESCRIPTORS:
            descriptor_seconds[name] += radial_seconds
        for name in SOED_COMPACT_GLOBAL_DESCRIPTORS:
            descriptor_seconds[name] += compact_global_seconds
        output: dict[str, np.ndarray] = {}
        for descriptor_name, values in local.items():
            started = time.perf_counter()
            array = np.asarray(values, dtype=np.float64)
            mean = np.mean(array, axis=0)
            blocks = [normalize_block(mean), normalize_block(np.std(array, axis=0))]
            if descriptor_name in (
                *SOED_CHEMICAL_DESCRIPTORS,
                *SOED_LOCAL_CHEMICAL_DESCRIPTORS,
            ):
                blocks.extend(self._local_chemical_blocks(sample, array, mean))
            if descriptor_name in SOED_CHEMICAL_DESCRIPTORS:
                blocks.append(normalize_block(self._global_features(sample, grids)))
            elif descriptor_name in SOED_LOCAL_CHEMICAL_DESCRIPTORS:
                blocks.append(normalize_block(self._composition_features_raw(sample)))
                if descriptor_name in SOED_RADIAL_ELECTRONIC_DESCRIPTORS:
                    blocks.append(normalize_block(radial_features))
                if descriptor_name in SOED_COMPACT_GLOBAL_DESCRIPTORS:
                    blocks.append(normalize_block(compact_global))
            output[descriptor_name] = normalize_block(np.concatenate(blocks)).astype(
                np.float32
            )
            descriptor_seconds[descriptor_name] += time.perf_counter() - started
        if not return_timings:
            return output
        return output, descriptor_seconds

    def local_feature_count(self, descriptor_name: str) -> int:
        channels = len(self.channel_sets[descriptor_name])
        if descriptor_name in SOED_SEPARABLE_DESCRIPTORS:
            single = int((SOED_L_MAX + 1) * SOED_N_MAX * (SOED_N_MAX + 1) // 2)
            return channels * single
        k_value = channels * SOED_N_MAX
        return int((SOED_L_MAX + 1) * k_value * (k_value + 1) // 2)

    def feature_count(self, descriptor_name: str) -> int:
        local_count = self.local_feature_count(descriptor_name)
        if descriptor_name in SOED_CHEMICAL_DESCRIPTORS:
            return (
                2 * local_count
                + len(SOED_CHEMICAL_PROPERTIES) * SOED_CHEMICAL_PROJECTION_SIZE
                + 4 * SOED_BLOCK_PROJECTION_SIZE
                + self.global_feature_count
            )
        if descriptor_name in SOED_LOCAL_CHEMICAL_DESCRIPTORS:
            count = (
                2 * local_count
                + len(SOED_CHEMICAL_PROPERTIES) * SOED_CHEMICAL_PROJECTION_SIZE
                + 4 * SOED_BLOCK_PROJECTION_SIZE
                + self.composition_feature_count
            )
            if descriptor_name in SOED_RADIAL_ELECTRONIC_DESCRIPTORS:
                count += self.electronic_radial_feature_count
            if descriptor_name in SOED_COMPACT_GLOBAL_DESCRIPTORS:
                count += self.compact_global_feature_count
            return count
        return 2 * local_count


class PeriodicPaperSOAP:
    def __init__(self, overlap_basis: PeriodicSOEDFamily, alpha: float) -> None:
        self.overlap_basis = overlap_basis
        self.alpha = float(alpha)

    @staticmethod
    def _deposit_atoms(
        shape: tuple[int, int, int], fractional_positions: np.ndarray
    ) -> np.ndarray:
        grid = np.zeros(shape, dtype=np.float64)
        scaled = np.mod(fractional_positions, 1.0) * np.asarray(shape, dtype=np.float64)
        lower = np.floor(scaled).astype(int)
        fraction = scaled - lower
        for dx in (0, 1):
            for dy in (0, 1):
                for dz in (0, 1):
                    corner = np.asarray((dx, dy, dz), dtype=int)
                    indices = np.mod(
                        lower + corner[None, :], np.asarray(shape, dtype=int)
                    )
                    weights = np.prod(
                        np.where(corner[None, :] == 1, fraction, 1.0 - fraction), axis=1
                    )
                    np.add.at(
                        grid, tuple(indices[:, axis] for axis in range(3)), weights
                    )
        return grid

    def _gaussian_kernel(
        self, shape: tuple[int, int, int], cell: np.ndarray
    ) -> np.ndarray:
        frequencies = [np.fft.fftfreq(size) * size for size in shape]
        hkl = np.stack(np.meshgrid(*frequencies, indexing="ij"), axis=-1)
        reciprocal = (
            2.0
            * np.pi
            * np.einsum("...i,ij->...j", hkl, np.linalg.inv(cell).T, optimize=True)
        )
        reciprocal_squared = np.einsum(
            "...i,...i->...", reciprocal, reciprocal, optimize=True
        )
        return np.exp(-reciprocal_squared / (4.0 * self.alpha))

    def _periodic_scaffold_density(self, sample: PeriodicDensity) -> np.ndarray:
        shape = tuple(map(int, sample.density.shape))
        atomic_grid = self._deposit_atoms(shape, sample.frac_positions)
        kernel = self._gaussian_kernel(shape, sample.cell)
        scaffold = np.fft.ifftn(np.fft.fftn(atomic_grid) * kernel).real
        scale = float(np.mean(np.abs(scaffold)))
        if scale <= 1e-14:
            raise ValueError(f"All-zero scaffold density for {sample.identifier}")
        return scaffold / scale

    def create(
        self,
        sample: PeriodicDensity,
        return_local: bool = False,
    ) -> np.ndarray | tuple[np.ndarray, np.ndarray]:
        scaffold = self._periodic_scaffold_density(sample)
        inverse_cell = np.linalg.inv(sample.cell)
        local = []
        for center in sample.frac_positions:
            sampled = self.overlap_basis._sample_grid(scaffold, center, inverse_cell)[
                None, :, :
            ]
            local.append(self.overlap_basis._power_spectrum(sampled))
        array = np.asarray(local, dtype=np.float64)
        blocks = [normalize_block(np.mean(array, axis=0))]
        if SOAP_INCLUDE_LOCAL_STD:
            blocks.append(normalize_block(np.std(array, axis=0)))
        descriptor = normalize_block(np.concatenate(blocks)).astype(np.float32)
        if return_local:
            return descriptor, array.astype(np.float32)
        return descriptor

    def feature_count(self) -> int:
        local_count = int((SOED_L_MAX + 1) * SOED_N_MAX * (SOED_N_MAX + 1) // 2)
        return local_count * (2 if SOAP_INCLUDE_LOCAL_STD else 1)


def representation_names() -> list[str]:
    return [
        *SOAP_REPRESENTATIONS,
        *SOAP_COMPOSITION_REPRESENTATIONS,
        *SOAP_LOCAL_CHEMICAL_REPRESENTATIONS,
        *SOED_CHANNEL_SETS,
        COMPOSITION_ONLY_NAME,
        DENSITY_COMPOSITION_NAME,
    ]


def active_atomic_numbers(
    frame: pd.DataFrame, store: StructureStore
) -> tuple[int, ...]:
    numbers: set[int] = set()
    for structure_index in frame["structure_index"].to_numpy(dtype=np.int64):
        atom_slice = store.atom_slice(int(structure_index))
        numbers.update(map(int, store.atomic_numbers[atom_slice]))
    if not numbers:
        raise ValueError("No elements remain after data-quality filtering.")
    return tuple(sorted(numbers))


def descriptor_cache_payload(
    frame: pd.DataFrame,
    root: Path,
    store: StructureStore,
    atomic_numbers: tuple[int, ...],
    dimensions: dict[str, int],
) -> dict[str, Any]:
    material_ids = frame["material_id"].astype(str).tolist()
    density_digest = hashlib.sha256()
    for material_id in material_ids:
        path = root / "charge_density" / f"{material_id}.npy"
        stat = path.stat()
        density_digest.update(
            f"{material_id}|{stat.st_size}|{stat.st_mtime_ns}\n".encode("utf-8")
        )
    return {
        "schema_version": DESCRIPTOR_CACHE_SCHEMA_VERSION,
        "implementation_revision": "periodic-soap-soed-v12-r1",
        "structure_dataset_digest": store.digest,
        "n_samples": len(frame),
        "material_id_sha256": hashlib.sha256(
            "\n".join(material_ids).encode("utf-8")
        ).hexdigest(),
        "structure_index_sha256": hashlib.sha256(
            frame["structure_index"].to_numpy(dtype=np.int64).tobytes()
        ).hexdigest(),
        "charge_density_stat_sha256": density_digest.hexdigest(),
        "active_atomic_numbers": atomic_numbers,
        "dimensions": dimensions,
        "soap": {
            "alpha_grid": SOAP_ALPHA_GRID,
            "include_local_std": SOAP_INCLUDE_LOCAL_STD,
            "atom_deposition": SOAP_ATOM_DEPOSITION,
        },
        "soed": {
            "r_cut": SOED_R_CUT,
            "n_max": SOED_N_MAX,
            "l_max": SOED_L_MAX,
            "n_radial": SOED_N_RADIAL,
            "n_angular": SOED_N_ANGULAR,
            "interpolation_order": SOED_INTERPOLATION_ORDER,
            "density_normalization": SOED_DENSITY_NORMALIZATION,
            "contrast_smoothing_fraction": SOED_CONTRAST_SMOOTHING_FRACTION,
            "reciprocal_grid": SOED_RECIPROCAL_GRID,
            "reciprocal_bins": SOED_RECIPROCAL_BINS,
            "chemical_projection_size": SOED_CHEMICAL_PROJECTION_SIZE,
            "block_projection_size": SOED_BLOCK_PROJECTION_SIZE,
            "channel_sets": SOED_CHANNEL_SETS,
            "chemical_descriptors": SOED_CHEMICAL_DESCRIPTORS,
            "chemical_properties": SOED_CHEMICAL_PROPERTIES,
            "local_chemical_descriptors": SOED_LOCAL_CHEMICAL_DESCRIPTORS,
            "radial_electronic_descriptors": SOED_RADIAL_ELECTRONIC_DESCRIPTORS,
            "compact_global_descriptors": SOED_COMPACT_GLOBAL_DESCRIPTORS,
            "electronic_radial_feature_count": PeriodicSOEDFamily.electronic_radial_feature_count,
            "compact_global_feature_count": PeriodicSOEDFamily.compact_global_feature_count,
            "separable_descriptors": SOED_SEPARABLE_DESCRIPTORS,
            "electron_nuclear_alpha": SOED_ELECTRON_NUCLEAR_ALPHA,
        },
    }


def descriptor_cache_fingerprint(payload: dict[str, Any]) -> str:
    serialized = json.dumps(
        payload,
        sort_keys=True,
        separators=(",", ":"),
        default=json_default,
    ).encode("utf-8")
    return hashlib.sha256(serialized).hexdigest()[:24]


def atomic_save_npy(path: Path, values: np.ndarray) -> None:
    ensure_dir(path.parent)
    temporary = path.with_name(path.name + ".tmp")
    with temporary.open("wb") as handle:
        np.save(handle, values, allow_pickle=False)
    os.replace(temporary, path)


def atomic_save_frame(frame: pd.DataFrame, path: Path, separator: str = ",") -> None:
    ensure_dir(path.parent)
    temporary = path.with_name(path.name + ".tmp")
    frame.to_csv(temporary, sep=separator, index=False)
    os.replace(temporary, path)


def atomic_save_json(payload: dict[str, Any], path: Path) -> None:
    ensure_dir(path.parent)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, ensure_ascii=False, default=json_default),
        encoding="utf-8",
    )
    os.replace(temporary, path)


def load_descriptor_cache(
    cache_dir: Path,
    fingerprint: str,
    dimensions: dict[str, int],
    n_samples: int,
    logger: logging.Logger,
) -> tuple[dict[str, np.ndarray], np.ndarray, pd.DataFrame] | None:
    manifest_path = cache_dir / "manifest.json"
    valid_path = cache_dir / "valid.npy"
    summary_path = cache_dir / "feature_summary.csv"
    array_paths = {name: cache_dir / f"{name}.npy" for name in dimensions}
    required = (manifest_path, valid_path, summary_path, *array_paths.values())
    if not all(path.exists() for path in required):
        return None
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        expected_shapes = {
            name: [int(n_samples), int(dimension)]
            for name, dimension in dimensions.items()
        }
        if (
            manifest.get("status") != "complete"
            or manifest.get("fingerprint") != fingerprint
            or manifest.get("shapes") != expected_shapes
            or manifest.get("dtype") != "float32"
        ):
            logger.warning(
                "Descriptor cache manifest mismatch; features will be recomputed: %s",
                cache_dir,
            )
            return None
        mmap_mode = "r" if DESCRIPTOR_CACHE_MMAP else None
        arrays = {
            name: np.load(path, mmap_mode=mmap_mode, allow_pickle=False)
            for name, path in array_paths.items()
        }
        valid = np.load(valid_path, mmap_mode=mmap_mode, allow_pickle=False)
        if valid.shape != (n_samples,) or valid.dtype != np.bool_:
            raise ValueError(
                f"Invalid cached validity vector: shape={valid.shape}, dtype={valid.dtype}"
            )
        for name, array in arrays.items():
            if array.shape != tuple(expected_shapes[name]) or array.dtype != np.float32:
                raise ValueError(
                    f"Invalid cached array {name}: shape={array.shape}, dtype={array.dtype}"
                )
        summary = pd.read_csv(summary_path)
        if set(summary["representation"].astype(str)) != set(dimensions):
            raise ValueError("Cached descriptor summary does not match representations")
        logger.info(
            "Descriptor cache hit: fingerprint=%s samples=%d valid=%d mmap=%s path=%s",
            fingerprint,
            n_samples,
            int(np.asarray(valid).sum()),
            DESCRIPTOR_CACHE_MMAP,
            cache_dir,
        )
        return arrays, valid, summary
    except Exception as exc:
        logger.warning(
            "Descriptor cache could not be loaded (%r); features will be recomputed.",
            exc,
        )
        return None


def write_descriptor_cache(
    cache_dir: Path,
    fingerprint: str,
    identity: dict[str, Any],
    arrays: dict[str, np.ndarray],
    valid: np.ndarray,
    summary: pd.DataFrame,
    timing: pd.DataFrame,
    logger: logging.Logger,
) -> None:
    try:
        ensure_dir(cache_dir)
        atomic_save_json(
            {
                "status": "writing",
                "fingerprint": fingerprint,
                "started_at_unix": time.time(),
            },
            cache_dir / "manifest.json",
        )
        required_bytes = int(
            sum(array.nbytes for array in arrays.values()) + valid.nbytes
        )
        free_bytes = int(shutil.disk_usage(cache_dir).free)
        if free_bytes < int(1.05 * required_bytes):
            raise OSError(
                f"Insufficient free space for descriptor cache: required about "
                f"{required_bytes / 2**30:.2f} GiB, free {free_bytes / 2**30:.2f} GiB"
            )
        logger.info(
            "Writing descriptor cache: approximately %.2f GiB; free space %.2f GiB",
            required_bytes / 2**30,
            free_bytes / 2**30,
        )
        for name, array in arrays.items():
            logger.info("Caching descriptor %s with shape=%s", name, array.shape)
            atomic_save_npy(
                cache_dir / f"{name}.npy", np.asarray(array, dtype=np.float32)
            )
        atomic_save_npy(cache_dir / "valid.npy", np.asarray(valid, dtype=np.bool_))
        atomic_save_frame(summary, cache_dir / "feature_summary.csv")
        atomic_save_frame(timing, cache_dir / "feature_timing.dat", separator="\t")
        manifest = {
            "status": "complete",
            "fingerprint": fingerprint,
            "created_at_unix": time.time(),
            "dtype": "float32",
            "shapes": {
                name: list(map(int, array.shape)) for name, array in arrays.items()
            },
            "valid_samples": int(valid.sum()),
            "identity": identity,
        }
        atomic_save_json(manifest, cache_dir / "manifest.json")
        logger.info(
            "Descriptor cache created: fingerprint=%s valid=%d/%d path=%s",
            fingerprint,
            int(valid.sum()),
            len(valid),
            cache_dir,
        )
    except Exception as exc:
        logger.warning(
            "Descriptor cache write failed (%r). Current run will continue with in-memory features.",
            exc,
        )


def generate_features(
    frame: pd.DataFrame,
    root: Path,
    output: Path,
    store: StructureStore,
    logger: logging.Logger,
) -> tuple[dict[str, np.ndarray], np.ndarray, pd.DataFrame]:
    atomic_numbers = active_atomic_numbers(frame, store)
    logger.info(
        "Initializing descriptor chemistry for %d active elements: %s",
        len(atomic_numbers),
        ",".join(map(str, atomic_numbers)),
    )
    soed = PeriodicSOEDFamily(SOED_CHANNEL_SETS, active_atomic_numbers=atomic_numbers)
    soap_models = {
        name: PeriodicPaperSOAP(soed, alpha)
        for name, alpha in SOAP_ALPHA_BY_NAME.items()
    }
    dimensions = {name: soap.feature_count() for name, soap in soap_models.items()}
    local_chemical_count = (
        next(iter(soap_models.values())).feature_count()
        + len(SOED_CHEMICAL_PROPERTIES) * SOED_CHEMICAL_PROJECTION_SIZE
        + 4 * SOED_BLOCK_PROJECTION_SIZE
        + soed.composition_feature_count
    )
    dimensions.update(
        {name: local_chemical_count for name in SOAP_LOCAL_CHEMICAL_REPRESENTATIONS}
    )
    dimensions.update({name: soed.feature_count(name) for name in SOED_CHANNEL_SETS})
    dimensions[COMPOSITION_ONLY_NAME] = soed.composition_feature_count
    diagnostics_dir = ensure_dir(output / "descriptor_diagnostics")
    logger.info(
        "Preparing descriptor-cache identity from samples, density files, and parameters."
    )
    cache_identity = descriptor_cache_payload(
        frame, root, store, atomic_numbers, dimensions
    )
    cache_fingerprint = descriptor_cache_fingerprint(cache_identity)
    descriptor_cache_dir = (
        store.cache_dir.parent / "features" / f"descriptor_{cache_fingerprint}"
    )
    if REUSE_DESCRIPTOR_CACHE and not FORCE_RECOMPUTE_FEATURES:
        cached = load_descriptor_cache(
            descriptor_cache_dir,
            cache_fingerprint,
            dimensions,
            len(frame),
            logger,
        )
        if cached is not None:
            cached_arrays, cached_valid, cached_summary = cached
            cached_summary.to_csv(diagnostics_dir / "feature_summary.csv", index=False)
            save_dat(cached_summary, diagnostics_dir / "feature_summary.dat")
            cached_timing = descriptor_cache_dir / "feature_timing.dat"
            if cached_timing.exists():
                shutil.copy2(cached_timing, diagnostics_dir / "feature_timing.dat")
            atomic_save_json(
                {
                    "status": "hit",
                    "fingerprint": cache_fingerprint,
                    "cache_dir": str(descriptor_cache_dir),
                    "mmap": DESCRIPTOR_CACHE_MMAP,
                },
                diagnostics_dir / "descriptor_cache_status.json",
            )
            (diagnostics_dir / "cache_policy.txt").write_text(
                "Parsed structures and SOAP/SOED features are cached. Descriptor cache identity includes sample order, density-file metadata, active elements, numerical parameters, and cache schema.\n",
                encoding="utf-8",
            )
            return cached_arrays, cached_valid, cached_summary
    elif FORCE_RECOMPUTE_FEATURES:
        logger.info("Descriptor cache bypassed because force recomputation is enabled.")
    arrays = {
        name: np.full((len(frame), dimension), np.nan, dtype=np.float32)
        for name, dimension in dimensions.items()
    }
    completed = np.zeros(len(frame), dtype=bool)
    timing_rows: list[dict[str, Any]] = []
    error_rows: list[dict[str, Any]] = []

    def compute_one(index: int):
        timings: dict[str, float] = {}
        try:
            sample = row_to_sample(index, frame.iloc[index], root, store)
            soap_values = {}
            for name, soap in soap_models.items():
                started = time.perf_counter()
                base_descriptor, local_array = soap.create(sample, return_local=True)
                base_elapsed = time.perf_counter() - started
                soap_values[name] = base_descriptor
                timings[name] = base_elapsed
                chemical_name = SOAP_LOCAL_CHEMICAL_BY_SOURCE[name]
                started = time.perf_counter()
                soap_values[chemical_name] = soed.local_chemical_descriptor(
                    sample, local_array
                )
                timings[chemical_name] = base_elapsed + time.perf_counter() - started
            soed_values, soed_timings = soed.create_all(sample, return_timings=True)
            timings.update(soed_timings)
            started = time.perf_counter()
            composition = soed.composition_features(sample)
            timings[COMPOSITION_ONLY_NAME] = time.perf_counter() - started
            return (
                index,
                {
                    **soap_values,
                    **soed_values,
                    COMPOSITION_ONLY_NAME: composition,
                },
                timings,
                None,
            )
        except Exception:
            return index, None, timings, traceback.format_exc()

    logger.info(
        "Descriptor generation starts: n=%d workers=%d dimensions=%s cache_fingerprint=%s",
        len(frame),
        DESCRIPTOR_WORKERS,
        dimensions,
        cache_fingerprint,
    )
    for offset in range(0, len(frame), DESCRIPTOR_CHUNK_SIZE):
        batch = np.arange(
            offset, min(offset + DESCRIPTOR_CHUNK_SIZE, len(frame)), dtype=int
        )
        with ThreadPoolExecutor(max_workers=DESCRIPTOR_WORKERS) as executor:
            futures = [executor.submit(compute_one, int(index)) for index in batch]
            for future in as_completed(futures):
                index, values, elapsed_by_representation, error = future.result()
                material_id = str(frame.iloc[index]["material_id"])
                timing_rows.extend(
                    {
                        "index": index,
                        "material_id": material_id,
                        "representation": name,
                        "seconds": elapsed,
                    }
                    for name, elapsed in elapsed_by_representation.items()
                )
                if error is not None or values is None:
                    error_rows.append(
                        {"index": index, "material_id": material_id, "error": error}
                    )
                    logger.error("Descriptor generation failed for %s", material_id)
                    continue
                for name, vector in values.items():
                    if vector.shape != (dimensions[name],):
                        raise RuntimeError(
                            f"Unexpected {name} shape for {material_id}: {vector.shape}"
                        )
                    arrays[name][index] = vector
                completed[index] = True
        logger.info(
            "Descriptor progress: %d/%d (%.1f%%)",
            int(completed.sum()),
            len(frame),
            100.0 * float(completed.mean()),
        )

    valid = completed & np.logical_and.reduce(
        [np.all(np.isfinite(array), axis=1) for array in arrays.values()]
    )
    timing_frame = pd.DataFrame(timing_rows)
    save_dat(timing_frame, diagnostics_dir / "feature_timing.dat")
    if error_rows:
        save_dat(pd.DataFrame(error_rows), diagnostics_dir / "feature_errors.dat")
    summary_rows = []
    for name, dimension in dimensions.items():
        values = (
            timing_frame.loc[timing_frame["representation"] == name, "seconds"]
            if len(timing_frame)
            else pd.Series(dtype=float)
        )
        if name in SOAP_REPRESENTATIONS:
            family = "Paper-inspired periodic SOAP"
        elif name in SOAP_LOCAL_CHEMICAL_REPRESENTATIONS:
            family = "Local-chemistry-matched periodic SOAP"
        elif name == COMPOSITION_ONLY_NAME:
            family = "Composition-only control"
        elif name in SOED_ENHANCED_CANDIDATES:
            family = "Electronic SOED candidate"
        else:
            family = "SOED ablation"
        summary_rows.append(
            {
                "representation": name,
                "family": family,
                "feature_dimension": dimension,
                "valid_samples": int(valid.sum()),
                "seconds_per_structure": (
                    float(values.mean()) if len(values) else np.nan
                ),
            }
        )
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(diagnostics_dir / "feature_summary.csv", index=False)
    save_dat(summary, diagnostics_dir / "feature_summary.dat")
    if REUSE_DESCRIPTOR_CACHE:
        write_descriptor_cache(
            descriptor_cache_dir,
            cache_fingerprint,
            cache_identity,
            arrays,
            valid,
            summary,
            timing_frame,
            logger,
        )
    atomic_save_json(
        {
            "status": "recomputed",
            "fingerprint": cache_fingerprint,
            "cache_dir": str(descriptor_cache_dir),
            "cache_enabled": REUSE_DESCRIPTOR_CACHE,
            "force_recompute": FORCE_RECOMPUTE_FEATURES,
        },
        diagnostics_dir / "descriptor_cache_status.json",
    )
    (diagnostics_dir / "cache_policy.txt").write_text(
        "Parsed structures and SOAP/SOED features are cached. Descriptor cache identity includes sample order, density-file metadata, active elements, numerical parameters, and cache schema.\n",
        encoding="utf-8",
    )
    logger.info("Feature dimensions: %s", dimensions)
    return arrays, valid, summary


def regression_metrics(y_true: np.ndarray, y_pred: np.ndarray) -> dict[str, float]:
    y_true = np.asarray(y_true, dtype=float)
    y_pred = np.asarray(y_pred, dtype=float)
    result = {
        "mae": float(mean_absolute_error(y_true, y_pred)),
        "rmse": float(np.sqrt(mean_squared_error(y_true, y_pred))),
        "r2": float(r2_score(y_true, y_pred)) if len(y_true) >= 2 else math.nan,
    }
    mask = y_true >= TAIL_THRESHOLD_EV
    result["n_tail"] = int(mask.sum())
    if mask.any():
        residual = y_pred[mask] - y_true[mask]
        result.update(
            {
                "tail_mae": float(np.mean(np.abs(residual))),
                "tail_rmse": float(np.sqrt(np.mean(residual**2))),
                "tail_bias": float(np.mean(residual)),
            }
        )
    else:
        result.update({"tail_mae": np.nan, "tail_rmse": np.nan, "tail_bias": np.nan})
    return result


def regression_sample_weights(
    y_value: np.ndarray,
    mode: str | bool = "tail_weighted",
) -> np.ndarray:
    y_value = np.asarray(y_value, dtype=float)
    if isinstance(mode, bool):
        mode = "tail_weighted" if mode else "unweighted"
    if mode == "unweighted":
        return np.ones(len(y_value), dtype=np.float32)
    if mode == "midgap_weighted":
        distance = (y_value - MIDGAP_CENTER_EV) / max(MIDGAP_WIDTH_EV, 1e-12)
        return (1.0 + MIDGAP_WEIGHT_MULTIPLIER * np.exp(-0.5 * distance**2)).astype(
            np.float32
        )
    if mode != "tail_weighted":
        raise ValueError(f"Unknown regression weighting mode: {mode}")
    if not USE_TAIL_SAMPLE_WEIGHTS:
        return np.ones(len(y_value), dtype=np.float32)
    denominator = max(TAIL_FULL_WEIGHT_EV - TAIL_THRESHOLD_EV, 1e-12)
    fraction = np.clip((y_value - TAIL_THRESHOLD_EV) / denominator, 0.0, 1.0)
    return (1.0 + (TAIL_WEIGHT_MULTIPLIER - 1.0) * fraction).astype(np.float32)


def validation_blend_objective(
    y_true: np.ndarray,
    prediction_matrix: np.ndarray,
    weights: np.ndarray,
    l2_strength: float,
    protect_midgap: bool,
) -> float:
    prediction = prediction_matrix @ weights
    residual = prediction - y_true
    value = float(np.mean(residual**2))
    if protect_midgap:
        midgap = (y_true >= MIDGAP_LOWER_EV) & (y_true < MIDGAP_UPPER_EV)
        if midgap.any():
            mid_residual = residual[midgap]
            value += MIDGAP_OBJECTIVE_WEIGHT * float(np.mean(mid_residual**2))
            excess = np.maximum(
                np.abs(mid_residual) - CATASTROPHIC_ERROR_THRESHOLD_EV,
                0.0,
            )
            value += CATASTROPHIC_OBJECTIVE_WEIGHT * float(np.mean(excess**2))
    if l2_strength > 0.0:
        equal = np.full(len(weights), 1.0 / len(weights), dtype=float)
        value += l2_strength * float(np.sum((weights - equal) ** 2))
    return value


def validation_simplex_weights(
    y_true: np.ndarray,
    predictions: Sequence[np.ndarray],
    l2_strength: float,
    protect_midgap: bool,
) -> np.ndarray:
    matrix = np.column_stack([np.asarray(value, dtype=float) for value in predictions])
    count = matrix.shape[1]
    initial = np.full(count, 1.0 / count, dtype=float)
    result = minimize(
        lambda value: validation_blend_objective(
            np.asarray(y_true, dtype=float),
            matrix,
            np.asarray(value, dtype=float),
            l2_strength,
            protect_midgap,
        ),
        initial,
        method="SLSQP",
        bounds=[(0.0, 1.0)] * count,
        constraints=[{"type": "eq", "fun": lambda value: float(np.sum(value) - 1.0)}],
        options={"ftol": 1e-12, "maxiter": 1000, "disp": False},
    )
    if not result.success or not np.all(np.isfinite(result.x)):
        return initial
    weights = np.clip(np.asarray(result.x, dtype=float), 0.0, 1.0)
    total = float(weights.sum())
    return weights / total if total > 0.0 else initial


def classification_metrics(
    y_true: np.ndarray, probability: np.ndarray, threshold: float
) -> dict[str, float]:
    y_true = np.asarray(y_true, dtype=np.int8)
    probability = np.clip(np.asarray(probability, dtype=float), 1e-7, 1.0 - 1e-7)
    prediction = (probability >= threshold).astype(np.int8)
    bins = np.clip(np.digitize(probability, np.linspace(0.1, 0.9, 9)), 0, 9)
    ece = 0.0
    for bin_index in range(10):
        mask = bins == bin_index
        if mask.any():
            ece += float(np.mean(mask)) * abs(
                float(np.mean(probability[mask])) - float(np.mean(y_true[mask]))
            )
    return {
        "auroc": float(roc_auc_score(y_true, probability)),
        "auprc": float(average_precision_score(y_true, probability)),
        "f1": float(f1_score(y_true, prediction, zero_division=0)),
        "balanced_accuracy": float(balanced_accuracy_score(y_true, prediction)),
        "mcc": float(matthews_corrcoef(y_true, prediction)),
        "brier": float(brier_score_loss(y_true, probability)),
        "logloss": float(log_loss(y_true, probability, labels=[0, 1])),
        "ece_10bin": float(ece),
    }


def select_classifier_threshold(y_true: np.ndarray, probability: np.ndarray) -> float:
    best_threshold = 0.5
    best_score = -math.inf
    for threshold in CLASSIFIER_THRESHOLD_GRID:
        score = classification_metrics(y_true, probability, threshold)[
            CLASSIFIER_THRESHOLD_METRIC
        ]
        if score > best_score:
            best_score = score
            best_threshold = float(threshold)
    return best_threshold


def subsample_indices(count: int, maximum: int, seed: int) -> np.ndarray:
    if maximum <= 0 or count <= maximum:
        return np.arange(count)
    return np.sort(
        np.random.default_rng(seed).choice(count, size=maximum, replace=False)
    )


def fit_xgb_regressor(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_valid: np.ndarray,
    y_valid: np.ndarray,
    params: dict[str, Any],
    use_gpu: bool,
    train_weights: np.ndarray,
    valid_weights: np.ndarray,
    logger: logging.Logger,
    log_structure: bool = True,
):
    from xgboost import XGBRegressor

    model = XGBRegressor(
        objective="reg:squarederror",
        eval_metric="rmse",
        random_state=RANDOM_SEED,
        n_jobs=max(1, os.cpu_count() or 1),
        tree_method="hist",
        device=f"cuda:{GPU_ID}" if use_gpu else "cpu",
        early_stopping_rounds=EARLY_STOPPING_ROUNDS,
        **params,
    )
    model.fit(
        x_train,
        y_train,
        sample_weight=train_weights,
        eval_set=[(x_train, y_train), (x_valid, y_valid)],
        sample_weight_eval_set=[train_weights, valid_weights],
        verbose=False,
    )
    evaluation = model.evals_result()
    history = pd.DataFrame(
        {
            "iteration": np.arange(1, len(evaluation["validation_0"]["rmse"]) + 1),
            "train_rmse": evaluation["validation_0"]["rmse"],
            "valid_rmse": evaluation["validation_1"]["rmse"],
        }
    )
    structure = json.dumps(model.get_params(), indent=2, default=json_default)
    if log_structure:
        logger.info("Model structure xgboost:\n%s", structure)
    return model, history, structure


def fit_xgb_classifier(
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_valid: np.ndarray,
    y_valid: np.ndarray,
    params: dict[str, Any],
    use_gpu: bool,
    logger: logging.Logger,
    log_structure: bool = True,
):
    from xgboost import XGBClassifier

    positive = max(1, int(np.sum(y_train == 1)))
    negative = max(1, int(np.sum(y_train == 0)))
    model = XGBClassifier(
        objective="binary:logistic",
        eval_metric="logloss",
        random_state=RANDOM_SEED,
        n_jobs=max(1, os.cpu_count() or 1),
        tree_method="hist",
        device=f"cuda:{GPU_ID}" if use_gpu else "cpu",
        early_stopping_rounds=EARLY_STOPPING_ROUNDS,
        scale_pos_weight=negative / positive if USE_BALANCED_CLASS_WEIGHTS else 1.0,
        **params,
    )
    model.fit(
        x_train,
        y_train,
        eval_set=[(x_train, y_train), (x_valid, y_valid)],
        verbose=False,
    )
    evaluation = model.evals_result()
    history = pd.DataFrame(
        {
            "iteration": np.arange(1, len(evaluation["validation_0"]["logloss"]) + 1),
            "train_logloss": evaluation["validation_0"]["logloss"],
            "valid_logloss": evaluation["validation_1"]["logloss"],
        }
    )
    structure = json.dumps(model.get_params(), indent=2, default=json_default)
    if log_structure:
        logger.info("Classifier structure xgboost:\n%s", structure)
    return model, history, structure


def xgb_predict(model: Any, values: np.ndarray) -> np.ndarray:
    import xgboost as xgb

    matrix = xgb.DMatrix(np.asarray(values, dtype=np.float32))
    best_iteration = getattr(model, "best_iteration", None)
    iteration_range = (
        (0, int(best_iteration) + 1) if best_iteration is not None else (0, 0)
    )
    return np.asarray(
        model.get_booster().predict(matrix, iteration_range=iteration_range),
        dtype=float,
    )


def tune_regressor(
    representation: str,
    task_name: str,
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_valid: np.ndarray,
    y_valid: np.ndarray,
    hardware: dict[str, Any],
    output: Path,
    logger: logging.Logger,
    weighting_mode: str,
) -> dict[str, Any]:
    params = dict(XGBOOST_PARAMS)
    if not ENABLE_OPTUNA or representation not in OPTUNA_REPRESENTATIONS:
        return params
    import optuna

    train_index = subsample_indices(len(y_train), OPTUNA_MAX_TRAIN_SAMPLES, RANDOM_SEED)
    valid_index = subsample_indices(
        len(y_valid), max(2000, OPTUNA_MAX_TRAIN_SAMPLES // 4), RANDOM_SEED + 1
    )
    train_weights = regression_sample_weights(y_train[train_index], weighting_mode)
    valid_weights = regression_sample_weights(y_valid[valid_index], weighting_mode)

    def objective(trial: optuna.Trial) -> float:
        trial_params = {
            "n_estimators": 4000,
            "max_depth": trial.suggest_int("max_depth", 3, 8),
            "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.08, log=True),
            "min_child_weight": trial.suggest_float(
                "min_child_weight", 3.0, 30.0, log=True
            ),
            "subsample": trial.suggest_float("subsample", 0.60, 1.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.50, 1.0),
            "reg_alpha": trial.suggest_float("reg_alpha", 1e-3, 2.0, log=True),
            "reg_lambda": trial.suggest_float("reg_lambda", 1.0, 30.0, log=True),
            "max_bin": trial.suggest_categorical("max_bin", [128, 256, 512]),
        }
        try:
            model, _, _ = fit_xgb_regressor(
                x_train[train_index],
                y_train[train_index],
                x_valid[valid_index],
                y_valid[valid_index],
                trial_params,
                bool(hardware.get("gpu_available")),
                train_weights,
                valid_weights,
                logger,
                False,
            )
        except Exception:
            model, _, _ = fit_xgb_regressor(
                x_train[train_index],
                y_train[train_index],
                x_valid[valid_index],
                y_valid[valid_index],
                trial_params,
                False,
                train_weights,
                valid_weights,
                logger,
                False,
            )
        prediction = np.maximum(xgb_predict(model, x_valid[valid_index]), 0.0)
        return float(np.sqrt(mean_squared_error(y_valid[valid_index], prediction)))

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(
        direction="minimize", sampler=optuna.samplers.TPESampler(seed=RANDOM_SEED)
    )
    study.optimize(objective, n_trials=OPTUNA_TRIALS, timeout=OPTUNA_TIMEOUT_SECONDS)
    params.update(study.best_params)
    optuna_dir = ensure_dir(output / "optuna")
    study.trials_dataframe().to_csv(
        optuna_dir / f"{task_name}_{representation}.csv", index=False
    )
    (optuna_dir / f"{task_name}_{representation}_best.json").write_text(
        json.dumps({"value": study.best_value, "params": params}, indent=2),
        encoding="utf-8",
    )
    logger.info(
        "Optuna %s %s best RMSE=%.6g params=%s",
        task_name,
        representation,
        study.best_value,
        study.best_params,
    )
    return params


def tune_classifier(
    representation: str,
    x_train: np.ndarray,
    y_train: np.ndarray,
    x_valid: np.ndarray,
    y_valid: np.ndarray,
    hardware: dict[str, Any],
    output: Path,
    logger: logging.Logger,
) -> dict[str, Any]:
    params = dict(XGBOOST_PARAMS)
    if not ENABLE_OPTUNA or representation not in OPTUNA_REPRESENTATIONS:
        return params
    import optuna

    train_index = subsample_indices(
        len(y_train), OPTUNA_MAX_TRAIN_SAMPLES, RANDOM_SEED + 2
    )
    valid_index = subsample_indices(
        len(y_valid), max(2000, OPTUNA_MAX_TRAIN_SAMPLES // 4), RANDOM_SEED + 3
    )

    def objective(trial: optuna.Trial) -> float:
        trial_params = {
            "n_estimators": 4000,
            "max_depth": trial.suggest_int("max_depth", 3, 8),
            "learning_rate": trial.suggest_float("learning_rate", 0.01, 0.08, log=True),
            "min_child_weight": trial.suggest_float(
                "min_child_weight", 3.0, 30.0, log=True
            ),
            "subsample": trial.suggest_float("subsample", 0.60, 1.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.50, 1.0),
            "reg_alpha": trial.suggest_float("reg_alpha", 1e-3, 2.0, log=True),
            "reg_lambda": trial.suggest_float("reg_lambda", 1.0, 30.0, log=True),
            "max_bin": trial.suggest_categorical("max_bin", [128, 256, 512]),
        }
        try:
            model, _, _ = fit_xgb_classifier(
                x_train[train_index],
                y_train[train_index],
                x_valid[valid_index],
                y_valid[valid_index],
                trial_params,
                bool(hardware.get("gpu_available")),
                logger,
                False,
            )
        except Exception:
            model, _, _ = fit_xgb_classifier(
                x_train[train_index],
                y_train[train_index],
                x_valid[valid_index],
                y_valid[valid_index],
                trial_params,
                False,
                logger,
                False,
            )
        probability = xgb_predict(model, x_valid[valid_index])
        return float(log_loss(y_valid[valid_index], probability, labels=[0, 1]))

    optuna.logging.set_verbosity(optuna.logging.WARNING)
    study = optuna.create_study(
        direction="minimize", sampler=optuna.samplers.TPESampler(seed=RANDOM_SEED + 1)
    )
    study.optimize(objective, n_trials=OPTUNA_TRIALS, timeout=OPTUNA_TIMEOUT_SECONDS)
    params.update(study.best_params)
    optuna_dir = ensure_dir(output / "optuna")
    study.trials_dataframe().to_csv(
        optuna_dir / f"classifier_{representation}.csv", index=False
    )
    return params


def fit_one_regression(
    representation: str,
    target: str,
    task_name: str,
    x_sets: dict[str, np.ndarray],
    y_sets: dict[str, np.ndarray],
    ids_sets: dict[str, np.ndarray],
    hardware: dict[str, Any],
    output: Path,
    logger: logging.Logger,
    weighting_mode: str = "tail_weighted",
) -> dict[str, Any]:
    run_dir = ensure_dir(
        output / "models" / target / task_name / representation / MODEL_NAME
    )
    sample_weights = {
        split: regression_sample_weights(values, weighting_mode)
        for split, values in y_sets.items()
    }
    params = tune_regressor(
        representation,
        task_name,
        x_sets["train"],
        y_sets["train"],
        x_sets["validation"],
        y_sets["validation"],
        hardware,
        output,
        logger,
        weighting_mode,
    )
    started = time.perf_counter()
    used_gpu = bool(hardware.get("gpu_available"))
    try:
        model, history, structure = fit_xgb_regressor(
            x_sets["train"],
            y_sets["train"],
            x_sets["validation"],
            y_sets["validation"],
            params,
            used_gpu,
            sample_weights["train"],
            sample_weights["validation"],
            logger,
        )
    except Exception as exc:
        if not used_gpu:
            raise
        logger.warning("XGBoost GPU regression failed (%r); retrying on CPU", exc)
        used_gpu = False
        model, history, structure = fit_xgb_regressor(
            x_sets["train"],
            y_sets["train"],
            x_sets["validation"],
            y_sets["validation"],
            params,
            False,
            sample_weights["train"],
            sample_weights["validation"],
            logger,
        )
    predictions = {
        split: xgb_predict(model, values) for split, values in x_sets.items()
    }
    if CLIP_NEGATIVE_GAP_PREDICTIONS:
        predictions = {
            split: np.maximum(value, 0.0) for split, value in predictions.items()
        }
    elapsed = time.perf_counter() - started
    joblib.dump(model, run_dir / "model.joblib")
    history.to_csv(run_dir / "rmse_history.csv", index=False)
    save_dat(history, run_dir / "rmse_history.dat")
    (run_dir / "model_structure.txt").write_text(structure, encoding="utf-8")
    (run_dir / "parameters.json").write_text(
        json.dumps(
            {
                **params,
                "weighting_mode": weighting_mode,
                "tail_sample_weights": bool(
                    USE_TAIL_SAMPLE_WEIGHTS and weighting_mode == "tail_weighted"
                ),
                "midgap_sample_weights": bool(weighting_mode == "midgap_weighted"),
                "tail_threshold_ev": TAIL_THRESHOLD_EV,
                "tail_weight_multiplier": TAIL_WEIGHT_MULTIPLIER,
                "midgap_interval_ev": [MIDGAP_LOWER_EV, MIDGAP_UPPER_EV],
                "midgap_weight_multiplier": MIDGAP_WEIGHT_MULTIPLIER,
                "device": "cuda" if used_gpu else "cpu",
            },
            indent=2,
            default=json_default,
        ),
        encoding="utf-8",
    )
    metric_rows = []
    prediction_frames = []
    best_iteration = getattr(model, "best_iteration", np.nan)
    for split in ("train", "validation", "test"):
        score = regression_metrics(y_sets[split], predictions[split])
        metric_rows.append(
            {
                "target": target,
                "task": task_name,
                "representation": representation,
                "model": MODEL_NAME,
                "split": split,
                **score,
                "n_samples": len(y_sets[split]),
                "feature_dimension": x_sets[split].shape[1],
                "training_seconds": elapsed,
                "used_gpu": used_gpu,
                "best_iteration": best_iteration,
                "weighting_mode": weighting_mode,
            }
        )
        prediction_frames.append(
            pd.DataFrame(
                {
                    "material_id": ids_sets[split],
                    "split": split,
                    "y_true": y_sets[split],
                    "y_pred": predictions[split],
                    "residual": predictions[split] - y_sets[split],
                    "sample_weight": sample_weights[split],
                    "is_tail": (y_sets[split] >= TAIL_THRESHOLD_EV).astype(np.int8),
                }
            )
        )
    metric_frame = pd.DataFrame(metric_rows)
    prediction_frame = pd.concat(prediction_frames, ignore_index=True)
    metric_frame.to_csv(run_dir / "metrics.csv", index=False)
    save_dat(metric_frame, run_dir / "metrics.dat")
    prediction_frame.to_csv(run_dir / "predictions.csv", index=False)
    save_dat(prediction_frame, run_dir / "predictions.dat")
    logger.info(
        "%s | %s | %s | xgboost | %.1fs | train/val/test RMSE = %.6g / %.6g / %.6g",
        target,
        task_name,
        representation,
        elapsed,
        *[
            float(metric_frame.loc[metric_frame["split"] == split, "rmse"].iloc[0])
            for split in ("train", "validation", "test")
        ],
    )
    for row in metric_rows:
        logger.info(
            "METRICS target=%s task=%s representation=%s model=xgboost split=%s n=%d MAE=%.8g RMSE=%.8g R2=%.8g tail_n=%d tail_RMSE=%.8g tail_bias=%.8g",
            row["target"],
            row["task"],
            row["representation"],
            row["split"],
            row["n_samples"],
            row["mae"],
            row["rmse"],
            row["r2"],
            row["n_tail"],
            row["tail_rmse"],
            row["tail_bias"],
        )
    return {
        "metrics": metric_frame,
        "predictions": prediction_frame,
        "history": history,
        "model_object": model,
        "used_gpu": used_gpu,
        "weighting_mode": weighting_mode,
    }


def result_validation_rmse(result: dict[str, Any]) -> float:
    row = result["metrics"].loc[result["metrics"]["split"] == "validation", "rmse"]
    return float(row.iloc[0])


def materialize_selected_regression(
    source: dict[str, Any],
    source_representation: str,
    selected_representation: str,
    selected_weighting: str,
    target: str,
    output: Path,
) -> dict[str, Any]:
    soap_source = SOAP_LOCAL_CHEMICAL_SOURCE.get(
        source_representation,
        SOAP_COMPOSITION_SOURCE.get(source_representation, source_representation),
    )
    soap_alpha = SOAP_ALPHA_BY_NAME.get(soap_source)
    metrics = source["metrics"].copy()
    if "source_task" not in metrics:
        metrics["source_task"] = metrics["task"]
    metrics["task"] = "direct_regression"
    metrics["source_representation"] = source_representation
    metrics["representation"] = selected_representation
    metrics["selected_weighting"] = selected_weighting
    if soap_alpha is not None:
        metrics["soap_alpha"] = soap_alpha

    predictions = source["predictions"].copy()
    predictions["source_representation"] = source_representation
    predictions["selected_weighting"] = selected_weighting
    if soap_alpha is not None:
        predictions["soap_alpha"] = soap_alpha

    history = source["history"].copy()
    history["source_representation"] = source_representation
    history["selected_weighting"] = selected_weighting
    if soap_alpha is not None:
        history["soap_alpha"] = soap_alpha

    blend_fraction = source.get("tail_weighted_fraction", np.nan)
    if np.isfinite(blend_fraction):
        metrics["tail_weighted_fraction"] = float(blend_fraction)
        predictions["tail_weighted_fraction"] = float(blend_fraction)
        history["tail_weighted_fraction"] = float(blend_fraction)

    run_dir = ensure_dir(
        output
        / "models"
        / target
        / "direct_regression"
        / selected_representation
        / MODEL_NAME
    )
    metrics.to_csv(run_dir / "metrics.csv", index=False)
    save_dat(metrics, run_dir / "metrics.dat")
    predictions.to_csv(run_dir / "predictions.csv", index=False)
    save_dat(predictions, run_dir / "predictions.dat")
    history.to_csv(run_dir / "rmse_history.csv", index=False)
    save_dat(history, run_dir / "rmse_history.dat")
    (run_dir / "selection.json").write_text(
        json.dumps(
            {
                "criterion": "minimum formula-disjoint validation RMSE",
                "source_representation": source_representation,
                "selected_weighting": selected_weighting,
                "source_task": str(source["metrics"]["task"].iloc[0]),
                "tail_weighted_fraction": (
                    float(blend_fraction) if np.isfinite(blend_fraction) else None
                ),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return {
        **source,
        "metrics": metrics,
        "predictions": predictions,
        "history": history,
        "source_representation": source_representation,
        "selected_weighting": selected_weighting,
        "tail_weighted_fraction": blend_fraction,
    }


def materialize_validation_weight_blend(
    weighted: dict[str, Any],
    unweighted: dict[str, Any],
    representation: str,
    target: str,
    output: Path,
) -> dict[str, Any]:
    weighted_predictions = (
        weighted["predictions"]
        .sort_values(["split", "material_id"])
        .reset_index(drop=True)
    )
    unweighted_predictions = (
        unweighted["predictions"]
        .sort_values(["split", "material_id"])
        .reset_index(drop=True)
    )
    if not np.array_equal(
        weighted_predictions[["split", "material_id"]].to_numpy(),
        unweighted_predictions[["split", "material_id"]].to_numpy(),
    ):
        raise RuntimeError(
            f"Weighted/unweighted prediction mismatch for {representation}"
        )
    if not np.allclose(
        weighted_predictions["y_true"], unweighted_predictions["y_true"]
    ):
        raise RuntimeError(f"Weighted/unweighted targets mismatch for {representation}")
    validation_mask = weighted_predictions["split"].to_numpy() == "validation"
    y_validation = weighted_predictions.loc[validation_mask, "y_true"].to_numpy(
        dtype=float
    )
    weighted_validation = weighted_predictions.loc[validation_mask, "y_pred"].to_numpy(
        dtype=float
    )
    unweighted_validation = unweighted_predictions.loc[
        validation_mask, "y_pred"
    ].to_numpy(dtype=float)
    blend_weight = min(
        WEIGHT_BLEND_GRID,
        key=lambda value: float(
            np.sqrt(
                np.mean(
                    (
                        value * weighted_validation
                        + (1.0 - value) * unweighted_validation
                        - y_validation
                    )
                    ** 2
                )
            )
        ),
    )
    predictions = weighted_predictions.copy()
    predictions["y_pred"] = blend_weight * weighted_predictions["y_pred"].to_numpy(
        dtype=float
    ) + (1.0 - blend_weight) * unweighted_predictions["y_pred"].to_numpy(dtype=float)
    predictions["residual"] = predictions["y_pred"] - predictions["y_true"]
    predictions["source_representation"] = representation
    predictions["selected_weighting"] = "validation_blend"
    predictions["tail_weighted_fraction"] = blend_weight
    weighted_metrics = weighted["metrics"].set_index("split")
    unweighted_metrics = unweighted["metrics"].set_index("split")
    metric_rows: list[dict[str, Any]] = []
    for split in ("train", "validation", "test"):
        subset = predictions.loc[predictions["split"] == split]
        score = regression_metrics(subset["y_true"], subset["y_pred"])
        metric_rows.append(
            {
                "target": target,
                "task": "direct_regression",
                "source_task": "validation_blend_of_weighted_and_unweighted",
                "representation": representation,
                "source_representation": representation,
                "model": MODEL_NAME,
                "split": split,
                **score,
                "n_samples": len(subset),
                "feature_dimension": int(
                    weighted_metrics.loc[split, "feature_dimension"]
                ),
                "training_seconds": float(
                    weighted_metrics.loc[split, "training_seconds"]
                    + unweighted_metrics.loc[split, "training_seconds"]
                ),
                "used_gpu": bool(
                    weighted_metrics.loc[split, "used_gpu"]
                    or unweighted_metrics.loc[split, "used_gpu"]
                ),
                "best_iteration": np.nan,
                "selected_weighting": "validation_blend",
                "tail_weighted_fraction": blend_weight,
            }
        )
    metrics = pd.DataFrame(metric_rows)
    history_source = weighted if blend_weight >= 0.5 else unweighted
    history = history_source["history"].copy()
    history["source_representation"] = representation
    history["selected_weighting"] = "validation_blend"
    history["tail_weighted_fraction"] = blend_weight
    run_dir = ensure_dir(
        output / "models" / target / "direct_regression" / representation / MODEL_NAME
    )
    metrics.to_csv(run_dir / "metrics.csv", index=False)
    save_dat(metrics, run_dir / "metrics.dat")
    predictions.to_csv(run_dir / "predictions.csv", index=False)
    save_dat(predictions, run_dir / "predictions.dat")
    history.to_csv(run_dir / "rmse_history.csv", index=False)
    save_dat(history, run_dir / "rmse_history.dat")
    (run_dir / "selection.json").write_text(
        json.dumps(
            {
                "criterion": "minimum validation RMSE over convex prediction blend",
                "tail_weighted_fraction": blend_weight,
                "grid": WEIGHT_BLEND_GRID,
            },
            indent=2,
            default=json_default,
        ),
        encoding="utf-8",
    )
    return {
        "metrics": metrics,
        "predictions": predictions,
        "history": history,
        "model_object": None,
        "used_gpu": bool(metrics["used_gpu"].iloc[0]),
        "source_representation": representation,
        "selected_weighting": "validation_blend",
        "tail_weighted_fraction": blend_weight,
    }


def materialize_validation_expert_blend(
    experts: dict[str, dict[str, Any]],
    representation: str,
    target: str,
    output: Path,
    protect_midgap: bool,
) -> dict[str, Any]:
    names = tuple(experts)
    ordered = {
        name: experts[name]["predictions"]
        .sort_values(["split", "material_id"])
        .reset_index(drop=True)
        for name in names
    }
    reference = ordered[names[0]]
    for name in names[1:]:
        candidate = ordered[name]
        if not np.array_equal(
            reference[["split", "material_id"]].to_numpy(),
            candidate[["split", "material_id"]].to_numpy(),
        ) or not np.allclose(reference["y_true"], candidate["y_true"]):
            raise RuntimeError(
                f"Expert prediction mismatch for {representation}: {name}"
            )
    validation_mask = reference["split"].to_numpy() == "validation"
    y_validation = reference.loc[validation_mask, "y_true"].to_numpy(dtype=float)
    weights = validation_simplex_weights(
        y_validation,
        [
            ordered[name].loc[validation_mask, "y_pred"].to_numpy(dtype=float)
            for name in names
        ],
        EXPERT_BLEND_L2,
        protect_midgap,
    )
    predictions = reference.copy()
    predictions["y_pred"] = sum(
        weight * ordered[name]["y_pred"].to_numpy(dtype=float)
        for name, weight in zip(names, weights)
    )
    predictions["residual"] = predictions["y_pred"] - predictions["y_true"]
    predictions["source_representation"] = representation
    predictions["selected_weighting"] = "validation_expert_blend"
    for name, weight in zip(names, weights):
        predictions[f"expert_weight_{name}"] = float(weight)

    metric_rows: list[dict[str, Any]] = []
    total_seconds = sum(
        float(experts[name]["metrics"]["training_seconds"].iloc[0]) for name in names
    )
    feature_dimension = int(experts[names[0]]["metrics"]["feature_dimension"].iloc[0])
    for split in ("train", "validation", "test"):
        subset = predictions.loc[predictions["split"] == split]
        score = regression_metrics(subset["y_true"], subset["y_pred"])
        row = {
            "target": target,
            "task": "direct_regression",
            "source_task": "validation_blend_of_regression_experts",
            "representation": representation,
            "source_representation": representation,
            "model": MODEL_NAME,
            "split": split,
            **score,
            "n_samples": len(subset),
            "feature_dimension": feature_dimension,
            "training_seconds": total_seconds,
            "used_gpu": any(bool(experts[name]["used_gpu"]) for name in names),
            "best_iteration": np.nan,
            "selected_weighting": "validation_expert_blend",
            "tail_weighted_fraction": float(
                weights[names.index("tail_weighted")]
                if "tail_weighted" in names
                else 0.0
            ),
            "midgap_weighted_fraction": float(
                weights[names.index("midgap_weighted")]
                if "midgap_weighted" in names
                else 0.0
            ),
        }
        for name, weight in zip(names, weights):
            row[f"expert_weight_{name}"] = float(weight)
        metric_rows.append(row)
    metrics = pd.DataFrame(metric_rows)
    history_source_name = names[int(np.argmax(weights))]
    history = experts[history_source_name]["history"].copy()
    history["source_representation"] = representation
    history["selected_weighting"] = "validation_expert_blend"
    history["history_source_expert"] = history_source_name
    for name, weight in zip(names, weights):
        history[f"expert_weight_{name}"] = float(weight)

    run_dir = ensure_dir(
        output / "models" / target / "direct_regression" / representation / MODEL_NAME
    )
    metrics.to_csv(run_dir / "metrics.csv", index=False)
    save_dat(metrics, run_dir / "metrics.dat")
    predictions.to_csv(run_dir / "predictions.csv", index=False)
    save_dat(predictions, run_dir / "predictions.dat")
    history.to_csv(run_dir / "rmse_history.csv", index=False)
    save_dat(history, run_dir / "rmse_history.dat")
    selection = {
        "criterion": "validation-only convex expert blend",
        "protect_midgap": protect_midgap,
        "midgap_interval_ev": [MIDGAP_LOWER_EV, MIDGAP_UPPER_EV],
        "midgap_objective_weight": MIDGAP_OBJECTIVE_WEIGHT if protect_midgap else 0.0,
        "catastrophic_error_threshold_ev": CATASTROPHIC_ERROR_THRESHOLD_EV,
        "catastrophic_objective_weight": (
            CATASTROPHIC_OBJECTIVE_WEIGHT if protect_midgap else 0.0
        ),
        "l2_strength": EXPERT_BLEND_L2,
        "expert_weights": dict(zip(names, map(float, weights))),
    }
    (run_dir / "selection.json").write_text(
        json.dumps(selection, indent=2, default=json_default), encoding="utf-8"
    )
    return {
        "metrics": metrics,
        "predictions": predictions,
        "history": history,
        "model_object": None,
        "component_models": {name: experts[name].get("model_object") for name in names},
        "expert_results": experts,
        "expert_weights": dict(zip(names, map(float, weights))),
        "used_gpu": bool(metrics["used_gpu"].iloc[0]),
        "source_representation": representation,
        "selected_weighting": "validation_expert_blend",
        "tail_weighted_fraction": selection["expert_weights"].get("tail_weighted", 0.0),
        "midgap_weighted_fraction": selection["expert_weights"].get(
            "midgap_weighted", 0.0
        ),
    }


def materialize_soed_candidate_ensemble(
    candidates: dict[str, dict[str, Any]],
    target: str,
    output: Path,
) -> dict[str, Any]:
    names = tuple(name for name in SOED_ENHANCED_CANDIDATES if name in candidates)
    ordered = {
        name: candidates[name]["predictions"]
        .sort_values(["split", "material_id"])
        .reset_index(drop=True)
        for name in names
    }
    reference = ordered[names[0]]
    for name in names[1:]:
        candidate = ordered[name]
        if not np.array_equal(
            reference[["split", "material_id"]].to_numpy(),
            candidate[["split", "material_id"]].to_numpy(),
        ) or not np.allclose(reference["y_true"], candidate["y_true"]):
            raise RuntimeError(f"SOED candidate prediction mismatch: {name}")
    validation_mask = reference["split"].to_numpy() == "validation"
    y_validation = reference.loc[validation_mask, "y_true"].to_numpy(dtype=float)
    weights = validation_simplex_weights(
        y_validation,
        [
            ordered[name].loc[validation_mask, "y_pred"].to_numpy(dtype=float)
            for name in names
        ],
        SOED_ENSEMBLE_L2,
        True,
    )
    if SOED_ENSEMBLE_MIN_WEIGHT > 0.0:
        weights = np.maximum(weights, SOED_ENSEMBLE_MIN_WEIGHT)
        weights /= weights.sum()
    predictions = reference.copy()
    predictions["y_pred"] = sum(
        weight * ordered[name]["y_pred"].to_numpy(dtype=float)
        for name, weight in zip(names, weights)
    )
    predictions["residual"] = predictions["y_pred"] - predictions["y_true"]
    predictions["source_representation"] = "prediction_ensemble"
    predictions["selected_weighting"] = "validation_soed_ensemble"
    for name, weight in zip(names, weights):
        predictions[f"candidate_weight_{name}"] = float(weight)

    maximum_dimension = max(
        int(candidates[name]["metrics"]["feature_dimension"].iloc[0]) for name in names
    )
    total_seconds = sum(
        float(candidates[name]["metrics"]["training_seconds"].iloc[0]) for name in names
    )
    metric_rows: list[dict[str, Any]] = []
    for split in ("train", "validation", "test"):
        subset = predictions.loc[predictions["split"] == split]
        row = {
            "target": target,
            "task": "direct_regression",
            "source_task": "validation_ensemble_of_soed_candidates",
            "representation": PRIMARY_SOED_NAME,
            "source_representation": "prediction_ensemble",
            "model": MODEL_NAME,
            "split": split,
            **regression_metrics(subset["y_true"], subset["y_pred"]),
            "n_samples": len(subset),
            "feature_dimension": maximum_dimension,
            "training_seconds": total_seconds,
            "used_gpu": any(bool(candidates[name]["used_gpu"]) for name in names),
            "best_iteration": np.nan,
            "selected_weighting": "validation_soed_ensemble",
            "ensemble_members": len(names),
        }
        for name, weight in zip(names, weights):
            row[f"candidate_weight_{name}"] = float(weight)
        metric_rows.append(row)
    metrics = pd.DataFrame(metric_rows)
    history_source = names[int(np.argmax(weights))]
    history = candidates[history_source]["history"].copy()
    history["source_representation"] = history_source
    history["selected_weighting"] = "validation_soed_ensemble"
    history["history_source_candidate"] = history_source
    for name, weight in zip(names, weights):
        history[f"candidate_weight_{name}"] = float(weight)
    run_dir = ensure_dir(
        output
        / "models"
        / target
        / "direct_regression"
        / PRIMARY_SOED_NAME
        / MODEL_NAME
    )
    metrics.to_csv(run_dir / "metrics.csv", index=False)
    save_dat(metrics, run_dir / "metrics.dat")
    predictions.to_csv(run_dir / "predictions.csv", index=False)
    save_dat(predictions, run_dir / "predictions.dat")
    history.to_csv(run_dir / "rmse_history.csv", index=False)
    save_dat(history, run_dir / "rmse_history.dat")
    selection = {
        "criterion": "validation-only regularized convex SOED candidate ensemble",
        "candidate_weights": dict(zip(names, map(float, weights))),
        "l2_strength": SOED_ENSEMBLE_L2,
        "midgap_protection": True,
        "test_set_used_for_selection": False,
    }
    (run_dir / "selection.json").write_text(
        json.dumps(selection, indent=2), encoding="utf-8"
    )
    save_dat(
        pd.DataFrame(
            {
                "candidate": names,
                "weight": weights,
                "label": [SOED_CANDIDATE_LABELS[name] for name in names],
            }
        ),
        run_dir / "candidate_ensemble_weights.dat",
    )
    return {
        "metrics": metrics,
        "predictions": predictions,
        "history": history,
        "model_object": None,
        "candidate_results": candidates,
        "candidate_weights": selection["candidate_weights"],
        "used_gpu": bool(metrics["used_gpu"].iloc[0]),
        "source_representation": "prediction_ensemble",
        "selected_weighting": "validation_soed_ensemble",
    }


def fit_one_classifier(
    representation: str,
    target: str,
    x_sets: dict[str, np.ndarray],
    y_sets: dict[str, np.ndarray],
    ids_sets: dict[str, np.ndarray],
    hardware: dict[str, Any],
    output: Path,
    logger: logging.Logger,
) -> dict[str, Any]:
    run_dir = ensure_dir(
        output
        / "models"
        / target
        / "metal_nonmetal_classification"
        / representation
        / MODEL_NAME
    )
    params = tune_classifier(
        representation,
        x_sets["train"],
        y_sets["train"],
        x_sets["validation"],
        y_sets["validation"],
        hardware,
        output,
        logger,
    )
    started = time.perf_counter()
    used_gpu = bool(hardware.get("gpu_available"))
    try:
        model, history, structure = fit_xgb_classifier(
            x_sets["train"],
            y_sets["train"],
            x_sets["validation"],
            y_sets["validation"],
            params,
            used_gpu,
            logger,
        )
    except Exception as exc:
        if not used_gpu:
            raise
        logger.warning("XGBoost GPU classifier failed (%r); retrying on CPU", exc)
        used_gpu = False
        model, history, structure = fit_xgb_classifier(
            x_sets["train"],
            y_sets["train"],
            x_sets["validation"],
            y_sets["validation"],
            params,
            False,
            logger,
        )
    probabilities = {
        split: xgb_predict(model, values) for split, values in x_sets.items()
    }
    elapsed = time.perf_counter() - started
    threshold = select_classifier_threshold(
        y_sets["validation"], probabilities["validation"]
    )
    joblib.dump(model, run_dir / "model.joblib")
    history.to_csv(run_dir / "logloss_history.csv", index=False)
    save_dat(history, run_dir / "logloss_history.dat")
    (run_dir / "model_structure.txt").write_text(structure, encoding="utf-8")
    metric_rows = []
    prediction_frames = []
    best_iteration = getattr(model, "best_iteration", np.nan)
    for split in ("train", "validation", "test"):
        score = classification_metrics(y_sets[split], probabilities[split], threshold)
        metric_rows.append(
            {
                "target": target,
                "task": "metal_nonmetal_classification",
                "representation": representation,
                "model": MODEL_NAME,
                "split": split,
                **score,
                "n_samples": len(y_sets[split]),
                "feature_dimension": x_sets[split].shape[1],
                "training_seconds": elapsed,
                "used_gpu": used_gpu,
                "best_iteration": best_iteration,
                "gap_threshold_ev": BAND_GAP_ZERO_THRESHOLD_EV,
                "probability_threshold": threshold,
            }
        )
        prediction_frames.append(
            pd.DataFrame(
                {
                    "material_id": ids_sets[split],
                    "split": split,
                    "y_true_nonmetal": y_sets[split],
                    "p_nonmetal": probabilities[split],
                    "y_pred_nonmetal": (probabilities[split] >= threshold).astype(
                        np.int8
                    ),
                }
            )
        )
    metric_frame = pd.DataFrame(metric_rows)
    prediction_frame = pd.concat(prediction_frames, ignore_index=True)
    metric_frame.to_csv(run_dir / "metrics.csv", index=False)
    save_dat(metric_frame, run_dir / "metrics.dat")
    prediction_frame.to_csv(run_dir / "predictions.csv", index=False)
    save_dat(prediction_frame, run_dir / "predictions.dat")
    for row in metric_rows:
        logger.info(
            "CLASSIFICATION target=%s representation=%s model=xgboost split=%s n=%d AUROC=%.8g AUPRC=%.8g F1=%.8g balanced_accuracy=%.8g MCC=%.8g Brier=%.8g ECE=%.8g threshold=%.4f",
            row["target"],
            row["representation"],
            row["split"],
            row["n_samples"],
            row["auroc"],
            row["auprc"],
            row["f1"],
            row["balanced_accuracy"],
            row["mcc"],
            row["brier"],
            row["ece_10bin"],
            threshold,
        )
    return {
        "metrics": metric_frame,
        "predictions": prediction_frame,
        "probabilities": probabilities,
        "decision_threshold": threshold,
        "history": history,
    }


REPRESENTATION_LABELS = {
    PRIMARY_SOAP_NAME: "Validation-selected periodic SOAP",
    COMPOSITION_ONLY_NAME: "Composition-only control",
    PRIMARY_CHEMICAL_SOAP_NAME: "Periodic SOAP + composition",
    PRIMARY_LOCAL_CHEMICAL_SOAP_NAME: "Validation-selected local-chemical SOAP",
    DENSITY_COMPOSITION_NAME: r"SOED ($\rho$) + composition",
    "psoed_density": r"SOED ($\rho$)",
    MATCHED_CHEMICAL_SOED_NAME: r"Local-chemical SOED ($\rho$)",
    "psoed_log_density": r"SOED ($\log\rho$)",
    "psoed_density_contrast": "SOED (density contrast)",
    "psoed_density_gradient": r"SOED ($\rho,|\nabla\rho|$)",
    "psoed_physics_multichannel": "Multichannel SOED",
    LEGACY_CHEMICAL_SOED_NAME: "Cross-channel chemical SOED",
    COMPACT_MULTICHANNEL_SOED_NAME: "Compact six-channel SOED",
    PRIMARY_SOED_NAME: "Validation-selected electronic SOED",
    **SOED_CANDIDATE_LABELS,
}


def representation_label(name: str) -> str:
    if name in SOAP_ALPHA_BY_NAME:
        return rf"Periodic SOAP ($\alpha={SOAP_ALPHA_BY_NAME[name]:g}$)"
    if name in SOAP_COMPOSITION_SOURCE:
        source = SOAP_COMPOSITION_SOURCE[name]
        return rf"Periodic SOAP + composition ($\alpha={SOAP_ALPHA_BY_NAME[source]:g}$)"
    if name in SOAP_LOCAL_CHEMICAL_SOURCE:
        source = SOAP_LOCAL_CHEMICAL_SOURCE[name]
        return rf"Local-chemical SOAP ($\alpha={SOAP_ALPHA_BY_NAME[source]:g}$)"
    return REPRESENTATION_LABELS.get(name, name.replace("_", " "))


def representation_color(name: str) -> str:
    if name == COMPOSITION_ONLY_NAME:
        return "tab:gray"
    if (
        name in SOAP_REPRESENTATIONS
        or name in SOAP_COMPOSITION_REPRESENTATIONS
        or name
        in (
            PRIMARY_SOAP_NAME,
            PRIMARY_CHEMICAL_SOAP_NAME,
        )
    ):
        return "tab:purple"
    if (
        name in SOAP_LOCAL_CHEMICAL_REPRESENTATIONS
        or name == PRIMARY_LOCAL_CHEMICAL_SOAP_NAME
    ):
        return "tab:pink"
    colors = {
        MATCHED_SOED_NAME: "tab:orange",
        MATCHED_CHEMICAL_SOED_NAME: "tab:cyan",
        DENSITY_COMPOSITION_NAME: "tab:orange",
        "psoed_log_density": "tab:green",
        "psoed_density_contrast": "tab:cyan",
        "psoed_density_gradient": "tab:red",
        "psoed_physics_multichannel": "tab:blue",
        LEGACY_CHEMICAL_SOED_NAME: "tab:olive",
        "psoed_density_local_chemical_global": "tab:blue",
        "psoed_density_local_chemical_radial": "tab:green",
        "psoed_density_local_chemical_radial_global": "tab:red",
        COMPACT_MULTICHANNEL_SOED_NAME: "tab:olive",
        PRIMARY_SOED_NAME: "tab:brown",
    }
    return colors.get(name, "tab:gray")


def configure_matplotlib() -> None:
    import matplotlib as mpl

    mpl.use("Agg")
    mpl.rcParams["font.family"] = "sans-serif"
    mpl.rcParams["font.sans-serif"] = list(FONT_FAMILY_PRIORITY)
    mpl.rcParams["axes.unicode_minus"] = False
    mpl.rcParams["axes.linewidth"] = 1.0
    mpl.rcParams["axes.spines.top"] = False
    mpl.rcParams["axes.spines.right"] = False
    mpl.rcParams["legend.frameon"] = False
    mpl.rcParams["xtick.labelsize"] = TICK_LABEL_FONT_SIZE
    mpl.rcParams["ytick.labelsize"] = TICK_LABEL_FONT_SIZE


def save_jpg(fig, path: Path) -> None:
    ensure_dir(path.parent)
    fig.patch.set_facecolor("white")
    fig.savefig(
        path, dpi=PLOT_DPI, format="jpg", bbox_inches="tight", facecolor="white"
    )


def create_workflow_files(output: Path) -> None:
    configure_matplotlib()
    import matplotlib.pyplot as plt
    from matplotlib.colors import to_hex
    from matplotlib.patches import FancyBboxPatch

    main_dir = ensure_dir(output / "main_figures")
    ratio = ":".join(
        f"{100 * value:g}" for value in (TRAIN_RATIO, VALID_RATIO, TEST_RATIO)
    )
    nodes = [
        (
            "A",
            "MP-20-Charge\nCIF + $\\rho_e(\\mathbf{r})$",
            0.07,
            0.50,
            TAB_COLORS["train"],
        ),
        ("Q", "Quality control\nPauling $X$ defined", 0.23, 0.50, "tab:olive"),
        (
            "B",
            f"Fresh formula split\nseed {RANDOM_SEED}; {ratio}",
            0.39,
            0.50,
            TAB_COLORS["validation"],
        ),
        (
            "C",
            "Local-chemical SOAP\nvalidation-selected $\\alpha$",
            0.56,
            0.68,
            representation_color(PRIMARY_LOCAL_CHEMICAL_SOAP_NAME),
        ),
        (
            "D",
            "Electronic SOED\ncandidate family\n$\\rho$ + radial/global",
            0.56,
            0.32,
            representation_color(MATCHED_CHEMICAL_SOED_NAME),
        ),
        ("G", "XGBoost experts\nearly stopping", 0.73, 0.50, TAB_COLORS["direct"]),
        (
            "F",
            "Validation-only\nselection + ensemble",
            0.86,
            0.32,
            representation_color(PRIMARY_SOED_NAME),
        ),
        ("H", "Frozen outer test\n5 fresh seeds", 0.94, 0.62, TAB_COLORS["test"]),
    ]
    edges = [
        ("A", "Q"),
        ("Q", "B"),
        ("B", "C"),
        ("B", "D"),
        ("C", "G"),
        ("D", "G"),
        ("G", "F"),
        ("F", "H"),
    ]
    lookup = {node[0]: node for node in nodes}
    fig, ax = plt.subplots(figsize=(9.0, 3.6))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    width, height = 0.14, 0.15
    for node_id, label, x_value, y_value, color in nodes:
        ax.add_patch(
            FancyBboxPatch(
                (x_value - width / 2, y_value - height / 2),
                width,
                height,
                boxstyle="round,pad=0.012,rounding_size=0.015",
                facecolor=color,
                edgecolor="tab:gray",
                linewidth=1.1,
                alpha=0.90,
            )
        )
        ax.text(
            x_value,
            y_value,
            label,
            ha="center",
            va="center",
            fontsize=6.4,
            color="white",
            fontweight="bold",
        )
    for source, target in edges:
        curvature = 0.0
        source_x, source_y = lookup[source][2], lookup[source][3]
        target_x, target_y = lookup[target][2], lookup[target][3]
        delta_x, delta_y = target_x - source_x, target_y - source_y
        source_scale = min(
            width / (2.0 * max(abs(delta_x), 1e-12)),
            height / (2.0 * max(abs(delta_y), 1e-12)),
        )
        target_scale = source_scale
        start_point = (
            source_x + source_scale * delta_x,
            source_y + source_scale * delta_y,
        )
        end_point = (
            target_x - target_scale * delta_x,
            target_y - target_scale * delta_y,
        )
        ax.annotate(
            "",
            xy=end_point,
            xytext=start_point,
            arrowprops=dict(
                arrowstyle="-|>",
                lw=1.3,
                color="tab:gray",
                shrinkA=0,
                shrinkB=0,
                connectionstyle=f"arc3,rad={curvature}",
            ),
        )
    ax.set_title(
        "Frozen-Confirmatory Electronic SOED Evaluation", fontsize=TITLE_FONT_SIZE - 1
    )
    save_jpg(fig, main_dir / "Figure01_workflow.jpg")
    plt.close(fig)
    rows = [
        {
            "type": "node",
            "source": node[0],
            "target": "",
            "label": node[1],
            "x": node[2],
            "y": node[3],
            "color": node[4],
        }
        for node in nodes
    ] + [
        {
            "type": "edge",
            "source": source,
            "target": target,
            "label": "",
            "x": "",
            "y": "",
            "color": "tab:gray",
        }
        for source, target in edges
    ]
    save_dat(pd.DataFrame(rows), main_dir / "Figure01_workflow.dat")
    mermaid = f"""flowchart LR
    A["<b>MP-20-Charge</b><br/>CIF + <i>ρ</i><sub>e</sub>(<b><i>r</i></b>)"]
    Q["<b>Quality control</b><br/>defined Pauling electronegativity"]
    B["<b>Fresh formula-grouped split</b><br/>seed {RANDOM_SEED}; {ratio}"]
    C["<b>Local-chemical SOAP</b><br/>validation-selected <i>α</i>"]
    D["<b>Electronic SOED candidates</b><br/><i>ρ</i> + radial/global"]
    G["<b>XGBoost experts</b><br/>early stopping"]
    F["<b>Validation-only selection</b><br/>SOED candidate ensemble"]
    H["<b>Frozen outer test</b><br/>five fresh seeds"]
    A --> Q
    Q --> B
    B --> C
    B --> D
    C --> G
    D --> G
    G --> F
    F --> H
"""
    (main_dir / "Figure01_workflow_mermaid.txt").write_text(mermaid, encoding="utf-8")
    mxfile = ET.Element("mxfile", host="app.diagrams.net")
    diagram = ET.SubElement(mxfile, "diagram", id="psoed", name="Page-1")
    model = ET.SubElement(
        diagram,
        "mxGraphModel",
        dx="1200",
        dy="800",
        grid="1",
        page="1",
        pageWidth="1169",
        pageHeight="827",
    )
    root = ET.SubElement(model, "root")
    ET.SubElement(root, "mxCell", id="0")
    ET.SubElement(root, "mxCell", id="1", parent="0")
    html_labels = {
        "A": "<b>MP-20-Charge</b><br>CIF + <i>ρ</i><sub>e</sub>(<b><i>r</i></b>)",
        "Q": "<b>Quality control</b><br>defined Pauling electronegativity",
        "B": f"<b>Fresh formula-grouped split</b><br>seed {RANDOM_SEED}; {ratio}",
        "C": "<b>Local-chemical SOAP</b><br>validation-selected <i>α</i>",
        "D": "<b>Electronic SOED candidates</b><br><i>ρ</i> + radial/global",
        "G": "<b>XGBoost experts</b><br>early stopping",
        "F": "<b>Validation-only selection</b><br>SOED candidate ensemble",
        "H": "<b>Frozen outer test</b><br>five fresh seeds",
    }
    for node_id, _, x_value, y_value, color in nodes:
        cell = ET.SubElement(
            root,
            "mxCell",
            id=node_id,
            value=html_labels[node_id],
            style=f"rounded=1;whiteSpace=wrap;html=1;fillColor={to_hex(color)};fontColor=#FFFFFF;strokeColor=#333333;fontSize=14;",
            vertex="1",
            parent="1",
        )
        ET.SubElement(
            cell,
            "mxGeometry",
            x=str(int(x_value * 950)),
            y=str(int((1 - y_value) * 600)),
            width="185",
            height="75",
            **{"as": "geometry"},
        )
    for edge_index, (source, target) in enumerate(edges, start=1):
        cell = ET.SubElement(
            root,
            "mxCell",
            id=f"edge{edge_index}",
            style="edgeStyle=orthogonalEdgeStyle;html=1;endArrow=block;strokeWidth=2;",
            edge="1",
            parent="1",
            source=source,
            target=target,
        )
        ET.SubElement(cell, "mxGeometry", relative="1", **{"as": "geometry"})
    ET.ElementTree(mxfile).write(
        main_dir / "Figure01_workflow.drawio", encoding="utf-8", xml_declaration=True
    )
    (output / "methodology_notes.txt").write_text(
        "Periodic SOAP implementation\n"
        "- No DScribe dependency is used.\n"
        "- The scaffold density follows the equal-width Gaussian construction in Eqs. 13-14 of Gugler and Reiher, JCTC 2022.\n"
        "- Periodicity is imposed by reciprocal-space Gaussian convolution on the native charge-density grid.\n"
        "- SOAP and the matched density-only SOED use identical centers, radial basis functions, angular quadrature, spherical-harmonic cutoff, invariant spectrum, local mean/std pooling, and descriptor dimension.\n"
        "- SOAP is a baseline only; no SOAP-SOED feature fusion is formed.\n"
        "- The primary fair comparison augments SOAP and density SOED with identical local element-property contrast projections, s/p/d/f block projections, and a 65-dimensional composition block.\n"
        "- The v13 SOED family starts from the local-chemical density SOED and tests compact, density-specific additions: atom-centered radial electron distributions, charge/nuclear imbalance summaries, and global/reciprocal density statistics.\n"
        "- A predefined 1-2 eV expert is trained for both the fair local-chemical SOAP family and the SOED candidates; expert weights use validation data only.\n"
        "- The primary SOED is a regularized convex prediction ensemble of predefined SOED candidates; it contains no SOAP prediction or SOAP feature.\n"
        "- Composition-only and simple composition-appended controls are retained in the SI; SOAP and SOED descriptors are never concatenated.\n"
        "- SOAP Gaussian width, regression-expert weights, and SOED candidate weights are selected only from formula-disjoint validation data.\n"
        f"- The confirmatory main split uses fresh seed {RANDOM_SEED}; five fresh grouped outer seeds repeat every validation choice after the algorithm is frozen. Legacy seed {EXPLORATORY_LEGACY_SEED} is not reused.\n"
        "- Statistical intervals resample reduced-formula groups and gap-bin p-values are adjusted by Benjamini-Hochberg.\n"
        "- Ensemble-weighted XGBoost gain is reported by physically defined feature block as an attribution diagnostic; because tree gain can favor large blocks, it is not interpreted as causal proof.\n"
        "- A fixed invariant vector is used for XGBoost; this is a periodic feature adaptation of the paper's normalized pairwise kernel, not the authors' original analytic kernel code.\n\n"
        "Data-quality and cache policy\n"
        "- Structures containing any element with undefined Pauling electronegativity are excluded before splitting.\n"
        "- Raw CIF and charge-density files are never deleted.\n"
        "- Parsed structures and SOAP/SOED features are cached with parameter-aware fingerprints.\n"
        "- Descriptor times are estimated as standalone postprocessing costs; shared intermediates are not used to make one candidate appear artificially faster.\n"
        "- Changing samples, density-file metadata, descriptor parameters, active elements, or cache schema automatically creates a different descriptor cache.\n",
        encoding="utf-8",
    )


def validation_selected_rows(metrics_frame: pd.DataFrame, task: str) -> pd.DataFrame:
    validation = metrics_frame.loc[
        (metrics_frame["task"] == task) & (metrics_frame["split"] == "validation")
    ].copy()
    test = metrics_frame.loc[
        (metrics_frame["task"] == task) & (metrics_frame["split"] == "test")
    ].copy()
    if validation.empty:
        return pd.DataFrame()
    validation = validation.rename(
        columns={
            "mae": "validation_mae",
            "rmse": "validation_rmse",
            "r2": "validation_r2",
            "tail_mae": "validation_tail_mae",
            "tail_rmse": "validation_tail_rmse",
            "tail_bias": "validation_tail_bias",
        }
    )
    test = test.rename(
        columns={
            "mae": "test_mae",
            "rmse": "test_rmse",
            "r2": "test_r2",
            "tail_mae": "test_tail_mae",
            "tail_rmse": "test_tail_rmse",
            "tail_bias": "test_tail_bias",
        }
    )
    metadata_columns = [
        column
        for column in validation.columns
        if column
        in {
            "selected_weighting",
            "source_representation",
            "soap_alpha",
            "tail_weighted_fraction",
            "midgap_weighted_fraction",
            "ensemble_members",
        }
        or column.startswith("expert_weight_")
        or column.startswith("candidate_weight_")
    ]
    validation_columns = [
        "representation",
        "model",
        "validation_mae",
        "validation_rmse",
        "validation_r2",
        "validation_tail_mae",
        "validation_tail_rmse",
        "validation_tail_bias",
        *metadata_columns,
    ]
    test_columns = [
        "representation",
        "model",
        "test_mae",
        "test_rmse",
        "test_r2",
        "test_tail_mae",
        "test_tail_rmse",
        "test_tail_bias",
        "feature_dimension",
        "training_seconds",
        "used_gpu",
        "best_iteration",
    ]
    return validation[validation_columns].merge(
        test[test_columns], on=["representation", "model"], how="left"
    )


def bootstrap_distribution(
    proposed_predictions: pd.DataFrame,
    baseline_predictions: pd.DataFrame,
    proposed_name: str,
    baseline_name: str,
    comparison_label: str,
    group_lookup: pd.DataFrame,
    analysis_scope: str = "all_test",
    minimum_gap_ev: float | None = None,
    maximum_gap_ev: float | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    proposed = proposed_predictions.loc[
        proposed_predictions["split"] == "test", ["material_id", "y_true", "y_pred"]
    ].rename(columns={"y_pred": "y_pred_soed"})
    baseline = baseline_predictions.loc[
        baseline_predictions["split"] == "test", ["material_id", "y_true", "y_pred"]
    ].rename(columns={"y_true": "y_true_baseline", "y_pred": "y_pred_soap"})
    paired = proposed.merge(
        baseline, on="material_id", how="inner", validate="one_to_one"
    )
    paired = paired.merge(
        group_lookup[["material_id", "reduced_formula"]].drop_duplicates("material_id"),
        on="material_id",
        how="left",
        validate="one_to_one",
    )
    if minimum_gap_ev is not None:
        paired = paired.loc[paired["y_true"] >= minimum_gap_ev].copy()
    if maximum_gap_ev is not None:
        paired = paired.loc[paired["y_true"] < maximum_gap_ev].copy()
    if paired.empty:
        raise RuntimeError(
            f"No paired test samples for analysis scope {analysis_scope}"
        )
    if paired["reduced_formula"].isna().any():
        raise RuntimeError("Missing reduced-formula group in paired bootstrap")
    if not np.allclose(paired["y_true"], paired["y_true_baseline"]):
        raise RuntimeError("SOAP and SOED test targets are not aligned")
    y_true = paired["y_true"].to_numpy(dtype=float)
    y_soed = paired["y_pred_soed"].to_numpy(dtype=float)
    y_soap = paired["y_pred_soap"].to_numpy(dtype=float)
    rmse_soed = float(np.sqrt(np.mean((y_soed - y_true) ** 2)))
    rmse_soap = float(np.sqrt(np.mean((y_soap - y_true) ** 2)))
    paired["squared_error_soed"] = (y_soed - y_true) ** 2
    paired["squared_error_soap"] = (y_soap - y_true) ** 2
    paired["absolute_error_soed"] = np.abs(y_soed - y_true)
    paired["absolute_error_soap"] = np.abs(y_soap - y_true)
    grouped = (
        paired.groupby("reduced_formula", sort=False)
        .agg(
            n_materials=("material_id", "size"),
            squared_error_soed=("squared_error_soed", "sum"),
            squared_error_soap=("squared_error_soap", "sum"),
            absolute_error_soed=("absolute_error_soed", "mean"),
            absolute_error_soap=("absolute_error_soap", "mean"),
        )
        .reset_index()
    )
    n_groups = len(grouped)
    group_n = grouped["n_materials"].to_numpy(dtype=float)
    group_sse_soed = grouped["squared_error_soed"].to_numpy(dtype=float)
    group_sse_soap = grouped["squared_error_soap"].to_numpy(dtype=float)
    rng = np.random.default_rng(RANDOM_SEED)
    deltas = np.empty(BOOTSTRAP_REPEATS, dtype=float)
    for start in range(0, BOOTSTRAP_REPEATS, 250):
        count = min(250, BOOTSTRAP_REPEATS - start)
        draws = rng.integers(0, n_groups, size=(count, n_groups))
        multiplicities = np.asarray(
            [np.bincount(row, minlength=n_groups) for row in draws], dtype=float
        )
        denominators = multiplicities @ group_n
        rmse_soap_draw = np.sqrt((multiplicities @ group_sse_soap) / denominators)
        rmse_soed_draw = np.sqrt((multiplicities @ group_sse_soed) / denominators)
        deltas[start : start + count] = rmse_soap_draw - rmse_soed_draw
    tail = (1.0 - BOOTSTRAP_CONFIDENCE) / 2.0
    lower, upper = np.quantile(deltas, (tail, 1.0 - tail))
    p_bootstrap = float((1 + np.sum(deltas <= 0.0)) / (BOOTSTRAP_REPEATS + 1))
    absolute_soed = paired["absolute_error_soed"].to_numpy(dtype=float)
    absolute_soap = paired["absolute_error_soap"].to_numpy(dtype=float)
    try:
        test = wilcoxon(
            grouped["absolute_error_soed"],
            grouped["absolute_error_soap"],
            alternative="less",
            zero_method="zsplit",
        )
        statistic, wilcoxon_p = float(test.statistic), float(test.pvalue)
    except ValueError:
        statistic, wilcoxon_p = np.nan, np.nan
    relative_gain = (rmse_soap - rmse_soed) / rmse_soap
    paired["absolute_error_soed"] = absolute_soed
    paired["absolute_error_soap"] = absolute_soap
    paired["squared_error_difference_soap_minus_soed"] = (y_soap - y_true) ** 2 - (
        y_soed - y_true
    ) ** 2
    paired.insert(0, "comparison", comparison_label)
    paired.insert(1, "analysis_scope", analysis_scope)
    paired.insert(2, "soed_representation", proposed_name)
    distribution = pd.DataFrame(
        {
            "comparison": comparison_label,
            "analysis_scope": analysis_scope,
            "soed_representation": proposed_name,
            "bootstrap_index": np.arange(1, BOOTSTRAP_REPEATS + 1),
            "rmse_gain_soap_minus_soed_ev": deltas,
            "bootstrap_unit": "reduced_formula",
            "n_groups": n_groups,
        }
    )
    summary = pd.DataFrame(
        [
            {
                "comparison": comparison_label,
                "analysis_scope": analysis_scope,
                "soed_representation": proposed_name,
                "soap_representation": baseline_name,
                "n_test": len(y_true),
                "n_test_groups": n_groups,
                "bootstrap_unit": "reduced_formula",
                "soap_rmse_ev": rmse_soap,
                "soed_rmse_ev": rmse_soed,
                "rmse_gain_ev": rmse_soap - rmse_soed,
                "relative_rmse_gain": relative_gain,
                "bootstrap_ci_lower_ev": lower,
                "bootstrap_ci_upper_ev": upper,
                "bootstrap_one_sided_p": p_bootstrap,
                "wilcoxon_statistic": statistic,
                "wilcoxon_one_sided_p": wilcoxon_p,
            }
        ]
    )
    return paired, distribution, summary


def stratified_errors(prediction_frames: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for method, prediction_frame in prediction_frames.items():
        data = prediction_frame.loc[prediction_frame["split"] == "test"].copy()
        data["gap_bin"] = pd.cut(
            data["y_true"],
            bins=TAIL_METRIC_BINS_EV,
            labels=TAIL_METRIC_LABELS,
            include_lowest=True,
        )
        for gap_bin, subset in data.groupby("gap_bin", observed=False):
            if subset.empty:
                continue
            score = regression_metrics(subset["y_true"], subset["y_pred"])
            rows.append(
                {
                    "method": method,
                    "gap_bin": str(gap_bin),
                    "n_samples": len(subset),
                    "bias_ev": float(subset["residual"].mean()),
                    "mae_ev": score["mae"],
                    "rmse_ev": score["rmse"],
                }
            )
    return pd.DataFrame(rows)


def subset_regression_metrics(
    prediction_frames: dict[str, pd.DataFrame],
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    for representation in MAIN_REPRESENTATIONS:
        if representation not in prediction_frames:
            continue
        data = prediction_frames[representation]
        for split in ("train", "validation", "test"):
            split_data = data.loc[data["split"] == split]
            subsets = {
                "all": np.ones(len(split_data), dtype=bool),
                "positive_gap_gt_0p01_ev": (
                    split_data["y_true"].to_numpy(dtype=float)
                    > BAND_GAP_ZERO_THRESHOLD_EV
                ),
                "high_gap_ge_3_ev": (
                    split_data["y_true"].to_numpy(dtype=float) >= TAIL_THRESHOLD_EV
                ),
            }
            for subset_name, mask in subsets.items():
                subset = split_data.loc[mask]
                if subset.empty:
                    continue
                score = regression_metrics(subset["y_true"], subset["y_pred"])
                rows.append(
                    {
                        "representation": representation,
                        "method": representation_label(representation),
                        "split": split,
                        "subset": subset_name,
                        "n_samples": len(subset),
                        "mae_ev": score["mae"],
                        "rmse_ev": score["rmse"],
                        "r2": score["r2"],
                        "bias_ev": float(subset["residual"].mean()),
                    }
                )
    return pd.DataFrame(rows)


def save_physical_interpretation(
    direct_results: dict[str, dict[str, Any]],
    candidate_weights: dict[str, float],
    output: Path,
    target: str,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    base_name = MATCHED_CHEMICAL_SOED_NAME
    base_validation = result_validation_rmse(direct_results[base_name])
    base_test = float(
        direct_results[base_name]["metrics"]
        .loc[direct_results[base_name]["metrics"]["split"] == "test", "rmse"]
        .iloc[0]
    )
    mechanism_map = {
        "psoed_density_local_chemical": (
            "local density and chemical contrast",
            "Separates similar geometric environments by electron accumulation, nuclear-charge mismatch, and local element-property contrast.",
        ),
        "psoed_density_local_chemical_global": (
            "global and reciprocal density",
            "Adds long-range density modulation and cell-scale electronic inhomogeneity that local pooling can average away.",
        ),
        "psoed_density_local_chemical_radial": (
            "radial charge and charge transfer",
            "Resolves shell-wise electron redistribution around atoms and electron-minus-nuclear imbalance associated with bonding and ionicity.",
        ),
        "psoed_density_local_chemical_radial_global": (
            "radial plus global electronic structure",
            "Combines short-range charge transfer with long-range density organization while retaining separate prediction models.",
        ),
    }
    rows = []
    for name in SOED_ENHANCED_CANDIDATES:
        result = direct_results[name]
        validation_rmse = result_validation_rmse(result)
        test_rmse = float(
            result["metrics"].loc[result["metrics"]["split"] == "test", "rmse"].iloc[0]
        )
        block, hypothesis = mechanism_map[name]
        rows.append(
            {
                "representation": name,
                "electronic_block": block,
                "physical_hypothesis": hypothesis,
                "ensemble_weight": float(candidate_weights.get(name, 0.0)),
                "validation_rmse_ev": validation_rmse,
                "test_rmse_ev": test_rmse,
                "validation_rmse_change_vs_local_density_ev": validation_rmse
                - base_validation,
                "test_rmse_change_vs_local_density_ev": test_rmse - base_test,
                "selection_role": "validation-controlled candidate",
            }
        )
    mechanism = pd.DataFrame(rows)

    soap = direct_results[PRIMARY_LOCAL_CHEMICAL_SOAP_NAME]["predictions"]
    soed = direct_results[PRIMARY_SOED_NAME]["predictions"]
    paired = soap[["material_id", "split", "y_true", "y_pred"]].merge(
        soed[["material_id", "split", "y_true", "y_pred"]],
        on=["material_id", "split"],
        suffixes=("_soap", "_soed"),
        validate="one_to_one",
    )
    paired = paired.loc[
        (paired["split"] == "test")
        & (paired["y_true_soap"] >= MIDGAP_LOWER_EV)
        & (paired["y_true_soap"] < MIDGAP_UPPER_EV)
    ].copy()
    paired["absolute_error_soap_ev"] = np.abs(
        paired["y_pred_soap"] - paired["y_true_soap"]
    )
    paired["absolute_error_soed_ev"] = np.abs(
        paired["y_pred_soed"] - paired["y_true_soed"]
    )
    paired["absolute_error_change_soap_minus_soed_ev"] = (
        paired["absolute_error_soap_ev"] - paired["absolute_error_soed_ev"]
    )
    paired["soed_catastrophic_error"] = (
        paired["absolute_error_soed_ev"] >= CATASTROPHIC_ERROR_THRESHOLD_EV
    )
    paired = paired.sort_values("absolute_error_soed_ev", ascending=False)

    local_dimension = int(2 * (SOED_L_MAX + 1) * SOED_N_MAX * (SOED_N_MAX + 1) // 2)
    chemical_dimension = len(SOED_CHEMICAL_PROPERTIES) * SOED_CHEMICAL_PROJECTION_SIZE
    orbital_dimension = 4 * SOED_BLOCK_PROJECTION_SIZE
    composition_dimension = PeriodicSOEDFamily.composition_feature_count

    def block_ranges(name: str) -> list[tuple[str, int, int]]:
        start = 0
        ranges = [("density invariants", start, start + local_dimension)]
        start += local_dimension
        ranges.append(("local chemical contrast", start, start + chemical_dimension))
        start += chemical_dimension
        ranges.append(("orbital s/p/d/f contrast", start, start + orbital_dimension))
        start += orbital_dimension
        ranges.append(("composition", start, start + composition_dimension))
        start += composition_dimension
        if name in SOED_RADIAL_ELECTRONIC_DESCRIPTORS:
            ranges.append(
                (
                    "radial charge transfer",
                    start,
                    start + PeriodicSOEDFamily.electronic_radial_feature_count,
                )
            )
            start += PeriodicSOEDFamily.electronic_radial_feature_count
        if name in SOED_COMPACT_GLOBAL_DESCRIPTORS:
            ranges.append(
                (
                    "global/reciprocal density",
                    start,
                    start + PeriodicSOEDFamily.compact_global_feature_count,
                )
            )
        return ranges

    importance_rows: list[dict[str, Any]] = []
    for candidate_name in SOED_ENHANCED_CANDIDATES:
        candidate_result = direct_results[candidate_name]
        candidate_weight = float(candidate_weights.get(candidate_name, 0.0))
        expert_results = candidate_result.get("expert_results", {})
        expert_weights = candidate_result.get("expert_weights", {})
        for expert_name, expert_result in expert_results.items():
            model = expert_result.get("model_object")
            if model is None:
                continue
            raw = model.get_booster().get_score(importance_type="gain")
            gain = {
                int(key[1:]): float(value)
                for key, value in raw.items()
                if key.startswith("f") and key[1:].isdigit()
            }
            total_gain = float(sum(gain.values()))
            if total_gain <= 0.0:
                continue
            expert_weight = float(expert_weights.get(expert_name, 0.0))
            for block_name, lower, upper in block_ranges(candidate_name):
                block_gain = (
                    sum(
                        value for index, value in gain.items() if lower <= index < upper
                    )
                    / total_gain
                )
                importance_rows.append(
                    {
                        "candidate": candidate_name,
                        "expert": expert_name,
                        "feature_block": block_name,
                        "candidate_weight": candidate_weight,
                        "expert_weight": expert_weight,
                        "within_model_gain_fraction": block_gain,
                        "ensemble_weighted_gain": (
                            candidate_weight * expert_weight * block_gain
                        ),
                    }
                )
    importance_detail = pd.DataFrame(importance_rows)
    if importance_detail.empty:
        block_importance = pd.DataFrame(
            columns=["feature_block", "ensemble_weighted_gain", "normalized_importance"]
        )
    else:
        block_importance = (
            importance_detail.groupby("feature_block", as_index=False)[
                "ensemble_weighted_gain"
            ]
            .sum()
            .sort_values("ensemble_weighted_gain", ascending=False)
        )
        total = float(block_importance["ensemble_weighted_gain"].sum())
        block_importance["normalized_importance"] = (
            block_importance["ensemble_weighted_gain"] / total if total > 0.0 else 0.0
        )

    interpretation_dir = ensure_dir(output / "physical_interpretation")
    mechanism.to_csv(
        interpretation_dir / f"{target}_electronic_block_evidence.csv", index=False
    )
    save_dat(mechanism, interpretation_dir / f"{target}_electronic_block_evidence.dat")
    paired.to_csv(
        interpretation_dir / f"{target}_midgap_case_analysis.csv", index=False
    )
    save_dat(paired, interpretation_dir / f"{target}_midgap_case_analysis.dat")
    importance_detail.to_csv(
        interpretation_dir / f"{target}_feature_block_gain_detail.csv", index=False
    )
    save_dat(
        importance_detail,
        interpretation_dir / f"{target}_feature_block_gain_detail.dat",
    )
    block_importance.to_csv(
        interpretation_dir / f"{target}_feature_block_importance.csv", index=False
    )
    save_dat(
        block_importance,
        interpretation_dir / f"{target}_feature_block_importance.dat",
    )
    return mechanism, paired, block_importance


def plot_main_figures(
    frame: pd.DataFrame,
    target: str,
    metrics_frame: pd.DataFrame,
    predictions: dict[str, pd.DataFrame],
    histories: dict[str, pd.DataFrame],
    feature_summary: pd.DataFrame,
    paired: pd.DataFrame,
    bootstrap: pd.DataFrame,
    comparison: pd.DataFrame,
    robustness: pd.DataFrame,
    mechanism: pd.DataFrame,
    block_importance: pd.DataFrame,
    output: Path,
) -> pd.DataFrame:
    configure_matplotlib()
    import matplotlib.pyplot as plt

    main_dir = ensure_dir(output / "main_figures")
    data = frame[["material_id", "split", target, "is_nonmetal"]].dropna().copy()
    save_dat(data, main_dir / f"Figure02_{target}_dataset.dat")
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.5))
    splits = ("train", "validation", "test")
    zero = [
        int(((data["split"] == split) & (data["is_nonmetal"] == 0)).sum())
        for split in splits
    ]
    positive = [
        int(((data["split"] == split) & (data["is_nonmetal"] == 1)).sum())
        for split in splits
    ]
    x_value = np.arange(3)
    axes[0].bar(x_value, zero, color="tab:gray", label=r"$E_g\leq0.01$ eV")
    axes[0].bar(
        x_value, positive, bottom=zero, color="tab:blue", label=r"$E_g>0.01$ eV"
    )
    axes[0].set_xticks(x_value, [value.capitalize() for value in splits])
    axes[0].set_ylabel("Number of materials", fontsize=AXIS_LABEL_FONT_SIZE)
    ratio_label = ":".join(
        f"{100 * value:g}" for value in (TRAIN_RATIO, VALID_RATIO, TEST_RATIO)
    )
    axes[0].set_title(
        f"(a) Formula-grouped {ratio_label} split", fontsize=TITLE_FONT_SIZE
    )
    axes[0].legend(fontsize=LEGEND_FONT_SIZE, loc="upper right")
    bins = np.linspace(
        BAND_GAP_ZERO_THRESHOLD_EV, max(8.0, float(data[target].max())), 50
    )
    for split in splits:
        values = data.loc[
            (data["split"] == split) & (data[target] > BAND_GAP_ZERO_THRESHOLD_EV),
            target,
        ]
        axes[1].hist(
            values,
            bins=bins,
            density=True,
            histtype="step",
            lw=LINE_WIDTH,
            color=TAB_COLORS[split],
            label=f"{split} (n={len(values)})",
        )
    axes[1].set_xlabel(r"DFT band gap, $E_g$ (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
    axes[1].set_ylabel("Probability density", fontsize=AXIS_LABEL_FONT_SIZE)
    axes[1].set_title("(b) Positive-gap distribution", fontsize=TITLE_FONT_SIZE)
    axes[1].legend(fontsize=LEGEND_FONT_SIZE, loc="upper right")
    for axis in axes:
        axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
    fig.tight_layout()
    save_jpg(fig, main_dir / f"Figure02_{target}_dataset.jpg")
    plt.close(fig)

    selected_all = validation_selected_rows(metrics_frame, "direct_regression")
    selected = selected_all.loc[
        selected_all["representation"].isin(MAIN_REPRESENTATIONS)
    ].copy()
    selected["plot_order"] = selected["representation"].map(
        {name: index for index, name in enumerate(MAIN_REPRESENTATIONS)}
    )
    selected = selected.sort_values("plot_order").drop(columns="plot_order")
    save_dat(selected, main_dir / f"Figure03_{target}_descriptor_benchmark.dat")
    from matplotlib.patches import Patch

    fig, axes = plt.subplots(
        1, 2, figsize=(11.6, 4.8), gridspec_kw={"width_ratios": (1.15, 1.0)}
    )
    benchmark_short_labels = {
        COMPOSITION_ONLY_NAME: "Composition",
        PRIMARY_SOAP_NAME: "Density SOAP",
        MATCHED_SOED_NAME: r"Matched SOED ($\rho$)",
        PRIMARY_LOCAL_CHEMICAL_SOAP_NAME: "Fair SOAP",
        PRIMARY_SOED_NAME: "SOED ensemble",
    }
    labels = [benchmark_short_labels[value] for value in selected["representation"]]
    colors = [representation_color(value) for value in selected["representation"]]
    order = np.arange(len(selected))
    axes[0].barh(
        order - 0.18,
        selected["validation_rmse"],
        height=0.34,
        color="white",
        edgecolor=colors,
        hatch="////",
        linewidth=1.2,
    )
    axes[0].barh(
        order + 0.18, selected["test_rmse"], height=0.34, color=colors, edgecolor=colors
    )
    axes[0].set_yticks(order, labels)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("RMSE (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
    axes[0].set_title("(a) Grouped validation and test", fontsize=TITLE_FONT_SIZE)
    axes[0].legend(
        handles=(
            Patch(
                facecolor="white",
                edgecolor="tab:gray",
                hatch="////",
                label="Validation",
            ),
            Patch(facecolor="tab:gray", label="Test"),
        ),
        fontsize=LEGEND_FONT_SIZE,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.16),
        ncol=2,
    )
    for index, (_, row) in enumerate(selected.iterrows()):
        axes[1].scatter(
            row["test_rmse"],
            row["test_r2"],
            s=62,
            color=colors[index],
            label=benchmark_short_labels[row["representation"]],
            zorder=3,
        )
    axes[1].set_xlabel("Test RMSE (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
    axes[1].set_ylabel(r"Test $R^2$", fontsize=AXIS_LABEL_FONT_SIZE)
    axes[1].set_title(r"(b) Test RMSE vs. $R^2$", fontsize=TITLE_FONT_SIZE)
    axes[1].legend(
        fontsize=LEGEND_FONT_SIZE - 1, loc="center left", bbox_to_anchor=(1.02, 0.5)
    )
    for axis in axes:
        axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
    fig.subplots_adjust(left=0.17, right=0.78, bottom=0.19, top=0.86, wspace=0.38)
    save_jpg(fig, main_dir / f"Figure03_{target}_descriptor_benchmark.jpg")
    plt.close(fig)

    proposed_test = predictions[PRIMARY_SOED_NAME].loc[
        predictions[PRIMARY_SOED_NAME]["split"] == "test"
    ]
    stratified = stratified_errors(
        {
            "Local-chemical SOAP": predictions[PRIMARY_LOCAL_CHEMICAL_SOAP_NAME],
            "SOED ensemble": predictions[PRIMARY_SOED_NAME],
        }
    )
    save_dat(proposed_test, main_dir / f"Figure04_{target}_proposed_soed_test.dat")
    save_dat(stratified, main_dir / f"Figure04_{target}_stratified_errors.dat")
    score = regression_metrics(proposed_test["y_true"], proposed_test["y_pred"])
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 3.7))
    upper = 1.04 * float(
        max(proposed_test["y_true"].max(), proposed_test["y_pred"].max())
    )
    axes[0].scatter(
        proposed_test["y_true"],
        proposed_test["y_pred"],
        s=10,
        alpha=0.20,
        color=representation_color(PRIMARY_SOED_NAME),
        edgecolors="none",
        rasterized=True,
    )
    axes[0].plot([0, upper], [0, upper], "--", color="tab:gray", lw=1.2)
    axes[0].set_xlim(0, upper)
    axes[0].set_ylim(0, upper)
    axes[0].set_aspect("equal", adjustable="box")
    axes[0].set_xlabel(r"DFT $E_g$ (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
    axes[0].set_ylabel(r"Predicted $E_g$ (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
    axes[0].set_title("(a) SOED ensemble parity", fontsize=TITLE_FONT_SIZE)
    axes[0].text(
        0.04,
        0.96,
        f"RMSE={score['rmse']:.3f} eV\nMAE={score['mae']:.3f} eV\n$R^2$={score['r2']:.3f}",
        transform=axes[0].transAxes,
        va="top",
        fontsize=ANNOTATION_FONT_SIZE,
    )
    method_colors = {
        "Local-chemical SOAP": representation_color(PRIMARY_LOCAL_CHEMICAL_SOAP_NAME),
        "SOED ensemble": representation_color(PRIMARY_SOED_NAME),
    }
    for method, subset in stratified.groupby("method", sort=False):
        axes[1].plot(
            np.arange(len(subset)),
            subset["rmse_ev"],
            "o-",
            lw=LINE_WIDTH,
            color=method_colors[method],
            label=method,
        )
    first = stratified.loc[stratified["method"] == "Local-chemical SOAP"]
    gap_labels = [str(row.gap_bin) for row in first.itertuples(index=False)]
    axes[1].set_xticks(np.arange(len(first)), gap_labels, rotation=22, ha="right")
    axes[1].set_xlabel("DFT band-gap interval (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
    axes[1].set_ylabel("Test RMSE (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
    axes[1].set_title("(b) Gap-resolved RMSE", fontsize=TITLE_FONT_SIZE)
    upper_axis = axes[1].get_ylim()[1]
    for index, row in enumerate(first.itertuples(index=False)):
        axes[1].text(
            index,
            upper_axis * 0.98,
            f"n={int(row.n_samples)}",
            rotation=90,
            ha="center",
            va="top",
            fontsize=ANNOTATION_FONT_SIZE - 1,
            color="tab:gray",
        )
    handles, labels = axes[1].get_legend_handles_labels()
    for axis in axes:
        axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
    fig.legend(
        handles,
        labels,
        fontsize=LEGEND_FONT_SIZE,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.01),
        ncol=2,
    )
    fig.tight_layout(rect=(0.0, 0.11, 1.0, 1.0))
    save_jpg(fig, main_dir / f"Figure04_{target}_proposed_soed.jpg")
    plt.close(fig)

    save_dat(paired, main_dir / f"Figure05_{target}_paired_errors.dat")
    save_dat(bootstrap, main_dir / f"Figure05_{target}_bootstrap.dat")
    save_dat(comparison, main_dir / f"Figure05_{target}_comparison.dat")
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 3.9))
    panels = (
        (
            "Density-only match",
            MATCHED_SOED_NAME,
            "all_test",
            representation_color(MATCHED_SOED_NAME),
        ),
        (
            "Selected electronic SOED",
            PRIMARY_SOED_NAME,
            "all_test",
            representation_color(PRIMARY_SOED_NAME),
        ),
        (
            r"Selected SOED, $E_g\geq3$ eV",
            PRIMARY_SOED_NAME,
            "high_gap_ge_3_ev",
            "tab:red",
        ),
    )
    for panel_index, (panel_label, representation, scope, color) in enumerate(panels):
        axis = axes[panel_index]
        samples = bootstrap.loc[
            (bootstrap["soed_representation"] == representation)
            & (bootstrap["analysis_scope"] == scope),
            "rmse_gain_soap_minus_soed_ev",
        ]
        stats = comparison.loc[
            (comparison["soed_representation"] == representation)
            & (comparison["analysis_scope"] == scope)
        ].iloc[0]
        axis.hist(samples, bins=55, density=True, color=color, alpha=0.82)
        axis.axvline(0.0, color="tab:gray", linestyle="--", lw=1.2)
        axis.axvline(
            stats["rmse_gain_ev"], color="tab:blue", lw=LINE_WIDTH, label="Observed"
        )
        axis.axvspan(
            stats["bootstrap_ci_lower_ev"],
            stats["bootstrap_ci_upper_ev"],
            color="tab:orange",
            alpha=0.22,
            label="95% CI",
        )
        axis.set_xlabel("Baseline − proposed RMSE (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
        axis.set_title(
            f"({chr(97 + panel_index)}) {panel_label}", fontsize=TITLE_FONT_SIZE
        )
        axis.text(
            0.04,
            0.96,
            f"gain={stats['rmse_gain_ev']:.3f} eV\nrelative={100*stats['relative_rmse_gain']:.1f}%\n95% CI [{stats['bootstrap_ci_lower_ev']:.3f}, {stats['bootstrap_ci_upper_ev']:.3f}]",
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=ANNOTATION_FONT_SIZE,
        )
        axis.legend(fontsize=LEGEND_FONT_SIZE, loc="upper right")
    axes[0].set_ylabel("Group-bootstrap density", fontsize=AXIS_LABEL_FONT_SIZE)
    for axis in axes:
        axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
    fig.tight_layout()
    save_jpg(fig, main_dir / f"Figure05_{target}_statistical_comparison.jpg")
    plt.close(fig)

    candidate_order = (PRIMARY_LOCAL_CHEMICAL_SOAP_NAME, *SOED_ENHANCED_CANDIDATES)
    ablation = selected_all.loc[
        selected_all["representation"].isin(candidate_order)
    ].copy()
    ablation["plot_order"] = ablation["representation"].map(
        {name: index for index, name in enumerate(candidate_order)}
    )
    ablation = ablation.sort_values("plot_order").drop(columns="plot_order")
    ablation = ablation.merge(
        feature_summary[["representation", "feature_dimension"]].drop_duplicates(
            "representation"
        ),
        on="representation",
        how="left",
        suffixes=("", "_summary"),
    )
    ensemble_weight_map = dict(
        zip(mechanism["representation"], mechanism["ensemble_weight"])
    )
    ablation["ensemble_weight"] = (
        ablation["representation"].map(ensemble_weight_map).fillna(0.0)
    )
    ablation["validation_selected"] = ablation["ensemble_weight"].eq(
        ablation["ensemble_weight"].max()
    ) & ablation["representation"].isin(SOED_ENHANCED_CANDIDATES)
    soap_rmse = float(
        ablation.loc[
            ablation["representation"] == PRIMARY_LOCAL_CHEMICAL_SOAP_NAME,
            "test_rmse",
        ].iloc[0]
    )
    ablation["relative_gain_percent"] = (
        (soap_rmse - ablation["test_rmse"]) / soap_rmse * 100.0
    )
    save_dat(ablation, main_dir / f"Figure06_{target}_electronic_ablation.dat")
    fig, axes = plt.subplots(1, 3, figsize=(14.8, 4.6))
    x_value = np.arange(len(ablation))
    width = 0.36
    colors = [representation_color(value) for value in ablation["representation"]]
    axes[0].bar(
        x_value - width / 2,
        ablation["validation_rmse"],
        width=width,
        facecolor="white",
        edgecolor=colors,
        hatch="////",
        linewidth=1.2,
        label="Validation",
    )
    axes[0].bar(
        x_value + width / 2,
        ablation["test_rmse"],
        width=width,
        color=colors,
        label="Test",
    )
    short_labels = {
        PRIMARY_LOCAL_CHEMICAL_SOAP_NAME: "Fair SOAP",
        "psoed_density_local_chemical": r"SOED $\rho$",
        "psoed_density_local_chemical_global": "+ global",
        "psoed_density_local_chemical_radial": "+ radial charge",
        "psoed_density_local_chemical_radial_global": "+ radial + global",
    }
    axes[0].set_xticks(
        x_value,
        [short_labels[value] for value in ablation["representation"]],
        rotation=18,
        ha="right",
    )
    axes[0].set_ylabel("RMSE (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
    axes[0].set_title("(a) Electronic-feature ablation", fontsize=TITLE_FONT_SIZE)
    axes[0].legend(fontsize=LEGEND_FONT_SIZE, frameon=False)
    ensemble_rows = mechanism.set_index("representation").reindex(
        SOED_ENHANCED_CANDIDATES
    )
    ensemble_colors = [representation_color(name) for name in SOED_ENHANCED_CANDIDATES]
    ensemble_labels = [
        r"Local $\rho$",
        "+ global",
        "+ radial",
        "+ radial + global",
    ]
    axes[1].bar(
        np.arange(len(ensemble_rows)),
        ensemble_rows["ensemble_weight"],
        color=ensemble_colors,
    )
    axes[1].set_xticks(
        np.arange(len(ensemble_rows)), ensemble_labels, rotation=18, ha="right"
    )
    axes[1].set_ylim(
        0.0, max(1.0, 1.12 * float(ensemble_rows["ensemble_weight"].max()))
    )
    axes[1].set_ylabel(
        "Validation-derived ensemble weight", fontsize=AXIS_LABEL_FONT_SIZE
    )
    axes[1].set_title("(b) SOED-only candidate ensemble", fontsize=TITLE_FONT_SIZE)
    for index, value in enumerate(ensemble_rows["ensemble_weight"]):
        axes[1].text(
            index,
            float(value) + 0.02,
            f"{float(value):.2f}",
            ha="center",
            va="bottom",
            fontsize=ANNOTATION_FONT_SIZE,
        )
    if not block_importance.empty:
        block_plot = block_importance.sort_values(
            "normalized_importance", ascending=True
        )
        block_colors = [
            TAB_PALETTE[index % len(TAB_PALETTE)] for index in range(len(block_plot))
        ]
        axes[2].barh(
            np.arange(len(block_plot)),
            block_plot["normalized_importance"],
            color=block_colors,
        )
        axes[2].set_yticks(np.arange(len(block_plot)), block_plot["feature_block"])
        axes[2].set_xlabel(
            "Ensemble-weighted gain fraction", fontsize=AXIS_LABEL_FONT_SIZE
        )
    axes[2].set_title("(c) XGBoost block importance", fontsize=TITLE_FONT_SIZE)
    for axis in axes:
        axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
    fig.subplots_adjust(left=0.07, right=0.98, bottom=0.22, top=0.87, wspace=0.40)
    save_jpg(fig, main_dir / f"Figure06_{target}_efficiency_ablation.jpg")
    plt.close(fig)
    return stratified


def plot_si_figures(
    frame: pd.DataFrame,
    target: str,
    metrics_frame: pd.DataFrame,
    predictions: dict[str, pd.DataFrame],
    histories: dict[str, pd.DataFrame],
    classification_predictions: pd.DataFrame | None,
    gating_summary: pd.DataFrame,
    feature_summary: pd.DataFrame,
    stratified: pd.DataFrame,
    selected_soap_source: str,
    selected_local_chemical_soap_source: str,
    selected_soed_source: str,
    robustness: pd.DataFrame,
    bootstrap: pd.DataFrame,
    comparison: pd.DataFrame,
    mechanism: pd.DataFrame,
    midgap_cases: pd.DataFrame,
    block_importance: pd.DataFrame,
    output: Path,
) -> None:
    configure_matplotlib()
    import matplotlib.pyplot as plt

    si_dir = ensure_dir(output / "si_figures")
    distribution = frame[["material_id", "split", target]].dropna().copy()
    save_dat(distribution, si_dir / f"FigureS01_{target}_split_cdf.dat")
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    for split in ("train", "validation", "test"):
        values = np.sort(
            distribution.loc[distribution["split"] == split, target].to_numpy(
                dtype=float
            )
        )
        ax.plot(
            values,
            np.arange(1, len(values) + 1) / len(values),
            color=TAB_COLORS[split],
            lw=LINE_WIDTH,
            label=split,
        )
    ax.set_xlabel(r"DFT $E_g$ (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
    ax.set_ylabel("Empirical cumulative probability", fontsize=AXIS_LABEL_FONT_SIZE)
    ax.set_title(
        "Band-gap distributions after grouped splitting", fontsize=TITLE_FONT_SIZE
    )
    ax.legend(fontsize=LEGEND_FONT_SIZE, frameon=False)
    ax.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
    fig.tight_layout()
    save_jpg(fig, si_dir / f"FigureS01_{target}_split_cdf.jpg")
    plt.close(fig)

    proposed = predictions[PRIMARY_SOED_NAME]
    save_dat(proposed, si_dir / f"FigureS02_{target}_proposed_parity_splits.dat")
    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.5))
    for axis, split in zip(axes, ("train", "validation", "test")):
        subset = proposed.loc[proposed["split"] == split]
        score = regression_metrics(subset["y_true"], subset["y_pred"])
        upper = 1.04 * float(max(subset["y_true"].max(), subset["y_pred"].max()))
        axis.scatter(
            subset["y_true"],
            subset["y_pred"],
            s=8,
            alpha=0.25,
            edgecolors="none",
            color=TAB_COLORS[split],
        )
        axis.plot([0, upper], [0, upper], "--", color="tab:gray", lw=1.1)
        axis.set_xlim(0, upper)
        axis.set_ylim(0, upper)
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlabel(r"DFT $E_g$ (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
        axis.set_title(
            f"{split}: RMSE={score['rmse']:.3f}, $R^2$={score['r2']:.3f}",
            fontsize=TITLE_FONT_SIZE,
        )
        axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
    axes[0].set_ylabel(r"Predicted $E_g$ (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
    fig.tight_layout()
    save_jpg(fig, si_dir / f"FigureS02_{target}_proposed_parity_splits.jpg")
    plt.close(fig)
    for split in ("train", "validation", "test"):
        subset = proposed.loc[proposed["split"] == split].copy()
        save_dat(subset, si_dir / f"FigureS02_{target}_parity_{split}.dat")
        score = regression_metrics(subset["y_true"], subset["y_pred"])
        upper = 1.04 * float(max(subset["y_true"].max(), subset["y_pred"].max()))
        fig, ax = plt.subplots(figsize=(5.2, 5.0))
        ax.scatter(
            subset["y_true"],
            subset["y_pred"],
            s=MARKER_SIZE,
            alpha=0.30,
            edgecolors="none",
            color=TAB_COLORS[split],
            label=split.capitalize(),
        )
        ax.plot([0, upper], [0, upper], "--", color="tab:gray", lw=1.2)
        ax.set_xlim(0, upper)
        ax.set_ylim(0, upper)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlabel(r"DFT $E_g$ (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
        ax.set_ylabel(r"Predicted $E_g$ (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
        ax.set_title(f"{split.capitalize()} parity", fontsize=TITLE_FONT_SIZE)
        ax.text(
            0.04,
            0.96,
            f"RMSE={score['rmse']:.3f} eV\nMAE={score['mae']:.3f} eV\n$R^2$={score['r2']:.3f}",
            transform=ax.transAxes,
            va="top",
            fontsize=ANNOTATION_FONT_SIZE,
        )
        ax.legend(fontsize=LEGEND_FONT_SIZE, frameon=False)
        ax.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
        fig.tight_layout()
        save_jpg(fig, si_dir / f"FigureS02_{target}_parity_{split}.jpg")
        plt.close(fig)
    save_dat(proposed, si_dir / f"FigureS02_{target}_parity_combined.dat")
    upper = 1.04 * float(max(proposed["y_true"].max(), proposed["y_pred"].max()))
    fig, ax = plt.subplots(figsize=(5.4, 5.1))
    for split in ("train", "validation", "test"):
        subset = proposed.loc[proposed["split"] == split]
        ax.scatter(
            subset["y_true"],
            subset["y_pred"],
            s=MARKER_SIZE,
            alpha=0.25,
            edgecolors="none",
            color=TAB_COLORS[split],
            label=split.capitalize(),
        )
    ax.plot([0, upper], [0, upper], "--", color="tab:gray", lw=1.2)
    ax.set_xlim(0, upper)
    ax.set_ylim(0, upper)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"DFT $E_g$ (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
    ax.set_ylabel(r"Predicted $E_g$ (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
    ax.set_title("Combined parity", fontsize=TITLE_FONT_SIZE)
    ax.legend(fontsize=LEGEND_FONT_SIZE, frameon=False)
    ax.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
    fig.tight_layout()
    save_jpg(fig, si_dir / f"FigureS02_{target}_parity_combined.jpg")
    plt.close(fig)

    residual_data = pd.concat(
        [
            predictions[name].assign(method=representation_label(name))
            for name in (PRIMARY_LOCAL_CHEMICAL_SOAP_NAME, PRIMARY_SOED_NAME)
        ],
        ignore_index=True,
    )
    save_dat(residual_data, si_dir / f"FigureS03_{target}_residuals.dat")
    fig, axes = plt.subplots(1, 3, figsize=(14.0, 4.3))
    for axis, split in zip(axes, ("train", "validation", "test")):
        for name in (PRIMARY_LOCAL_CHEMICAL_SOAP_NAME, PRIMARY_SOED_NAME):
            color = representation_color(name)
            values = predictions[name].loc[
                predictions[name]["split"] == split, "residual"
            ]
            axis.hist(
                values,
                bins=70,
                density=True,
                histtype="step",
                lw=LINE_WIDTH,
                color=color,
                label=representation_label(name),
            )
        axis.axvline(0.0, color="tab:gray", linestyle="--", lw=1.0)
        axis.set_xlabel(r"Residual, $\hat E_g-E_g$ (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
        axis.set_title(split.capitalize(), fontsize=TITLE_FONT_SIZE)
        axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
    axes[0].set_ylabel("Probability density", fontsize=AXIS_LABEL_FONT_SIZE)
    axes[-1].legend(fontsize=LEGEND_FONT_SIZE, frameon=False)
    fig.tight_layout()
    save_jpg(fig, si_dir / f"FigureS03_{target}_residuals.jpg")
    plt.close(fig)

    learning_data = []
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5))
    for axis, name in zip(axes, (PRIMARY_LOCAL_CHEMICAL_SOAP_NAME, PRIMARY_SOED_NAME)):
        history = histories[name].copy()
        history["representation"] = name
        learning_data.append(history)
        axis.plot(
            history["iteration"],
            history["train_rmse"],
            color=TAB_COLORS["train"],
            lw=LINE_WIDTH,
            label="Train",
        )
        axis.plot(
            history["iteration"],
            history["valid_rmse"],
            color=TAB_COLORS["validation"],
            lw=LINE_WIDTH,
            label="Validation",
        )
        best_index = history["valid_rmse"].idxmin()
        axis.axvline(
            history.loc[best_index, "iteration"],
            color="tab:gray",
            linestyle="--",
            lw=1.1,
        )
        axis.set_xlabel("Boosting iteration", fontsize=AXIS_LABEL_FONT_SIZE)
        weighting = str(
            history.get("selected_weighting", pd.Series(["unweighted"])).iloc[0]
        )
        objective_label = (
            "Tail-weighted RMSE (eV)" if weighting == "tail_weighted" else "RMSE (eV)"
        )
        axis.set_ylabel(objective_label, fontsize=AXIS_LABEL_FONT_SIZE)
        axis.set_title(representation_label(name), fontsize=TITLE_FONT_SIZE)
        axis.legend(fontsize=LEGEND_FONT_SIZE, frameon=False)
        axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
    save_dat(
        pd.concat(learning_data, ignore_index=True),
        si_dir / f"FigureS04_{target}_learning_curves.dat",
    )
    fig.tight_layout()
    save_jpg(fig, si_dir / f"FigureS04_{target}_learning_curves.jpg")
    plt.close(fig)

    save_dat(gating_summary, si_dir / f"FigureS05_{target}_gating.dat")
    if not gating_summary.empty:
        fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))
        x_value = np.arange(len(gating_summary))
        colors = [
            TAB_PALETTE[index % len(TAB_PALETTE)]
            for index in range(len(gating_summary))
        ]
        for axis, metric, title in (
            (axes[0], "test_mae", "(a) Test MAE"),
            (axes[1], "test_rmse", "(b) Test RMSE"),
        ):
            axis.bar(x_value, gating_summary[metric], color=colors)
            axis.set_xticks(
                x_value, gating_summary["method_label"], rotation=18, ha="right"
            )
            axis.set_ylabel("Error (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
            axis.set_title(title, fontsize=TITLE_FONT_SIZE)
            axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
        fig.tight_layout()
        save_jpg(fig, si_dir / f"FigureS05_{target}_gating.jpg")
        plt.close(fig)

    if classification_predictions is not None:
        test = classification_predictions.loc[
            classification_predictions["split"] == "test"
        ]
        fpr, tpr, roc_thresholds = roc_curve(
            test["y_true_nonmetal"], test["p_nonmetal"]
        )
        precision, recall, pr_thresholds = precision_recall_curve(
            test["y_true_nonmetal"], test["p_nonmetal"]
        )
        observed, predicted = calibration_curve(
            test["y_true_nonmetal"], test["p_nonmetal"], n_bins=10, strategy="quantile"
        )
        save_dat(
            pd.DataFrame({"fpr": fpr, "tpr": tpr, "threshold": roc_thresholds}),
            si_dir / f"FigureS06_{target}_roc.dat",
        )
        save_dat(
            pd.DataFrame(
                {
                    "recall": recall,
                    "precision": precision,
                    "threshold": np.r_[pr_thresholds, np.nan],
                }
            ),
            si_dir / f"FigureS06_{target}_pr.dat",
        )
        save_dat(
            pd.DataFrame(
                {"predicted_probability": predicted, "observed_fraction": observed}
            ),
            si_dir / f"FigureS06_{target}_calibration.dat",
        )
        fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.3))
        axes[0].plot(fpr, tpr, color="tab:blue", lw=LINE_WIDTH)
        axes[0].plot([0, 1], [0, 1], "--", color="tab:gray")
        axes[0].set_xlabel("False-positive rate", fontsize=AXIS_LABEL_FONT_SIZE)
        axes[0].set_ylabel("True-positive rate", fontsize=AXIS_LABEL_FONT_SIZE)
        axes[0].set_title("(a) ROC", fontsize=TITLE_FONT_SIZE)
        axes[1].plot(recall, precision, color="tab:orange", lw=LINE_WIDTH)
        axes[1].set_xlabel("Recall", fontsize=AXIS_LABEL_FONT_SIZE)
        axes[1].set_ylabel("Precision", fontsize=AXIS_LABEL_FONT_SIZE)
        axes[1].set_title("(b) Precision–recall", fontsize=TITLE_FONT_SIZE)
        axes[2].plot(predicted, observed, "o-", color="tab:green", lw=LINE_WIDTH)
        axes[2].plot([0, 1], [0, 1], "--", color="tab:gray")
        axes[2].set_xlabel("Predicted probability", fontsize=AXIS_LABEL_FONT_SIZE)
        axes[2].set_ylabel("Observed fraction", fontsize=AXIS_LABEL_FONT_SIZE)
        axes[2].set_title("(c) Calibration", fontsize=TITLE_FONT_SIZE)
        for axis in axes:
            axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
        fig.tight_layout()
        save_jpg(fig, si_dir / f"FigureS06_{target}_classifier.jpg")
        plt.close(fig)

    weighted = validation_selected_rows(
        metrics_frame, "direct_regression_tail_weighted"
    )
    unweighted = validation_selected_rows(metrics_frame, "direct_regression_unweighted")
    tail_ablation = (
        weighted.merge(
            unweighted,
            on="representation",
            suffixes=("_weighted", "_unweighted"),
            how="inner",
        )
        if not unweighted.empty
        else pd.DataFrame()
    )
    if not tail_ablation.empty:
        tail_ablation = tail_ablation.loc[
            tail_ablation["representation"].isin(
                (selected_local_chemical_soap_source, selected_soed_source)
            )
        ].copy()
    save_dat(tail_ablation, si_dir / f"FigureS07_{target}_tail_weighting.dat")
    if not tail_ablation.empty:
        fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5))
        x_value = np.arange(len(tail_ablation))
        width = 0.36
        for axis, metric_name, title in (
            (axes[0], "test_rmse", "(a) Overall RMSE"),
            (axes[1], "test_tail_rmse", "(b) High-gap RMSE"),
        ):
            axis.bar(
                x_value - width / 2,
                tail_ablation[f"{metric_name}_unweighted"],
                width=width,
                color="tab:gray",
                label="Unweighted",
            )
            axis.bar(
                x_value + width / 2,
                tail_ablation[f"{metric_name}_weighted"],
                width=width,
                color="tab:red",
                label="Tail weighted",
            )
            axis.set_xticks(
                x_value,
                [
                    representation_label(value)
                    for value in tail_ablation["representation"]
                ],
                rotation=25,
                ha="right",
            )
            axis.set_ylabel("RMSE (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
            axis.set_title(title, fontsize=TITLE_FONT_SIZE)
            axis.legend(fontsize=LEGEND_FONT_SIZE, frameon=False)
            axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
        fig.tight_layout()
        save_jpg(fig, si_dir / f"FigureS07_{target}_tail_weighting.jpg")
        plt.close(fig)

    save_dat(feature_summary, si_dir / "FigureS08_descriptor_size_cost.dat")
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.3))
    labels = [
        representation_label(value) for value in feature_summary["representation"]
    ]
    descriptor_colors = [
        TAB_PALETTE[index % len(TAB_PALETTE)] for index in range(len(feature_summary))
    ]
    y_value = np.arange(len(labels))
    axes[0].barh(y_value, feature_summary["feature_dimension"], color=descriptor_colors)
    axes[0].set_xscale("log")
    axes[0].set_yticks(y_value, labels)
    axes[0].invert_yaxis()
    axes[0].set_xlabel("Feature dimension", fontsize=AXIS_LABEL_FONT_SIZE)
    axes[0].set_title("(a) Descriptor size", fontsize=TITLE_FONT_SIZE)
    axes[1].barh(
        y_value, feature_summary["seconds_per_structure"], color=descriptor_colors
    )
    axes[1].set_yticks(y_value, labels)
    axes[1].invert_yaxis()
    axes[1].set_xlabel(
        "Estimated standalone time (s structure$^{-1}$)", fontsize=AXIS_LABEL_FONT_SIZE
    )
    axes[1].set_title("(b) Descriptor cost (DFT excluded)", fontsize=TITLE_FONT_SIZE)
    for axis in axes:
        axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
    fig.tight_layout()
    save_jpg(fig, si_dir / "FigureS08_descriptor_size_cost.jpg")
    plt.close(fig)
    save_dat(stratified, si_dir / f"FigureS09_{target}_full_stratified_errors.dat")
    if not stratified.empty:
        methods = list(dict.fromkeys(stratified["method"].astype(str)))
        bins = list(dict.fromkeys(stratified["gap_bin"].astype(str)))
        method_colors = {
            "Local-chemical SOAP": representation_color(
                PRIMARY_LOCAL_CHEMICAL_SOAP_NAME
            ),
            "SOED ensemble": representation_color(PRIMARY_SOED_NAME),
        }
        fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5))
        x_value = np.arange(len(bins), dtype=float)
        width = 0.72 / max(len(methods), 1)
        for method_index, method in enumerate(methods):
            subset = (
                stratified.loc[stratified["method"] == method]
                .set_index("gap_bin")
                .reindex(bins)
            )
            offset = (method_index - (len(methods) - 1) / 2.0) * width
            axes[0].bar(
                x_value + offset,
                subset["mae_ev"],
                width=width,
                color=method_colors[method],
                label=method,
            )
            axes[1].bar(
                x_value + offset,
                subset["bias_ev"],
                width=width,
                color=method_colors[method],
                label=method,
            )
        for axis, ylabel, title in (
            (axes[0], "Test MAE (eV)", "(a) Absolute error"),
            (axes[1], "Mean signed error (eV)", "(b) Prediction bias"),
        ):
            axis.set_xticks(x_value, bins, rotation=28, ha="right")
            axis.set_xlabel("DFT band-gap interval (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
            axis.set_ylabel(ylabel, fontsize=AXIS_LABEL_FONT_SIZE)
            axis.set_title(title, fontsize=TITLE_FONT_SIZE)
            axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
        axes[1].axhline(0.0, color="tab:gray", linestyle="--", lw=1.0)
        handles, labels = axes[0].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.01),
            ncol=max(1, len(methods)),
            fontsize=LEGEND_FONT_SIZE,
        )
        fig.tight_layout(rect=(0.0, 0.12, 1.0, 1.0))
        save_jpg(fig, si_dir / f"FigureS09_{target}_full_stratified_errors.jpg")
        plt.close(fig)

    alpha_rows = validation_selected_rows(metrics_frame, "direct_regression")
    alpha_rows = alpha_rows.loc[
        alpha_rows["representation"].isin(
            (*SOAP_REPRESENTATIONS, *SOAP_LOCAL_CHEMICAL_REPRESENTATIONS)
        )
    ].copy()
    alpha_rows["soap_family"] = np.where(
        alpha_rows["representation"].isin(SOAP_LOCAL_CHEMICAL_REPRESENTATIONS),
        "Local-chemical SOAP",
        "Density-scaffold SOAP",
    )
    alpha_rows["soap_source"] = alpha_rows["representation"].map(
        lambda name: SOAP_LOCAL_CHEMICAL_SOURCE.get(name, name)
    )
    alpha_rows["soap_alpha"] = alpha_rows["soap_source"].map(SOAP_ALPHA_BY_NAME)
    alpha_rows = alpha_rows.sort_values(["soap_family", "soap_alpha"])
    save_dat(alpha_rows, si_dir / f"FigureS10_{target}_soap_alpha_sensitivity.dat")
    if not alpha_rows.empty:
        fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4), sharey=True)
        family_specs = (
            ("Density-scaffold SOAP", selected_soap_source, "tab:purple"),
            ("Local-chemical SOAP", selected_local_chemical_soap_source, "tab:pink"),
        )
        for axis, (family, selected_source, color) in zip(axes, family_specs):
            subset = alpha_rows.loc[alpha_rows["soap_family"] == family]
            axis.plot(
                subset["soap_alpha"],
                subset["validation_rmse"],
                "o-",
                color=color,
                lw=LINE_WIDTH,
                label="Validation",
            )
            axis.plot(
                subset["soap_alpha"],
                subset["test_rmse"],
                "s--",
                color="tab:blue",
                lw=LINE_WIDTH,
                label="Test (diagnostic)",
            )
            selected_row = subset.loc[subset["representation"] == selected_source].iloc[
                0
            ]
            axis.axvline(
                float(selected_row["soap_alpha"]),
                color="tab:gray",
                linestyle=":",
                lw=1.2,
                label="Validation-selected",
            )
            axis.set_xscale("log", base=2)
            axis.set_xticks(
                list(SOAP_ALPHA_GRID), [f"{value:g}" for value in SOAP_ALPHA_GRID]
            )
            axis.set_xlabel(
                r"Gaussian width, $\alpha$ ($\mathrm{\AA}^{-2}$)",
                fontsize=AXIS_LABEL_FONT_SIZE,
            )
            axis.set_title(family, fontsize=TITLE_FONT_SIZE)
            axis.legend(fontsize=LEGEND_FONT_SIZE - 1, loc="best")
            axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
        axes[0].set_ylabel("RMSE (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
        fig.tight_layout()
        save_jpg(fig, si_dir / f"FigureS10_{target}_soap_alpha_sensitivity.jpg")
        plt.close(fig)

    save_dat(robustness, si_dir / f"FigureS11_{target}_repeated_group_splits.dat")
    if not robustness.empty:
        fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.5))
        method_order = (
            PRIMARY_LOCAL_CHEMICAL_SOAP_NAME,
            PRIMARY_SOED_NAME,
        )
        colors = tuple(representation_color(name) for name in method_order)
        for method, color in zip(method_order, colors):
            subset = robustness.loc[robustness["representation"] == method].sort_values(
                "split_seed"
            )
            robust_label = (
                "Validation-selected local-chemical SOAP"
                if method == PRIMARY_LOCAL_CHEMICAL_SOAP_NAME
                else "Validation-selected electronic SOED"
            )
            axes[0].plot(
                subset["split_seed"],
                subset["test_rmse"],
                "o-",
                color=color,
                lw=LINE_WIDTH,
                label=robust_label,
            )
        pivot = robustness.pivot(
            index="split_seed", columns="representation", values="test_rmse"
        )
        for method in (PRIMARY_SOED_NAME,):
            color = representation_color(method)
            gain = (
                100.0
                * (pivot[PRIMARY_LOCAL_CHEMICAL_SOAP_NAME] - pivot[method])
                / pivot[PRIMARY_LOCAL_CHEMICAL_SOAP_NAME]
            )
            axes[1].plot(
                gain.index,
                gain.values,
                "o-",
                color=color,
                lw=LINE_WIDTH,
                label="Validation-selected electronic SOED",
            )
        axes[0].set_ylabel("Test RMSE (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
        axes[0].set_title("(a) Absolute performance", fontsize=TITLE_FONT_SIZE)
        axes[1].axhline(0.0, color="tab:gray", linestyle="--", lw=1.0)
        axes[1].set_ylabel(
            "RMSE improvement over local-chemical SOAP (%)",
            fontsize=AXIS_LABEL_FONT_SIZE,
        )
        axes[1].set_title("(b) Paired improvement", fontsize=TITLE_FONT_SIZE)
        for axis in axes:
            axis.set_xlabel("Formula-group split seed", fontsize=AXIS_LABEL_FONT_SIZE)
            axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
        handles_left, labels_left = axes[0].get_legend_handles_labels()
        handles_right, labels_right = axes[1].get_legend_handles_labels()
        legend_entries = dict(
            zip(labels_left + labels_right, handles_left + handles_right)
        )
        fig.legend(
            list(legend_entries.values()),
            list(legend_entries.keys()),
            fontsize=LEGEND_FONT_SIZE - 1,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.01),
            ncol=2,
        )
        fig.tight_layout(rect=(0.0, 0.15, 1.0, 1.0))
        save_jpg(fig, si_dir / f"FigureS11_{target}_repeated_group_splits.jpg")
        plt.close(fig)

    fairness_order = [COMPOSITION_ONLY_NAME, PRIMARY_LOCAL_CHEMICAL_SOAP_NAME]
    if selected_soed_source != MATCHED_CHEMICAL_SOED_NAME:
        fairness_order.append(MATCHED_CHEMICAL_SOED_NAME)
    fairness_order.append(PRIMARY_SOED_NAME)
    fairness_order = tuple(fairness_order)
    fairness = validation_selected_rows(metrics_frame, "direct_regression")
    fairness = fairness.loc[fairness["representation"].isin(fairness_order)].copy()
    fairness["plot_order"] = fairness["representation"].map(
        {name: index for index, name in enumerate(fairness_order)}
    )
    fairness = fairness.sort_values("plot_order").drop(columns="plot_order")
    fair_bootstrap = bootstrap.loc[
        (bootstrap["soed_representation"] == PRIMARY_SOED_NAME)
        & (bootstrap["analysis_scope"] == "all_test")
    ].copy()
    fair_comparison = comparison.loc[
        (comparison["soed_representation"] == PRIMARY_SOED_NAME)
        & (comparison["analysis_scope"] == "all_test")
    ].copy()
    save_dat(fairness, si_dir / f"FigureS12_{target}_chemistry_fairness_benchmark.dat")
    save_dat(
        fair_bootstrap, si_dir / f"FigureS12_{target}_chemistry_fairness_bootstrap.dat"
    )
    if not fairness.empty and not fair_comparison.empty:
        fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.5))
        x_value = np.arange(len(fairness))
        width = 0.36
        fair_colors = [
            representation_color(value) for value in fairness["representation"]
        ]
        axes[0].bar(
            x_value - width / 2,
            fairness["validation_rmse"],
            width=width,
            facecolor="white",
            edgecolor=fair_colors,
            hatch="////",
            linewidth=1.2,
            label="Validation",
        )
        axes[0].bar(
            x_value + width / 2,
            fairness["test_rmse"],
            width=width,
            color=fair_colors,
            label="Test",
        )
        axes[0].set_xticks(
            x_value,
            [representation_label(value) for value in fairness["representation"]],
            rotation=24,
            ha="right",
        )
        axes[0].set_ylabel("RMSE (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
        axes[0].set_title(
            "(a) Local-chemistry-controlled benchmark", fontsize=TITLE_FONT_SIZE
        )
        axes[0].legend(fontsize=LEGEND_FONT_SIZE, loc="upper right")
        stats = fair_comparison.iloc[0]
        samples = fair_bootstrap["rmse_gain_soap_minus_soed_ev"]
        axes[1].hist(
            samples,
            bins=55,
            density=True,
            color=representation_color(PRIMARY_SOED_NAME),
            alpha=0.82,
        )
        axes[1].axvline(0.0, color="tab:gray", linestyle="--", lw=1.2)
        axes[1].axvline(
            stats["rmse_gain_ev"], color="tab:blue", lw=LINE_WIDTH, label="Observed"
        )
        axes[1].axvspan(
            stats["bootstrap_ci_lower_ev"],
            stats["bootstrap_ci_upper_ev"],
            color="tab:green",
            alpha=0.22,
            label="95% CI",
        )
        axes[1].set_xlabel(
            "RMSE gain (local-chemical SOAP − SOED; eV)",
            fontsize=AXIS_LABEL_FONT_SIZE,
        )
        axes[1].set_ylabel("Group-bootstrap density", fontsize=AXIS_LABEL_FONT_SIZE)
        axes[1].set_title("(b) Fair paired inference", fontsize=TITLE_FONT_SIZE)
        axes[1].legend(fontsize=LEGEND_FONT_SIZE, loc="upper right")
        for axis in axes:
            axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
        fig.tight_layout()
        save_jpg(fig, si_dir / f"FigureS12_{target}_chemistry_fairness.jpg")
        plt.close(fig)

    gap_order = (
        "gap_le_0p01_ev",
        "gap_0p01_to_0p5_ev",
        "gap_0p5_to_1_ev",
        "gap_1_to_2_ev",
        "gap_2_to_3_ev",
        "gap_3_to_6_ev",
        "gap_gt_6_ev",
    )
    gap_labels = (
        r"$E_g\leq0.01$",
        "0.01–0.5",
        "0.5–1",
        "1–2",
        "2–3",
        "3–6",
        r"$E_g>6$",
    )
    gap_inference = comparison.loc[
        (comparison["soed_representation"] == PRIMARY_SOED_NAME)
        & comparison["analysis_scope"].isin(gap_order)
    ].copy()
    gap_inference["plot_order"] = gap_inference["analysis_scope"].map(
        {name: index for index, name in enumerate(gap_order)}
    )
    gap_inference = gap_inference.sort_values("plot_order").drop(columns="plot_order")
    save_dat(gap_inference, si_dir / f"FigureS13_{target}_gap_bin_inference.dat")
    if len(gap_inference) == len(gap_order):
        fig, ax = plt.subplots(figsize=(7.2, 4.8))
        y_value = np.arange(len(gap_inference))
        gain = gap_inference["rmse_gain_ev"].to_numpy(dtype=float)
        lower = gap_inference["bootstrap_ci_lower_ev"].to_numpy(dtype=float)
        upper = gap_inference["bootstrap_ci_upper_ev"].to_numpy(dtype=float)
        significant = gap_inference["bootstrap_p_fdr_bh"].to_numpy(dtype=float) < 0.05
        for index in range(len(gap_inference)):
            ax.errorbar(
                gain[index],
                y_value[index],
                xerr=np.asarray(
                    [[gain[index] - lower[index]], [upper[index] - gain[index]]]
                ),
                fmt="o",
                color="tab:blue" if significant[index] else "tab:gray",
                markerfacecolor="tab:blue" if significant[index] else "white",
                markeredgecolor="tab:blue" if significant[index] else "tab:gray",
                capsize=3,
                lw=1.3,
            )
            q_value = float(gap_inference.iloc[index]["bootstrap_p_fdr_bh"])
            ax.annotate(
                f"q={q_value:.3f}",
                (upper[index], y_value[index]),
                xytext=(5, 0),
                textcoords="offset points",
                va="center",
                fontsize=ANNOTATION_FONT_SIZE - 1,
                color="tab:blue" if significant[index] else "tab:gray",
            )
        ax.axvline(0.0, color="tab:red", linestyle="--", lw=1.0)
        ax.set_yticks(y_value, gap_labels)
        ax.invert_yaxis()
        ax.set_xlabel(
            "RMSE gain (local-chemical SOAP − selected SOED; eV)",
            fontsize=AXIS_LABEL_FONT_SIZE,
        )
        ax.set_ylabel("DFT band-gap interval (eV)", fontsize=AXIS_LABEL_FONT_SIZE)
        ax.set_title(
            "Gap-resolved grouped-bootstrap inference", fontsize=TITLE_FONT_SIZE
        )
        from matplotlib.lines import Line2D

        ax.legend(
            handles=(
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="none",
                    markerfacecolor="tab:blue",
                    markeredgecolor="tab:blue",
                    label=r"FDR $q<0.05$",
                ),
                Line2D(
                    [0],
                    [0],
                    marker="o",
                    linestyle="none",
                    markerfacecolor="white",
                    markeredgecolor="tab:gray",
                    label=r"FDR $q\geq0.05$",
                ),
            ),
            fontsize=LEGEND_FONT_SIZE - 1,
            loc="lower right",
        )
        ax.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
        fig.tight_layout()
        save_jpg(fig, si_dir / f"FigureS13_{target}_gap_bin_inference.jpg")
        plt.close(fig)

    save_dat(mechanism, si_dir / f"FigureS14_{target}_electronic_mechanisms.dat")
    if not mechanism.empty:
        fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))
        short = (r"Local $\rho$", "+ global", "+ radial", "+ radial + global")
        x_value = np.arange(len(mechanism))
        colors = [representation_color(value) for value in mechanism["representation"]]
        axes[0].bar(x_value, mechanism["ensemble_weight"], color=colors)
        axes[0].set_ylabel(
            "Validation-derived ensemble weight", fontsize=AXIS_LABEL_FONT_SIZE
        )
        axes[0].set_title("(a) Electronic-channel use", fontsize=TITLE_FONT_SIZE)
        axes[1].bar(
            x_value - 0.18,
            -mechanism["validation_rmse_change_vs_local_density_ev"],
            width=0.36,
            facecolor="white",
            edgecolor=colors,
            hatch="////",
            linewidth=1.2,
            label="Validation",
        )
        axes[1].bar(
            x_value + 0.18,
            -mechanism["test_rmse_change_vs_local_density_ev"],
            width=0.36,
            color=colors,
            label="Test (diagnostic)",
        )
        axes[1].axhline(0.0, color="tab:gray", linestyle="--", lw=1.0)
        axes[1].set_ylabel(
            r"RMSE gain vs. local $\rho$ SOED (eV)", fontsize=AXIS_LABEL_FONT_SIZE
        )
        axes[1].set_title("(b) Controlled block ablation", fontsize=TITLE_FONT_SIZE)
        axes[1].legend(fontsize=LEGEND_FONT_SIZE - 1, loc="best")
        for axis in axes:
            axis.set_xticks(x_value, short, rotation=18, ha="right")
            axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
        fig.subplots_adjust(left=0.09, right=0.98, bottom=0.23, top=0.87, wspace=0.30)
        save_jpg(fig, si_dir / f"FigureS14_{target}_electronic_mechanisms.jpg")
        plt.close(fig)

    save_dat(midgap_cases, si_dir / f"FigureS15_{target}_midgap_extremes.dat")
    if not midgap_cases.empty:
        fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.4))
        axes[0].scatter(
            midgap_cases["absolute_error_soap_ev"],
            midgap_cases["absolute_error_soed_ev"],
            s=18,
            alpha=0.55,
            color="tab:brown",
            edgecolors="none",
        )
        upper = 1.05 * float(
            max(
                midgap_cases["absolute_error_soap_ev"].max(),
                midgap_cases["absolute_error_soed_ev"].max(),
            )
        )
        axes[0].plot([0, upper], [0, upper], "--", color="tab:gray", lw=1.0)
        axes[0].set_xlim(0, upper)
        axes[0].set_ylim(0, upper)
        axes[0].set_xlabel(
            "Fair SOAP absolute error (eV)", fontsize=AXIS_LABEL_FONT_SIZE
        )
        axes[0].set_ylabel(
            "SOED ensemble absolute error (eV)", fontsize=AXIS_LABEL_FONT_SIZE
        )
        axes[0].set_title("(a) 1–2 eV paired errors", fontsize=TITLE_FONT_SIZE)
        changes = midgap_cases["absolute_error_change_soap_minus_soed_ev"]
        axes[1].hist(changes, bins=35, color="tab:orange", alpha=0.85)
        axes[1].axvline(0.0, color="tab:gray", linestyle="--", lw=1.0)
        axes[1].set_xlabel(
            "Absolute-error gain, SOAP − SOED (eV)", fontsize=AXIS_LABEL_FONT_SIZE
        )
        axes[1].set_ylabel("Number of test materials", fontsize=AXIS_LABEL_FONT_SIZE)
        axes[1].set_title("(b) Mid-gap gain distribution", fontsize=TITLE_FONT_SIZE)
        for axis in axes:
            axis.tick_params(labelsize=TICK_LABEL_FONT_SIZE)
        fig.tight_layout()
        save_jpg(fig, si_dir / f"FigureS15_{target}_midgap_extremes.jpg")
        plt.close(fig)


def save_tables(
    metrics_frame: pd.DataFrame,
    classification_frame: pd.DataFrame,
    gating_summary: pd.DataFrame,
    stratified: pd.DataFrame,
    feature_summary: pd.DataFrame,
    paired: pd.DataFrame,
    bootstrap: pd.DataFrame,
    comparison: pd.DataFrame,
    robustness: pd.DataFrame,
    robustness_summary: pd.DataFrame,
    subset_metrics: pd.DataFrame,
    mechanism: pd.DataFrame,
    midgap_cases: pd.DataFrame,
    block_importance: pd.DataFrame,
    target: str,
    output: Path,
) -> None:
    table_dir = ensure_dir(output / "tables")
    feature_info = feature_summary.rename(
        columns={"feature_dimension": "descriptor_dimension"}
    )
    main = validation_selected_rows(metrics_frame, "direct_regression").merge(
        feature_info, on="representation", how="left"
    )
    main = main.loc[main["representation"].isin(MAIN_REPRESENTATIONS)].copy()
    main["method"] = main["representation"].map(representation_label)
    main["table_order"] = main["representation"].map(
        {name: index for index, name in enumerate(MAIN_REPRESENTATIONS)}
    )
    main = main.sort_values("table_order").drop(columns="table_order")
    display_columns = [
        "method",
        "selected_weighting",
        "tail_weighted_fraction",
        "descriptor_dimension",
        "validation_rmse",
        "test_mae",
        "test_rmse",
        "test_r2",
        "test_tail_rmse",
        "seconds_per_structure",
    ]
    display = main[[column for column in display_columns if column in main]].copy()
    display.to_csv(table_dir / f"Table1_{target}_primary_benchmark.csv", index=False)
    save_dat(display, table_dir / f"Table1_{target}_primary_benchmark.dat")
    latex_methods = {
        representation_label(PRIMARY_SOAP_NAME): representation_label(
            PRIMARY_SOAP_NAME
        ),
        representation_label(
            MATCHED_SOED_NAME
        ): r"Matched periodic SOED ($\rho$, 210-D)",
        representation_label(
            PRIMARY_LOCAL_CHEMICAL_SOAP_NAME
        ): r"Local-chemical periodic SOAP",
        representation_label(
            MATCHED_CHEMICAL_SOED_NAME
        ): r"Local-chemical periodic SOED ($\rho$)",
        representation_label("psoed_log_density"): r"Periodic SOED ($\log\rho$)",
        representation_label(
            "psoed_density_gradient"
        ): r"Periodic SOED ($\rho,|\nabla\rho|$)",
        representation_label(
            "psoed_physics_multichannel"
        ): r"Physics-multichannel periodic SOED",
        representation_label(PRIMARY_SOED_NAME): r"Validation-selected electronic SOED",
    }
    latex_lines = [
        r"\begin{table*}[t]",
        r"\caption{Frozen-confirmatory periodic SOAP/SOED benchmark after undefined-Pauling-element filtering. Local-chemical SOAP and density SOED use identical element-property contrast, orbital-block, composition, and invariant pooling blocks. The electronic SOED family adds compact radial charge-transfer and global density features without SOAP--SOED fusion. SOAP width, regression-expert weights, and SOED candidate-ensemble weights were determined only from formula-disjoint validation data; the held-out test set was not used for selection. Dimension denotes the largest single ensemble member. The high-gap subset has $E_g\geq3$ eV.}",
        rf"\label{{tab:{target}_primary}}",
        r"\centering",
        r"\small",
        r"\setlength{\tabcolsep}{4.5pt}",
        r"\begin{tabular}{llrrrrrrr}",
        r"\toprule",
        r"Representation & Weighting & Dimension & Val. RMSE & Test MAE & Test RMSE & Test $R^2$ & High-gap RMSE & Time \\",
        r" & & & (eV) & (eV) & (eV) &  & (eV) & (s structure$^{-1}$) \\",
        r"\midrule",
    ]
    for row in display.itertuples(index=False):
        method = latex_methods.get(str(row.method), str(row.method).replace("_", r"\_"))
        latex_lines.append(
            f"{method} & {str(row.selected_weighting).replace('_', ' ')} & {int(row.descriptor_dimension):d} & {row.validation_rmse:.3f} & "
            f"{row.test_mae:.3f} & {row.test_rmse:.3f} & {row.test_r2:.3f} & "
            f"{row.test_tail_rmse:.3f} & {row.seconds_per_structure:.3f} " + r"\\"
        )
    latex_lines.extend(
        [
            r"\bottomrule",
            r"\end{tabular}",
            r"\end{table*}",
        ]
    )
    (table_dir / f"Table1_{target}_primary_benchmark.tex").write_text(
        "\n".join(latex_lines) + "\n", encoding="utf-8"
    )
    quality_dir = output / "data_quality"
    filter_summary_path = quality_dir / "undefined_pauling_filter_summary.csv"
    excluded_path = quality_dir / "excluded_undefined_pauling_elements.csv"
    alpha_sensitivity = validation_selected_rows(metrics_frame, "direct_regression")
    alpha_sensitivity = alpha_sensitivity.loc[
        alpha_sensitivity["representation"].isin(
            (*SOAP_REPRESENTATIONS, *SOAP_LOCAL_CHEMICAL_REPRESENTATIONS)
        )
    ].copy()
    alpha_sensitivity["soap_family"] = np.where(
        alpha_sensitivity["representation"].isin(SOAP_LOCAL_CHEMICAL_REPRESENTATIONS),
        "local_chemical",
        "density_scaffold",
    )
    alpha_sensitivity["soap_source"] = alpha_sensitivity["representation"].map(
        lambda name: SOAP_LOCAL_CHEMICAL_SOURCE.get(name, name)
    )
    alpha_sensitivity["soap_alpha"] = alpha_sensitivity["soap_source"].map(
        SOAP_ALPHA_BY_NAME
    )
    primary_row = main.loc[main["representation"] == PRIMARY_SOED_NAME]
    selected_soed_source = (
        str(primary_row["source_representation"].iloc[0])
        if not primary_row.empty and "source_representation" in primary_row
        else MATCHED_CHEMICAL_SOED_NAME
    )
    fairness_order = [COMPOSITION_ONLY_NAME, PRIMARY_LOCAL_CHEMICAL_SOAP_NAME]
    if selected_soed_source != MATCHED_CHEMICAL_SOED_NAME:
        fairness_order.append(MATCHED_CHEMICAL_SOED_NAME)
    fairness_order.append(PRIMARY_SOED_NAME)
    fairness_order = tuple(fairness_order)
    fairness = validation_selected_rows(metrics_frame, "direct_regression")
    fairness = fairness.loc[fairness["representation"].isin(fairness_order)].copy()
    fairness["table_order"] = fairness["representation"].map(
        {name: index for index, name in enumerate(fairness_order)}
    )
    fairness = fairness.sort_values("table_order").drop(columns="table_order")
    fairness_comparison = comparison.loc[
        (comparison["soed_representation"] == PRIMARY_SOED_NAME)
        & (comparison["analysis_scope"] == "all_test")
    ].copy()
    gap_inference = comparison.loc[
        (comparison["soed_representation"] == PRIMARY_SOED_NAME)
        & comparison["analysis_scope"].str.startswith("gap_")
    ].copy()
    exports = {
        "TableS01_undefined_pauling_filter_summary": (
            pd.read_csv(filter_summary_path)
            if filter_summary_path.exists()
            else pd.DataFrame()
        ),
        "TableS02_excluded_structures": (
            pd.read_csv(excluded_path) if excluded_path.exists() else pd.DataFrame()
        ),
        f"TableS03_{target}_all_regression_metrics": metrics_frame,
        f"TableS04_{target}_classification_metrics": classification_frame,
        f"TableS05_{target}_gating_summary": gating_summary,
        f"TableS06_{target}_stratified_errors": stratified,
        "TableS07_descriptor_summary": feature_summary,
        f"TableS08_{target}_paired_material_errors": paired,
        f"TableS09_{target}_bootstrap_distribution": bootstrap,
        f"TableS10_{target}_statistical_comparison": comparison,
        f"TableS11_{target}_soap_alpha_sensitivity": alpha_sensitivity,
        f"TableS12_{target}_repeated_group_splits": robustness,
        f"TableS13_{target}_repeated_group_split_summary": robustness_summary,
        f"TableS14_{target}_chemistry_fairness_benchmark": fairness,
        f"TableS15_{target}_chemistry_fairness_inference": fairness_comparison,
        f"TableS16_{target}_positive_and_high_gap_metrics": subset_metrics,
        f"TableS17_{target}_gap_bin_inference": gap_inference,
        f"TableS18_{target}_electronic_mechanism_evidence": mechanism,
        f"TableS19_{target}_midgap_case_analysis": midgap_cases,
        f"TableS20_{target}_feature_block_importance": block_importance,
    }
    for name, frame in exports.items():
        frame.to_csv(table_dir / f"{name}.csv", index=False)
        save_dat(frame, table_dir / f"{name}.dat")


def build_x_sets(
    representation: str,
    arrays: dict[str, np.ndarray],
    split_indices: dict[str, np.ndarray],
) -> dict[str, np.ndarray]:
    if representation in SOAP_COMPOSITION_SOURCE:
        components = (SOAP_COMPOSITION_SOURCE[representation], COMPOSITION_ONLY_NAME)
    elif representation == DENSITY_COMPOSITION_NAME:
        components = (MATCHED_SOED_NAME, COMPOSITION_ONLY_NAME)
    else:
        components = (representation,)
    return {
        split: np.ascontiguousarray(
            np.concatenate(
                [
                    np.asarray(arrays[name][indices], dtype=np.float32)
                    for name in components
                ],
                axis=1,
            )
            if len(components) > 1
            else np.asarray(arrays[components[0]][indices], dtype=np.float32)
        )
        for split, indices in split_indices.items()
    }


def split_labels_for_seed(frame: pd.DataFrame, seed: int) -> np.ndarray:
    indices = np.arange(len(frame))
    groups = frame["reduced_formula"].to_numpy(dtype=str)
    strata = frame["band_gap_stratum"].to_numpy(dtype=str)
    train_valid, test = best_group_holdout(indices, groups, strata, TEST_RATIO, seed)
    train, valid = best_group_holdout(
        train_valid,
        groups[train_valid],
        strata[train_valid],
        VALID_RATIO / (TRAIN_RATIO + VALID_RATIO),
        seed + 1,
    )
    labels = np.full(len(frame), "", dtype=object)
    labels[train] = "train"
    labels[valid] = "validation"
    labels[test] = "test"
    group_sets = {
        name: set(groups[labels == name]) for name in ("train", "validation", "test")
    }
    if any(
        group_sets[left] & group_sets[right]
        for left, right in (
            ("train", "validation"),
            ("train", "test"),
            ("validation", "test"),
        )
    ):
        raise RuntimeError(f"Reduced-formula leakage in robustness split seed={seed}")
    return labels.astype(str)


def fit_robustness_candidate(
    x_sets: dict[str, np.ndarray],
    y_sets: dict[str, np.ndarray],
    weighting_mode: str,
    hardware: dict[str, Any],
    logger: logging.Logger,
) -> dict[str, Any]:
    weights = {
        split: regression_sample_weights(values, weighting_mode)
        for split, values in y_sets.items()
    }
    use_gpu = bool(hardware.get("gpu_available"))
    try:
        model, _, _ = fit_xgb_regressor(
            x_sets["train"],
            y_sets["train"],
            x_sets["validation"],
            y_sets["validation"],
            dict(XGBOOST_PARAMS),
            use_gpu,
            weights["train"],
            weights["validation"],
            logger,
            False,
        )
    except Exception as exc:
        if not use_gpu:
            raise
        logger.warning("Robustness XGBoost GPU fit failed (%r); retrying on CPU", exc)
        model, _, _ = fit_xgb_regressor(
            x_sets["train"],
            y_sets["train"],
            x_sets["validation"],
            y_sets["validation"],
            dict(XGBOOST_PARAMS),
            False,
            weights["train"],
            weights["validation"],
            logger,
            False,
        )
        use_gpu = False
    scores = {}
    predictions = {}
    for split in ("validation", "test"):
        prediction = xgb_predict(model, x_sets[split])
        if CLIP_NEGATIVE_GAP_PREDICTIONS:
            prediction = np.maximum(prediction, 0.0)
        predictions[split] = prediction
        scores[split] = regression_metrics(y_sets[split], prediction)
    return {
        "scores": scores,
        "predictions": predictions,
        "best_iteration": getattr(model, "best_iteration", np.nan),
        "used_gpu": use_gpu,
    }


def run_repeated_split_robustness(
    frame: pd.DataFrame,
    arrays: dict[str, np.ndarray],
    descriptor_valid: np.ndarray,
    target: str,
    hardware: dict[str, Any],
    output: Path,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not RUN_REPEATED_SPLIT_ROBUSTNESS:
        return pd.DataFrame(), pd.DataFrame()
    target_values = frame[target].to_numpy(dtype=float)
    valid_target = descriptor_valid & np.isfinite(target_values)
    rows = []
    candidates = (
        *SOAP_REPRESENTATIONS,
        *SOAP_LOCAL_CHEMICAL_REPRESENTATIONS,
        MATCHED_SOED_NAME,
        *SOED_ENHANCED_CANDIDATES,
    )
    for seed in ROBUSTNESS_SPLIT_SEEDS:
        labels = split_labels_for_seed(frame, int(seed))
        split_indices = {
            split: np.flatnonzero(valid_target & (labels == split))
            for split in ("train", "validation", "test")
        }
        y_sets = {split: target_values[index] for split, index in split_indices.items()}
        selected_by_representation = {}
        for representation in candidates:
            x_sets = build_x_sets(representation, arrays, split_indices)
            alternatives = {}
            weighting_modes = ["tail_weighted", "unweighted"]
            if USE_MIDGAP_EXPERT and representation in MIDGAP_EXPERT_REPRESENTATIONS:
                weighting_modes.append("midgap_weighted")
            for weighting in weighting_modes:
                alternatives[weighting] = fit_robustness_candidate(
                    x_sets, y_sets, weighting, hardware, logger
                )
            if USE_VALIDATION_WEIGHT_BLEND:
                alternative_names = tuple(alternatives)
                blend_weights = validation_simplex_weights(
                    y_sets["validation"],
                    [
                        alternatives[name]["predictions"]["validation"]
                        for name in alternative_names
                    ],
                    EXPERT_BLEND_L2,
                    representation in MIDGAP_EXPERT_REPRESENTATIONS,
                )
                blended_predictions = {
                    split: sum(
                        weight * alternatives[name]["predictions"][split]
                        for name, weight in zip(alternative_names, blend_weights)
                    )
                    for split in ("validation", "test")
                }
                blended_result = {
                    "scores": {
                        split: regression_metrics(
                            y_sets[split], blended_predictions[split]
                        )
                        for split in ("validation", "test")
                    },
                    "predictions": blended_predictions,
                    "best_iteration": np.nan,
                    "used_gpu": bool(
                        alternatives["tail_weighted"]["used_gpu"]
                        or alternatives["unweighted"]["used_gpu"]
                    ),
                    "tail_weighted_fraction": float(
                        blend_weights[alternative_names.index("tail_weighted")]
                    ),
                    "midgap_weighted_fraction": float(
                        blend_weights[alternative_names.index("midgap_weighted")]
                        if "midgap_weighted" in alternative_names
                        else 0.0
                    ),
                    "expert_weights": dict(
                        zip(alternative_names, map(float, blend_weights))
                    ),
                }
                selected_by_representation[representation] = (
                    "validation_blend",
                    blended_result,
                )
            else:
                selected_weighting = min(
                    alternatives,
                    key=lambda name: alternatives[name]["scores"]["validation"]["rmse"],
                )
                selected_by_representation[representation] = (
                    selected_weighting,
                    alternatives[selected_weighting],
                )
        selected_soap_source = min(
            SOAP_REPRESENTATIONS,
            key=lambda name: selected_by_representation[name][1]["scores"][
                "validation"
            ]["rmse"],
        )
        selected_local_chemical_soap_source = min(
            SOAP_LOCAL_CHEMICAL_REPRESENTATIONS,
            key=lambda name: selected_by_representation[name][1]["scores"][
                "validation"
            ]["rmse"],
        )
        selected_soed_source = min(
            SOED_ENHANCED_CANDIDATES,
            key=lambda name: selected_by_representation[name][1]["scores"][
                "validation"
            ]["rmse"],
        )
        soed_candidate_names = tuple(SOED_ENHANCED_CANDIDATES)
        soed_ensemble_weights = validation_simplex_weights(
            y_sets["validation"],
            [
                selected_by_representation[name][1]["predictions"]["validation"]
                for name in soed_candidate_names
            ],
            SOED_ENSEMBLE_L2,
            True,
        )
        soed_ensemble_predictions = {
            split: sum(
                weight * selected_by_representation[name][1]["predictions"][split]
                for name, weight in zip(soed_candidate_names, soed_ensemble_weights)
            )
            for split in ("validation", "test")
        }
        soed_ensemble_result = {
            "scores": {
                split: regression_metrics(
                    y_sets[split], soed_ensemble_predictions[split]
                )
                for split in ("validation", "test")
            },
            "predictions": soed_ensemble_predictions,
            "best_iteration": np.nan,
            "used_gpu": any(
                selected_by_representation[name][1]["used_gpu"]
                for name in soed_candidate_names
            ),
            "tail_weighted_fraction": np.nan,
            "candidate_weights": dict(
                zip(soed_candidate_names, map(float, soed_ensemble_weights))
            ),
        }
        reporting = {
            PRIMARY_SOAP_NAME: selected_soap_source,
            MATCHED_SOED_NAME: MATCHED_SOED_NAME,
            PRIMARY_LOCAL_CHEMICAL_SOAP_NAME: selected_local_chemical_soap_source,
            PRIMARY_SOED_NAME: "prediction_ensemble",
        }
        for representation, source in reporting.items():
            if representation == PRIMARY_SOED_NAME and USE_SOED_CANDIDATE_ENSEMBLE:
                weighting, result = "validation_soed_ensemble", soed_ensemble_result
            else:
                weighting, result = selected_by_representation[
                    selected_soed_source if source == "prediction_ensemble" else source
                ]
            row = {
                "split_seed": int(seed),
                "representation": representation,
                "source_representation": source,
                "soap_alpha": SOAP_ALPHA_BY_NAME.get(
                    SOAP_LOCAL_CHEMICAL_SOURCE.get(source, source), np.nan
                ),
                "selected_weighting": weighting,
                "best_iteration": result["best_iteration"],
                "tail_weighted_fraction": result.get("tail_weighted_fraction", np.nan),
                "midgap_weighted_fraction": result.get(
                    "midgap_weighted_fraction", np.nan
                ),
                "used_gpu": result["used_gpu"],
                "n_train": len(split_indices["train"]),
                "n_validation": len(split_indices["validation"]),
                "n_test": len(split_indices["test"]),
                "n_train_groups": frame.loc[
                    split_indices["train"], "reduced_formula"
                ].nunique(),
                "n_validation_groups": frame.loc[
                    split_indices["validation"], "reduced_formula"
                ].nunique(),
                "n_test_groups": frame.loc[
                    split_indices["test"], "reduced_formula"
                ].nunique(),
            }
            for candidate_name, candidate_weight in result.get(
                "candidate_weights", {}
            ).items():
                row[f"candidate_weight_{candidate_name}"] = candidate_weight
            for split in ("validation", "test"):
                row.update(
                    {
                        f"{split}_{key}": value
                        for key, value in result["scores"][split].items()
                    }
                )
            rows.append(row)
        logger.info(
            "Robustness split seed=%d completed; selected SOAP alpha=%g local-chemical SOAP alpha=%g electronic SOED=%s",
            seed,
            SOAP_ALPHA_BY_NAME[selected_soap_source],
            SOAP_ALPHA_BY_NAME[
                SOAP_LOCAL_CHEMICAL_SOURCE[selected_local_chemical_soap_source]
            ],
            selected_soed_source,
        )
    details = pd.DataFrame(rows)
    summary = (
        details.groupby("representation", sort=False)
        .agg(
            n_seeds=("split_seed", "nunique"),
            test_rmse_mean=("test_rmse", "mean"),
            test_rmse_std=("test_rmse", "std"),
            test_rmse_min=("test_rmse", "min"),
            test_rmse_max=("test_rmse", "max"),
            test_mae_mean=("test_mae", "mean"),
            test_r2_mean=("test_r2", "mean"),
            test_tail_rmse_mean=("test_tail_rmse", "mean"),
        )
        .reset_index()
    )
    results_dir = ensure_dir(output / "results")
    details.to_csv(results_dir / f"{target}_repeated_group_splits.csv", index=False)
    save_dat(details, results_dir / f"{target}_repeated_group_splits.dat")
    summary.to_csv(
        results_dir / f"{target}_repeated_group_splits_summary.csv", index=False
    )
    save_dat(summary, results_dir / f"{target}_repeated_group_splits_summary.dat")
    return details, summary


def save_gate_result(
    task_name: str,
    representation: str,
    target: str,
    predictions: dict[str, np.ndarray],
    probabilities: dict[str, np.ndarray],
    positive_predictions: dict[str, np.ndarray],
    y_sets: dict[str, np.ndarray],
    class_sets: dict[str, np.ndarray],
    ids_sets: dict[str, np.ndarray],
    threshold: float,
    feature_dimension: int,
    training_seconds: float,
    used_gpu: bool,
    output: Path,
) -> dict[str, Any]:
    metric_rows = []
    prediction_frames = []
    for split in ("train", "validation", "test"):
        score = regression_metrics(y_sets[split], predictions[split])
        metric_rows.append(
            {
                "target": target,
                "task": task_name,
                "representation": representation,
                "model": MODEL_NAME,
                "split": split,
                **score,
                "n_samples": len(y_sets[split]),
                "feature_dimension": feature_dimension,
                "training_seconds": training_seconds,
                "used_gpu": used_gpu,
                "best_iteration": np.nan,
                "gap_threshold_ev": BAND_GAP_ZERO_THRESHOLD_EV,
                "probability_threshold": (
                    threshold if task_name == "hard_hurdle" else np.nan
                ),
            }
        )
        prediction_frames.append(
            pd.DataFrame(
                {
                    "material_id": ids_sets[split],
                    "split": split,
                    "y_true": y_sets[split],
                    "y_true_nonmetal": class_sets[split],
                    "p_nonmetal": probabilities[split],
                    "positive_gap_prediction": positive_predictions[split],
                    "y_pred": predictions[split],
                    "residual": predictions[split] - y_sets[split],
                }
            )
        )
    metric_frame = pd.DataFrame(metric_rows)
    prediction_frame = pd.concat(prediction_frames, ignore_index=True)
    run_dir = ensure_dir(
        output / "models" / target / task_name / representation / MODEL_NAME
    )
    metric_frame.to_csv(run_dir / "metrics.csv", index=False)
    save_dat(metric_frame, run_dir / "metrics.dat")
    prediction_frame.to_csv(run_dir / "predictions.csv", index=False)
    save_dat(prediction_frame, run_dir / "predictions.dat")
    return {"metrics": metric_frame, "predictions": prediction_frame}


def build_gating_summary(metrics_frame: pd.DataFrame) -> pd.DataFrame:
    labels = {
        "direct_regression": "Direct regression",
        "soft_gate": "Soft gate",
        "hard_hurdle": "Hard hurdle",
    }
    rows = []
    for task, method_label in labels.items():
        subset = validation_selected_rows(
            metrics_frame.loc[metrics_frame["representation"] == PRIMARY_SOED_NAME],
            task,
        )
        if subset.empty:
            continue
        row = subset.iloc[0].to_dict()
        row["task"] = task
        row["method_label"] = method_label
        rows.append(row)
    return pd.DataFrame(rows)


def run_target(
    target: str,
    frame: pd.DataFrame,
    arrays: dict[str, np.ndarray],
    descriptor_valid: np.ndarray,
    feature_summary: pd.DataFrame,
    hardware: dict[str, Any],
    output: Path,
    logger: logging.Logger,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    target_values = frame[target].to_numpy(dtype=float)
    target_valid = descriptor_valid & np.isfinite(target_values)
    split_indices = {
        split: np.flatnonzero(target_valid & (frame["split"].to_numpy() == split))
        for split in ("train", "validation", "test")
    }
    if any(len(indices) == 0 for indices in split_indices.values()):
        raise RuntimeError(f"Empty split for {target}")
    y_sets = {split: target_values[indices] for split, indices in split_indices.items()}
    class_sets = {
        split: (values > BAND_GAP_ZERO_THRESHOLD_EV).astype(np.int8)
        for split, values in y_sets.items()
    }
    material_ids = frame["material_id"].to_numpy(dtype=str)
    ids_sets = {
        split: material_ids[indices] for split, indices in split_indices.items()
    }
    metric_frames: list[pd.DataFrame] = []
    classification_frames: list[pd.DataFrame] = []
    direct_results: dict[str, dict[str, Any]] = {}
    classifier_results: dict[str, dict[str, Any]] = {}
    gate_results: dict[str, dict[str, dict[str, Any]]] = {}

    for representation in representation_names():
        logger.info(
            "Training target=%s representation=%s model=xgboost", target, representation
        )
        x_sets = build_x_sets(representation, arrays, split_indices)
        if not all(np.all(np.isfinite(values)) for values in x_sets.values()):
            raise ValueError(f"Non-finite features in {representation}")
        if RUN_DIRECT_REGRESSION:
            weighted = fit_one_regression(
                representation,
                target,
                "direct_regression_tail_weighted",
                x_sets,
                y_sets,
                ids_sets,
                hardware,
                output,
                logger,
                weighting_mode="tail_weighted",
            )
            unweighted = fit_one_regression(
                representation,
                target,
                "direct_regression_unweighted",
                x_sets,
                y_sets,
                ids_sets,
                hardware,
                output,
                logger,
                weighting_mode="unweighted",
            )
            metric_frames.extend((weighted["metrics"], unweighted["metrics"]))
            alternatives = {"tail_weighted": weighted, "unweighted": unweighted}
            if USE_MIDGAP_EXPERT and representation in MIDGAP_EXPERT_REPRESENTATIONS:
                midgap_weighted = fit_one_regression(
                    representation,
                    target,
                    "direct_regression_midgap_weighted",
                    x_sets,
                    y_sets,
                    ids_sets,
                    hardware,
                    output,
                    logger,
                    weighting_mode="midgap_weighted",
                )
                alternatives["midgap_weighted"] = midgap_weighted
                metric_frames.append(midgap_weighted["metrics"])
            if USE_VALIDATION_WEIGHT_BLEND:
                direct = materialize_validation_expert_blend(
                    alternatives,
                    representation,
                    target,
                    output,
                    protect_midgap=(representation in MIDGAP_EXPERT_REPRESENTATIONS),
                )
                selected_weighting = direct["selected_weighting"]
            else:
                selected_weighting = min(
                    alternatives,
                    key=lambda name: result_validation_rmse(alternatives[name]),
                )
                direct = materialize_selected_regression(
                    alternatives[selected_weighting],
                    representation,
                    representation,
                    selected_weighting,
                    target,
                    output,
                )
            direct_results[representation] = direct
            metric_frames.append(direct["metrics"])
            logger.info(
                "Validation selected weighting=%s for %s (RMSE=%.6g)",
                selected_weighting,
                representation,
                result_validation_rmse(direct),
            )
        if not RUN_HURDLE_AUXILIARY or representation not in HURDLE_REPRESENTATIONS:
            continue
        classifier = fit_one_classifier(
            representation,
            target,
            x_sets,
            class_sets,
            ids_sets,
            hardware,
            output,
            logger,
        )
        classifier_results[representation] = classifier
        classification_frames.append(classifier["metrics"])
        positive_masks = {
            split: y_sets[split] > BAND_GAP_ZERO_THRESHOLD_EV for split in y_sets
        }
        positive = fit_one_regression(
            representation,
            target,
            "positive_regression",
            {split: x_sets[split][positive_masks[split]] for split in x_sets},
            {split: y_sets[split][positive_masks[split]] for split in y_sets},
            {split: ids_sets[split][positive_masks[split]] for split in ids_sets},
            hardware,
            output,
            logger,
            weighting_mode="tail_weighted",
        )
        metric_frames.append(positive["metrics"])
        positive_predictions = {
            split: np.maximum(xgb_predict(positive["model_object"], values), 0.0)
            for split, values in x_sets.items()
        }
        gate_predictions: dict[str, dict[str, np.ndarray]] = {}
        if RUN_SOFT_GATE:
            gate_predictions["soft_gate"] = {
                split: classifier["probabilities"][split] * positive_predictions[split]
                for split in x_sets
            }
        if RUN_HARD_GATE:
            gate_predictions["hard_hurdle"] = {
                split: np.where(
                    classifier["probabilities"][split]
                    >= classifier["decision_threshold"],
                    positive_predictions[split],
                    0.0,
                )
                for split in x_sets
            }
        gate_results[representation] = {}
        total_seconds = float(
            classifier["metrics"]["training_seconds"].iloc[0]
            + positive["metrics"]["training_seconds"].iloc[0]
        )
        used_gpu = bool(
            classifier["metrics"]["used_gpu"].iloc[0]
            or positive["metrics"]["used_gpu"].iloc[0]
        )
        for task_name, prediction_values in gate_predictions.items():
            result = save_gate_result(
                task_name,
                representation,
                target,
                prediction_values,
                classifier["probabilities"],
                positive_predictions,
                y_sets,
                class_sets,
                ids_sets,
                classifier["decision_threshold"],
                x_sets["train"].shape[1],
                total_seconds,
                used_gpu,
                output,
            )
            gate_results[representation][task_name] = result
            metric_frames.append(result["metrics"])
            gate_validation_rmse = float(
                result["metrics"]
                .loc[result["metrics"]["split"] == "validation", "rmse"]
                .iloc[0]
            )
            gate_test_rmse = float(
                result["metrics"]
                .loc[result["metrics"]["split"] == "test", "rmse"]
                .iloc[0]
            )
            logger.info(
                "%s target=%s representation=%s validation_RMSE=%.6g test_RMSE=%.6g",
                task_name.upper(),
                target,
                representation,
                gate_validation_rmse,
                gate_test_rmse,
            )

    selected_soap_source = min(
        SOAP_REPRESENTATIONS,
        key=lambda name: result_validation_rmse(direct_results[name]),
    )
    soap_source = direct_results[selected_soap_source]
    direct_results[PRIMARY_SOAP_NAME] = materialize_selected_regression(
        soap_source,
        selected_soap_source,
        PRIMARY_SOAP_NAME,
        soap_source["selected_weighting"],
        target,
        output,
    )
    metric_frames.append(direct_results[PRIMARY_SOAP_NAME]["metrics"])
    selected_alpha = SOAP_ALPHA_BY_NAME[selected_soap_source]
    REPRESENTATION_LABELS[PRIMARY_SOAP_NAME] = (
        rf"Periodic SOAP ($\alpha={selected_alpha:g}$)"
    )
    logger.info(
        "Validation selected SOAP alpha=%g source=%s weighting=%s",
        selected_alpha,
        selected_soap_source,
        soap_source["selected_weighting"],
    )

    selected_chemical_soap_source = min(
        SOAP_COMPOSITION_REPRESENTATIONS,
        key=lambda name: result_validation_rmse(direct_results[name]),
    )
    chemical_soap_source = direct_results[selected_chemical_soap_source]
    direct_results[PRIMARY_CHEMICAL_SOAP_NAME] = materialize_selected_regression(
        chemical_soap_source,
        selected_chemical_soap_source,
        PRIMARY_CHEMICAL_SOAP_NAME,
        chemical_soap_source["selected_weighting"],
        target,
        output,
    )
    metric_frames.append(direct_results[PRIMARY_CHEMICAL_SOAP_NAME]["metrics"])
    selected_chemical_source = SOAP_COMPOSITION_SOURCE[selected_chemical_soap_source]
    selected_chemical_alpha = SOAP_ALPHA_BY_NAME[selected_chemical_source]
    REPRESENTATION_LABELS[PRIMARY_CHEMICAL_SOAP_NAME] = (
        rf"Periodic SOAP + composition ($\alpha={selected_chemical_alpha:g}$)"
    )
    logger.info(
        "Validation selected composition-matched SOAP alpha=%g source=%s weighting=%s",
        selected_chemical_alpha,
        selected_chemical_soap_source,
        chemical_soap_source["selected_weighting"],
    )

    selected_local_chemical_soap_source = min(
        SOAP_LOCAL_CHEMICAL_REPRESENTATIONS,
        key=lambda name: result_validation_rmse(direct_results[name]),
    )
    local_chemical_soap_source = direct_results[selected_local_chemical_soap_source]
    direct_results[PRIMARY_LOCAL_CHEMICAL_SOAP_NAME] = materialize_selected_regression(
        local_chemical_soap_source,
        selected_local_chemical_soap_source,
        PRIMARY_LOCAL_CHEMICAL_SOAP_NAME,
        local_chemical_soap_source["selected_weighting"],
        target,
        output,
    )
    metric_frames.append(direct_results[PRIMARY_LOCAL_CHEMICAL_SOAP_NAME]["metrics"])
    selected_local_source = SOAP_LOCAL_CHEMICAL_SOURCE[
        selected_local_chemical_soap_source
    ]
    selected_local_alpha = SOAP_ALPHA_BY_NAME[selected_local_source]
    REPRESENTATION_LABELS[PRIMARY_LOCAL_CHEMICAL_SOAP_NAME] = (
        rf"Local-chemical SOAP ($\alpha={selected_local_alpha:g}$)"
    )
    logger.info(
        "Validation selected local-chemical SOAP alpha=%g source=%s weighting=%s",
        selected_local_alpha,
        selected_local_chemical_soap_source,
        local_chemical_soap_source["selected_weighting"],
    )

    selected_soed_source = min(
        SOED_ENHANCED_CANDIDATES,
        key=lambda name: result_validation_rmse(direct_results[name]),
    )
    selected_soed = direct_results[selected_soed_source]
    if USE_SOED_CANDIDATE_ENSEMBLE:
        direct_results[PRIMARY_SOED_NAME] = materialize_soed_candidate_ensemble(
            {name: direct_results[name] for name in SOED_ENHANCED_CANDIDATES},
            target,
            output,
        )
    else:
        direct_results[PRIMARY_SOED_NAME] = materialize_selected_regression(
            selected_soed,
            selected_soed_source,
            PRIMARY_SOED_NAME,
            selected_soed["selected_weighting"],
            target,
            output,
        )
    metric_frames.append(direct_results[PRIMARY_SOED_NAME]["metrics"])
    REPRESENTATION_LABELS[PRIMARY_SOED_NAME] = (
        "Validation-selected electronic SOED ensemble"
        if USE_SOED_CANDIDATE_ENSEMBLE
        else "Validation-selected " + SOED_CANDIDATE_LABELS[selected_soed_source]
    )
    logger.info(
        "Validation selected electronic SOED source=%s weighting=%s; primary=%s RMSE=%.6g",
        selected_soed_source,
        selected_soed["selected_weighting"],
        "candidate_ensemble" if USE_SOED_CANDIDATE_ENSEMBLE else selected_soed_source,
        result_validation_rmse(direct_results[PRIMARY_SOED_NAME]),
    )
    if USE_SOED_CANDIDATE_ENSEMBLE:
        logger.info(
            "SOED_ENSEMBLE_WEIGHTS %s",
            json.dumps(
                direct_results[PRIMARY_SOED_NAME].get("candidate_weights", {}),
                sort_keys=True,
            ),
        )

    if selected_soed_source in classifier_results:
        classifier_source = classifier_results[selected_soed_source]
        classifier_metrics = classifier_source["metrics"].copy()
        classifier_metrics["source_representation"] = selected_soed_source
        classifier_metrics["representation"] = PRIMARY_SOED_NAME
        classifier_predictions = classifier_source["predictions"].copy()
        classifier_predictions["source_representation"] = selected_soed_source
        classifier_results[PRIMARY_SOED_NAME] = {
            **classifier_source,
            "metrics": classifier_metrics,
            "predictions": classifier_predictions,
        }
        classification_frames.append(classifier_metrics)
        classifier_dir = ensure_dir(
            output
            / "models"
            / target
            / "metal_nonmetal_classification"
            / PRIMARY_SOED_NAME
            / MODEL_NAME
        )
        classifier_metrics.to_csv(classifier_dir / "metrics.csv", index=False)
        save_dat(classifier_metrics, classifier_dir / "metrics.dat")
        classifier_predictions.to_csv(classifier_dir / "predictions.csv", index=False)
        save_dat(classifier_predictions, classifier_dir / "predictions.dat")

    if selected_soed_source in gate_results:
        gate_results[PRIMARY_SOED_NAME] = {}
        for task_name, source_gate in gate_results[selected_soed_source].items():
            gate_metrics = source_gate["metrics"].copy()
            gate_metrics["source_representation"] = selected_soed_source
            gate_metrics["representation"] = PRIMARY_SOED_NAME
            gate_predictions = source_gate["predictions"].copy()
            gate_predictions["source_representation"] = selected_soed_source
            alias_gate = {
                **source_gate,
                "metrics": gate_metrics,
                "predictions": gate_predictions,
            }
            gate_results[PRIMARY_SOED_NAME][task_name] = alias_gate
            metric_frames.append(gate_metrics)
            gate_dir = ensure_dir(
                output / "models" / target / task_name / PRIMARY_SOED_NAME / MODEL_NAME
            )
            gate_metrics.to_csv(gate_dir / "metrics.csv", index=False)
            save_dat(gate_metrics, gate_dir / "metrics.dat")
            gate_predictions.to_csv(gate_dir / "predictions.csv", index=False)
            save_dat(gate_predictions, gate_dir / "predictions.dat")

    required = {
        PRIMARY_SOAP_NAME,
        PRIMARY_CHEMICAL_SOAP_NAME,
        PRIMARY_LOCAL_CHEMICAL_SOAP_NAME,
        COMPOSITION_ONLY_NAME,
        MATCHED_SOED_NAME,
        MATCHED_CHEMICAL_SOED_NAME,
        DENSITY_COMPOSITION_NAME,
        PRIMARY_SOED_NAME,
    }
    if not required.issubset(direct_results):
        raise RuntimeError(
            f"Missing primary direct results: {required - set(direct_results)}"
        )
    all_metrics = pd.concat(metric_frames, ignore_index=True)
    all_classification = (
        pd.concat(classification_frames, ignore_index=True)
        if classification_frames
        else pd.DataFrame()
    )
    results_dir = ensure_dir(output / "results")
    all_metrics.to_csv(results_dir / f"{target}_all_metrics.csv", index=False)
    save_dat(all_metrics, results_dir / f"{target}_all_metrics.dat")
    all_classification.to_csv(
        results_dir / f"{target}_classification_metrics.csv", index=False
    )
    save_dat(all_classification, results_dir / f"{target}_classification_metrics.dat")

    robustness, robustness_summary = run_repeated_split_robustness(
        frame, arrays, descriptor_valid, target, hardware, output, logger
    )

    paired_matched, bootstrap_matched, comparison_matched = bootstrap_distribution(
        direct_results[MATCHED_SOED_NAME]["predictions"],
        direct_results[PRIMARY_SOAP_NAME]["predictions"],
        MATCHED_SOED_NAME,
        PRIMARY_SOAP_NAME,
        "dimension-matched density SOED versus SOAP",
        frame,
    )
    paired_local_matched, bootstrap_local_matched, comparison_local_matched = (
        bootstrap_distribution(
            direct_results[MATCHED_CHEMICAL_SOED_NAME]["predictions"],
            direct_results[PRIMARY_LOCAL_CHEMICAL_SOAP_NAME]["predictions"],
            MATCHED_CHEMICAL_SOED_NAME,
            PRIMARY_LOCAL_CHEMICAL_SOAP_NAME,
            "local-chemistry-matched density SOED versus SOAP",
            frame,
        )
    )
    paired_final, bootstrap_final, comparison_final = bootstrap_distribution(
        direct_results[PRIMARY_SOED_NAME]["predictions"],
        direct_results[PRIMARY_LOCAL_CHEMICAL_SOAP_NAME]["predictions"],
        PRIMARY_SOED_NAME,
        PRIMARY_LOCAL_CHEMICAL_SOAP_NAME,
        "validation-selected electronic SOED versus local-chemical SOAP",
        frame,
    )
    paired_fair, bootstrap_fair, comparison_fair = bootstrap_distribution(
        direct_results[DENSITY_COMPOSITION_NAME]["predictions"],
        direct_results[PRIMARY_CHEMICAL_SOAP_NAME]["predictions"],
        DENSITY_COMPOSITION_NAME,
        PRIMARY_CHEMICAL_SOAP_NAME,
        "composition-matched density SOED versus SOAP",
        frame,
    )
    paired_final_tail, bootstrap_final_tail, comparison_final_tail = (
        bootstrap_distribution(
            direct_results[PRIMARY_SOED_NAME]["predictions"],
            direct_results[PRIMARY_LOCAL_CHEMICAL_SOAP_NAME]["predictions"],
            PRIMARY_SOED_NAME,
            PRIMARY_LOCAL_CHEMICAL_SOAP_NAME,
            "validation-selected electronic SOED versus local-chemical SOAP, high-gap subset",
            frame,
            analysis_scope="high_gap_ge_3_ev",
            minimum_gap_ev=TAIL_THRESHOLD_EV,
        )
    )
    gap_bootstrap_outputs: list[tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]] = []
    gap_scopes = (
        ("gap_le_0p01_ev", None, np.nextafter(BAND_GAP_ZERO_THRESHOLD_EV, np.inf)),
        ("gap_0p01_to_0p5_ev", np.nextafter(BAND_GAP_ZERO_THRESHOLD_EV, np.inf), 0.5),
        ("gap_0p5_to_1_ev", 0.5, 1.0),
        ("gap_1_to_2_ev", 1.0, 2.0),
        ("gap_2_to_3_ev", 2.0, 3.0),
        ("gap_3_to_6_ev", 3.0, 6.0),
        ("gap_gt_6_ev", 6.0, None),
    )
    for scope, minimum_gap, maximum_gap in gap_scopes:
        gap_bootstrap_outputs.append(
            bootstrap_distribution(
                direct_results[PRIMARY_SOED_NAME]["predictions"],
                direct_results[PRIMARY_LOCAL_CHEMICAL_SOAP_NAME]["predictions"],
                PRIMARY_SOED_NAME,
                PRIMARY_LOCAL_CHEMICAL_SOAP_NAME,
                f"validation-selected electronic SOED versus local-chemical SOAP, {scope}",
                frame,
                analysis_scope=scope,
                minimum_gap_ev=minimum_gap,
                maximum_gap_ev=maximum_gap,
            )
        )
    paired = pd.concat(
        (
            paired_matched,
            paired_local_matched,
            paired_final,
            paired_fair,
            paired_final_tail,
            *(value[0] for value in gap_bootstrap_outputs),
        ),
        ignore_index=True,
    )
    bootstrap = pd.concat(
        (
            bootstrap_matched,
            bootstrap_local_matched,
            bootstrap_final,
            bootstrap_fair,
            bootstrap_final_tail,
            *(value[1] for value in gap_bootstrap_outputs),
        ),
        ignore_index=True,
    )
    comparison = pd.concat(
        (
            comparison_matched,
            comparison_local_matched,
            comparison_final,
            comparison_fair,
            comparison_final_tail,
            *(value[2] for value in gap_bootstrap_outputs),
        ),
        ignore_index=True,
    )
    gap_mask = comparison["analysis_scope"].str.startswith("gap_")
    comparison["bootstrap_p_fdr_bh"] = np.nan
    if gap_mask.any():
        p_values = comparison.loc[gap_mask, "bootstrap_one_sided_p"].to_numpy(
            dtype=float
        )
        order = np.argsort(p_values)
        adjusted_sorted = np.minimum.accumulate(
            (p_values[order] * len(p_values) / np.arange(1, len(p_values) + 1))[::-1]
        )[::-1]
        adjusted = np.empty_like(adjusted_sorted)
        adjusted[order] = np.clip(adjusted_sorted, 0.0, 1.0)
        comparison.loc[gap_mask, "bootstrap_p_fdr_bh"] = adjusted
    comparison.to_csv(
        results_dir / f"{target}_soed_vs_soap_statistical_comparison.csv", index=False
    )
    save_dat(
        comparison, results_dir / f"{target}_soed_vs_soap_statistical_comparison.dat"
    )
    for row in comparison.to_dict(orient="records"):
        logger.info("SOED_VS_SOAP %s", json.dumps(row, default=json_default))

    selected_feature = feature_summary.loc[
        feature_summary["representation"] == selected_soap_source
    ].copy()
    selected_feature["source_representation"] = selected_soap_source
    selected_feature["representation"] = PRIMARY_SOAP_NAME
    selected_feature["soap_alpha"] = selected_alpha
    composition_feature = feature_summary.loc[
        feature_summary["representation"] == COMPOSITION_ONLY_NAME
    ].copy()
    selected_chemical_feature = feature_summary.loc[
        feature_summary["representation"] == selected_chemical_source
    ].copy()
    selected_chemical_feature["source_representation"] = selected_chemical_soap_source
    selected_chemical_feature["representation"] = PRIMARY_CHEMICAL_SOAP_NAME
    selected_chemical_feature["soap_alpha"] = selected_chemical_alpha
    selected_chemical_feature["feature_dimension"] += int(
        composition_feature["feature_dimension"].iloc[0]
    )
    selected_chemical_feature["seconds_per_structure"] += float(
        composition_feature["seconds_per_structure"].iloc[0]
    )
    selected_local_chemical_feature = feature_summary.loc[
        feature_summary["representation"] == selected_local_chemical_soap_source
    ].copy()
    selected_local_chemical_feature["source_representation"] = (
        selected_local_chemical_soap_source
    )
    selected_local_chemical_feature["representation"] = PRIMARY_LOCAL_CHEMICAL_SOAP_NAME
    selected_local_chemical_feature["soap_alpha"] = selected_local_alpha
    selected_soed_feature = (
        feature_summary.loc[
            feature_summary["representation"].isin(SOED_ENHANCED_CANDIDATES)
        ]
        .sort_values("feature_dimension")
        .tail(1)
        .copy()
    )
    selected_soed_feature["source_representation"] = (
        "prediction_ensemble" if USE_SOED_CANDIDATE_ENSEMBLE else selected_soed_source
    )
    selected_soed_feature["representation"] = PRIMARY_SOED_NAME
    selected_soed_feature["ensemble_members"] = (
        len(SOED_ENHANCED_CANDIDATES) if USE_SOED_CANDIDATE_ENSEMBLE else 1
    )
    density_composition_feature = feature_summary.loc[
        feature_summary["representation"] == MATCHED_SOED_NAME
    ].copy()
    density_composition_feature["source_representation"] = MATCHED_SOED_NAME
    density_composition_feature["representation"] = DENSITY_COMPOSITION_NAME
    density_composition_feature["feature_dimension"] += int(
        composition_feature["feature_dimension"].iloc[0]
    )
    density_composition_feature["seconds_per_structure"] += float(
        composition_feature["seconds_per_structure"].iloc[0]
    )
    feature_summary_plot = pd.concat(
        (
            selected_feature,
            feature_summary.loc[
                feature_summary["representation"].isin(SOED_CHANNEL_SETS)
            ],
            composition_feature,
            selected_chemical_feature,
            selected_local_chemical_feature,
            selected_soed_feature,
            density_composition_feature,
        ),
        ignore_index=True,
    )
    direct_predictions = {
        name: result["predictions"] for name, result in direct_results.items()
    }
    direct_histories = {
        name: result["history"] for name, result in direct_results.items()
    }
    mechanism, midgap_cases, block_importance = save_physical_interpretation(
        direct_results,
        direct_results[PRIMARY_SOED_NAME].get("candidate_weights", {}),
        output,
        target,
    )
    subset_metrics = subset_regression_metrics(direct_predictions)
    stratified = plot_main_figures(
        frame.loc[target_valid].copy(),
        target,
        all_metrics,
        direct_predictions,
        direct_histories,
        feature_summary_plot,
        paired,
        bootstrap,
        comparison,
        robustness,
        mechanism,
        block_importance,
        output,
    )
    gating_summary = build_gating_summary(all_metrics)
    classification_predictions = classifier_results.get(PRIMARY_SOED_NAME, {}).get(
        "predictions"
    )
    plot_si_figures(
        frame.loc[target_valid].copy(),
        target,
        all_metrics,
        direct_predictions,
        direct_histories,
        classification_predictions,
        gating_summary,
        feature_summary_plot,
        stratified,
        selected_soap_source,
        selected_local_chemical_soap_source,
        selected_soed_source,
        robustness,
        bootstrap,
        comparison,
        mechanism,
        midgap_cases,
        block_importance,
        output,
    )
    save_tables(
        all_metrics,
        all_classification,
        gating_summary,
        stratified,
        feature_summary_plot,
        paired,
        bootstrap,
        comparison,
        robustness,
        robustness_summary,
        subset_metrics,
        mechanism,
        midgap_cases,
        block_importance,
        target,
        output,
    )
    return all_metrics, comparison


def main() -> int:
    global FORCE_RECOMPUTE_STRUCTURES, FORCE_RECOMPUTE_FEATURES
    args = parse_args()
    if args.force_structures:
        FORCE_RECOMPUTE_STRUCTURES = True
        FORCE_RECOMPUTE_FEATURES = True
    if args.force_features:
        FORCE_RECOMPUTE_FEATURES = True
    output = ensure_dir(args.output.expanduser().resolve())
    logger = setup_logger(output)
    logger.info("Script build: %s | file=%s", SCRIPT_BUILD, Path(__file__).resolve())
    started = time.perf_counter()
    seed_everything(RANDOM_SEED)
    create_workflow_files(output)
    if args.workflow_only:
        logger.info("Workflow files written to %s", output / "main_figures")
        return 0
    config = config_snapshot()
    config["dataset"]["max_samples"] = args.limit
    config["cache"]["root"] = str(args.cache.expanduser().resolve())
    (output / "run_config.json").write_text(
        json.dumps(config, indent=2, ensure_ascii=False, default=json_default),
        encoding="utf-8",
    )
    protocol_payload = {
        "script_build": SCRIPT_BUILD,
        "status": (
            "frozen_confirmatory" if CONFIRMATORY_PROTOCOL_FROZEN else "exploratory"
        ),
        "main_seed": RANDOM_SEED,
        "outer_seeds": ROBUSTNESS_SPLIT_SEEDS,
        "legacy_seed_excluded": EXPLORATORY_LEGACY_SEED,
        "selection_data": "training and formula-disjoint validation only",
        "test_feedback_allowed": False,
        "config_sha256": hashlib.sha256(
            json.dumps(config, sort_keys=True, default=json_default).encode("utf-8")
        ).hexdigest(),
    }
    (output / "confirmatory_protocol_lock.json").write_text(
        json.dumps(protocol_payload, indent=2), encoding="utf-8"
    )
    logger.info(
        "Main parameters:\n%s",
        json.dumps(config, indent=2, ensure_ascii=False, default=json_default),
    )
    hardware = detect_device(logger)
    root = resolve_dataset_root(args.root)
    frame, structure_store = prepare_structure_cache(
        load_metadata(root, args.limit, logger), args.cache, logger
    )
    frame = exclude_undefined_pauling_structures(frame, structure_store, output, logger)
    if len(frame) < 10:
        raise RuntimeError(
            "Too few structures remain after undefined-Pauling filtering"
        )
    frame = assign_group_stratified_splits(frame, output, logger)
    arrays, descriptor_valid, feature_summary = generate_features(
        frame, root, output, structure_store, logger
    )
    all_metrics = []
    comparison_results = []
    for target in TARGETS:
        target_metrics, comparison = run_target(
            target,
            frame,
            arrays,
            descriptor_valid,
            feature_summary,
            hardware,
            output,
            logger,
        )
        all_metrics.append(target_metrics)
        comparison_results.append(comparison)
    combined = pd.concat(all_metrics, ignore_index=True)
    combined_comparison = pd.concat(comparison_results, ignore_index=True)
    combined.to_csv(output / "results" / "all_targets_metrics.csv", index=False)
    save_dat(combined, output / "results" / "all_targets_metrics.dat")
    combined_comparison.to_csv(
        output / "results" / "all_targets_statistical_comparisons.csv", index=False
    )
    save_dat(
        combined_comparison,
        output / "results" / "all_targets_statistical_comparisons.dat",
    )
    elapsed = time.perf_counter() - started
    primary_test = combined.loc[
        (combined["task"] == "direct_regression")
        & (combined["split"] == "test")
        & (combined["representation"] == PRIMARY_SOED_NAME)
    ].iloc[0]
    summary = {
        "status": "completed",
        "runtime_seconds": elapsed,
        "runtime_hours": elapsed / 3600.0,
        "dataset_root": str(root),
        "output_root": str(output),
        "samples": len(frame),
        "descriptor_valid_samples": int(descriptor_valid.sum()),
        "split_counts": frame["split"].value_counts().to_dict(),
        "model": MODEL_NAME,
        "confirmatory_protocol": protocol_payload,
        "primary_representation": PRIMARY_SOED_NAME,
        "primary_test": primary_test.to_dict(),
        "soed_vs_soap": combined_comparison.to_dict(orient="records"),
        "hardware": hardware,
    }
    (output / "run_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, default=json_default),
        encoding="utf-8",
    )
    logger.info("Run completed in %.1f seconds (%.3f hours)", elapsed, elapsed / 3600.0)
    logger.info("Outputs: %s", output)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print(
            "Interrupted by user; completed structure and descriptor caches remain reusable.",
            file=sys.stderr,
        )
        raise

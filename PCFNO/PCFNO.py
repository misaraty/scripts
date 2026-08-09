import argparse
import atexit
import faulthandler
import gc
import csv
import hashlib
import json
import logging
import math
import os
import platform
import random
import shutil
import sys
import time
import tempfile
from contextlib import nullcontext
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Sequence, Tuple

os.environ.setdefault("TORCH_SHOW_CPP_STACKTRACES", "1")
os.environ.setdefault("MALLOC_ARENA_MAX", "2")
os.environ.setdefault("CUDA_MODULE_LOADING", "LAZY")
os.environ.setdefault("OMP_NUM_THREADS", str(max(1, min(4, os.cpu_count() or 1))))
os.environ.setdefault("MKL_NUM_THREADS", str(max(1, min(4, os.cpu_count() or 1))))
os.environ.setdefault(
    "PYTORCH_ALLOC_CONF",
    "max_split_size_mb:128,garbage_collection_threshold:0.80",
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.ndimage import binary_dilation, binary_erosion, distance_transform_edt
from scipy.signal import detrend, welch
from scipy.stats import wilcoxon
from torch.utils.data import DataLoader, Dataset

CONFIG: Dict[str, Any] = {
    "data_path": "./cylinder",
    "dataset_type": "prop",
    "output_dir": "./PCFNO_PROP_speed_PoF_results_v15_fixed_blend_Re24_180_900_seed44_45_46",
    "resolution": 64,
    "batch_size": 8,
    "epochs": 300,
    "lr": 1e-3,
    "weight_decay": 1e-6,
    "width": 64,
    "modes": 16,
    "fno_layers": 4,
    "padding": 24,
    "time_samples_per_case": 64,
    "spectrum_time_samples": 256,
    "data_delta_time": 0.001,
    "train_fraction": 0.8,
    "val_fraction": 0.1,
    "split_seed": 42,
    "split_mode": "grouped_re",
    "split_profile": "confirmation",
    "split_profiles": {
        "development": {
            "test_reynolds": [20.0, 160.0, 1000.0],
            "validation_reynolds": [12.0, 60.0, 240.0, 500.0, 800.0],
            "excluded_reynolds": [],
        },
        "confirmation": {
            "test_reynolds": [24.0, 180.0, 900.0],
            "validation_reynolds": [16.0, 80.0, 320.0, 500.0, 800.0],
            "excluded_reynolds": [],
        },
    },
    "test_reynolds": [24.0, 180.0, 900.0],
    "validation_reynolds": [16.0, 80.0, 320.0, 500.0, 800.0],
    "excluded_reynolds": [],
    "re_group_decimals": 8,
    "training_seeds": [44, 45, 46],
    "representative_re_targets": [24.0, 180.0, 900.0],
    "num_workers": 0,
    "pin_memory": False,
    "use_amp": False,
    "amp_dtype": "auto",
    "amp_init_scale": 1024.0,
    "amp_growth_factor": 2.0,
    "amp_backoff_factor": 0.5,
    "amp_growth_interval": 2000,
    "max_nonfinite_batches_per_epoch": 20,
    "spectral_micro_batch_size": 8,
    "cuda_maintenance_interval": 5,
    "cufft_plan_cache_max_size": 16,
    "checkpoint_interval_epochs": 5,
    "save_history_every_epoch": True,
    "enable_native_faulthandler": True,
    "log_resource_usage": True,
    "early_stopping_patience": 45,
    "early_stopping_min_delta": 1e-6,
    "minimum_epochs_before_early_stopping": 120,
    "validation_wake_weight": 0.02,
    "validation_boundary_weight": 0.001,
    "validation_initial_weight": 0.01,
    "validation_high_re_weight": 0.02,
    "validation_high_re_threshold_scaled": 0.80,
    "physics_warmup_fraction": 0.05,
    "physics_warmup_epochs": 5,
    "pcfno_data_pretrain_fraction": 0.0,
    "lambda_relative": 0.03,
    "lambda_gradient": 0.005,
    "lambda_boundary": 0.001,
    "lambda_wake": 0.0,
    "lambda_non_degradation": 0.0,
    "lambda_residual_regularization": 0.002,
    "lambda_uncertainty": 0.003,
    "lambda_calibration": 0.00002,
    "positive_output_mode": "smooth_relu",
    "positive_output_epsilon": 1e-4,
    "enforce_positive_output": True,
    "initial_time_threshold_scaled": 0.02,
    "initial_sample_weight": 1.10,
    "gradient_exclusion_pixels": 1,
    "pcfno_warm_start_from_fno": True,
    "pcfno_freeze_backbone": True,
    "pcfno_adapter_lr": 5.0e-4,
    "adapter_width": 32,
    "pcfno_max_epochs": 24,
    "pcfno_minimum_epochs_before_early_stopping": 12,
    "pcfno_early_stopping_patience": 8,
    "adapter_correction_limit": 0.20,
    "adapter_focus_floor": 0.20,
    "adapter_boundary_scale_over_d": 0.75,
    "adapter_wake_x_smoothing": 0.30,
    "adapter_wake_y_smoothing": 0.30,
    "fixed_residual_alpha": 0.50,
    "unblended_residual_alpha": 1.0,
    "uncertainty_refinement_epochs": 12,
    "uncertainty_refinement_lr": 1.0e-3,
    "uncertainty_refinement_patience": 4,
    "uncertainty_refinement_min_delta": 1.0e-6,
    "uncertainty_refinement_log_weight": 0.01,
    "strict_re_matching": True,
    "runtime_warmup": 10,
    "runtime_repeats": 30,
    "cfd_reference_seconds_per_frame": 1.18,
    "save_jpg_dpi": 600,
    "resume_existing": True,
    "force_retrain": False,
    "uncertainty_subsample": 200000,
    "posthoc_sigma_calibration": True,
    "sigma_scale_min": 0.25,
    "sigma_scale_max": 20.0,
    "hard_zero_solid": True,
    "wake_x_min_over_d": 0.5,
    "wake_x_max_over_d": 8.0,
    "wake_half_height_over_d": 2.0,
    "near_wall_distance_over_d": 0.5,
    "bootstrap_repeats": 5000,
    "confidence_level": 0.95,
    "uncertainty_plot_percentile": 99.5,
    "uncertainty_plot_gridsize": 50,
    "max_train_batches": None,
    "max_eval_batches": None,
    "device": "cuda" if torch.cuda.is_available() else "cpu",
}

INPUT_XD = 0
INPUT_YD = 1
INPUT_SOLID = 2
INPUT_SDF = 3
INPUT_BOUNDARY = 4
INPUT_RE = 5
INPUT_UIN = 6
INPUT_TIME = 7
INPUT_CHANNELS = 8
OUTPUT_CHANNELS = 1

BASELINE_NAMES = ["UNet", "DeepONet", "FNO", "PCFNO"]
ABLATION_NAMES = [
    "PCFNO_without_SDF_adapter",
    "PCFNO_without_physics_focus",
    "PCFNO_without_physics_losses",
    "PCFNO_without_residual_blending",
]
ALL_EXPERIMENTS = BASELINE_NAMES + ABLATION_NAMES

MODEL_COLORS = {
    "UNet": "#0072B2",
    "DeepONet": "#E69F00",
    "FNO": "#009E73",
    "PCFNO": "#D55E00",
    "PCFNO_without_SDF_adapter": "#56B4E9",
    "PCFNO_without_physics_focus": "#CC79A7",
    "PCFNO_without_physics_losses": "#E69F00",
    "PCFNO_without_residual_blending": "#000000",
}

DISPLAY_NAMES = {
    "UNet": "U-Net",
    "DeepONet": "DeepONet",
    "FNO": "FNO",
    "PCFNO": "PCFNO",
    "PCFNO_without_SDF_adapter": "w/o SDF adapter",
    "PCFNO_without_physics_focus": "w/o physics focus",
    "PCFNO_without_physics_losses": "w/o physics losses",
    "PCFNO_without_residual_blending": "w/o residual blending",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Speed-magnitude PCFNO study for Physics of Fluids"
    )
    parser.add_argument("--data-path", type=str, default=None)
    parser.add_argument("--output-dir", type=str, default=None)
    parser.add_argument("--epochs", type=int, default=None)
    parser.add_argument("--batch-size", type=int, default=None)
    parser.add_argument("--resolution", type=int, default=None)
    parser.add_argument("--time-samples", type=int, default=None)
    parser.add_argument("--device", type=str, default=None)
    parser.add_argument(
        "--split-profile",
        type=str,
        choices=["development", "confirmation"],
        default=None,
        help="Grouped-Re split profile. confirmation is the default.",
    )
    parser.add_argument("--seeds", type=int, nargs="+", default=None)
    parser.add_argument("--models", type=str, nargs="+", default=None)
    parser.add_argument("--smoke-test", action="store_true")
    parser.add_argument("--skip-plots", action="store_true")
    parser.add_argument("--force-retrain", action="store_true")
    parser.add_argument("--no-resume", action="store_true")
    parser.add_argument(
        "--spectral-micro-batch-size",
        type=int,
        default=None,
        help="Micro-batch size for FNO/PCFNO. The effective loader batch remains unchanged.",
    )
    return parser.parse_args()


def apply_cli_overrides(args: argparse.Namespace) -> None:
    mapping = {
        "data_path": args.data_path,
        "output_dir": args.output_dir,
        "epochs": args.epochs,
        "batch_size": args.batch_size,
        "resolution": args.resolution,
        "time_samples_per_case": args.time_samples,
        "device": args.device,
        "split_profile": args.split_profile,
        "spectral_micro_batch_size": args.spectral_micro_batch_size,
    }
    for key, value in mapping.items():
        if value is not None:
            CONFIG[key] = value
    if args.seeds:
        CONFIG["training_seeds"] = [int(seed) for seed in args.seeds]
    if args.force_retrain:
        CONFIG["force_retrain"] = True
    if args.no_resume:
        CONFIG["resume_existing"] = False
    if args.models:
        unknown = [name for name in args.models if name not in ALL_EXPERIMENTS]
        if unknown:
            raise ValueError(f"Unknown model names: {unknown}")
        CONFIG["selected_models"] = list(args.models)
    else:
        CONFIG["selected_models"] = list(ALL_EXPERIMENTS)
    if args.smoke_test:
        CONFIG.update(
            {
                "resolution": min(int(CONFIG["resolution"]), 16),
                "batch_size": min(int(CONFIG["batch_size"]), 2),
                "epochs": 1,
                "pcfno_max_epochs": 1,
                "pcfno_minimum_epochs_before_early_stopping": 1,
                "pcfno_early_stopping_patience": 1,
                "uncertainty_refinement_epochs": 1,
                "uncertainty_refinement_patience": 1,
                "physics_warmup_epochs": 1,
                "width": 8,
                "modes": 4,
                "fno_layers": 1,
                "padding": 4,
                "time_samples_per_case": 2,
                "spectrum_time_samples": 8,
                "training_seeds": [42],
                "use_amp": False,
                "early_stopping_patience": 2,
                "max_train_batches": 1,
                "max_eval_batches": 1,
                "runtime_warmup": 1,
                "runtime_repeats": 2,
                "spectral_micro_batch_size": 1,
                "checkpoint_interval_epochs": 1,
                "save_jpg_dpi": 120,
                "force_retrain": True,
            }
        )
    resolve_derived_config()


def resolve_derived_config() -> None:
    profile_name = str(CONFIG.get("split_profile", "confirmation"))
    profiles = CONFIG.get("split_profiles", {})
    if profile_name not in profiles:
        raise ValueError(
            f"Unknown split_profile={profile_name!r}; available={sorted(profiles)}"
        )
    profile = profiles[profile_name]
    CONFIG["test_reynolds"] = [float(value) for value in profile["test_reynolds"]]
    CONFIG["validation_reynolds"] = [
        float(value) for value in profile["validation_reynolds"]
    ]
    CONFIG["excluded_reynolds"] = [
        float(value) for value in profile.get("excluded_reynolds", [])
    ]
    CONFIG["representative_re_targets"] = list(CONFIG["test_reynolds"])

    fraction = float(CONFIG.get("physics_warmup_fraction", 0.10))
    configured = CONFIG.get("physics_warmup_epochs")
    if configured is None:
        CONFIG["physics_warmup_epochs"] = max(
            1, int(round(float(CONFIG["epochs"]) * fraction))
        )
    else:
        CONFIG["physics_warmup_epochs"] = max(1, int(configured))
    CONFIG["pcfno_data_pretrain_fraction"] = 0.0


FAULT_LOG_HANDLE: Optional[Any] = None


def setup_native_diagnostics(output_dir: Path, logger: logging.Logger) -> None:
    global FAULT_LOG_HANDLE
    if bool(CONFIG.get("enable_native_faulthandler", True)):
        fault_path = output_dir / "native_crash_trace.log"
        FAULT_LOG_HANDLE = fault_path.open("a", buffering=1, encoding="utf-8")
        faulthandler.enable(file=FAULT_LOG_HANDLE, all_threads=True)
        atexit.register(lambda: FAULT_LOG_HANDLE and FAULT_LOG_HANDLE.flush())
        logger.info("Native faulthandler enabled: %s", fault_path)
    try:
        torch.set_num_threads(max(1, min(4, os.cpu_count() or 1)))
        torch.set_num_interop_threads(1)
    except RuntimeError:
        pass
    if torch.cuda.is_available() and str(CONFIG["device"]).startswith("cuda"):
        try:
            torch.backends.cuda.cufft_plan_cache.max_size = int(
                CONFIG.get("cufft_plan_cache_max_size", 16)
            )
            logger.info(
                "cuFFT plan cache max_size=%d",
                int(CONFIG.get("cufft_plan_cache_max_size", 16)),
            )
        except Exception as exc:
            logger.warning("Could not limit cuFFT plan cache: %s", exc)


def process_rss_mb() -> float:
    try:
        with open("/proc/self/status", "r", encoding="utf-8") as handle:
            for line in handle:
                if line.startswith("VmRSS:"):
                    return float(line.split()[1]) / 1024.0
    except OSError:
        pass
    return float("nan")


def resource_snapshot() -> Dict[str, float]:
    row = {"cpu_rss_mb": process_rss_mb()}
    if torch.cuda.is_available() and str(CONFIG["device"]).startswith("cuda"):
        device = torch.device(CONFIG["device"])
        row.update(
            {
                "cuda_allocated_mb": torch.cuda.memory_allocated(device) / 2**20,
                "cuda_reserved_mb": torch.cuda.memory_reserved(device) / 2**20,
                "cuda_max_allocated_mb": torch.cuda.max_memory_allocated(device)
                / 2**20,
            }
        )
    else:
        row.update(
            {
                "cuda_allocated_mb": 0.0,
                "cuda_reserved_mb": 0.0,
                "cuda_max_allocated_mb": 0.0,
            }
        )
    return row


def cuda_maintenance(epoch: int, logger: logging.Logger) -> None:
    interval = max(int(CONFIG.get("cuda_maintenance_interval", 0)), 0)
    if not torch.cuda.is_available() or not str(CONFIG["device"]).startswith("cuda"):
        gc.collect()
        return
    if interval <= 0 or (epoch + 1) % interval != 0:
        return
    try:
        torch.cuda.synchronize()
        before = resource_snapshot()
        torch.cuda.empty_cache()
        try:
            torch.backends.cuda.cufft_plan_cache.clear()
        except Exception:
            pass
        gc.collect()
        after = resource_snapshot()
        logger.info(
            "CUDA maintenance epoch %d | allocated %.1f->%.1f MB | reserved %.1f->%.1f MB | RSS %.1f MB",
            epoch + 1,
            before["cuda_allocated_mb"],
            after["cuda_allocated_mb"],
            before["cuda_reserved_mb"],
            after["cuda_reserved_mb"],
            after["cpu_rss_mb"],
        )
    except Exception as exc:
        logger.warning("CUDA maintenance failed at epoch %d: %s", epoch + 1, exc)


def tree_to_cpu(value: Any) -> Any:
    if isinstance(value, torch.Tensor):
        return value.detach().cpu()
    if isinstance(value, dict):
        return {key: tree_to_cpu(item) for key, item in value.items()}
    if isinstance(value, list):
        return [tree_to_cpu(item) for item in value]
    if isinstance(value, tuple):
        return tuple(tree_to_cpu(item) for item in value)
    return value


def safe_torch_load(path: Path, map_location: Any) -> Any:
    try:
        return torch.load(path, map_location=map_location, weights_only=False)
    except TypeError:
        return torch.load(path, map_location=map_location)


def atomic_torch_save(payload: Dict[str, Any], destination: Path) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary = destination.with_suffix(destination.suffix + ".tmp")
    torch.save(tree_to_cpu(payload), temporary)
    os.replace(temporary, destination)
    gc.collect()


def last_checkpoint_path(output_dir: Path, name: str, seed: int) -> Path:
    return output_dir / "checkpoints" / name / f"seed_{seed}_last.pt"


def rng_state_payload() -> Dict[str, Any]:
    payload: Dict[str, Any] = {
        "python_rng_state": random.getstate(),
        "numpy_rng_state": np.random.get_state(),
        "torch_rng_state": torch.get_rng_state(),
    }
    if torch.cuda.is_available():
        payload["cuda_rng_state_all"] = torch.cuda.get_rng_state_all()
    return payload


def restore_rng_state(checkpoint: Dict[str, Any]) -> None:
    try:
        if "python_rng_state" in checkpoint:
            random.setstate(checkpoint["python_rng_state"])
        if "numpy_rng_state" in checkpoint:
            np.random.set_state(checkpoint["numpy_rng_state"])
        if "torch_rng_state" in checkpoint:
            torch.set_rng_state(checkpoint["torch_rng_state"])
        if torch.cuda.is_available() and "cuda_rng_state_all" in checkpoint:
            torch.cuda.set_rng_state_all(checkpoint["cuda_rng_state_all"])
    except Exception:
        pass


def spectral_micro_batch_size(model: nn.Module, full_batch_size: int) -> int:
    if model_has_complex_parameters(model):
        configured = int(CONFIG.get("spectral_micro_batch_size", full_batch_size))
        return max(1, min(configured, full_batch_size))
    return full_batch_size


def set_seed(seed: int) -> None:
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def setup_logging(output_dir: Path) -> logging.Logger:
    output_dir.mkdir(parents=True, exist_ok=True)
    logger = logging.getLogger("PCFNO_SPEED")
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    formatter = logging.Formatter("%(asctime)s | %(levelname)s | %(message)s")
    stream_handler = logging.StreamHandler(sys.stdout)
    stream_handler.setFormatter(formatter)
    file_handler = logging.FileHandler(
        output_dir / "run.log", mode="w", encoding="utf-8"
    )
    file_handler.setFormatter(formatter)
    logger.addHandler(stream_handler)
    logger.addHandler(file_handler)
    return logger


def choose_font() -> str:
    for candidate in ["Arial", "Liberation Sans", "DejaVu Sans"]:
        try:
            fm.findfont(candidate, fallback_to_default=False)
            return candidate
        except ValueError:
            continue
    return "DejaVu Sans"


def configure_matplotlib() -> str:
    chosen_font = choose_font()
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "mathtext.fontset": "dejavusans",
            "font.size": 9.0,
            "axes.labelsize": 9.5,
            "axes.titlesize": 10.0,
            "legend.fontsize": 8.2,
            "xtick.labelsize": 8.5,
            "ytick.labelsize": 8.5,
            "axes.linewidth": 0.8,
            "lines.linewidth": 1.6,
            "lines.markersize": 4.0,
            "figure.dpi": 150,
            "savefig.dpi": int(CONFIG["save_jpg_dpi"]),
            "savefig.transparent": False,
            "axes.unicode_minus": False,
            "axes.grid": False,
            "legend.frameon": False,
        }
    )
    return chosen_font


def save_figure(fig: plt.Figure, base_path: Path) -> None:
    base_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        base_path.with_suffix(".jpg"),
        dpi=int(CONFIG["save_jpg_dpi"]),
        bbox_inches="tight",
        facecolor="white",
        pil_kwargs={"quality": 95, "subsampling": 0},
    )
    plt.close(fig)


def save_dat(
    path: Path,
    rows: Sequence[Dict[str, Any]],
    fieldnames: Optional[Sequence[str]] = None,
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    if fieldnames is None:
        ordered: List[str] = []
        for row in rows:
            for key in row.keys():
                if key not in ordered:
                    ordered.append(key)
        fieldnames = ordered
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(fieldnames),
            delimiter="\t",
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def read_dat(path: Path) -> List[Dict[str, Any]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    rows: List[Dict[str, Any]] = []
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for raw in reader:
            row: Dict[str, Any] = {}
            for key, value in raw.items():
                if value is None or value == "":
                    row[key] = value
                    continue
                try:
                    numeric = float(value)
                    row[key] = int(numeric) if numeric.is_integer() else numeric
                except ValueError:
                    row[key] = value
            rows.append(row)
    return rows


def save_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(data, handle, indent=2, ensure_ascii=False)


def cleanup_figure_formats(output_dir: Path) -> None:
    for extension in ["*.png", "*.pdf", "*.jpeg", "*.svg", "*.eps"]:
        for path in output_dir.rglob(extension):
            path.unlink(missing_ok=True)


def unzip_cfdbench(path: Path) -> None:
    import zipfile

    zip_info = {"bc.zip": "bc", "geo.zip": "geo", "prop.zip": "prop"}
    for zip_name, folder_name in zip_info.items():
        zip_path = path / zip_name
        extract_path = path / folder_name
        if extract_path.exists() or not zip_path.exists():
            continue
        extract_path.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zip_path, "r") as archive:
            archive.extractall(extract_path)


def safe_case_number(case_name: str) -> int:
    digits = "".join(ch for ch in case_name if ch.isdigit())
    return int(digits) if digits else 0


@dataclass(frozen=True)
class CaseInfo:
    name: str
    directory: Path
    num_frames: int
    height: int
    width: int
    params: Dict[str, float]
    reynolds: float
    inlet_velocity: float
    diameter: float


def load_case_params(case_dir: Path) -> Dict[str, float]:
    path = case_dir / "case.json"
    if not path.exists():
        raise FileNotFoundError(path)
    with path.open("r", encoding="utf-8") as handle:
        raw = json.load(handle)
    return {str(key): float(value) for key, value in raw.items()}


def discover_cases(data_root: Path, dataset_type: str) -> List[CaseInfo]:
    data_path = data_root / dataset_type
    if not data_path.exists():
        raise RuntimeError(f"Dataset directory not found: {data_path}")
    case_dirs = sorted(
        [p for p in data_path.iterdir() if p.is_dir() and p.name.startswith("case")],
        key=lambda p: safe_case_number(p.name),
    )
    cases: List[CaseInfo] = []
    for case_dir in case_dirs:
        speed_path = case_dir / "u.npy"
        if not speed_path.exists():
            raise FileNotFoundError(f"Missing first-channel u.npy in {case_dir}")
        shape = np.load(speed_path, mmap_mode="r").shape
        if len(shape) == 2:
            num_frames, height, width = 1, int(shape[0]), int(shape[1])
        elif len(shape) == 3:
            num_frames, height, width = int(shape[0]), int(shape[1]), int(shape[2])
        else:
            raise ValueError(f"Expected 2D or 3D u.npy in {case_dir}, got {shape}")
        params = load_case_params(case_dir)
        density = params["density"]
        viscosity = params["viscosity"]
        inlet_velocity = params["vel_in"]
        diameter = 2.0 * params["radius"]
        reynolds = density * inlet_velocity * diameter / max(viscosity, 1e-12)
        cases.append(
            CaseInfo(
                name=case_dir.name,
                directory=case_dir,
                num_frames=num_frames,
                height=height,
                width=width,
                params=params,
                reynolds=reynolds,
                inlet_velocity=inlet_velocity,
                diameter=diameter,
            )
        )
    if not cases:
        raise RuntimeError(f"No cases found in {data_path}")
    return cases


def canonical_re(value: float) -> float:
    return round(float(value), int(CONFIG.get("re_group_decimals", 8)))


def group_cases_by_re(cases: Sequence[CaseInfo]) -> Dict[float, List[CaseInfo]]:
    groups: Dict[float, List[CaseInfo]] = {}
    for case in cases:
        groups.setdefault(canonical_re(case.reynolds), []).append(case)
    for value in groups:
        groups[value] = sorted(
            groups[value], key=lambda item: safe_case_number(item.name)
        )
    return dict(sorted(groups.items()))


def exact_re_values(
    available: Sequence[float],
    targets: Sequence[float],
    excluded: Optional[set[float]] = None,
    label: str = "requested",
) -> List[float]:
    available_values = {canonical_re(value) for value in available}
    excluded_values = (
        set() if excluded is None else {canonical_re(value) for value in excluded}
    )
    requested = [canonical_re(value) for value in targets]
    if len(requested) != len(set(requested)):
        raise ValueError(f"Duplicate Reynolds groups in {label}: {requested}")
    missing = [value for value in requested if value not in available_values]
    if missing:
        raise ValueError(
            f"Missing Reynolds groups for {label}: {missing}; available={sorted(available_values)}"
        )
    blocked = [value for value in requested if value in excluded_values]
    if blocked:
        raise ValueError(
            f"Reynolds groups for {label} overlap excluded groups: {blocked}"
        )
    return requested


def nearest_unique_re_values(
    available: Sequence[float],
    targets: Sequence[float],
    excluded: Optional[set[float]] = None,
) -> List[float]:
    excluded = set() if excluded is None else set(excluded)
    chosen: List[float] = []
    for target in targets:
        candidates = [
            value
            for value in available
            if value not in excluded and value not in chosen
        ]
        if not candidates:
            break
        selected = min(
            candidates, key=lambda value: (abs(value - float(target)), value)
        )
        chosen.append(selected)
    return chosen


def representative_case_for_group(cases: Sequence[CaseInfo]) -> CaseInfo:
    ordered = sorted(cases, key=lambda item: safe_case_number(item.name))
    return ordered[len(ordered) // 2]


def nearest_unique_cases(
    cases: Sequence[CaseInfo], targets: Sequence[float]
) -> List[CaseInfo]:
    groups = group_cases_by_re(cases)
    if bool(CONFIG.get("strict_re_matching", True)):
        values = exact_re_values(list(groups), targets, label="representative")
    else:
        values = nearest_unique_re_values(list(groups), targets)
    return [representative_case_for_group(groups[value]) for value in values]


def split_cases(
    cases: Sequence[CaseInfo],
) -> Tuple[List[CaseInfo], List[CaseInfo], List[CaseInfo], List[CaseInfo]]:
    groups = group_cases_by_re(cases)
    all_re_values = list(groups)
    if len(all_re_values) < 8:
        raise ValueError("At least eight distinct Reynolds numbers are required")

    excluded_targets = CONFIG.get("excluded_reynolds") or []
    if bool(CONFIG.get("strict_re_matching", True)):
        excluded_re = (
            exact_re_values(all_re_values, excluded_targets, label="excluded")
            if excluded_targets
            else []
        )
    else:
        excluded_re = nearest_unique_re_values(all_re_values, excluded_targets)
    active_re_values = [
        value for value in all_re_values if value not in set(excluded_re)
    ]
    if len(active_re_values) < 8:
        raise ValueError(
            "The selected split profile leaves fewer than eight active Reynolds groups"
        )

    test_targets = CONFIG.get("test_reynolds", CONFIG["representative_re_targets"])
    if bool(CONFIG.get("strict_re_matching", True)):
        test_re = exact_re_values(active_re_values, test_targets, label="test")
    else:
        test_re = nearest_unique_re_values(active_re_values, test_targets)
    if len(test_re) < 3:
        raise ValueError("Could not select three distinct held-out Reynolds groups")

    validation_targets = CONFIG.get("validation_reynolds") or []
    if validation_targets and bool(CONFIG.get("strict_re_matching", True)):
        val_re = exact_re_values(
            active_re_values,
            validation_targets,
            excluded=set(test_re),
            label="validation",
        )
    else:
        val_re = nearest_unique_re_values(
            active_re_values, validation_targets, excluded=set(test_re)
        )
    if not val_re:
        remaining = [value for value in active_re_values if value not in test_re]
        count = max(
            1,
            int(round(len(active_re_values) * float(CONFIG.get("val_fraction", 0.1)))),
        )
        quantiles = np.linspace(0.15, 0.85, count)
        target_values = [float(np.quantile(remaining, q)) for q in quantiles]
        val_re = nearest_unique_re_values(remaining, target_values)

    train_re = [
        value for value in active_re_values if value not in set(test_re) | set(val_re)
    ]
    if not train_re or not val_re:
        raise RuntimeError(
            "Grouped-Re split produced an empty training or validation set"
        )

    train_cases = [case for value in train_re for case in groups[value]]
    val_cases = [case for value in val_re for case in groups[value]]
    test_cases = [case for value in test_re for case in groups[value]]
    representative = [representative_case_for_group(groups[value]) for value in test_re]

    train_set, val_set, test_set = set(train_re), set(val_re), set(test_re)
    if train_set & val_set or train_set & test_set or val_set & test_set:
        raise RuntimeError("Reynolds-number leakage detected after grouped split")
    if (train_set | val_set | test_set) & set(excluded_re):
        raise RuntimeError("Excluded Reynolds groups leaked into an active split")

    CONFIG["resolved_train_reynolds"] = train_re
    CONFIG["resolved_validation_reynolds"] = val_re
    CONFIG["resolved_test_reynolds"] = test_re
    CONFIG["resolved_excluded_reynolds"] = excluded_re
    CONFIG["representative_re_targets"] = test_re
    return train_cases, val_cases, test_cases, representative


def select_frame_indices(num_frames: int, count: Optional[int]) -> np.ndarray:
    if count is None or count <= 0 or count >= num_frames:
        return np.arange(num_frames, dtype=np.int64)
    return np.unique(np.linspace(0, num_frames - 1, int(count), dtype=np.int64))


def compute_time_scale(cases: Sequence[CaseInfo]) -> float:
    maximum = 0.0
    for case in cases:
        indices = select_frame_indices(case.num_frames, CONFIG["time_samples_per_case"])
        if indices.size == 0:
            continue
        final_time = float(indices[-1]) * float(CONFIG["data_delta_time"])
        time_star = final_time * case.inlet_velocity / max(case.diameter, 1e-12)
        maximum = max(maximum, time_star)
    return max(maximum, 1.0)


def resize_frame(array: np.ndarray, resolution: int) -> np.ndarray:
    if array.shape == (resolution, resolution):
        return np.asarray(array, dtype=np.float32)
    tensor = torch.from_numpy(np.array(array, dtype=np.float32, copy=True))[None, None]
    resized = F.interpolate(
        tensor, size=(resolution, resolution), mode="bilinear", align_corners=False
    )
    return resized[0, 0].numpy().astype(np.float32)


class CFDBenchSpeedDataset(Dataset):
    def __init__(
        self,
        cases: Sequence[CaseInfo],
        time_scale: float,
        resolution: int,
        time_samples_per_case: Optional[int],
        re_min: float,
        re_max: float,
    ):
        self.cases = list(cases)
        self.time_scale = float(time_scale)
        self.resolution = int(resolution)
        self.time_samples_per_case = time_samples_per_case
        self.re_min = float(re_min)
        self.re_max = float(re_max)
        self.samples: List[Tuple[int, int]] = []
        for case_index, case in enumerate(self.cases):
            for frame_index in select_frame_indices(
                case.num_frames, self.time_samples_per_case
            ):
                self.samples.append((case_index, int(frame_index)))
        if not self.samples:
            raise RuntimeError("Dataset contains no selected frames")
        self._geometry_cache: Dict[str, Tuple[np.ndarray, ...]] = {}
        self._speed_cache: Dict[str, np.ndarray] = {}

    def __len__(self) -> int:
        return len(self.samples)

    def _speed_array(self, case: CaseInfo) -> np.ndarray:
        if case.name not in self._speed_cache:
            self._speed_cache[case.name] = np.load(
                case.directory / "u.npy", mmap_mode="r"
            )
        return self._speed_cache[case.name]

    def _geometry(self, case: CaseInfo) -> Tuple[np.ndarray, ...]:
        if case.name in self._geometry_cache:
            return self._geometry_cache[case.name]
        params = case.params
        resolution = self.resolution
        x_min, x_max = params["x_min"], params["x_max"]
        y_min, y_max = params["y_min"], params["y_max"]
        dx = (x_max - x_min) / resolution
        dy = (y_max - y_min) / resolution
        x_coord = x_min + (np.arange(resolution, dtype=np.float32) + 0.5) * dx
        y_coord = y_min + (np.arange(resolution, dtype=np.float32) + 0.5) * dy
        x_grid, y_grid = np.meshgrid(x_coord, y_coord)
        center_x = params.get("center_x", 0.0)
        center_y = params.get("center_y", 0.0)
        radius = params["radius"]
        diameter = max(2.0 * radius, 1e-12)
        x_over_d = (x_grid - center_x) / diameter
        y_over_d = (y_grid - center_y) / diameter
        solid = (x_grid - center_x) ** 2 + (y_grid - center_y) ** 2 <= radius**2
        fluid = ~solid
        dist_fluid = distance_transform_edt(fluid, sampling=(dy, dx))
        dist_solid = distance_transform_edt(solid, sampling=(dy, dx))
        sdf_over_d = np.clip((dist_fluid - dist_solid) / diameter, -4.0, 4.0) / 4.0
        cylinder_wall_fluid = binary_dilation(solid) & (~solid)
        no_slip = cylinder_wall_fluid.copy()
        no_slip[0, :] = True
        no_slip[-1, :] = True
        boundary = no_slip.copy()
        boundary[:, 0] = True
        boundary[:, -1] = True
        geometry = (
            x_over_d.astype(np.float32),
            y_over_d.astype(np.float32),
            solid.astype(np.float32),
            sdf_over_d.astype(np.float32),
            boundary.astype(np.float32),
            no_slip.astype(np.float32),
        )
        self._geometry_cache[case.name] = geometry
        return geometry

    def __getitem__(
        self, index: int
    ) -> Tuple[torch.Tensor, torch.Tensor, torch.Tensor, Dict[str, Any]]:
        case_index, frame_index = self.samples[index]
        case = self.cases[case_index]
        speed_all = self._speed_array(case)
        if speed_all.ndim == 2:
            speed = np.asarray(speed_all, dtype=np.float32)
        else:
            speed = np.asarray(speed_all[frame_index], dtype=np.float32)
        speed = resize_frame(speed, self.resolution)
        speed = np.maximum(speed, 0.0)
        x_over_d, y_over_d, solid, sdf, boundary, no_slip = self._geometry(case)
        if bool(CONFIG.get("hard_zero_solid", True)):
            speed = speed.copy()
            speed[solid > 0.5] = 0.0
        physical_time = frame_index * float(CONFIG["data_delta_time"])
        time_star = physical_time * case.inlet_velocity / max(case.diameter, 1e-12)
        time_scaled = time_star / self.time_scale
        log_re = math.log10(max(case.reynolds, 1e-8))
        log_min = math.log10(max(self.re_min, 1e-8))
        log_max = math.log10(max(self.re_max, self.re_min + 1e-8))
        re_scaled = (log_re - log_min) / max(log_max - log_min, 1e-12)
        re_field = np.full_like(speed, re_scaled, dtype=np.float32)
        uin_field = np.full_like(speed, case.inlet_velocity, dtype=np.float32)
        time_field = np.full_like(speed, time_scaled, dtype=np.float32)
        inputs = np.stack(
            [
                x_over_d,
                y_over_d,
                solid,
                sdf,
                boundary,
                re_field,
                uin_field,
                time_field,
            ],
            axis=0,
        ).astype(np.float32)
        target = (speed / max(case.inlet_velocity, 1e-12))[None].astype(np.float32)
        metadata: Dict[str, Any] = {
            "case": case.name,
            "case_number": safe_case_number(case.name),
            "frame_index": frame_index,
            "num_frames": case.num_frames,
            "time": physical_time,
            "time_star": time_star,
            "time_scaled": time_scaled,
            "Re": case.reynolds,
            "Uin": case.inlet_velocity,
            "diameter": case.diameter,
            "radius": case.diameter / 2.0,
            "center_x": case.params.get("center_x", 0.0),
            "center_y": case.params.get("center_y", 0.0),
        }
        return (
            torch.from_numpy(inputs),
            torch.from_numpy(target),
            torch.from_numpy(no_slip[None]),
            metadata,
        )


def make_loader(dataset: Dataset, shuffle: bool, seed: int) -> DataLoader:
    generator = torch.Generator().manual_seed(int(seed))
    return DataLoader(
        dataset,
        batch_size=int(CONFIG["batch_size"]),
        shuffle=shuffle,
        generator=generator,
        num_workers=int(CONFIG["num_workers"]),
        pin_memory=bool(CONFIG["pin_memory"] and torch.cuda.is_available()),
        persistent_workers=bool(CONFIG["num_workers"] > 0),
    )


class SpectralConv2d(nn.Module):
    def __init__(self, in_channels: int, out_channels: int, modes1: int, modes2: int):
        super().__init__()
        self.in_channels = in_channels
        self.out_channels = out_channels
        self.modes1 = modes1
        self.modes2 = modes2
        scale = 1.0 / max(in_channels * out_channels, 1)
        self.weights1 = nn.Parameter(
            scale
            * torch.randn(in_channels, out_channels, modes1, modes2, dtype=torch.cfloat)
        )
        self.weights2 = nn.Parameter(
            scale
            * torch.randn(in_channels, out_channels, modes1, modes2, dtype=torch.cfloat)
        )

    @staticmethod
    def complex_multiply(inputs: torch.Tensor, weights: torch.Tensor) -> torch.Tensor:
        return torch.einsum("bixy,ioxy->boxy", inputs, weights)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        batch_size, _, height, width = x.shape
        max_modes1 = min(self.modes1, height // 2)
        max_modes2 = min(self.modes2, width // 2 + 1)
        output_dtype = x.dtype
        with torch.autocast(device_type=x.device.type, enabled=False):
            x_float = x.float()
            x_ft = torch.fft.rfft2(x_float)
            out_ft = torch.zeros(
                batch_size,
                self.out_channels,
                height,
                width // 2 + 1,
                dtype=torch.cfloat,
                device=x.device,
            )
            out_ft[:, :, :max_modes1, :max_modes2] = self.complex_multiply(
                x_ft[:, :, :max_modes1, :max_modes2],
                self.weights1[:, :, :max_modes1, :max_modes2],
            )
            out_ft[:, :, -max_modes1:, :max_modes2] = self.complex_multiply(
                x_ft[:, :, -max_modes1:, :max_modes2],
                self.weights2[:, :, :max_modes1, :max_modes2],
            )
            result = torch.fft.irfft2(out_ft, s=(height, width))
        return result.to(dtype=output_dtype)


class FNOBlock(nn.Module):
    def __init__(self, width: int, modes: int):
        super().__init__()
        groups = 8 if width % 8 == 0 else 1
        self.spectral = SpectralConv2d(width, width, modes, modes)
        self.pointwise = nn.Conv2d(width, width, 1)
        self.norm = nn.GroupNorm(groups, width)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return F.gelu(self.norm(x + self.spectral(x) + self.pointwise(x)))


class FNOCore(nn.Module):
    def __init__(self, modes: int, width: int, layers: int):
        super().__init__()
        self.width = width
        self.lift = nn.Conv2d(INPUT_CHANNELS, width, 1)
        self.blocks = nn.ModuleList([FNOBlock(width, modes) for _ in range(layers)])
        self.decoder = nn.Sequential(
            nn.Conv2d(width, 128, 1),
            nn.GELU(),
            nn.Conv2d(128, OUTPUT_CHANNELS, 1),
        )

    def encode(self, inputs: torch.Tensor) -> torch.Tensor:
        x = self.lift(inputs)
        padding = int(CONFIG["padding"])
        if padding > 0:
            x = F.pad(x, [0, padding, 0, padding], mode="replicate")
        for block in self.blocks:
            x = block(x)
        if padding > 0:
            x = x[:, :, :-padding, :-padding]
        return x

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        return self.decoder(self.encode(inputs))


class FNO2d(nn.Module):
    def __init__(self, modes: int, width: int, layers: int):
        super().__init__()
        self.core = FNOCore(modes, width, layers)

    def forward(self, inputs: torch.Tensor) -> torch.Tensor:
        return self.core(inputs)


class DoubleConv(nn.Module):
    def __init__(self, in_channels: int, out_channels: int):
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv2d(in_channels, out_channels, 3, padding=1),
            nn.BatchNorm2d(out_channels),
            nn.GELU(),
            nn.Conv2d(out_channels, out_channels, 3, padding=1),
            nn.BatchNorm2d(out_channels),
            nn.GELU(),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.net(x)


class UNet(nn.Module):
    def __init__(self):
        super().__init__()
        self.down1 = DoubleConv(INPUT_CHANNELS, 64)
        self.pool1 = nn.MaxPool2d(2)
        self.down2 = DoubleConv(64, 128)
        self.pool2 = nn.MaxPool2d(2)
        self.middle = DoubleConv(128, 256)
        self.up2 = nn.ConvTranspose2d(256, 128, 2, stride=2)
        self.conv2 = DoubleConv(256, 128)
        self.up1 = nn.ConvTranspose2d(128, 64, 2, stride=2)
        self.conv1 = DoubleConv(128, 64)
        self.output = nn.Conv2d(64, OUTPUT_CHANNELS, 1)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        x1 = self.down1(x)
        x2 = self.down2(self.pool1(x1))
        x3 = self.middle(self.pool2(x2))
        x = self.up2(x3)
        if x.shape[-2:] != x2.shape[-2:]:
            x = F.interpolate(
                x, size=x2.shape[-2:], mode="bilinear", align_corners=False
            )
        x = self.conv2(torch.cat([x, x2], dim=1))
        x = self.up1(x)
        if x.shape[-2:] != x1.shape[-2:]:
            x = F.interpolate(
                x, size=x1.shape[-2:], mode="bilinear", align_corners=False
            )
        return self.output(self.conv1(torch.cat([x, x1], dim=1)))


class DeepONetBranch(nn.Module):
    def __init__(self, hidden_dim: int = 128):
        super().__init__()
        self.encoder = nn.Sequential(
            nn.Conv2d(INPUT_CHANNELS, 32, 3, padding=1),
            nn.GELU(),
            nn.Conv2d(32, 64, 3, stride=2, padding=1),
            nn.GELU(),
            nn.AdaptiveAvgPool2d((8, 8)),
            nn.Flatten(),
            nn.Linear(64 * 8 * 8, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        return self.encoder(x)


class DeepONetTrunk(nn.Module):
    def __init__(self, hidden_dim: int = 128):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(2, hidden_dim),
            nn.GELU(),
            nn.Linear(hidden_dim, hidden_dim),
        )

    def forward(self, coordinates: torch.Tensor) -> torch.Tensor:
        return self.net(coordinates)


class DeepONet(nn.Module):
    def __init__(self, hidden_dim: int = 128):
        super().__init__()
        self.branch = DeepONetBranch(hidden_dim)
        self.trunk = DeepONetTrunk(hidden_dim)
        self.output = nn.Linear(hidden_dim, OUTPUT_CHANNELS)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        batch, _, height, width = x.shape
        branch_feature = self.branch(x).unsqueeze(1)
        coordinates = x[:, :2].permute(0, 2, 3, 1).reshape(batch, height * width, 2)
        trunk_feature = self.trunk(coordinates)
        features = branch_feature * trunk_feature
        outputs = self.output(features).reshape(batch, height, width, OUTPUT_CHANNELS)
        return outputs.permute(0, 3, 1, 2)


class SDFEncoder(nn.Module):
    def __init__(self, out_channels: int):
        super().__init__()
        hidden = max(out_channels // 2, 8)
        self.net = nn.Sequential(
            nn.Conv2d(1, hidden, 3, padding=1),
            nn.GELU(),
            nn.Conv2d(hidden, out_channels, 3, padding=1),
            nn.GELU(),
        )

    def forward(self, sdf: torch.Tensor) -> torch.Tensor:
        return self.net(sdf)


class PhysicsGuidedResidualAdapter(nn.Module):
    def __init__(self, backbone_width: int, adapter_width: int):
        super().__init__()
        self.input_projection = nn.Conv2d(
            backbone_width + INPUT_CHANNELS + adapter_width,
            adapter_width,
            1,
        )
        groups = 8 if adapter_width % 8 == 0 else 1
        self.local_path = nn.Sequential(
            nn.Conv2d(
                adapter_width,
                adapter_width,
                3,
                padding=1,
                groups=adapter_width,
            ),
            nn.GELU(),
            nn.Conv2d(adapter_width, adapter_width, 1),
            nn.GroupNorm(groups, adapter_width),
            nn.GELU(),
        )
        self.output_head = nn.Conv2d(adapter_width, OUTPUT_CHANNELS, 1)
        nn.init.zeros_(self.output_head.weight)
        nn.init.zeros_(self.output_head.bias)

    def forward(
        self,
        backbone_feature: torch.Tensor,
        adapter_inputs: torch.Tensor,
        sdf_feature: torch.Tensor,
    ) -> Tuple[torch.Tensor, torch.Tensor]:
        feature = self.input_projection(
            torch.cat(
                [backbone_feature.detach(), adapter_inputs, sdf_feature],
                dim=1,
            )
        )
        feature = feature + self.local_path(feature)
        residual = self.output_head(feature)
        return residual, feature


class PCFNO(nn.Module):

    def __init__(
        self,
        modes: int,
        width: int,
        layers: int,
        use_sdf_adapter: bool = True,
        use_physics_focus: bool = True,
        residual_alpha: float = 0.5,
        uncertainty: bool = True,
    ):
        super().__init__()
        self.core = FNOCore(modes, width, layers)
        self.use_sdf_adapter = bool(use_sdf_adapter)
        self.use_physics_focus = bool(use_physics_focus)
        self.residual_alpha = float(residual_alpha)
        self.uncertainty = bool(uncertainty)
        self.correction_enabled = True
        self.backbone_is_frozen = False

        adapter_width = int(CONFIG.get("adapter_width", 32))
        self.sdf_encoder = SDFEncoder(adapter_width)
        self.adapter = PhysicsGuidedResidualAdapter(
            backbone_width=width,
            adapter_width=adapter_width,
        )
        self.uncertainty_head = nn.Sequential(
            nn.Conv2d(adapter_width, adapter_width, 3, padding=1),
            nn.GELU(),
            nn.Conv2d(adapter_width, OUTPUT_CHANNELS, 1),
        )
        self.register_buffer(
            "inference_alpha",
            torch.tensor(float(self.residual_alpha)),
        )
        self.register_buffer(
            "warm_started_from_fno", torch.tensor(0, dtype=torch.uint8)
        )

    def load_fno_backbone(self, fno_state_dict: Dict[str, torch.Tensor]) -> None:
        core_state = {
            key[len("core.") :]: value
            for key, value in fno_state_dict.items()
            if key.startswith("core.")
        }
        if not core_state:
            raise RuntimeError("The supplied FNO checkpoint has no core.* parameters")
        self.core.load_state_dict(core_state, strict=True)
        self.warm_started_from_fno.fill_(1)

    def freeze_backbone(self) -> None:
        for parameter in self.core.parameters():
            parameter.requires_grad_(False)
        self.backbone_is_frozen = True

    def physics_focus(self, inputs: torch.Tensor) -> torch.Tensor:
        if not self.use_physics_focus:
            return torch.ones_like(inputs[:, :1])

        sdf_over_d = torch.abs(inputs[:, INPUT_SDF : INPUT_SDF + 1]) * 4.0
        boundary_scale = max(
            float(CONFIG.get("adapter_boundary_scale_over_d", 0.75)), 1e-6
        )
        boundary_focus = torch.exp(-sdf_over_d / boundary_scale)

        x = inputs[:, INPUT_XD : INPUT_XD + 1]
        y = inputs[:, INPUT_YD : INPUT_YD + 1]
        x_smooth = max(float(CONFIG.get("adapter_wake_x_smoothing", 0.30)), 1e-6)
        y_smooth = max(float(CONFIG.get("adapter_wake_y_smoothing", 0.30)), 1e-6)
        wake_left = torch.sigmoid((x - float(CONFIG["wake_x_min_over_d"])) / x_smooth)
        wake_right = torch.sigmoid((float(CONFIG["wake_x_max_over_d"]) - x) / x_smooth)
        wake_vertical = torch.sigmoid(
            (float(CONFIG["wake_half_height_over_d"]) - torch.abs(y)) / y_smooth
        )
        time_activity = 0.20 + 0.80 * torch.sigmoid(
            (inputs[:, INPUT_TIME : INPUT_TIME + 1] - 0.08) / 0.04
        )
        re_activity = 0.25 + 0.75 * torch.sigmoid(
            (inputs[:, INPUT_RE : INPUT_RE + 1] - 0.25) / 0.12
        )
        wake_focus = (
            wake_left * wake_right * wake_vertical * time_activity * re_activity
        )

        floor = float(CONFIG.get("adapter_focus_floor", 0.20))
        combined = torch.maximum(boundary_focus, wake_focus).clamp(0.0, 1.0)
        return floor + (1.0 - floor) * combined

    def adapter_inputs(self, inputs: torch.Tensor) -> Tuple[torch.Tensor, torch.Tensor]:
        adapter_width = int(CONFIG.get("adapter_width", 32))
        if not self.use_sdf_adapter:
            result = inputs.clone()
            result[:, INPUT_SDF] = 0.0
            sdf_feature = torch.zeros(
                inputs.shape[0],
                adapter_width,
                inputs.shape[-2],
                inputs.shape[-1],
                dtype=inputs.dtype,
                device=inputs.device,
            )
            return result, sdf_feature

        sdf_feature = self.sdf_encoder(inputs[:, INPUT_SDF : INPUT_SDF + 1])
        return inputs, sdf_feature

    def positive_components(
        self,
        inputs: torch.Tensor,
        alpha_override: Optional[float] = None,
        return_uncertainty: bool = False,
    ) -> Tuple[torch.Tensor, Optional[torch.Tensor], Dict[str, torch.Tensor]]:
        if self.backbone_is_frozen:
            with torch.no_grad():
                backbone_feature = self.core.encode(inputs)
                base_raw = self.core.decoder(backbone_feature)
        else:
            backbone_feature = self.core.encode(inputs)
            base_raw = self.core.decoder(backbone_feature)
        base_positive = positive_speed(base_raw)

        if not self.correction_enabled:
            components = {
                "base_positive": base_positive,
                "corrected_positive": base_positive,
                "positive_residual": torch.zeros_like(base_positive),
                "raw_residual": torch.zeros_like(base_positive),
                "focus": torch.ones_like(base_positive),
            }
            return base_positive, None, components

        adapter_inputs, sdf_feature = self.adapter_inputs(inputs)
        raw_residual, adapter_feature = self.adapter(
            backbone_feature, adapter_inputs, sdf_feature
        )
        focus = self.physics_focus(inputs)
        correction_limit = float(CONFIG.get("adapter_correction_limit", 0.25))
        bounded_residual = correction_limit * torch.tanh(raw_residual) * focus
        corrected_positive = positive_speed(base_raw + bounded_residual)
        positive_residual = corrected_positive - base_positive

        if alpha_override is None:
            alpha_value = float(self.inference_alpha.item())
        else:
            alpha_value = float(alpha_override)
        alpha_value = min(max(alpha_value, 0.0), 1.0)
        positive_prediction = base_positive + alpha_value * positive_residual

        sigma: Optional[torch.Tensor]
        if return_uncertainty and self.uncertainty:
            sigma = F.softplus(self.uncertainty_head(adapter_feature.detach())) + 1e-5
        else:
            sigma = None

        components = {
            "base_positive": base_positive,
            "corrected_positive": corrected_positive,
            "positive_residual": positive_residual,
            "raw_residual": bounded_residual,
            "focus": focus,
        }
        return positive_prediction, sigma, components

    def forward(
        self,
        inputs: torch.Tensor,
        return_uncertainty: bool = False,
        alpha_override: Optional[float] = None,
    ) -> Any:
        positive_prediction, sigma, _ = self.positive_components(
            inputs,
            alpha_override=alpha_override,
            return_uncertainty=return_uncertainty,
        )
        if return_uncertainty and self.uncertainty:
            return positive_prediction, sigma
        return positive_prediction


def model_factory(name: str) -> nn.Module:
    if name == "UNet":
        return UNet()
    if name == "DeepONet":
        return DeepONet()
    if name == "FNO":
        return FNO2d(CONFIG["modes"], CONFIG["width"], CONFIG["fno_layers"])

    common = {
        "modes": CONFIG["modes"],
        "width": CONFIG["width"],
        "layers": CONFIG["fno_layers"],
    }
    flags = {
        "use_sdf_adapter": name != "PCFNO_without_SDF_adapter",
        "use_physics_focus": name != "PCFNO_without_physics_focus",
        "residual_alpha": (
            float(CONFIG.get("unblended_residual_alpha", 1.0))
            if name == "PCFNO_without_residual_blending"
            else float(CONFIG.get("fixed_residual_alpha", 0.5))
        ),
        "uncertainty": True,
    }
    if name.startswith("PCFNO"):
        return PCFNO(**common, **flags)
    raise KeyError(f"Unknown experiment: {name}")


def loss_configuration(name: str) -> Dict[str, float]:
    is_pcfno = name.startswith("PCFNO")
    result = {
        "lambda_relative": float(CONFIG["lambda_relative"]),
        "lambda_gradient": float(CONFIG["lambda_gradient"]) if is_pcfno else 0.0,
        "lambda_boundary": float(CONFIG["lambda_boundary"]) if is_pcfno else 0.0,
        "lambda_wake": 0.0,
        "lambda_non_degradation": (
            float(CONFIG["lambda_non_degradation"]) if is_pcfno else 0.0
        ),
        "lambda_residual_regularization": (
            float(CONFIG["lambda_residual_regularization"]) if is_pcfno else 0.0
        ),
        "lambda_uncertainty": (
            float(CONFIG["lambda_uncertainty"]) if is_pcfno else 0.0
        ),
        "lambda_calibration": (
            float(CONFIG["lambda_calibration"]) if is_pcfno else 0.0
        ),
    }
    if name == "PCFNO_without_physics_losses":
        result["lambda_gradient"] = 0.0
        result["lambda_boundary"] = 0.0
    return result


def fluid_mask(inputs: torch.Tensor) -> torch.Tensor:
    return inputs[:, INPUT_SOLID : INPUT_SOLID + 1] < 0.5


def wake_mask(inputs: torch.Tensor) -> torch.Tensor:
    fluid = fluid_mask(inputs)
    x = inputs[:, INPUT_XD : INPUT_XD + 1]
    y = inputs[:, INPUT_YD : INPUT_YD + 1]
    return (
        fluid
        & (x >= float(CONFIG["wake_x_min_over_d"]))
        & (x <= float(CONFIG["wake_x_max_over_d"]))
        & (torch.abs(y) <= float(CONFIG["wake_half_height_over_d"]))
    )


def valid_gradient_mask(inputs: torch.Tensor) -> torch.Tensor:
    solid = inputs[:, INPUT_SOLID : INPUT_SOLID + 1] > 0.5
    pixels = max(int(CONFIG["gradient_exclusion_pixels"]), 0)
    invalid = solid.float()
    for _ in range(pixels):
        invalid = F.max_pool2d(invalid, kernel_size=3, stride=1, padding=1)
    valid = invalid < 0.5
    valid[:, :, 0, :] = False
    valid[:, :, -1, :] = False
    valid[:, :, :, 0] = False
    valid[:, :, :, -1] = False
    return valid


def spatial_gradients(
    field: torch.Tensor, inputs: torch.Tensor
) -> Tuple[torch.Tensor, torch.Tensor]:
    height, width = field.shape[-2:]
    x_coord = inputs[:, INPUT_XD, 0, :]
    y_coord = inputs[:, INPUT_YD, :, 0]
    dx = (
        ((x_coord[:, -1] - x_coord[:, 0]) / max(width - 1, 1))
        .abs()
        .clamp_min(1e-12)
        .view(-1, 1, 1)
    )
    dy = (
        ((y_coord[:, -1] - y_coord[:, 0]) / max(height - 1, 1))
        .abs()
        .clamp_min(1e-12)
        .view(-1, 1, 1)
    )
    df_dy, df_dx = torch.gradient(field, dim=(1, 2))
    return df_dx / dx, df_dy / dy


def supervised_terms(
    prediction: torch.Tensor,
    target: torch.Tensor,
    inputs: torch.Tensor,
    apply_training_time_weight: bool = False,
) -> Tuple[torch.Tensor, torch.Tensor]:
    mask = fluid_mask(inputs)
    sample_mse: List[torch.Tensor] = []
    sample_relative: List[torch.Tensor] = []
    sample_weights: List[torch.Tensor] = []
    threshold = float(CONFIG.get("initial_time_threshold_scaled", 0.02))
    initial_weight = float(CONFIG.get("initial_sample_weight", 1.0))

    for index in range(prediction.shape[0]):
        valid = mask[index, 0]
        squared_error = (prediction[index, 0][valid] - target[index, 0][valid]).square()
        sample_mse.append(torch.mean(squared_error))
        numerator = torch.sum(squared_error)
        denominator = torch.sum(target[index, 0][valid] ** 2).clamp_min(1e-12)
        sample_relative.append(torch.sqrt(numerator / denominator))

        if apply_training_time_weight:
            time_value = inputs[index, INPUT_TIME, 0, 0]
            weight = torch.where(
                time_value <= threshold,
                time_value.new_tensor(initial_weight),
                time_value.new_tensor(1.0),
            )
        else:
            weight = prediction.new_tensor(1.0)
        sample_weights.append(weight)

    weights = torch.stack(sample_weights)
    weights = weights / weights.sum().clamp_min(1e-12)
    mse = torch.sum(weights * torch.stack(sample_mse))
    relative = torch.sum(weights * torch.stack(sample_relative))
    return mse, relative


def wake_matching_loss(
    prediction: torch.Tensor, target: torch.Tensor, inputs: torch.Tensor
) -> torch.Tensor:
    valid = wake_mask(inputs)
    if not torch.any(valid):
        return prediction.new_zeros(())
    error = prediction[valid] - target[valid]
    scale = torch.mean(target[valid] ** 2).detach().clamp_min(1e-4)
    return torch.mean(error**2) / scale


def gradient_matching_loss(
    prediction: torch.Tensor, target: torch.Tensor, inputs: torch.Tensor
) -> torch.Tensor:
    pred_dx, pred_dy = spatial_gradients(prediction[:, 0], inputs)
    true_dx, true_dy = spatial_gradients(target[:, 0], inputs)
    valid = valid_gradient_mask(inputs)[:, 0]
    if not torch.any(valid):
        return prediction.new_zeros(())
    return torch.mean(
        (pred_dx[valid] - true_dx[valid]) ** 2 + (pred_dy[valid] - true_dy[valid]) ** 2
    )


def conservative_non_degradation_loss(
    prediction: torch.Tensor,
    base_prediction: torch.Tensor,
    target: torch.Tensor,
    inputs: torch.Tensor,
) -> torch.Tensor:
    mask = fluid_mask(inputs)
    penalties: List[torch.Tensor] = []
    for index in range(prediction.shape[0]):
        valid = mask[index, 0]
        pred_mse = torch.mean(
            (prediction[index, 0][valid] - target[index, 0][valid]).square()
        )
        with torch.no_grad():
            base_mse = torch.mean(
                (base_prediction[index, 0][valid] - target[index, 0][valid]).square()
            )
        penalties.append(F.relu(pred_mse - base_mse))
    if not penalties:
        return prediction.new_zeros(())
    return torch.mean(torch.stack(penalties))


def residual_regularization_loss(
    positive_residual: torch.Tensor, inputs: torch.Tensor
) -> torch.Tensor:
    valid = fluid_mask(inputs)
    if not torch.any(valid):
        return positive_residual.new_zeros(())
    return torch.mean(positive_residual[valid].square())


def positive_speed(raw_prediction: torch.Tensor) -> torch.Tensor:
    if not bool(CONFIG.get("enforce_positive_output", True)):
        return raw_prediction
    mode = str(CONFIG.get("positive_output_mode", "smooth_relu")).lower()
    if mode == "softplus":
        beta = float(CONFIG.get("positive_output_beta", 6.0))
        return F.softplus(raw_prediction, beta=beta)
    if mode != "smooth_relu":
        raise ValueError(f"Unsupported positive_output_mode={mode!r}")
    epsilon = float(CONFIG.get("positive_output_epsilon", 1e-4))
    return 0.5 * (
        raw_prediction + torch.sqrt(raw_prediction.square() + epsilon * epsilon)
    )


def immersed_boundary_values(
    positive_prediction: torch.Tensor, inputs: torch.Tensor
) -> torch.Tensor:
    phi = inputs[:, INPUT_SDF : INPUT_SDF + 1]
    solid = inputs[:, INPUT_SOLID : INPUT_SOLID + 1] > 0.5
    values: List[torch.Tensor] = []

    def edge_values(
        p1: torch.Tensor,
        p2: torch.Tensor,
        phi1: torch.Tensor,
        phi2: torch.Tensor,
        s1: torch.Tensor,
        s2: torch.Tensor,
    ) -> None:
        crossing = s1 != s2
        if not torch.any(crossing):
            return
        a1 = torch.abs(phi1)
        a2 = torch.abs(phi2)
        interpolated = (a2 * p1 + a1 * p2) / (a1 + a2).clamp_min(1e-8)
        values.append(interpolated[crossing])

    edge_values(
        positive_prediction[:, :, :, :-1],
        positive_prediction[:, :, :, 1:],
        phi[:, :, :, :-1],
        phi[:, :, :, 1:],
        solid[:, :, :, :-1],
        solid[:, :, :, 1:],
    )
    edge_values(
        positive_prediction[:, :, :-1, :],
        positive_prediction[:, :, 1:, :],
        phi[:, :, :-1, :],
        phi[:, :, 1:, :],
        solid[:, :, :-1, :],
        solid[:, :, 1:, :],
    )
    if not values:
        return positive_prediction.new_zeros((0,))
    return torch.cat(values)


def immersed_boundary_loss(
    positive_prediction: torch.Tensor, inputs: torch.Tensor
) -> torch.Tensor:
    boundary_values = immersed_boundary_values(positive_prediction, inputs)
    if boundary_values.numel() == 0:
        return positive_prediction.new_zeros(())
    return torch.mean(boundary_values**2)


def uncertainty_regression_loss(
    sigma: torch.Tensor,
    prediction: torch.Tensor,
    target: torch.Tensor,
    inputs: torch.Tensor,
) -> torch.Tensor:
    mask = fluid_mask(inputs)
    absolute_error = torch.abs(prediction - target).detach()
    return torch.mean((sigma[mask] - absolute_error[mask]) ** 2)


def calibration_loss(
    sigma: torch.Tensor,
    prediction: torch.Tensor,
    target: torch.Tensor,
    inputs: torch.Tensor,
) -> torch.Tensor:
    mask = fluid_mask(inputs)
    absolute_error = torch.abs(prediction - target).detach().clamp_min(1e-5)
    sigma_valid = sigma.clamp_min(1e-5)
    return torch.mean(
        (torch.log(sigma_valid[mask]) - torch.log(absolute_error[mask])) ** 2
    )


def apply_output_constraints(
    positive_prediction: torch.Tensor, inputs: torch.Tensor
) -> torch.Tensor:
    if not bool(CONFIG.get("hard_zero_solid", True)):
        return positive_prediction
    fluid = (inputs[:, INPUT_SOLID : INPUT_SOLID + 1] < 0.5).to(
        positive_prediction.dtype
    )
    return positive_prediction * fluid


def forward_raw_with_sigma(
    model: nn.Module, inputs: torch.Tensor
) -> Tuple[torch.Tensor, Optional[torch.Tensor]]:
    if isinstance(model, PCFNO):
        positive_raw, sigma, _ = model.positive_components(
            inputs, return_uncertainty=model.uncertainty
        )
        return positive_raw, sigma
    return model(inputs), None


def forward_with_sigma(
    model: nn.Module,
    inputs: torch.Tensor,
    return_positive_raw: bool = False,
    return_components: bool = False,
    alpha_override: Optional[float] = None,
) -> Any:
    components: Optional[Dict[str, torch.Tensor]] = None

    if isinstance(model, PCFNO):
        positive_raw, sigma, components = model.positive_components(
            inputs,
            alpha_override=alpha_override,
            return_uncertainty=model.uncertainty,
        )
    else:
        raw_prediction = model(inputs)
        positive_raw = positive_speed(raw_prediction)
        sigma = None

    prediction = apply_output_constraints(positive_raw, inputs)
    if sigma is not None and bool(CONFIG.get("hard_zero_solid", True)):
        fluid = (inputs[:, INPUT_SOLID : INPUT_SOLID + 1] < 0.5).to(sigma.dtype)
        sigma = sigma * fluid

    if return_components and return_positive_raw:
        return prediction, sigma, positive_raw, components
    if return_components:
        return prediction, sigma, components
    if return_positive_raw:
        return prediction, sigma, positive_raw
    return prediction, sigma


def model_has_complex_parameters(model: nn.Module) -> bool:
    return any(torch.is_complex(parameter) for parameter in model.parameters())


def model_uses_amp(model: nn.Module) -> bool:
    return bool(
        CONFIG["use_amp"]
        and str(CONFIG["device"]).startswith("cuda")
        and torch.cuda.is_available()
        and not model_has_complex_parameters(model)
    )


def model_amp_dtype(model: nn.Module) -> Optional[torch.dtype]:
    if not model_uses_amp(model):
        return None
    requested = str(CONFIG.get("amp_dtype", "auto")).lower()
    if requested == "float32":
        return None
    if requested in {"bfloat16", "bf16"}:
        if hasattr(torch.cuda, "is_bf16_supported") and torch.cuda.is_bf16_supported():
            return torch.bfloat16
        return torch.float16
    if requested in {"float16", "fp16"}:
        return torch.float16
    if hasattr(torch.cuda, "is_bf16_supported") and torch.cuda.is_bf16_supported():
        return torch.bfloat16
    return torch.float16


def model_autocast_context(model: nn.Module) -> Any:
    dtype = model_amp_dtype(model)
    return (
        torch.autocast(device_type="cuda", dtype=dtype)
        if dtype is not None
        else nullcontext()
    )


def create_grad_scaler(model: nn.Module) -> Any:
    dtype = model_amp_dtype(model)
    enabled = dtype == torch.float16
    kwargs = {
        "enabled": enabled,
        "init_scale": float(CONFIG["amp_init_scale"]),
        "growth_factor": float(CONFIG["amp_growth_factor"]),
        "backoff_factor": float(CONFIG["amp_backoff_factor"]),
        "growth_interval": int(CONFIG["amp_growth_interval"]),
    }
    try:
        return torch.amp.GradScaler("cuda", **kwargs)
    except (AttributeError, TypeError):
        return torch.cuda.amp.GradScaler(**kwargs)


def precision_mode_name(model: nn.Module) -> str:
    dtype = model_amp_dtype(model)
    if dtype == torch.bfloat16:
        return "BF16 autocast"
    if dtype == torch.float16:
        return "FP16 autocast with GradScaler"
    return "FP32"


def smooth_curriculum(epoch: int, warmup_epochs: int) -> float:
    if warmup_epochs <= 0:
        return 1.0
    progress = min(max(float(epoch + 1) / float(warmup_epochs), 0.0), 1.0)
    return 0.5 - 0.5 * math.cos(math.pi * progress)


def train_epoch(
    model: nn.Module,
    loader: DataLoader,
    optimizer: torch.optim.Optimizer,
    scaler: Any,
    epoch: int,
    loss_cfg: Dict[str, float],
) -> Dict[str, float]:
    model.train()
    if isinstance(model, PCFNO) and bool(CONFIG.get("pcfno_freeze_backbone", True)):
        model.core.eval()

    keys = [
        "total",
        "data_mse",
        "relative",
        "gradient",
        "boundary",
        "non_degradation",
        "residual_regularization",
        "uncertainty",
        "calibration",
        "weighted_relative",
        "weighted_gradient",
        "weighted_boundary",
        "weighted_non_degradation",
        "weighted_residual_regularization",
        "weighted_uncertainty",
        "weighted_calibration",
    ]
    totals = {key: 0.0 for key in keys}
    successful_batches = 0
    attempted_batches = 0
    nonfinite_batches = 0

    if isinstance(model, PCFNO):
        model.correction_enabled = True
        warmup = smooth_curriculum(epoch, int(CONFIG["physics_warmup_epochs"]))
    else:
        warmup = 0.0

    for batch_index, (inputs_cpu, target_cpu, _, _) in enumerate(loader):
        if CONFIG["max_train_batches"] is not None and batch_index >= int(
            CONFIG["max_train_batches"]
        ):
            break
        attempted_batches += 1
        full_batch_size = int(inputs_cpu.shape[0])
        micro_size = spectral_micro_batch_size(model, full_batch_size)
        optimizer.zero_grad(set_to_none=True)
        batch_values = {key: 0.0 for key in keys}
        batch_failed = False

        for start_index in range(0, full_batch_size, micro_size):
            stop_index = min(start_index + micro_size, full_batch_size)
            fraction = float(stop_index - start_index) / float(full_batch_size)
            inputs = inputs_cpu[start_index:stop_index].to(
                CONFIG["device"], non_blocking=True
            )
            target = target_cpu[start_index:stop_index].to(
                CONFIG["device"], non_blocking=True
            )

            with model_autocast_context(model):
                if isinstance(model, PCFNO):
                    prediction, sigma, positive_raw, components = forward_with_sigma(
                        model,
                        inputs,
                        return_positive_raw=True,
                        return_components=True,
                    )
                    assert components is not None
                    base_prediction = apply_output_constraints(
                        components["base_positive"], inputs
                    )
                    positive_residual = components["positive_residual"]
                else:
                    prediction, sigma, positive_raw = forward_with_sigma(
                        model, inputs, return_positive_raw=True
                    )
                    base_prediction = prediction.detach()
                    positive_residual = torch.zeros_like(prediction)

                data_mse, relative = supervised_terms(
                    prediction, target, inputs, apply_training_time_weight=True
                )
                gradient = (
                    gradient_matching_loss(prediction, target, inputs)
                    if loss_cfg["lambda_gradient"] > 0
                    else prediction.new_zeros(())
                )
                boundary = (
                    immersed_boundary_loss(positive_raw, inputs)
                    if loss_cfg["lambda_boundary"] > 0
                    else prediction.new_zeros(())
                )
                non_degradation = (
                    conservative_non_degradation_loss(
                        prediction, base_prediction, target, inputs
                    )
                    if loss_cfg["lambda_non_degradation"] > 0
                    else prediction.new_zeros(())
                )
                residual_regularization = (
                    residual_regularization_loss(positive_residual, inputs)
                    if loss_cfg["lambda_residual_regularization"] > 0
                    else prediction.new_zeros(())
                )

                if sigma is not None and loss_cfg["lambda_uncertainty"] > 0:
                    uncertainty = uncertainty_regression_loss(
                        sigma, prediction, target, inputs
                    )
                    calibration = calibration_loss(sigma, prediction, target, inputs)
                else:
                    uncertainty = prediction.new_zeros(())
                    calibration = prediction.new_zeros(())

                weighted_relative = loss_cfg["lambda_relative"] * relative
                weighted_gradient = warmup * loss_cfg["lambda_gradient"] * gradient
                weighted_boundary = warmup * loss_cfg["lambda_boundary"] * boundary
                weighted_non_degradation = (
                    warmup * loss_cfg["lambda_non_degradation"] * non_degradation
                )
                weighted_residual_regularization = (
                    warmup
                    * loss_cfg["lambda_residual_regularization"]
                    * residual_regularization
                )
                weighted_uncertainty = (
                    warmup * loss_cfg["lambda_uncertainty"] * uncertainty
                )
                weighted_calibration = (
                    warmup * loss_cfg["lambda_calibration"] * calibration
                )

                total_loss = (
                    data_mse
                    + weighted_relative
                    + weighted_gradient
                    + weighted_boundary
                    + weighted_non_degradation
                    + weighted_residual_regularization
                    + weighted_uncertainty
                    + weighted_calibration
                )

            if not torch.isfinite(total_loss):
                batch_failed = True
                break

            scaler.scale(total_loss * fraction).backward()
            values = {
                "total": total_loss,
                "data_mse": data_mse,
                "relative": relative,
                "gradient": gradient,
                "boundary": boundary,
                "non_degradation": non_degradation,
                "residual_regularization": residual_regularization,
                "uncertainty": uncertainty,
                "calibration": calibration,
                "weighted_relative": weighted_relative,
                "weighted_gradient": weighted_gradient,
                "weighted_boundary": weighted_boundary,
                "weighted_non_degradation": weighted_non_degradation,
                "weighted_residual_regularization": weighted_residual_regularization,
                "weighted_uncertainty": weighted_uncertainty,
                "weighted_calibration": weighted_calibration,
            }
            for key, value in values.items():
                batch_values[key] += float(value.detach()) * fraction
            del (
                inputs,
                target,
                prediction,
                sigma,
                positive_raw,
                base_prediction,
                positive_residual,
            )

        if batch_failed:
            nonfinite_batches += 1
            optimizer.zero_grad(set_to_none=True)
            if nonfinite_batches > int(CONFIG["max_nonfinite_batches_per_epoch"]):
                raise FloatingPointError(
                    f"Too many non-finite losses in epoch {epoch + 1}"
                )
            continue

        if scaler.is_enabled():
            scaler.unscale_(optimizer)
        trainable_parameters = [
            parameter for parameter in model.parameters() if parameter.requires_grad
        ]
        grad_norm = torch.nn.utils.clip_grad_norm_(
            trainable_parameters,
            max_norm=5.0,
            error_if_nonfinite=False,
            foreach=False,
        )
        if not bool(torch.isfinite(torch.as_tensor(grad_norm)).item()):
            nonfinite_batches += 1
            optimizer.zero_grad(set_to_none=True)
            if scaler.is_enabled():
                scaler.update()
            if nonfinite_batches > int(CONFIG["max_nonfinite_batches_per_epoch"]):
                raise FloatingPointError(
                    f"Too many non-finite gradients in epoch {epoch + 1}"
                )
            continue

        scaler.step(optimizer)
        scaler.update()
        for key in keys:
            totals[key] += batch_values[key]
        successful_batches += 1

    if attempted_batches == 0:
        raise RuntimeError("Training loader produced no batches")
    if successful_batches == 0:
        raise FloatingPointError(f"All {attempted_batches} batches were skipped")

    stats = {key: value / successful_batches for key, value in totals.items()}
    stats.update(
        {
            "successful_batches": float(successful_batches),
            "attempted_batches": float(attempted_batches),
            "nonfinite_batches": float(nonfinite_batches),
            "skipped_fraction": nonfinite_batches / max(attempted_batches, 1),
            "amp_scale": float(scaler.get_scale()) if scaler.is_enabled() else 1.0,
            "physics_ramp": float(warmup),
            "physics_warmup_epochs": float(CONFIG["physics_warmup_epochs"]),
            "spectral_micro_batch_size": float(
                spectral_micro_batch_size(model, int(CONFIG["batch_size"]))
            ),
        }
    )
    return stats


def validate_epoch(model: nn.Module, loader: DataLoader) -> Dict[str, float]:
    model.eval()
    threshold = float(CONFIG.get("initial_time_threshold_scaled", 0.02))

    blend_alpha = float("nan")
    blend_alpha_raw = float("nan")
    blend_validation_mse_gain = float("nan")
    if isinstance(model, PCFNO):
        blend_alpha = float(model.inference_alpha.item())
        blend_alpha_raw = blend_alpha
        blend_validation_mse_gain = 0.0

    mse_total = 0.0
    rel_total = 0.0
    wake_mse_total = 0.0
    boundary_mse_total = 0.0
    initial_mse_total = 0.0
    initial_samples = 0
    high_re_mse_total = 0.0
    high_re_samples = 0
    samples = 0

    with torch.no_grad():
        for batch_index, (inputs_cpu, target_cpu, _, _) in enumerate(loader):
            if CONFIG["max_eval_batches"] is not None and batch_index >= int(
                CONFIG["max_eval_batches"]
            ):
                break
            full_batch_size = int(inputs_cpu.shape[0])
            micro_size = spectral_micro_batch_size(model, full_batch_size)
            for start_index in range(0, full_batch_size, micro_size):
                stop_index = min(start_index + micro_size, full_batch_size)
                current = stop_index - start_index
                inputs = inputs_cpu[start_index:stop_index].to(
                    CONFIG["device"], non_blocking=True
                )
                target = target_cpu[start_index:stop_index].to(
                    CONFIG["device"], non_blocking=True
                )

                with model_autocast_context(model):
                    prediction, _, positive_raw = forward_with_sigma(
                        model, inputs, return_positive_raw=True
                    )
                    mse, relative = supervised_terms(prediction, target, inputs)
                    wake_valid = wake_mask(inputs)
                    wake_mse = (
                        torch.mean((prediction[wake_valid] - target[wake_valid]) ** 2)
                        if torch.any(wake_valid)
                        else mse.new_zeros(())
                    )
                    boundary_values = immersed_boundary_values(positive_raw, inputs)
                    boundary_mse = (
                        torch.mean(boundary_values.square())
                        if boundary_values.numel()
                        else mse.new_zeros(())
                    )
                    initial_batch = inputs[:, INPUT_TIME, 0, 0] <= threshold
                    if torch.any(initial_batch):
                        initial_mask = fluid_mask(inputs[initial_batch])
                        initial_error = (
                            prediction[initial_batch][initial_mask]
                            - target[initial_batch][initial_mask]
                        )
                        initial_mse = torch.mean(initial_error.square())
                        initial_count = int(torch.sum(initial_batch).item())
                    else:
                        initial_mse = mse.new_zeros(())
                        initial_count = 0
                    high_re_batch = inputs[:, INPUT_RE, 0, 0] >= float(
                        CONFIG.get("validation_high_re_threshold_scaled", 0.80)
                    )
                    if torch.any(high_re_batch):
                        high_re_mask = fluid_mask(inputs[high_re_batch])
                        high_re_error = (
                            prediction[high_re_batch][high_re_mask]
                            - target[high_re_batch][high_re_mask]
                        )
                        high_re_mse = torch.mean(high_re_error.square())
                        high_re_count = int(torch.sum(high_re_batch).item())
                    else:
                        high_re_mse = mse.new_zeros(())
                        high_re_count = 0

                if not torch.isfinite(mse):
                    raise FloatingPointError(
                        f"Non-finite validation loss at batch {batch_index + 1}"
                    )
                mse_total += float(mse) * current
                rel_total += float(relative) * current
                wake_mse_total += float(wake_mse) * current
                boundary_mse_total += float(boundary_mse) * current
                initial_mse_total += float(initial_mse) * initial_count
                initial_samples += initial_count
                high_re_mse_total += float(high_re_mse) * high_re_count
                high_re_samples += high_re_count
                samples += current
                del inputs, target, prediction, positive_raw

    if samples == 0:
        raise RuntimeError("Validation loader produced no samples")

    mse = mse_total / samples
    relative = rel_total / samples
    wake_mse = wake_mse_total / samples
    boundary_mse = boundary_mse_total / samples
    initial_mse = initial_mse_total / max(initial_samples, 1)
    high_re_mse = high_re_mse_total / max(high_re_samples, 1)
    objective = (
        mse
        + float(CONFIG.get("validation_wake_weight", 0.05)) * wake_mse
        + float(CONFIG.get("validation_boundary_weight", 0.002)) * boundary_mse
        + float(CONFIG.get("validation_initial_weight", 0.02)) * initial_mse
        + float(CONFIG.get("validation_high_re_weight", 0.0)) * high_re_mse
    )
    return {
        "mse": mse,
        "relative_l2": relative,
        "wake_mse": wake_mse,
        "boundary_mse": boundary_mse,
        "initial_mse": initial_mse,
        "high_re_mse": high_re_mse,
        "objective": objective,
        "blend_alpha": blend_alpha,
        "blend_alpha_raw": blend_alpha_raw,
        "blend_validation_mse_gain": blend_validation_mse_gain,
    }


def uncertainty_head_epoch(
    model: PCFNO,
    loader: DataLoader,
    optimizer: Optional[torch.optim.Optimizer],
) -> float:
    model.eval()
    model.uncertainty_head.train(optimizer is not None)
    total_loss = 0.0
    total_samples = 0
    log_weight = float(CONFIG.get("uncertainty_refinement_log_weight", 0.01))
    for batch_index, (inputs_cpu, target_cpu, _, _) in enumerate(loader):
        if CONFIG["max_eval_batches"] is not None and batch_index >= int(
            CONFIG["max_eval_batches"]
        ):
            break
        full_batch_size = int(inputs_cpu.shape[0])
        micro_size = spectral_micro_batch_size(model, full_batch_size)
        if optimizer is not None:
            optimizer.zero_grad(set_to_none=True)
        for start_index in range(0, full_batch_size, micro_size):
            stop_index = min(start_index + micro_size, full_batch_size)
            current = stop_index - start_index
            inputs = inputs_cpu[start_index:stop_index].to(
                CONFIG["device"], non_blocking=True
            )
            target = target_cpu[start_index:stop_index].to(
                CONFIG["device"], non_blocking=True
            )
            prediction, sigma = forward_with_sigma(model, inputs)
            if sigma is None:
                raise RuntimeError("PCFNO uncertainty head returned no sigma")
            valid = fluid_mask(inputs)
            absolute_error = torch.abs(prediction - target).detach().clamp_min(1e-5)
            sigma_valid = sigma.clamp_min(1e-5)
            regression = torch.mean(
                (sigma_valid[valid] - absolute_error[valid]).square()
            )
            logarithmic = torch.mean(
                (
                    torch.log(sigma_valid[valid]) - torch.log(absolute_error[valid])
                ).square()
            )
            loss = regression + log_weight * logarithmic
            if optimizer is not None:
                fraction = float(current) / float(full_batch_size)
                (loss * fraction).backward()
            total_loss += float(loss.detach()) * current
            total_samples += current
            del inputs, target, prediction, sigma, absolute_error, sigma_valid
        if optimizer is not None:
            torch.nn.utils.clip_grad_norm_(model.uncertainty_head.parameters(), 1.0)
            optimizer.step()
    if total_samples == 0:
        raise RuntimeError("Uncertainty refinement loader produced no samples")
    return total_loss / total_samples


def refine_uncertainty_head(
    model: PCFNO,
    train_loader: DataLoader,
    val_loader: DataLoader,
    name: str,
    seed: int,
    logger: logging.Logger,
) -> Dict[str, Any]:
    requested = max(int(CONFIG.get("uncertainty_refinement_epochs", 0)), 0)
    if requested == 0 or not model.uncertainty:
        return {
            "uncertainty_refinement_epochs": 0,
            "uncertainty_refinement_best_epoch": 0,
            "uncertainty_refinement_best_val_loss": float("nan"),
            "uncertainty_refinement_time_seconds": 0.0,
        }
    for parameter in model.parameters():
        parameter.requires_grad_(False)
    for parameter in model.uncertainty_head.parameters():
        parameter.requires_grad_(True)
    optimizer = torch.optim.AdamW(
        model.uncertainty_head.parameters(),
        lr=float(CONFIG.get("uncertainty_refinement_lr", 1e-3)),
        weight_decay=float(CONFIG.get("weight_decay", 0.0)),
    )
    best_loss = float("inf")
    best_epoch = 0
    best_state: Optional[Dict[str, Any]] = None
    patience = 0
    patience_limit = max(int(CONFIG.get("uncertainty_refinement_patience", 4)), 1)
    min_delta = float(CONFIG.get("uncertainty_refinement_min_delta", 1e-6))
    start = time.perf_counter()
    completed = 0
    for epoch in range(requested):
        train_loss = uncertainty_head_epoch(model, train_loader, optimizer)
        with torch.no_grad():
            val_loss = uncertainty_head_epoch(model, val_loader, None)
        completed = epoch + 1
        logger.info(
            "%s seed %d | uncertainty refinement %d/%d | train %.4e | val %.4e",
            name,
            seed,
            completed,
            requested,
            train_loss,
            val_loss,
        )
        if val_loss < best_loss - min_delta:
            best_loss = val_loss
            best_epoch = completed
            best_state = tree_to_cpu(model.uncertainty_head.state_dict())
            patience = 0
        else:
            patience += 1
        if patience >= patience_limit:
            break
    if best_state is not None:
        model.uncertainty_head.load_state_dict(best_state)
    model.eval()
    elapsed = time.perf_counter() - start
    return {
        "uncertainty_refinement_epochs": completed,
        "uncertainty_refinement_best_epoch": best_epoch,
        "uncertainty_refinement_best_val_loss": best_loss,
        "uncertainty_refinement_time_seconds": elapsed,
    }


def training_signature(name: str) -> str:
    keys = [
        "resolution",
        "epochs",
        "lr",
        "weight_decay",
        "width",
        "modes",
        "fno_layers",
        "padding",
        "time_samples_per_case",
        "split_mode",
        "resolved_train_reynolds",
        "resolved_validation_reynolds",
        "resolved_test_reynolds",
        "resolved_excluded_reynolds",
        "split_profile",
        "physics_warmup_fraction",
        "physics_warmup_epochs",
        "lambda_relative",
        "lambda_gradient",
        "lambda_boundary",
        "lambda_non_degradation",
        "lambda_residual_regularization",
        "lambda_uncertainty",
        "lambda_calibration",
        "positive_output_mode",
        "positive_output_epsilon",
        "enforce_positive_output",
        "hard_zero_solid",
        "spectral_micro_batch_size",
        "validation_wake_weight",
        "validation_boundary_weight",
        "validation_initial_weight",
        "validation_high_re_weight",
        "validation_high_re_threshold_scaled",
        "initial_time_threshold_scaled",
        "initial_sample_weight",
        "pcfno_warm_start_from_fno",
        "pcfno_freeze_backbone",
        "pcfno_adapter_lr",
        "adapter_width",
        "pcfno_max_epochs",
        "pcfno_minimum_epochs_before_early_stopping",
        "pcfno_early_stopping_patience",
        "adapter_correction_limit",
        "adapter_focus_floor",
        "adapter_boundary_scale_over_d",
        "fixed_residual_alpha",
        "unblended_residual_alpha",
        "uncertainty_refinement_epochs",
        "uncertainty_refinement_lr",
        "uncertainty_refinement_patience",
        "uncertainty_refinement_min_delta",
        "uncertainty_refinement_log_weight",
        "uncertainty_plot_percentile",
        "uncertainty_plot_gridsize",
        "strict_re_matching",
        "use_amp",
        "amp_dtype",
    ]
    payload = {
        "model": name,
        "implementation_version": "v15_fixed_blend_decoupled_uncertainty_pcfno",
    }
    payload.update({key: CONFIG.get(key) for key in keys})
    text_value = json.dumps(payload, sort_keys=True, default=str).encode("utf-8")
    return hashlib.sha256(text_value).hexdigest()


def model_requested_epochs(name: str) -> int:
    if name.startswith("PCFNO"):
        return int(CONFIG.get("pcfno_max_epochs", CONFIG["epochs"]))
    return int(CONFIG["epochs"])


def model_minimum_epochs(name: str) -> int:
    if name.startswith("PCFNO"):
        return int(
            CONFIG.get(
                "pcfno_minimum_epochs_before_early_stopping",
                CONFIG.get("minimum_epochs_before_early_stopping", 0),
            )
        )
    return int(CONFIG.get("minimum_epochs_before_early_stopping", 0))


def model_early_stopping_patience(name: str) -> int:
    if name.startswith("PCFNO"):
        return int(
            CONFIG.get(
                "pcfno_early_stopping_patience",
                CONFIG.get("early_stopping_patience", 0),
            )
        )
    return int(CONFIG.get("early_stopping_patience", 0))


def checkpoint_path(output_dir: Path, name: str, seed: int) -> Path:
    return output_dir / "checkpoints" / name / f"seed_{seed}.pt"


def history_path(output_dir: Path, name: str, seed: int) -> Path:
    return output_dir / "data" / "histories" / f"training_history_{name}_seed{seed}.dat"


def train_model(
    name: str,
    seed: int,
    train_dataset: Dataset,
    val_dataset: Dataset,
    output_dir: Path,
    logger: logging.Logger,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    ckpt_path = checkpoint_path(output_dir, name, seed)
    last_path = last_checkpoint_path(output_dir, name, seed)
    hist_path = history_path(output_dir, name, seed)
    ckpt_path.parent.mkdir(parents=True, exist_ok=True)
    existing_history = read_dat(hist_path)

    if (
        bool(CONFIG["resume_existing"])
        and not bool(CONFIG["force_retrain"])
        and ckpt_path.exists()
        and existing_history
    ):
        checkpoint = safe_torch_load(ckpt_path, map_location="cpu")
        completed_epochs = int(existing_history[-1].get("epoch", len(existing_history)))
        signature_matches = checkpoint.get("training_signature") == training_signature(
            name
        )
        training_completed = bool(checkpoint.get("training_completed", False))
        requested_epochs = model_requested_epochs(name)
        if signature_matches and (
            training_completed or completed_epochs >= requested_epochs
        ):
            runtime = {
                "model": name,
                "seed": seed,
                "training_time_seconds": float(
                    existing_history[-1].get("training_time_seconds", float("nan"))
                ),
                "epochs_completed": len(existing_history),
                "best_epoch": int(checkpoint.get("best_epoch", len(existing_history))),
                "best_val_mse": float(checkpoint.get("best_val_mse", float("nan"))),
                "best_val_objective": float(
                    checkpoint.get(
                        "best_val_objective",
                        checkpoint.get("best_val_mse", float("nan")),
                    )
                ),
                "precision_mode": checkpoint.get("precision_mode", "unknown"),
                "stop_reason": checkpoint.get("stop_reason", "completed"),
                "fno_pretraining_time_seconds": float(
                    checkpoint.get("warm_start_metadata", {}).get(
                        "fno_pretraining_time_seconds", 0.0
                    )
                ),
                "effective_training_time_seconds": float(
                    existing_history[-1].get("training_time_seconds", float("nan"))
                )
                + float(
                    checkpoint.get("warm_start_metadata", {}).get(
                        "fno_pretraining_time_seconds", 0.0
                    )
                ),
                "uncertainty_refinement_epochs": int(
                    checkpoint.get("uncertainty_refinement_metadata", {}).get(
                        "uncertainty_refinement_epochs", 0
                    )
                ),
                "uncertainty_refinement_best_epoch": int(
                    checkpoint.get("uncertainty_refinement_metadata", {}).get(
                        "uncertainty_refinement_best_epoch", 0
                    )
                ),
                "uncertainty_refinement_best_val_loss": float(
                    checkpoint.get("uncertainty_refinement_metadata", {}).get(
                        "uncertainty_refinement_best_val_loss", float("nan")
                    )
                ),
                "uncertainty_refinement_time_seconds": float(
                    checkpoint.get("uncertainty_refinement_metadata", {}).get(
                        "uncertainty_refinement_time_seconds", 0.0
                    )
                ),
            }
            logger.info("Skipping %s seed %d: completed checkpoint found", name, seed)
            return existing_history, runtime

    set_seed(seed)
    train_loader = make_loader(train_dataset, shuffle=True, seed=seed)
    val_loader = make_loader(val_dataset, shuffle=False, seed=seed)
    model = model_factory(name).to(CONFIG["device"])
    loss_cfg = loss_configuration(name)

    warm_start_metadata: Dict[str, Any] = {
        "warm_started_from_fno": False,
        "fno_checkpoint": "",
        "fno_pretraining_time_seconds": 0.0,
    }
    if isinstance(model, PCFNO) and bool(CONFIG.get("pcfno_warm_start_from_fno", True)):
        fno_path = checkpoint_path(output_dir, "FNO", seed)
        if not fno_path.exists():
            raise RuntimeError(
                f"{name} seed {seed} requires the matching FNO checkpoint: {fno_path}. "
                "Train FNO first or include FNO in --models."
            )
        fno_checkpoint = safe_torch_load(fno_path, map_location=CONFIG["device"])
        model.load_fno_backbone(fno_checkpoint["state_dict"])
        if bool(CONFIG.get("pcfno_freeze_backbone", True)):
            model.freeze_backbone()
        fno_history = read_dat(history_path(output_dir, "FNO", seed))
        fno_pretraining_time = (
            float(fno_history[-1].get("training_time_seconds", 0.0))
            if fno_history
            else 0.0
        )
        warm_start_metadata = {
            "warm_started_from_fno": True,
            "fno_checkpoint": str(fno_path),
            "fno_pretraining_time_seconds": fno_pretraining_time,
        }
        logger.info(
            "%s seed %d warm-started from FNO checkpoint %s; backbone_frozen=%s",
            name,
            seed,
            fno_path,
            bool(CONFIG.get("pcfno_freeze_backbone", True)),
        )

    trainable_parameters = [
        parameter for parameter in model.parameters() if parameter.requires_grad
    ]
    if not trainable_parameters:
        raise RuntimeError(f"No trainable parameters remain for {name}")
    model_lr = (
        float(CONFIG.get("pcfno_adapter_lr", CONFIG["lr"]))
        if isinstance(model, PCFNO)
        else float(CONFIG["lr"])
    )
    optimizer = torch.optim.AdamW(
        trainable_parameters,
        lr=model_lr,
        weight_decay=float(CONFIG["weight_decay"]),
    )
    requested_epochs = model_requested_epochs(name)
    scheduler = torch.optim.lr_scheduler.CosineAnnealingLR(
        optimizer,
        T_max=max(requested_epochs, 1),
        eta_min=model_lr * 0.02,
    )
    scaler = create_grad_scaler(model)
    precision = precision_mode_name(model)
    best_objective = float("inf")
    best_val_mse = float("inf")
    best_epoch = 0
    patience = 0
    history: List[Dict[str, Any]] = []
    start_epoch = 0
    previous_elapsed = 0.0
    stop_reason = "max_epochs"

    if (
        bool(CONFIG["resume_existing"])
        and not bool(CONFIG["force_retrain"])
        and last_path.exists()
    ):
        recovery = safe_torch_load(last_path, map_location=CONFIG["device"])
        if recovery.get("training_signature") == training_signature(name):
            model.load_state_dict(recovery["state_dict"])
            optimizer.load_state_dict(recovery["optimizer_state_dict"])
            scheduler.load_state_dict(recovery["scheduler_state_dict"])
            if scaler.is_enabled() and recovery.get("scaler_state_dict"):
                scaler.load_state_dict(recovery["scaler_state_dict"])
            best_objective = float(
                recovery.get(
                    "best_val_objective", recovery.get("best_val_mse", best_objective)
                )
            )
            best_val_mse = float(recovery.get("best_val_mse", best_val_mse))
            best_epoch = int(recovery.get("best_epoch", 0))
            patience = int(recovery.get("patience", 0))
            start_epoch = int(recovery.get("epoch", 0))
            previous_elapsed = float(recovery.get("training_time_seconds", 0.0))
            history = [
                row
                for row in existing_history
                if int(row.get("epoch", 0)) <= start_epoch
            ]
            restore_rng_state(recovery)
            if (
                getattr(train_loader, "generator", None) is not None
                and recovery.get("loader_generator_state") is not None
            ):
                train_loader.generator.set_state(recovery["loader_generator_state"])
            logger.info(
                "Recovered %s seed %d from epoch %d (best epoch %d, best objective %.4e)",
                name,
                seed,
                start_epoch,
                best_epoch,
                best_objective,
            )
        else:
            logger.warning(
                "Ignoring incompatible last checkpoint for %s seed %d", name, seed
            )

    logger.info(
        "%s seed %d | precision=%s | complex=%s | micro_batch=%d/%d | start_epoch=%d | loss=%s",
        name,
        seed,
        precision,
        model_has_complex_parameters(model),
        spectral_micro_batch_size(model, int(CONFIG["batch_size"])),
        int(CONFIG["batch_size"]),
        start_epoch,
        loss_cfg,
    )

    start_time = time.perf_counter()
    min_delta = float(CONFIG.get("early_stopping_min_delta", 0.0))
    for epoch in range(start_epoch, requested_epochs):
        train_stats = train_epoch(
            model, train_loader, optimizer, scaler, epoch, loss_cfg
        )
        val_stats = validate_epoch(model, val_loader)
        scheduler.step()
        resources = resource_snapshot()
        row: Dict[str, Any] = {
            "model": name,
            "seed": seed,
            "epoch": epoch + 1,
            "learning_rate": optimizer.param_groups[0]["lr"],
            "val_mse": val_stats["mse"],
            "val_relative_l2": val_stats["relative_l2"],
            "val_wake_mse": val_stats["wake_mse"],
            "val_boundary_mse": val_stats["boundary_mse"],
            "val_initial_mse": val_stats["initial_mse"],
            "val_high_re_mse": val_stats["high_re_mse"],
            "val_objective": val_stats["objective"],
            "val_blend_alpha": val_stats.get("blend_alpha", float("nan")),
            "val_blend_alpha_raw": val_stats.get("blend_alpha_raw", float("nan")),
            "val_blend_validation_mse_gain": val_stats.get(
                "blend_validation_mse_gain", float("nan")
            ),
        }
        row.update({f"train_{key}": value for key, value in train_stats.items()})
        row.update(resources)
        history.append(row)
        logger.info(
            "%s seed %d | epoch %d/%d | train %.4e | val %.4e | wake %.4e | boundary %.4e | initial %.4e | highRe %.4e | obj %.4e | rel %.4e | alpha %.3f | skipped %d/%d | GPU alloc/res %.0f/%.0f MB | RSS %.0f MB",
            name,
            seed,
            epoch + 1,
            requested_epochs,
            train_stats["total"],
            val_stats["mse"],
            val_stats["wake_mse"],
            val_stats["boundary_mse"],
            val_stats["initial_mse"],
            val_stats["high_re_mse"],
            val_stats["objective"],
            val_stats["relative_l2"],
            float(val_stats.get("blend_alpha", float("nan"))),
            int(train_stats["nonfinite_batches"]),
            int(train_stats["attempted_batches"]),
            resources["cuda_allocated_mb"],
            resources["cuda_reserved_mb"],
            resources["cpu_rss_mb"],
        )

        improved = val_stats["objective"] < best_objective - min_delta
        if improved:
            best_objective = val_stats["objective"]
            best_val_mse = val_stats["mse"]
            best_epoch = epoch + 1
            patience = 0
            atomic_torch_save(
                {
                    "model_name": name,
                    "seed": seed,
                    "state_dict": model.state_dict(),
                    "best_val_mse": best_val_mse,
                    "best_val_objective": best_objective,
                    "best_epoch": best_epoch,
                    "requested_epochs": requested_epochs,
                    "precision_mode": precision,
                    "training_signature": training_signature(name),
                    "config": CONFIG,
                    "warm_start_metadata": warm_start_metadata,
                },
                ckpt_path,
            )
        else:
            patience += 1

        elapsed_now = previous_elapsed + (time.perf_counter() - start_time)
        for saved_row in history:
            saved_row["training_time_seconds"] = elapsed_now
            saved_row["best_epoch"] = best_epoch
        if bool(CONFIG.get("save_history_every_epoch", True)):
            save_dat(hist_path, history)

        interval = max(int(CONFIG.get("checkpoint_interval_epochs", 5)), 1)
        can_stop = (epoch + 1) >= model_minimum_epochs(name)
        if (epoch + 1) % interval == 0 or (
            can_stop and patience >= model_early_stopping_patience(name)
        ):
            recovery_payload: Dict[str, Any] = {
                "model_name": name,
                "seed": seed,
                "epoch": epoch + 1,
                "state_dict": model.state_dict(),
                "optimizer_state_dict": optimizer.state_dict(),
                "scheduler_state_dict": scheduler.state_dict(),
                "scaler_state_dict": scaler.state_dict() if scaler.is_enabled() else {},
                "best_val_mse": best_val_mse,
                "best_val_objective": best_objective,
                "best_epoch": best_epoch,
                "patience": patience,
                "requested_epochs": requested_epochs,
                "precision_mode": precision,
                "training_signature": training_signature(name),
                "training_time_seconds": elapsed_now,
                "loader_generator_state": (
                    train_loader.generator.get_state()
                    if getattr(train_loader, "generator", None) is not None
                    else None
                ),
                "config": CONFIG,
                "warm_start_metadata": warm_start_metadata,
            }
            recovery_payload.update(rng_state_payload())
            atomic_torch_save(recovery_payload, last_path)

        cuda_maintenance(epoch, logger)
        if can_stop and patience >= model_early_stopping_patience(name):
            stop_reason = "early_stopping"
            logger.info("Early stopping %s seed %d at epoch %d", name, seed, epoch + 1)
            break

    elapsed = previous_elapsed + (time.perf_counter() - start_time)
    if not ckpt_path.exists():
        raise RuntimeError(f"No checkpoint created for {name}, seed {seed}")
    best_payload = safe_torch_load(ckpt_path, map_location=CONFIG["device"])
    uncertainty_refinement_metadata: Dict[str, Any] = {
        "uncertainty_refinement_epochs": 0,
        "uncertainty_refinement_best_epoch": 0,
        "uncertainty_refinement_best_val_loss": float("nan"),
        "uncertainty_refinement_time_seconds": 0.0,
    }
    if isinstance(model, PCFNO):
        model.load_state_dict(best_payload["state_dict"])
        uncertainty_refinement_metadata = refine_uncertainty_head(
            model, train_loader, val_loader, name, seed, logger
        )
        elapsed += float(
            uncertainty_refinement_metadata["uncertainty_refinement_time_seconds"]
        )
        best_payload["state_dict"] = tree_to_cpu(model.state_dict())
    for row in history:
        row["training_time_seconds"] = elapsed
        row["best_epoch"] = best_epoch
        row["training_completed"] = True
        row["stop_reason"] = stop_reason
        row.update(uncertainty_refinement_metadata)
    save_dat(hist_path, history)
    final_epoch = int(history[-1]["epoch"]) if history else int(start_epoch)
    final_recovery_payload: Dict[str, Any] = {
        "model_name": name,
        "seed": seed,
        "epoch": final_epoch,
        "state_dict": model.state_dict(),
        "optimizer_state_dict": optimizer.state_dict(),
        "scheduler_state_dict": scheduler.state_dict(),
        "scaler_state_dict": scaler.state_dict() if scaler.is_enabled() else {},
        "best_val_mse": best_val_mse,
        "best_val_objective": best_objective,
        "best_epoch": best_epoch,
        "patience": patience,
        "requested_epochs": requested_epochs,
        "precision_mode": precision,
        "training_signature": training_signature(name),
        "training_time_seconds": elapsed,
        "loader_generator_state": (
            train_loader.generator.get_state()
            if getattr(train_loader, "generator", None) is not None
            else None
        ),
        "config": CONFIG,
        "warm_start_metadata": warm_start_metadata,
        "uncertainty_refinement_metadata": uncertainty_refinement_metadata,
        "training_completed": True,
        "stop_reason": stop_reason,
    }
    final_recovery_payload.update(rng_state_payload())
    atomic_torch_save(final_recovery_payload, last_path)
    best_payload = tree_to_cpu(best_payload)
    best_payload["training_completed"] = True
    best_payload["stop_reason"] = stop_reason
    best_payload["epochs_completed"] = final_epoch
    best_payload["training_time_seconds"] = elapsed
    best_payload["uncertainty_refinement_metadata"] = uncertainty_refinement_metadata
    atomic_torch_save(best_payload, ckpt_path)
    runtime = {
        "model": name,
        "seed": seed,
        "training_time_seconds": elapsed,
        "epochs_completed": len(history),
        "best_epoch": best_epoch,
        "best_val_mse": best_val_mse,
        "best_val_objective": best_objective,
        "precision_mode": precision,
        "stop_reason": stop_reason,
        "fno_pretraining_time_seconds": float(
            warm_start_metadata.get("fno_pretraining_time_seconds", 0.0)
        ),
        "effective_training_time_seconds": elapsed
        + float(warm_start_metadata.get("fno_pretraining_time_seconds", 0.0)),
        "uncertainty_refinement_epochs": int(
            uncertainty_refinement_metadata["uncertainty_refinement_epochs"]
        ),
        "uncertainty_refinement_best_epoch": int(
            uncertainty_refinement_metadata["uncertainty_refinement_best_epoch"]
        ),
        "uncertainty_refinement_best_val_loss": float(
            uncertainty_refinement_metadata["uncertainty_refinement_best_val_loss"]
        ),
        "uncertainty_refinement_time_seconds": float(
            uncertainty_refinement_metadata["uncertainty_refinement_time_seconds"]
        ),
    }
    del model, optimizer, scheduler, scaler, train_loader, val_loader
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
        try:
            torch.backends.cuda.cufft_plan_cache.clear()
        except Exception:
            pass
    return history, runtime


def load_trained_model(name: str, seed: int, output_dir: Path) -> nn.Module:
    path = checkpoint_path(output_dir, name, seed)
    if not path.exists():
        raise FileNotFoundError(path)
    set_seed(seed)
    model = model_factory(name).to(CONFIG["device"])
    checkpoint = safe_torch_load(path, map_location=CONFIG["device"])
    model.load_state_dict(checkpoint["state_dict"])
    if isinstance(model, PCFNO):
        model.correction_enabled = True
    model.eval()
    return model


def metadata_value(metadata: Dict[str, Any], key: str, index: int) -> Any:
    value = metadata[key]
    if isinstance(value, torch.Tensor):
        item = value[index].detach().cpu().item()
        return int(item) if not torch.is_floating_point(value) else float(item)
    if isinstance(value, np.ndarray):
        return value[index].item()
    if isinstance(value, (list, tuple)):
        return value[index]
    return value


def pearson_from_arrays(x: np.ndarray, y: np.ndarray) -> float:
    if x.size < 2 or y.size < 2 or np.std(x) < 1e-12 or np.std(y) < 1e-12:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def benchmark_inference(model: nn.Module, loader: DataLoader) -> Dict[str, float]:
    try:
        inputs, _, _, _ = next(iter(loader))
    except StopIteration:
        return {
            "Inference_ms_per_frame": float("nan"),
            "Inference_ms_per_frame_std": float("nan"),
            "Throughput_frames_per_second": float("nan"),
            "Runtime_batch_size": 0,
        }
    inputs = inputs.to(CONFIG["device"], non_blocking=True)
    warmup = int(CONFIG["runtime_warmup"])
    repeats = int(CONFIG["runtime_repeats"])
    with torch.no_grad():
        for _ in range(warmup):
            with model_autocast_context(model):
                _ = model(inputs)
        if torch.cuda.is_available() and str(CONFIG["device"]).startswith("cuda"):
            torch.cuda.synchronize()
        timings: List[float] = []
        for _ in range(repeats):
            if torch.cuda.is_available() and str(CONFIG["device"]).startswith("cuda"):
                torch.cuda.synchronize()
            start = time.perf_counter()
            with model_autocast_context(model):
                _ = model(inputs)
            if torch.cuda.is_available() and str(CONFIG["device"]).startswith("cuda"):
                torch.cuda.synchronize()
            timings.append(time.perf_counter() - start)
    per_frame_ms = (
        np.asarray(timings, dtype=np.float64) * 1000.0 / max(inputs.shape[0], 1)
    )
    mean_ms = float(np.mean(per_frame_ms))
    return {
        "Inference_ms_per_frame": mean_ms,
        "Inference_ms_per_frame_std": (
            float(np.std(per_frame_ms, ddof=1)) if len(per_frame_ms) > 1 else 0.0
        ),
        "Throughput_frames_per_second": 1000.0 / max(mean_ms, 1e-12),
        "Runtime_batch_size": int(inputs.shape[0]),
    }


def sigma_scale_path(output_dir: Path, name: str, seed: int) -> Path:
    return output_dir / "data" / "uncertainty" / f"sigma_scale_{name}_seed{seed}.dat"


def load_sigma_scale(output_dir: Path, name: str, seed: int) -> float:
    path = sigma_scale_path(output_dir, name, seed)
    rows = read_dat(path)
    if rows:
        return float(rows[0].get("sigma_scale", 1.0))
    return 1.0


def fit_uncertainty_scale(
    name: str,
    seed: int,
    val_dataset: Dataset,
    output_dir: Path,
    logger: logging.Logger,
) -> Dict[str, Any]:
    model = load_trained_model(name, seed, output_dir)
    if not isinstance(model, PCFNO) or not model.uncertainty:
        del model
        return {
            "model": name,
            "seed": seed,
            "sigma_scale": 1.0,
            "calibration_enabled": False,
        }
    loader = make_loader(val_dataset, shuffle=False, seed=seed)
    ratio_sq_sum = 0.0
    ratio_samples: List[np.ndarray] = []
    count = 0
    raw_cover1 = 0
    raw_cover2 = 0
    raw_nll_sum = 0.0
    model.eval()
    with torch.no_grad():
        for batch_index, (inputs, target, _, _) in enumerate(loader):
            if CONFIG["max_eval_batches"] is not None and batch_index >= int(
                CONFIG["max_eval_batches"]
            ):
                break
            inputs = inputs.to(CONFIG["device"], non_blocking=True)
            target = target.to(CONFIG["device"], non_blocking=True)
            with model_autocast_context(model):
                prediction, sigma = forward_with_sigma(model, inputs)
            if sigma is None:
                continue
            valid = fluid_mask(inputs)
            error = torch.abs(prediction - target)[valid]
            sigma_valid = sigma[valid].clamp_min(1e-6)
            ratio = error / sigma_valid
            ratio_sq_sum += torch.sum(ratio**2).item()
            remaining = int(CONFIG["uncertainty_subsample"]) - sum(
                values.size for values in ratio_samples
            )
            if remaining > 0:
                ratio_np = ratio.detach().float().cpu().numpy().ravel()
                if ratio_np.size > remaining:
                    rng = np.random.default_rng(seed + batch_index)
                    chosen = rng.choice(ratio_np.size, remaining, replace=False)
                    ratio_np = ratio_np[chosen]
                ratio_samples.append(ratio_np)
            count += error.numel()
            raw_cover1 += int(torch.sum(error <= sigma_valid).item())
            raw_cover2 += int(torch.sum(error <= 2.0 * sigma_valid).item())
            raw_nll_sum += torch.sum(
                0.5 * (error / sigma_valid) ** 2 + torch.log(sigma_valid)
            ).item()
    nll_scale = math.sqrt(ratio_sq_sum / max(count, 1))
    ratios = np.concatenate(ratio_samples) if ratio_samples else np.array([])
    if ratios.size:
        multipliers = np.linspace(0.35, 2.75, 25)
        nominal = np.asarray(
            [math.erf(value / math.sqrt(2.0)) for value in multipliers]
        )
        coverage_scales = np.asarray(
            [
                np.quantile(ratios, probability) / multiplier
                for probability, multiplier in zip(nominal, multipliers)
                if 0.02 < probability < 0.995
            ],
            dtype=np.float64,
        )
        scale = float(np.median(coverage_scales))
        calibration_method = "median_multi_coverage_quantile_scale"
    else:
        scale = float(nll_scale)
        calibration_method = "gaussian_nll_scale_fallback"
    scale = float(
        np.clip(
            scale,
            float(CONFIG["sigma_scale_min"]),
            float(CONFIG["sigma_scale_max"]),
        )
    )
    calibrated_cover1 = 0
    calibrated_cover2 = 0
    calibrated_nll_sum = 0.0
    if count > 0:
        loader = make_loader(val_dataset, shuffle=False, seed=seed)
        with torch.no_grad():
            for batch_index, (inputs, target, _, _) in enumerate(loader):
                if CONFIG["max_eval_batches"] is not None and batch_index >= int(
                    CONFIG["max_eval_batches"]
                ):
                    break
                inputs = inputs.to(CONFIG["device"], non_blocking=True)
                target = target.to(CONFIG["device"], non_blocking=True)
                with model_autocast_context(model):
                    prediction, sigma = forward_with_sigma(model, inputs)
                if sigma is None:
                    continue
                valid = fluid_mask(inputs)
                error = torch.abs(prediction - target)[valid]
                sigma_valid = (sigma[valid] * scale).clamp_min(1e-6)
                calibrated_cover1 += int(torch.sum(error <= sigma_valid).item())
                calibrated_cover2 += int(torch.sum(error <= 2.0 * sigma_valid).item())
                calibrated_nll_sum += torch.sum(
                    0.5 * (error / sigma_valid) ** 2 + torch.log(sigma_valid)
                ).item()
    row = {
        "model": name,
        "seed": seed,
        "sigma_scale": scale,
        "nll_optimal_scale": nll_scale,
        "calibration_method": calibration_method,
        "calibration_enabled": bool(CONFIG.get("posthoc_sigma_calibration", True)),
        "validation_points": count,
        "raw_coverage_1sigma": raw_cover1 / max(count, 1),
        "raw_coverage_2sigma": raw_cover2 / max(count, 1),
        "calibrated_coverage_1sigma": calibrated_cover1 / max(count, 1),
        "calibrated_coverage_2sigma": calibrated_cover2 / max(count, 1),
        "raw_gaussian_nll": raw_nll_sum / max(count, 1),
        "calibrated_gaussian_nll": calibrated_nll_sum / max(count, 1),
    }
    save_dat(sigma_scale_path(output_dir, name, seed), [row])
    logger.info(
        "%s seed %d validation sigma scale=%.4f | raw coverage=(%.3f, %.3f) | calibrated=(%.3f, %.3f)",
        name,
        seed,
        scale,
        row["raw_coverage_1sigma"],
        row["raw_coverage_2sigma"],
        row["calibrated_coverage_1sigma"],
        row["calibrated_coverage_2sigma"],
    )
    del model
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return row


def region_masks(inputs: torch.Tensor) -> Dict[str, torch.Tensor]:
    fluid = fluid_mask(inputs)
    x = inputs[:, INPUT_XD : INPUT_XD + 1]
    y = inputs[:, INPUT_YD : INPUT_YD + 1]
    wake = (
        fluid
        & (x >= float(CONFIG["wake_x_min_over_d"]))
        & (x <= float(CONFIG["wake_x_max_over_d"]))
        & (torch.abs(y) <= float(CONFIG["wake_half_height_over_d"]))
    )
    distance_over_d = torch.clamp(4.0 * inputs[:, INPUT_SDF : INPUT_SDF + 1], min=0.0)
    near_wall = (
        fluid
        & (distance_over_d > 0.0)
        & (distance_over_d <= float(CONFIG["near_wall_distance_over_d"]))
    )
    return {"Global": fluid, "Wake": wake, "Near_wall": near_wall}


def tensor_region_metrics(
    prediction: torch.Tensor,
    target: torch.Tensor,
    mask: torch.Tensor,
) -> Dict[str, float]:
    if not torch.any(mask):
        return {
            "count": 0.0,
            "sq": 0.0,
            "abs": 0.0,
            "target_sq": 0.0,
            "target_sum": 0.0,
            "Relative_L2": float("nan"),
            "RMSE": float("nan"),
            "MAE": float("nan"),
            "R2": float("nan"),
        }
    error = prediction[mask] - target[mask]
    true = target[mask]
    count = float(error.numel())
    sq = float(torch.sum(error**2).item())
    absolute = float(torch.sum(torch.abs(error)).item())
    target_sq = float(torch.sum(true**2).item())
    target_sum = float(torch.sum(true).item())
    variance_sum = max(target_sq - target_sum**2 / max(count, 1.0), 1e-12)
    return {
        "count": count,
        "sq": sq,
        "abs": absolute,
        "target_sq": target_sq,
        "target_sum": target_sum,
        "Relative_L2": math.sqrt(sq / max(target_sq, 1e-12)),
        "RMSE": math.sqrt(sq / max(count, 1.0)),
        "MAE": absolute / max(count, 1.0),
        "R2": 1.0 - sq / variance_sum,
    }


def evaluate_model(
    name: str,
    seed: int,
    test_dataset: Dataset,
    output_dir: Path,
    logger: logging.Logger,
) -> Tuple[Dict[str, Any], List[Dict[str, Any]]]:
    loader = make_loader(test_dataset, shuffle=False, seed=seed)
    model = load_trained_model(name, seed, output_dir)
    parameter_count = sum(parameter.numel() for parameter in model.parameters())
    runtime = benchmark_inference(model, loader)
    sigma_scale = (
        load_sigma_scale(output_dir, name, seed)
        if bool(CONFIG.get("posthoc_sigma_calibration", True))
        else 1.0
    )
    per_sample: List[Dict[str, Any]] = []
    region_sums: Dict[str, Dict[str, float]] = {
        region: {
            "sq": 0.0,
            "abs": 0.0,
            "target_sq": 0.0,
            "target_sum": 0.0,
            "count": 0.0,
        }
        for region in ["Global", "Wake", "Near_wall"]
    }
    sums = {
        "sq_physical": 0.0,
        "abs_physical": 0.0,
        "gradient_sq": 0.0,
        "gradient_count": 0,
        "boundary_sq": 0.0,
        "boundary_count": 0,
        "negative_count": 0,
        "fluid_count": 0,
        "sigma_raw_sum": 0.0,
        "sigma_cal_sum": 0.0,
        "unc_count": 0,
        "cover1_raw": 0,
        "cover2_raw": 0,
        "cover1_cal": 0,
        "cover2_cal": 0,
        "nll_raw_sum": 0.0,
        "nll_cal_sum": 0.0,
    }
    sampled_errors: List[np.ndarray] = []
    sampled_sigmas_raw: List[np.ndarray] = []
    sampled_sigmas_cal: List[np.ndarray] = []
    model.eval()
    with torch.no_grad():
        for batch_index, (inputs, target, wall_mask, metadata) in enumerate(loader):
            if CONFIG["max_eval_batches"] is not None and batch_index >= int(
                CONFIG["max_eval_batches"]
            ):
                break
            inputs = inputs.to(CONFIG["device"], non_blocking=True)
            target = target.to(CONFIG["device"], non_blocking=True)
            wall_mask = wall_mask.to(CONFIG["device"], non_blocking=True)
            with model_autocast_context(model):
                prediction, sigma_raw, positive_raw = forward_with_sigma(
                    model, inputs, return_positive_raw=True
                )
            sigma_cal = sigma_raw * sigma_scale if sigma_raw is not None else None
            masks = region_masks(inputs)
            error_dim = prediction - target
            pred_dx, pred_dy = spatial_gradients(prediction[:, 0], inputs)
            true_dx, true_dy = spatial_gradients(target[:, 0], inputs)
            grad_valid = valid_gradient_mask(inputs)[:, 0]
            for batch_item in range(inputs.shape[0]):
                row: Dict[str, Any] = {
                    "model": name,
                    "seed": seed,
                    "case": metadata_value(metadata, "case", batch_item),
                    "case_number": metadata_value(metadata, "case_number", batch_item),
                    "frame_index": metadata_value(metadata, "frame_index", batch_item),
                    "time": metadata_value(metadata, "time", batch_item),
                    "time_star": metadata_value(metadata, "time_star", batch_item),
                    "Re": metadata_value(metadata, "Re", batch_item),
                    "Uin": metadata_value(metadata, "Uin", batch_item),
                }
                for region_name, region_mask in masks.items():
                    metrics = tensor_region_metrics(
                        prediction[batch_item : batch_item + 1],
                        target[batch_item : batch_item + 1],
                        region_mask[batch_item : batch_item + 1],
                    )
                    prefix = "" if region_name == "Global" else f"{region_name}_"
                    row[f"{prefix}Relative_L2"] = metrics["Relative_L2"]
                    row[f"{prefix}RMSE_over_Uin"] = metrics["RMSE"]
                    row[f"{prefix}MAE_over_Uin"] = metrics["MAE"]
                    row[f"{prefix}R2"] = metrics["R2"]
                uin = float(row["Uin"])
                row["RMSE_m_per_s"] = float(row["RMSE_over_Uin"]) * uin
                row["MAE_m_per_s"] = float(row["MAE_over_Uin"]) * uin
                gv = grad_valid[batch_item]
                row["Gradient_RMSE"] = torch.sqrt(
                    torch.mean(
                        (pred_dx[batch_item][gv] - true_dx[batch_item][gv]) ** 2
                        + (pred_dy[batch_item][gv] - true_dy[batch_item][gv]) ** 2
                    )
                ).item()
                boundary_values = immersed_boundary_values(
                    positive_raw[batch_item : batch_item + 1],
                    inputs[batch_item : batch_item + 1],
                )
                boundary_rmse = (
                    torch.sqrt(torch.mean(boundary_values**2)).item()
                    if boundary_values.numel()
                    else float("nan")
                )
                row["Boundary_RMSE_over_Uin"] = boundary_rmse
                row["Wall_RMSE_over_Uin"] = boundary_rmse
                global_valid = masks["Global"][batch_item, 0]
                pred_valid = prediction[batch_item, 0][global_valid]
                row["Negative_fraction"] = torch.mean((pred_valid < 0).float()).item()
                if sigma_raw is not None and sigma_cal is not None:
                    raw_item = sigma_raw[batch_item, 0][global_valid].clamp_min(1e-6)
                    cal_item = sigma_cal[batch_item, 0][global_valid].clamp_min(1e-6)
                    abs_error = torch.abs(error_dim[batch_item, 0][global_valid])
                    row.update(
                        {
                            "Sigma_scale": sigma_scale,
                            "Mean_sigma_raw_over_Uin": torch.mean(raw_item).item(),
                            "Mean_sigma_over_Uin": torch.mean(cal_item).item(),
                            "Coverage_1sigma_raw": torch.mean(
                                (abs_error <= raw_item).float()
                            ).item(),
                            "Coverage_2sigma_raw": torch.mean(
                                (abs_error <= 2.0 * raw_item).float()
                            ).item(),
                            "Coverage_1sigma": torch.mean(
                                (abs_error <= cal_item).float()
                            ).item(),
                            "Coverage_2sigma": torch.mean(
                                (abs_error <= 2.0 * cal_item).float()
                            ).item(),
                            "Gaussian_NLL_raw": torch.mean(
                                0.5
                                * (error_dim[batch_item, 0][global_valid] / raw_item)
                                ** 2
                                + torch.log(raw_item)
                            ).item(),
                            "Gaussian_NLL": torch.mean(
                                0.5
                                * (error_dim[batch_item, 0][global_valid] / cal_item)
                                ** 2
                                + torch.log(cal_item)
                            ).item(),
                        }
                    )
                per_sample.append(row)
            for region_name, region_mask in masks.items():
                metrics = tensor_region_metrics(prediction, target, region_mask)
                for key in ["sq", "abs", "target_sq", "target_sum", "count"]:
                    region_sums[region_name][key] += float(metrics[key])
            global_mask = masks["Global"]
            uin_field = inputs[:, INPUT_UIN : INPUT_UIN + 1]
            physical_error = ((prediction - target) * uin_field)[global_mask]
            sums["sq_physical"] += torch.sum(physical_error**2).item()
            sums["abs_physical"] += torch.sum(torch.abs(physical_error)).item()
            gradient_error = ((pred_dx - true_dx) ** 2 + (pred_dy - true_dy) ** 2)[
                grad_valid
            ]
            sums["gradient_sq"] += torch.sum(gradient_error).item()
            sums["gradient_count"] += gradient_error.numel()
            boundary_values = immersed_boundary_values(positive_raw, inputs)
            sums["boundary_sq"] += torch.sum(boundary_values**2).item()
            sums["boundary_count"] += int(boundary_values.numel())
            sums["negative_count"] += int(
                torch.sum((prediction[global_mask] < 0)).item()
            )
            sums["fluid_count"] += int(torch.sum(global_mask).item())
            if sigma_raw is not None and sigma_cal is not None:
                raw_valid = sigma_raw[global_mask].clamp_min(1e-6)
                cal_valid = sigma_cal[global_mask].clamp_min(1e-6)
                signed_error = error_dim[global_mask]
                abs_error_valid = torch.abs(signed_error)
                sums["sigma_raw_sum"] += torch.sum(raw_valid).item()
                sums["sigma_cal_sum"] += torch.sum(cal_valid).item()
                sums["unc_count"] += raw_valid.numel()
                sums["cover1_raw"] += int(
                    torch.sum(abs_error_valid <= raw_valid).item()
                )
                sums["cover2_raw"] += int(
                    torch.sum(abs_error_valid <= 2.0 * raw_valid).item()
                )
                sums["cover1_cal"] += int(
                    torch.sum(abs_error_valid <= cal_valid).item()
                )
                sums["cover2_cal"] += int(
                    torch.sum(abs_error_valid <= 2.0 * cal_valid).item()
                )
                sums["nll_raw_sum"] += torch.sum(
                    0.5 * (signed_error / raw_valid) ** 2 + torch.log(raw_valid)
                ).item()
                sums["nll_cal_sum"] += torch.sum(
                    0.5 * (signed_error / cal_valid) ** 2 + torch.log(cal_valid)
                ).item()
                remaining = int(CONFIG["uncertainty_subsample"]) - sum(
                    array.size for array in sampled_errors
                )
                if remaining > 0:
                    error_np = abs_error_valid.detach().float().cpu().numpy().ravel()
                    raw_np = raw_valid.detach().float().cpu().numpy().ravel()
                    cal_np = cal_valid.detach().float().cpu().numpy().ravel()
                    if error_np.size > remaining:
                        rng = np.random.default_rng(seed + batch_index)
                        selected = rng.choice(error_np.size, remaining, replace=False)
                        error_np = error_np[selected]
                        raw_np = raw_np[selected]
                        cal_np = cal_np[selected]
                    sampled_errors.append(error_np)
                    sampled_sigmas_raw.append(raw_np)
                    sampled_sigmas_cal.append(cal_np)
    global_values = region_sums["Global"]
    global_count = max(global_values["count"], 1.0)
    global_variance = max(
        global_values["target_sq"] - global_values["target_sum"] ** 2 / global_count,
        1e-12,
    )
    summary: Dict[str, Any] = {
        "model": name,
        "seed": seed,
        "MSE_dimensionless": global_values["sq"] / global_count,
        "RMSE_over_Uin": math.sqrt(global_values["sq"] / global_count),
        "MAE_over_Uin": global_values["abs"] / global_count,
        "RMSE_m_per_s": math.sqrt(sums["sq_physical"] / global_count),
        "MAE_m_per_s": sums["abs_physical"] / global_count,
        "Relative_L2": math.sqrt(
            global_values["sq"] / max(global_values["target_sq"], 1e-12)
        ),
        "R2": 1.0 - global_values["sq"] / global_variance,
        "Gradient_RMSE": math.sqrt(
            sums["gradient_sq"] / max(int(sums["gradient_count"]), 1)
        ),
        "Boundary_RMSE_over_Uin": math.sqrt(
            sums["boundary_sq"] / max(int(sums["boundary_count"]), 1)
        ),
        "Wall_RMSE_over_Uin": math.sqrt(
            sums["boundary_sq"] / max(int(sums["boundary_count"]), 1)
        ),
        "Negative_fraction": sums["negative_count"] / max(int(sums["fluid_count"]), 1),
        "Mean_sample_Relative_L2": float(
            np.mean([row["Relative_L2"] for row in per_sample])
        ),
        "Std_sample_Relative_L2": float(
            np.std([row["Relative_L2"] for row in per_sample])
        ),
        "Parameters": parameter_count,
        "Evaluated_frames": len(per_sample),
    }
    for region_name in ["Wake", "Near_wall"]:
        values = region_sums[region_name]
        count = max(values["count"], 1.0)
        variance = max(values["target_sq"] - values["target_sum"] ** 2 / count, 1e-12)
        summary[f"{region_name}_Relative_L2"] = math.sqrt(
            values["sq"] / max(values["target_sq"], 1e-12)
        )
        summary[f"{region_name}_RMSE_over_Uin"] = math.sqrt(values["sq"] / count)
        summary[f"{region_name}_MAE_over_Uin"] = values["abs"] / count
        summary[f"{region_name}_R2"] = 1.0 - values["sq"] / variance
        summary[f"{region_name}_points"] = values["count"]
    summary.update(runtime)
    summary["Reference_speedup_vs_CFD"] = (
        float(CONFIG["cfd_reference_seconds_per_frame"])
        * 1000.0
        / max(summary["Inference_ms_per_frame"], 1e-12)
    )
    if sums["unc_count"] > 0:
        errors = np.concatenate(sampled_errors) if sampled_errors else np.array([])
        sigmas_raw = (
            np.concatenate(sampled_sigmas_raw) if sampled_sigmas_raw else np.array([])
        )
        sigmas_cal = (
            np.concatenate(sampled_sigmas_cal) if sampled_sigmas_cal else np.array([])
        )
        count_unc = max(int(sums["unc_count"]), 1)
        summary.update(
            {
                "Sigma_scale": sigma_scale,
                "Mean_sigma_raw_over_Uin": sums["sigma_raw_sum"] / count_unc,
                "Mean_sigma_over_Uin": sums["sigma_cal_sum"] / count_unc,
                "Gaussian_NLL_raw": sums["nll_raw_sum"] / count_unc,
                "Gaussian_NLL": sums["nll_cal_sum"] / count_unc,
                "Coverage_1sigma_raw": sums["cover1_raw"] / count_unc,
                "Coverage_2sigma_raw": sums["cover2_raw"] / count_unc,
                "Coverage_1sigma": sums["cover1_cal"] / count_unc,
                "Coverage_2sigma": sums["cover2_cal"] / count_unc,
                "Uncertainty_error_correlation_raw": pearson_from_arrays(
                    sigmas_raw, errors
                ),
                "Uncertainty_error_correlation": pearson_from_arrays(
                    sigmas_cal, errors
                ),
                "Calibration_error_raw": abs(sums["cover1_raw"] / count_unc - 0.6827)
                + abs(sums["cover2_raw"] / count_unc - 0.9545),
                "Calibration_error": abs(sums["cover1_cal"] / count_unc - 0.6827)
                + abs(sums["cover2_cal"] / count_unc - 0.9545),
            }
        )
        uncertainty_dir = output_dir / "data" / "uncertainty"
        uncertainty_dir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            uncertainty_dir / f"uncertainty_{name}_seed{seed}.npz",
            error=errors,
            sigma_raw=sigmas_raw,
            sigma_calibrated=sigmas_cal,
            sigma_scale=np.asarray([sigma_scale], dtype=np.float64),
        )
    save_dat(
        output_dir / "data" / "per_sample" / f"per_sample_{name}_seed{seed}.dat",
        per_sample,
    )
    logger.info("%s seed %d metrics: %s", name, seed, json.dumps(summary, indent=2))
    del model
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    return summary, per_sample


def aggregate_model_summaries(
    seed_rows: Sequence[Dict[str, Any]]
) -> List[Dict[str, Any]]:
    grouped: Dict[str, List[Dict[str, Any]]] = {}
    for row in seed_rows:
        grouped.setdefault(str(row["model"]), []).append(row)
    result: List[Dict[str, Any]] = []
    for model_name, rows in grouped.items():
        aggregate: Dict[str, Any] = {"model": model_name, "n_seeds": len(rows)}
        numeric_keys: List[str] = []
        for row in rows:
            for key, value in row.items():
                if key in {"model", "seed"}:
                    continue
                if (
                    isinstance(value, (int, float, np.integer, np.floating))
                    and key not in numeric_keys
                ):
                    numeric_keys.append(key)
        for key in numeric_keys:
            values = np.asarray(
                [
                    float(row[key])
                    for row in rows
                    if key in row and np.isfinite(float(row[key]))
                ],
                dtype=np.float64,
            )
            if values.size:
                aggregate[f"{key}_mean"] = float(np.mean(values))
                aggregate[f"{key}_std"] = (
                    float(np.std(values, ddof=1)) if values.size > 1 else 0.0
                )
                aggregate[f"{key}_median"] = float(np.median(values))
        result.append(aggregate)
    order = {name: index for index, name in enumerate(ALL_EXPERIMENTS)}
    return sorted(result, key=lambda row: order.get(str(row["model"]), 999))


def rows_for_names(
    rows: Sequence[Dict[str, Any]], names: Sequence[str]
) -> List[Dict[str, Any]]:
    mapping = {str(row["model"]): row for row in rows}
    return [mapping[name] for name in names if name in mapping]


def best_seed_for_model(name: str, seeds: Sequence[int], output_dir: Path) -> int:
    candidates: List[Tuple[float, int]] = []
    for seed in seeds:
        path = checkpoint_path(output_dir, name, seed)
        if path.exists():
            checkpoint = safe_torch_load(path, map_location="cpu")
            candidates.append(
                (
                    float(
                        checkpoint.get(
                            "best_val_objective",
                            checkpoint.get("best_val_mse", float("inf")),
                        )
                    ),
                    int(seed),
                )
            )
    if not candidates:
        raise RuntimeError(f"No checkpoints found for {name}")
    return min(candidates)[1]


def dataset_indices_for_case(
    dataset: CFDBenchSpeedDataset, case_name: str
) -> List[int]:
    return [
        index
        for index, (case_index, _) in enumerate(dataset.samples)
        if dataset.cases[case_index].name == case_name
    ]


def predict_indices(
    model: nn.Module,
    dataset: CFDBenchSpeedDataset,
    indices: Sequence[int],
    sigma_scale: float = 1.0,
) -> List[Dict[str, Any]]:
    records: List[Dict[str, Any]] = []
    model.eval()
    with torch.no_grad():
        for index in indices:
            inputs, target, wall_mask, metadata = dataset[index]
            inputs_batch = inputs[None].to(CONFIG["device"])
            with model_autocast_context(model):
                prediction, sigma = forward_with_sigma(model, inputs_batch)
            if sigma is not None:
                sigma = sigma * float(sigma_scale)
            records.append(
                {
                    "inputs": inputs.numpy(),
                    "target": target.numpy(),
                    "prediction": prediction[0].detach().float().cpu().numpy(),
                    "sigma": (
                        sigma[0].detach().float().cpu().numpy()
                        if sigma is not None
                        else None
                    ),
                    "sigma_scale": float(sigma_scale),
                    "wall_mask": wall_mask.numpy(),
                    "metadata": metadata,
                }
            )
    return records


def plot_ready_field(field: np.ndarray, inputs: np.ndarray) -> np.ndarray:
    result = np.asarray(field, dtype=np.float64).copy()
    result[inputs[INPUT_SOLID] > 0.5] = np.nan
    return result


def solid_aware_cmap(name: str) -> Any:
    cmap = plt.get_cmap(name).copy()
    cmap.set_bad("#D9D9D9")
    return cmap


def contour_solid(ax: plt.Axes, inputs: np.ndarray) -> None:
    cylinder = Circle(
        (0.0, 0.0),
        radius=0.5,
        facecolor="#D9D9D9",
        edgecolor="black",
        linewidth=0.75,
        zorder=5,
        antialiased=True,
    )
    ax.add_patch(cylinder)
    ax.set_aspect("equal", adjustable="box")


def field_extent(inputs: np.ndarray) -> Tuple[float, float, float, float]:
    x = inputs[INPUT_XD, 0, :]
    y = inputs[INPUT_YD, :, 0]
    return float(x.min()), float(x.max()), float(y.min()), float(y.max())


def plot_training_convergence(
    histories: Dict[Tuple[str, int], List[Dict[str, Any]]], output_dir: Path
) -> None:
    figure_rows: List[Dict[str, Any]] = []
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35), constrained_layout=True)
    panels = [
        (axes[0], BASELINE_NAMES, "Baseline models"),
        (axes[1], ["PCFNO"] + ABLATION_NAMES, "PCFNO ablations"),
    ]
    aggregated: Dict[str, Tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for panel_index, (ax, names, title) in enumerate(panels):
        for name in names:
            seed_histories = [
                rows for (model, _), rows in histories.items() if model == name and rows
            ]
            if not seed_histories:
                continue
            seed_maxima = [
                max(int(row["epoch"]) for row in rows) for rows in seed_histories
            ]
            common_max_epoch = min(seed_maxima)
            epoch_values: List[int] = []
            means: List[float] = []
            stds: List[float] = []
            counts: List[int] = []
            for epoch in range(1, common_max_epoch + 1):
                values = [
                    float(row["val_mse"])
                    for rows in seed_histories
                    for row in rows
                    if int(row["epoch"]) == epoch
                ]
                if len(values) != len(seed_histories):
                    continue
                epoch_values.append(epoch)
                means.append(float(np.mean(values)))
                stds.append(float(np.std(values, ddof=1)) if len(values) > 1 else 0.0)
                counts.append(len(values))
                figure_rows.append(
                    {
                        "record_type": "aggregate_common_epoch",
                        "panel": panel_index,
                        "model": name,
                        "epoch": epoch,
                        "val_mse_mean": means[-1],
                        "val_mse_std": stds[-1],
                        "contributing_seeds": len(values),
                    }
                )
            if not epoch_values:
                continue
            epochs_array = np.asarray(epoch_values, dtype=np.float64)
            means_array = np.asarray(means, dtype=np.float64)
            stds_array = np.asarray(stds, dtype=np.float64)
            counts_array = np.maximum(np.asarray(counts, dtype=np.float64), 1.0)
            aggregated[f"{panel_index}:{name}"] = (
                epochs_array,
                means_array,
                stds_array,
            )
            color = MODEL_COLORS[name]
            ax.plot(epochs_array, means_array, label=DISPLAY_NAMES[name], color=color)
            half_width = 1.96 * stds_array / np.sqrt(counts_array)
            lower = np.maximum(means_array - half_width, means_array * 0.20)
            upper = means_array + half_width
            ax.fill_between(
                epochs_array,
                lower,
                upper,
                color=color,
                alpha=0.13,
                linewidth=0,
            )
        ax.set_yscale("log")
        ax.set_xlabel("Epoch")
        ax.set_ylabel("Validation MSE")
        ax.set_title(title)
        ax.grid(True, which="both", alpha=0.2, linewidth=0.6)
        ax.legend(ncol=2 if len(names) > 4 else 1, loc="best")
    pcfno_curve = aggregated.get("0:PCFNO")
    if pcfno_curve is not None and pcfno_curve[0].size > 1:
        inset = axes[0].inset_axes([0.53, 0.08, 0.43, 0.29])
        epochs_array, means_array, _ = pcfno_curve
        baseline_value = max(float(means_array[0]), 1e-15)
        relative_change = 100.0 * (means_array / baseline_value - 1.0)
        inset.plot(
            epochs_array, relative_change, color=MODEL_COLORS["PCFNO"], linewidth=1.25
        )
        best_index = int(np.argmin(means_array))
        inset.plot(
            epochs_array[best_index],
            relative_change[best_index],
            marker="o",
            markersize=3.2,
            color=MODEL_COLORS["PCFNO"],
        )
        inset.axhline(0.0, color="#000000", linewidth=0.65, linestyle="--")
        inset.set_xlabel("Epoch", fontsize=6.8)
        inset.set_ylabel("PCFNO change (%)", fontsize=6.8)
        inset.tick_params(labelsize=6.4)
        inset.grid(True, alpha=0.18, linewidth=0.45)
    for (name, seed), rows in histories.items():
        for row in rows:
            raw = dict(row)
            raw["record_type"] = "raw_seed"
            figure_rows.append(raw)
    stem = "fig01_training_convergence"
    save_dat(output_dir / "figure_data" / f"{stem}.dat", figure_rows)
    save_figure(fig, output_dir / "figures" / stem)


def plot_pcfno_loss_decomposition(
    histories: Dict[Tuple[str, int], List[Dict[str, Any]]], output_dir: Path
) -> None:
    rows = [rows for (name, _), rows in histories.items() if name == "PCFNO" and rows]
    if not rows:
        return
    components = [
        ("train_data_mse", "Data MSE", "#0072B2"),
        ("train_weighted_relative", "Relative term", "#E69F00"),
        ("train_weighted_gradient", "Gradient term", "#009E73"),
        ("train_weighted_boundary", "Boundary term", "#D55E00"),
        (
            "train_weighted_residual_regularization",
            "Residual regularization",
            "#56B4E9",
        ),
        ("train_total", "Total loss", "#000000"),
    ]
    seed_maxima = [max(int(row["epoch"]) for row in history) for history in rows]
    max_epoch = min(seed_maxima)
    fig, ax = plt.subplots(figsize=(7.2, 4.05))
    fig.subplots_adjust(left=0.10, right=0.98, bottom=0.16, top=0.90)
    figure_rows: List[Dict[str, Any]] = []
    plotted_values: List[float] = []
    for key, label, color in components:
        epochs: List[int] = []
        means: List[float] = []
        stds: List[float] = []
        for epoch in range(1, max_epoch + 1):
            values = [
                float(row[key])
                for history in rows
                for row in history
                if int(row["epoch"]) == epoch and key in row
            ]
            if len(values) != len(rows):
                continue
            epochs.append(epoch)
            means.append(float(np.mean(values)))
            stds.append(float(np.std(values, ddof=1)) if len(values) > 1 else 0.0)
            figure_rows.append(
                {
                    "component": key,
                    "label": label,
                    "epoch": epoch,
                    "mean": means[-1],
                    "std": stds[-1],
                    "contributing_seeds": len(values),
                }
            )
        if not epochs:
            continue
        plot_values = np.maximum(np.asarray(means), 1e-12)
        plotted_values.extend(plot_values.tolist())
        ax.plot(epochs, plot_values, label=label, color=color)
        lower = np.maximum(plot_values - np.asarray(stds), 1e-12)
        upper = plot_values + np.asarray(stds)
        ax.fill_between(epochs, lower, upper, color=color, alpha=0.10, linewidth=0)
    warmup_epoch = int(CONFIG["physics_warmup_epochs"])
    ax.axvline(
        warmup_epoch,
        color="#F0E442",
        linestyle="--",
        linewidth=1.0,
        label=f"Curriculum end ({warmup_epoch} epochs)",
    )
    if plotted_values:
        positive_min = max(min(plotted_values), 1e-12)
        positive_max = max(plotted_values)
        ax.set_ylim(max(positive_min * 0.015, 1e-12), positive_max * 1.8)
    ax.set_yscale("log")
    ax.set_xlabel("Epoch")
    ax.set_ylabel("Weighted loss contribution")
    ax.set_title("PCFNO loss decomposition", pad=8)
    ax.grid(True, which="both", alpha=0.2, linewidth=0.6)
    ax.legend(
        loc="lower center",
        bbox_to_anchor=(0.5, 0.015),
        ncol=3,
        frameon=False,
        columnspacing=1.15,
        handlelength=2.0,
        fontsize=7.4,
    )
    stem = "fig02_pcfno_loss_decomposition"
    save_dat(output_dir / "figure_data" / f"{stem}.dat", figure_rows)
    save_figure(fig, output_dir / "figures" / stem)


def plot_metric_bars(
    rows: Sequence[Dict[str, Any]],
    metrics: Sequence[Tuple[str, str]],
    title: str,
    stem: str,
    output_dir: Path,
) -> None:
    if not rows:
        return
    fig, axes = plt.subplots(1, len(metrics), figsize=(7.2, 3.45))
    if len(metrics) == 1:
        axes = [axes]
    labels = [DISPLAY_NAMES[str(row["model"])] for row in rows]
    colors = [MODEL_COLORS[str(row["model"])] for row in rows]
    for ax, (metric, ylabel) in zip(axes, metrics):
        means = [float(row.get(f"{metric}_mean", float("nan"))) for row in rows]
        stds = [float(row.get(f"{metric}_std", 0.0)) for row in rows]
        positions = np.arange(len(rows))
        ax.bar(
            positions,
            means,
            yerr=stds,
            color=colors,
            edgecolor="black",
            linewidth=0.5,
            capsize=2.5,
            alpha=0.90,
        )
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.set_ylabel(ylabel)
        ax.grid(True, axis="y", alpha=0.2, linewidth=0.6)
        ax.margins(x=0.08)
    fig.suptitle(title, y=0.965, fontsize=11.0)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.27, top=0.80, wspace=0.34)
    save_dat(output_dir / "figure_data" / f"{stem}.dat", list(rows))
    save_figure(fig, output_dir / "figures" / stem)


def plot_ablation_bars(rows: Sequence[Dict[str, Any]], output_dir: Path) -> None:
    if not rows:
        return
    metrics = [
        ("Relative_L2", r"Relative $L_2$ error"),
        ("Gradient_RMSE", "Gradient RMSE"),
        ("Boundary_RMSE_over_Uin", r"Boundary RMSE/$U_{in}$"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 4.35), constrained_layout=True)
    labels = [DISPLAY_NAMES[str(row["model"])] for row in rows]
    positions = np.arange(len(rows))
    colors = [MODEL_COLORS[str(row["model"])] for row in rows]
    for index, (ax, (metric, xlabel)) in enumerate(zip(axes, metrics)):
        means = [float(row.get(f"{metric}_mean", float("nan"))) for row in rows]
        stds = [float(row.get(f"{metric}_std", 0.0)) for row in rows]
        ax.barh(
            positions,
            means,
            xerr=stds,
            color=colors,
            edgecolor="black",
            linewidth=0.5,
            capsize=2.2,
            alpha=0.90,
        )
        ax.set_xlabel(xlabel)
        ax.grid(True, axis="x", alpha=0.2, linewidth=0.6)
        ax.invert_yaxis()
        if index == 0:
            ax.set_yticks(positions)
            ax.set_yticklabels(labels)
        else:
            ax.set_yticks(positions)
            ax.set_yticklabels([])
    fig.suptitle("PCFNO ablation study", y=1.04)
    stem = "fig04_ablation_metrics"
    save_dat(output_dir / "figure_data" / f"{stem}.dat", list(rows))
    save_figure(fig, output_dir / "figures" / stem)


def plot_runtime_comparison(rows: Sequence[Dict[str, Any]], output_dir: Path) -> None:
    if not rows:
        return
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.25), constrained_layout=True)
    labels = [DISPLAY_NAMES[str(row["model"])] for row in rows]
    colors = [MODEL_COLORS[str(row["model"])] for row in rows]
    positions = np.arange(len(rows))
    train_means = [
        float(
            row.get(
                "effective_training_time_seconds_mean",
                row.get("training_time_seconds_mean", float("nan")),
            )
        )
        / 60.0
        for row in rows
    ]
    train_stds = [
        float(row.get("training_time_seconds_std", 0.0)) / 60.0 for row in rows
    ]
    infer_means = [
        float(row.get("Inference_ms_per_frame_mean", float("nan"))) for row in rows
    ]
    infer_stds = [float(row.get("Inference_ms_per_frame_std", 0.0)) for row in rows]
    axes[0].bar(
        positions,
        train_means,
        yerr=train_stds,
        color=colors,
        edgecolor="black",
        linewidth=0.5,
        capsize=2.5,
    )
    axes[0].set_ylabel("Training time (min)")
    axes[1].bar(
        positions,
        infer_means,
        yerr=infer_stds,
        color=colors,
        edgecolor="black",
        linewidth=0.5,
        capsize=2.5,
    )
    axes[1].set_ylabel("Inference time (ms frame$^{-1}$)")
    for ax in axes:
        ax.set_xticks(positions)
        ax.set_xticklabels(labels, rotation=35, ha="right")
        ax.grid(True, axis="y", alpha=0.2, linewidth=0.6)
    stem = "fig05_runtime_comparison"
    save_dat(output_dir / "figure_data" / f"{stem}.dat", list(rows))
    save_figure(fig, output_dir / "figures" / stem)


def plot_baseline_field_comparison(
    model_records: Dict[str, Dict[str, Any]], output_dir: Path
) -> None:
    if "PCFNO" not in model_records:
        return
    reference = model_records["PCFNO"]
    inputs = reference["inputs"]
    target_raw = reference["target"][0]
    fluid = inputs[INPUT_SOLID] < 0.5
    extent = field_extent(inputs)
    vmax = max(float(np.percentile(target_raw[fluid], 99.5)), 1e-6)
    names = [name for name in BASELINE_NAMES if name in model_records]
    error_raw = {
        name: np.abs(model_records[name]["prediction"][0] - target_raw)
        for name in names
    }
    error_values = np.concatenate([error_raw[name][fluid] for name in names])
    error_vmax = max(float(np.percentile(error_values, 99.0)), 1e-6)
    fig = plt.figure(figsize=(9.0, 4.75))
    grid = fig.add_gridspec(
        4,
        5,
        height_ratios=[1.0, 0.075, 1.0, 0.075],
        left=0.055,
        right=0.985,
        bottom=0.09,
        top=0.86,
        wspace=0.18,
        hspace=0.38,
    )
    top_axes = [fig.add_subplot(grid[0, column]) for column in range(5)]
    top_fields = [("CFD", target_raw)] + [
        (DISPLAY_NAMES[name], model_records[name]["prediction"][0]) for name in names
    ]
    speed_image = None
    for column, (title, field) in enumerate(top_fields):
        ax = top_axes[column]
        speed_image = ax.imshow(
            plot_ready_field(field, inputs),
            origin="lower",
            extent=extent,
            cmap=solid_aware_cmap("viridis"),
            vmin=0,
            vmax=vmax,
            aspect="equal",
        )
        contour_solid(ax, inputs)
        ax.set_title(title, pad=4)
        ax.set_xlabel(r"$(x-x_c)/D$", labelpad=2)
        if column == 0:
            ax.set_ylabel(r"$(y-y_c)/D$")
        else:
            ax.set_yticklabels([])
    cax_speed = fig.add_subplot(grid[1, :])
    speed_bar = fig.colorbar(
        speed_image,
        cax=cax_speed,
        orientation="horizontal",
    )
    speed_bar.set_label(r"$|\mathbf{U}|/U_{in}$", labelpad=2)
    cax_speed.xaxis.set_label_position("top")
    blank = fig.add_subplot(grid[2, 0])
    blank.axis("off")
    error_axes = [fig.add_subplot(grid[2, column]) for column in range(1, 5)]
    error_image = None
    rows: List[Dict[str, Any]] = []
    x = inputs[INPUT_XD]
    y = inputs[INPUT_YD]
    for ax, name in zip(error_axes, names):
        error = error_raw[name]
        error_image = ax.imshow(
            plot_ready_field(error, inputs),
            origin="lower",
            extent=extent,
            cmap=solid_aware_cmap("magma"),
            vmin=0,
            vmax=error_vmax,
            aspect="equal",
        )
        contour_solid(ax, inputs)
        ax.set_title(f"{DISPLAY_NAMES[name]} absolute error", pad=4)
        ax.set_xlabel(r"$(x-x_c)/D$", labelpad=2)
        ax.set_yticklabels([])
        for xd, yd, truth, pred, err, solid_value in zip(
            x.ravel(),
            y.ravel(),
            target_raw.ravel(),
            model_records[name]["prediction"][0].ravel(),
            error.ravel(),
            inputs[INPUT_SOLID].ravel(),
        ):
            rows.append(
                {
                    "model": name,
                    "x_over_D": xd,
                    "y_over_D": yd,
                    "solid": solid_value,
                    "CFD_speed_over_Uin": truth,
                    "prediction_over_Uin": pred,
                    "absolute_error_over_Uin": err,
                }
            )
    cax_error = fig.add_subplot(grid[3, 1:])
    fig.colorbar(
        error_image,
        cax=cax_error,
        orientation="horizontal",
        label=r"Absolute error in $|\mathbf{U}|/U_{in}$",
    )
    metadata = reference["metadata"]
    fig.suptitle(
        f"{metadata['case']}: Re={metadata['Re']:.1f}, $t^*$={metadata['time_star']:.2f}",
        y=0.965,
        fontsize=11.0,
    )
    stem = "fig06_baseline_field_comparison"
    save_dat(output_dir / "figure_data" / f"{stem}.dat", rows)
    save_figure(fig, output_dir / "figures" / stem)


def plot_three_reynolds_cases(
    records: Sequence[Dict[str, Any]], output_dir: Path
) -> None:
    if not records:
        return
    fluid_values = np.concatenate(
        [record["target"][0][record["inputs"][INPUT_SOLID] < 0.5] for record in records]
    )
    vmax = max(float(np.percentile(fluid_values, 99.5)), 1e-6)
    errors = [
        np.abs(record["prediction"][0] - record["target"][0]) for record in records
    ]
    error_values = np.concatenate(
        [
            error[record["inputs"][INPUT_SOLID] < 0.5]
            for error, record in zip(errors, records)
        ]
    )
    error_vmax = max(float(np.percentile(error_values, 99.0)), 1e-6)
    nrows = len(records)
    fig = plt.figure(figsize=(7.2, 1.62 * nrows + 0.92))
    grid = fig.add_gridspec(
        nrows + 2,
        3,
        height_ratios=[1.0] * nrows + [0.12, 0.08],
        left=0.09,
        right=0.985,
        bottom=0.07,
        top=0.93,
        wspace=0.14,
        hspace=0.14,
    )
    axes = np.asarray(
        [[fig.add_subplot(grid[row, col]) for col in range(3)] for row in range(nrows)]
    )
    rows: List[Dict[str, Any]] = []
    speed_image = None
    error_image = None
    for row_index, record in enumerate(records):
        inputs = record["inputs"]
        target = record["target"][0]
        pred = record["prediction"][0]
        error = np.abs(pred - target)
        extent = field_extent(inputs)
        fields = [target, pred, error]
        titles = ["CFD", "PCFNO", "Absolute error"]
        cmaps = ["viridis", "viridis", "magma"]
        limits = [(0, vmax), (0, vmax), (0, error_vmax)]
        for col, (field, title, cmap, limits_pair) in enumerate(
            zip(fields, titles, cmaps, limits)
        ):
            image = axes[row_index, col].imshow(
                plot_ready_field(field, inputs),
                origin="lower",
                extent=extent,
                cmap=solid_aware_cmap(cmap),
                vmin=limits_pair[0],
                vmax=limits_pair[1],
                aspect="equal",
            )
            if col < 2:
                speed_image = image
            else:
                error_image = image
            contour_solid(axes[row_index, col], inputs)
            axes[row_index, col].set_title(title if row_index == 0 else "", pad=4)
            if row_index == nrows - 1:
                axes[row_index, col].set_xlabel(r"$(x-x_c)/D$", labelpad=2)
            else:
                axes[row_index, col].set_xticklabels([])
            if col == 0:
                axes[row_index, col].set_ylabel(
                    rf"Re={record['metadata']['Re']:.1f}" + "\n" + r"$(y-y_c)/D$"
                )
            else:
                axes[row_index, col].set_yticklabels([])
        x = inputs[INPUT_XD]
        y = inputs[INPUT_YD]
        for xd, yd, truth, estimate, err, solid_value in zip(
            x.ravel(),
            y.ravel(),
            target.ravel(),
            pred.ravel(),
            error.ravel(),
            inputs[INPUT_SOLID].ravel(),
        ):
            rows.append(
                {
                    "case": record["metadata"]["case"],
                    "Re": record["metadata"]["Re"],
                    "time_star": record["metadata"]["time_star"],
                    "x_over_D": xd,
                    "y_over_D": yd,
                    "solid": solid_value,
                    "CFD_speed_over_Uin": truth,
                    "PCFNO_speed_over_Uin": estimate,
                    "absolute_error_over_Uin": err,
                }
            )
    cbar_grid = grid[nrows + 1, :].subgridspec(1, 3, wspace=0.28)
    cax_speed = fig.add_subplot(cbar_grid[0, :2])
    cax_error = fig.add_subplot(cbar_grid[0, 2])
    fig.colorbar(
        speed_image,
        cax=cax_speed,
        orientation="horizontal",
        label=r"$|\mathbf{U}|/U_{in}$",
    )
    fig.colorbar(
        error_image,
        cax=cax_error,
        orientation="horizontal",
        label="Absolute error",
    )
    stem = "fig07_three_reynolds_cases"
    save_dat(output_dir / "figure_data" / f"{stem}.dat", rows)
    save_figure(fig, output_dir / "figures" / stem)


def plot_time_evolution(records: Sequence[Dict[str, Any]], output_dir: Path) -> None:
    if not records:
        return
    fluid_values = np.concatenate(
        [record["target"][0][record["inputs"][INPUT_SOLID] < 0.5] for record in records]
    )
    vmax = max(float(np.percentile(fluid_values, 99.5)), 1e-6)
    errors = [
        np.abs(record["prediction"][0] - record["target"][0]) for record in records
    ]
    error_values = np.concatenate(
        [
            error[record["inputs"][INPUT_SOLID] < 0.5]
            for error, record in zip(errors, records)
        ]
    )
    error_vmax = max(float(np.percentile(error_values, 99.0)), 1e-6)
    nrows = len(records)
    fig = plt.figure(figsize=(7.2, 1.62 * nrows + 0.92))
    grid = fig.add_gridspec(
        nrows + 2,
        3,
        height_ratios=[1.0] * nrows + [0.12, 0.08],
        left=0.09,
        right=0.985,
        bottom=0.07,
        top=0.93,
        wspace=0.14,
        hspace=0.14,
    )
    axes = np.asarray(
        [[fig.add_subplot(grid[row, col]) for col in range(3)] for row in range(nrows)]
    )
    rows: List[Dict[str, Any]] = []
    speed_image = None
    error_image = None
    for row_index, record in enumerate(records):
        inputs = record["inputs"]
        target = record["target"][0]
        pred = record["prediction"][0]
        error = np.abs(pred - target)
        extent = field_extent(inputs)
        for col, (field, cmap, vmin, vmax_local) in enumerate(
            [
                (target, "viridis", 0, vmax),
                (pred, "viridis", 0, vmax),
                (error, "magma", 0, error_vmax),
            ]
        ):
            image = axes[row_index, col].imshow(
                plot_ready_field(field, inputs),
                origin="lower",
                extent=extent,
                cmap=solid_aware_cmap(cmap),
                vmin=vmin,
                vmax=vmax_local,
                aspect="equal",
            )
            if col < 2:
                speed_image = image
            else:
                error_image = image
            contour_solid(axes[row_index, col], inputs)
            if row_index == nrows - 1:
                axes[row_index, col].set_xlabel(r"$(x-x_c)/D$", labelpad=2)
            else:
                axes[row_index, col].set_xticklabels([])
            if col == 0:
                axes[row_index, col].set_ylabel(
                    rf"$t^*$={record['metadata']['time_star']:.2f}"
                    + "\n"
                    + r"$(y-y_c)/D$"
                )
            else:
                axes[row_index, col].set_yticklabels([])
        if row_index == 0:
            for col, title in enumerate(["CFD", "PCFNO", "Absolute error"]):
                axes[row_index, col].set_title(title, pad=4)
        x = inputs[INPUT_XD]
        y = inputs[INPUT_YD]
        for xd, yd, truth, estimate, err, solid_value in zip(
            x.ravel(),
            y.ravel(),
            target.ravel(),
            pred.ravel(),
            error.ravel(),
            inputs[INPUT_SOLID].ravel(),
        ):
            rows.append(
                {
                    "case": record["metadata"]["case"],
                    "Re": record["metadata"]["Re"],
                    "time_star": record["metadata"]["time_star"],
                    "x_over_D": xd,
                    "y_over_D": yd,
                    "solid": solid_value,
                    "CFD_speed_over_Uin": truth,
                    "PCFNO_speed_over_Uin": estimate,
                    "absolute_error_over_Uin": err,
                }
            )
    cbar_grid = grid[nrows + 1, :].subgridspec(1, 3, wspace=0.28)
    cax_speed = fig.add_subplot(cbar_grid[0, :2])
    cax_error = fig.add_subplot(cbar_grid[0, 2])
    fig.colorbar(
        speed_image,
        cax=cax_speed,
        orientation="horizontal",
        label=r"$|\mathbf{U}|/U_{in}$",
    )
    fig.colorbar(
        error_image,
        cax=cax_error,
        orientation="horizontal",
        label="Absolute error",
    )
    stem = "fig08_time_evolution"
    save_dat(output_dir / "figure_data" / f"{stem}.dat", rows)
    save_figure(fig, output_dir / "figures" / stem)


def grouped_curve(
    rows: Sequence[Dict[str, Any]],
    key: str,
    value_key: str,
    unit_key: str = "case",
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    grouped_units: Dict[float, Dict[str, List[float]]] = {}
    for row in rows:
        x_value = float(row[key])
        unit = str(row.get(unit_key, row.get("Re", "unknown")))
        grouped_units.setdefault(x_value, {}).setdefault(unit, []).append(
            float(row[value_key])
        )

    x_values = np.asarray(sorted(grouped_units), dtype=np.float64)
    means: List[float] = []
    lowers: List[float] = []
    uppers: List[float] = []
    counts: List[int] = []
    repeats = max(int(CONFIG.get("bootstrap_repeats", 5000)), 200)
    confidence = float(CONFIG.get("confidence_level", 0.95))
    alpha = 0.5 * (1.0 - confidence)

    for group_index, x_value in enumerate(x_values):
        unit_values = np.asarray(
            [
                np.mean(values)
                for values in grouped_units[float(x_value)].values()
                if values
            ],
            dtype=np.float64,
        )
        means.append(float(np.mean(unit_values)))
        counts.append(int(unit_values.size))
        if unit_values.size <= 1:
            lowers.append(means[-1])
            uppers.append(means[-1])
            continue
        rng = np.random.default_rng(
            int(CONFIG["split_seed"]) + 7919 * (group_index + 1)
        )
        sampled_indices = rng.integers(
            0, unit_values.size, size=(repeats, unit_values.size)
        )
        bootstrap_means = np.mean(unit_values[sampled_indices], axis=1)
        lowers.append(float(np.quantile(bootstrap_means, alpha)))
        uppers.append(float(np.quantile(bootstrap_means, 1.0 - alpha)))

    return (
        x_values,
        np.asarray(means),
        np.asarray(lowers),
        np.asarray(uppers),
        np.asarray(counts, dtype=np.int64),
    )


def plot_error_vs_time_and_re(
    per_sample_by_model_seed: Dict[Tuple[str, int], List[Dict[str, Any]]],
    output_dir: Path,
) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.25), constrained_layout=True)
    figure_rows: List[Dict[str, Any]] = []
    for name in ["FNO", "PCFNO"]:
        rows = [
            row
            for (model, _), model_rows in per_sample_by_model_seed.items()
            if model == name
            for row in model_rows
        ]
        if not rows:
            continue
        for ax, key, xlabel in [
            (axes[0], "time_star", r"Dimensionless time $t^*=tU_{in}/D$"),
            (axes[1], "Re", r"Reynolds number $Re$"),
        ]:
            x, mean, lower, upper, n_units = grouped_curve(
                rows, key, "Relative_L2", unit_key="case"
            )
            color = MODEL_COLORS[name]
            ax.plot(x, mean, marker="o", label=DISPLAY_NAMES[name], color=color)
            ax.fill_between(
                x,
                np.maximum(lower, 0),
                upper,
                color=color,
                alpha=0.16,
                linewidth=0,
            )
            ax.set_xlabel(xlabel)
            ax.set_ylabel(r"Relative $L_2$ error")
            ax.grid(True, alpha=0.2, linewidth=0.6)
            for xv, mv, lv, uv, nv in zip(x, mean, lower, upper, n_units):
                figure_rows.append(
                    {
                        "model": name,
                        "group_variable": key,
                        "group_value": xv,
                        "Relative_L2_mean": mv,
                        "CI_lower": lv,
                        "CI_upper": uv,
                        "independent_case_units": int(nv),
                        "confidence_level": float(CONFIG.get("confidence_level", 0.95)),
                    }
                )
    for ax in axes:
        ax.legend()
        ax.text(
            0.98,
            0.03,
            "Shading: 95% case-cluster bootstrap CI",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=7.0,
            color="tab:gray",
        )
    stem = "fig09_error_vs_time_and_re"
    save_dat(output_dir / "figure_data" / f"{stem}.dat", figure_rows)
    save_figure(fig, output_dir / "figures" / stem)


def plot_uncertainty(output_dir: Path, seeds: Sequence[int]) -> None:
    error_arrays: List[np.ndarray] = []
    raw_arrays: List[np.ndarray] = []
    calibrated_arrays: List[np.ndarray] = []
    for seed in seeds:
        path = output_dir / "data" / "uncertainty" / f"uncertainty_PCFNO_seed{seed}.npz"
        if path.exists():
            data = np.load(path)
            error_arrays.append(data["error"])
            raw_arrays.append(data["sigma_raw"])
            calibrated_arrays.append(data["sigma_calibrated"])
    if not error_arrays:
        return
    error = np.concatenate(error_arrays)
    sigma_raw = np.concatenate(raw_arrays)
    sigma_cal = np.concatenate(calibrated_arrays)
    maximum = min(error.size, 30000)
    rng = np.random.default_rng(int(CONFIG["split_seed"]))
    selected = (
        rng.choice(error.size, maximum, replace=False)
        if error.size > maximum
        else np.arange(error.size)
    )
    selected_error = error[selected]
    selected_sigma = sigma_cal[selected]
    percentile = float(CONFIG.get("uncertainty_plot_percentile", 99.5))
    limit = max(
        float(np.percentile(selected_sigma, percentile)),
        float(np.percentile(selected_error, percentile)),
        1e-8,
    )
    visible = (
        np.isfinite(selected_sigma)
        & np.isfinite(selected_error)
        & (selected_sigma >= 0.0)
        & (selected_error >= 0.0)
        & (selected_sigma <= limit)
        & (selected_error <= limit)
    )
    visible_sigma = selected_sigma[visible]
    visible_error = selected_error[visible]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.35), constrained_layout=True)
    density = axes[0].hexbin(
        visible_sigma,
        visible_error,
        gridsize=int(CONFIG.get("uncertainty_plot_gridsize", 50)),
        mincnt=1,
        bins="log",
        cmap="viridis",
        extent=(0.0, limit, 0.0, limit),
        linewidths=0.0,
    )
    colorbar = fig.colorbar(density, ax=axes[0], pad=0.02, fraction=0.05)
    colorbar.set_label("Count per hexagon", fontsize=7.6)
    colorbar.ax.tick_params(labelsize=6.8)
    axes[0].plot([0, limit], [0, limit], color="#000000", linestyle="--", lw=0.9)
    bin_edges = np.linspace(0.0, limit, 16)
    bin_centers: List[float] = []
    medians: List[float] = []
    lower_quantiles: List[float] = []
    upper_quantiles: List[float] = []
    for left, right in zip(bin_edges[:-1], bin_edges[1:]):
        mask = (visible_sigma >= left) & (visible_sigma < right)
        if np.count_nonzero(mask) < 20:
            continue
        values = visible_error[mask]
        bin_centers.append(0.5 * (left + right))
        medians.append(float(np.median(values)))
        lower_quantiles.append(float(np.quantile(values, 0.10)))
        upper_quantiles.append(float(np.quantile(values, 0.90)))
    if bin_centers:
        axes[0].fill_between(
            bin_centers,
            lower_quantiles,
            upper_quantiles,
            color="#D9D9D9",
            alpha=0.45,
            linewidth=0.0,
            label="10–90% interval",
        )
        axes[0].plot(
            bin_centers,
            medians,
            color="#D55E00",
            linewidth=1.25,
            label="Binned median",
        )
        axes[0].legend(loc="upper left", fontsize=7.0)
    axes[0].set_xlim(0, limit)
    axes[0].set_ylim(0, limit)
    axes[0].set_xlabel(r"Calibrated $\sigma/U_{in}$")
    axes[0].set_ylabel(r"Absolute error in $|\mathbf{U}|/U_{in}$")
    multipliers = np.linspace(0.25, 3.0, 40)
    theoretical = np.asarray(
        [math.erf(multiplier / math.sqrt(2.0)) for multiplier in multipliers]
    )
    empirical_raw = np.asarray(
        [np.mean(error <= multiplier * sigma_raw) for multiplier in multipliers]
    )
    empirical_cal = np.asarray(
        [np.mean(error <= multiplier * sigma_cal) for multiplier in multipliers]
    )
    axes[1].plot(
        theoretical,
        empirical_raw,
        marker="o",
        markersize=2.8,
        linewidth=1.1,
        label="Raw",
        color="#000000",
    )
    axes[1].plot(
        theoretical,
        empirical_cal,
        marker="o",
        markersize=2.8,
        linewidth=1.5,
        label="Validation-calibrated",
        color="#D55E00",
    )
    axes[1].plot([0, 1], [0, 1], color="#0072B2", linestyle="--", lw=0.9, label="Ideal")
    axes[1].set_xlabel("Nominal Gaussian coverage")
    axes[1].set_ylabel("Empirical coverage")
    axes[1].set_xlim(0, 1)
    axes[1].set_ylim(0, 1)
    axes[1].grid(True, alpha=0.2, linewidth=0.6)
    axes[1].legend(loc="upper left")
    rows: List[Dict[str, Any]] = []
    for sigma_value, error_value in zip(selected_sigma, selected_error):
        rows.append(
            {
                "record_type": "scatter_calibrated",
                "sigma_over_Uin": sigma_value,
                "absolute_error_over_Uin": error_value,
                "visible_in_density_panel": bool(
                    np.isfinite(sigma_value)
                    and np.isfinite(error_value)
                    and 0.0 <= sigma_value <= limit
                    and 0.0 <= error_value <= limit
                ),
                "plot_percentile": percentile,
                "plot_limit": limit,
            }
        )
    for center, median, lower, upper in zip(
        bin_centers, medians, lower_quantiles, upper_quantiles
    ):
        rows.append(
            {
                "record_type": "binned_error",
                "sigma_over_Uin": center,
                "median_absolute_error": median,
                "q10_absolute_error": lower,
                "q90_absolute_error": upper,
            }
        )
    for multiplier, nominal, raw, calibrated in zip(
        multipliers, theoretical, empirical_raw, empirical_cal
    ):
        rows.append(
            {
                "record_type": "calibration",
                "sigma_multiplier": multiplier,
                "nominal_coverage": nominal,
                "raw_empirical_coverage": raw,
                "calibrated_empirical_coverage": calibrated,
            }
        )
    stem = "fig10_uncertainty_calibration"
    save_dat(output_dir / "figure_data" / f"{stem}.dat", rows)
    save_figure(fig, output_dir / "figures" / stem)


def nearest_grid_index(
    inputs: np.ndarray, x_over_d: float, y_over_d: float
) -> Tuple[int, int]:
    x = inputs[INPUT_XD, 0, :]
    y = inputs[INPUT_YD, :, 0]
    return int(np.argmin(np.abs(y - y_over_d))), int(np.argmin(np.abs(x - x_over_d)))


def dominant_dimensionless_frequency(
    times: np.ndarray, signal: np.ndarray, diameter: float, uin: float
) -> Tuple[float, np.ndarray, np.ndarray]:
    if times.size < 16:
        return float("nan"), np.array([]), np.array([])
    start = times.size // 2
    time_segment = times[start:]
    signal_segment = detrend(signal[start:], type="linear")
    dt = float(np.median(np.diff(time_segment)))
    if dt <= 0:
        return float("nan"), np.array([]), np.array([])
    nperseg = min(256, max(16, signal_segment.size // 2))
    frequencies, power = welch(
        signal_segment, fs=1.0 / dt, nperseg=nperseg, detrend="linear"
    )
    positive = frequencies > 0
    frequencies = frequencies[positive]
    power = power[positive]
    if power.size == 0:
        return float("nan"), frequencies, power
    dominant = float(frequencies[int(np.argmax(power))] * diameter / max(uin, 1e-12))
    return dominant, frequencies * diameter / max(uin, 1e-12), power


def plot_probe_and_spectrum(
    records: Sequence[Dict[str, Any]], output_dir: Path
) -> None:
    if not records:
        return
    records = sorted(records, key=lambda record: float(record["metadata"]["time"]))
    reference = records[0]
    diameter = float(reference["metadata"]["diameter"])
    uin = float(reference["metadata"]["Uin"])
    iy, ix = nearest_grid_index(reference["inputs"], 2.0, 0.5)
    times = np.asarray([float(record["metadata"]["time"]) for record in records])
    time_star = np.asarray(
        [float(record["metadata"]["time_star"]) for record in records]
    )
    true_signal = np.asarray([record["target"][0, iy, ix] for record in records])
    pred_signal = np.asarray([record["prediction"][0, iy, ix] for record in records])
    true_frequency, true_axis, true_power = dominant_dimensionless_frequency(
        times, true_signal, diameter, uin
    )
    pred_frequency, pred_axis, pred_power = dominant_dimensionless_frequency(
        times, pred_signal, diameter, uin
    )
    true_power_n = (
        true_power / max(float(np.max(true_power)), 1e-12)
        if true_power.size
        else true_power
    )
    pred_power_n = (
        pred_power / max(float(np.max(pred_power)), 1e-12)
        if pred_power.size
        else pred_power
    )
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.25), constrained_layout=True)
    axes[0].plot(time_star, true_signal, label="CFD", color="tab:blue")
    axes[0].plot(
        time_star, pred_signal, "--", label="PCFNO", color=MODEL_COLORS["PCFNO"]
    )
    axes[0].set_xlabel(r"Dimensionless time $t^*=tU_{in}/D$")
    axes[0].set_ylabel(r"Probe $|\mathbf{U}|/U_{in}$")
    axes[0].grid(True, alpha=0.2, linewidth=0.6)
    axes[0].legend()
    if true_axis.size:
        axes[1].plot(
            true_axis,
            true_power_n,
            label=f"CFD, $St_s$={true_frequency:.3f}",
            color="tab:blue",
        )
        axes[1].plot(
            pred_axis,
            pred_power_n,
            "--",
            label=f"PCFNO, $St_s$={pred_frequency:.3f}",
            color=MODEL_COLORS["PCFNO"],
        )
    axes[1].set_xlabel(r"Dimensionless frequency $St_s=fD/U_{in}$")
    axes[1].set_ylabel("Normalized spectral power")
    axes[1].grid(True, alpha=0.2, linewidth=0.6)
    if axes[1].lines:
        axes[1].legend()
    rows: List[Dict[str, Any]] = []
    for t, ts, truth, pred in zip(times, time_star, true_signal, pred_signal):
        rows.append(
            {
                "record_type": "time_series",
                "time_s": t,
                "time_star": ts,
                "CFD_speed_over_Uin": truth,
                "PCFNO_speed_over_Uin": pred,
                "probe_x_over_D": 2.0,
                "probe_y_over_D": 0.5,
            }
        )
    for index in range(max(len(true_axis), len(pred_axis))):
        rows.append(
            {
                "record_type": "spectrum",
                "CFD_St_axis": (
                    true_axis[index] if index < len(true_axis) else float("nan")
                ),
                "CFD_normalized_power": (
                    true_power_n[index] if index < len(true_power_n) else float("nan")
                ),
                "PCFNO_St_axis": (
                    pred_axis[index] if index < len(pred_axis) else float("nan")
                ),
                "PCFNO_normalized_power": (
                    pred_power_n[index] if index < len(pred_power_n) else float("nan")
                ),
            }
        )
    stem = "fig11_speed_probe_spectrum"
    save_dat(output_dir / "figure_data" / f"{stem}.dat", rows)
    save_figure(fig, output_dir / "figures" / stem)
    save_dat(
        output_dir / "tables" / "table09_speed_probe_frequency.dat",
        [
            {
                "case": reference["metadata"]["case"],
                "Re": reference["metadata"]["Re"],
                "CFD_speed_probe_St": true_frequency,
                "PCFNO_speed_probe_St": pred_frequency,
                "absolute_difference": (
                    abs(pred_frequency - true_frequency)
                    if np.isfinite(true_frequency) and np.isfinite(pred_frequency)
                    else float("nan")
                ),
                "note": "Speed-magnitude probe frequency; it may emphasize harmonics of vortex shedding.",
            }
        ],
    )


def plot_speed_profiles(record: Dict[str, Any], output_dir: Path) -> None:
    inputs = record["inputs"]
    target = record["target"][0]
    prediction = record["prediction"][0]
    x = inputs[INPUT_XD, 0, :]
    y = inputs[INPUT_YD, :, 0]
    locations = [1.5, 2.5, 3.5]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 3.15), sharey=True)
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.17, top=0.78, wspace=0.11)
    rows: List[Dict[str, Any]] = []
    for ax, location in zip(axes, locations):
        ix = int(np.argmin(np.abs(x - location)))
        ax.plot(target[:, ix], y, label="CFD", color="tab:blue")
        ax.plot(prediction[:, ix], y, "--", label="PCFNO", color=MODEL_COLORS["PCFNO"])
        ax.set_title(rf"$(x-x_c)/D={x[ix]:.2f}$", pad=6)
        ax.set_xlabel(r"$|\mathbf{U}|/U_{in}$")
        ax.grid(True, alpha=0.2, linewidth=0.6)
        for y_value, truth, pred in zip(y, target[:, ix], prediction[:, ix]):
            rows.append(
                {
                    "x_over_D": x[ix],
                    "y_over_D": y_value,
                    "CFD_speed_over_Uin": truth,
                    "PCFNO_speed_over_Uin": pred,
                }
            )
    axes[0].set_ylabel(r"$(y-y_c)/D$")
    axes[0].legend(loc="center right")
    metadata = record["metadata"]
    fig.suptitle(
        f"{metadata['case']}: Re={metadata['Re']:.1f}, $t^*$={metadata['time_star']:.2f}",
        y=0.965,
        fontsize=10.5,
    )
    stem = "fig12_speed_profiles"
    save_dat(output_dir / "figure_data" / f"{stem}.dat", rows)
    save_figure(fig, output_dir / "figures" / stem)


def aggregate_three_re_metrics(
    records: Sequence[Dict[str, Any]], model_name: str, seed: int
) -> List[Dict[str, Any]]:
    rows: List[Dict[str, Any]] = []
    for record in records:
        target = record["target"][0]
        prediction = record["prediction"][0]
        fluid = record["inputs"][INPUT_SOLID] < 0.5
        error = prediction[fluid] - target[fluid]
        true = target[fluid]
        rows.append(
            {
                "model": model_name,
                "seed": seed,
                "case": record["metadata"]["case"],
                "Re": record["metadata"]["Re"],
                "time_star": record["metadata"]["time_star"],
                "Relative_L2": float(
                    np.linalg.norm(error) / max(np.linalg.norm(true), 1e-12)
                ),
                "RMSE_over_Uin": float(np.sqrt(np.mean(error**2))),
                "MAE_over_Uin": float(np.mean(np.abs(error))),
            }
        )
    return rows


def paired_fno_pcfno_statistics(
    per_sample_by_model_seed: Dict[Tuple[str, int], List[Dict[str, Any]]]
) -> List[Dict[str, Any]]:
    output_rows: List[Dict[str, Any]] = []
    common_seeds = sorted(
        {
            seed
            for model, seed in per_sample_by_model_seed
            if model == "FNO" and ("PCFNO", seed) in per_sample_by_model_seed
        }
    )
    metrics = ["Relative_L2", "Wake_Relative_L2", "Near_wall_Relative_L2"]
    repeats = max(int(CONFIG.get("bootstrap_repeats", 5000)), 200)
    confidence = float(CONFIG.get("confidence_level", 0.95))
    alpha = 0.5 * (1.0 - confidence)

    def summarize_pairs(
        scope: str,
        analysis_unit: str,
        metric: str,
        pairs: Sequence[Tuple[float, float]],
        rng_seed: int,
    ) -> None:
        if not pairs:
            return
        fno_values = np.asarray([pair[0] for pair in pairs], dtype=np.float64)
        pcfno_values = np.asarray([pair[1] for pair in pairs], dtype=np.float64)
        differences = fno_values - pcfno_values
        if len(differences) >= 5 and not np.allclose(differences, 0.0):
            try:
                p_value = float(
                    wilcoxon(
                        differences,
                        alternative="greater",
                        zero_method="wilcox",
                    ).pvalue
                )
            except ValueError:
                p_value = float("nan")
        else:
            p_value = float("nan")

        if len(differences) > 1:
            rng = np.random.default_rng(rng_seed)
            indices = rng.integers(
                0, len(differences), size=(repeats, len(differences))
            )
            bootstrap_means = np.mean(differences[indices], axis=1)
            ci_lower = float(np.quantile(bootstrap_means, alpha))
            ci_upper = float(np.quantile(bootstrap_means, 1.0 - alpha))
        else:
            ci_lower = ci_upper = float(differences[0])

        output_rows.append(
            {
                "scope": scope,
                "analysis_unit": analysis_unit,
                "metric": metric,
                "paired_units": len(pairs),
                "FNO_mean": float(np.mean(fno_values)),
                "PCFNO_mean": float(np.mean(pcfno_values)),
                "mean_paired_improvement": float(np.mean(differences)),
                "improvement_CI_lower": ci_lower,
                "improvement_CI_upper": ci_upper,
                "confidence_level": confidence,
                "relative_improvement_percent": 100.0
                * float(np.mean(differences))
                / max(float(np.mean(fno_values)), 1e-12),
                "PCFNO_win_fraction": float(np.mean(differences > 0)),
                "wilcoxon_one_sided_p": p_value,
            }
        )

    all_frame_pairs: Dict[str, Dict[Tuple[str, int], List[Tuple[float, float]]]] = {
        metric: {} for metric in metrics
    }
    for seed in common_seeds:
        fno_map = {
            (str(row["case"]), int(row["frame_index"])): row
            for row in per_sample_by_model_seed[("FNO", seed)]
        }
        pcfno_map = {
            (str(row["case"]), int(row["frame_index"])): row
            for row in per_sample_by_model_seed[("PCFNO", seed)]
        }
        keys = sorted(set(fno_map) & set(pcfno_map))
        for metric_index, metric in enumerate(metrics):
            frame_pairs: List[Tuple[float, float]] = []
            case_groups: Dict[str, List[Tuple[float, float]]] = {}
            re_groups: Dict[float, List[Tuple[float, float]]] = {}
            for key in keys:
                fno_row, pcfno_row = fno_map[key], pcfno_map[key]
                if metric not in fno_row or metric not in pcfno_row:
                    continue
                pair = (float(fno_row[metric]), float(pcfno_row[metric]))
                if not np.all(np.isfinite(pair)):
                    continue
                frame_pairs.append(pair)
                case_groups.setdefault(str(key[0]), []).append(pair)
                re_groups.setdefault(canonical_re(float(fno_row["Re"])), []).append(
                    pair
                )
                all_frame_pairs[metric].setdefault(key, []).append(pair)

            summarize_pairs(
                f"seed_{seed}",
                "frame",
                metric,
                frame_pairs,
                int(CONFIG["split_seed"]) + seed * 101 + metric_index,
            )
            case_pairs = [
                (
                    float(np.mean([pair[0] for pair in values])),
                    float(np.mean([pair[1] for pair in values])),
                )
                for values in case_groups.values()
            ]
            summarize_pairs(
                f"seed_{seed}",
                "case",
                metric,
                case_pairs,
                int(CONFIG["split_seed"]) + seed * 211 + metric_index,
            )
            re_pairs = [
                (
                    float(np.mean([pair[0] for pair in values])),
                    float(np.mean([pair[1] for pair in values])),
                )
                for values in re_groups.values()
            ]
            summarize_pairs(
                f"seed_{seed}",
                "Re_group",
                metric,
                re_pairs,
                int(CONFIG["split_seed"]) + seed * 307 + metric_index,
            )

    for metric_index, metric in enumerate(metrics):
        frame_pair_map = all_frame_pairs[metric]
        averaged_frames: Dict[Tuple[str, int], Tuple[float, float]] = {
            key: (
                float(np.mean([pair[0] for pair in values])),
                float(np.mean([pair[1] for pair in values])),
            )
            for key, values in frame_pair_map.items()
        }
        summarize_pairs(
            "seed_averaged",
            "frame",
            metric,
            list(averaged_frames.values()),
            int(CONFIG["split_seed"]) + 10000 + metric_index,
        )

        case_groups: Dict[str, List[Tuple[float, float]]] = {}
        re_groups: Dict[float, List[Tuple[float, float]]] = {}
        reference_rows = (
            per_sample_by_model_seed[("FNO", common_seeds[0])] if common_seeds else []
        )
        re_lookup = {
            (str(row["case"]), int(row["frame_index"])): canonical_re(float(row["Re"]))
            for row in reference_rows
        }
        for key, pair in averaged_frames.items():
            case_groups.setdefault(key[0], []).append(pair)
            if key in re_lookup:
                re_groups.setdefault(re_lookup[key], []).append(pair)

        case_pairs = [
            (
                float(np.mean([pair[0] for pair in values])),
                float(np.mean([pair[1] for pair in values])),
            )
            for values in case_groups.values()
        ]
        summarize_pairs(
            "seed_averaged",
            "case",
            metric,
            case_pairs,
            int(CONFIG["split_seed"]) + 20000 + metric_index,
        )

        re_pairs = [
            (
                float(np.mean([pair[0] for pair in values])),
                float(np.mean([pair[1] for pair in values])),
            )
            for values in re_groups.values()
        ]
        summarize_pairs(
            "seed_averaged",
            "Re_group",
            metric,
            re_pairs,
            int(CONFIG["split_seed"]) + 30000 + metric_index,
        )

    return output_rows


def generate_publication_outputs(
    test_dataset: CFDBenchSpeedDataset,
    representative_cases: Sequence[CaseInfo],
    histories: Dict[Tuple[str, int], List[Dict[str, Any]]],
    aggregate_summaries: List[Dict[str, Any]],
    per_sample_by_model_seed: Dict[Tuple[str, int], List[Dict[str, Any]]],
    output_dir: Path,
    logger: logging.Logger,
) -> None:
    seeds = [int(seed) for seed in CONFIG["training_seeds"]]
    plot_training_convergence(histories, output_dir)
    plot_pcfno_loss_decomposition(histories, output_dir)
    plot_metric_bars(
        rows_for_names(aggregate_summaries, BASELINE_NAMES),
        [
            ("Relative_L2", r"Global relative $L_2$"),
            ("Wake_Relative_L2", r"Wake relative $L_2$"),
            ("Near_wall_Relative_L2", r"Near-wall relative $L_2$"),
        ],
        "Baseline comparison",
        "fig03_baseline_metrics",
        output_dir,
    )
    plot_ablation_bars(
        rows_for_names(aggregate_summaries, ["PCFNO"] + ABLATION_NAMES),
        output_dir,
    )
    plot_runtime_comparison(
        rows_for_names(aggregate_summaries, BASELINE_NAMES), output_dir
    )
    chosen_seeds = {
        name: best_seed_for_model(name, seeds, output_dir)
        for name in BASELINE_NAMES
        if checkpoint_path(output_dir, name, seeds[0]).parent.exists()
    }
    save_dat(
        output_dir / "data" / "figure_seed_selection.dat",
        [{"model": name, "seed": seed} for name, seed in chosen_seeds.items()],
    )
    middle_case = sorted(representative_cases, key=lambda case: case.reynolds)[
        len(representative_cases) // 2
    ]
    middle_indices = dataset_indices_for_case(test_dataset, middle_case.name)
    representative_index = middle_indices[len(middle_indices) // 2]
    model_records: Dict[str, Dict[str, Any]] = {}
    for name in BASELINE_NAMES:
        if name not in chosen_seeds:
            continue
        model = load_trained_model(name, chosen_seeds[name], output_dir)
        model_records[name] = predict_indices(
            model, test_dataset, [representative_index]
        )[0]
        del model
        if torch.cuda.is_available():
            torch.cuda.empty_cache()
    plot_baseline_field_comparison(model_records, output_dir)
    pcfno_seed = chosen_seeds.get("PCFNO", seeds[0])
    pcfno_sigma_scale = load_sigma_scale(output_dir, "PCFNO", pcfno_seed)
    pcfno = load_trained_model("PCFNO", pcfno_seed, output_dir)
    three_records: List[Dict[str, Any]] = []
    for case in sorted(representative_cases, key=lambda item: item.reynolds):
        indices = dataset_indices_for_case(test_dataset, case.name)
        selected_index = indices[int(round(0.75 * (len(indices) - 1)))]
        three_records.extend(
            predict_indices(
                pcfno,
                test_dataset,
                [selected_index],
                sigma_scale=pcfno_sigma_scale,
            )
        )
    plot_three_reynolds_cases(three_records, output_dir)
    save_dat(
        output_dir / "tables" / "table08_three_reynolds_metrics.dat",
        aggregate_three_re_metrics(three_records, "PCFNO", pcfno_seed),
    )
    evolution_indices = [
        middle_indices[index]
        for index in np.linspace(0, len(middle_indices) - 1, 3, dtype=int)
    ]
    evolution_records = predict_indices(
        pcfno,
        test_dataset,
        evolution_indices,
        sigma_scale=pcfno_sigma_scale,
    )
    plot_time_evolution(evolution_records, output_dir)
    plot_error_vs_time_and_re(per_sample_by_model_seed, output_dir)
    plot_uncertainty(output_dir, seeds)
    spectrum_dataset = CFDBenchSpeedDataset(
        [middle_case],
        test_dataset.time_scale,
        test_dataset.resolution,
        int(CONFIG["spectrum_time_samples"]),
        test_dataset.re_min,
        test_dataset.re_max,
    )
    spectrum_records = predict_indices(
        pcfno,
        spectrum_dataset,
        list(range(len(spectrum_dataset))),
        sigma_scale=pcfno_sigma_scale,
    )
    plot_probe_and_spectrum(spectrum_records, output_dir)
    plot_speed_profiles(model_records["PCFNO"], output_dir)
    del pcfno
    if torch.cuda.is_available():
        torch.cuda.empty_cache()
    descriptions = {
        "fig01_training_convergence": "Seed-averaged convergence of baselines and ablations",
        "fig02_pcfno_loss_decomposition": "PCFNO weighted loss contributions",
        "fig03_baseline_metrics": "Global, wake-region, and near-wall accuracy metrics",
        "fig04_ablation_metrics": "PCFNO ablation metrics",
        "fig05_runtime_comparison": "Training and inference times",
        "fig06_baseline_field_comparison": "Speed-magnitude fields and errors",
        "fig07_three_reynolds_cases": "PCFNO results for three fully held-out Reynolds-number groups",
        "fig08_time_evolution": "Speed-magnitude evolution for a representative case",
        "fig09_error_vs_time_and_re": "FNO and PCFNO error versus time and Reynolds number",
        "fig10_uncertainty_calibration": "Raw and validation-calibrated PCFNO uncertainty",
        "fig11_speed_probe_spectrum": "Speed-magnitude probe signal and spectrum",
        "fig12_speed_profiles": "Downstream speed-magnitude profiles",
    }
    manifest: List[Dict[str, Any]] = []
    for stem, description in descriptions.items():
        image = output_dir / "figures" / f"{stem}.jpg"
        data = output_dir / "figure_data" / f"{stem}.dat"
        if image.exists():
            manifest.append(
                {
                    "figure": image.name,
                    "data_file": data.name,
                    "data_exists": data.exists(),
                    "description": description,
                }
            )
    save_dat(output_dir / "figure_data" / "figure_manifest.dat", manifest)


def save_dataset_and_config_tables(
    train_cases: Sequence[CaseInfo],
    val_cases: Sequence[CaseInfo],
    test_cases: Sequence[CaseInfo],
    excluded_cases: Sequence[CaseInfo],
    representative_cases: Sequence[CaseInfo],
    train_dataset: CFDBenchSpeedDataset,
    val_dataset: CFDBenchSpeedDataset,
    test_dataset: CFDBenchSpeedDataset,
    time_scale: float,
    output_dir: Path,
    chosen_font: str,
) -> None:
    split_rows: List[Dict[str, Any]] = []
    split_re_sets: Dict[str, set[float]] = {}
    for split, cases, dataset in [
        ("train", train_cases, train_dataset),
        ("validation", val_cases, val_dataset),
        ("test", test_cases, test_dataset),
    ]:
        re_values = sorted({canonical_re(case.reynolds) for case in cases})
        split_re_sets[split] = set(re_values)
        split_rows.append(
            {
                "split": split,
                "split_profile": CONFIG.get("split_profile", "confirmation"),
                "cases": len(cases),
                "unique_Re_groups": len(re_values),
                "selected_frames": len(dataset),
                "Re_min": min(re_values),
                "Re_max": max(re_values),
                "Re_values": ",".join(f"{value:g}" for value in re_values),
            }
        )
    excluded_values = sorted(
        float(value) for value in CONFIG.get("resolved_excluded_reynolds", [])
    )
    if excluded_values:
        split_rows.append(
            {
                "split": "excluded",
                "split_profile": CONFIG.get("split_profile", "confirmation"),
                "cases": len(excluded_cases),
                "unique_Re_groups": len(excluded_values),
                "selected_frames": 0,
                "Re_min": min(excluded_values),
                "Re_max": max(excluded_values),
                "Re_values": ",".join(f"{value:g}" for value in excluded_values),
            }
        )

    overlap_rows = [
        {
            "comparison": "train_validation",
            "overlap": sorted(split_re_sets["train"] & split_re_sets["validation"]),
        },
        {
            "comparison": "train_test",
            "overlap": sorted(split_re_sets["train"] & split_re_sets["test"]),
        },
        {
            "comparison": "validation_test",
            "overlap": sorted(split_re_sets["validation"] & split_re_sets["test"]),
        },
        {
            "comparison": "active_excluded",
            "overlap": sorted(
                (
                    split_re_sets["train"]
                    | split_re_sets["validation"]
                    | split_re_sets["test"]
                )
                & set(excluded_values)
            ),
        },
    ]
    if any(row["overlap"] for row in overlap_rows):
        raise RuntimeError(f"Reynolds leakage detected: {overlap_rows}")
    for row in overlap_rows:
        row["overlap"] = "none"
    save_dat(output_dir / "tables" / "table01_dataset_split.dat", split_rows)
    save_dat(output_dir / "data" / "reynolds_split_audit.dat", overlap_rows)

    hyper_rows = [{"parameter": key, "value": value} for key, value in CONFIG.items()]
    hyper_rows.extend(
        [
            {"parameter": "time_star_scale", "value": time_scale},
            {"parameter": "actual_plot_font", "value": chosen_font},
            {
                "parameter": "target_definition",
                "value": "first public channel u.npy interpreted as speed magnitude and normalized by Uin",
            },
            {
                "parameter": "split_definition",
                "value": (
                    "complete Reynolds-number groups are disjoint across train, "
                    "validation, and test; confirmation evaluates unseen "
                    "Reynolds-number groups within the covered parameter range"
                ),
            },
            {
                "parameter": "boundary_constraint",
                "value": "zero-speed interpolation at the SDF zero level across solid-fluid grid edges",
            },
            {
                "parameter": "pcfno_geometry_fusion",
                "value": (
                    "local physics-guided residual adapter using frozen FNO features, "
                    "physical input fields, and encoded SDF features"
                ),
            },
            {
                "parameter": "statistical_unit",
                "value": "seed-averaged physical frames, cases, and Reynolds groups with case-cluster bootstrap confidence intervals",
            },
        ]
    )
    save_dat(output_dir / "tables" / "table02_hyperparameters.dat", hyper_rows)
    representative_names = {case.name for case in representative_cases}
    split_map = {case.name: "train" for case in train_cases}
    split_map.update({case.name: "validation" for case in val_cases})
    split_map.update({case.name: "test" for case in test_cases})
    split_map.update({case.name: "excluded" for case in excluded_cases})
    case_rows: List[Dict[str, Any]] = []
    for case in (
        list(train_cases) + list(val_cases) + list(test_cases) + list(excluded_cases)
    ):
        case_rows.append(
            {
                "case": case.name,
                "split": split_map[case.name],
                "representative_Re_case": case.name in representative_names,
                "frames": case.num_frames,
                "Re": case.reynolds,
                "Re_group": canonical_re(case.reynolds),
                "Uin": case.inlet_velocity,
                "diameter": case.diameter,
                "density": case.params["density"],
                "viscosity": case.params["viscosity"],
            }
        )
    save_dat(output_dir / "data" / "case_split_and_parameters.dat", case_rows)
    save_dat(
        output_dir / "tables" / "table03_representative_reynolds_cases.dat",
        [
            {
                "case": case.name,
                "Re": case.reynolds,
                "target_Re": canonical_re(case.reynolds),
                "split": "test",
            }
            for case in representative_cases
        ],
    )
    save_json(output_dir / "config.json", CONFIG)
    system_rows = [
        {"item": "python", "value": sys.version.replace("\n", " ")},
        {"item": "platform", "value": platform.platform()},
        {"item": "torch", "value": torch.__version__},
        {"item": "device", "value": CONFIG["device"]},
        {"item": "cuda_available", "value": torch.cuda.is_available()},
        {
            "item": "cuda_device",
            "value": (
                torch.cuda.get_device_name(0) if torch.cuda.is_available() else "NA"
            ),
        },
    ]
    save_dat(output_dir / "data" / "system_information.dat", system_rows)


def audit_first_channel(cases: Sequence[CaseInfo], output_dir: Path) -> None:
    rows: List[Dict[str, Any]] = []
    selected = nearest_unique_cases(cases, CONFIG["representative_re_targets"])
    for case in selected:
        array = np.load(case.directory / "u.npy", mmap_mode="r")
        indices = select_frame_indices(case.num_frames, 3)
        for index in indices:
            frame = np.asarray(
                array if array.ndim == 2 else array[int(index)], dtype=np.float64
            )
            rows.append(
                {
                    "case": case.name,
                    "Re": case.reynolds,
                    "frame_index": int(index),
                    "minimum": float(np.min(frame)),
                    "maximum": float(np.max(frame)),
                    "mean": float(np.mean(frame)),
                    "negative_fraction": float(np.mean(frame < -1e-8)),
                }
            )
    save_dat(output_dir / "data" / "first_channel_speed_audit.dat", rows)


def main() -> None:
    args = parse_args()
    apply_cli_overrides(args)
    chosen_font = configure_matplotlib()
    output_dir = Path(CONFIG["output_dir"]).resolve()
    if output_dir.exists() and args.smoke_test:
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    cleanup_figure_formats(output_dir)
    logger = setup_logging(output_dir)
    setup_native_diagnostics(output_dir, logger)
    logger.info("Configuration: %s", json.dumps(CONFIG, indent=2))
    logger.info("Plot font selected: %s", chosen_font)
    logger.info(
        "Grouped-Re split profile: %s | requested test=%s | excluded=%s",
        CONFIG.get("split_profile"),
        CONFIG.get("test_reynolds"),
        CONFIG.get("excluded_reynolds"),
    )
    if CONFIG.get("split_profile") == "development":
        logger.warning(
            "The development split is intended for reproducibility checks; "
            "use the confirmation profile for the reported experiments."
        )
    data_root = Path(CONFIG["data_path"]).resolve()
    unzip_cfdbench(data_root)
    if CONFIG["dataset_type"] != "prop":
        raise ValueError(
            "This script is intentionally restricted to CFDBench cylinder PROP"
        )
    cases = discover_cases(data_root, CONFIG["dataset_type"])
    native_shapes = sorted({(case.height, case.width) for case in cases})
    if any(
        shape != (int(CONFIG["resolution"]), int(CONFIG["resolution"]))
        for shape in native_shapes
    ):
        logger.warning(
            "Requested resolution=%d differs from native data shapes=%s. Values are interpolated; 64 is recommended for the primary PoF experiment.",
            int(CONFIG["resolution"]),
            native_shapes,
        )
    logger.info(
        "Physics curriculum: cosine ramp over %d epochs (%.1f%% of %d PCFNO epochs)",
        int(CONFIG["physics_warmup_epochs"]),
        100.0
        * float(CONFIG["physics_warmup_epochs"])
        / max(int(CONFIG["pcfno_max_epochs"]), 1),
        int(CONFIG["pcfno_max_epochs"]),
    )
    audit_first_channel(cases, output_dir)
    train_cases, val_cases, test_cases, representative_cases = split_cases(cases)
    excluded_re_set = set(CONFIG.get("resolved_excluded_reynolds", []))
    excluded_cases = [
        case for case in cases if canonical_re(case.reynolds) in excluded_re_set
    ]
    re_min = min(case.reynolds for case in train_cases)
    re_max = max(case.reynolds for case in train_cases)
    time_scale = compute_time_scale(train_cases)
    train_dataset = CFDBenchSpeedDataset(
        train_cases,
        time_scale,
        CONFIG["resolution"],
        CONFIG["time_samples_per_case"],
        re_min,
        re_max,
    )
    val_dataset = CFDBenchSpeedDataset(
        val_cases,
        time_scale,
        CONFIG["resolution"],
        CONFIG["time_samples_per_case"],
        re_min,
        re_max,
    )
    test_dataset = CFDBenchSpeedDataset(
        test_cases,
        time_scale,
        CONFIG["resolution"],
        CONFIG["time_samples_per_case"],
        re_min,
        re_max,
    )
    save_dataset_and_config_tables(
        train_cases,
        val_cases,
        test_cases,
        excluded_cases,
        representative_cases,
        train_dataset,
        val_dataset,
        test_dataset,
        time_scale,
        output_dir,
        chosen_font,
    )
    logger.info(
        "Cases: train=%d, validation=%d, test=%d",
        len(train_cases),
        len(val_cases),
        len(test_cases),
    )
    logger.info(
        "Representative cases from fully held-out Re groups: %s",
        [(case.name, case.reynolds) for case in representative_cases],
    )
    selected_models = [
        name for name in CONFIG["selected_models"] if name in ALL_EXPERIMENTS
    ]
    if (
        any(name.startswith("PCFNO") for name in selected_models)
        and "FNO" not in selected_models
    ):
        selected_models = ["FNO"] + selected_models
        logger.info(
            "Automatically prepended FNO because PCFNO uses the matching FNO checkpoint."
        )
    seeds = [int(seed) for seed in CONFIG["training_seeds"]]
    histories: Dict[Tuple[str, int], List[Dict[str, Any]]] = {}
    training_runtime_rows: List[Dict[str, Any]] = []
    for name in selected_models:
        for seed in seeds:
            history, runtime = train_model(
                name, seed, train_dataset, val_dataset, output_dir, logger
            )
            histories[(name, seed)] = history
            training_runtime_rows.append(runtime)
    save_dat(
        output_dir / "data" / "training_runtime_by_seed.dat", training_runtime_rows
    )
    calibration_rows: List[Dict[str, Any]] = []
    for name in selected_models:
        if not name.startswith("PCFNO"):
            continue
        for seed in seeds:
            calibration_rows.append(
                fit_uncertainty_scale(name, seed, val_dataset, output_dir, logger)
            )
    save_dat(
        output_dir / "tables" / "table07_sigma_calibration_by_seed.dat",
        calibration_rows,
    )
    seed_summaries: List[Dict[str, Any]] = []
    per_sample_by_model_seed: Dict[Tuple[str, int], List[Dict[str, Any]]] = {}
    for name in selected_models:
        for seed in seeds:
            summary, per_sample = evaluate_model(
                name, seed, test_dataset, output_dir, logger
            )
            training_row = next(
                row
                for row in training_runtime_rows
                if row["model"] == name and int(row["seed"]) == seed
            )
            summary["training_time_seconds"] = training_row["training_time_seconds"]
            summary["fno_pretraining_time_seconds"] = training_row.get(
                "fno_pretraining_time_seconds", 0.0
            )
            summary["effective_training_time_seconds"] = training_row.get(
                "effective_training_time_seconds",
                training_row["training_time_seconds"],
            )
            summary["epochs_completed"] = training_row["epochs_completed"]
            summary["best_epoch"] = training_row["best_epoch"]
            summary["training_time_per_epoch_seconds"] = float(
                training_row["training_time_seconds"]
            ) / max(int(training_row["epochs_completed"]), 1)
            seed_summaries.append(summary)
            per_sample_by_model_seed[(name, seed)] = per_sample
    aggregate_summaries = aggregate_model_summaries(seed_summaries)
    save_dat(output_dir / "data" / "all_model_metrics_by_seed.dat", seed_summaries)
    save_dat(
        output_dir / "data" / "all_model_metrics_aggregate.dat", aggregate_summaries
    )
    save_dat(
        output_dir / "tables" / "table04_baseline_metrics.dat",
        rows_for_names(aggregate_summaries, BASELINE_NAMES),
    )
    save_dat(
        output_dir / "tables" / "table05_ablation_metrics.dat",
        rows_for_names(aggregate_summaries, ["PCFNO"] + ABLATION_NAMES),
    )
    runtime_fields = [
        "model",
        "n_seeds",
        "training_time_seconds_mean",
        "training_time_seconds_std",
        "fno_pretraining_time_seconds_mean",
        "effective_training_time_seconds_mean",
        "effective_training_time_seconds_std",
        "epochs_completed_mean",
        "best_epoch_mean",
        "best_val_objective_mean",
        "training_time_per_epoch_seconds_mean",
        "training_time_per_epoch_seconds_std",
        "Inference_ms_per_frame_mean",
        "Inference_ms_per_frame_std",
        "Throughput_frames_per_second_mean",
        "Reference_speedup_vs_CFD_mean",
        "Parameters_mean",
    ]
    save_dat(
        output_dir / "tables" / "table06_model_runtime.dat",
        aggregate_summaries,
        runtime_fields,
    )
    uncertainty_rows = [
        row for row in aggregate_summaries if "Mean_sigma_over_Uin_mean" in row
    ]
    save_dat(
        output_dir / "tables" / "table07_uncertainty_metrics.dat", uncertainty_rows
    )
    diagnostic_rows: List[Dict[str, Any]] = []
    aggregate_map = {row["model"]: row for row in aggregate_summaries}
    if "FNO" in aggregate_map and "PCFNO" in aggregate_map:
        fno_error = float(aggregate_map["FNO"]["Relative_L2_mean"])
        pcfno_error = float(aggregate_map["PCFNO"]["Relative_L2_mean"])
        diagnostic_rows.append(
            {
                "comparison": "PCFNO_vs_FNO",
                "FNO_Relative_L2": fno_error,
                "PCFNO_Relative_L2": pcfno_error,
                "PCFNO_relative_improvement_percent": 100.0
                * (fno_error - pcfno_error)
                / max(fno_error, 1e-12),
                "PCFNO_outperforms_FNO": pcfno_error < fno_error,
                "FNO_Relative_L2_median": float(
                    aggregate_map["FNO"].get("Relative_L2_median", float("nan"))
                ),
                "PCFNO_Relative_L2_median": float(
                    aggregate_map["PCFNO"].get("Relative_L2_median", float("nan"))
                ),
            }
        )
        if pcfno_error >= fno_error:
            logger.warning(
                "PCFNO does not outperform FNO after multi-seed evaluation; do not claim an accuracy advantage without further tuning."
            )
    save_dat(output_dir / "data" / "pcfno_vs_fno_diagnostic.dat", diagnostic_rows)
    paired_rows = paired_fno_pcfno_statistics(per_sample_by_model_seed)
    save_dat(
        output_dir / "tables" / "table10_paired_fno_pcfno_statistics.dat",
        paired_rows,
    )
    save_dat(
        output_dir / "tables" / "table11_case_re_statistical_summary.dat",
        [
            row
            for row in paired_rows
            if row.get("scope") == "seed_averaged"
            and row.get("analysis_unit") in {"case", "Re_group"}
        ],
    )
    if (
        not args.skip_plots
        and all(name in selected_models for name in BASELINE_NAMES)
        and "PCFNO" in selected_models
    ):
        generate_publication_outputs(
            test_dataset,
            representative_cases,
            histories,
            aggregate_summaries,
            per_sample_by_model_seed,
            output_dir,
            logger,
        )
    logger.info("All experiments completed. Results saved to %s", output_dir)


if __name__ == "__main__":
    main()

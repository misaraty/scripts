## [中文版本](https://www.misaraty.com/2026-10-06_mcrystalframer/)

## mCrystalFramer

`mCrystalFramer` is a standalone pure-PyTorch CrystalFramer-style implementation for crystal band-gap prediction from CIF structures.

The code retains the main ideas used by CrystalFramer, including periodic-image attention, distance-dependent attention bias, Gaussian distance encoding, attention-head-specific dynamic max frames, frame-based angular encoding, Transformer-style residual blocks, crystal-level mean pooling, and a scalar regression head. It uses `pymatgen` for CIF parsing and does not depend on PyTorch Geometric, DGL, JARVIS, CuPy, custom CUDA kernels, or the original CrystalFramer repository.

This standalone version uses finite periodic-image summation implemented with native PyTorch tensor operations. CUDA is still used automatically when available, but its numerical behavior and runtime can differ from those of the original fused CuPy/CUDA implementation.

## Requirements

```bash
pip install torch pymatgen numpy pandas openpyxl scikit-learn matplotlib tqdm
```

Install Optuna only when `USE_OPTUNA = True`:

```bash
pip install optuna
```

## Data

Prepare the following files in the script directory:

```text
mCrystalFramer_v2.py
data.xlsx
cif/
|-- 1.cif
|-- 2.cif
|-- 3.cif
```

The first two columns of `data.xlsx` are used:

| cif | bandgap |
| ---: | ---: |
| 1 | 1.23 |
| 2 | 0.87 |
| 3 | 2.15 |

The value `1` in the first column corresponds to `./cif/1.cif`. Integer values, integer-valued floating-point values such as `1.0`, names without the `.cif` suffix, and complete CIF filenames are accepted.

## Usage

```bash
python mCrystalFramer_v2.py
```

The script automatically validates the CIF files, performs a fixed 80/10/10 train/validation/test split with `SEED = 42`, standardizes the target using training-set statistics, trains with `MSELoss`, applies validation-RMSE-based learning-rate scheduling and early stopping, reloads the best checkpoint, evaluates all three splits, and generates the figures and data files. CUDA and BF16 automatic mixed precision are used when available; otherwise, the script runs on CPU.

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. With the default settings, the output directory is:

```text
CrystalFramer_v2/
|-- CrystalFramer_best.pt
|-- figure/
|-- dat/
|-- table/
|-- log/
|-- split/
|-- cache/
```

The outputs include train/validation/test MAE, RMSE, and R2; individual and combined parity plots with their data; the train/validation RMSE curve and its data; the reusable fixed data split; the best checkpoint; and the complete training log. The checkpoint stores the model configuration, target normalization statistics, best epoch, best validation RMSE, and graph configuration.

The script also provides `load_trained_model()` and `predict_cifs()` for loading the saved checkpoint and predicting additional CIF structures from another Python program.

## Memory note

Periodic attention and frame-based angular encoding can require substantial GPU memory for structures with many atoms. If CUDA memory is insufficient, reduce the parameters in the following order:

```python
BATCH_SIZE = 1
ANGLE_BASIS_DIM = 24
NUM_HEADS = 4
DISTANCE_BASIS_DIM = 32
```

Keep `LATTICE_RANGE = 1` whenever possible because setting it to zero removes noncentral periodic images. `ANGLE_BASIS_DIM` must be divisible by 3, and `MODEL_DIM` must be divisible by `NUM_HEADS`.

## Citation

Original CrystalFramer reference:

* [Ito Y, Taniai T, Igarashi R, et al. Rethinking the role of frames for SE (3)-invariant crystal structure modeling. International Conference on Learning Representations. 2025, 2025: 8655-8676.](https://proceedings.iclr.cc/paper_files/paper/2025/hash/187d94b3c93343f0e925b5cf729eadd5-Abstract-Conference.html)

This work:

To be added after the paper is officially published.

## [中文版本](https://www.misaraty.com/2026-10-06_mcrystalformer/)

## mCrystalformer

`mCrystalformer` is a standalone pure-PyTorch implementation of Crystalformer for crystal band-gap prediction from CIF structures.

The code retains the central Crystalformer design: fully connected periodic attention, real-space and reciprocal-space periodic encoding, query-dependent Gaussian attention, radial value encoding, alternating Latticeformer encoder blocks, T-Fixup initialization, crystal-level pooling, and a scalar regression head. `pymatgen` is used for CIF parsing and primitive-cell standardization.

This version does not depend on PyTorch Geometric, DGL, JARVIS, CuPy, `pytorch-pfn-extras`, custom CUDA source files, JSON parameter files, or the original Crystalformer repository. When CUDA is available, the standard PyTorch tensor operations still run on the GPU. The pure-PyTorch implementation may be slower than the fused CUDA-kernel implementation in the original repository.

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
mCrystalformer_v2.py
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

The first column must contain integer-compatible identifiers. The value `1` is converted to `1.cif` and matched to `./cif/1.cif`. The second column contains the band gap in eV.

Invalid targets, missing CIF files, unparsable structures, invalid occupancies, and singular lattices are skipped and recorded in the training log. At least 20 valid structures are required.

## Usage

```bash
python mCrystalformer_v2.py
```

The script automatically performs CIF validation and caching, primitive-cell standardization, a fixed 80/10/10 train/validation/test split with `SEED = 42`, target normalization fitted only on the training set, model training, validation-RMSE early stopping, checkpoint reloading, final evaluation, and plotting. CUDA is used when available; otherwise, the script runs on CPU.

The default model uses four 128-dimensional Latticeformer blocks with eight attention heads. `DOMAIN = "real-reci"` alternates real-space and reciprocal-space periodic attention between blocks. Training uses MSE loss, AdamW, inverse-square-root learning-rate decay, and a maximum of 300 epochs. Optional Optuna optimization is controlled by `USE_OPTUNA`.

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. With the default settings, the output directory is:

```text
Crystalformer_v2/
|-- Crystalformer_best.pt
|-- figure/
|   |-- Crystalformer_parity_train.jpg
|   |-- Crystalformer_parity_val.jpg
|   |-- Crystalformer_parity_test.jpg
|   |-- Crystalformer_parity_all.jpg
|   `-- Crystalformer_rmse_curve.jpg
|-- dat/
|-- table/
|-- log/
|-- split/
`-- cache/
```

The outputs include train/validation/test MAE, RMSE, and R2; individual and combined parity plots with tab-separated prediction data; the RMSE learning curve and data; the reusable fixed split; the parsed-structure cache; the best checkpoint selected by validation RMSE; and the complete training log. Set a new `RUN_VERSION` to create a separate result directory without overwriting earlier runs.

## Citation

Original Crystalformer reference:

* [Taniai T, Igarashi R, Suzuki Y, et al. Crystalformer: Infinitely connected attention for periodic structure encoding. International Conference on Learning Representations. 2024, 2024: 45083-45105](https://proceedings.iclr.cc/paper_files/paper/2024/hash/c428adf74782c2092d254329b6b02482-Abstract-Conference.html)

This work:

To be added after the paper is officially published.

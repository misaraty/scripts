## [中文版本](https://www.misaraty.com/2026-10-06_mm3gnet/)

## mM3GNet

`mM3GNet` is a standalone PyTorch/PyTorch Geometric implementation of M3GNet for crystal band-gap prediction from CIF structures.

The code retains the main M3GNet components, including the periodic radius graph, atomic embeddings, spherical Bessel radial basis, angular three-body basis, explicit three-body interactions, gated residual bond and atom updates, weighted atom readout, and regression head. It uses `pymatgen` for CIF parsing and periodic graph construction and does not depend on TensorFlow or the original M3GNet repository.

## Requirements

```bash
pip install torch torch-geometric pymatgen numpy pandas openpyxl scikit-learn matplotlib tqdm
```

Install Optuna only when `USE_OPTUNA = True`:

```bash
pip install optuna
```

## Data

Prepare the following files in the script directory:

```text
mM3GNet_v1.py
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

The value `1` in the first column corresponds to `./cif/1.cif`. Numeric values such as `1.0` are also converted to `1.cif`.

## Usage

```bash
python mM3GNet_v1.py
```

The script automatically performs CIF validation, periodic graph and triplet construction, a fixed 80/10/10 train/validation/test split, model training, early stopping, evaluation, and plotting. CUDA is used when available; otherwise, the script runs on CPU. BF16 is enabled automatically on supported CUDA devices.

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. With the default settings, the output directory is:

```text
M3GNet_v1/
|-- M3GNet_best.pt
|-- figure/
|-- dat/
|-- table/
|-- log/
|-- split/
|-- cache/
```

The outputs include train/validation/test MAE, RMSE, and R2; parity plots and data; the RMSE learning curve and data; the fixed data split; cached crystal graphs; the best checkpoint; and the complete training log.

## Citation

Original M3GNet reference:

* [Chen C, Ong S P. A universal graph deep learning interatomic potential for the periodic table. Nature Computational Science, 2022, 2(11): 718-728.](https://www.nature.com/articles/s43588-022-00349-3)

This work:

To be added after the paper is officially published.

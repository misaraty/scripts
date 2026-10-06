## [中文版本](https://www.misaraty.com/2026-10-06_mpotnet/)

## mPotNet

`mPotNet` is a standalone PyTorch/PyTorch Geometric implementation of PotNet for crystal band-gap prediction from CIF structures.

The code retains the central PotNet design: a local periodic neighbor graph is combined with a complete atom-pair graph carrying periodic Coulomb, London-dispersion, and Pauli-potential features. The potential features and local distances are expanded with radial basis functions and processed by gated message-passing layers, global mean pooling, and a scalar regression head.

The script uses `pymatgen` for CIF parsing and periodic graph construction. The periodic summation routines are implemented directly with NumPy and SciPy, so the code does not depend on JARVIS, the original PotNet repository, repository-internal modules, or separately compiled Cython/GSL extensions.

## Requirements

```bash
pip install torch torch-geometric pymatgen scipy numpy pandas openpyxl scikit-learn matplotlib tqdm
```

Install Optuna only when `USE_OPTUNA = True`:

```bash
pip install optuna
```

## Data

Prepare the following files in the script directory:

```text
mPotNet_v1.py
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

The value `1` in the first column corresponds to `./cif/1.cif`. A `.cif` suffix is added automatically when it is absent. Invalid, missing, duplicated, non-finite, or disordered samples are skipped and recorded in the log.

## Usage

```bash
python mPotNet_v1.py
```

The script automatically performs CIF validation, periodic-potential preprocessing and caching, a fixed 80/10/10 train/validation/test split, target normalization fitted only on the training set, model training, early stopping, evaluation, and plotting. CUDA is used when available; otherwise, the script runs on CPU. BF16 automatic mixed precision is enabled by default on CUDA.

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. With the default settings, the output directory is:

```text
PotNet_v1/
|-- PotNet_best.pt
|-- figure/
|-- dat/
|-- table/
|-- log/
|-- split/
|-- cache/
```

The outputs include train/validation/test MAE, RMSE, and R2; individual and combined parity plots and data; the RMSE learning curve and data; the reusable fixed split; cached PotNet graphs; the best checkpoint selected by validation RMSE; and the complete training log.

The functions `load_trained_model()` and `predict_cifs()` can be imported by another Python program to load the saved checkpoint and predict additional CIF structures.

## Memory Note

PotNet constructs a complete long-range graph with approximately `N^2` atom-pair edges per structure. GPU memory therefore grows rapidly with the number of atoms. For large unit cells, reduce `BATCH_SIZE` first, followed by `HIDDEN_DIM` and `POTENTIAL_RBF_BINS`. Keep `USE_AMP = True` and `AMP_DTYPE = "bfloat16"` on supported NVIDIA GPUs.

## Citation

Original PotNet reference:

* [Lin Y, Yan K, Luo Y, et al. Efficient approximations of complete interatomic potentials for crystal property prediction. International conference on machine learning. PMLR, 2023: 21260-21287](https://proceedings.mlr.press/v202/lin23m.html)

This work:

To be added after the paper is officially published.

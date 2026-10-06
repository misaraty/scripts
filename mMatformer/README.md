## [中文版本](https://www.misaraty.com/2026-10-06_mmatformer/)

## mMatformer

`mMatformer` is a standalone PyTorch/PyTorch Geometric implementation of Matformer for crystal band-gap prediction from CIF structures.

The code retains the main Matformer components, including the periodic crystal multigraph, embedded 92-dimensional element descriptors, radial-basis distance encoding, multi-head gated attention, residual node updates, graph-level mean pooling, and regression head. It uses `pymatgen` for CIF parsing and periodic graph construction and does not depend on DGL, JARVIS, `torch_scatter`, or the original Matformer repository.

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
mMatformer_v1.py
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

The value `1` in the first column corresponds to `./cif/1.cif`.

## Usage

```bash
python mMatformer_v1.py
```

The script automatically performs CIF validation, periodic graph construction and caching, a fixed 80/10/10 train/validation/test split, target standardization, model training, early stopping, evaluation, and plotting. CUDA is used when available; otherwise, the script runs on CPU.

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. With the default settings, the output directory is:

```text
Matformer_v1/
|-- Matformer_best.pt
|-- figure/
|-- dat/
|-- table/
|-- log/
|-- split/
|-- cache/
```

The outputs include train/validation/test MAE, RMSE, and R2; individual and combined parity plots with their data; the RMSE learning curve and its data; the fixed data split; the best checkpoint selected by validation RMSE; and the complete training log.

## Citation

Original Matformer reference:

* [Yan K, Liu Y, Lin Y, et al. Periodic graph transformers for crystal material property prediction. Advances in Neural Information Processing Systems, 2022, 35: 15066-15080](https://proceedings.neurips.cc/paper_files/paper/2022/hash/6145c70a4a4bf353a31ac5496a72a72d-Abstract-Conference.html)

This work:

To be added after the paper is officially published.
## [中文版本](待补充)

## mALIGNN

`mALIGNN` is a standalone PyTorch/PyTorch Geometric implementation of ALIGNN for crystal band-gap prediction from CIF structures.

The code retains the main ALIGNN components, including the periodic crystal graph, bond-angle line graph, distance and angle encoding, edge-gated message passing, graph pooling, and regression head. It uses `pymatgen` for CIF parsing and graph construction and does not depend on DGL, JARVIS, or the original ALIGNN repository.

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
mALIGNN_v2.py
data.xlsx
cif/
|-- 1.cif
|-- 2.cif
`-- 3.cif
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
python mALIGNN_v2.py
```

The script automatically performs CIF validation, a fixed 80/10/10 train/validation/test split, model training, early stopping, evaluation, and plotting. CUDA is used when available; otherwise, the script runs on CPU.

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. With the default settings, the output directory is:

```text
ALIGNN_v2/
|-- ALIGNN_best.pt
|-- figure/
|-- dat/
|-- table/
|-- log/
|-- split/
`-- cache/
```

The outputs include train/validation/test MAE, RMSE, and R2; parity plots and data; the RMSE learning curve and data; the fixed data split; the best checkpoint; and the complete training log.

## Citation

Original ALIGNN reference:

* [Choudhary K, DeCost B. Atomistic line graph neural network for improved materials property predictions. npj Computational Materials, 2021, 7(1): 185](https://www.nature.com/articles/s41524-021-00650-1)

This work:

To be added after the paper is officially published.

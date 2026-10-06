## [中文版本](https://www.misaraty.com/2026-10-06_mdensegnn/)

## mDenseGNN

`mDenseGNN` is a standalone PyTorch/PyTorch Geometric implementation of a DenseGNN-style crystal graph neural network for band-gap regression from CIF structures.

The code retains the principal DenseGNN mechanisms used by this implementation: dense connections between message-passing layers, hierarchical node-edge-graph residual updates, learned atomic embeddings augmented with periodic-table properties, periodic Voronoi crystal graphs, distance and Voronoi ridge-area edge features, graph-level pooling, and a scalar regression head. It uses `pymatgen` and SciPy for CIF parsing and periodic graph construction and does not depend on the original DenseGNN repository or its internal modules.

This standalone version is restricted to single-target CIF-to-band-gap regression. Force-field training, classification, multitask learning, Matbench downloading, cross-validation, and ablation workflows are not included.

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
mDenseGNN_v1.py
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

The first column must contain numeric CIF identifiers. The script converts each value with `int(raw_name)` and appends `.cif`; therefore, the value `1` corresponds to `./cif/1.cif`. The second column contains the band gap in eV. Empty entries, nonfinite targets, missing CIF files, duplicate identifiers, invalid CIF files, and disordered or partially occupied structures are skipped and recorded in the log.

The input locations and target metadata can be changed near the beginning of the script:

```python
EXCEL_PATH = "./data.xlsx"
CIF_DIR = "./cif"
TARGET_NAME = "bandgap"
TARGET_UNIT = "eV"
```

## Usage

```bash
python mDenseGNN_v1.py
```

The script automatically performs CIF validation, periodic Voronoi graph construction and caching, a fixed 80/10/10 train/validation/test split with `seed=42`, target standardization, model training, early stopping, evaluation, and plotting. CUDA is used when available; otherwise, the script runs on CPU.

Training uses `nn.MSELoss()` and AdamW. Validation RMSE is monitored by `ReduceLROnPlateau`, and the checkpoint with the lowest validation RMSE is saved. With the default configuration, the model uses five DenseGNN layers with 128-dimensional hidden features, a Voronoi ridge-area threshold of 0.10, batch size 32, a maximum of 300 epochs, and early-stopping patience of 50 epochs.

Set `USE_OPTUNA = True` to enable the optional Optuna search over hidden dimension, depth, batch size, learning rate, and weight decay. The fixed train/validation/test split is reused during optimization.

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. With the default settings, the output directory is:

```text
DenseGNN_v1/
|-- DenseGNN_best.pt
|-- figure/
|   |-- DenseGNN_parity_train.jpg
|   |-- DenseGNN_parity_val.jpg
|   |-- DenseGNN_parity_test.jpg
|   |-- DenseGNN_parity_all.jpg
|   `-- DenseGNN_rmse_curve.jpg
|-- dat/
|-- table/
|   `-- DenseGNN_metrics.dat
|-- log/
|   `-- DenseGNN_training.log
|-- split/
|   `-- DenseGNN_split.csv
`-- cache/
    `-- DenseGNN_graphs/
```

The outputs include train/validation/test MAE, RMSE, and R2; separate and combined parity plots with their `.dat` files; the training/validation RMSE curve and its `.dat` file; the fixed data split; cached crystal graphs; the best checkpoint; and the complete training log. Figures are saved as 600 dpi JPG files.

The saved checkpoint can also be loaded programmatically with `load_trained_model()`, and `predict_cifs()` can be used to predict band gaps for new CIF files.

## Citation

Original DenseGNN reference:

* [Du H, Wang J, Hui J, et al. DenseGNN: universal and scalable deeper graph neural networks for high-performance property prediction in crystals and molecules. npj Computational Materials, 2024, 10(1): 292.](https://www.nature.com/articles/s41524-024-01444-x)

Original DenseGNN repository:

* [https://github.com/dhw059/DenseGNN](https://github.com/dhw059/DenseGNN)

This work:

To be added after the paper is officially published.

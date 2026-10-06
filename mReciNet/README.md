## [中文版本](https://www.misaraty.com/2026-10-06_mrecinet/)

## mReciNet

`mReciNet` is a standalone PyTorch/PyTorch Geometric implementation of ReciNet for crystal band-gap prediction from CIF structures.

The code combines a periodic local crystal graph with reciprocal-space long-range modeling. It retains the main ReciNet components, including gated local message passing, reciprocal-vector generation from each crystal lattice, reciprocal structure-factor updates, layerwise fusion of local and long-range representations, graph-level mean pooling, and a scalar regression head. CIF parsing and periodic graph construction are performed with `pymatgen`; the script does not depend on JARVIS, the original ReciNet repository, YAML configuration files, or repository-internal modules.

The default model uses 92-dimensional occupancy-weighted atomic-number features, a 4.0 Å local cutoff, at most 16 neighbors per atom, the 16 shortest nonzero reciprocal vectors, four local/reciprocal interaction layers, a hidden dimension of 256, and a reciprocal down-projection dimension of 64.

## Requirements

```bash
pip install torch torch-geometric pymatgen numpy pandas openpyxl scikit-learn matplotlib tqdm
```

`torch-scatter` is optional. When it is available, the script uses `torch_scatter.scatter_add`; otherwise, it automatically falls back to PyTorch `index_add_`.

Install Optuna only when `USE_OPTUNA = True`:

```bash
pip install optuna
```

## Data

Prepare the following files in the script directory:

```text
mReciNet_v1.py
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

The first-column value is converted to an integer CIF identifier. For example, `1` corresponds to `./cif/1.cif`. Invalid targets, missing CIF files, structures without valid periodic neighbors, and CIF parsing failures are skipped and recorded in the training log.

## Usage

```bash
python mReciNet_v1.py
```

The script automatically performs CIF validation and graph caching, a fixed 80/10/10 train/validation/test split with seed 42, target normalization fitted only on the training set, model training with MSE loss, validation-RMSE early stopping, checkpoint reloading, final evaluation, and plotting. CUDA is used when available; otherwise, the script runs on CPU.

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. With the default settings, the output directory is:

```text
ReciNet_v1/
|-- ReciNet_best.pt
|-- figure/
|-- dat/
|-- table/
|-- log/
|-- split/
|-- cache/
```

The outputs include train/validation/test MAE, RMSE, and R2; individual and combined parity plots with source data; the RMSE learning curve and source data; the reusable fixed split; cached crystal graphs; the best checkpoint; and the complete training log. The checkpoint stores the model configuration, graph configuration, target normalization, best epoch, validation RMSE, and model weights. The script also provides `load_trained_model()` and `predict_cifs()` for prediction on new CIF files.

Set `USE_OPTUNA = True` to optimize the batch size, learning rate, weight decay, hidden dimension, number of interaction layers, reciprocal down-projection dimension, and dropout using only the training and validation sets. The test set is not used during hyperparameter selection.

## Citation

Original ReciNet reference:

* [Nie J, Xiao P, Ji K, et al. ReciNet: Reciprocal space-aware long-range modeling for crystalline property prediction. arXiv preprint arXiv:2502.02748, 2025](https://arxiv.org/abs/2502.02748)

This work:

To be added after the paper is officially published.

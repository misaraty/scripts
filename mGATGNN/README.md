## [中文版本](https://www.misaraty.com/2026-10-06_mgatgnn/)

## mGATGNN

`mGATGNN` is a standalone PyTorch/PyTorch Geometric implementation of GATGNN for crystal band-gap prediction from CIF structures.

The code retains the main GATGNN components, including periodic crystal-graph construction, Gaussian distance expansion, augmented multi-head graph attention (AGAT), composition-guided global attention, graph pooling, and the regression head. It uses `pymatgen` for CIF parsing and periodic-neighbor construction and does not depend on the original GATGNN repository.

The original 92-dimensional CGCNN element descriptors are embedded directly in the script. Therefore, an external `atom_init.json` file and the CGCNN repository are not required.

## Model

The default model uses 12 periodic neighbors per atom, 41 Gaussian distance features, three AGAT layers, four attention heads, 64 hidden features, and a 103-dimensional elemental-composition vector for global attention. The common graph and model settings can be edited near the beginning of the script.

Training uses MSE loss, AdamW, target standardization, validation-RMSE model selection, `ReduceLROnPlateau`, early stopping, gradient clipping, and BF16 automatic mixed precision on supported CUDA devices. Optuna hyperparameter optimization is available as an optional switch.

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
mGATGNN_v1.py
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

The value `1` in the first column corresponds to `./cif/1.cif`. The `.cif` suffix is added automatically when it is omitted.

## Usage

```bash
python mGATGNN_v1.py
```

The script automatically performs CIF validation, graph caching, a fixed 80/10/10 train/validation/test split, model training, early stopping, evaluation, and plotting. CUDA is used when available; otherwise, the script runs on CPU.

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. With the default settings, the output directory is:

```text
GATGNN_v1/
|-- GATGNN_best.pt
|-- figure/
|-- dat/
|-- table/
|-- log/
|-- split/
|-- cache/
```

The outputs include train/validation/test MAE, RMSE, and R2; individual and combined parity plots and their data; the RMSE learning curve and its data; the fixed data split; the best checkpoint; cached crystal graphs; and the complete training log.

The saved checkpoint contains the model state, model and graph configurations, target normalization, embedded atom descriptors, best epoch, best validation RMSE, target metadata, and training parameters. The script also provides `load_trained_model()` and `predict_cifs()` for prediction on new CIF structures.

## Citation

Original GATGNN reference:

* [Louis S Y, Zhao Y, Nasiri A, et al. Graph convolutional neural networks with global attention for improved materials property prediction. Physical Chemistry Chemical Physics, 2020, 22(32): 18141-18148](https://pubs.rsc.org/cp/article-abstract/22/32/18141/679461/Graph-convolutional-neural-networks-with-global)

This work:

To be added after the paper is officially published.

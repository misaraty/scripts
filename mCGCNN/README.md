## [中文版本](https://www.misaraty.com/2026-02-28_mcgcnn/)

## mCGCNN

`mCGCNN` is a standalone PyTorch implementation of the Crystal Graph Convolutional Neural Network (CGCNN) for crystal band-gap prediction from CIF structures.

The code retains the main CGCNN components, including the original 92-dimensional elemental descriptors, periodic crystal-neighbor construction, Gaussian distance expansion, gated crystal graph convolution, crystal-level mean pooling, and a scalar regression head. The original `atom_init.json` descriptors are compressed and embedded directly in the script, so no external feature file or original CGCNN repository is required. `pymatgen` is used for CIF parsing and periodic-neighbor construction.

## Model and Training

The default configuration uses:

| Setting | Default |
| --- | ---: |
| Cutoff radius | 8.0 Å |
| Maximum neighbors per atom | 12 |
| Gaussian distance step | 0.2 Å |
| Atom hidden dimension | 64 |
| CGCNN convolution layers | 3 |
| Crystal hidden dimension | 128 |
| Batch size | 64 |
| Maximum epochs | 300 |
| Learning rate | 1.0 × 10⁻³ |
| Weight decay | 1.0 × 10⁻⁵ |
| Early-stopping patience | 40 epochs |

The script uses `MSELoss` for optimization, AdamW, `ReduceLROnPlateau`, and validation RMSE for early stopping and best-checkpoint selection. Target normalization is fitted only on the training set. CUDA and automatic mixed precision are used when available; otherwise, the script runs on CPU.

## Requirements

```bash
pip install torch pymatgen numpy pandas openpyxl scikit-learn matplotlib
```

Install Optuna only when `USE_OPTUNA = True`:

```bash
pip install optuna
```

PyTorch Geometric, DGL, JARVIS, and the original CGCNN package are not required.

## Data

Prepare the following files in the script directory:

```text
mCGCNN_v2.py
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

The value `1` in the first column corresponds to `./cif/1.cif`. CIF names may be supplied with or without the `.cif` suffix. Integer-like spreadsheet values such as `1.0` are normalized to `1.cif`.

Only ordered structures are accepted. Missing CIF files, invalid targets, unreadable structures, disordered structures, and partially occupied structures are skipped and recorded in the training log.

## Usage

```bash
python mCGCNN_v2.py
```

The script automatically performs CIF validation and graph caching, a fixed 80/10/10 train/validation/test split with `SEED = 42`, model training, early stopping, evaluation, and plotting. A compatible existing split file is reused to keep repeated runs consistent.

Optuna is disabled by default. To optimize the CGCNN dimensions, number of convolution layers, number of hidden layers, learning rate, and weight decay, set:

```python
USE_OPTUNA = True
```

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. With the default settings, the output directory is:

```text
CGCNN_v2/
|-- CGCNN_best.pt
|-- figure/
|   |-- CGCNN_parity_train.jpg
|   |-- CGCNN_parity_val.jpg
|   |-- CGCNN_parity_test.jpg
|   |-- CGCNN_parity_all.jpg
|   `-- CGCNN_rmse_curve.jpg
|-- dat/
|-- table/
|   `-- CGCNN_metrics.dat
|-- log/
|   `-- CGCNN_training.log
|-- split/
|   `-- CGCNN_split.csv
`-- cache/
```

The checkpoint contains the model state, model and graph configurations, target normalization statistics, best epoch, validation RMSE, training history, target metadata, and embedded atomic descriptors. The script also provides `load_trained_model()` and `predict_cifs()` for prediction on new CIF structures.

## v1

This script implements a Crystal Graph Convolutional Neural Network (`CGCNN`) training framework in `PyTorch` for predicting target properties of crystalline materials, such as the superconducting critical temperature `Tc`. It reads structure IDs and labels from `data.xlsx`, loads the corresponding structure files from the `cif` directory, and constructs different types of atom or edge level features according to `CIFDATA_VARIANT`, including multiple encoding schemes for the electronegativity difference `Δχ`. The workflow supports either a `train/test` split or `KFold` cross validation, integrates `EarlyStopping` and learning rate scheduling, and reports `MAE`, `RMSE`, and `R2` metrics while saving the result files. When `USE_OPTUNA` is enabled, the script automatically performs hyperparameter optimization and saves the model checkpoints and logs for each trial.

Run: `python mCGCNN_v1.py`

* `CIFDATA_VARIANT`: Feature construction mode (`origin` for baseline features, `atom` for adding atomic electronegativity features, `edge` for adding edge electronegativity difference features).

* `DELTA_EN_FEAT_MODE`: Encoding scheme for `Δχ` features (`raw`, `poly`, `rbf`, `fourier`, or `all`).

* `n_folds`: Number of cross-validation folds (an integer for `KFold`, or `'none'` for a single split).

* `train_ratio` / `test_ratio`: Proportions for splitting the dataset into training and testing subsets.

* `batch_size`: Number of samples per gradient update.

* `lr`: Initial learning rate.

* `epochs`: Maximum number of training epochs.

* `patience`: Patience parameter for `EarlyStopping`.

* `atom_fea_len`: Embedding dimension for atom features.

* `h_fea_len`: Hidden dimension of the fully connected layers.

* `n_conv`: Number of graph convolution layers.

* `n_h`: Number of fully connected hidden layers.

* `USE_OPTUNA`: Whether to enable automatic hyperparameter optimization with `Optuna`.

* `OPTUNA_TRIALS`: Number of trials for the `Optuna` search.

## Citation

Original CGCNN reference:

* [Xie T, Grossman J C. Crystal graph convolutional neural networks for an accurate and interpretable prediction of material properties. Physical review letters, 2018, 120(14): 145301](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.120.145301)

This work:

v2:

To be added after the paper is officially published.

v1:

* [Zhang, Z.*; Liu, Y.; Liu, J.; Zhang, W.; Xiong, Q. Electronegativity Informed Graph Neural Networks for Superconducting Temperature Prediction with Generative Crystal Validation. Inorg. Chem. 2026, 65, 9625-9632](https://pubs.acs.org/doi/full/10.1021/acs.inorgchem.6c01169)

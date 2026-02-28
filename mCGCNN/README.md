## [中文版本](https://www.misaraty.com/2026-02-28_mcgcnn/)

## Modified CGCNN (mCGCNN)

This script implements a Crystal Graph Convolutional Neural Network (`CGCNN`) training framework in `PyTorch` for predicting target properties of crystalline materials, such as the superconducting critical temperature `Tc`. It reads structure IDs and labels from `data.xlsx`, loads the corresponding structure files from the `cif` directory, and constructs different types of atom or edge level features according to `CIFDATA_VARIANT`, including multiple encoding schemes for the electronegativity difference `Δχ`. The workflow supports either a `train/test` split or `KFold` cross validation, integrates `EarlyStopping` and learning rate scheduling, and reports `MAE`, `RMSE`, and `R2` metrics while saving the result files. When `USE_OPTUNA` is enabled, the script automatically performs hyperparameter optimization and saves the model checkpoints and logs for each trial.

## Usage

Run: `python mCGCNN_v22.py`

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

To be added after the paper is officially published.
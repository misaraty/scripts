## [中文版本](https://www.misaraty.com/2026-10-06_mchemprop/)

## mChemprop

`mChemprop` is a standalone PyTorch implementation of a Chemprop-style directed message passing neural network (D-MPNN) for molecular property regression from SMILES.

The code includes RDKit molecular parsing, Chemprop-style atom and bond features, directed bond message passing with reverse-edge exclusion, molecular aggregation, and a feed-forward regression head. It does not depend on the official Chemprop package or any internal module from the original repository.

## Requirements

```bash
pip install torch rdkit numpy pandas openpyxl scikit-learn matplotlib
```

Install Optuna only when `USE_OPTUNA = True`:

```bash
pip install optuna
```

## Data

Prepare the following files in the script directory:

```text
mChemprop_v2.py
data.xlsx
```

The first two columns of `data.xlsx` are used. The first row must contain column headers:

| SMILES | target |
| --- | ---: |
| CCO | 1.23 |
| CC(=O)O | 0.87 |
| c1ccccc1 | 2.15 |

The first column contains SMILES strings, and the second column contains numeric regression targets. Invalid SMILES, empty values, NaN values, and infinite targets are skipped and recorded in the training log. At least 10 valid molecules are required.

Set the property name and unit at the beginning of the script:

```python
TARGET_NAME = "target"
TARGET_UNIT = "unit"
```

## Usage

```bash
python mChemprop_v2.py
```

The script automatically validates and canonicalizes SMILES, performs a fixed 80/10/10 train/validation/test split, standardizes targets using training data only, trains the model with MSE loss, applies validation-RMSE early stopping, reloads the best checkpoint, evaluates all three subsets, and generates plots and data tables. CUDA and automatic mixed precision are used when available; otherwise, the script runs on CPU.

The saved split is reused when it matches the current dataset. Optional Optuna optimization uses only the training and validation sets and is disabled by default.

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. With the default settings, the output directory is:

```text
Chemprop_v2/
|-- Chemprop_best.pt
|-- figure/
|   |-- Chemprop_parity_train.jpg
|   |-- Chemprop_parity_val.jpg
|   |-- Chemprop_parity_test.jpg
|   |-- Chemprop_parity_all.jpg
|   `-- Chemprop_rmse_curve.jpg
|-- dat/
|-- table/
|-- log/
`-- split/
```

The outputs include train/validation/test MAE, RMSE, and R2; parity plots and their DAT files; the RMSE learning curve and its DAT file; the reusable data split; the best checkpoint; and the complete training log. When Optuna is enabled, its best parameters are saved in `table/Chemprop_optuna_best_params.json`.

The trained checkpoint can also be loaded from another Python program with `load_trained_model()`, and new molecules can be predicted with `predict_smiles()`.

## Citation

Chemprop v2:

* [Graff D E, Morgan N K, Burns J W, et al. Chemprop v2: An efficient, modular machine learning package for chemical property prediction[J]. Journal of Chemical Information and Modeling, 2026, 66(1): 28-33](https://pubs.acs.org/jcisd8/article-abstract/66/1/28/5080274/Chemprop-v2-An-Efficient-Modular-Machine-Learning)

D-MPNN theory:

* [Yang K, Swanson K, Jin W, et al. Analyzing learned molecular representations for property prediction[J]. Journal of chemical information and modeling, 2019, 59(8): 3370](https://pubs.acs.org/jcisd8/article/59/8/3370/863796/Analyzing-Learned-Molecular-Representations-for)

This work:

To be added after the paper is officially published.

## [中文版本](https://www.misaraty.com/2026-10-06_mchemprop/)

## mChemprop v3

`mChemprop_v3.py` is a standalone PyTorch implementation of a Chemprop-style directed message passing neural network (D-MPNN) for single-task and multi-task molecular property regression from SMILES.

The script includes RDKit molecular parsing, Chemprop v2-style atom and bond features, directed bond message passing with reverse-edge exclusion, molecular aggregation, and a feed-forward regression head. It does not import the official Chemprop package or any internal module from the original repository.

The single-task and multi-task modes use the same D-MPNN architecture. A single-task dataset produces one output, while a multi-task dataset produces one output per target column and shares the molecular representation across all tasks.

## Main features

- Automatic single-task or multi-task regression from one Excel file
- First Excel column used as SMILES; all remaining columns used as regression targets
- Partial target missingness supported through masked MSE
- Independent target standardization fitted on the training set for every task
- Fixed and reusable 80/10/10 train/validation/test split with seed 42
- Validation-only early stopping, checkpoint selection, scheduling, and optional Optuna optimization
- Per-task MAE, RMSE, and R2 for the training, validation, and test sets
- Per-task parity plots, prediction tables, and RMSE curves
- Automatic CUDA selection and automatic mixed precision when CUDA is available
- Importable `load_trained_model()` and `predict_smiles()` functions
- Version-isolated output directories controlled by `RUN_VERSION`

## Requirements

```bash
pip install torch rdkit numpy pandas openpyxl scikit-learn matplotlib
```

Install Optuna only when `USE_OPTUNA = True`:

```bash
pip install optuna
```

## Data

Place the script and Excel file in the same directory:

```text
mChemprop_v3.py
data.xlsx
```

The first row of `data.xlsx` must contain unique, non-empty column names. The first column contains SMILES strings, and every remaining column is treated as a continuous regression task.

Single-task example:

| SMILES | Activity |
| --- | ---: |
| CCO | 1.23 |
| CC(=O)O | 0.87 |
| c1ccccc1 | 2.15 |

Multi-task example:

| SMILES | MIC_Ecoli | MIC_Saureus | Cytotoxicity |
| --- | ---: | ---: | ---: |
| CCO | 8.0 | 16.0 | 120.0 |
| CCN | 4.0 |  | 85.0 |
| c1ccccc1 |  | 32.0 | 96.0 |

Missing task labels may be left empty or stored as `NaN`. A molecule is retained when at least one target is valid. Invalid SMILES and rows without any valid target are skipped and recorded in the log. Each task must have at least one valid label in both the training and validation sets. At least 10 valid molecules are required for the fixed split.

Optional units are assigned by exact Excel column name near the beginning of the script:

```python
TARGET_UNITS = {
    "MIC_Ecoli": "ug/mL",
    "MIC_Saureus": "ug/mL",
    "Cytotoxicity": "uM",
}
```

Leave `TARGET_UNITS = {}` when units are not needed.

## Training and model selection

For every target, the mean and standard deviation are fitted using only its available training labels. Training uses masked MSE in the per-task standardized target space, so missing labels do not contribute to the loss.

At the end of each epoch, RMSE is calculated separately for every task in its original unit. For multi-task runs, the arithmetic mean of the per-task validation RMSE values is used for:

- early stopping;
- best-epoch and checkpoint selection;
- learning-rate scheduling;
- the Optuna objective.

The test set is not used during optimization or model selection. Because the macro RMSE directly averages per-task RMSE values, targets with very different physical units or numerical scales should be interpreted carefully.

## Usage

```bash
python mChemprop_v3.py
```

The script validates and canonicalizes SMILES, creates or reuses the fixed split, trains the D-MPNN, saves the best checkpoint, reloads that checkpoint, evaluates all three subsets, and generates the complete report. CUDA is selected automatically when available; otherwise, the script runs on CPU.

Optuna is disabled by default:

```python
USE_OPTUNA = False
```

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. The default directory is:

```text
Chemprop_v3/
|-- Chemprop_best.pt
|-- figure/
|-- dat/
|-- table/
|-- log/
`-- split/
```

For a single-task dataset, the main plot names remain concise:

```text
figure/Chemprop_parity_train.jpg
figure/Chemprop_parity_val.jpg
figure/Chemprop_parity_test.jpg
figure/Chemprop_parity_all.jpg
figure/Chemprop_rmse_curve.jpg
```

For a multi-task dataset, each task receives its own parity plots, DAT files, and RMSE curve. For example:

```text
figure/Chemprop_MIC_Ecoli_parity_train.jpg
figure/Chemprop_MIC_Ecoli_parity_val.jpg
figure/Chemprop_MIC_Ecoli_parity_test.jpg
figure/Chemprop_MIC_Ecoli_parity_all.jpg
figure/Chemprop_MIC_Ecoli_rmse_curve.jpg
```

The principal tabular outputs are:

```text
dat/Chemprop_parity_all_tasks.dat
dat/Chemprop_rmse_curve.dat
table/Chemprop_metrics.dat
log/Chemprop_training.log
split/Chemprop_split.csv
```

`Chemprop_metrics.dat` contains per-task metrics and a `__macro__` summary row for each split. `Chemprop_rmse_curve.dat` contains macro-averaged and per-task training/validation RMSE columns. When Optuna is enabled, the best parameters are saved as `table/Chemprop_optuna_best_params.json`.

## Prediction from Python

The best checkpoint stores the model configuration, task names, task units, number of tasks, target-scaling parameters, feature dimensions, seed, run version, and best validation result.

New molecules can be predicted from another Python program:

```python
from mChemprop_v3 import predict_smiles

predictions = predict_smiles(
    ["CCO", "CC(=O)O", "c1ccccc1"],
    "./Chemprop_v3/Chemprop_best.pt",
)
print(predictions)
```

The returned table contains one `prediction::<task_name>` column for every regression task. Invalid SMILES are retained in the returned table with missing predictions.

## Historical releases

`v2` was the standalone single-task release. It read the first Excel column as SMILES and the second column as one continuous target, then generated one checkpoint, one set of parity plots, one RMSE curve, and one set of metrics. v3 retains this single-task workflow while extending the same D-MPNN implementation to multiple targets and partially missing labels.

## Citation

Original Chemprop reference:

* [Graff D E, Morgan N K, Burns J W, et al. Chemprop v2: An efficient, modular machine learning package for chemical property prediction[J]. Journal of Chemical Information and Modeling, 2026, 66(1): 28-33](https://pubs.acs.org/jcisd8/article-abstract/66/1/28/5080274/Chemprop-v2-An-Efficient-Modular-Machine-Learning)

* [Yang K, Swanson K, Jin W, et al. Analyzing learned molecular representations for property prediction[J]. Journal of chemical information and modeling, 2019, 59(8): 3370](https://pubs.acs.org/jcisd8/article/59/8/3370/863796/Analyzing-Learned-Molecular-Representations-for)

This work:

To be added after the paper is officially published.

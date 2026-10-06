## [中文版本](https://www.misaraty.com/2026-10-06_mmodnet/)

## mMODNet

`mMODNet` is a standalone PyTorch implementation of the Material Optimal Descriptor Network (MODNet) for crystal band-gap prediction from CIF structures.

The code uses `pymatgen` to parse CIF files and `matminer` to generate composition and structure descriptors. It retains the main MODNet workflow: material descriptor generation, normalized-mutual-information relevance-redundancy (RR) feature selection, descriptor preprocessing, hierarchical dense representation blocks, and scalar regression. The script does not depend on TensorFlow, the original `modnet` package, or repository-internal modules.

The default featurization includes BandCenter, ElementFraction, Magpie and elemental-property statistics, stoichiometry, transition-metal fraction, valence-orbital fractions, density, Ewald energy, global symmetry, and structural-complexity descriptors. Feature imputation, variance filtering, RR ranking, and standardization are fitted using the training set only.

## Requirements

```bash
pip install torch pymatgen matminer numpy pandas openpyxl scikit-learn matplotlib
```

Install Optuna only when `USE_OPTUNA = True`:

```bash
pip install optuna
```

## Data

Prepare the following files in the script directory:

```text
mMODNet_v1.py
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

The value `1` in the first column corresponds to `./cif/1.cif`. The current script converts the first-column value with `int(raw_name)`, so CIF identifiers must be numeric.

## Usage

```bash
python mMODNet_v1.py
```

The script automatically performs CIF validation, descriptor generation and caching, a fixed 80/10/10 train/validation/test split, training-set-only RR feature selection, model training, early stopping, evaluation, and plotting. CUDA and BF16 automatic mixed precision are used when available; otherwise, the script runs on CPU.

The default configuration selects 128 descriptors from at most 256 RR candidates and uses hierarchical dense blocks with dimensions `(256, 128)`, `(128, 64)`, and `(64, 32)`. Training uses MSE loss, while validation RMSE controls learning-rate scheduling, early stopping, and best-checkpoint selection.

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. With the default settings, the output directory is:

```text
MODNet_v1/
|-- MODNet_best.pt
|-- figure/
|-- dat/
|-- table/
|-- log/
|-- split/
|-- cache/
```

The outputs include train/validation/test MAE, RMSE, and R2; separate and combined parity plots and data; the RMSE learning curve and data; the RR-selected descriptor ranking; the fixed data split; the descriptor cache; the best checkpoint; and the complete training log. The checkpoint stores the model configuration and preprocessing information required by the included `load_trained_model()` and `predict_cifs()` functions.

## Citation

Original MODNet reference:

* [De Breuck P P, Hautier G, Rignanese G M. Materials property prediction for limited datasets enabled by feature selection and joint learning with MODNet. npj computational materials, 2021, 7(1): 83](https://www.nature.com/articles/s41524-021-00552-2)

This work:

To be added after the paper is officially published.

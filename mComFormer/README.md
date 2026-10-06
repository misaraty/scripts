## [中文版本](https://www.misaraty.com/2026-10-06_mcomformer/)

## mComFormer

`mComFormer` is a standalone PyTorch/PyTorch Geometric implementation of iComFormer and eComFormer for crystal band-gap prediction from CIF structures. Both models are included in one Python script and can be selected through `MODEL_NAME`.

The code uses `pymatgen` for CIF parsing and periodic graph construction and does not depend on DGL, JARVIS, or the original ComFormer repository.

## Models

| Model | Geometric representation | Main components |
| --- | --- | --- |
| iComFormer | SE(3)-invariant | Interatomic distances, three canonical lattice reference vectors, angle encoding, four node-wise ComFormer layers, and one edge update layer |
| eComFormer | SO(3)-equivariant | Interatomic displacement vectors, spherical harmonics, three node-wise ComFormer layers, and one equivariant tensor-product update layer |

Both models use a learned 92-dimensional atomic embedding, radial basis distance encoding, graph-level mean pooling, and a scalar regression head.

## Requirements

Install the common dependencies:

```bash
pip install torch torch-geometric pymatgen numpy pandas openpyxl scikit-learn matplotlib tqdm
```

eComFormer additionally requires e3nn:

```bash
pip install e3nn
```

Install Optuna only when `USE_OPTUNA = True`:

```bash
pip install optuna
```

## Data

Prepare the following files in the script directory:

```text
mComFormer_v1.py
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

The value `1` in the first column is converted to `1.cif` and corresponds to `./cif/1.cif`. CIF structures containing disordered or partially occupied sites are skipped and recorded in the training log.

## Model Selection

Select one of the two models at the beginning of `mComFormer_v1.py`:

```python
MODEL_NAME = "iComFormer"
```

or:

```python
MODEL_NAME = "eComFormer"
```

The run version can be changed independently:

```python
RUN_VERSION = "v1"
```

## Usage

```bash
python mComFormer_v1.py
```

The script automatically performs CIF validation, periodic graph construction and caching, a fixed 80/10/10 train/validation/test split, target normalization based only on the training set, model training, early stopping, evaluation, and plotting. CUDA is used when available; otherwise, the script runs on CPU.

Training uses MSE loss, AdamW, and ReduceLROnPlateau. The best checkpoint is selected according to validation RMSE. The saved split is reused when the same valid CIF set is detected.

## Outputs

All results are saved in a directory determined by `MODEL_NAME` and `RUN_VERSION`. For example:

```text
iComFormer_v1/
|-- iComFormer_best.pt
|-- figure/
|-- dat/
|-- table/
|-- log/
|-- split/
|-- cache/
```

or:

```text
eComFormer_v1/
|-- eComFormer_best.pt
|-- figure/
|-- dat/
|-- table/
|-- log/
|-- split/
|-- cache/
```

The outputs include train/validation/test MAE, RMSE, and R2; individual and combined parity plots with data files; the RMSE learning curve and data; the fixed data split; cached crystal graphs; the best checkpoint; optional Optuna results; and the complete training log.

The script also provides `load_trained_model()` and `predict_cifs()` for loading a saved checkpoint and predicting band gaps for new CIF structures.

## Citation

Original ComFormer reference:

* [Yan K, Fu C, Qian X, et al. Complete and efficient graph transformers for crystal material property prediction. International Conference on Learning Representations. 2024, 2024: 2564-2590.](https://proceedings.iclr.cc/paper_files/paper/2024/hash/0ab51646ca369140c3c3ece011b66587-Abstract-Conference.html)

This work:

To be added after the paper is officially published.

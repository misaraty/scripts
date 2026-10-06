## [中文版本](https://www.misaraty.com/2026-10-06_mschnet/)

## mSchNet

`mSchNet` is a standalone PyTorch implementation of SchNet for crystal band-gap prediction from CIF structures.

The code retains the main SchNet components, including atomic-number embeddings, Gaussian radial basis expansion, continuous-filter convolution, shifted-softplus activations, residual interaction blocks, atom-wise outputs, and crystal-level pooling. Periodic neighbors are constructed directly from CIF structures with `pymatgen`. The current configuration uses a smooth cosine cutoff, which can be disabled through `USE_COSINE_CUTOFF` for comparison with a hard cutoff. The script does not depend on TensorFlow, PyTorch Geometric, ASE databases, or the original SchNet repository.

## Requirements

```bash
pip install torch pymatgen numpy pandas openpyxl scikit-learn matplotlib
```

Install Optuna only when `USE_OPTUNA = True`:

```bash
pip install optuna
```

## Data

Prepare the following files in the script directory:

```text
mSchNet_v1.py
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

The value `1` in the first column corresponds to `./cif/1.cif`.

## Usage

```bash
python mSchNet_v1.py
```

The script automatically performs CIF validation, periodic graph construction, graph caching, a fixed 80/10/10 train/validation/test split, target normalization based only on the training set, model training, early stopping, evaluation, and plotting. CUDA is used when available; otherwise, the script runs on CPU. BF16 automatic mixed precision is enabled by default on CUDA devices.

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. With the current default settings, the output directory is:

```text
SchNet_v1_USE_COSINE_CUTOFF/
|-- SchNet_best.pt
|-- figure/
|-- dat/
|-- table/
|-- log/
|-- split/
|-- cache/
```

The outputs include train/validation/test MAE, RMSE, and R2; individual and combined parity plots and data; the train/validation RMSE learning curve and data; the reusable fixed data split; the cached crystal graphs; the best checkpoint selected by validation RMSE; and the complete training log.

## Citation

Original SchNet reference:

* [Schütt K, Kindermans P J, Sauceda Felix H E, et al. Schnet: A continuous-filter convolutional neural network for modeling quantum interactions. Advances in neural information processing systems, 2017, 30](https://proceedings.neurips.cc/paper_files/paper/2017/hash/303ed4c69846ab36c2904d3ba8573050-Abstract.html)

This work:

To be added after the paper is officially published.

## [中文版本](https://www.misaraty.com/2026-10-06_mcongn/)

`mcoNGN_v1.py` is a standalone PyTorch/PyTorch Geometric program for crystal-property regression from CIF structures. It contains both **coGN** and **coNGN**; select the model by changing the single `MODEL_NAME` option at the top of the file.

The default target is the band gap in eV. The target name, unit, data paths, graph construction, network size, training schedule, plotting style, and optional Optuna search can all be changed in the user-configuration section.

## Models

| Setting | coGN | coNGN |
| --- | --- | --- |
| `MODEL_NAME` | `"coGN"` | `"coNGN"` |
| Crystal graph | Adaptive periodic 24-nearest-neighbor graph | Adaptive periodic Voronoi graph |
| Edge inputs | Interatomic-distance RBF | Distance RBF + Voronoi ridge-area RBF |
| Message passing | Connectivity-optimized graph-network blocks | Connectivity-optimized blocks with nested line-graph updates |
| Angular information | Not used | Bond-angle RBF on the line graph |
| Default hidden width | 128 | 160 |
| Default output head | Linear | Two-layer nonlinear MLP |

Both variants use learned atomic embeddings together with fixed elemental descriptors, five processing blocks, global mean pooling, and a scalar regression output.

## Requirements

Python 3.10 or later is recommended. Install a PyTorch build suitable for the local CPU/CUDA environment, followed by the remaining packages:

```bash
pip install torch torch-geometric pymatgen numpy pandas openpyxl scikit-learn matplotlib tqdm
```

Optuna is required only when `USE_OPTUNA = True`:

```bash
pip install optuna
```

## Data

Place `data.xlsx`, the script, and the `cif` directory in the same working directory:

```text
.
├── mcoNGN_v1.py
├── data.xlsx
└── cif
    ├── 1.cif
    ├── 2.cif
    └── ...
```

Only the first two columns of `data.xlsx` are read. The first column must contain an integer CIF identifier; for example, `1` is mapped to `./cif/1.cif`. The second column contains the regression target.

| CIF identifier | bandgap |
| ---: | ---: |
| 1 | 1.237 |
| 2 | 2.104 |
| ... | ... |

The program removes rows with invalid targets, missing/unreadable CIFs, or structures that cannot produce a valid graph. At least 10 valid graphs must remain before training can continue.

## Usage

Choose one model in the user-configuration section:

```python
MODEL_NAME = "coGN"
# or
MODEL_NAME = "coNGN"
```

Then run:

```bash
python mcoNGN_v1.py
```

By default, the code uses a deterministic 80%/10%/10% train/validation/test split (`SEED = 42`), standardizes targets with training-set statistics, trains with AdamW and MSE loss, reduces the learning rate on a validation plateau, and applies early stopping using validation RMSE. CUDA is used automatically when available; automatic mixed precision is enabled on CUDA.

The split file and graph cache are reused on later runs. When changing the dataset or graph-related settings, use a new `RUN_VERSION` or move the corresponding old output directory first.

## Output

The selected model controls the output directory:

```text
coGN_v1/                         # or coNGN_v1/
├── coGN_best.pt                 # selected-model checkpoint
├── cache/coGN_graphs/           # cached graphs
├── dat/
│   ├── coGN_parity_train.dat
│   ├── coGN_parity_val.dat
│   ├── coGN_parity_test.dat
│   ├── coGN_parity_all.dat
│   └── coGN_rmse_curve.dat
├── figure/
│   ├── coGN_parity_train.jpg
│   ├── coGN_parity_val.jpg
│   ├── coGN_parity_test.jpg
│   ├── coGN_parity_all.jpg
│   └── coGN_rmse_curve.jpg
├── log/coGN_training.log
├── split/coGN_split.csv
└── table/coGN_metrics.dat
```

When `MODEL_NAME = "coNGN"`, every `coGN` prefix in the example above is replaced with `coNGN`. If Optuna is enabled, `table/<MODEL_NAME>_optuna_best_params.json` is also written.

The checkpoint includes the trained weights, model and graph configurations, target-normalization statistics, best epoch, best validation RMSE, seed, and training parameters. The supplied `predict_cifs()` function can load this checkpoint and predict new CIF files.

## Citation

Original ALIGNN reference:

* [Ruff R, Reiser P, Stühmer J, et al. Connectivity optimized nested line graph networks for crystal structures. Digital Discovery, 2024, 3(3): 594-601.](https://pubs.rsc.org/dd/article/3/3/594/846022/Connectivity-optimized-nested-line-graph-networks)

This work:

To be added after the paper is officially published.

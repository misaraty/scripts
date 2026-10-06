## [中文版本](https://www.misaraty.com/2026-02-27_mmegnet/)

## mMEGNet

`mMEGNet` is a standalone PyTorch/PyTorch Geometric implementation of MEGNet for crystal band-gap prediction from CIF structures.

The code retains the main MEGNet components, including the periodic crystal graph, atomic embeddings, Gaussian distance expansion, edge/node/global state updates, residual MEGNet blocks, dual Set2Set pooling, and regression head. It uses `pymatgen` for CIF parsing and periodic graph construction and does not depend on TensorFlow, Keras, or the original MEGNet repository.

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
mMEGNet_v2.py
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
python mMEGNet_v2.py
```

The script automatically performs CIF validation, periodic graph caching, a fixed 80/10/10 train/validation/test split, target normalization using the training set, model training, early stopping, evaluation, and plotting. CUDA is used when available; otherwise, the script runs on CPU. BF16 automatic mixed precision is enabled on supported CUDA devices.

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. With the default settings, the output directory is:

```text
MEGNet_v2/
|-- MEGNet_best.pt
|-- figure/
|-- dat/
|-- table/
|-- log/
|-- split/
|-- cache/
```

The outputs include train/validation/test MAE, RMSE, and R2; parity plots and data; the RMSE learning curve and data; the fixed data split; cached crystal graphs; the best checkpoint; and the complete training log.

## v1

This script implements a MEGNet graph neural network pipeline to predict the superconducting critical temperature Tc of crystalline materials. It reads structure IDs and Tc labels from data.xlsx, loads the corresponding structure files from the cif directory, and optionally adds electronegativity mean and standard deviation as global features. The data are split using either a train test scheme or KFold cross validation, and the model is trained with EarlyStopping and learning rate scheduling. The script outputs RMSE, MAE, and R2 metrics, saves prediction scatter plots, and in megnet_tuned mode uses Optuna to optimize network widths and learning rate while automatically saving logs and model files.

The execution commands are `python mMEGNet_v14.py`.

* `TARGET_COL`: Name of the prediction target column (e.g., `'tc'` for superconducting critical temperature).

* `method`: Training mode selection (`'megnet_default'` for fixed architecture, `'megnet_tuned'` for Optuna hyperparameter optimization).

* `batch_size`: Number of samples per gradient update affecting convergence stability and memory usage.

* `lr`: Initial learning rate determining optimization step size.

* `train_ratio` / `test_ratio`: Proportions for dataset splitting between training and testing.

* `USE_EN_GLOBAL`: Whether to include electronegativity-based global state features.

* `n1, n2, n3`: Control the width (capacity) of the MEGNet network; larger values increase model complexity but may risk overfitting.

* `epochs`: Maximum number of training epochs; sets the upper limit of training iterations.

* `n_folds`: Number of cross-validation folds; controls robustness of model evaluation (`5` recommended, or `'none'` for single split).


> [!NOTE]
> Replace the `get_atom_features` function in `~\anaconda3\Lib\site-packages\megnet\data\graph.py` with the following implementation:
>
> ```python
>     @staticmethod
>     def get_atom_features(structure) -> List[Any]:
>         z_list = []
>         for site in structure:
>             try:
>                 if getattr(site, "is_ordered", True):
>                     z = int(site.specie.Z)
>                 else:
>                     items = list(site.species.items())  # [(Element, occ), ...]
>                     items.sort(key=lambda kv: (float(kv[1]), getattr(kv[0], "Z", 0)), reverse=True)
>                     z = int(getattr(items[0][0], "Z", 0))
>             except Exception:
>                 tok = str(site.species_string).split()[0]
>                 try:
>                     from pymatgen.core.periodic_table import Element
>                     z = int(Element(tok).Z)
>                 except Exception:
>                     z = 0
>             z_list.append(z)
>         return np.array(z_list, dtype="int32").tolist()
> ```
>
> This modification prevents errors when reading CIF files containing partially occupied (disordered) atomic sites by selecting the species with the highest occupancy for each site.

## Citation

Original MEGNet reference:

* [Chen C, Ye W, Zuo Y, et al. Graph networks as a universal machine learning framework for molecules and crystals. Chemistry of Materials, 2019, 31(9): 3564-3572](https://pubs.acs.org/cmatex/article-abstract/31/9/3564/1288793/Graph-Networks-as-a-Universal-Machine-Learning)

This work:

v2:

To be added after the paper is officially published.

v1:

* [Zhang, Z.*; Liu, Y.; Liu, J.; Zhang, W.; Xiong, Q. Electronegativity Informed Graph Neural Networks for Superconducting Temperature Prediction with Generative Crystal Validation. Inorg. Chem. 2026, 65, 9625-9632](https://pubs.acs.org/doi/full/10.1021/acs.inorgchem.6c01169)
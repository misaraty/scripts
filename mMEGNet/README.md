## [中文版本](https://www.misaraty.com/2026-02-20_imod/)

## modified MEGNet (mMEGNet)

This script implements a MEGNet graph neural network pipeline to predict the superconducting critical temperature Tc of crystalline materials. It reads structure IDs and Tc labels from data.xlsx, loads the corresponding structure files from the cif directory, and optionally adds electronegativity mean and standard deviation as global features. The data are split using either a train test scheme or KFold cross validation, and the model is trained with EarlyStopping and learning rate scheduling. The script outputs RMSE, MAE, and R2 metrics, saves prediction scatter plots, and in megnet_tuned mode uses Optuna to optimize network widths and learning rate while automatically saving logs and model files.

## Usage

* The execution commands are `python mMEGNet_v14.py`.

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

To be added after the paper is officially published.
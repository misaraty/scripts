## [中文版本](https://www.misaraty.com/2026-10-06_flycrysnet/)

## FlyCrysNet

`FlyCrysNet` is a standalone PyTorch/PyTorch Geometric model for crystal band-gap prediction from CIF structures. It combines a periodic crystal graph encoder with sparse latent message passing derived from the complete *Drosophila* male central nervous system (MaleCNS) connectome.

The crystal encoder uses atomic embeddings, radial-basis distance encoding, four edge-gated graph-convolution layers, and global mean pooling. The pooled crystal representation is projected onto 256 latent nodes with eight channels per node. A MaleCNS-derived directed weighted topology then performs adaptively gated sparse propagation. The topology gate and synaptic-weight exponent are learnable, and attention pooling plus a global skip pathway produces the final band-gap prediction.

The script also includes four necessary topology controls: self-only, random topology, degree-preserving randomization, and unweighted MaleCNS topology. All variants use the same fixed data split and training settings, while their best checkpoints are selected independently by validation RMSE.

## Requirements

```bash
pip install torch torch-geometric pymatgen numpy pandas openpyxl scikit-learn matplotlib tqdm pyarrow
```

Install Optuna only when `USE_OPTUNA = True`:

```bash
pip install optuna
```

## Data

Prepare the following files and directories in the same directory as the script:

```text
FlyCrysNet_v3.py
data.xlsx
cif/
|-- 1.cif
|-- 2.cif
|-- 3.cif
fly_connectome/
|-- connectome-weights-male-cns-v1.0-minconf-0.5.feather
```

The `connectome-weights-male-cns-v1.0-minconf-0.5.feather` file is approximately 1.1 GB. Keep its filename unchanged and place it in `./fly_connectome/`.

The first two columns of `data.xlsx` are used:

| cif | bandgap |
| ---: | ---: |
| 1 | 1.23 |
| 2 | 0.87 |
| 3 | 2.15 |

The script converts each value in the first column to an integer filename. For example, `1` corresponds to `./cif/1.cif`.

## MaleCNS topology

On the first run, the script scans the Feather file twice. It first ranks neurons by weighted degree and selects the top 256 neurons, then extracts their directed weighted induced subgraph. Positive synaptic weights are transformed with `log1p`, self-loops are added, and the resulting sparse topology is normalized and cached.

The cached topology is saved as:

```text
FlyCrysNet_v3/cache/MaleCNS_top256.npz
```

Subsequent runs reuse this cache when the source Feather file, its modification time, and `CONNECTOME_NODES` remain unchanged. The first run is therefore substantially slower than later runs.

## Usage

```bash
python FlyCrysNet_v3.py
```

The script automatically performs CIF validation, periodic graph construction, a fixed 80/10/10 train/validation/test split, target standardization, model training, early stopping, evaluation, plotting, and the topology ablations. CUDA with bfloat16 automatic mixed precision is used when available; otherwise, the script runs on CPU.

Set `RUN_ABLATIONS = False` to run only the full weighted MaleCNS model. Set `USE_OPTUNA = True` to enable the optional hyperparameter search.

All results are saved in the directory determined by `MODEL_NAME` and `RUN_VERSION`. With the default settings, the output directory is:

```text
FlyCrysNet_v3/
|-- FlyCrysNet_best.pt
|-- figure/
|-- dat/
|-- table/
|-- log/
|-- split/
|-- cache/
|   |-- MaleCNS_top256.npz
|   `-- crystal_graphs/
`-- ablation_checkpoints/
```

The outputs include train/validation/test MAE, RMSE, and R2; parity plots and their data; the RMSE learning curve and its data; the fixed data split; the best full-model checkpoint; ablation checkpoints; topology-control metrics, predictions, learning curves, and figures; cached crystal graphs and MaleCNS topology; and the complete training log.

## Citation

MaleCNS connectome reference:

* [Berg S, Beckett I R, Costa M, et al. Sexual dimorphism in the complete Drosophila male central nervous system connectome. Cell, 2026, 189(18): 5504-5526. e15](https://www.cell.com/cell/fulltext/S0092-8674(26)00942-6)

This work:

To be added after the paper is officially published.

## **[中文版本](https://www.misaraty.com/2025-12-19_matlab%E7%89%88%E6%9C%ACsisso/)**

## MATLAB Implementation of SISSO

A MATLAB implementation of the `Sure Independence Screening and Sparsifying Operator (SISSO)` method for symbolic regression based on compressed sensing.

The code is adapted from the official MATLAB repository
[SISSORegressor_MATLAB](https://github.com/NREL/SISSORegressor_MATLAB)
with additional modifications for descriptor construction, data splitting, and model evaluation.

## Design Philosophy

- explicit implementation of `Sure Independence Screening (SIS)` and `Sparsifying Operator (SO)` without black-box solvers

- single-file `.m` scripts for easy inspection and modification

- explicit evaluation of model generalization performance

## MATLAB Versions

### mSISSO_MATLAB_v1

- A minimal single-file MATLAB implementation derived from [SISSORegressor_MATLAB](https://github.com/NREL/SISSORegressor_MATLAB).

- Output:

```shell
Loaded data_v4.csv: 3624 samples, 8 base features.
Generated descriptor library: 3624 samples, 28747 features.
Split: train=2537, val=544, test=543
Searching for models up to 3 dimensions, considering 26 new features per iteration.
          RMSE            Model
1D:  0.320931  0.599 - 0.227 ((χ_X^2)/((χ_B'^(1/4))+eps))
2D:  0.197306  0.256 - 0.232 ((χ_X^2)/((χ_B'^(1/4))+eps)) + 0.448 (χ_A+χ_B''+(-r_A))
3D:  0.178237  0.065 - 0.232 ((χ_X^2)/((χ_B'^(1/4))+eps)) + 0.443 (χ_A+χ_B''+(-r_A)) + 0.133 ((1/(r_B'+eps))-log(r_B''+eps))

=== Metrics (Model Dim = 3) ===
Train: MAE=0.1212, RMSE=0.1782, R2=0.9455
Val  : MAE=0.1228, RMSE=0.1688, R2=0.9536
Test : MAE=0.1225, RMSE=0.1613, R2=0.9540
```

### mSISSO_MATLAB_v2

- An extended version of v1, also provided as a single `.m` file.

- Inherits the core SISSO workflow (standardization, SIS, and SO).

- Adds:

  - more flexible and robust descriptor construction

  - explicit Train / Validation / Test data splitting

  - standardized reporting of R2, MAE, and RMSE for generalization assessment

## Fortran Implementation of SISSO

- [SISSO](https://github.com/rouyang2017/SISSO)

  The official and most complete implementation, recommended for high-throughput applications.

## Python Implementation of SISSO

### pysisso

- [pysisso](https://github.com/Matgenix/pysisso)

### tutorial-compressed-sensing

- [Symbolic regression via compressed sensing: a tutorial (NOMAD Notebook)](https://nomad-lab.eu/prod/analytics/public/user/misaraty/notebooks/tutorials/compressed_sensing.ipynb)

- [tutorial-compressed-sensing (GitLab mirror by Luigi Sbailo)](https://gitlab.mpcdf.mpg.de/nomad-lab/ai-toolkit/tutorial-compressed-sensing)

## C++ Implementation of SISSO

- [sissopp](https://gitlab.com/sissopp_developers/sissopp)

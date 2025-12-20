## **[中文版本](https://www.misaraty.com/2025-12-19_matlab%E7%89%88%E6%9C%ACsisso/)**

## Modified MATLAB Implementation of SISSO (mSISSO_MATLAB)

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
Fitting 'small' data: 82 data points, 115 features.
Searching for models up to 3 dimemsions, considering 10 new features per iteration.
          RMSE            Model
1D:	0.296696	1.922 - 0.478 (r_p(A)+r_d(B)) 
2D:	0.218070	7.495 - 3.483 (r_p(A)+r_d(B)) + 0.392 (r_p(A)+r_d(B))^2 
3D:	0.193928	7.280 - 3.528 (r_p(A)+r_d(B)) + 0.405 (r_p(A)+r_d(B))^2 + 0.293 |r_s(A)-r_d(B)| 
 
Fitting 'big' data: 82 data points, 3391 features.
Searching for models up to 3 dimemsions, considering 26 new features per iteration.
          RMSE            Model
1D:	0.137310	-0.327 - 0.055 (IP(A)+IP(B))/r_p(A)^2 
2D:	0.100216	-0.145 + 0.114 |IP(B)-EA(B)|/r_p(A)^2 - 1.482 |r_s(A)-r_p(B)|/exp(r_s(A)) 
3D:	0.076428	-0.005 + 0.109 |IP(B)-EA(B)|/r_p(A)^2 - 1.766 |r_s(A)-r_p(B)|/exp(r_s(A)) - 6.032 |r_s(B)-r_p(B)|/(r_p(B)+r_d(A))^2 
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

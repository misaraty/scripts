## **[中文版本](https://www.misaraty.com/2025-12-21_xlrf/)**

## XLRF

Within a unified tree-ensemble framework, MATLAB-based XGBoost-/LightGBM-style gradient and random forest models (XLRF) are constructed and evaluated.

## XLRF_v1.0

* In the `xgboost-style` and `lightgbm-style` implementations, native XGBoost/LightGBM libraries are not used. Instead, MATLAB's  
  `fitrensemble (LSBoost)` combined with `decision-tree` templates is employed to emulate their behavior, mapping only high-level concepts such as learning rate, number of boosting iterations, tree depth/leaf count, and row/column subsampling.

* The random forest model is implemented using `TreeBagger`, which is MATLAB’s standard realization and is relatively close to the original random forest algorithm.

> [!NOTE]
> * When running MATLAB code, avoid using decimal points in `.m` script filenames (e.g., `XLRF_v1.0.m`). Use integer or underscore formats instead (e.g., `XLRF_v1.m` or `XLRF_v1_0.m`).
> * MATLAB allows only a single dot for the file extension; additional dots can cause script name parsing errors and lead to execution failure.

## XLRF_v1.1

* Compared with `v1.0`, model interpretability is computed using MATLAB's built-in importance measures:  
  `predictorImportance` for GBDT models, and `OOBPermutedPredictorDeltaError` for random forests by enabling  
  `OOBPrediction` and `OOBPredictorImportance`.

## XLRF_v1.2

* Compared with `v1.0`, custom permutation importance is introduced for interpretability analysis, evaluating feature relevance via test-set RMSE increases under repeated feature shuffling.

## XLRF_v2.0

* Building on `v1.0`, validation-based optimal boosting iteration selection is introduced.  
  Models are first trained with a large maximum number of iterations, and the optimal `NumLearningCycles` is then determined using a validation set, providing more robust control of overfitting.

## XLRF_v2.1

* Compared with `v2.0`, custom permutation importance is used for interpretability evaluation instead of built-in importance metrics.

## Summary

* Recommended versions: `XLRF_v1.2` > `XLRF_v2.1` > others.

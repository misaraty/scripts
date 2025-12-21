## **[中文版本](https://www.misaraty.com/2025-12-21_mgptips/)**

## mGPTIPS

A lightly modified implementation based on the symbolic regression software `GPTIPS`, released as `mGPTIPS`.

## Original GPTIPS

- [GPTIPS](https://sites.google.com/site/gptips4matlab/)

- Tutorial:

  [GPTIPS 2 – Tutorial 1 – “A Hard Problem for Standard GP”](https://sites.google.com/site/gptips4matlab/gptips-2-tutorial-examples/tutorial1)

## mGPTIPS

### Common Features

- Both implementations are based on a genetic programming symbolic regression framework, where analytical expressions are constructed from a function set and input variables, and optimized through population evolution, tournament selection, elite preservation, and explicit control of model complexity to obtain interpretable formulas.

### Improvements

- mGPTIPS introduces random Train / Validation / Test data splitting, applies z-score standardization driven by training-set statistics, incorporates joint constraints from validation performance and model complexity during model selection, and automatically saves predictions, evaluation metrics (R2, MAE, RMSE), model complexity, and parity plots, providing a more robust and reproducible workflow for practical data-driven studies.

- Based on the original `GPTIPS` distribution, only one existing `.m` file was modified and one additional `.m` file was added.

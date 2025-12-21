## **[中文版本](https://www.misaraty.com/2025-12-21_xlrf/)**

## RGNN

* This code implements a regression modeling approach based on a Residual Gated Neural Network, combined with K-fold cross-validation, exponential moving average, and ensemble learning to improve prediction stability and generalization performance.

## Features

* A deep neural network for data regression analysis is implemented, employing K-fold cross-validation to robustly evaluate model generalization performance. In each fold, the data are divided into training, validation, and test sets, and feature normalization is performed using statistics computed only from the training set to avoid data leakage. The model adopts a multilayer perceptron architecture with residual connections and gating mechanisms, which enhances feature representation through stacked nonlinear transformations while improving training stability via residual learning.

* In terms of training strategy, the framework integrates adaptive optimization, weight decay, gradient clipping, learning-rate scheduling, and early stopping to ensure numerical stability and mitigate overfitting. In addition, a small-scale ensemble is constructed within each fold, and an exponential moving average of network parameters is maintained, with ensemble-averaged predictions used as the final output. This design effectively reduces prediction variance and improves model robustness. The entire workflow automatically saves trained models and evaluation results and generates parity plots with consistent scales, forming a reproducible regression modeling framework.

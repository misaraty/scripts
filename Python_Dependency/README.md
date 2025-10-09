# Python Dependency

resolving dependency conflicts in packages

```bash
# Clear pip cache to avoid version conflicts
pip cache purge

# Basic scientific stack for pymatgen and spinney
pip install numpy==1.23.5
pip install pandas==2.1.4

# Crystal graph networks: MEGNet and M3GNet
pip install megnet
pip install m3gnet --no-deps

# AMSET: ab initio transport
pip install numba==0.58.1 spglib==2.2.0

# Graph neural networks: Spektral
pip install spektral --no-deps

# TensorFlow GNN (Graph Neural Network for TF)
pip install tensorflow-gnn --no-deps

# GPU support (CUDA 11.x)
pip install tensorflow-gpu==2.10.1
pip install pandas==1.5.3

# Atomic descriptor generation: DScribe
pip install dscribe==2.1.0

# Materials Project API
pip install mp-api==0.40.0

# JAX-based differentiable DFT (jax_dft)
pip install jax==0.4.23 jaxlib==0.4.23 numpy==1.23.5

# Time series deep learning (tsai)
pip install tsai fastcore==1.5.29 numpy==1.23.5
from tsai.all import *

# CGCNN
pip uninstall scikit-learn -y
pip install scikit-learn==1.4.0

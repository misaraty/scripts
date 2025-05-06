## Fewest Switch Surface Hopping — Tully Model 1

This repository contains two versions of a simple Fewest Switch Surface Hopping (FSSH) simulation for Tully Model 1 based on [fewest-switch-surface-hopping](https://github.com/Paul-St-Young/fewest-switch-surface-hopping).

## **[中文版本](https://www.misaraty.com/2025-05-06_%E6%9C%80%E5%B0%8F%E9%9D%A2%E8%B7%B3%E8%B7%83%E6%B3%95/)**

## 📁 version1

### Files:
- `model.nb`: Wolfram Mathematica notebook
- `fssh_v1.py`: Python script

### How to run:
1. Open and run `model.nb` in Mathematica to generate:
   - `Tully-model1-e1.dat`
   - `Tully-model1-e2.dat`
   - `Tully-model1-d12.dat`
   - `potential.pdf`

2. Then run the Python script:
   ```bash
   python fssh_v1.py
   ```
   It loads the `.dat` files and performs surface hopping simulation, generating:
   - `single.jpg`: potential and coupling plot
   - `model-1-prob.jpg`: hopping statistics

## 📁 version2

### File:
- `fssh_v2.py`: Fully integrated Python version

### How to run:
```bash
python fssh_v2.py
```
This script:
- Computes electronic surfaces and couplings directly
- Runs the FSSH simulation
- Generates the same output figures:
  - `single.jpg`
  - `model-1-prob.jpg`

## Reference
- [Tully J C. Molecular dynamics with electronic transitions[J]. The Journal of Chemical Physics, 1990, 93(2): 1061-1071.](https://xuv.scs.illinois.edu/chem540/GroupProjects/Benke_Hammer_SourcePaper1.pdf)
- [mudslide](https://github.com/smparker/mudslide)
- [FSSH](https://github.com/binggu56/FSSH)
- [Jain A, Sindhu A. Pedagogical overview of the fewest switches surface hopping method[J]. ACS omega, 2022, 7(50): 45810-45824.](https://pubs.acs.org/doi/10.1021/acsomega.2c04843)

## [中文版本](https://www.misaraty.com/2026-09-26_periodic_soed/)

## Periodic SOED

Periodic SOED is a periodic electron-density-based descriptor for representing inorganic crystals and predicting materials properties. It extends the smooth overlap of electron densities (SOED) concept from molecular systems to periodic crystals by sampling valence electron density around atomic centers under periodic boundary conditions and constructing rotationally invariant features using radial functions and spherical harmonics.

This repository provides the implementation of Periodic SOED for crystal band gap prediction, together with a periodic scaffold SOAP baseline, chemistry-matched comparisons, reduced-formula-grouped data splitting, XGBoost regression, statistical analysis, and figure regeneration. The benchmark is based on the MP-20-Charge dataset.

## Usage

Download the MP-20-Charge dataset from:

https://figshare.com/ndownloader/files/58973161

Extract the downloaded archive as a folder named `MP-20-Charge` in the same directory as `Periodic_SOED_v13.py`.

Then run:

```bash
python Periodic_SOED_v13.py
```

The script performs data preparation, reduced-formula-grouped splitting, periodic scaffold SOAP and Periodic SOED descriptor generation, XGBoost regression, validation-based model selection, statistical evaluation, robustness analysis, and export of the numerical results used in the manuscript.

## Citation

To be added after the paper is officially published.
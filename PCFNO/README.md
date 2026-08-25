## [中文版本](https://www.misaraty.com/2026-08-25_%E6%99%B6%E4%BD%93%E9%9A%8F%E6%9C%BA%E5%90%88%E9%87%91%E5%92%8C%E7%82%B9%E7%BC%BA%E9%99%B7%E6%A8%A1%E5%9E%8B%E7%94%9F%E6%88%90%E5%B7%A5%E5%85%B7/)

## PCFNO

PCFNO is a physically constrained Fourier neural operator for nonautoregressive prediction of normalized speed magnitude in cylinder wakes across Reynolds numbers. The model combines a frozen FNO backbone with a lightweight residual correction using geometry and boundary information.

<div align="center">
  <img src="./workflow.jpg" width="80%"/><br>
  Computational workflow of PCFNO
</div>

## Usage

Create a `cylinder` folder in the project directory.

Download `bc.zip`, `geo.zip`, and `prop.zip` from the [CFDBench cylinder dataset](https://huggingface.co/datasets/chen-yingfa/CFDBench/tree/main/cylinder).

Place the three ZIP files in `./cylinder/`. The script automatically extracts them when needed.

Run:

```bash
python PCFNO.py
```

The reported experiments use the cylinder PROP subset.

## Citation

To be added after the paper is officially published.

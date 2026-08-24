## [中文版本](https://www.misaraty.com/2026-08-21_paw%E8%B5%9D%E5%8A%BF%E8%87%AA%E5%8A%A8%E7%94%9F%E6%88%90%E5%B7%A5%E5%85%B7/)

# PAWGen

PAWGen is a lightweight Python tool for automatically generating VASP POTCAR files from POSCAR structures using the VASP 6.5.1 PAW potential library.

## Features

* Automatic POTCAR generation from POSCAR

* Supports VASP 6.5.1 PAW PBE and GW potentials

* Supports four potential selection modes:

  * `default`

    * Recommended PBE PAW potentials

  * `basic`

    * Basic PBE elemental potentials

  * `default_gw`

    * Recommended GW PAW potentials

  * `basic_gw`

    * Basic GW PAW potentials

* Automatically reports:

  * Selected potentials
  
  * ENMAX
  
  * Valence electrons
  
  * Total NELECT
  
  * Estimated VBM
  
  * Recommended ENCUT

## Installation

Download:

```
potpaw_PBE.64.tgz
```

from:

```
https://vasp.at/vasp-portal/
```

Extract:

```bash
tar -xzf potpaw_PBE.64.tgz
```

Place `PAWGen.py` into the extracted PAW potential directory:

```
potpaw6.5.1/
├── PAWGen.py
├── C/
├── Pb_d/
├── I/
└── ...
```

The directory should contain the corresponding PAW potential folders with `POTCAR` files.

## Configuration

Edit the following variables in `PAWGen.py`.

### Linux

```python
poscar_path = "./"
pot_path_root = "/home/user/potpaw6.5.1/"
```

### Windows

```python
poscar_path = "C:/Users/user/Desktop/"
pot_path_root = "C:/vasp/potpaw6.5.1/"
```

The folder specified by `poscar_path` should contain:

```
POSCAR
```

The generated file will be:

```
POTCAR
```

## Usage

### Default mode

Recommended PBE potentials:

```bash
python PAWGen.py default
```

or simply:

```bash
python PAWGen.py
```

### Basic mode

Basic PBE elemental potentials:

```bash
python PAWGen.py basic
```

### Default GW mode

Recommended GW potentials:

```bash
python PAWGen.py default_gw
```

### Basic GW mode

Basic GW potentials:

```bash
python PAWGen.py basic_gw
```

## Example Output

```bash
Mode: default

Element  Potential     ENMAX (eV)   Valence
-------------------------------------------------
C        C                400.000         4
N        N                400.000         5
H        H                250.000         1
Pb       Pb_d             237.835        14
I        I                175.647         7

Maximum ENMAX: 400 eV
Recommended ENCUT: 520 eV

NELECT: 432
VBM: 216
```

> [!NOTE]
> VASP software and PAW potentials must be obtained through official VASP channels.
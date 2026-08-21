## [中文版本]()

## PAWGen

PAWGen is a lightweight Python tool for automatically generating VASP
POTCAR files from POSCAR structures using the VASP 6.5.1 PAW potential
library.

## Features

-   Automatic POTCAR generation from POSCAR
-   Supports VASP 6.5.1 PAW PBE potentials
-   Supports three modes:
    -   default: recommended PBE PAW potentials
    -   basic: basic elemental potentials
    -   gw: GW potentials when available
-   Reports:
    -   selected potentials
    -   ENMAX
    -   valence electrons
    -   NELECT
    -   VBM
    -   recommended ENCUT

## Installation

Download:

    potpaw_PBE.64.tgz

from:

    https://vasp.at/vasp-portal/

Extract:

``` bash
tar -xzf potpaw_PBE.64.tgz
```

Place PAWGen.py in the extracted directory:

    potpaw6.5.1/
    ├── PAWGen.py
    ├── Cu/
    ├── In_d/
    └── S/

## Configuration

Linux:

``` python
poscar_path = "./"
pot_path_root = "/home/user/potpaw6.5.1/"
```

Windows:

``` python
poscar_path = "C:/Users/user/Desktop/"
pot_path_root = "C:/vasp/potpaw6.5.1/"
```

## Usage

Default:

``` bash
python PAWGen.py
```

Recommended PBE potentials:

``` bash
python PAWGen.py default
```

Basic potentials:

``` bash
python PAWGen.py basic
```

GW potentials:

``` bash
python PAWGen.py gw
```

## Example Output

``` bash
Mode: default

Element  Potential        ENMAX(eV)   Valence
C       C              400         4
N       N              400         5
H       H              250         1
Sn      Sn_d           241         14
I       I              176         7

Maximum ENMAX: 400 eV
Recommended ENCUT: 520 eV

NELECT: 864
VBM: 432
```

> VASP software and PAW potentials must be obtained through official VASP
channels.

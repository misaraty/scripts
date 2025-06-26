## **[中文版本](https://www.misaraty.com/2025-06-26_md%E9%94%AE%E9%95%BF%E9%94%AE%E8%A7%92%E5%88%86%E6%9E%90/)**

## MD Simulation

* **POSCAR**

Use [CaTiO₃ mp-5827](https://legacy.materialsproject.org/materials/mp-5827/) as an example. First perform a high-precision SCF calculation to generate the structure input for MD.

* **INCAR**

See the file: `INCAR` for recommended settings for NVE simulations at 300 K with van der Waals correction.

> [!NOTE]
> You can use `SMASS = -1` for pre-heating and then switch to `SMASS = -3` for equilibrium sampling.

## Generating pfiles

### Method 1: Using OUTCAR + CONTCAR

* Required files: `OUTCAR` and `CONTCAR` after MD simulation.

> [!NOTE]
> You can also use `OUTCAR + POSCAR` by slightly modifying the script.

* Use the script: `md_outcar_pfiles_v2.2.py`

This script reads coordinates from `OUTCAR` and creates individual structure files (`p0001`, `p0002`, ...) in a new folder `./pfiles`.

> [!NOTE]
> The generated `pXXXX` files are in **Cartesian coordinates**.

### Method 2: Using XDATCAR

* Required file: `XDATCAR` after MD simulation.

* Use the script: `md_xdatcar_pfiles_v1.1.py`

This script reads each MD frame and writes it to separate files in `./pfiles`.

> [!NOTE]
> The generated `pXXXX` files are in **fractional coordinates**.

## Bond Length Analysis

> [!NOTE]
> The script automatically detects and supports both Cartesian and Direct coordinates.

* Use the script: `bond_length_time_v7.2.py`  
This script calculates bond lengths between specified atom pairs over time.

> [!NOTE]
> - Set `num = 1000` for the number of MD steps.
> - Set `name = [[2, 4], [2, 5]]` to specify atom pairs (1-based indexing).
> - Works for any lattice system, including triclinic.

* Output file: `bond_length_time.dat`

* Example plot:

<div align="center">
  <img src="./bond_length.png" width="50%"/>
</div>

## Bond Angle Analysis

* Use the script: `angle_time_v7.2.py`  
This script computes time-dependent angles between triplets of atoms.

> [!NOTE]
> - Set `num = 1000` for the number of MD steps.
> - Set `name = [[3, 2, 5], [4, 2, 5]]` for analyzing angles ∠A–B–C (B is the vertex).
> - Atom indices are 1-based. Multiple angle triplets are supported.
> - Supports arbitrary crystal systems.

* Output file: `angle_time.dat`

* Example plot:

<div align="center">
  <img src="./bond_angle.png" width="50%"/>
</div>

## **[中文版本](https://www.misaraty.com/2025-05-21_epw%E8%AE%A1%E7%AE%97%E8%B6%85%E5%AF%BC%E8%BD%AC%E5%8F%98%E6%B8%A9%E5%BA%A6/)**

## Directory Structure

```shell
EPW Input Files Only/
├── job.sh                  # Job submission script
├── plot.py                 # Script for plotting Tc from EPW output
├── pseudo/                 # Pseudopotential directory
│   ├── B.pz-vbc.UPF
│   ├── Mg.pz-n-vbc.UPF
│   ├── Nb_ONCV_PBE-1.2.upf
│   └── pb_s.UPF
├── phonon/                 # SCF and phonon calculations
│   ├── scf.in              # Input file for SCF calculation
│   └── ph.in               # Input file for phonon calculation
└── epw/                    # NSCF and EPW calculations
    ├── nscf.in             # Input file for NSCF calculation
    ├── epw1.in             # First EPW run (α²F, λ, gap)
    └── epw2.in             # Second EPW run (Tc calculation)
```

## Tutorial

### Step 1: SCF Calculation (phonon/scf.in)
```bash
cd phonon
pw.x < scf.in > scf.out
```
Get self-consistent charge density.

### Step 2: Phonon Calculation (phonon/ph.in)
```bash
ph.x < ph.in > ph.out
```
Compute phonon dynamical matrices and electron-phonon interactions.

### Step 3: NSCF Calculation (epw/nscf.in)
```
cd ../epw
pw.x < nscf.in > nscf.out
```
Provide wavefunctions for Wannier interpolation.

### Step 4: First EPW Run (epw/epw1.in)
```bash
epw.x < epw1.in > epw1.out
```
Compute electron-phonon matrix elements, α²F(ω), λ, DOS, gap Δ(ω).

### Step 5: Second EPW Run (epw/epw2.in)
```bash
epw.x < epw2.in > epw2.out
```
Read α²F and solve linearized Eliashberg equation to get Tc.

### Step 6: Plotting (plot.py)
```bash
python plot.py
```
Plot of temperature vs. maximum eigenvalue; Tc where eigenvalue = 1.

<div align="center">
  <img src="./EPW%20Full%20Calculation/epw/tc.jpg" width="60%"/>
</div>

> [!NOTE]
> Modify: `/opt/ohpc/pub/apps/q-e-qe-7.3/EPW/bin/pp.py`，
> Replace the lines:
> ```python
> Enter the number of irr. q-points
> user_input = input(
>         'Enter the prefix used for PH calculations (e.g. diam)\n')
> prefix = str(user_input)
> ```
> with:
> ```python
> import sys
> prefix = sys.argv[1]
> ```

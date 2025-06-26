## **[中文版本](https://www.misaraty.com/2025-06-26_md%E9%94%AE%E9%95%BF%E9%94%AE%E8%A7%92%E5%88%86%E6%9E%90/)**

## **MD Simulation**

* **POSCAR**

Using [CaTiO₃ mp-5827](https://legacy.materialsproject.org/materials/mp-5827/) as an example, perform a high-precision SCF calculation to obtain a `POSCAR` file for `MD` input.

* **INCAR**

```shell
PREC = Accurate
SMASS = -3   #micro canonical ensemble
LREAL = A   #real space
IBRION = 0   #molecular dynamics
NBLOCK = 1   #update XDATCAR every x steps
TEBEG = 300   #start temperature
TEEND = 300   #final temperature
ISIF = 2   #2, ions change; 3, shape and ions change
ISYM = 0   #symmetry
NSW = 1000   #max ionic steps
POTIM = 1   #step width scaling
EDIFFG = -0.01   #ionic relaxation
EDIFF = 1e-5   #electronic SC-loop
LCHARG = .F.   #not save CHGCAR CHG
LWAVE = .F.   #not save WAVECAR
ISMEAR = 0   #gaussian smearing
SIGMA = 0.05   #the width of the smearing in eV
ALGO = Normal   #electronic minimisation algorithm
IVDW = 12   #DFT-D3 method with Becke-Jonson damping, van der Waals
```

> [!NOTE]
> If needed, use `SMASS = -1` for initial thermalization, then switch to `SMASS = -3` for temperature control.

## **Generating pfiles**

### **Method 1: Based on OUTCAR + CONTCAR**

Required: `OUTCAR` and `CONTCAR` from an `MD` run.

> [!NOTE]
> You may also use `OUTCAR + POSCAR` with slight modifications to the script.

Use the script `md_outcar_pfiles_v2.2.py` to generate pfiles.

```python
import os
import re
import math
import shutil

def print_intro():
    print(
        "Convert OUTCAR from MD simulation to multiple pfiles. Supports generating 0-9999 pfiles.\n\n"
        
        "Usage:\n"
        "Place OUTCAR and CONTCAR from MD simulation in the same directory.\n"
        "The script will generate multiple pfiles in this directory, matching the number of MD steps (NSW in INCAR).\n\n"
        
        'copyright by misaraty (misaraty@163.com) last update: 2025-06-26\n'
    )

def extract_nsw(outcar_path):
    # Extract NSW value from OUTCAR
    with open(outcar_path, 'r') as f:
        for line in f:
            if '   NSW' in line:
                parts = line.split()
                return int(parts[2])
    raise ValueError("NSW keyword not found in OUTCAR.")

def get_atom_count(contcar_path):
    # Calculate total number of atoms from CONTCAR
    with open(contcar_path, 'r') as f:
        lines = f.readlines()
    atom_counts = list(map(int, lines[6].split()))
    return sum(atom_counts)

def extract_forces(outcar_path, atom_count):
    # Extract atomic positions from all MD steps in OUTCAR
    forces = []
    with open(outcar_path, 'r') as f:
        lines = f.readlines()

    for idx, line in enumerate(lines):
        if 'POSITION                                       TOTAL-FORCE' in line:
            block_start = idx + 2
            block = lines[block_start:block_start + atom_count]
            for atom_line in block:
                coords = atom_line.split()[:3]
                forces.append("  " + "  ".join(coords))
    return forces

def write_pfiles(nsw, atom_count, contcar_path, forces, output_dir="pfiles"):
    # Write pfiles based on extracted coordinates
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    # Read first 7 lines of CONTCAR (lattice & atom info)
    with open(contcar_path, 'r') as f:
        contcar_lines = f.readlines()[:7]

    for step in range(nsw):
        filename = os.path.join(output_dir, f"p{step+1:04d}")
        with open(filename, 'w') as f:
            f.writelines(contcar_lines)
            f.write("C\n")  # atom label (can be adjusted if necessary)
            start = step * atom_count
            end = start + atom_count
            f.writelines(line + '\n' for line in forces[start:end])

def main():
    # Set working directory to script location and execute conversion
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    outcar = 'OUTCAR'
    contcar = 'CONTCAR'

    print_intro()
    nsw = extract_nsw(outcar)
    atom_count = get_atom_count(contcar)
    forces = extract_forces(outcar, atom_count)
    write_pfiles(nsw, atom_count, contcar, forces)
    print(f"{nsw} pfiles have been created.\n")

if __name__ == "__main__":
    main()
```

A `./pfiles` folder will be created with `p0001`, `p0002`, ..., etc.

> [!NOTE]
> This method generates `pXXXX` files in fractional coordinates.





















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

## **Bond Length Analysis**

> [!NOTE]
> The script automatically detects and supports both Cartesian and fractional coordinates.

* Use the script `bond_length_time_v7.2.py` to calculate time-dependent bond lengths.

```python
import os
import re
import shutil

def print_intro():
    print(
        "Convert XDATCAR to multiple pfiles\n\n"
        
        "Usage:\n"
        "Place the XDATCAR file in the current directory.\n"
        "This script will automatically create a 'pfiles' directory and output one file per step.\n\n"
        
        'copyright by misaraty (misaraty@163.com) last update: 2025-06-26\n'
    )

def read_header_lines(xdatcar_path):
    # Read the first 7 lines from XDATCAR: system name, scaling factor, 3x3 lattice, element symbols, and atom counts
    with open(xdatcar_path, 'r') as f:
        header = [next(f) for _ in range(7)]
    return header

def count_atoms(header_lines):
    # Count the total number of atoms from the 7th line of the header
    atom_counts = list(map(int, header_lines[6].split()))
    return sum(atom_counts)

def parse_xdatcar_blocks(xdatcar_path, atom_count):
    """
    Parse XDATCAR and return a list of snapshots.
    Each snapshot is a list of atomic coordinates for one MD step.
    """
    snapshots = []
    with open(xdatcar_path, 'r') as f:
        lines = f.readlines()

    config_indices = [i for i, line in enumerate(lines) if line.startswith('Direct configuration=')]

    for idx in range(len(config_indices)):
        start = config_indices[idx] + 1
        end = start + atom_count
        coords = lines[start:end]
        coords = ['  '.join(line.strip().split()) + '\n' for line in coords]
        snapshots.append(coords)

    return snapshots

def write_pfiles(header_lines, snapshots, output_dir='pfiles'):
    # Write each snapshot to a separate file in the specified output directory
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir)

    for i, coords in enumerate(snapshots, 1):
        fname = os.path.join(output_dir, f"p{i:04d}")
        with open(fname, 'w') as f:
            f.writelines(header_lines)
            f.write("Direct\n")
            f.writelines(coords)

def main():
    print_intro()
    xdatcar = 'XDATCAR'
    if not os.path.exists(xdatcar):
        raise FileNotFoundError("XDATCAR file not found in the current directory.")

    header = read_header_lines(xdatcar)
    atom_count = count_atoms(header)
    snapshots = parse_xdatcar_blocks(xdatcar, atom_count)
    write_pfiles(header, snapshots)
    print(f"Generated {len(snapshots)} pfiles in the './pfiles' directory.")

if __name__ == '__main__':
    main()
```

> [!NOTE]
> `num = 1000` should match the `NSW` value in `INCAR`.
> 
> `name = [[2, 4], [2, 5]]` defines atom pairs (index starts from 1).
> 
> Supports all lattice types, including non-orthogonal cells.

* Output file: `bond_length_time.dat`.























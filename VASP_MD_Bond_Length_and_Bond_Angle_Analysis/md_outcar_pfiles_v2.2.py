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

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

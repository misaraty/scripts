import os
os.chdir(os.path.split(os.path.realpath(__file__))[0])
print('copyright by Zhaosheng Zhang (misaraty@163.com)\n' + 'last update: 2025-07-22\n')

# ======== Parameter Settings ========
zmin = 0.0  # Minimum fractional z-coordinate to fix
zmax = 0.16  # Maximum fractional z-coordinate to fix
input_path = "POSCAR"
output_path = "POSCAR_sd"
# ====================================

def normalise_frac(z):
    """Normalize a fractional coordinate to the range [0, 1)."""
    while z < 0.0:
        z += 1.0
    while z >= 1.0:
        z -= 1.0
    return z

def process_poscar(lines, zmin, zmax):
    """
    Process POSCAR assuming:
    - No 'Selective Dynamics' present
    - Line 7 is 'Direct'
    - Coordinates start from line 8
    """
    header = lines[:8]
    natoms_list = list(map(int, header[6].split()))
    total_atoms = sum(natoms_list)

    coord_mode_idx = 7            # Line index of 'Direct'
    coord_start_idx = 8           # Coordinates start from line 8 (index 8)
    coord_end_idx = coord_start_idx + total_atoms
    coords = lines[coord_start_idx:coord_end_idx]

    processed = []
    count_fixed = 0

    for line in coords:
        parts = line.strip().split()
        if len(parts) < 3:
            continue
        try:
            z = normalise_frac(float(parts[2]))
        except ValueError:
            continue
        if zmin <= z <= zmax:
            tag = ["F", "F", "F"]
            count_fixed += 1
        else:
            tag = ["T", "T", "T"]
        # Preserve original formatting
        new_line = f"{'  '.join(parts[:3])}  {tag[0]} {tag[1]} {tag[2]}\n"
        processed.append(new_line)

    # Compose output lines
    new_lines = header[:7]
    new_lines.append("Selective Dynamics\n")
    new_lines.append(lines[coord_mode_idx])  # 'Direct'
    new_lines.extend(processed)
    if coord_end_idx < len(lines):
        new_lines.extend(lines[coord_end_idx:])  # Keep comments or blank lines

    # Print stats
    percent = 100.0 * count_fixed / total_atoms
    print(f"✅ Process completed.")
    print(f"Total atoms: {total_atoms}")
    print(f"Fixed atoms: {count_fixed}")
    print(f"Fixed fraction: {percent:.2f}%")

    return new_lines

# === Main execution ===
with open(input_path, 'r') as f:
    lines = f.readlines()

new_lines = process_poscar(lines, zmin, zmax)

with open(output_path, 'w') as f:
    f.writelines(new_lines)

print(f"Modified POSCAR saved to: {output_path}")

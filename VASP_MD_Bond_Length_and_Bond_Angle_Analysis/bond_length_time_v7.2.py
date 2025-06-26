import os
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# Set current working directory to the script location
os.chdir(os.path.dirname(os.path.abspath(__file__)))
print('copyright by misaraty (misaraty@163.com)\n' + 'last update: 2025-06-26\n')

# Parameters
num = 1000  # Number of files
name = [[2, 4], [2, 5]]  # Atom pairs (indices start from 1)

# Read lattice information and coordinate type from p0001
with open('./pfiles/p0001') as f:
    lines = f.readlines()

# Read 3×3 lattice matrix
lattice = np.array([
    list(map(float, lines[2].split())),
    list(map(float, lines[3].split())),
    list(map(float, lines[4].split()))
])

# Calculate inverse lattice matrix
inv_lattice = np.linalg.inv(lattice)

# Check if coordinates are in Direct mode
is_direct = lines[7].strip().startswith('Direct')

# Total number of atoms
natoms = sum(map(int, lines[6].split()))
coord_start = 8  # Starting line of coordinates

def get_cart_coords(filepath):
    # Read coordinates from pfile and convert to Cartesian if needed
    with open(filepath) as f:
        lines = f.readlines()
    coords = np.array([
        list(map(float, lines[coord_start + i].split()))
        for i in range(natoms)
    ])
    if is_direct:
        coords = coords @ lattice  # Convert from Direct to Cartesian
    return coords

def minimum_image(vec, lattice, inv_lattice):
    # Apply minimum image convention to a displacement vector
    frac = vec @ inv_lattice           # Convert to fractional
    frac -= np.round(frac)            # Apply periodic boundary condition
    return frac @ lattice             # Convert back to Cartesian

# Main loop to compute distances
list_sum = []

for pair in name:
    idx_a = pair[0] - 1
    idx_b = pair[1] - 1
    dist_list = []

    for i in range(1, num + 1):
        coords = get_cart_coords(f'./pfiles/p{i:04d}')
        delta = coords[idx_a] - coords[idx_b]
        delta = minimum_image(delta, lattice, inv_lattice)
        dist = np.linalg.norm(delta)
        dist_list.append(dist)

    list_sum.append(dist_list)

# Save results to file
np.savetxt('bond_length_time.dat', np.array(list_sum).T, header=str(name))

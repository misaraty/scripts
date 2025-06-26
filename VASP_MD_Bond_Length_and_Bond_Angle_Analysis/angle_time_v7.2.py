import os
import numpy as np
import warnings
warnings.filterwarnings('ignore')

# Set the path to the current script directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))
print('copyright by misaraty (misaraty@163.com)\n' + 'last update: 2025-06-26\n')

# Parameter settings
num = 1000  # Number of frames
name = [[3, 2, 5], [4, 2, 5]]  # Angles defined by triplets A-B-C, B is the vertex

# Read lattice info and coordinate type from p0001
with open('./pfiles/p0001') as f:
    lines = f.readlines()

# Lattice matrix and its inverse
lattice = np.array([
    list(map(float, lines[2].split())),
    list(map(float, lines[3].split())),
    list(map(float, lines[4].split()))
])
inv_lattice = np.linalg.inv(lattice)

# Check if coordinates are in Direct mode
is_direct = lines[7].strip().startswith('Direct')

# Total number of atoms
natoms = sum(map(int, lines[6].split()))
coord_start = 8  # Line where atomic coordinates start

# Function to read coordinates and convert to Cartesian if needed
def get_cart_coords(filepath):
    with open(filepath) as f:
        lines = f.readlines()
    coords = np.array([
        list(map(float, lines[coord_start + i].split()))
        for i in range(natoms)
    ])
    if is_direct:
        coords = coords @ lattice
    return coords

# Apply the minimum image convention
def minimum_image(vec, lattice, inv_lattice):
    frac = vec @ inv_lattice
    frac -= np.round(frac)
    return frac @ lattice

# Main loop to compute angles over all frames
list_sum = []

for triplet in name:
    idx_a, idx_b, idx_c = [x - 1 for x in triplet]
    angle_list = []

    for i in range(1, num + 1):
        coords = get_cart_coords(f'./pfiles/p{i:04d}')
        posA = coords[idx_a]
        posB = coords[idx_b]
        posC = coords[idx_c]

        vecBA = minimum_image(posA - posB, lattice, inv_lattice)
        vecBC = minimum_image(posC - posB, lattice, inv_lattice)

        # Calculate angle in degrees
        dot_product = np.dot(vecBA, vecBC)
        norm_BA = np.linalg.norm(vecBA)
        norm_BC = np.linalg.norm(vecBC)
        cos_theta = np.clip(dot_product / (norm_BA * norm_BC), -1.0, 1.0)
        angle_deg = np.degrees(np.arccos(cos_theta))

        angle_list.append(angle_deg)

    list_sum.append(angle_list)

# Save angle vs. time data
np.savetxt('angle_time.dat', np.array(list_sum).T, header=str(name))

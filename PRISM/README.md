## [中文版本]()

## PRISM: Pymatgen Random-alloy and point-defect Integrated Structure Maker

PRISM uses only `pymatgen` to generate reproducible supercells, SQS-like random alloys, vacancies, substitutions, and tetrahedral/octahedral interstitial structures from CIF or POSCAR files. 

PRISM is a general structure generator for periodic materials calculations. By editing a single configuration section at the beginning of the Python file, users can read an ordered CIF, POSCAR, or CONTCAR, construct a supercell, perform count/fraction/percentage substitutions, rank SQS-like alloy candidates, and build point-defect structures. Random-alloy ranking combines short-range pair statistics, local composition uniformity, and dopant separation. Multiple vacancies or substitutions can be randomly enumerated and ranked by their minimum periodic separation. Interstitial candidates can be identified using tetrahedral, octahedral, or largest-void geometry. Every structure is written as both CIF and POSCAR, while `manifest.csv` records provenance, composition, random seed, scores, site indices, coordinates, and defect distances.

## Quick start

1. Place the input structure next to `prism_structure_generator.py`.

2. Edit the `USER CONFIGURATION` section at the beginning of the script.

3. Run:

```bash
python prism_structure_generator.py
```

The input, output directory, and diagonal supercell multipliers can also be overridden temporarily:

```bash
python prism_structure_generator.py \
    -i input.cif \
    -o generated_structures \
    --supercell 2 2 2
```

To read a POSCAR:

```bash
python prism_structure_generator.py -i POSCAR
```

Show command-line help:

```bash
python prism_structure_generator.py --help
```

## Basic configuration

```python
INPUT_STRUCTURE = "input.cif"
OUTPUT_DIRECTORY = "generated_structures"
SUPERCELL_MATRIX = [2, 2, 2]
RANDOM_SEED = 2026
```

- `INPUT_STRUCTURE`: input CIF, POSCAR, or CONTCAR.

- `OUTPUT_DIRECTORY`: output directory.

- `SUPERCELL_MATRIX`: either `[a, b, c]` or a general 3 x 3 integer matrix.

- `RANDOM_SEED`: fixed seed for reproducibility.

All relative paths are resolved from the directory containing the script.

## SQS-like random alloys

### Percentage substitution

The following rule replaces 25% of the Te sites in the supercell with Se:

```python
GENERATE_RANDOM_ALLOY = True

RANDOM_REPLACEMENTS = [
    {"source": "Te", "target": "Se", "percent": 25},
]
```

Each rule must define exactly one of `count`, `fraction`, and `percent`:

```python
{"source": "Te", "target": "Se", "count": 8}
{"source": "Te", "target": "Se", "fraction": 0.25}
{"source": "Te", "target": "Se", "percent": 25}
```

Fractional counts use conventional half-up rounding:

```text
replacement count = int(parent-site count * fraction + 0.5)
```

For example, 25% of 30 Te sites is 7.5, which is rounded to 8 replacements.

### Ranking parameters

```python
SQS_TRIALS = 5000
SQS_KEEP = 3
SQS_PAIR_SHELLS = 6
SQS_SHELL_TOLERANCE = 0.05
SQS_LOCAL_NEIGHBORS = 8

SQS_PAIR_WEIGHT = 1.0
SQS_LOCAL_UNIFORMITY_WEIGHT = 2.0
SQS_SEPARATION_WEIGHT = 2.0
```

The combined score contains three components:

- `pair`: deviation of pair-shell statistics from an ideal random alloy.

- `local`: spatial fluctuations in local substitution concentration.

- `separation`: clustering penalty and periodic dopant dispersion.

For stronger emphasis on random-alloy statistics, use:

```python
SQS_PAIR_WEIGHT = 1.0
SQS_LOCAL_UNIFORMITY_WEIGHT = 0.5
SQS_SEPARATION_WEIGHT = 0.5
```

For stronger visual and spatial uniformity, use the default `1.0, 2.0, 2.0` or increase the final two weights.

> PRISM produces SQS-like structures, not rigorous multi-body SQS structures as defined by ATAT/mcsqs. If rigorous SQS quality is central to a publication, validate the selected structure with ATAT/mcsqs or icet.

## Defect parent

```python
DEFECT_BASE = "best_sqs"
```

- `best_sqs`: independently generate each defect from the best SQS-like structure.

- `supercell`: independently generate each defect from the unsubstituted supercell.

Each defect rule is applied independently to the same parent. PRISM does not automatically chain vacancy, substitution, and interstitial operations into one structure.

## Vacancies

```python
GENERATE_VACANCIES = True

VACANCY_DEFECTS = [
    {
        "name": "V_Te",
        "element": "Te",
        "count": 1,
        "selection": "random_separated",
        "structures": 3,
        "trials": 2000,
        "min_pair_distance": 0.0,
    },
]
```

Available `selection` modes:

- `random`: select random sites without comparing defect distances.

- `random_separated`: exhaust a small search space, or randomly enumerate a large search space, and rank configurations by decreasing minimum periodic separation. Recommended for multiple defects of the same type.

- `separated`: start near the cell center and use greedy farthest-point selection.

- `center`: select sites nearest the cell center. Useful only for visual placement of a single defect and not recommended for multiple defects.

A single vacancy has no vacancy-vacancy distance. In a low-symmetry structure or random alloy, generate several single-vacancy structures because nominally identical species can have different local chemical environments.

## Substitutions

The following example substitutes two Te sites with P and retains the three configurations with the largest periodic P-P separation:

```python
GENERATE_SUBSTITUTIONS = True

SUBSTITUTION_DEFECTS = [
    {
        "name": "P_Te",
        "source": "Te",
        "target": "P",
        "count": 2,
        "selection": "random_separated",
        "structures": 3,
        "trials": 5000,
        "min_pair_distance": 0.0,
    },
]
```

To require a P-P distance of at least 6 angstrom:

```python
"min_pair_distance": 6.0
```

The minimum periodic distance is included in the filename:

```text
03_substitution_P_Te_01_dmin11.301A.cif
POSCAR_03_substitution_P_Te_01_dmin11.301A
```

Configurations with the same `dmin` can still have different local alloy environments. For production calculations, relax all retained configurations and compare their total energies. If dopant clustering or binding is of interest, also construct short- and intermediate-distance configurations.

## Interstitials

An interstitial atom should not be placed at an arbitrary random coordinate. PRISM first identifies geometric candidates and then writes structures.

### Tetrahedral interstitial

```python
{
    "name": "Cd_i_tet",
    "element": "Cd",
    "site_type": "tetrahedral",
    "count": 1,
    "structures": 3,
    "grid": [12, 12, 12],
    "min_distance": 1.20,
    "tetra_distance_tolerance": 0.35,
    "tetra_angle_tolerance": 0.35,
    "fifth_neighbor_gap": 0.10,
    "dedup_distance": 0.70,
}
```

### Octahedral interstitial

```python
{
    "name": "Cd_i_oct",
    "element": "Cd",
    "site_type": "octahedral",
    "count": 1,
    "structures": 3,
    "grid": [12, 12, 12],
    "min_distance": 1.20,
    "octa_distance_tolerance": 0.35,
    "octa_angle_tolerance": 0.35,
    "seventh_neighbor_gap": 0.10,
    "dedup_distance": 0.70,
}
```

### Largest void

If the crystal does not contain a standard tetrahedral or octahedral void, use:

```python
"site_type": "largest_void"
```

Tetrahedral and octahedral labels are geometric candidates based on a grid search, near-neighbor distance similarity, and angular deviation from ideal polyhedra. They do not prove energetic stability. Relax every interstitial structure with first-principles calculations or a reliable potential and verify whether the atom migrates to another site.

## Outputs

Typical output tree:

```text
generated_structures/
|-- 00_supercell.cif
|-- POSCAR_00_supercell
|-- 01_sqs_rank01_score....cif
|-- POSCAR_01_sqs_rank01_score...
|-- 02_vacancy_....cif
|-- 03_substitution_....cif
|-- 04_interstitial_....cif
|-- interstitial_candidates_....csv
`-- manifest.csv
```

`manifest.csv` records:

- structure name and category;

- composition and atom count;

- random seed;

- combined and component SQS-like scores;

- defect site indices and fractional coordinates;

- periodic vacancy/substitution pair distances;

- interstitial type, coordinates, void radius, and neighbor geometry.

## Important notes

1. The input must be ordered. PRISM rejects a disordered CIF with partial occupancies.

2. Every distance uses periodic boundary conditions. Display at least a `2 x 2 x 2` replication when inspecting a structure in VESTA.

3. The cell center has no special physical meaning in a periodic bulk structure; `center` is only a visual convenience.

4. Maximizing defect separation is appropriate for dilute, approximately independent defects. Defect-complex studies must also include close configurations.

5. PRISM generates initial structures only. It does not calculate defect formation energies, charge states, finite-size corrections, or chemical potentials.

6. CIF and POSCAR outputs retain the same site order; POSCAR writing does not intentionally sort the structure.

7. Relax, inspect, and compare the energies of all candidate structures before using them in a publication.

## Citation

To be added after the paper is officially published.

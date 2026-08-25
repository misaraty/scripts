"""
PRISM: Pymatgen Random-alloy and point-defect Integrated Structure Maker.

Main capabilities:
1. Read a CIF, POSCAR, or CONTCAR and build a supercell.
2. Enumerate random occupations and select SQS-like configurations using
   pair-correlation, local-uniformity, and dopant-separation scores.
3. Independently generate vacancy, substitution, and interstitial structures
   from either the pristine supercell or the best SQS-like configuration.
4. Write every structure as both CIF and VASP POSCAR and create manifest.csv.

PRISM does not call ATAT/mcsqs. Its SQS-like mode ranks randomly generated
configurations using approximate short-range and spatial metrics. For work that
requires a rigorous multi-body SQS, validate the final structure with ATAT/mcsqs
or icet.

Quick start:
1. Edit the USER CONFIGURATION section below.
2. Place the input CIF or POSCAR in the same directory as this script.
3. Run: python prism_structure_generator.py

The input, output directory, and diagonal supercell can also be overridden:
    python prism_structure_generator.py -i other.cif \
        -o other_results --supercell 3 3 2
"""

import argparse
import csv
import json
import math
import os
import random
import re
import sys
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path


os.chdir(os.path.split(os.path.realpath(__file__))[0])

try:
    from pymatgen.core import Structure
    from pymatgen.io.cif import CifWriter
    from pymatgen.io.vasp import Poscar
except ImportError:
    raise SystemExit("pymatgen is required. Install it with: pip install pymatgen")


# =============================================================================
# USER CONFIGURATION: normally, only edit this section.
# =============================================================================

# Input structure: CIF, POSCAR, or CONTCAR.
INPUT_STRUCTURE = "input.cif"

# Output directory.
OUTPUT_DIRECTORY = "generated_structures"

# Supercell matrix. [2, 2, 2] doubles the cell along a, b, and c.
# A general 3 x 3 integer matrix is also accepted, e.g.
# [[2, 0, 0], [0, 2, 0], [0, 0, 2]].
SUPERCELL_MATRIX = [2, 2, 2]

# A fixed seed makes every random operation reproducible.
RANDOM_SEED = 2026


# ---------------------- SQS-like random substitution ----------------------

# Enable random substitution and SQS-like ranking.
GENERATE_RANDOM_ALLOY = True

# Exactly one of count, fraction, and percent must be given in each rule.
# The example replaces 25% of the Te sites in the supercell with Se.
# Half-up rounding is used: number = int(number_of_Te_sites * 25 / 100 + 0.5).
# Example: 30 Te sites * 25% = 7.5, so 8 Te sites are replaced.
RANDOM_REPLACEMENTS = [
    {"source": "Te", "target": "Se", "percent": 25},
]

# Number of random occupation trials. At least 5000 trials are recommended
# when replacing 25% of 32 parent sites. More trials expand the search but do
# not by themselves improve the discriminating power of the score.
SQS_TRIALS = 5000

# Number of lowest-scoring SQS-like configurations to save.
SQS_KEEP = 3

# Number of pair-distance shells included in the SQS-like score.
# Three shells often produce many ties; 5-8 shells are usually preferable.
SQS_PAIR_SHELLS = 6

# Pair distances within this tolerance (angstrom) are assigned to one shell.
# A 0.20-angstrom tolerance can merge distinct shells in distorted or
# low-symmetry structures; 0.05 angstrom is a practical starting value.
SQS_SHELL_TOLERANCE = 0.05

# Number of same-parent-sublattice neighbors used to evaluate local uniformity.
SQS_LOCAL_NEIGHBORS = 8

# Weights of the three SQS-like score components:
# pair: favors the short-range pair statistics of an ideal random alloy.
# local: suppresses large spatial fluctuations in substitution concentration.
# separation: suppresses dopant clustering and favors visually dispersed sites.
# For a more random-alloy-like structure, use 1.0, 0.5, 0.5.
# For stronger spatial uniformity, use the default 1.0, 2.0, 2.0.
SQS_PAIR_WEIGHT = 1.0
SQS_LOCAL_UNIFORMITY_WEIGHT = 2.0
SQS_SEPARATION_WEIGHT = 2.0


# ----------------------------- Defect parent -----------------------------

# "best_sqs": independently generate defects from the best SQS-like structure.
# "supercell": independently generate defects from the unsubstituted supercell.
DEFECT_BASE = "best_sqs"


# -------------------------------- Vacancy --------------------------------

GENERATE_VACANCIES = False

# element: species to remove; count: vacancies per generated structure.
# Available selection modes:
#   "random": random sites without comparing defect-defect distances.
#   "random_separated": random enumeration, or exhaustive enumeration for a
#       small search space, followed by maximum-distance ranking (recommended).
#   "separated": greedy farthest-point selection starting near the cell center.
#   "center": sites nearest the cell center; not recommended for many defects.
# structures: number of distinct vacancy structures to save.
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


# ------------------------------ Substitution ------------------------------

GENERATE_SUBSTITUTIONS = True

# Replace two Te atoms with P. Up to 5000 candidate combinations are examined,
# and the three configurations with the largest periodic P-P separation are saved.
SUBSTITUTION_DEFECTS = [
    {
        "name": "P_Te",
        "source": "Te",
        "target": "P",
        "count": 2,
        "selection": "random_separated",
        "structures": 3,
        "trials": 5000,
        # 0.0 disables the hard cutoff and ranks candidates only by distance.
        # Set this to 6.0 to require a minimum P-P distance of 6 angstrom.
        "min_pair_distance": 0.0,
    },
]


# ------------------------------ Interstitial ------------------------------

GENERATE_INTERSTITIALS = False

# site_type="tetrahedral": search for tetrahedral interstitial sites.
# site_type="octahedral": search for octahedral interstitial sites.
# site_type="largest_void": search for the largest void without fixing geometry.
# Both tetrahedral and octahedral rules below run when GENERATE_INTERSTITIALS=True.
INTERSTITIAL_DEFECTS = [
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
    },
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
    },
]

# =============================================================================
# END OF USER CONFIGURATION
# =============================================================================


def validate_ordered_structure(structure):
    """Require an ordered input because PRISM creates explicit occupations."""
    disordered = [
        index + 1 for index, site in enumerate(structure) if not site.is_ordered
    ]
    if disordered:
        preview = ", ".join(map(str, disordered[:10]))
        raise ValueError(
            f"The input contains partial or disordered sites, including site(s) {preview}. "
            "Create an ordered parent structure before running PRISM."
        )


def element_symbol(site):
    """Return the element symbol of an ordered site, including oxidized Species."""
    return site.specie.symbol


def composition_text(structure):
    """Return an element-count string suitable for the output manifest."""
    counts = Counter(element_symbol(site) for site in structure)
    return " ".join(f"{element}{counts[element]}" for element in sorted(counts))


def safe_name(text):
    """Convert user-supplied text into a safe filename component."""
    value = re.sub(r"[^A-Za-z0-9_.+-]+", "_", str(text).strip())
    return value.strip("_") or "structure"


def write_structure_pair(structure, output_dir, stem):
    """Write one structure as CIF and as a POSCAR without reordering sites."""
    stem = safe_name(stem)
    cif_path = output_dir / f"{stem}.cif"
    poscar_path = output_dir / f"POSCAR_{stem}"
    CifWriter(structure, symprec=None).write_file(str(cif_path))
    Poscar(structure, sort_structure=False, comment=stem).write_file(str(poscar_path))
    return cif_path, poscar_path


def normalize_supercell_matrix(value):
    """Accept either [a, b, c] or a general 3 x 3 integer scaling matrix."""
    if (
        isinstance(value, list)
        and len(value) == 3
        and all(isinstance(item, (int, float)) for item in value)
    ):
        matrix = [int(round(item)) for item in value]
        if any(item <= 0 for item in matrix):
            raise ValueError("Diagonal supercell multipliers must be positive integers")
        return matrix
    if (
        isinstance(value, list)
        and len(value) == 3
        and all(isinstance(row, list) and len(row) == 3 for row in value)
    ):
        return [[int(round(item)) for item in row] for row in value]
    raise ValueError("supercell must be [a, b, c] or a 3 x 3 integer matrix")


def number_from_rule(rule, pool_size):
    """Resolve count, fraction, or percent using half-up rounding."""
    has_count = "count" in rule
    has_fraction = "fraction" in rule
    has_percent = "percent" in rule
    if sum((has_count, has_fraction, has_percent)) != 1:
        raise ValueError(
            "Each replacement rule must define exactly one of count, fraction, or percent"
        )
    if has_count:
        number = int(rule["count"])
    else:
        fraction = (
            float(rule["fraction"]) if has_fraction else float(rule["percent"]) / 100.0
        )
        if not 0.0 <= fraction <= 1.0:
            raise ValueError("fraction must be 0-1 and percent must be 0-100")
        # Do not use Python round(), which applies bankers' rounding (2.5 -> 2).
        number = int(pool_size * fraction + 0.5)
    if number < 0:
        raise ValueError("The number of replacements cannot be negative")
    return number


def prepare_alloy_groups(supercell, replacement_rules):
    """Build fixed parent sublattices and prevent overlapping replacements."""
    initial_symbols = [element_symbol(site) for site in supercell]
    grouped_rules = defaultdict(list)
    for rule_index, raw_rule in enumerate(replacement_rules):
        rule = dict(raw_rule)
        rule["_rule_index"] = rule_index
        source = str(rule["source"])
        target = str(rule["target"])
        if source == target:
            raise ValueError(
                f"Replacement rule {rule_index + 1} has identical source and target"
            )
        grouped_rules[source].append(rule)

    groups = []
    for source, rules in grouped_rules.items():
        pool = [
            index for index, symbol in enumerate(initial_symbols) if symbol == source
        ]
        if not pool:
            raise ValueError(f"The supercell contains no parent species {source}")
        prepared = []
        total = 0
        for rule in rules:
            count = number_from_rule(rule, len(pool))
            total += count
            prepared.append({"target": str(rule["target"]), "count": count})
        if total > len(pool):
            raise ValueError(
                f"Parent species {source} has {len(pool)} sites, but {total} replacements were requested"
            )
        groups.append({"source": source, "pool": pool, "rules": prepared})
    return groups


def replacement_count_messages(supercell, replacement_rules):
    """Report percentage conversion and the final rounded replacement count."""
    initial_symbols = [element_symbol(site) for site in supercell]
    messages = []
    for rule in replacement_rules:
        source = str(rule["source"])
        target = str(rule["target"])
        pool_size = initial_symbols.count(source)
        if pool_size == 0:
            raise ValueError(f"The supercell contains no parent species {source}")
        number = number_from_rule(rule, pool_size)
        if "percent" in rule:
            description = f"{float(rule['percent']):g}%"
        elif "fraction" in rule:
            description = f"fraction={float(rule['fraction']):g}"
        else:
            description = f"count={int(rule['count'])}"
        messages.append(
            f"{source} -> {target}: {pool_size} parent sites, {description}, "
            f"rounded replacement count = {number}"
        )
    return messages


def build_pair_shells(structure, indices, shell_count, tolerance):
    """Group same-parent-sublattice pairs into periodic distance shells."""
    pairs = []
    for left_pos, left in enumerate(indices):
        for right in indices[left_pos + 1 :]:
            distance = float(structure.get_distance(left, right))
            if distance > 1.0e-8:
                pairs.append((distance, left, right))
    pairs.sort(key=lambda item: item[0])

    shells = []
    centers = []
    for distance, left, right in pairs:
        matched = None
        for shell_index, center in enumerate(centers):
            if abs(distance - center) <= tolerance:
                matched = shell_index
                break
        if matched is None:
            if len(shells) >= shell_count:
                continue
            centers.append(distance)
            shells.append([])
            matched = len(shells) - 1
        shells[matched].append((left, right))
    return shells, centers


def pair_correlation_score(structure, alloy_groups, prepared_shells):
    """Score pair-shell deviations; lower values better match a random alloy."""
    total_score = 0.0
    used_shells = 0
    for group, shells in zip(alloy_groups, prepared_shells):
        pool = group["pool"]
        labels = [element_symbol(structure[index]) for index in pool]
        concentrations = {
            key: value / len(pool) for key, value in Counter(labels).items()
        }
        species = sorted(concentrations)

        for shell in shells:
            if not shell:
                continue
            observed_counts = Counter()
            for left, right in shell:
                pair = tuple(
                    sorted(
                        (
                            element_symbol(structure[left]),
                            element_symbol(structure[right]),
                        )
                    )
                )
                observed_counts[pair] += 1
            pair_number = len(shell)
            shell_score = 0.0
            for first_pos, first in enumerate(species):
                for second in species[first_pos:]:
                    pair = (first, second)
                    if first == second:
                        expected_probability = concentrations[first] ** 2
                    else:
                        expected_probability = (
                            2.0 * concentrations[first] * concentrations[second]
                        )
                    observed_probability = observed_counts[pair] / pair_number
                    shell_score += (observed_probability - expected_probability) ** 2
            total_score += shell_score
            used_shells += 1
    return total_score / used_shells if used_shells else 0.0


def local_uniformity_score(structure, alloy_groups, neighbor_count):
    """Score local composition fluctuations; lower values are more uniform."""
    group_scores = []
    for group in alloy_groups:
        pool = group["pool"]
        if len(pool) <= 1:
            continue
        labels = {index: element_symbol(structure[index]) for index in pool}
        concentrations = {
            key: value / len(pool) for key, value in Counter(labels.values()).items()
        }
        species = sorted(concentrations)
        neighbors_used = min(max(1, int(neighbor_count)), len(pool) - 1)
        squared_errors = []

        for center in pool:
            nearest = sorted(
                (float(structure.get_distance(center, other)), other)
                for other in pool
                if other != center
            )[:neighbors_used]
            local_counts = Counter(labels[index] for _, index in nearest)
            for symbol in species:
                local_fraction = local_counts[symbol] / neighbors_used
                squared_errors.append((local_fraction - concentrations[symbol]) ** 2)
        if squared_errors:
            group_scores.append(sum(squared_errors) / len(squared_errors))
    return sum(group_scores) / len(group_scores) if group_scores else 0.0


def dopant_separation_score(structure, alloy_groups, shell_centers):
    """Score periodic dopant separation; lower values indicate less clustering."""
    penalties = []
    minimum_distances = []
    mean_nearest_distances = []

    for group_index, group in enumerate(alloy_groups):
        source = group["source"]
        pool = group["pool"]
        labels = {index: element_symbol(structure[index]) for index in pool}
        targets = sorted(set(labels.values()) - {source})
        centers = shell_centers[group_index]
        reference_distance = centers[0] if centers else 1.0

        for target in targets:
            target_indices = [index for index in pool if labels[index] == target]
            if len(target_indices) < 2:
                continue
            nearest_for_each = []
            all_pair_distances = []
            for left_pos, left in enumerate(target_indices):
                distances = [
                    float(structure.get_distance(left, right))
                    for right in target_indices
                    if right != left
                ]
                nearest_for_each.append(min(distances))
                for right in target_indices[left_pos + 1 :]:
                    all_pair_distances.append(
                        float(structure.get_distance(left, right))
                    )

            minimum_distance = min(all_pair_distances)
            mean_nearest = sum(nearest_for_each) / len(nearest_for_each)
            minimum_distances.append(minimum_distance)
            mean_nearest_distances.append(mean_nearest)

            # Consider both the closest pair and the mean nearest-dopant distance.
            penalties.append(
                0.5 * reference_distance / max(minimum_distance, 1.0e-12)
                + 0.5 * reference_distance / max(mean_nearest, 1.0e-12)
            )

    return {
        "score": sum(penalties) / len(penalties) if penalties else 0.0,
        "minimum_distance_A": min(minimum_distances) if minimum_distances else None,
        "mean_nearest_distance_A": (
            sum(mean_nearest_distances) / len(mean_nearest_distances)
            if mean_nearest_distances
            else None
        ),
    }


def alloy_score_components(
    structure, alloy_groups, prepared_shells, shell_centers, local_neighbors
):
    """Return pair, local-uniformity, and dopant-separation score components."""
    separation = dopant_separation_score(structure, alloy_groups, shell_centers)
    return {
        "pair": pair_correlation_score(structure, alloy_groups, prepared_shells),
        "local": local_uniformity_score(structure, alloy_groups, local_neighbors),
        "separation": separation["score"],
        "minimum_dopant_distance_A": separation["minimum_distance_A"],
        "mean_nearest_dopant_distance_A": separation["mean_nearest_distance_A"],
    }


def normalize_candidate_components(candidates, weights):
    """Normalize and weight components so numerical scale does not dominate."""
    component_names = ("pair", "local", "separation")
    extrema = {}
    for name in component_names:
        values = [candidate["components"][name] for candidate in candidates]
        extrema[name] = (min(values), max(values))

    for candidate in candidates:
        normalized = {}
        total = 0.0
        weight_sum = 0.0
        for name in component_names:
            low, high = extrema[name]
            # Each component is a nonnegative error or penalty. Scale it to
            # 0-1 by the candidate-set maximum without subtracting the minimum,
            # so normalization alone does not force the best score to zero.
            if high <= 1.0e-15:
                value = 0.0
            else:
                value = candidate["components"][name] / high
            normalized[name] = value
            weight = max(0.0, float(weights[name]))
            total += weight * value
            weight_sum += weight
        if weight_sum <= 0.0:
            raise ValueError("The three SQS-like score weights cannot all be zero")
        candidate["normalized_components"] = normalized
        candidate["score"] = total / weight_sum
    return extrema


def generate_sqs_like(supercell, sqs_config, seed):
    """Generate random candidates and retain the best combined-score structures."""
    rules = sqs_config.get("replacements", [])
    if not sqs_config.get("enabled", False) or not rules:
        return []

    groups = prepare_alloy_groups(supercell, rules)
    shell_count = max(1, int(sqs_config.get("pair_shells", 3)))
    tolerance = float(sqs_config.get("shell_tolerance", 0.20))
    prepared_shells = []
    shell_centers = []
    for group in groups:
        shells, centers = build_pair_shells(
            supercell, group["pool"], shell_count, tolerance
        )
        prepared_shells.append(shells)
        shell_centers.append(centers)

    trials = max(1, int(sqs_config.get("trials", 500)))
    keep = max(1, int(sqs_config.get("keep", 3)))
    local_neighbors = max(1, int(sqs_config.get("local_neighbors", 8)))
    weights = {
        "pair": float(sqs_config.get("pair_weight", 1.0)),
        "local": float(sqs_config.get("local_weight", 2.0)),
        "separation": float(sqs_config.get("separation_weight", 2.0)),
    }
    unique_candidates = {}

    for trial in range(trials):
        rng = random.Random(seed + trial * 104729)
        candidate = supercell.copy()
        signature_parts = []
        for group in groups:
            shuffled = list(group["pool"])
            rng.shuffle(shuffled)
            cursor = 0
            assigned = []
            for rule in group["rules"]:
                selected = sorted(shuffled[cursor : cursor + rule["count"]])
                cursor += rule["count"]
                for index in selected:
                    candidate.replace(index, rule["target"])
                assigned.append((rule["target"], tuple(selected)))
            signature_parts.append((group["source"], tuple(assigned)))
        signature = tuple(signature_parts)
        if signature in unique_candidates:
            continue
        components = alloy_score_components(
            candidate, groups, prepared_shells, shell_centers, local_neighbors
        )
        unique_candidates[signature] = {
            "trial": trial,
            "structure": candidate,
            "components": components,
        }

    candidates = list(unique_candidates.values())
    if not candidates:
        raise RuntimeError("No SQS-like candidate structure could be generated")
    component_ranges = normalize_candidate_components(candidates, weights)
    ranked = sorted(
        candidates,
        key=lambda item: (
            item["score"],
            item["components"]["pair"],
            item["components"]["local"],
            item["components"]["separation"],
            item["trial"],
        ),
    )[:keep]
    if not ranked:
        raise RuntimeError("No SQS-like candidate structure could be generated")
    return [
        {
            "structure": item["structure"],
            "score": item["score"],
            "trial": item["trial"],
            "components": item["components"],
            "normalized_components": item["normalized_components"],
            "component_ranges": component_ranges,
            "weights": weights,
            "shell_centers": shell_centers,
        }
        for item in ranked
    ]


def distance_to_cell_center(structure, index):
    """Return the shortest periodic distance from a site to the cell center."""
    distance, _ = structure.lattice.get_distance_and_image(
        structure[index].frac_coords, [0.5, 0.5, 0.5]
    )
    return float(distance)


def choose_defect_indices(structure, candidates, count, selection, rng, variant):
    """Choose defect sites randomly, near the center, or greedily separated."""
    if count <= 0:
        raise ValueError("Defect count must be a positive integer")
    if count > len(candidates):
        raise ValueError(
            f"Only {len(candidates)} candidate sites are available; cannot select {count}"
        )

    if selection == "random":
        return sorted(rng.sample(candidates, count))

    center_order = sorted(
        candidates, key=lambda index: distance_to_cell_center(structure, index)
    )
    if selection == "center":
        start = variant % len(center_order)
        rotated = center_order[start:] + center_order[:start]
        return sorted(rotated[:count])

    if selection == "separated":
        # Start near the cell center, then repeatedly add the farthest site.
        selected = [center_order[variant % len(center_order)]]
        remaining = set(candidates) - set(selected)
        while len(selected) < count:
            best = max(
                remaining,
                key=lambda index: min(
                    float(structure.get_distance(index, chosen)) for chosen in selected
                ),
            )
            selected.append(best)
            remaining.remove(best)
        return sorted(selected)

    raise ValueError("selection must be random, random_separated, center, or separated")


def defect_distance_statistics(structure, indices):
    """Return periodic vacancy/substitution pair distances in angstrom."""
    if len(indices) < 2:
        return {
            "minimum_pair_distance_A": None,
            "mean_pair_distance_A": None,
            "all_pair_distances_A": [],
        }
    distances = [
        float(structure.get_distance(left, right))
        for left, right in combinations(indices, 2)
    ]
    return {
        "minimum_pair_distance_A": min(distances),
        "mean_pair_distance_A": sum(distances) / len(distances),
        "all_pair_distances_A": sorted(distances),
    }


def ranked_random_defect_selections(structure, candidates, rule, seed, rule_index):
    """Enumerate candidate sets and rank them by maximum periodic separation."""
    count = int(rule.get("count", 1))
    number = max(1, int(rule.get("structures", 1)))
    trials = max(number, int(rule.get("trials", 2000)))
    minimum_required = max(0.0, float(rule.get("min_pair_distance", 0.0)))
    if count <= 0:
        raise ValueError("Defect count must be a positive integer")
    if count > len(candidates):
        raise ValueError(
            f"Only {len(candidates)} candidate sites are available; cannot select {count}"
        )

    rng = random.Random(seed + rule_index * 100003)

    # A single defect has no defect-defect distance; select distinct sites
    # reproducibly at random.
    if count == 1:
        shuffled = list(candidates)
        rng.shuffle(shuffled)
        return [[index] for index in shuffled[: min(number, len(shuffled))]]

    total_combinations = math.comb(len(candidates), count)
    sampled = set()

    # Exhaust a small combination space to find the true maximum separation;
    # otherwise sample random combinations.
    if total_combinations <= trials and total_combinations <= 200000:
        sampled.update(combinations(candidates, count))
    else:
        attempts = 0
        max_attempts = max(trials * 20, 1000)
        while len(sampled) < trials and attempts < max_attempts:
            sampled.add(tuple(sorted(rng.sample(candidates, count))))
            attempts += 1

    ranked = []
    rejected = 0
    for indices in sampled:
        stats = defect_distance_statistics(structure, indices)
        minimum_distance = stats["minimum_pair_distance_A"]
        if minimum_distance is not None and minimum_distance < minimum_required:
            rejected += 1
            continue
        ranked.append(
            (
                -float(minimum_distance),
                -float(stats["mean_pair_distance_A"]),
                rng.random(),
                tuple(indices),
            )
        )
    ranked.sort()

    if not ranked:
        raise RuntimeError(
            f"No defect set satisfies min_pair_distance={minimum_required:g} angstrom. "
            "Lower the cutoff, enlarge the supercell, or increase trials."
        )
    if len(ranked) < number:
        print(
            f"Warning: only {len(ranked)} sets satisfy the distance requirement, "
            f"fewer than the requested {number}; {rejected} sets were rejected.",
            file=sys.stderr,
        )
    return [list(item[3]) for item in ranked[:number]]


def unique_defect_selections(structure, candidates, rule, seed, rule_index):
    """Generate reproducible, preferably unique site sets for one defect rule."""
    count = int(rule.get("count", 1))
    number = max(1, int(rule.get("structures", 1)))
    selection = str(rule.get("selection", "random")).lower()
    if selection == "random_separated":
        return ranked_random_defect_selections(
            structure, candidates, rule, seed, rule_index
        )
    results = []
    seen = set()
    attempts = 0
    max_attempts = max(100, number * 30)
    while len(results) < number and attempts < max_attempts:
        rng = random.Random(seed + rule_index * 100003 + attempts * 9176)
        indices = choose_defect_indices(
            structure, candidates, count, selection, rng, attempts
        )
        signature = tuple(indices)
        if signature not in seen:
            seen.add(signature)
            results.append(indices)
        attempts += 1
    if len(results) < number:
        print(
            f"Warning: only {len(results)} unique configurations were generated, "
            f"fewer than the requested {number}.",
            file=sys.stderr,
        )
    return results


def vector_dot(first, second):
    return sum(float(a) * float(b) for a, b in zip(first, second))


def vector_norm(vector):
    return math.sqrt(max(0.0, vector_dot(vector, vector)))


def nearest_site_data(structure, frac_coords, number=5):
    """Return periodic distances, indices, and nearest-image vectors."""
    data = []
    for index, site in enumerate(structure):
        distance, image = structure.lattice.get_distance_and_image(
            frac_coords, site.frac_coords
        )
        vector_frac = [
            float(site.frac_coords[axis])
            + float(image[axis])
            - float(frac_coords[axis])
            for axis in range(3)
        ]
        vector_cart = structure.lattice.get_cartesian_coords(vector_frac)
        data.append((float(distance), index, vector_frac, vector_cart))
    data.sort(key=lambda item: item[0])
    return data[:number]


def tetrahedral_quality(structure, frac_coords, config):
    """Test whether four nearly equidistant neighbors form a tetrahedron."""
    nearest = nearest_site_data(structure, frac_coords, number=5)
    if len(nearest) < 5:
        return None
    first_four = nearest[:4]
    distances = [item[0] for item in first_four]
    mean_distance = sum(distances) / 4.0
    min_distance = float(config.get("min_distance", 1.0))
    distance_tolerance = float(config.get("tetra_distance_tolerance", 0.35))
    angle_tolerance = float(config.get("tetra_angle_tolerance", 0.35))
    fifth_gap_limit = float(config.get("fifth_neighbor_gap", 0.10))

    if min(distances) < min_distance:
        return None
    distance_spread = max(distances) - min(distances)
    if distance_spread > distance_tolerance:
        return None
    fifth_gap = nearest[4][0] - max(distances)
    if fifth_gap < fifth_gap_limit:
        return None

    # Vectors from a regular tetrahedron center to two vertices have cosine -1/3.
    angle_errors = []
    vectors = [item[3] for item in first_four]
    for left in range(4):
        for right in range(left + 1, 4):
            denominator = vector_norm(vectors[left]) * vector_norm(vectors[right])
            if denominator <= 1.0e-12:
                return None
            cosine = vector_dot(vectors[left], vectors[right]) / denominator
            angle_errors.append((cosine + 1.0 / 3.0) ** 2)
    angle_rms = math.sqrt(sum(angle_errors) / len(angle_errors))
    if angle_rms > angle_tolerance:
        return None

    # Lower is better; the final term weakly favors a larger void.
    score = (
        distance_spread / max(mean_distance, 1.0e-8) + angle_rms - 0.02 * mean_distance
    )
    return {
        "frac_coords": [float(value) % 1.0 for value in frac_coords],
        "score": float(score),
        "radius": float(min(distances)),
        "coordination_type": "tetrahedral",
        "neighbor_indices": [item[1] for item in first_four],
        "neighbor_distances": distances,
        "angle_rms": float(angle_rms),
    }


def refine_tetrahedral_point(structure, frac_coords):
    """Refine a grid point using the centroid of four nearest periodic images."""
    nearest = nearest_site_data(structure, frac_coords, number=4)
    if len(nearest) < 4:
        return [float(value) % 1.0 for value in frac_coords]
    unwrapped_vertices = []
    for _, _, vector_frac, _ in nearest:
        unwrapped_vertices.append(
            [float(frac_coords[axis]) + vector_frac[axis] for axis in range(3)]
        )
    return [
        sum(vertex[axis] for vertex in unwrapped_vertices) / 4.0 % 1.0
        for axis in range(3)
    ]


def octahedral_quality(structure, frac_coords, config):
    """Test whether six nearly equidistant neighbors form an octahedron."""
    nearest = nearest_site_data(structure, frac_coords, number=7)
    if len(nearest) < 7:
        return None
    first_six = nearest[:6]
    distances = [item[0] for item in first_six]
    mean_distance = sum(distances) / 6.0
    min_distance = float(config.get("min_distance", 1.0))
    distance_tolerance = float(config.get("octa_distance_tolerance", 0.35))
    angle_tolerance = float(config.get("octa_angle_tolerance", 0.35))
    seventh_gap_limit = float(config.get("seventh_neighbor_gap", 0.10))

    if min(distances) < min_distance:
        return None
    distance_spread = max(distances) - min(distances)
    if distance_spread > distance_tolerance:
        return None
    seventh_gap = nearest[6][0] - max(distances)
    if seventh_gap < seventh_gap_limit:
        return None

    # An ideal octahedron has 15 vertex-vector angles: three are 180 degrees
    # and twelve are 90 degrees. The three lowest sorted cosines should be
    # close to -1, and the remaining twelve should be close to zero.
    cosines = []
    vectors = [item[3] for item in first_six]
    for left in range(6):
        for right in range(left + 1, 6):
            denominator = vector_norm(vectors[left]) * vector_norm(vectors[right])
            if denominator <= 1.0e-12:
                return None
            cosine = vector_dot(vectors[left], vectors[right]) / denominator
            cosines.append(max(-1.0, min(1.0, cosine)))
    cosines.sort()
    errors = [(value + 1.0) ** 2 for value in cosines[:3]]
    errors.extend(value**2 for value in cosines[3:])
    angle_rms = math.sqrt(sum(errors) / len(errors))
    if angle_rms > angle_tolerance:
        return None

    score = (
        distance_spread / max(mean_distance, 1.0e-8) + angle_rms - 0.02 * mean_distance
    )
    return {
        "frac_coords": [float(value) % 1.0 for value in frac_coords],
        "score": float(score),
        "radius": float(min(distances)),
        "coordination_type": "octahedral",
        "neighbor_indices": [item[1] for item in first_six],
        "neighbor_distances": distances,
        "angle_rms": float(angle_rms),
    }


def refine_octahedral_point(structure, frac_coords):
    """Refine a grid point using the centroid of six nearest periodic images."""
    nearest = nearest_site_data(structure, frac_coords, number=6)
    if len(nearest) < 6:
        return [float(value) % 1.0 for value in frac_coords]
    unwrapped_vertices = []
    for _, _, vector_frac, _ in nearest:
        unwrapped_vertices.append(
            [float(frac_coords[axis]) + vector_frac[axis] for axis in range(3)]
        )
    return [
        sum(vertex[axis] for vertex in unwrapped_vertices) / 6.0 % 1.0
        for axis in range(3)
    ]


def periodic_point_distance(structure, first, second):
    distance, _ = structure.lattice.get_distance_and_image(first, second)
    return float(distance)


def find_interstitial_candidates(structure, config):
    """Search periodic tetrahedral, octahedral, or largest-void sites."""
    grid = config.get("grid", [12, 12, 12])
    if not isinstance(grid, list) or len(grid) != 3:
        raise ValueError("Interstitial grid must be [nx, ny, nz]")
    nx, ny, nz = [int(value) for value in grid]
    if min(nx, ny, nz) < 2:
        raise ValueError("Each interstitial grid dimension must be at least 2")
    site_type = str(config.get("site_type", "tetrahedral")).lower()
    if site_type not in {"tetrahedral", "octahedral", "largest_void"}:
        raise ValueError("site_type must be tetrahedral, octahedral, or largest_void")

    raw = []
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                point = [ix / nx, iy / ny, iz / nz]
                if site_type == "tetrahedral":
                    refined = refine_tetrahedral_point(structure, point)
                    quality = tetrahedral_quality(structure, refined, config)
                    if quality is not None:
                        raw.append(quality)
                elif site_type == "octahedral":
                    refined = refine_octahedral_point(structure, point)
                    quality = octahedral_quality(structure, refined, config)
                    if quality is not None:
                        raw.append(quality)
                else:
                    nearest = nearest_site_data(structure, point, number=1)
                    radius = nearest[0][0]
                    if radius >= float(config.get("min_distance", 1.0)):
                        raw.append(
                            {
                                "frac_coords": point,
                                "score": -radius,
                                "radius": radius,
                                "coordination_type": "largest_void",
                                "neighbor_indices": [nearest[0][1]],
                                "neighbor_distances": [radius],
                                "angle_rms": 0.0,
                            }
                        )

    raw.sort(key=lambda item: (item["score"], -item["radius"]))
    dedup_distance = float(config.get("dedup_distance", 0.70))
    unique = []
    for candidate in raw:
        if all(
            periodic_point_distance(
                structure, candidate["frac_coords"], accepted["frac_coords"]
            )
            >= dedup_distance
            for accepted in unique
        ):
            unique.append(candidate)
    if not unique:
        raise RuntimeError(
            f"No acceptable {site_type} interstitial site was found. Increase the grid "
            "and the corresponding distance/angle tolerances, reduce min_distance and "
            "the fifth/seventh-neighbor gap, or use site_type='largest_void'."
        )
    return unique


def select_interstitial_sets(structure, candidates, count, number):
    """Greedily select high-quality interstitial sites that remain separated."""
    if count <= 0:
        raise ValueError("Interstitial count must be a positive integer")
    if count > len(candidates):
        raise ValueError(
            f"Only {len(candidates)} interstitial candidates were found; cannot add {count} atoms"
        )
    results = []
    seen = set()
    for start in range(min(len(candidates), max(number * 5, number))):
        selected = [start]
        remaining = set(range(len(candidates))) - set(selected)
        while len(selected) < count:
            best = max(
                remaining,
                key=lambda index: (
                    min(
                        periodic_point_distance(
                            structure,
                            candidates[index]["frac_coords"],
                            candidates[chosen]["frac_coords"],
                        )
                        for chosen in selected
                    ),
                    -candidates[index]["score"],
                ),
            )
            selected.append(best)
            remaining.remove(best)
        signature = tuple(sorted(selected))
        if signature not in seen:
            seen.add(signature)
            results.append(selected)
        if len(results) >= number:
            break
    return results


def add_manifest_row(rows, stem, category, structure, seed, details, score=""):
    rows.append(
        {
            "stem": stem,
            "category": category,
            "composition": composition_text(structure),
            "atoms": len(structure),
            "score": score,
            "seed": seed,
            "details": details,
        }
    )


def generate_vacancies(base, rules, output_dir, seed, manifest_rows):
    """Apply each vacancy rule independently to the same parent structure."""
    for rule_index, rule in enumerate(rules):
        element = str(rule["element"])
        candidates = [
            index for index, site in enumerate(base) if element_symbol(site) == element
        ]
        if not candidates:
            raise ValueError(f"Vacancy species {element} is absent from the parent")
        selections = unique_defect_selections(base, candidates, rule, seed, rule_index)
        name = safe_name(rule.get("name", f"V_{element}"))
        for variant, indices in enumerate(selections, start=1):
            structure = base.copy()
            distance_stats = defect_distance_statistics(base, indices)
            coords = [
                list(map(float, structure[index].frac_coords)) for index in indices
            ]
            structure.remove_sites(sorted(indices, reverse=True))
            minimum_distance = distance_stats["minimum_pair_distance_A"]
            distance_tag = (
                "" if minimum_distance is None else f"_dmin{minimum_distance:.3f}A"
            )
            stem = f"02_vacancy_{name}_{variant:02d}{distance_tag}"
            write_structure_pair(structure, output_dir, stem)
            details = json.dumps(
                {
                    "selection": rule.get("selection", "random"),
                    "removed_indices_1based": [i + 1 for i in indices],
                    "frac_coords": coords,
                    **distance_stats,
                },
                ensure_ascii=False,
            )
            add_manifest_row(manifest_rows, stem, "vacancy", structure, seed, details)


def generate_substitutions(base, rules, output_dir, seed, manifest_rows):
    """Apply each substitution rule independently to the same parent structure."""
    for rule_index, rule in enumerate(rules):
        source = str(rule["source"])
        target = str(rule["target"])
        candidates = [
            index for index, site in enumerate(base) if element_symbol(site) == source
        ]
        if not candidates:
            raise ValueError(f"Substitution parent species {source} is absent")
        selections = unique_defect_selections(
            base, candidates, rule, seed, rule_index + 1000
        )
        name = safe_name(rule.get("name", f"{target}_{source}"))
        for variant, indices in enumerate(selections, start=1):
            structure = base.copy()
            distance_stats = defect_distance_statistics(base, indices)
            coords = [
                list(map(float, structure[index].frac_coords)) for index in indices
            ]
            for index in indices:
                structure.replace(index, target)
            minimum_distance = distance_stats["minimum_pair_distance_A"]
            distance_tag = (
                "" if minimum_distance is None else f"_dmin{minimum_distance:.3f}A"
            )
            stem = f"03_substitution_{name}_{variant:02d}{distance_tag}"
            write_structure_pair(structure, output_dir, stem)
            details = json.dumps(
                {
                    "source": source,
                    "target": target,
                    "selection": rule.get("selection", "random"),
                    "indices_1based": [i + 1 for i in indices],
                    "frac_coords": coords,
                    **distance_stats,
                },
                ensure_ascii=False,
            )
            add_manifest_row(
                manifest_rows, stem, "substitution", structure, seed, details
            )


def generate_interstitials(base, rules, output_dir, seed, manifest_rows):
    """Generate interstitials and save candidate coordinates and quality metrics."""
    for rule_index, rule in enumerate(rules):
        element = str(rule["element"])
        name = safe_name(rule.get("name", f"{element}_i"))
        candidates = find_interstitial_candidates(base, rule)
        candidate_path = output_dir / f"interstitial_candidates_{name}.csv"
        with candidate_path.open("w", newline="", encoding="utf-8-sig") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=[
                    "rank",
                    "coordination_type",
                    "frac_x",
                    "frac_y",
                    "frac_z",
                    "score",
                    "void_radius_A",
                    "angle_rms",
                    "neighbor_indices_1based",
                    "neighbor_distances_A",
                ],
            )
            writer.writeheader()
            for rank, candidate in enumerate(candidates, start=1):
                writer.writerow(
                    {
                        "rank": rank,
                        "coordination_type": candidate["coordination_type"],
                        "frac_x": candidate["frac_coords"][0],
                        "frac_y": candidate["frac_coords"][1],
                        "frac_z": candidate["frac_coords"][2],
                        "score": candidate["score"],
                        "void_radius_A": candidate["radius"],
                        "angle_rms": candidate["angle_rms"],
                        "neighbor_indices_1based": " ".join(
                            str(i + 1) for i in candidate["neighbor_indices"]
                        ),
                        "neighbor_distances_A": " ".join(
                            f"{value:.8f}" for value in candidate["neighbor_distances"]
                        ),
                    }
                )

        count = int(rule.get("count", 1))
        number = max(1, int(rule.get("structures", 1)))
        selections = select_interstitial_sets(base, candidates, count, number)
        for variant, selected in enumerate(selections, start=1):
            structure = base.copy()
            chosen = [candidates[index] for index in selected]
            for candidate in chosen:
                structure.append(
                    element,
                    candidate["frac_coords"],
                    coords_are_cartesian=False,
                    validate_proximity=False,
                )
            stem = f"04_interstitial_{name}_{variant:02d}"
            write_structure_pair(structure, output_dir, stem)
            details = json.dumps(
                {
                    "site_type": rule.get("site_type", "tetrahedral"),
                    "candidate_ranks": [index + 1 for index in selected],
                    "coordination_type": [item["coordination_type"] for item in chosen],
                    "frac_coords": [item["frac_coords"] for item in chosen],
                    "void_radius_A": [item["radius"] for item in chosen],
                },
                ensure_ascii=False,
            )
            add_manifest_row(
                manifest_rows, stem, "interstitial", structure, seed, details
            )


def write_manifest(output_dir, rows):
    """Save structure provenance, composition, scores, and defect sites."""
    path = output_dir / "manifest.csv"
    with path.open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "stem",
                "category",
                "composition",
                "atoms",
                "score",
                "seed",
                "details",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "PRISM: supercell, SQS-like alloy, vacancy, substitution, and "
            "interstitial structure generator"
        )
    )
    parser.add_argument(
        "-i",
        "--input",
        help=f"Override the configured input structure (current: {INPUT_STRUCTURE})",
    )
    parser.add_argument(
        "-o",
        "--output",
        help=f"Override the configured output directory (current: {OUTPUT_DIRECTORY})",
    )
    parser.add_argument(
        "--supercell",
        nargs=3,
        type=int,
        metavar=("A", "B", "C"),
        help="Override diagonal supercell multipliers, e.g. --supercell 2 2 2",
    )
    return parser.parse_args()


def main():
    args = parse_arguments()

    # Assemble the user-facing parameters into internal configuration mappings.
    matrix_parameter = (
        args.supercell if args.supercell is not None else SUPERCELL_MATRIX
    )
    sqs_config = {
        "enabled": GENERATE_RANDOM_ALLOY,
        "trials": SQS_TRIALS,
        "keep": SQS_KEEP,
        "pair_shells": SQS_PAIR_SHELLS,
        "shell_tolerance": SQS_SHELL_TOLERANCE,
        "local_neighbors": SQS_LOCAL_NEIGHBORS,
        "pair_weight": SQS_PAIR_WEIGHT,
        "local_weight": SQS_LOCAL_UNIFORMITY_WEIGHT,
        "separation_weight": SQS_SEPARATION_WEIGHT,
        "replacements": RANDOM_REPLACEMENTS,
    }
    defects_config = {
        "base": DEFECT_BASE,
        "vacancies": VACANCY_DEFECTS if GENERATE_VACANCIES else [],
        "substitutions": SUBSTITUTION_DEFECTS if GENERATE_SUBSTITUTIONS else [],
        "interstitials": INTERSTITIAL_DEFECTS if GENERATE_INTERSTITIALS else [],
    }

    seed = int(RANDOM_SEED)
    input_path = Path(args.input if args.input is not None else INPUT_STRUCTURE)
    if not input_path.is_file():
        raise FileNotFoundError(
            f"Input structure not found: {input_path}. Edit INPUT_STRUCTURE or "
            "provide a file with -i."
        )
    output_dir = Path(args.output if args.output is not None else OUTPUT_DIRECTORY)
    output_dir.mkdir(parents=True, exist_ok=True)

    structure = Structure.from_file(str(input_path))
    validate_ordered_structure(structure)
    supercell = structure.copy()
    matrix = normalize_supercell_matrix(matrix_parameter)
    supercell.make_supercell(matrix)
    manifest_rows = []

    # Always save the complete unsubstituted, defect-free supercell first.
    supercell_stem = "00_supercell"
    write_structure_pair(supercell, output_dir, supercell_stem)
    add_manifest_row(
        manifest_rows,
        supercell_stem,
        "supercell",
        supercell,
        seed,
        json.dumps({"matrix": matrix}, ensure_ascii=False),
    )

    # Rank SQS-like structures from best to worst; rank01 is the default defect parent.
    if GENERATE_RANDOM_ALLOY:
        for message in replacement_count_messages(supercell, RANDOM_REPLACEMENTS):
            print(f"Random substitution: {message}")
    sqs_results = generate_sqs_like(supercell, sqs_config, seed)
    for rank, result in enumerate(sqs_results, start=1):
        # Scientific notation prevents a small score from appearing as 0.00000000.
        stem = f"01_sqs_rank{rank:02d}_score{result['score']:.6e}"
        write_structure_pair(result["structure"], output_dir, stem)
        details = json.dumps(
            {
                "trial": result["trial"],
                "raw_score_components": result["components"],
                "normalized_score_components": result["normalized_components"],
                "component_ranges_in_all_trials": result["component_ranges"],
                "score_weights": result["weights"],
                "pair_shell_centers_A": result["shell_centers"],
            },
            ensure_ascii=False,
        )
        add_manifest_row(
            manifest_rows,
            stem,
            "sqs_like",
            result["structure"],
            seed,
            details,
            f"{result['score']:.12g}",
        )

    base_choice = str(defects_config.get("base", "best_sqs")).lower()
    if base_choice == "best_sqs" and sqs_results:
        defect_base = sqs_results[0]["structure"]
        base_description = "best_sqs"
    elif base_choice in {"best_sqs", "supercell"}:
        defect_base = supercell
        base_description = "supercell"
    else:
        raise ValueError("defects.base must be best_sqs or supercell")

    generate_vacancies(
        defect_base,
        defects_config.get("vacancies", []),
        output_dir,
        seed,
        manifest_rows,
    )
    generate_substitutions(
        defect_base,
        defects_config.get("substitutions", []),
        output_dir,
        seed,
        manifest_rows,
    )
    generate_interstitials(
        defect_base,
        defects_config.get("interstitials", []),
        output_dir,
        seed,
        manifest_rows,
    )
    write_manifest(output_dir, manifest_rows)

    print(f"Input structure: {input_path.resolve()}")
    print(f"Supercell matrix: {matrix}")
    print(f"Complete supercell: {composition_text(supercell)}, {len(supercell)} atoms")
    print(f"Number of SQS-like structures: {len(sqs_results)}")
    if sqs_results:
        best = sqs_results[0]
        components = best["components"]
        print(
            f"Best combined score: {best['score']:.12g} "
            "(normalized within the candidate set; lower is better)"
        )
        print(
            "Best raw components: "
            f"pair={components['pair']:.12g}, "
            f"local={components['local']:.12g}, "
            f"separation={components['separation']:.12g}"
        )
        print(
            "Best dopant distances: "
            f"minimum={components['minimum_dopant_distance_A']} angstrom, "
            f"mean_nearest={components['mean_nearest_dopant_distance_A']} angstrom"
        )
    print(f"Defect parent: {base_description}")
    print(f"Total number of structures: {len(manifest_rows)}")
    print(f"Output directory: {output_dir.resolve()}")


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"Error: {error}", file=sys.stderr)
        raise

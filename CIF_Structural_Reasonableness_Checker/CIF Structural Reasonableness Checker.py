import os
from pymatgen.io.cif import CifParser

os.chdir(os.path.split(os.path.realpath(__file__))[0])

CIF_FILE = "1550125.cif"

MIN_DIST_THRESH = 0.8
L_MIN = 1.0
L_MAX = 30.0
ANGLE_MIN = 20.0
ANGLE_MAX = 160.0
DENSITY_MIN = 1.0
DENSITY_MAX = 15.0
OCC_TOL = 0.05

COORD_CUTOFF = 3.0
MAX_CN = 12
CHARGE_TOL = 0.25


def check_min_distance(struct):
    try:
        dmat = struct.distance_matrix
        n = len(struct)
        min_dist = min(dmat[i][j] for i in range(n) for j in range(n) if i != j)
        if min_dist < MIN_DIST_THRESH:
            return False, f"min_dist {min_dist:.3f} < {MIN_DIST_THRESH}"
        return True, f"min_dist {min_dist:.3f}"
    except Exception as e:
        return False, f"distance_matrix error: {e}"


def check_lattice(struct):
    lat = struct.lattice
    a, b, c = lat.a, lat.b, lat.c
    alpha, beta, gamma = lat.alpha, lat.beta, lat.gamma
    msgs = []
    ok = True
    for name, v in [("a", a), ("b", b), ("c", c)]:
        if not (L_MIN <= v <= L_MAX):
            ok = False
            msgs.append(f"{name}={v:.3f} out of range [{L_MIN}, {L_MAX}]")
    for name, v in [("alpha", alpha), ("beta", beta), ("gamma", gamma)]:
        if not (ANGLE_MIN <= v <= ANGLE_MAX):
            ok = False
            msgs.append(f"{name}={v:.3f} out of range [{ANGLE_MIN}, {ANGLE_MAX}]")
    if ok:
        return True, f"a={a:.3f}, b={b:.3f}, c={c:.3f}, alpha={alpha:.2f}, beta={beta:.2f}, gamma={gamma:.2f}"
    return False, "; ".join(msgs)


def check_density(struct):
    try:
        rho = struct.density
        if not (DENSITY_MIN <= rho <= DENSITY_MAX):
            return False, f"density {rho:.3f} g/cm^3 out of range [{DENSITY_MIN}, {DENSITY_MAX}]"
        return True, f"density {rho:.3f} g/cm^3"
    except Exception as e:
        return False, f"density error: {e}"


def check_occupancies(struct):
    bad_sites = []
    for i, site in enumerate(struct.sites):
        occ = sum(site.species.values())
        if not (1.0 - OCC_TOL <= occ <= 1.0 + OCC_TOL):
            bad_sites.append((i, occ, site.species_string))
    if bad_sites:
        msg_parts = [f"site {idx} ({sp}): occ={occ:.3f}" for idx, occ, sp in bad_sites]
        return False, "occupancy not close to 1: " + "; ".join(msg_parts)
    return True, "all site occupancies close to 1.0"


def normalize_oxi_guess(g_raw):
    # convert keys to string for robust handling
    g = {}
    for k, v in g_raw.items():
        if hasattr(k, "symbol"):
            key = k.symbol
        else:
            key = str(k)
        g[key] = v
    return g


def check_oxidation_states(struct):
    comp = struct.composition
    try:
        guesses_raw = comp.oxi_state_guesses()
    except Exception as e:
        return False, f"oxi_state_guesses error: {e}"
    if not guesses_raw:
        return False, "no oxidation state guesses returned"
    g0 = normalize_oxi_guess(guesses_raw[0])
    msg = ", ".join(f"{el}:{val:.2f}" for el, val in g0.items())
    return True, "one oxidation-state guess: " + msg


def check_charge_neutrality(struct):
    comp = struct.composition
    try:
        guesses_raw = comp.oxi_state_guesses()
    except Exception as e:
        return False, f"oxi_state_guesses error: {e}"
    if not guesses_raw:
        return False, "no oxidation state guesses returned"

    guesses = [normalize_oxi_guess(g) for g in guesses_raw]
    el_amt = comp.get_el_amt_dict()

    best_charge = None
    best_guess = None
    for g in guesses:
        net = 0.0
        for el, amt in el_amt.items():
            if el in g:
                net += g[el] * amt
        if best_charge is None or abs(net) < abs(best_charge):
            best_charge = net
            best_guess = g

    if best_charge is None:
        return False, "could not compute net charge from guesses"

    if abs(best_charge) <= CHARGE_TOL:
        msg_guess = ", ".join(f"{el}:{val:.2f}" for el, val in best_guess.items())
        return True, f"net charge ~ {best_charge:.3f} (within ±{CHARGE_TOL}), guess: {msg_guess}"
    else:
        msg_guess = ", ".join(f"{el}:{val:.2f}" for el, val in best_guess.items())
        return False, f"net charge {best_charge:.3f} (>|{CHARGE_TOL}|), best guess: {msg_guess}"


def check_coordination(struct):
    per_site = []
    for i, site in enumerate(struct.sites):
        neighs = struct.get_neighbors(site, COORD_CUTOFF)
        cn = len(neighs)
        per_site.append((i, site.species_string, cn))

    bad_sites = [(i, sp, cn) for i, sp, cn in per_site if cn < 1 or cn > MAX_CN]
    if bad_sites:
        msg_parts = [f"site {i} ({sp}): CN={cn}" for i, sp, cn in bad_sites]
        return False, "abnormal coordination: " + "; ".join(msg_parts)

    by_elem = {}
    for _, sp, cn in per_site:
        by_elem.setdefault(sp, []).append(cn)
    summary = []
    for sp, cns in by_elem.items():
        summary.append(f"{sp}: min={min(cns)}, max={max(cns)}, avg={sum(cns)/len(cns):.2f}")
    return True, "coordination OK; " + "; ".join(summary)


def main():
    if not os.path.exists(CIF_FILE):
        print(f"File not found: {CIF_FILE}")
        return

    print(f"Checking CIF file: {CIF_FILE}")
    try:
        parser = CifParser(CIF_FILE)
        # use parse_structures to avoid FutureWarning
        structs = parser.parse_structures(primitive=True)
        if len(structs) == 0:
            print("PARSE FAIL: no structure parsed")
            return
        struct = structs[0]
    except Exception as e:
        print(f"PARSE FAIL: {e}")
        return

    results = []

    try:
        if struct.is_valid():
            results.append(("is_valid", True, "pymatgen is_valid() == True"))
        else:
            results.append(("is_valid", False, "pymatgen is_valid() == False"))
    except Exception as e:
        results.append(("is_valid", False, f"is_valid() error: {e}"))

    ok_min_dist, msg_min_dist = check_min_distance(struct)
    results.append(("min_distance", ok_min_dist, msg_min_dist))

    ok_lat, msg_lat = check_lattice(struct)
    results.append(("lattice", ok_lat, msg_lat))

    ok_rho, msg_rho = check_density(struct)
    results.append(("density", ok_rho, msg_rho))

    ok_occ, msg_occ = check_occupancies(struct)
    results.append(("occupancy", ok_occ, msg_occ))

    ok_oxi, msg_oxi = check_oxidation_states(struct)
    results.append(("oxidation_states", ok_oxi, msg_oxi))

    ok_charge, msg_charge = check_charge_neutrality(struct)
    results.append(("charge_neutrality", ok_charge, msg_charge))

    ok_coord, msg_coord = check_coordination(struct)
    results.append(("coordination", ok_coord, msg_coord))

    print("\n=== Check results ===")
    for name, ok, msg in results:
        status = "OK" if ok else "FAIL"
        print(f"{name:18s}: {status} - {msg}")

    all_ok = all(ok for _, ok, _ in results)
    print("\nOverall:", "REASONABLE" if all_ok else "HAS ISSUES")


if __name__ == "__main__":
    main()

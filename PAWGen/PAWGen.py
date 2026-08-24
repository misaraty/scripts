import math
import os
import sys

os.chdir(os.path.split(os.path.realpath(__file__))[0])

# https://vasp.at/wiki/Available_pseudopotentials

# POSCAR location
# windows
poscar_path = "C:/Users/misar/Desktop/"
pot_path_root = "./"
# linux
# poscar_path = "./"
# pot_path_root = "/home/misaraty/bin/potpaw6.5.1/"

# Available modes:
#   default    : recommended PBE potentials
#   basic      : simplest PBE potentials
#   default_gw : recommended GW potentials
#   basic_gw   : simplest GW potentials
potential_mode = "default"


def build_table(text):
    table = {}
    for line_number, line in enumerate(text.strip().splitlines(), 1):
        fields = line.split()
        if len(fields) != 4:
            raise ValueError(f"Invalid potential-table row {line_number}: {line}")
        symbol, potential, enmax, valence = fields
        if symbol in table:
            raise ValueError(f"Duplicate element in potential table: {symbol}")
        table[symbol] = (potential, float(enmax), int(valence))
    return table


PBE_DEFAULT = build_table(
    """
H H 250.000 1
He He 478.896 2
Li Li_sv 499.034 3
Be Be 247.543 2
B B 318.614 3
C C 400.000 4
N N 400.000 5
O O 400.000 6
F F 400.000 7
Ne Ne 343.606 8
Na Na_pv 259.561 7
Mg Mg 200.000 2
Al Al 240.300 3
Si Si 245.345 4
P P 255.040 5
S S 258.689 6
Cl Cl 262.472 7
Ar Ar 266.408 8
K K_sv 259.264 9
Ca Ca_sv 266.622 10
Sc Sc_sv 222.660 11
Ti Ti_sv 274.610 12
V V_sv 263.673 13
Cr Cr_pv 265.681 12
Mn Mn_pv 269.864 13
Fe Fe 267.882 8
Co Co 267.968 9
Ni Ni 269.532 10
Cu Cu 295.446 11
Zn Zn 276.723 12
Ga Ga_d 282.691 13
Ge Ge_d 310.294 14
As As 208.702 5
Se Se 211.555 6
Br Br 216.285 7
Kr Kr 185.331 8
Rb Rb_sv 220.112 9
Sr Sr_sv 229.353 10
Y Y_sv 202.626 11
Zr Zr_sv 229.898 12
Nb Nb_sv 293.235 13
Mo Mo_sv 242.676 14
Tc Tc_pv 263.523 13
Ru Ru_pv 240.049 14
Rh Rh_pv 247.408 15
Pd Pd 250.925 10
Ag Ag 249.844 11
Cd Cd 274.336 12
In In_d 239.211 13
Sn Sn_d 241.083 14
Sb Sb 172.069 5
Te Te 174.982 6
I I 175.647 7
Xe Xe 153.118 8
Cs Cs_sv 220.318 9
Ba Ba_sv 187.181 10
La La 219.292 11
Ce Ce 273.042 12
Pr Pr_3 181.719 11
Nd Nd_3 182.619 11
Pm Pm_3 176.959 11
Sm Sm_3 177.087 11
Eu Eu_2 99.328 8
Gd Gd_3 154.332 9
Tb Tb_3 155.613 9
Dy Dy_3 155.713 9
Ho Ho_3 154.137 9
Er Er_3 155.037 9
Tm Tm_3 149.221 9
Yb Yb_2 112.578 8
Lu Lu_3 154.992 9
Hf Hf_pv 220.334 10
Ta Ta_pv 223.667 11
W W_sv 223.057 14
Re Re 226.216 7
Os Os 228.022 8
Ir Ir 210.864 9
Pt Pt 230.283 10
Au Au 229.943 11
Hg Hg 233.204 12
Tl Tl_d 237.053 13
Pb Pb_d 237.835 14
Bi Bi_d 242.839 15
Po Po_d 264.565 16
At At 161.430 7
Rn Rn 151.497 8
Fr Fr_sv 214.540 9
Ra Ra_sv 237.367 10
Ac Ac 172.351 11
Th Th 247.306 12
Pa Pa 252.193 13
U U 252.502 14
Np Np 254.260 15
Pu Pu 254.353 16
Am Am 255.875 17
Cm Cm 257.953 18
"""
)

PBE_BASIC = build_table(
    """
H H 250.000 1
He He 478.896 2
Li Li 140.000 1
Be Be 247.543 2
B B 318.614 3
C C 400.000 4
N N 400.000 5
O O 400.000 6
F F 400.000 7
Ne Ne 343.606 8
Na Na 101.968 1
Mg Mg 200.000 2
Al Al 240.300 3
Si Si 245.345 4
P P 255.040 5
S S 258.689 6
Cl Cl 262.472 7
Ar Ar 266.408 8
K K_pv 116.731 7
Ca Ca_pv 119.559 8
Sc Sc 154.763 3
Ti Ti 178.330 4
V V 192.543 5
Cr Cr 227.080 6
Mn Mn 269.864 7
Fe Fe 267.882 8
Co Co 267.968 9
Ni Ni 269.532 10
Cu Cu 295.446 11
Zn Zn 276.723 12
Ga Ga 134.678 3
Ge Ge 173.807 4
As As 208.702 5
Se Se 211.555 6
Br Br 216.285 7
Kr Kr 185.331 8
Rb Rb_pv 121.882 7
Sr Sr_sv 229.353 10
Y Y_sv 202.626 11
Zr Zr_sv 229.898 12
Nb Nb_pv 208.608 11
Mo Mo 224.584 6
Tc Tc 228.694 7
Ru Ru 213.271 8
Rh Rh 228.996 9
Pd Pd 250.925 10
Ag Ag 249.844 11
Cd Cd 274.336 12
In In 95.934 3
Sn Sn 103.236 4
Sb Sb 172.069 5
Te Te 174.982 6
I I 175.647 7
Xe Xe 153.118 8
Cs Cs_sv 220.318 9
Ba Ba_sv 187.181 10
La La_s 136.530 9
Ce Ce_3 176.506 11
Pr Pr_3 181.719 11
Nd Nd_3 182.619 11
Pm Pm_3 176.959 11
Sm Sm_3 177.087 11
Eu Eu_2 99.328 8
Gd Gd_3 154.332 9
Tb Tb_3 155.613 9
Dy Dy_3 155.713 9
Ho Ho_3 154.137 9
Er Er_2 119.750 8
Tm Tm_3 149.221 9
Yb Yb_2 112.578 8
Lu Lu_3 154.992 9
Hf Hf 220.334 4
Ta Ta 223.667 5
W W 223.057 6
Re Re 226.216 7
Os Os 228.022 8
Ir Ir 210.864 9
Pt Pt 230.283 10
Au Au 229.943 11
Hg Hg 233.204 12
Tl Tl 90.140 3
Pb Pb 97.973 4
Bi Bi 105.037 5
Po Po 159.707 6
At At 161.430 7
Rn Rn 151.497 8
Fr Fr_sv 214.540 9
Ra Ra_sv 237.367 10
Ac Ac 172.351 11
Th Th_s 169.363 10
Pa Pa_s 193.466 11
U U 252.502 14
Np Np 254.260 15
Pu Pu 254.353 16
Am Am 255.875 17
Cm Cm 257.953 18
Cf Cf 414.614 20
"""
)

GW_DEFAULT = build_table(
    """
H H_GW 300.000 1
He He_GW 405.780 2
Li Li_sv_GW 433.699 3
Be Be_sv_GW 537.454 4
B B_GW 318.614 3
C C_GW 413.992 4
N N_GW 420.902 5
O O_GW 414.635 6
F F_GW 487.698 7
Ne Ne_GW 432.275 8
Na Na_sv_GW 372.853 9
Mg Mg_sv_GW 429.893 10
Al Al_GW 240.300 3
Si Si_GW 245.345 4
P P_GW 255.040 5
S S_GW 258.689 6
Cl Cl_GW 262.472 7
Ar Ar_GW 290.599 8
K K_sv_GW 248.998 9
Ca Ca_sv_GW 281.430 10
Sc Sc_sv_GW 378.961 11
Ti Ti_sv_GW 383.774 12
V V_sv_GW 382.321 13
Cr Cr_sv_GW 384.932 14
Mn Mn_sv_GW 384.627 15
Fe Fe_sv_GW 387.837 16
Co Co_sv_GW 387.491 17
Ni Ni_sv_GW 389.645 18
Cu Cu_sv_GW 467.331 19
Zn Zn_sv_GW 401.665 20
Ga Ga_d_GW 404.602 13
Ge Ge_d_GW 375.434 14
As As_GW 208.702 5
Se Se_GW 211.555 6
Br Br_GW 216.285 7
Kr Kr_GW 252.232 8
Rb Rb_sv_GW 221.197 9
Sr Sr_sv_GW 224.817 10
Y Y_sv_GW 339.758 11
Zr Zr_sv_GW 346.364 12
Nb Nb_sv_GW 353.872 13
Mo Mo_sv_GW 344.914 14
Tc Tc_sv_GW 351.044 15
Ru Ru_sv_GW 348.106 16
Rh Rh_sv_GW 351.206 17
Pd Pd_sv_GW 356.093 18
Ag Ag_sv_GW 354.430 19
Cd Cd_sv_GW 361.806 20
In In_d_GW 278.624 13
Sn Sn_d_GW 260.066 14
Sb Sb_d_GW 263.100 15
Te Te_GW 174.982 6
I I_GW 175.647 7
Xe Xe_GW 179.547 8
Cs Cs_sv_GW 198.101 9
Ba Ba_sv_GW 267.020 10
La La_GW 313.688 11
Ce Ce_GW 304.625 12
Hf Hf_sv_GW 309.037 12
Ta Ta_sv_GW 286.008 13
W W_sv_GW 317.132 14
Re Re_sv_GW 317.012 15
Os Os_sv_GW 319.773 16
Ir Ir_sv_GW 319.843 17
Pt Pt_sv_GW 323.669 18
Au Au_sv_GW 306.658 19
Hg Hg_sv_GW 312.028 20
Tl Tl_d_GW 237.053 15
Pb Pb_d_GW 237.809 16
Bi Bi_d_GW 261.876 17
Po Po_d_GW 267.847 18
At At_d_GW 266.251 17
Rn Rn_d_GW 267.347 18
"""
)

GW_BASIC = build_table(
    """
H H_GW 300.000 1
He He_GW 405.780 2
Li Li_GW 112.104 1
Be Be_GW 247.543 2
B B_GW 318.614 3
C C_GW 413.992 4
N N_GW 420.902 5
O O_GW 414.635 6
F F_GW 487.698 7
Ne Ne_GW 432.275 8
Na Na_sv_GW 372.853 9
Mg Mg_GW 126.143 2
Al Al_GW 240.300 3
Si Si_GW 245.345 4
P P_GW 255.040 5
S S_GW 258.689 6
Cl Cl_GW 262.472 7
Ar Ar_GW 290.599 8
K K_sv_GW 248.998 9
Ca Ca_sv_GW 281.430 10
Sc Sc_sv_GW 378.961 11
Ti Ti_sv_GW 383.774 12
V V_sv_GW 382.321 13
Cr Cr_sv_GW 384.932 14
Mn Mn_GW 278.466 7
Fe Fe_GW 321.007 8
Co Co_GW 323.400 9
Ni Ni_GW 357.323 10
Cu Cu_GW 417.039 11
Zn Zn_GW 328.191 12
Ga Ga_GW 134.678 3
Ge Ge_GW 173.807 4
As As_GW 208.702 5
Se Se_GW 211.555 6
Br Br_GW 216.285 7
Kr Kr_GW 252.232 8
Rb Rb_sv_GW 221.197 9
Sr Sr_sv_GW 224.817 10
Y Y_sv_GW 339.758 11
Zr Zr_sv_GW 346.364 12
Nb Nb_sv_GW 353.872 13
Mo Mo_sv_GW 344.914 14
Tc Tc_sv_GW 351.044 15
Ru Ru_sv_GW 348.106 16
Rh Rh_GW 247.408 9
Pd Pd_GW 250.925 10
Ag Ag_GW 249.844 11
Cd Cd_GW 254.045 12
In In_d_GW 278.624 13
Sn Sn_d_GW 260.066 14
Sb Sb_GW 172.069 5
Te Te_GW 174.982 6
I I_GW 175.647 7
Xe Xe_GW 179.547 8
Cs Cs_sv_GW 198.101 9
Ba Ba_sv_GW 267.020 10
La La_GW 313.688 11
Ce Ce_GW 304.625 12
Hf Hf_sv_GW 309.037 12
Ta Ta_sv_GW 286.008 13
W W_sv_GW 317.132 14
Re Re_sv_GW 317.012 15
Os Os_sv_GW 319.773 16
Ir Ir_sv_GW 319.843 17
Pt Pt_GW 248.716 10
Au Au_GW 248.344 11
Hg Hg_sv_GW 312.028 20
Tl Tl_d_GW 237.053 15
Pb Pb_d_GW 237.809 16
Bi Bi_GW 146.530 5
Po Po_d_GW 267.847 18
At At_d_GW 266.251 17
Rn Rn_d_GW 267.347 18
"""
)

POTENTIALS = {
    "default": PBE_DEFAULT,
    "basic": PBE_BASIC,
    "default_gw": GW_DEFAULT,
    "basic_gw": GW_BASIC,
}


def read_poscar(path):
    with open(path, "r", encoding="utf-8") as f:
        lines = f.readlines()

    if len(lines) < 7:
        raise ValueError(f"Invalid POSCAR: {path} has fewer than 7 lines")

    elements = lines[5].split()
    try:
        counts = [int(value) for value in lines[6].split()]
    except ValueError as exc:
        raise ValueError(
            "This script requires a VASP 5/6 POSCAR with element symbols on "
            "line 6 and atom counts on line 7"
        ) from exc

    if not elements or len(elements) != len(counts):
        raise ValueError(
            "The numbers of element symbols and atom counts in POSCAR do not match"
        )
    if any(count <= 0 for count in counts):
        raise ValueError("All atom counts in POSCAR must be positive integers")

    return elements, counts


def select_potential(symbol, mode):
    table = POTENTIALS[mode]
    if symbol not in table:
        if mode.endswith("_gw"):
            raise ValueError(
                f"No VASP 6.5.1 GW potential is listed for element {symbol} "
                f"in mode {mode}; PBE fallback is intentionally disabled"
            )
        raise ValueError(
            f"No VASP 6.5.1 potential information is listed for element {symbol} "
            f"in mode {mode}"
        )
    return table[symbol]


def resolve_potcars(elements, mode):
    selected = []
    missing = []

    for symbol in elements:
        potential, enmax, valence = select_potential(symbol, mode)
        potcar = os.path.join(pot_path_root, potential, "POTCAR")
        selected.append((symbol, potential, enmax, valence, potcar))
        if not os.path.isfile(potcar):
            missing.append(potcar)

    if missing:
        details = "\n".join(f"  {path}" for path in missing)
        raise FileNotFoundError(f"Missing POTCAR file(s):\n{details}")

    return selected


def write_potcar(path, selected):
    temporary_path = path + ".tmp"
    try:
        with open(temporary_path, "w", encoding="utf-8", newline="\n") as fout:
            for _, _, _, _, source_path in selected:
                with open(source_path, "r", encoding="utf-8", errors="replace") as fin:
                    content = fin.read()
                fout.write(content)
                if content and not content.endswith("\n"):
                    fout.write("\n")
        os.replace(temporary_path, path)
    except Exception:
        if os.path.exists(temporary_path):
            os.remove(temporary_path)
        raise


def main():
    mode = potential_mode
    if len(sys.argv) > 1:
        mode = sys.argv[1].lower()

    if mode not in POTENTIALS:
        available = ", ".join(POTENTIALS)
        raise ValueError(f"Unknown mode: {mode}. Available modes: {available}")

    poscar_file = os.path.join(poscar_path, "POSCAR")
    output_file = os.path.join(poscar_path, "POTCAR")
    elements, counts = read_poscar(poscar_file)
    selected = resolve_potcars(elements, mode)

    total_valence = 0
    max_enmax = 0.0

    print(f"Mode: {mode}\n")
    print(f"{'Element':<9}{'Potential':<12}{'ENMAX (eV)':>12}{'Valence':>10}")
    print("-" * 49)

    for (symbol, potential, enmax, valence, _), count in zip(selected, counts):
        print(f"{symbol:<9}{potential:<12}{enmax:>12.3f}{valence:>10d}")
        total_valence += valence * count
        max_enmax = max(max_enmax, enmax)

    write_potcar(output_file, selected)
    recommended_encut = int(math.ceil(max_enmax * 1.3 / 10.0) * 10)

    print(f"\nMaximum ENMAX: {max_enmax:.0f} eV")
    print(f"Recommended ENCUT: {recommended_encut} eV\n")

    print(f"NELECT: {total_valence}")
    print("VBM:", total_valence // 2)


if __name__ == "__main__":
    main()

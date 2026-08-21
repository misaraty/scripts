import os
import re
import sys

os.chdir(os.path.split(os.path.realpath(__file__))[0])

# POSCAR location
# windows
poscar_path = "C:/Users/misar/Desktop/"
pot_path_root = "./"
# linux
# poscar_path = "./"
# pot_path_root = "/home/misaraty/bin/potpaw6.5.1/"

# mode:
#   default : recommended PBE potentials
#   basic   : simplest elemental potentials
#   gw      : GW potentials when available
potential_mode = "default"
if len(sys.argv) > 1:
    input_mode = sys.argv[1].lower()
    if input_mode in ["default", "basic", "gw"]:
        potential_mode = input_mode
    else:
        raise ValueError("Unknown mode. Available modes: default, basic, gw")

print("Mode:", potential_mode, "\n")

# VASP 6.5.1 recommended PBE PAW potentials
# Source:
# https://www.vasp.at/wiki/index.php/Available_PAW_potentials

elements_data = [
    "H  H	250	1          ",
    "He  He	479	2          ",
    "Li  Li_sv	499	3      ",
    "Be  Be	248	2          ",
    "B  B	319	3          ",
    "C  C	400	4          ",
    "N  N	400	5          ",
    "O  O	400	6          ",
    "F  F	400	7          ",
    "Ne  Ne	344	8          ",
    "Na  Na_pv	260	7      ",
    "Mg  Mg	200	2          ",
    "Al  Al	240	3          ",
    "Si  Si	245	4          ",
    "P  P	255	5          ",
    "S  S	259	6          ",
    "Cl  Cl	262	7          ",
    "Ar  Ar	266	8          ",
    "K   K_sv	259	9      ",
    "Ca  Ca_sv	267	10     ",
    "Sc  Sc_sv	223	11     ",
    "Ti  Ti_sv	275	12     ",
    "V   V_sv	264	13     ",
    "Cr  Cr_pv	266	12     ",
    "Mn  Mn_pv	270	13     ",
    "Fe  Fe	268	8          ",
    "Co  Co	268	9          ",
    "Ni  Ni	270	10         ",
    "Cu  Cu	295	11         ",
    "Zn  Zn	277	12         ",
    "Ga  Ga_d	283	13     ",
    "Ge  Ge_d	310	14     ",
    "As  As	209	5          ",
    "Se  Se	212	6          ",
    "Br  Br	216	7          ",
    "Kr  Kr	185	8          ",
    "Rb  Rb_sv	220	9      ",
    "Sr  Sr_sv	229	10     ",
    "Y   Y_sv	203	11     ",
    "Zr  Zr_sv	230	12     ",
    "Nb  Nb_sv	293	13     ",
    "Mo  Mo_sv	243	14     ",
    "Tc  Tc_pv	264	13     ",
    "Ru  Ru_pv	240	14     ",
    "Rh  Rh_pv	247	15     ",
    "Pd  Pd	251	10         ",
    "Ag  Ag	250	11         ",
    "Cd  Cd	274	12         ",
    "In  In_d	239	13     ",
    "Sn  Sn_d	241	14     ",
    "Sb  Sb	172	5          ",
    "Te  Te	175	6          ",
    "I  I	176	7          ",
    "Xe  Xe	153	8          ",
    "Cs  Cs_sv	220	9      ",
    "Ba  Ba_sv	187	10     ",
    "La  La	219	11         ",
    "Ce  Ce	273	12         ",
    "Pr  Pr_3	182	11     ",
    "Nd  Nd_3	183	11     ",
    "Pm  Pm_3	177	11     ",
    "Sm  Sm_3	177	11     ",
    "Eu  Eu_2	99	8      ",
    "Gd  Gd_3	154	9      ",
    "Tb  Tb_3	156	9      ",
    "Dy  Dy_3	156	9      ",
    "Ho  Ho_3	154	9      ",
    "Er  Er_3	155	9      ",
    "Tm  Tm_3	149	9      ",
    "Yb  Yb_2	113	8      ",
    "Lu  Lu_3	155	9      ",
    "Hf  Hf_pv	220	10     ",
    "Ta  Ta_pv	224	11     ",
    "W   W_sv	223	14     ",
    "Re  Re	226	7          ",
    "Os  Os	228	8          ",
    "Ir  Ir	211	9          ",
    "Pt  Pt	230	10         ",
    "Au  Au	230	11         ",
    "Hg  Hg	233	12         ",
    "Tl  Tl_d	237	13     ",
    "Pb  Pb_d	238	14     ",
    "Bi  Bi_d	243	15     ",
    "Po  Po_d	265	16     ",
    "At  At	161 7     ",
    "Rn  Rn	152	8          ",
    "Fr  Fr_sv	215	9      ",
    "Ra  Ra_sv	237	10     ",
    "Ac  Ac	172	11         ",
    "Th  Th	247	12         ",
    "Pa  Pa	252	13         ",
    "U  U	253	14         ",
    "Np  Np	254	15         ",
    "Pu  Pu	254	16         ",
    "Am  Am	256	17         ",
    "Cm  Cm	258	18         ",
]

gw_appendix = {
    "H": "H_GW",
    "C": "C_GW",
    "N": "N_GW",
    "O": "O_GW",
    "Si": "Si_GW",
    "S": "S_GW",
    "Cu": "Cu_GW",
    "In": "In_GW",
    "Ga": "Ga_GW",
    "Se": "Se_GW",
    "I": "I_GW",
    "Zn": "Zn_GW",
}

basic_appendix = {
    "H": "H",
    "C": "C",
    "N": "N",
    "O": "O",
    "S": "S",
    "Cu": "Cu",
    "In": "In",
    "Ga": "Ga",
    "Se": "Se",
    "I": "I",
    "Zn": "Zn",
}


def read_poscar(path):
    with open(path, "r") as f:
        lines = f.readlines()

    elements = lines[5].split()
    counts = list(map(int, lines[6].split()))
    return elements, counts


def parse_default(symbol):
    for item in elements_data:
        data = item.split()
        if data[0] == symbol:
            return data[1], int(float(data[2])), int(data[3])
    raise ValueError(f"No potential information for {symbol}")


def select_potential(symbol):
    default, cutoff, valence = parse_default(symbol)

    if potential_mode == "gw":
        appendix = gw_appendix.get(symbol, default)
    elif potential_mode == "basic":
        appendix = basic_appendix.get(symbol, default)
    else:
        appendix = default

    return appendix, cutoff, valence


elements, counts = read_poscar(poscar_path + "POSCAR")

total_valence = 0
max_enmax = 0

print("Element  Potential        ENMAX(eV)   Valence")

with open(poscar_path + "POTCAR", "w") as fout:

    for symbol, count in zip(elements, counts):

        appendix, cutoff, valence = select_potential(symbol)

        pot_path = pot_path_root + appendix + "/POTCAR"

        if not os.path.exists(pot_path):
            raise FileNotFoundError(f"Missing POTCAR directory: {appendix}")

        print(f"{symbol:<8}{appendix:<15}{cutoff:<12}{valence}")

        with open(pot_path, "r") as fin:
            fout.write(fin.read())

        total_valence += valence * count
        if cutoff > max_enmax:
            max_enmax = cutoff

recommended_encut = int(((max_enmax * 1.3) + 9) // 10 * 10)

print("\nMaximum ENMAX:", max_enmax, "eV")
print("Recommended ENCUT:", recommended_encut, "eV")

print("\nNELECT:", total_valence)
print("VBM:", total_valence // 2)

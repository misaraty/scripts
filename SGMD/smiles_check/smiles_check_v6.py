import os, re, io, sys
import pandas as pd

try:
    os.chdir(os.path.split(os.path.realpath(__file__))[0])
except NameError:
    pass

from pathlib import Path
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolDescriptors, Descriptors, Lipinski, Crippen
from rdkit.Chem.Draw import rdMolDraw2D
from PIL import Image, ImageDraw, ImageFont

DM_CSV_PATH = "Dm.csv"
DM_SMILES_COL = "SMILES"
ENABLE_DM_FILTER = True

ENABLE_RING_FILTER = True
MIN_MAX_RING_SIZE = 12


LOG_FILE = open("smiles_check.log", "w", encoding="utf-8")


def log_print(msg):
    print(msg)
    LOG_FILE.write(msg + "\n")
    LOG_FILE.flush()


BASE_DIR = Path("runs/mol_chemprop_multi/outputs")
CSV_PATH = str(BASE_DIR / "generated_with_pred_MERGED_R5.csv")
SMILES_COL = "smiles"
D_COL = "D"

OUT_DIR = Path("./pic")
OUT_DIR.mkdir(parents=True, exist_ok=True)

IMG_SIZE = (800, 800)
TEXT_H = 220
CANVAS_SIZE = (IMG_SIZE[0], IMG_SIZE[1] + TEXT_H)
DPI = (300, 300)

PAD = 16
LINE_GAP = 38
TEXT_Y_OFFSET = -40
LABEL_FONT_SIZE = 34

MULTILINE_SPACING = 4
BG_MODE = "white"


def get_font(size=28):
    for p in [
        "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
        "C:/Windows/Fonts/arial.ttf",
        "/System/Library/Fonts/Supplemental/Arial.ttf",
    ]:
        if os.path.exists(p):
            return ImageFont.truetype(p, size=size)
    return ImageFont.load_default()


def normalize_smiles(smiles):
    """转换为 Canonical SMILES，用于跨文件一致匹配"""
    if not isinstance(smiles, str):
        return None
    smiles = smiles.strip()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    return Chem.MolToSmiles(mol, canonical=True)


def load_dm_smiles_set(csv_path=DM_CSV_PATH, smiles_col=DM_SMILES_COL):
    """读取 Dm.csv 并返回 canonical smiles 的 set；读不到就返回 None（表示跳过筛选）"""
    if not os.path.exists(csv_path):
        log_print(f"[warn] {csv_path} not found -> Dm filter skipped")
        return None
    try:
        df_dm = pd.read_csv(csv_path)
    except Exception as e:
        log_print(f"[warn] read {csv_path} failed -> Dm filter skipped ({e})")
        return None
    if smiles_col not in df_dm.columns:
        log_print(f"[warn] {csv_path} missing column {smiles_col} -> Dm filter skipped")
        return None

    dm_norm = df_dm[smiles_col].apply(normalize_smiles)
    dm_set = set(dm_norm.dropna())
    log_print(f"[info] Dm filter loaded: {len(dm_set)} unique canonical SMILES")
    return dm_set


FONT_LABEL = get_font(LABEL_FONT_SIZE)

# Reasonableness validation rules
CHECK = {
    # Structural validity
    "require_neutral": True,  # formal charge == 0
    "require_single_fragment": True,  # Salt / multi-component structures not allowed
    "require_no_radicals": False,  # Radicals not allowed
    "require_sanitize": True,  # Must pass RDKit SanitizeMol
    "require_chiral_sanity": False,  # Optional: stricter check (usually unnecessary)
    # Element constraints
    "enable_atom_filter": False,
    "allowed_atoms": {"C", "N", "O", "H"},
    # Size / complexity constraints
    "min_atoms": 2,
    "max_atoms": 80,  # With explicit hydrogens the count becomes large, heavy atoms are more stable
    "min_heavy_atoms": 3,
    "max_heavy_atoms": 60,
    # Typical physicochemical ranges
    "mw_min": 30.0,
    "mw_max": 1200.0,
    "logp_min": -2.0,
    "logp_max": 8.0,
    "tpsa_min": 0.0,
    "tpsa_max": 300.0,
    "hbd_max": 10,
    "hba_max": 20,
    "rotb_max": 20,
    "rings_max": 12,
    "arom_rings_max": 8,
    # Optional advanced filtering
    "enable_pains": True,  # Enable if FilterCatalog is available, otherwise skip automatically
    "enable_sa_score": True,  # Enable if sascorer can be imported, otherwise skip automatically
    "sa_max": 7.5,  # 1 (easy to synthesize) ~ 10 (hard to synthesize), typical threshold 6–8
}

SA_SCORER = None
try:
    here = (
        os.path.dirname(os.path.abspath(__file__))
        if "__file__" in globals()
        else os.getcwd()
    )
    sa_dir = os.path.join(here, "SA_Score")
    if os.path.isdir(sa_dir) and sa_dir not in sys.path:
        sys.path.insert(0, sa_dir)

    import sascorer as _sascorer

    SA_SCORER = _sascorer
    log_print("[info] SA scorer available: ENABLED")
except Exception as e:
    log_print(f"[warn] SA scorer not available -> skipped ({e})")
    SA_SCORER = None


def sanitize_filename(s, maxlen=80):
    s = re.sub(r"[\\/:*?\"<>|\s]+", "_", str(s)).strip("._")
    reserved = {
        "CON",
        "PRN",
        "AUX",
        "NUL",
        "COM1",
        "COM2",
        "COM3",
        "COM4",
        "COM5",
        "COM6",
        "COM7",
        "COM8",
        "COM9",
        "LPT1",
        "LPT2",
        "LPT3",
        "LPT4",
        "LPT5",
        "LPT6",
        "LPT7",
        "LPT8",
        "LPT9",
    }
    if s.upper() in reserved:
        s = f"_{s}_"
    return s[:maxlen]


def text_width(draw, s, font):
    box = draw.textbbox((0, 0), s, font=font)
    return box[2] - box[0]


def line_height(draw, font):
    box = draw.textbbox((0, 0), "Ag", font=font)
    return box[3] - box[1]


def wrap_text_by_pixel(draw, text, font, max_width_px):
    text = str(text)
    lines = []
    cur = ""
    for ch in text:
        trial = cur + ch
        if text_width(draw, trial, font) <= max_width_px or cur == "":
            cur = trial
        else:
            lines.append(cur)
            cur = ch
    if cur:
        lines.append(cur)
    return lines


def mol_to_image(mol, size, bg_mode="white"):
    drawer = rdMolDraw2D.MolDraw2DCairo(size[0], size[1])
    if bg_mode == "transparent":
        drawer.drawOptions().clearBackground = False
    else:
        drawer.drawOptions().clearBackground = True
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    png = drawer.GetDrawingText()
    img = Image.open(io.BytesIO(png)).convert("RGBA")
    if bg_mode == "transparent":
        return img
    bg = Image.new("RGBA", img.size, (255, 255, 255, 255))
    img = Image.alpha_composite(bg, img).convert("RGB")
    return img


PAINS_FILTER = None
try:
    from rdkit.Chem import FilterCatalog
    from rdkit.Chem.FilterCatalog import FilterCatalogParams

    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)
    PAINS_FILTER = FilterCatalog.FilterCatalog(params)
    log_print("[info] PAINS FilterCatalog available: ENABLED")
except Exception as e:
    log_print(f"[warn] PAINS FilterCatalog not available -> skipped ({e})")
    PAINS_FILTER = None


def count_radical_electrons(mol):
    try:
        if hasattr(Descriptors, "NumRadicalElectrons"):
            return int(Descriptors.NumRadicalElectrons(mol))
    except Exception:
        pass
    try:
        return int(sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms()))
    except Exception:
        return 0


def calc_max_ring_size(mol):
    """返回最大环元数；无环返回0。"""
    rings = Chem.GetSymmSSSR(mol)
    if not rings:
        return 0
    return max(len(r) for r in rings)


def calc_basic_props(mol):
    heavy = mol.GetNumHeavyAtoms()
    atoms = mol.GetNumAtoms()
    rings = rdMolDescriptors.CalcNumRings(mol)
    arom_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    rotb = Lipinski.NumRotatableBonds(mol)
    hbd = Lipinski.NumHDonors(mol)
    hba = Lipinski.NumHAcceptors(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    formal_charge = Chem.rdmolops.GetFormalCharge(mol)
    radicals = count_radical_electrons(mol)
    frags = len(Chem.GetMolFrags(mol))
    return {
        "atoms": atoms,
        "heavy_atoms": heavy,
        "rings": rings,
        "arom_rings": arom_rings,
        "rotb": rotb,
        "hbd": hbd,
        "hba": hba,
        "tpsa": tpsa,
        "mw": mw,
        "logp": logp,
        "formal_charge": formal_charge,
        "radicals": radicals,
        "frags": frags,
    }


def check_allowed_atoms(mol, allowed_atoms):
    bad = set()
    for a in mol.GetAtoms():
        sym = a.GetSymbol()
        if sym not in allowed_atoms:
            bad.add(sym)
    return bad


def check_pains(mol):
    if PAINS_FILTER is None:
        return None
    entry = PAINS_FILTER.GetFirstMatch(mol)
    if entry is None:
        return False
    return True


def calc_sa_score(mol):
    if SA_SCORER is None:
        return None
    try:
        return float(SA_SCORER.calculateScore(mol))
    except Exception:
        return None


def validate_mol(mol, can_smi):
    reasons = []

    # 1) Force sanitize (valence, aromaticity, kekulization checks, etc.)
    if CHECK["require_sanitize"]:
        try:
            Chem.SanitizeMol(mol)
        except Exception as e:
            reasons.append(f"sanitize_fail:{type(e).__name__}")
            return False, reasons, None

    if "." in can_smi:
        return False, ["dot_in_smiles"], None

    props = calc_basic_props(mol)

    # 2) Neutral / radical / multi-component checks
    if CHECK["require_neutral"] and props["formal_charge"] != 0:
        reasons.append(f"charge!=0({props['formal_charge']})")
    if CHECK["require_no_radicals"] and props["radicals"] != 0:
        reasons.append(f"radicals!=0({props['radicals']})")
    if CHECK["require_single_fragment"] and props["frags"] != 1:
        reasons.append(f"fragments!=1({props['frags']})")

    # 3) Element whitelist
    if CHECK.get("enable_atom_filter", False):
        allowed = CHECK.get("allowed_atoms", set())
        if allowed:
            bad_atoms = check_allowed_atoms(mol, allowed)
            if bad_atoms:
                reasons.append(f"bad_atoms:{','.join(sorted(bad_atoms))}")

    # 4) Size / complexity / physicochemical thresholds
    if props["atoms"] < CHECK["min_atoms"] or props["atoms"] > CHECK["max_atoms"]:
        reasons.append(f"atoms_out({props['atoms']})")
    if (
        props["heavy_atoms"] < CHECK["min_heavy_atoms"]
        or props["heavy_atoms"] > CHECK["max_heavy_atoms"]
    ):
        reasons.append(f"heavy_out({props['heavy_atoms']})")

    if props["mw"] < CHECK["mw_min"] or props["mw"] > CHECK["mw_max"]:
        reasons.append(f"mw_out({props['mw']:.1f})")
    if props["logp"] < CHECK["logp_min"] or props["logp"] > CHECK["logp_max"]:
        reasons.append(f"logp_out({props['logp']:.2f})")
    if props["tpsa"] < CHECK["tpsa_min"] or props["tpsa"] > CHECK["tpsa_max"]:
        reasons.append(f"tpsa_out({props['tpsa']:.1f})")

    if props["hbd"] > CHECK["hbd_max"]:
        reasons.append(f"hbd_out({props['hbd']})")
    if props["hba"] > CHECK["hba_max"]:
        reasons.append(f"hba_out({props['hba']})")
    if props["rotb"] > CHECK["rotb_max"]:
        reasons.append(f"rotb_out({props['rotb']})")

    if props["rings"] > CHECK["rings_max"]:
        reasons.append(f"rings_out({props['rings']})")
    if props["arom_rings"] > CHECK["arom_rings_max"]:
        reasons.append(f"arom_rings_out({props['arom_rings']})")

    # 5) PAINS filter
    if CHECK["enable_pains"]:
        pains = check_pains(mol)
        if pains is True:
            reasons.append("PAINS_hit")
        elif pains is None:
            props["pains"] = "skipped"
        else:
            props["pains"] = "pass"

    # 6) SA score
    if CHECK["enable_sa_score"]:
        sa = calc_sa_score(mol)
        if sa is None:
            props["sa_score"] = "skipped"
        else:
            props["sa_score"] = sa
            if sa > CHECK["sa_max"]:
                reasons.append(f"sa_out({sa:.2f})")

    ok = len(reasons) == 0
    return ok, reasons, props


def draw_card(mol, can_smiles, formula, mw, d_val):
    AllChem.Compute2DCoords(mol)

    mol_img = mol_to_image(mol, IMG_SIZE, bg_mode=BG_MODE)

    canvas = Image.new("RGB", CANVAS_SIZE, (255, 255, 255))
    if mol_img.mode == "RGBA":
        canvas.paste(mol_img.convert("RGB"), (0, 0), mol_img.split()[-1])
    else:
        canvas.paste(mol_img, (0, 0))

    draw = ImageDraw.Draw(canvas)
    y = IMG_SIZE[1] + TEXT_Y_OFFSET

    prefix = "SMILES: "
    prefix_w = text_width(draw, prefix, FONT_LABEL)
    max_w = CANVAS_SIZE[0] - PAD * 2 - prefix_w
    smiles_lines = wrap_text_by_pixel(draw, can_smiles, FONT_LABEL, max_w)

    draw.text((PAD, y), prefix + smiles_lines[0], fill=(0, 0, 0), font=FONT_LABEL)
    lh = line_height(draw, FONT_LABEL)
    for line in smiles_lines[1:]:
        y += lh + MULTILINE_SPACING
        draw.text((PAD + prefix_w, y), line, fill=(0, 0, 0), font=FONT_LABEL)

    y += lh + 8
    draw.text((PAD, y), f"Formula: {formula}", fill=(0, 0, 0), font=FONT_LABEL)
    y += LINE_GAP
    draw.text((PAD, y), f"MW: {mw:.2f} g/mol", fill=(0, 0, 0), font=FONT_LABEL)
    y += LINE_GAP
    draw.text((PAD, y), f"Log P_eff: {d_val:.2f}", fill=(0, 0, 0), font=FONT_LABEL)

    return canvas


def main():
    if not os.path.exists(CSV_PATH):
        raise FileNotFoundError(CSV_PATH)

    df = pd.read_csv(CSV_PATH)
    if SMILES_COL not in df.columns or D_COL not in df.columns:
        raise ValueError(f"missing columns: {SMILES_COL},{D_COL}")

    df[D_COL] = pd.to_numeric(df[D_COL], errors="coerce")
    df = df.dropna(subset=[SMILES_COL, D_COL]).sort_values(by=D_COL, ascending=False)

    df = df[df[D_COL] > 9.4]
    log_print(f"[info] rows before checks: {len(df)}")

    dm_set = None
    if ENABLE_DM_FILTER:
        dm_set = load_dm_smiles_set()

    n_total = 0
    n_pass = 0
    n_fail = 0

    for _, row in df.iterrows():
        n_total += 1
        smi = str(row[SMILES_COL]).strip()
        d_val = float(row[D_COL])

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            n_fail += 1
            log_print(f"[skip] invalid SMILES parse: {smi}")
            continue

        can_smi = Chem.MolToSmiles(mol, canonical=True)

        if dm_set is not None and can_smi in dm_set:
            n_fail += 1
            log_print(f"[reject] {can_smi} | reasons=already_in_Dm.csv")
            continue

        if ENABLE_RING_FILTER:
            max_ring = calc_max_ring_size(mol)
            if max_ring < MIN_MAX_RING_SIZE:
                n_fail += 1
                log_print(
                    f"[reject] {can_smi} | reasons=max_ring<{MIN_MAX_RING_SIZE}({max_ring})"
                )
                continue

        ok, reasons, props = validate_mol(mol, can_smi)
        if not ok:
            n_fail += 1
            log_print(f"[reject] {can_smi} | reasons={';'.join(reasons)}")
            continue
        n_pass += 1

        formula = rdMolDescriptors.CalcMolFormula(mol)
        mw = Descriptors.MolWt(mol)

        try:
            card = draw_card(mol, can_smi, formula, mw, d_val)
            fpath = OUT_DIR / f"{sanitize_filename(can_smi)}__D_{d_val:.2f}.png"
            card.save(str(fpath), dpi=DPI)

            if props is None:
                props = {}
            sa_v = props.get("sa_score", "na")
            pains_v = props.get("pains", "na")
            log_print(
                f"[saved] {fpath} | "
                f"MW={props.get('mw', float('nan')):.2f} LogP={props.get('logp', float('nan')):.2f} "
                f"TPSA={props.get('tpsa', float('nan')):.1f} HBD={props.get('hbd','?')} HBA={props.get('hba','?')} "
                f"RotB={props.get('rotb','?')} Rings={props.get('rings','?')} AromRings={props.get('arom_rings','?')} "
                f"SA={sa_v} PAINS={pains_v}"
            )
        except Exception as e:
            n_fail += 1
            log_print(f"[error] {can_smi} | {e}")

    log_print(f"[summary] total={n_total} pass={n_pass} fail={n_fail}")
    log_print(
        "[summary] NOTE: PAINS/SA may be skipped depending on your RDKit installation."
    )


if __name__ == "__main__":
    main()
    LOG_FILE.close()

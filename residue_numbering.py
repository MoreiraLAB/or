from pathlib import Path
import pandas as pd
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parent
PDB_DIR = ROOT / "structural_complexes"
SUMMARY = ROOT / "summary"
OFFSET_FILE = SUMMARY / "pdb_numbering_offsets.csv"


def ensure_offsets():
    SUMMARY.mkdir(exist_ok=True)

    rows = []
    parser = PDBParser(QUIET=True)

    for pdb in sorted(PDB_DIR.glob("*.pdb")):
        structure = parser.get_structure(pdb.stem, str(pdb))
        for chain in structure[0]:
            nums = []
            for residue in chain:
                het, resseq, icode = residue.id
                if het == " ":
                    nums.append(resseq)

            if nums:
                first = min(nums)
                last = max(nums)
                rows.append({
                    "pdb": pdb.name,
                    "complex": pdb.stem,
                    "chain": chain.id,
                    "first_pdb_residue": first,
                    "last_pdb_residue": last,
                    "n_residues": len(set(nums)),
                    "offset_if_metadata_starts_at_1": first - 1,
                })

    df = pd.DataFrame(rows)
    df.to_csv(OFFSET_FILE, index=False)
    return df


def load_offsets():
    if not OFFSET_FILE.exists():
        return ensure_offsets()
    return pd.read_csv(OFFSET_FILE)


def pdb_to_metadata_residue(pdb_name, chain, pdb_residue):
    """
    Convert PDB numbering to metadata numbering.

    Receptor chain A usually has offset 0.
    Partner chain B may have offsets, e.g. 5, 6, 7, 8, 9, 25.
    """
    df = load_offsets()
    row = df[(df["pdb"] == Path(str(pdb_name)).name) & (df["chain"] == str(chain))]

    if row.empty:
        return int(pdb_residue)

    offset = int(row.iloc[0]["offset_if_metadata_starts_at_1"])
    return int(pdb_residue) - offset


def metadata_to_pdb_residue(pdb_name, chain, metadata_residue):
    df = load_offsets()
    row = df[(df["pdb"] == Path(str(pdb_name)).name) & (df["chain"] == str(chain))]

    if row.empty:
        return int(metadata_residue)

    offset = int(row.iloc[0]["offset_if_metadata_starts_at_1"])
    return int(metadata_residue) + offset


def in_interval(x, start, end, end_exclusive=False):
    x = int(x)
    start = int(start)
    end = int(end)

    if end_exclusive:
        return start <= x < end
    return start <= x <= end

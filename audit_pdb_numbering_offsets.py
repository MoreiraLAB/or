from pathlib import Path
import pandas as pd
from Bio.PDB import PDBParser

ROOT = Path(__file__).resolve().parent
PDB_DIR = ROOT / "structural_complexes"
OUT = ROOT / "summary"
OUT.mkdir(exist_ok=True)

def first_last_residue(pdb_file):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure(pdb_file.stem, str(pdb_file))
    rows = []
    for chain in s[0]:
        nums = []
        for r in chain:
            het, resseq, icode = r.id
            if het == " ":
                nums.append(resseq)
        if nums:
            rows.append({
                "pdb": pdb_file.name,
                "complex": pdb_file.stem,
                "chain": chain.id,
                "first_pdb_residue": min(nums),
                "last_pdb_residue": max(nums),
                "n_residues": len(set(nums)),
                "offset_if_metadata_starts_at_1": min(nums) - 1
            })
    return rows

all_rows = []
for pdb in sorted(PDB_DIR.glob("*.pdb")):
    all_rows.extend(first_last_residue(pdb))

df = pd.DataFrame(all_rows)
df.to_csv(OUT / "pdb_numbering_offsets.csv", index=False)

print(df.to_string(index=False))
print()
print("Written:", OUT / "pdb_numbering_offsets.csv")

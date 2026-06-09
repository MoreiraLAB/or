from pathlib import Path
import shutil
import pandas as pd
import argparse

CSV = Path("metadata/template_comparisons_numbering_opioids.csv")

FIXES = {
    "6DDF": {
        "ICL1": "97-100",
        "TM2": "101-131",
        "ECL3": "307-309",
        "TM7": "310-339",
    },
    "6U1N": {
        "ECL1": "88-91",
        "TM3": "92-127",
    },
    "6PWC": {
        "ICL1": "92-95",
        "TM2": "96-130",
    },
}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    args = ap.parse_args()

    df = pd.read_csv(CSV, sep=";")

    print("Columns:", list(df.columns))
    print()

    for template_key, domain_map in FIXES.items():
        mask = df["template"].astype(str).str.upper().str.contains(template_key)

        if not mask.any():
            print("NOT FOUND template:", template_key)
            continue

        for domain, new_value in domain_map.items():
            if domain not in df.columns:
                print("NOT FOUND domain column:", domain)
                continue

            for idx in df.index[mask]:
                old_value = df.at[idx, domain]
                print(f"row {idx}: {template_key} {domain}: {old_value} -> {new_value}")

                if args.apply:
                    df.at[idx, domain] = new_value

    if args.apply:
        backup = CSV.with_suffix(".csv.bak_before_template_fix")
        shutil.copy2(CSV, backup)
        df.to_csv(CSV, sep=";", index=False)
        print()
        print("Applied fixes.")
        print("Backup:", backup)
        print("Updated:", CSV)
    else:
        print()
        print("Dry run only. To apply:")
        print("python fix_template_comparisons_numbering.py --apply")

if __name__ == "__main__":
    main()

from pathlib import Path
import re
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "images"
OUT.mkdir(exist_ok=True)

AA_ORDER = list("ACDEFGHIKLMNPQRSTVWY")

AA_COLORS = {
    "A": "#8dd3c7", "C": "#ffffb3", "D": "#fb8072", "E": "#fb8072",
    "F": "#bebada", "G": "#80b1d3", "H": "#b3de69", "I": "#fccde5",
    "K": "#bc80bd", "L": "#ccebc5", "M": "#ffed6f", "N": "#fdb462",
    "P": "#d9d9d9", "Q": "#fdb462", "R": "#bc80bd", "S": "#b3de69",
    "T": "#b3de69", "V": "#fccde5", "W": "#bebada", "Y": "#bebada",
}

FAMILY_ORDER = ["Arr 6PWC", "Arr 6U1N", "Gi/o", "Gs", "Gq/11", "G12/13", "Other"]

def pick_col(df, names):
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    for c in df.columns:
        if any(n.lower() in c.lower() for n in names):
            return c
    return None

def partner_family(partner):
    p = str(partner)
    if "Arr" in p or "ARR" in p:
        if "6PWC" in p:
            return "Arr 6PWC"
        if "6U1N" in p:
            return "Arr 6U1N"
        return "Arr"
    if re.match(r"Gi|Go|Gob|Gz", p):
        return "Gi/o"
    if re.match(r"Gs", p):
        return "Gs"
    if re.match(r"Gq|G11|G14|G15", p):
        return "Gq/11"
    if re.match(r"G12", p):
        return "G12/13"
    return "Other"

def find_table():
    for p in [
        ROOT / "processed_results" / "interface_contacts.csv",
        ROOT / "processed_results" / "interface_residues.csv",
        ROOT / "summary" / "interface_contacts.csv",
    ]:
        if p.exists():
            return p
    raise FileNotFoundError("No interface table found.")

def load_data():
    df = pd.read_csv(find_table())

    complex_col = pick_col(df, ["complex"])
    receptor_col = pick_col(df, ["receptor", "DXR"])
    partner_col = pick_col(df, ["partner"])

    rec_aa = pick_col(df, ["receptor_aa", "rec_aa", "aa_receptor", "receptor_residue_aa"])
    par_aa = pick_col(df, ["partner_aa", "par_aa", "aa_partner", "partner_residue_aa"])

    if complex_col is None:
        df["complex"] = df[receptor_col].astype(str) + "-" + df[partner_col].astype(str)
        complex_col = "complex"

    if receptor_col is None:
        df["receptor"] = df[complex_col].astype(str).str.split("-").str[0]
        receptor_col = "receptor"

    if partner_col is None:
        df["partner"] = df[complex_col].astype(str).str.split("-", n=1).str[1]
        partner_col = "partner"

    if rec_aa is None or par_aa is None:
        raise ValueError(f"Could not identify amino-acid columns. Columns: {list(df.columns)}")

    rows = []
    for side, aa_col in [("Receptor side", rec_aa), ("Partner side", par_aa)]:
        d = df[[complex_col, receptor_col, partner_col, aa_col]].copy()
        d.columns = ["complex", "receptor", "partner", "aa"]
        d["side"] = side
        d["aa"] = d["aa"].astype(str).str.extract(r"([A-Z])", expand=False)
        d = d[d["aa"].isin(AA_ORDER)]
        d["family"] = d["partner"].map(partner_family)
        rows.append(d)

    return pd.concat(rows, ignore_index=True)

def order_complexes(d):
    order = []
    for fam in FAMILY_ORDER:
        order.extend(sorted(d.loc[d["family"] == fam, "complex"].unique()))
    order.extend(sorted(set(d["complex"]) - set(order)))
    return order

def plot_side(df, side, outname):
    d = df[df["side"] == side].copy()

    counts = d.groupby(["family", "complex", "aa"]).size().reset_index(name="n")
    counts["percent"] = 100 * counts["n"] / counts.groupby("complex")["n"].transform("sum")

    order = order_complexes(counts)
    pivot = counts.pivot_table(index="complex", columns="aa", values="percent", fill_value=0)
    pivot = pivot.reindex(order)

    for aa in AA_ORDER:
        if aa not in pivot.columns:
            pivot[aa] = 0
    pivot = pivot[AA_ORDER]

    fams = counts.drop_duplicates("complex").set_index("complex")["family"].reindex(order)

    fig, ax = plt.subplots(figsize=(22, 8))
    bottom = [0] * len(pivot)

    x = range(len(pivot))
    for aa in AA_ORDER:
        vals = pivot[aa].values
        ax.bar(
            x,
            vals,
            bottom=bottom,
            color=AA_COLORS[aa],
            label=aa,
            width=0.88,
            edgecolor="white",
            linewidth=0.25,
        )
        bottom = [a + b for a, b in zip(bottom, vals)]

    # Separators and family headers
    last = None
    for i, fam in enumerate(fams):
        if i > 0 and fam != last:
            ax.axvline(i - 0.5, color="black", lw=1.0, alpha=0.45)
        last = fam

    for fam in FAMILY_ORDER:
        idx = [i for i, f in enumerate(fams) if f == fam]
        if idx:
            ax.text(
                (min(idx) + max(idx)) / 2,
                103,
                fam,
                ha="center",
                va="bottom",
                fontsize=11,
                fontweight="bold",
            )

    ax.set_ylim(0, 112)
    ax.set_ylabel("Interface residue percentage (%)", fontsize=12)
    ax.set_title(f"Figure S2 paper-like interfacial amino-acid percentage — {side}", fontsize=15, fontweight="bold")
    ax.set_xticks(list(x))
    ax.set_xticklabels(pivot.index, rotation=90, fontsize=7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

    ax.legend(
        title="Residue",
        bbox_to_anchor=(1.005, 1),
        loc="upper left",
        frameon=False,
        ncol=1,
        fontsize=8,
    )

    fig.tight_layout()
    fig.savefig(OUT / f"{outname}.png", dpi=300)
    fig.savefig(OUT / f"{outname}.svg")
    plt.close(fig)

    print(OUT / f"{outname}.png")
    print(OUT / f"{outname}.svg")

def main():
    df = load_data()
    plot_side(df, "Receptor side", "figure_s2_paper_like_receptor")
    plot_side(df, "Partner side", "figure_s2_paper_like_partner")

if __name__ == "__main__":
    main()

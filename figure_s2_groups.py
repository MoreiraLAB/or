from pathlib import Path
import re
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "images"
OUT.mkdir(exist_ok=True)

AA_GROUPS = {
    "Nonpolar aliphatic": list("AVLIMG"),
    "Aromatic": list("FYW"),
    "Polar uncharged": list("STNQCP"),
    "Basic positive": list("KRH"),
    "Acidic negative": list("DE"),
}

AA_GROUP_COLORS = {
    "Nonpolar aliphatic": "#F4A261",
    "Aromatic": "#E76F51",
    "Polar uncharged": "#2A9D8F",
    "Basic positive": "#457B9D",
    "Acidic negative": "#9B5DE5",
}

GROUP_ORDER = [
    "Gi/o",
    "Gs",
    "Gq/11",
    "G12/13",
    "Arrestin 6PWC",
    "Arrestin 6U1N",
    "Other",
]

def aa_to_group(aa):
    aa = str(aa).strip().upper()
    for group, aas in AA_GROUPS.items():
        if aa in aas:
            return group
    return "Other"

def partner_family(partner):
    p = str(partner)
    if re.search(r"Arr|ARR", p):
        if "6PWC" in p:
            return "Arrestin 6PWC"
        if "6U1N" in p:
            return "Arrestin 6U1N"
        return "Arrestin"
    if re.match(r"Gi|Go|Gob|Gz", p):
        return "Gi/o"
    if re.match(r"Gs", p):
        return "Gs"
    if re.match(r"Gq|G11|G14|G15", p):
        return "Gq/11"
    if re.match(r"G12", p):
        return "G12/13"
    return "Other"

def find_contacts():
    candidates = [
        ROOT / "processed_results" / "interface_contacts.csv",
        ROOT / "processed_results" / "interface_residues.csv",
        ROOT / "summary" / "interface_contacts.csv",
    ]
    for p in candidates:
        if p.exists():
            return p
    raise FileNotFoundError("No interface contact/residue table found in processed_results/ or summary/.")

def pick_col(df, names):
    lower = {c.lower(): c for c in df.columns}
    for n in names:
        if n.lower() in lower:
            return lower[n.lower()]
    for c in df.columns:
        cl = c.lower()
        if any(n.lower() in cl for n in names):
            return c
    return None

def load_long_table():
    path = find_contacts()
    df = pd.read_csv(path)

    complex_col = pick_col(df, ["complex", "Complex"])
    receptor_col = pick_col(df, ["receptor", "DXR"])
    partner_col = pick_col(df, ["partner", "Partner"])

    rec_aa_col = pick_col(df, ["receptor_aa", "rec_aa", "aa_receptor", "residue_receptor", "receptor_residue_aa"])
    par_aa_col = pick_col(df, ["partner_aa", "par_aa", "aa_partner", "residue_partner", "partner_residue_aa"])

    if complex_col is None:
        if receptor_col and partner_col:
            df["complex"] = df[receptor_col].astype(str) + "-" + df[partner_col].astype(str)
            complex_col = "complex"
        else:
            raise ValueError("Could not identify complex/receptor/partner columns.")

    if receptor_col is None:
        df["receptor"] = df[complex_col].astype(str).str.split("-").str[0]
        receptor_col = "receptor"

    if partner_col is None:
        df["partner"] = df[complex_col].astype(str).str.split("-", n=1).str[1]
        partner_col = "partner"

    if rec_aa_col is None or par_aa_col is None:
        raise ValueError(
            "Could not identify receptor/partner amino-acid columns. "
            f"Columns available: {list(df.columns)}"
        )

    rows = []
    for side, aa_col in [("Receptor side", rec_aa_col), ("Partner side", par_aa_col)]:
        tmp = df[[complex_col, receptor_col, partner_col, aa_col]].copy()
        tmp.columns = ["complex", "receptor", "partner", "aa"]
        tmp["side"] = side
        tmp["aa"] = tmp["aa"].astype(str).str.extract(r"([A-Z])", expand=False)
        tmp = tmp.dropna(subset=["aa"])
        tmp["aa_group"] = tmp["aa"].map(aa_to_group)
        tmp["partner_family"] = tmp["partner"].map(partner_family)
        rows.append(tmp)

    return pd.concat(rows, ignore_index=True)

def make_percentages(df, side):
    d = df[df["side"] == side].copy()
    counts = (
        d.groupby(["partner_family", "complex", "aa_group"])
        .size()
        .reset_index(name="n")
    )
    totals = counts.groupby(["partner_family", "complex"])["n"].transform("sum")
    counts["percent"] = 100 * counts["n"] / totals
    return counts

def ordered_complexes(d):
    order = []
    for fam in GROUP_ORDER:
        vals = sorted(d.loc[d["partner_family"] == fam, "complex"].unique())
        order.extend(vals)
    remaining = sorted(set(d["complex"]) - set(order))
    return order + remaining

def plot_side(df, side, filename):
    d = make_percentages(df, side)
    order = ordered_complexes(d)
    pivot = (
        d.pivot_table(index="complex", columns="aa_group", values="percent", fill_value=0)
        .reindex(order)
    )

    for g in AA_GROUPS:
        if g not in pivot.columns:
            pivot[g] = 0
    pivot = pivot[list(AA_GROUPS.keys())]

    families = (
        d.drop_duplicates("complex")
        .set_index("complex")["partner_family"]
        .reindex(order)
    )

    fig, ax = plt.subplots(figsize=(18, 8))
    bottom = [0] * len(pivot)

    x = range(len(pivot))
    for group in pivot.columns:
        vals = pivot[group].values
        ax.bar(
            x,
            vals,
            bottom=bottom,
            label=group,
            color=AA_GROUP_COLORS[group],
            edgecolor="white",
            linewidth=0.4,
        )
        bottom = [a + b for a, b in zip(bottom, vals)]

    # family separators
    last_family = None
    for i, fam in enumerate(families):
        if i > 0 and fam != last_family:
            ax.axvline(i - 0.5, color="black", linewidth=0.8, alpha=0.35)
        last_family = fam

    # family labels above plot
    for fam in GROUP_ORDER:
        idx = [i for i, f in enumerate(families) if f == fam]
        if idx:
            mid = (min(idx) + max(idx)) / 2
            ax.text(mid, 104, fam, ha="center", va="bottom", fontsize=10, fontweight="bold")

    ax.set_ylim(0, 112)
    ax.set_ylabel("Interface amino-acid group percentage (%)", fontsize=12)
    ax.set_title(f"Figure S2-style interface composition — {side}", fontsize=15, fontweight="bold")
    ax.set_xticks(list(x))
    ax.set_xticklabels(pivot.index, rotation=90, fontsize=7)
    ax.legend(
        title="Amino-acid group",
        bbox_to_anchor=(1.01, 1),
        loc="upper left",
        frameon=False,
    )
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()

    png = OUT / f"{filename}.png"
    svg = OUT / f"{filename}.svg"
    fig.savefig(png, dpi=300)
    fig.savefig(svg)
    plt.close(fig)
    print(png)
    print(svg)

def main():
    df = load_long_table()
    plot_side(df, "Receptor side", "figure_s2_groups_receptor_paper_style")
    plot_side(df, "Partner side", "figure_s2_groups_partner_paper_style")

if __name__ == "__main__":
    main()

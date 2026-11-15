"""Ask whether cell-marker RNA accounts for the villous BMP score difference.

The BMP score is the mean of log2(median-ratio normalized count + 1) over the
fixed 27 Reactome genes. Marker scores use that same normalization and the
same OpenBio villous size factors. The adjusted group coefficient is the
EOPE+FGR minus isolated-FGR difference from a linear model with an HC3
interval, matching the Piekos sensitivity.
"""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from analyze import read_counts
from extensions import median_ratio_size_factors, welch_effect


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT.parents[3] / "openbio-datathon-2026" / "track-01" / "data"
SCORES = ROOT.parent / "python_results" / "bmp_donor_scores.csv"
OUT = ROOT / "results"
FIGURES = ROOT / "figures"

PANELS = {
    "endothelial": ["PECAM1", "VWF", "CDH5", "CLDN5"],
    "syncytiotrophoblast": ["CGA", "ERVW-1", "PSG1"],
    "perivascular": ["RGS5", "MCAM", "PDGFRB", "NOTCH3"],
    "receptors": ["BMPR1A", "BMPR1B", "BMPR2"],
    "antagonists": ["GREM2", "NOG"],
    "canonical_smads": ["SMAD1", "SMAD5", "SMAD9"],
}
MODELS = [
    ("unadjusted", []),
    ("endothelial", ["endothelial"]),
    ("syncytiotrophoblast", ["syncytiotrophoblast"]),
    ("both_markers", ["endothelial", "syncytiotrophoblast"]),
    ("age", ["ga_delivery_weeks"]),
    ("both_markers_and_age", ["endothelial", "syncytiotrophoblast", "ga_delivery_weeks"]),
    ("both_markers_and_perivascular", ["endothelial", "syncytiotrophoblast", "perivascular"]),
]
GENE_MARKERS = {
    "BMPR1A": "endothelial",
    "BMPR1B": "endothelial",
    "BMPR2": "endothelial",
    "SMAD1": "endothelial",
    "SMAD5": "endothelial",
    "SMAD9": "endothelial",
    "GREM2": "syncytiotrophoblast",
}
INK = "#182D3D"
MUTED = "#4D6271"
TEAL = "#005E68"
GRID = "#E2EBEE"


def hc3(outcome, design):
    """Group coefficient and an HC3 interval. Column 1 of design is the group indicator."""
    orthogonal, triangular = np.linalg.qr(design)
    beta = np.linalg.solve(triangular, orthogonal.T @ outcome)
    leverage = np.square(orthogonal).sum(axis=1)
    residuals = (outcome - design @ beta) / (1 - leverage)
    weights = np.linalg.solve(triangular, orthogonal.T)
    se = np.linalg.norm(weights[1] * residuals)
    half = stats.t.ppf(0.975, len(outcome) - design.shape[1]) * se
    return {"effect": beta[1], "ci_low": beta[1] - half, "ci_high": beta[1] + half,
            "df": len(outcome) - design.shape[1]}


def panel_score(log_counts, genes):
    missing = [gene for gene in genes if gene not in log_counts.index]
    if missing:
        raise ValueError(f"Marker genes absent from the villous counts: {missing}")
    return log_counts.loc[genes].mean(axis=0)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DATA, help="Organizer Track 1 data directory")
    parser.add_argument("--scores", type=Path, default=SCORES, help="BMP donor scores from run.py")
    args = parser.parse_args()
    counts = read_counts(args.data / "track01_counts_vt.csv")
    meta = pd.read_csv(args.data / "track01_metadata_vt.csv").set_index("sample_id")
    published = pd.read_csv(args.scores)
    published = published.loc[published.cohort.eq("OpenBio")].set_index("sample_id")
    factors = median_ratio_size_factors(counts)
    if not np.allclose(factors.loc[published.index], published.size_factor):
        raise ValueError("Recomputed size factors do not match the published BMP scores")
    log_counts = np.log2(counts.div(factors, axis=1) + 1)
    primary = meta.loc[meta.group.isin(["EOPE+FGR", "FGR"])].copy()
    primary["bmp27"] = published.loc[primary.index, "bmp27"]
    for name, genes in PANELS.items():
        primary[name] = panel_score(log_counts, genes).loc[primary.index]
    primary["combined"] = primary.group.eq("EOPE+FGR").astype(float)
    ages = primary.groupby("group").ga_delivery_weeks.agg(["min", "max"])
    lower, upper = ages["min"].max(), ages["max"].min()
    primary["in_shared_age"] = primary.ga_delivery_weeks.between(lower, upper)

    donor_columns = ["donor_id", "group", "ga_delivery_weeks", "in_shared_age", "bmp27", *PANELS]
    primary[donor_columns].reset_index().to_csv(OUT / "composition_donor_scores.csv", index=False)

    case = primary.group.eq("EOPE+FGR")
    shared = primary.loc[primary.in_shared_age]
    rows = []
    for column in ["bmp27", *PANELS]:
        for comparison, frame in (("openbio_primary", primary), ("openbio_overlap", shared)):
            arm = frame.group.eq("EOPE+FGR")
            rows.append({
                "comparison": comparison, "quantity": column, "method": "welch",
                "n_a": int(arm.sum()), "n_b": int((~arm).sum()),
                **welch_effect(frame.loc[arm, column], frame.loc[~arm, column]),
            })
    for name, covariates in MODELS:
        design = np.column_stack([np.ones(len(primary)), primary.combined, *[primary[c] for c in covariates]])
        rows.append({
            "comparison": "openbio_primary", "quantity": "bmp27", "method": f"hc3:{name}",
            "n_a": int(case.sum()), "n_b": int((~case).sum()),
            "covariates": ";".join(covariates),
            **hc3(primary.bmp27.to_numpy(), design),
        })
    effects = pd.DataFrame(rows)
    effects.to_csv(OUT / "composition_effects.csv", index=False)

    gene_rows = []
    for gene, marker in GENE_MARKERS.items():
        expression = log_counts.loc[gene, primary.index]
        design = np.column_stack([np.ones(len(primary)), primary.combined, primary[marker]])
        correlation = expression.corr(primary[marker])
        gene_rows.append({
            "gene": gene, "method": "welch", "marker_score": marker,
            "marker_pearson_r": correlation, **welch_effect(expression[case], expression[~case]),
        })
        gene_rows.append({
            "gene": gene, "method": "hc3", "marker_score": marker,
            "marker_pearson_r": correlation, **hc3(expression.to_numpy(), design),
        })
    pd.DataFrame(gene_rows).to_csv(OUT / "composition_gene_effects.csv", index=False)

    bmp_check = welch_effect(primary.loc[case, "bmp27"], primary.loc[~case, "bmp27"])["effect"]
    if abs(bmp_check - 0.273062367445176) > 1e-9:
        raise ValueError(f"BMP score difference changed: {bmp_check}")
    _figure(effects)
    _print(effects, primary)


def _figure(effects):
    order = [f"hc3:{name}" for name, _ in MODELS]
    labels = [
        "Unadjusted",
        "Endothelial markers",
        "Syncytiotrophoblast markers",
        "Both marker scores",
        "Delivery age",
        "Both markers and age",
        "Both markers and perivascular",
    ]
    table = effects.set_index("method").loc[order]
    y = np.arange(len(order))[::-1]
    fig, ax = plt.subplots(figsize=(11.2, 6.6))
    fig.subplots_adjust(left=.32, right=.96, top=.78, bottom=.24)
    ax.axvline(0, color="#C5D3D8", linewidth=1)
    ax.errorbar(
        table.effect, y,
        xerr=[table.effect - table.ci_low, table.ci_high - table.effect],
        fmt="o", color=TEAL, ecolor=TEAL, elinewidth=1.6, capsize=3, markersize=6,
    )
    ax.set_yticks(y, labels)
    ax.set_xlabel("EOPE+FGR minus isolated FGR, log₂ BMP score")
    ax.tick_params(length=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    ax.spines["left"].set_color(GRID)
    ax.spines["bottom"].set_color(GRID)
    ax.grid(axis="x", color=GRID)
    fig.text(.065, .95, "Cell-marker RNA does not shrink the villous BMP difference",
             fontsize=16, va="top", weight="bold", color=INK)
    fig.text(.065, .885,
             "OpenBio villi, 7 EOPE+FGR versus 5 isolated FGR  |  The point estimate stays near +0.3 and the intervals widen",
             fontsize=11, va="top", color=MUTED)
    fig.text(.065, .040,
             "Each row is the group coefficient for the same 27-gene score. Marker scores average log₂ normalized counts\n"
             "for the published endothelial or syncytiotrophoblast markers. Age is weeks at delivery. Intervals are HC3.",
             fontsize=9, va="bottom", color=MUTED, linespacing=1.35)
    fig.savefig(FIGURES / "composition.png", dpi=170, facecolor="white")
    plt.close(fig)


def _print(effects, primary):
    print(effects.to_string(index=False, float_format=lambda v: f"{v:.4f}"))
    print("correlations within 12")
    for column in ("endothelial", "syncytiotrophoblast", "perivascular", "ga_delivery_weeks"):
        r = np.corrcoef(primary.bmp27, primary[column])[0, 1]
        print(f"  bmp27 vs {column}: {r:.3f}")
    print("means")
    print(primary.groupby("group")[["bmp27", *PANELS]].mean().round(3))


if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    FIGURES.mkdir(parents=True, exist_ok=True)
    main()

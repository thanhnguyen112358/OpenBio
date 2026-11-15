"""Exploratory re-ranking of the genes that pass in fit A, with and without S24.

Two scores rank the passing genes:
  conservative  the end of the 95% interval nearest zero, |log2FC| - 1.96 * lfcSE
  weighted      0.7 x fold-change score + 0.3 x evidence score, each min-max scaled
                to 0-1 across the passing genes (a teammate's proposal)
Both are computed in fit A and in fit E, which is fit A without S24.
Neither score was prespecified. Run deseq2_fits.py first.
"""

import argparse
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import CASE, CONTROL, GROUP_COLORS, RESULTS_DIR, load_inputs, q1_cohorts
from deseq2_fits import BENCHMARK_GENES, run_fit


INFLUENTIAL_SAMPLE = "S24_VT"
FC_WEIGHT, EVIDENCE_WEIGHT = 0.7, 0.3
TOP_N = 10
METHODS = {
    "conservative": "Conservative fold change",
    "weighted": f"{FC_WEIGHT:.0%} fold change + {EVIDENCE_WEIGHT:.0%} evidence",
}
DESEQ_DIR = RESULTS_DIR / "deseq2"
OUT_DIR = RESULTS_DIR / "candidates"


def minmax(values):
    return (values - values.min()) / (values.max() - values.min())


def score(results):
    """Both scores and their ranks for the genes passing in one fit."""
    hits = results[results["passes"]].copy()
    abs_lfc = hits["log2FoldChange"].abs()
    evidence = -np.log10(hits["padj"])
    hits["conservative_abs_log2fc"] = abs_lfc - 1.96 * hits["lfcSE"]
    hits["weighted_score"] = FC_WEIGHT * minmax(abs_lfc) + EVIDENCE_WEIGHT * minmax(evidence)
    hits["rank_conservative"] = hits["conservative_abs_log2fc"].rank(ascending=False, method="min").astype(int)
    hits["rank_weighted"] = hits["weighted_score"].rank(ascending=False, method="min").astype(int)
    return hits


def top(ranked, method):
    return ranked.sort_values(f"rank_{method}").head(TOP_N)


def normalized_counts(counts, size_factors, cohort, genes):
    return counts.loc[genes, cohort.index] / size_factors[cohort.index]


def s24_position(normalized, cohort):
    """Where S24 sits between the other controls and the EOPE cases, gene by gene.

    On a log scale, 0 is the median of the other four controls and 1 is the EOPE median.
    """
    logged = np.log2(normalized + 1)
    cases = logged.loc[:, cohort["group"].eq(CASE)].median(axis=1)
    others = logged.loc[:, cohort["group"].eq(CONTROL) & (cohort.index != INFLUENTIAL_SAMPLE)].median(axis=1)
    return ((logged[INFLUENTIAL_SAMPLE] - others) / (cases - others)).round(2)


# ---------- plots ----------

def plot_rankings(ranked_by_fit, output_file):
    fig, axes = plt.subplots(1, len(ranked_by_fit), figsize=(7.5 * len(ranked_by_fit), 6))
    styles = {"conservative": ("tab:blue", "o", 80), "weighted": ("tab:orange", "s", 45)}
    for ax, (label, ranked) in zip(axes, ranked_by_fit.items()):
        y = -np.log10(ranked["padj"])
        ax.scatter(ranked["log2FoldChange"], y, s=6, color="lightgray",
                   label=f"{len(ranked)} genes passing in this fit")
        labelled = set()
        for method, (color, marker, size) in styles.items():
            chosen = top(ranked, method)
            labelled |= set(chosen.index)
            ax.scatter(chosen["log2FoldChange"], y[chosen.index], s=size, marker=marker, facecolors="none",
                       edgecolors=color, linewidths=1.3, label=f"Top {TOP_N}: {METHODS[method].lower()}")
        for gene in sorted(labelled):
            ax.annotate(gene, (ranked.loc[gene, "log2FoldChange"], y[gene]), fontsize=7,
                        xytext=(4, 3), textcoords="offset points")
        ax.axvline(0, lw=0.5, color="black")
        ax.set_xlabel(f"log2 Fold Change ({CASE} vs. {CONTROL})")
        ax.set_ylabel("Strength of Evidence (-log10 Adjusted p-value)")
        ax.set_title(label)
        ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=1, frameon=False, fontsize=8)
    fig.suptitle("Which Passing Genes Each Ranking Puts in Its Top 10")
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_candidate_expression(normalized, cohort, ranked, title, output_file):
    genes = list(normalized.index)
    ncols = 5
    nrows = math.ceil(len(genes) / ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 3.4 * nrows), squeeze=False)
    is_case = cohort["group"].eq(CASE)
    for ax, gene in zip(axes.flat, genes):
        values = normalized.loc[gene] + 1
        groups = [(CASE, values[is_case]), (CONTROL, values[~is_case].drop(INFLUENTIAL_SAMPLE))]
        for position, (group, shown) in enumerate(groups):
            x = position + np.linspace(-0.15, 0.15, len(shown))
            ax.scatter(x, shown, s=22, color=GROUP_COLORS[group])
        ax.scatter(1.25, values[INFLUENTIAL_SAMPLE], marker="*", s=150, color="crimson", zorder=4)
        ax.set_yscale("log")
        ax.set_xticks([0, 1], [CASE, CONTROL])
        ax.set_xlim(-0.5, 1.6)
        row = ranked.loc[gene]
        ax.set_title(f"{gene}\nconservative #{row['rank_conservative']}, weighted #{row['rank_weighted']}",
                     fontsize=8)
    for ax in axes.flat[len(genes):]:
        ax.axis("off")
    for ax in axes[:, 0]:
        ax.set_ylabel("Depth-Normalized Count + 1\n(Log Scale)", fontsize=8)
    fig.suptitle(f"{title}. S24 Is the Red Star Beside the Other Four Controls.", fontsize=11)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-cpus", type=int, default=1, help="Worker processes for PyDESeq2")
    args = parser.parse_args()

    counts, metadata, _ = load_inputs()
    full = q1_cohorts(metadata)[0]["eope_vs_iptb"]
    fit_a = pd.read_csv(DESEQ_DIR / "fit_A_full.csv", index_col=0)
    across = pd.read_csv(DESEQ_DIR / "fit_A_hits_across_fits.csv", index_col=0)
    size_factors = pd.read_csv(DESEQ_DIR / "cooks_by_sample_fit_A.csv", index_col=0)["size_factor"]
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    print("Fitting fit E: fit A without S24 ...", flush=True)
    fit_e, _, _ = run_fit(counts, full.drop(INFLUENTIAL_SAMPLE), "~group", args.n_cpus)
    fit_e.to_csv(OUT_DIR / "fit_E_without_S24.csv", float_format="%.6g")
    ranked_a, ranked_e = score(fit_a), score(fit_e)

    genes = ranked_a.index.union(ranked_e.index)
    position = s24_position(normalized_counts(counts, size_factors, full, genes), full)

    table = ranked_a[["baseMean", "log2FoldChange", "ci95_low", "ci95_high", "padj",
                      "conservative_abs_log2fc", "weighted_score", "rank_conservative", "rank_weighted",
                      "max_cooks_sample"]].copy()
    table["n_fits_passing_A_to_D"] = across["n_fits_passing"]
    without = fit_e.reindex(table.index)
    table["log2FoldChange_without_S24"] = without["log2FoldChange"]
    table["padj_without_S24"] = without["padj"]
    table["passes_without_S24"] = without["passes"].fillna(False).astype(bool)
    for method in METHODS:
        table[f"rank_{method}_without_S24"] = ranked_e[f"rank_{method}"].reindex(table.index).astype("Int64")
    table["s24_position"] = position
    table["is_benchmark_gene"] = table.index.isin(BENCHMARK_GENES)
    table.sort_values("rank_conservative").to_csv(OUT_DIR / "candidate_ranking.csv", float_format="%.4g")

    rows = []
    for fit_label, ranked in (("A_with_S24", ranked_a), ("E_without_S24", ranked_e)):
        for method in METHODS:
            for gene, row in top(ranked, method).iterrows():
                rows.append({
                    "fit": fit_label, "method": method, "rank": row[f"rank_{method}"], "gene": gene,
                    "baseMean": row["baseMean"], "log2FoldChange": row["log2FoldChange"],
                    "ci95_low": row["ci95_low"], "ci95_high": row["ci95_high"], "padj": row["padj"],
                    "s24_position": position[gene],
                    "passes_in_fit_A": bool(fit_a["passes"].get(gene, False)),
                    "n_fits_passing_A_to_D": across["n_fits_passing"].get(gene, 0),
                })
    tops = pd.DataFrame(rows)
    tops.to_csv(OUT_DIR / "top10_by_method.csv", index=False, float_format="%.4g")

    plot_rankings({"Fit A, all 13 pregnancies": ranked_a, "Fit E, without S24": ranked_e},
                  OUT_DIR / "ranking_comparison.png")
    top_genes = {name: list(dict.fromkeys(gene for method in METHODS for gene in top(ranked, method).index))
                 for name, ranked in (("A", ranked_a), ("E", ranked_e))}
    plot_candidate_expression(normalized_counts(counts, size_factors, full, top_genes["A"]), full, ranked_a,
                              "Top Candidate Genes in Each Pregnancy, Fit A",
                              OUT_DIR / "top_candidates_expression.png")
    only_without = [gene for gene in top_genes["E"] if gene not in top_genes["A"]]
    plot_candidate_expression(normalized_counts(counts, size_factors, full, only_without), full, ranked_e,
                              "Genes in the Top 10 Only Without S24 (Ranks Are From Fit E)",
                              OUT_DIR / "top_candidates_only_without_s24.png")

    print(f"Passing genes: fit A {len(ranked_a)}, fit E {len(ranked_e)}")
    print(tops.to_string(index=False))
    for method in METHODS:
        a, e = set(top(ranked_a, method).index), set(top(ranked_e, method).index)
        print(f"{method}: {len(a & e)} of {TOP_N} top genes shared with and without S24")
    print("S24 position across fit A passing genes:", table["s24_position"].describe().round(2).to_dict())


if __name__ == "__main__":
    main()

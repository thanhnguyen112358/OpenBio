"""Gene-level DESeq2 fits for Q1: EOPE without FGR vs. iPTB in placental villi.

Fits four raw-count models with PyDESeq2 and compares them:
  A  all 13 donors, ~ group                       (primary)
  B  all 13 donors, ~ delivery age + group        (delivery-age adjusted)
  C  shared delivery-age window, ~ group          (5 vs. 5 sensitivity)
  D  fit A without S22, ~ group                   (donor-influence check)
"""

import argparse
import json
import platform
from importlib.metadata import version

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats
from scipy.stats import f as f_dist

from common import (CASE, CONTROL, FLAGGED_SAMPLE, GROUP_COLORS, GROUP_ORDER, RESULTS_DIR,
                    load_inputs, q1_cohorts)


# Fixed on November 14, 2026, before any gene-level result was inspected.
PADJ_CUTOFF = 0.05
LFC_CUTOFF = 1.0          # |log2 fold change| >= 1, i.e. at least a 2-fold difference
BENCHMARK_GENES = ["FLT1", "ENG", "LEP", "INHBA"]    # prespecified in the Q1 plan
WELL_MEASURED_BASEMEAN = 10                          # organizer workflow's caution level
CUTOFF_GRID = [(0.01, 1.0), (0.05, 1.0), (0.10, 1.0), (0.05, 0.0), (0.05, 0.585), (0.05, 2.0)]

OUT_DIR = RESULTS_DIR / "deseq2"
FIT_LABELS = {
    "A_full": "A: all 13 donors",
    "B_age_adjusted": "B: 13 donors, delivery-age adjusted",
    "C_shared_age": "C: shared delivery-age window",
    "D_without_S22": "D: fit A without S22",
}
FIT_COLORS = {"A_full": "black", "B_age_adjusted": "tab:blue",
              "C_shared_age": "tab:green", "D_without_S22": "crimson"}


def run_fit(counts, cohort, design, n_cpus, contrast=(CASE, CONTROL)):
    """Fit DESeq2 on raw counts and return (per-gene results, per-sample Cook's summary).

    `contrast` is (numerator group, denominator group); the default is the Q1 comparison.
    """
    inference = DefaultInference(n_cpus=n_cpus)
    # `seq_run` (sequencing run, see batch_audit.py) is carried along only when the caller adds it.
    metadata = cohort[[c for c in ("group", "ga_delivery_weeks", "fetal_sex", "seq_run") if c in cohort]].copy()
    metadata["group"] = metadata["group"].astype(str)
    # Centering delivery age leaves the group effect unchanged and lets every gene's fit
    # converge. With raw weeks (about 26 to 33), 56 genes failed to converge in fit B.
    metadata["ga_delivery_weeks_centered"] = (
        metadata["ga_delivery_weeks"] - metadata["ga_delivery_weeks"].mean())
    dds = DeseqDataSet(
        counts=counts.loc[:, cohort.index].T.astype(int),   # PyDESeq2 wants samples as rows
        metadata=metadata,
        design=design,
        inference=inference,
        quiet=True,
    )
    dds.deseq2()
    stats = DeseqStats(dds, contrast=["group", *contrast], inference=inference, quiet=True)
    stats.summary()

    results = stats.results_df.copy()
    results.index.name = "gene"
    results["ci95_low"] = results["log2FoldChange"] - 1.96 * results["lfcSE"]
    results["ci95_high"] = results["log2FoldChange"] + 1.96 * results["lfcSE"]
    results["passes"] = passes(results, PADJ_CUTOFF, LFC_CUTOFF)

    # Cook's distance: how strongly one sample's count pulls a gene's fitted model.
    cooks = pd.DataFrame(dds.layers["cooks"], index=dds.obs_names, columns=dds.var_names)
    has_cooks = cooks.notna().any(axis=0)
    results["max_cooks"] = cooks.max(axis=0)
    results["max_cooks_sample"] = cooks.loc[:, has_cooks].idxmax(axis=0).reindex(results.index)
    results["outlier_count_replaced"] = dds.var["replaced"].astype(bool)
    results["pvalue_removed_by_cooks"] = dds.var["_pvalue_cooks_outlier"].astype(bool)
    results["lfc_converged"] = dds.var["_LFC_converged"]

    n_coefficients = dds.obsm["design_matrix"].shape[1]
    cooks_cutoff = f_dist.ppf(0.99, n_coefficients, len(cohort) - n_coefficients)
    tested = results["padj"].notna() | results["pvalue_removed_by_cooks"]
    top = results.loc[tested & has_cooks, "max_cooks_sample"].value_counts()
    by_sample = pd.DataFrame({
        "group": cohort["group"],
        # DESeq2's depth correction: each sample's counts are divided by this number.
        "size_factor": dds.obs["size_factors"].round(3),
        "genes_where_top_cooks_sample": top.reindex(cohort.index).fillna(0).astype(int),
        "genes_above_cooks_cutoff": (cooks.loc[:, has_cooks] > cooks_cutoff).sum(axis=1),
    })
    by_sample["pct_of_genes_where_top"] = (
        by_sample["genes_where_top_cooks_sample"] / by_sample["genes_where_top_cooks_sample"].sum() * 100
    ).round(1)
    return results.sort_values("padj"), by_sample, float(cooks_cutoff)


def passes(results, padj_cutoff, lfc_cutoff):
    return (results["padj"] < padj_cutoff) & (results["log2FoldChange"].abs() >= lfc_cutoff)


def fit_summary(name, cohort, design, results):
    hits = results[results["passes"]]
    return {
        "fit": name,
        "design": design,
        f"n_{CASE}": int(cohort["group"].eq(CASE).sum()),
        f"n_{CONTROL}": int(cohort["group"].eq(CONTROL).sum()),
        "genes_with_reads": int((results["baseMean"] > 0).sum()),
        "genes_with_padj": int(results["padj"].notna().sum()),
        "genes_passing": len(hits),
        f"higher_in_{CASE}": int((hits["log2FoldChange"] > 0).sum()),
        f"lower_in_{CASE}": int((hits["log2FoldChange"] < 0).sum()),
        "genes_padj_below_cutoff": int((results["padj"] < PADJ_CUTOFF).sum()),
        "genes_with_outlier_count_replaced": int(results["outlier_count_replaced"].sum()),
        "genes_pvalue_removed_by_cooks": int(results["pvalue_removed_by_cooks"].sum()),
        "genes_lfc_not_converged": int((results["lfc_converged"] == False).sum()),  # noqa: E712
    }


def cutoff_sensitivity(fits):
    rows = []
    for name, results in fits.items():
        for padj_cutoff, lfc_cutoff in CUTOFF_GRID:
            keep = passes(results, padj_cutoff, lfc_cutoff)
            rows.append({
                "fit": name, "padj_cutoff": padj_cutoff, "abs_log2fc_cutoff": lfc_cutoff,
                "fold_change_cutoff": round(2 ** lfc_cutoff, 2),
                "is_prespecified": (padj_cutoff, lfc_cutoff) == (PADJ_CUTOFF, LFC_CUTOFF),
                "genes_passing": int(keep.sum()),
                f"higher_in_{CASE}": int((keep & (results["log2FoldChange"] > 0)).sum()),
                f"lower_in_{CASE}": int((keep & (results["log2FoldChange"] < 0)).sum()),
            })
    return pd.DataFrame(rows)


def hits_across_fits(fits):
    """Follow the genes that pass in fit A through the three other fits."""
    primary = fits["A_full"]
    table = primary.loc[primary["passes"], ["baseMean", "log2FoldChange", "ci95_low", "ci95_high", "padj"]]
    table = table.add_suffix("_A")
    direction = np.sign(table["log2FoldChange_A"])
    for name in ("B_age_adjusted", "C_shared_age", "D_without_S22"):
        letter = name[0]
        other = fits[name].reindex(table.index)
        table[f"log2FoldChange_{letter}"] = other["log2FoldChange"]
        table[f"padj_{letter}"] = other["padj"]
        table[f"same_direction_{letter}"] = np.sign(other["log2FoldChange"]) == direction
        table[f"passes_{letter}"] = other["passes"].fillna(False).astype(bool)
    table["same_direction_all_fits"] = table.filter(like="same_direction_").all(axis=1)
    table["n_fits_passing"] = 1 + table.filter(like="passes_").sum(axis=1)
    table["is_benchmark_gene"] = table.index.isin(BENCHMARK_GENES)
    return table.sort_values(["n_fits_passing", "padj_A"], ascending=[False, True])


def influence_table(with_sample, without_sample):
    """Gene-by-gene change in the EOPE effect when S22 is left out (fit D minus fit A)."""
    both = with_sample.add_suffix("_A").join(without_sample.add_suffix("_D"), how="inner")
    both = both[both["log2FoldChange_A"].notna() & both["log2FoldChange_D"].notna()]
    table = both[["baseMean_A", "log2FoldChange_A", "ci95_low_A", "ci95_high_A", "padj_A", "passes_A",
                  "log2FoldChange_D", "ci95_low_D", "ci95_high_D", "padj_D", "passes_D",
                  "max_cooks_sample_A"]].copy()
    table["change_in_log2fc"] = table["log2FoldChange_D"] - table["log2FoldChange_A"]
    table["direction_flips"] = np.sign(table["log2FoldChange_D"]) != np.sign(table["log2FoldChange_A"])
    table["D_estimate_outside_A_ci95"] = ~table["log2FoldChange_D"].between(
        table["ci95_low_A"], table["ci95_high_A"])
    table["ci95_width_ratio_D_over_A"] = (
        (table["ci95_high_D"] - table["ci95_low_D"]) / (table["ci95_high_A"] - table["ci95_low_A"]))
    table["is_benchmark_gene"] = table.index.isin(BENCHMARK_GENES)
    return table


def influence_summary(table):
    hits_a = table[table["passes_A"]]
    well = table[table["baseMean_A"] >= WELL_MEASURED_BASEMEAN]
    return pd.Series({
        "genes_compared": len(table),
        "fit_A_hits": len(hits_a),
        "fit_A_hits_still_passing_without_S22": int(hits_a["passes_D"].sum()),
        "fit_A_hits_same_direction_without_S22": int((~hits_a["direction_flips"]).sum()),
        "fit_A_hits_with_D_estimate_outside_A_ci95": int(hits_a["D_estimate_outside_A_ci95"].sum()),
        "fit_A_hits_median_abs_change_in_log2fc": round(hits_a["change_in_log2fc"].abs().median(), 3),
        "fit_A_hits_max_abs_change_in_log2fc": round(hits_a["change_in_log2fc"].abs().max(), 3),
        "fit_A_hits_median_ci95_width_ratio": round(hits_a["ci95_width_ratio_D_over_A"].median(), 3),
        "fit_A_hits_where_S22_is_top_cooks_sample": int(hits_a["max_cooks_sample_A"].eq(FLAGGED_SAMPLE).sum()),
        "genes_passing_only_without_S22": int((table["passes_D"] & ~table["passes_A"]).sum()),
        "pearson_r_log2fc_well_measured_genes": round(
            well["log2FoldChange_A"].corr(well["log2FoldChange_D"]), 3),
    }, name="value")


def donor_expression(counts, metadata, genes, n_cpus):
    """Depth-normalized counts of a few genes in every villous sample, all five groups.

    Uses DESeq2 size factors fitted on all 30 samples, so the groups share one scale.
    """
    dds = DeseqDataSet(
        counts=counts.T.astype(int),
        metadata=metadata[["group"]].astype(str),
        design="~group",                      # required by the constructor; no model is fitted
        inference=DefaultInference(n_cpus=n_cpus),
        quiet=True,
    )
    dds.fit_size_factors()
    normalized = pd.DataFrame(dds.layers["normed_counts"], index=dds.obs_names, columns=dds.var_names)
    table = metadata[["group"]].copy()
    table["size_factor_all_30_samples"] = dds.obs["size_factors"].round(3)
    return table.join(normalized[genes].round(1))


# ---------- plots ----------

def plot_donor_expression(table, genes, output_file):
    fig, axes = plt.subplots(1, len(genes), figsize=(4.3 * len(genes), 5))
    for ax, gene in zip(axes, genes):
        for position, group in enumerate(GROUP_ORDER):
            values = table.loc[table["group"].eq(group), gene].sort_index()
            x = position + np.linspace(-0.25, 0.25, len(values))
            ax.scatter(x, values, s=22, color=GROUP_COLORS.get(group, "grey"), zorder=3)
            ax.hlines(values.median(), position - 0.35, position + 0.35, color="black", lw=1.2, zorder=4)
            for sample, x_sample in zip(values.index, x):
                if sample == FLAGGED_SAMPLE:
                    ax.scatter(x_sample, values[sample], s=120, facecolors="none", edgecolors="crimson", zorder=5)
                if group in (CASE, CONTROL):
                    ax.annotate(sample.replace("_VT", ""), (x_sample, values[sample]), fontsize=5,
                                alpha=0.7, xytext=(2, 3), textcoords="offset points")
        ax.set_yscale("log")
        ax.set_xticks(range(len(GROUP_ORDER)), GROUP_ORDER, fontsize=8)
        ax.set_xlabel("Clinical Group")
        ax.set_title(gene)
    axes[0].set_ylabel("Depth-Normalized Read Count (Log Scale)")
    fig.suptitle(f"Expression of the Benchmark Genes in Each of the {len(table)} Pregnancies. "
                 "Black Lines Are Group Medians. S22 Is Circled in Red.", fontsize=10)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_volcano(results, output_file):
    tested = results.dropna(subset=["padj"])
    hits = tested[tested["passes"]]
    y = lambda frame: -np.log10(frame["padj"].clip(lower=1e-300))  # noqa: E731
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    ax.scatter(tested["log2FoldChange"], y(tested), s=4, color="lightgray", label="Did not pass")
    ax.scatter(hits["log2FoldChange"], y(hits), s=7, color="crimson",
               label=f"{len(hits)} genes pass (padj < {PADJ_CUTOFF}, at least 2-fold)")
    ax.axhline(-np.log10(PADJ_CUTOFF), ls="--", lw=0.8, color="black")
    for cut in (-LFC_CUTOFF, LFC_CUTOFF):
        ax.axvline(cut, ls="--", lw=0.8, color="black")
    for gene in BENCHMARK_GENES:
        if gene in tested.index:
            point = (tested.loc[gene, "log2FoldChange"], y(tested.loc[[gene]]).iloc[0])
            ax.scatter(*point, s=40, facecolors="none", edgecolors="tab:blue")
            ax.annotate(gene, point, xytext=(5, 5), textcoords="offset points", fontsize=9, color="tab:blue")
    ax.set_xlabel(f"log2 Fold Change ({CASE} vs. {CONTROL}; +1 = Twice as High in {CASE})")
    ax.set_ylabel("Strength of Evidence (-log10 Adjusted p-value)")
    ax.set_title(f"Gene Activity Differences, {CASE} without FGR vs. {CONTROL} (Fit A, n=13)")
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_benchmarks(fits, output_file):
    fig, axes = plt.subplots(1, len(BENCHMARK_GENES), figsize=(3.2 * len(BENCHMARK_GENES), 4.2), sharey=True)
    for ax, gene in zip(axes, BENCHMARK_GENES):
        for position, (name, results) in enumerate(fits.items()):
            row = results.loc[gene]
            ax.errorbar(position, row["log2FoldChange"],
                        yerr=[[row["log2FoldChange"] - row["ci95_low"]], [row["ci95_high"] - row["log2FoldChange"]]],
                        fmt="o", color=FIT_COLORS[name], capsize=4,
                        markerfacecolor=FIT_COLORS[name] if row["passes"] else "white")
        ax.axhline(0, lw=0.8, color="black")
        ax.set_xticks(range(len(fits)))
        ax.set_xticklabels([name[0] for name in fits])
        ax.set_xlim(-0.6, len(fits) - 0.4)
        ax.set_title(gene)
        ax.set_xlabel("Fit")
    axes[0].set_ylabel(f"log2 Fold Change, {CASE} vs. {CONTROL}\n(Bars Show 95% Intervals)")
    fig.suptitle("Benchmark Genes Across the Four Fits. Filled Points Pass the Prespecified Cutoffs.\n"
                 + " · ".join(FIT_LABELS.values()), fontsize=10)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_lfc_comparison(fits, others, title, output_file):
    primary = fits["A_full"]
    fig, axes = plt.subplots(1, len(others), figsize=(6 * len(others), 5.5), squeeze=False)
    for ax, name in zip(axes[0], others):
        both = primary.add_suffix("_A").join(fits[name].add_suffix("_other"), how="inner")
        both = both[(both["baseMean_A"] >= WELL_MEASURED_BASEMEAN)
                    & both["log2FoldChange_A"].notna() & both["log2FoldChange_other"].notna()]
        hits = both[both["passes_A"]]
        r = both["log2FoldChange_A"].corr(both["log2FoldChange_other"])
        ax.scatter(both["log2FoldChange_A"], both["log2FoldChange_other"], s=3, color="lightgray",
                   label=f"{len(both):,} well-measured genes (r = {r:.2f})")
        ax.scatter(hits["log2FoldChange_A"], hits["log2FoldChange_other"], s=7, color="crimson",
                   label=f"{len(hits)} of them pass the cutoffs in fit A")
        limit = np.nanmax(np.abs(both[["log2FoldChange_A", "log2FoldChange_other"]].to_numpy())) * 1.05
        ax.plot([-limit, limit], [-limit, limit], lw=0.8, color="black")
        ax.axhline(0, lw=0.5, color="black")
        ax.axvline(0, lw=0.5, color="black")
        ax.set_xlim(-limit, limit)
        ax.set_ylim(-limit, limit)
        ax.set_xlabel(f"log2 Fold Change in Fit {FIT_LABELS['A_full']}")
        ax.set_ylabel(f"log2 Fold Change in Fit {FIT_LABELS[name]}")
        ax.legend(loc="upper left", fontsize=8)
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-cpus", type=int, default=1, help="Worker processes for PyDESeq2")
    args = parser.parse_args()

    counts, metadata, _ = load_inputs()
    cohorts, (low, high) = q1_cohorts(metadata)
    full, shared = cohorts["eope_vs_iptb"], cohorts["shared_age"]
    plan = {
        "A_full": (full, "~group"),
        "B_age_adjusted": (full, "~ga_delivery_weeks_centered + group"),
        "C_shared_age": (shared, "~group"),
        "D_without_S22": (full.drop(FLAGGED_SAMPLE), "~group"),
    }
    missing = [gene for gene in BENCHMARK_GENES if gene not in counts.index]
    if missing:
        raise ValueError(f"Benchmark genes missing from the count table: {missing}")
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fits, summaries, cooks_tables, cooks_cutoffs = {}, [], {}, {}
    for name, (cohort, design) in plan.items():
        print(f"Fitting {FIT_LABELS[name]} ({design}, n={len(cohort)}) ...", flush=True)
        results, by_sample, cooks_cutoff = run_fit(counts, cohort, design, args.n_cpus)
        results.to_csv(OUT_DIR / f"fit_{name}.csv", float_format="%.6g")
        fits[name], cooks_tables[name], cooks_cutoffs[name] = results, by_sample, cooks_cutoff
        summaries.append(fit_summary(name, cohort, design, results))

    summary = pd.DataFrame(summaries).set_index("fit")
    summary.to_csv(OUT_DIR / "summary.csv")
    cutoff_sensitivity(fits).to_csv(OUT_DIR / "cutoff_sensitivity.csv", index=False)
    hits = hits_across_fits(fits)
    hits.to_csv(OUT_DIR / "fit_A_hits_across_fits.csv")

    benchmark = pd.concat({name: results.loc[BENCHMARK_GENES] for name, results in fits.items()},
                          names=["fit", "gene"]).reset_index()
    benchmark = benchmark[["gene", "fit", "baseMean", "log2FoldChange", "lfcSE", "ci95_low", "ci95_high",
                           "padj", "passes", "max_cooks_sample"]].sort_values(["gene", "fit"])
    benchmark.to_csv(OUT_DIR / "benchmark_genes.csv", index=False)
    expression = donor_expression(counts, metadata, BENCHMARK_GENES, args.n_cpus)
    expression.to_csv(OUT_DIR / "benchmark_expression_by_donor.csv")

    influence = influence_table(fits["A_full"], fits["D_without_S22"])
    influence[influence["passes_A"] | influence["passes_D"] | influence["is_benchmark_gene"]].sort_values(
        "change_in_log2fc", key=abs, ascending=False).to_csv(OUT_DIR / "s22_influence.csv")
    s22_summary = influence_summary(influence)
    s22_summary.to_csv(OUT_DIR / "s22_influence_summary.csv", header=True)
    cooks_tables["A_full"].sort_values("genes_where_top_cooks_sample", ascending=False).to_csv(
        OUT_DIR / "cooks_by_sample_fit_A.csv")

    plot_volcano(fits["A_full"], OUT_DIR / "volcano_fit_A.png")
    plot_benchmarks(fits, OUT_DIR / "benchmark_genes.png")
    plot_donor_expression(expression, BENCHMARK_GENES, OUT_DIR / "benchmark_expression_by_donor.png")
    plot_lfc_comparison(fits, ["B_age_adjusted", "C_shared_age"],
                        "Does the EOPE Effect Hold After Accounting for Delivery Age?",
                        OUT_DIR / "lfc_fit_A_vs_age_checks.png")
    plot_lfc_comparison(fits, ["D_without_S22"], "Does S22 Drive the EOPE Effect?",
                        OUT_DIR / "lfc_fit_A_vs_without_S22.png")

    run_info = {
        "padj_cutoff": PADJ_CUTOFF, "abs_log2fc_cutoff": LFC_CUTOFF,
        "contrast": f"{CASE} vs. {CONTROL}",
        "shared_age_window_weeks": [float(low), float(high)],
        "samples": {name: list(cohort.index) for name, (cohort, _) in plan.items()},
        "cooks_cutoff": cooks_cutoffs,
        "python": platform.python_version(),
        "packages": {p: version(p) for p in ("pydeseq2", "pandas", "numpy", "scipy", "scikit-learn", "anndata")},
    }
    (OUT_DIR / "run_info.json").write_text(json.dumps(run_info, indent=2) + "\n")

    print(summary.T.to_string())
    print(benchmark.to_string(index=False))
    print(s22_summary.to_string())
    print(f"Fit A hits: {len(hits)}; same direction in all fits: {int(hits['same_direction_all_fits'].sum())}; "
          f"passing in all four fits: {int((hits['n_fits_passing'] == 4).sum())}")


if __name__ == "__main__":
    main()

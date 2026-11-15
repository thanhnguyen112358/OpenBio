"""Contextual comparisons for Q1 step 4 (exploratory; not the primary analysis).

  G  EOPE+FGR vs. EOPE   Does adding growth restriction change the EOPE pattern?
  H  FGR vs. iPTB        Is the pattern present with growth restriction alone?
  I  EOPE vs. TB         How does a term-birth control change the answer? (no delivery-age overlap)
  J  EOPE+FGR vs. iPTB   Is the primary pattern also seen in EOPE with growth restriction?

Each is an unadjusted `~ group` DESeq2 fit on raw counts of the two groups only, with the same
PyDESeq2 defaults and selection rule as fits A-D. Groups differ in delivery age (see the
summary), so these are directional context. The saved fit A is read, not refitted.

Run from the OpenBio root:
  .\\.venv\\Scripts\\python.exe "PCA/VST DESeq2/context_contrasts.py"
"""

import json
import platform
from importlib.metadata import version

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import CASE, CONTROL, DATA_DIR, RESULTS_DIR, load_inputs
from deseq2_fits import BENCHMARK_GENES, LFC_CUTOFF, PADJ_CUTOFF, WELL_MEASURED_BASEMEAN, run_fit
from pathway_analysis import FDR_CUTOFF, GMT, PRESPECIFIED, SHORT, ranking, read_gmt, run_gsea, symbol_map

OUT_DIR = RESULTS_DIR / "context_contrasts"
CONTRASTS = {
    "G_EOPE+FGR_vs_EOPE": ("EOPE+FGR", "EOPE"),
    "H_FGR_vs_iPTB": ("FGR", "iPTB"),
    "I_EOPE_vs_TB": ("EOPE", "TB"),
    "J_EOPE+FGR_vs_iPTB": ("EOPE+FGR", "iPTB"),
}


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    counts, metadata, _ = load_inputs()
    mapping = symbol_map(pd.read_csv(DATA_DIR / "track01_gene_annotation.csv"))
    fit_a = pd.read_csv(RESULTS_DIR / "deseq2" / "fit_A_full.csv", index_col=0)
    a_hits = fit_a.index[fit_a["passes"]]
    a_up, a_down = fit_a.index[fit_a["passes"] & (fit_a["log2FoldChange"] > 0)], \
        fit_a.index[fit_a["passes"] & (fit_a["log2FoldChange"] < 0)]

    summaries, benchmarks, pathways = [], [], []
    results_by_name = {"A_EOPE_vs_iPTB": fit_a}
    info = {"exploratory": True, "fits": {}}
    for name, (numerator, denominator) in CONTRASTS.items():
        cohort = metadata[metadata["group"].isin([numerator, denominator])]
        print(f"Fitting {name}: {numerator} ({(cohort.group == numerator).sum()}) vs. "
              f"{denominator} ({(cohort.group == denominator).sum()}) ...", flush=True)
        results, cooks, cutoff = run_fit(counts, cohort, "~group", 1, contrast=(numerator, denominator))
        results.to_csv(OUT_DIR / f"fit_{name}.csv", float_format="%.6g")
        cooks.to_csv(OUT_DIR / f"cooks_{name}.csv")
        results_by_name[name] = results
        info["fits"][name] = {"contrast": [numerator, denominator], "samples": list(cohort.index),
                              "cooks_cutoff": cutoff}

        hits = results[results["passes"]]
        ages = cohort.groupby("group")["ga_delivery_weeks"].agg(["min", "max"])
        joined = fit_a[["log2FoldChange", "baseMean"]].join(results[["log2FoldChange", "passes"]],
                                                            lsuffix="_A", rsuffix="_new", how="inner")
        well = joined[joined["baseMean"] >= WELL_MEASURED_BASEMEAN].dropna()
        on_hits = joined.loc[joined.index.intersection(a_hits)].dropna(subset=["log2FoldChange_new"])
        summaries.append({
            "fit": name, "numerator": numerator, "denominator": denominator,
            "n_numerator": int((cohort.group == numerator).sum()),
            "n_denominator": int((cohort.group == denominator).sum()),
            "delivery_weeks_numerator": f"{ages.loc[numerator, 'min']:.1f}-{ages.loc[numerator, 'max']:.1f}",
            "delivery_weeks_denominator": f"{ages.loc[denominator, 'min']:.1f}-{ages.loc[denominator, 'max']:.1f}",
            "genes_with_padj": int(results["padj"].notna().sum()),
            "genes_passing": len(hits),
            "higher_in_numerator": int((hits["log2FoldChange"] > 0).sum()),
            "lower_in_numerator": int((hits["log2FoldChange"] < 0).sum()),
            "genes_with_outlier_count_replaced": int(results["outlier_count_replaced"].sum()),
            "genes_pvalue_removed_by_cooks": int(results["pvalue_removed_by_cooks"].sum()),
            "genes_lfc_not_converged": int((results["lfc_converged"] == False).sum()),  # noqa: E712
            "A_up_hits_also_positive": int((on_hits.loc[on_hits.index.intersection(a_up), "log2FoldChange_new"] > 0).sum()),
            "A_up_hits": len(a_up),
            "A_down_hits_also_negative": int((on_hits.loc[on_hits.index.intersection(a_down), "log2FoldChange_new"] < 0).sum()),
            "A_down_hits": len(a_down),
            "A_hits_also_passing": int(on_hits["passes"].fillna(False).astype(bool).sum()),
            "median_log2fc_of_A_up_hits": round(on_hits.loc[on_hits.index.intersection(a_up), "log2FoldChange_new"].median(), 3),
            "median_log2fc_of_A_down_hits": round(on_hits.loc[on_hits.index.intersection(a_down), "log2FoldChange_new"].median(), 3),
            "pearson_r_with_A_well_measured": round(well["log2FoldChange_A"].corr(well["log2FoldChange_new"]), 3),
        })
        for gene in BENCHMARK_GENES:
            row = results.loc[gene]
            benchmarks.append({"fit": name, "gene": gene, **{k: row[k] for k in (
                "baseMean", "log2FoldChange", "lfcSE", "ci95_low", "ci95_high", "padj", "passes")}})

        print(f"  GSEA preranked on {name} ...", flush=True)
        rank, log = ranking(results, mapping)
        info["fits"][name]["ranking_log"] = log
        gsea, _ = run_gsea(rank)
        gsea.to_csv(OUT_DIR / f"hallmark_{name}.csv", index=False, float_format="%.6g")
        pathways.append(gsea.assign(fit=name))

    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT_DIR / "summary.csv", index=False)
    benchmark = pd.DataFrame(benchmarks)
    saved_a = fit_a.loc[BENCHMARK_GENES, ["baseMean", "log2FoldChange", "lfcSE", "ci95_low", "ci95_high",
                                          "padj", "passes"]].reset_index().assign(fit="A_EOPE_vs_iPTB")
    benchmark = pd.concat([saved_a, benchmark], ignore_index=True)
    benchmark.to_csv(OUT_DIR / "benchmark_genes.csv", index=False, float_format="%.6g")

    saved_gsea = pd.read_csv(RESULTS_DIR / "pathways" / "hallmark_A_full.csv").assign(fit="A_EOPE_vs_iPTB")
    allp = pd.concat([saved_gsea] + pathways, ignore_index=True)
    nes = allp.pivot(index="pathway", columns="fit", values="NES")
    nes = nes[["A_EOPE_vs_iPTB", *CONTRASTS]]
    nes.round(3).to_csv(OUT_DIR / "nes_all_contrasts.csv")
    pre = allp[allp["pathway"].isin(PRESPECIFIED)][["fit", "pathway", "NES", "nominal_p", "fdr_q",
                                                     "set_size_tested", "leading_edge_size"]]
    pre.to_csv(OUT_DIR / "prespecified_pathways.csv", index=False, float_format="%.6g")

    plot_context(benchmark, pre, OUT_DIR / "context_benchmarks_and_pathways.png")

    info.update({"python": platform.python_version(),
                 "packages": {p: version(p) for p in ("pydeseq2", "gseapy", "pandas", "numpy", "scipy")},
                 "selection_rule": {"padj": PADJ_CUTOFF, "abs_log2fc": LFC_CUTOFF},
                 "gene_sets": str(GMT.name)})
    (OUT_DIR / "run_info.json").write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")
    print(summary.T.to_string())
    print(benchmark[["fit", "gene", "log2FoldChange", "padj", "passes"]].to_string(index=False))
    print(pre.to_string(index=False))


def plot_context(benchmark, pre, output_file):
    fits = ["A_EOPE_vs_iPTB", *CONTRASTS]
    labels = ["A  EOPE\nvs iPTB", "G  EOPE+FGR\nvs EOPE", "H  FGR\nvs iPTB", "I  EOPE\nvs TB",
              "J  EOPE+FGR\nvs iPTB"]
    colors = ["black", "tab:purple", "tab:green", "tab:grey", "tab:red"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), gridspec_kw={"width_ratios": [4, 3]})
    ax = axes[0]
    width = 0.15
    for i, fit in enumerate(fits):
        sub = benchmark[benchmark["fit"] == fit].set_index("gene").reindex(BENCHMARK_GENES)
        x = np.arange(len(BENCHMARK_GENES)) + (i - 2) * width
        filled = sub["passes"].astype(bool).to_numpy()
        ax.errorbar(x, sub["log2FoldChange"], yerr=[sub["log2FoldChange"] - sub["ci95_low"],
                                                   sub["ci95_high"] - sub["log2FoldChange"]],
                    fmt="none", ecolor=colors[i], lw=1, capsize=2)
        ax.scatter(x[filled], sub["log2FoldChange"][filled], color=colors[i], s=36, zorder=3)
        ax.scatter(x[~filled], sub["log2FoldChange"][~filled], facecolors="white", edgecolors=colors[i],
                   s=36, zorder=3)
        ax.scatter([], [], color=colors[i], label=labels[i].replace("\n", " "))
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(range(len(BENCHMARK_GENES)), BENCHMARK_GENES)
    ax.set_ylabel("log2 fold change (numerator vs. denominator)\nbars = 95% intervals")
    ax.set_title("Benchmark genes. Filled = passes padj < 0.05 and |log2FC| ≥ 1", fontsize=10)
    ax.legend(fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.08), frameon=False)
    ax.grid(axis="y", lw=0.4, alpha=0.4)

    ax = axes[1]
    for i, fit in enumerate(fits):
        sub = pre[pre["fit"] == fit].set_index("pathway").reindex(PRESPECIFIED)
        x = np.arange(len(PRESPECIFIED)) + (i - 2) * width
        filled = (sub["fdr_q"] < FDR_CUTOFF).to_numpy()
        ax.scatter(x[filled], sub["NES"][filled], color=colors[i], s=40, zorder=3)
        ax.scatter(x[~filled], sub["NES"][~filled], facecolors="white", edgecolors=colors[i], s=40, zorder=3)
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(range(len(PRESPECIFIED)), [SHORT[p].replace(" ", "\n") for p in PRESPECIFIED])
    ax.set_ylabel("NES (+ = toward higher in numerator group)")
    ax.set_title(f"Prespecified Hallmark sets. Filled = FDR q < {FDR_CUTOFF}", fontsize=10)
    ax.grid(axis="y", lw=0.4, alpha=0.4)
    fig.suptitle("Context comparisons (exploratory). Groups differ in delivery age; "
                 "I (vs. term birth) has no delivery-age overlap.", fontsize=10)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

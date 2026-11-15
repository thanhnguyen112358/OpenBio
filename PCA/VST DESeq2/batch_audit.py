"""Sequencing-run audit for the Q1 comparison (exploratory).

Ackerman et al. (2023, Methods) state that the PTB and TB control specimens were isolated and
sequenced earlier, concurrently with samples S16, S20, S21 and S22. The run is not in the released
metadata, so it is inferred here from that sentence:

  earlier run : all iPTB (5), all TB (5), EOPE S16, S20, S21, S22
  later run   : EOPE S06, S07, S18, S19, all EOPE+FGR (7), all FGR (5)

Consequences: EOPE vs. iPTB is partly confounded with run; FGR or EOPE+FGR vs. either control is
fully confounded; a same-run contrast exists (4 EOPE vs. 5 iPTB).

This script
  1. tabulates run against group, read depth and the PCA axes,
  2. fits  K: earlier run only, ~ group (4 vs. 5)
           L: all 13, ~ seq_run + group
           M: EOPE only, ~ seq_run (4 later vs. 4 earlier): what run alone does, within one diagnosis
     and compares each with the saved fit A (not refitted),
  3. runs Hallmark GSEA on K, L, M with 1,000 permutations (fewer than the 10,000 used for A-F,
     for run time; treat q-values near the cutoff with care),
  4. shows the bulk Ackerman ratios and scores by group and run.

Run from the OpenBio root:
  .\\.venv\\Scripts\\python.exe "PCA/VST DESeq2/batch_audit.py"
"""

import json

import gseapy
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

import pathway_analysis as pa
from common import CASE, CONTROL, DATA_DIR, GROUP_ORDER, RESULTS_DIR, load_inputs
from deseq2_fits import BENCHMARK_GENES, WELL_MEASURED_BASEMEAN, run_fit

OUT_DIR = RESULTS_DIR / "batch_audit"
EARLIER_RUN_EOPE = ["S16_VT", "S20_VT", "S21_VT", "S22_VT"]
GSEA_PERMUTATIONS = 1000
SETS_OF_INTEREST = [*pa.PRESPECIFIED, "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MTORC1_SIGNALING",
                    "HALLMARK_INTERFERON_ALPHA_RESPONSE"]


def assign_run(metadata):
    earlier = metadata["group"].isin(["iPTB", "TB"]) | metadata.index.isin(EARLIER_RUN_EOPE)
    return pd.Series(np.where(earlier, "earlier", "later"), index=metadata.index, name="seq_run")


def gsea(results, mapping):
    rank, _ = pa.ranking(results, mapping)
    result = gseapy.prerank(rnk=rank["stat"], gene_sets=str(pa.GMT), permutation_num=GSEA_PERMUTATIONS,
                            min_size=pa.MIN_SIZE, max_size=pa.MAX_SIZE, weight=1, seed=pa.SEED, threads=4,
                            outdir=None, no_plot=True, verbose=False)
    table = result.res2d.rename(columns={"Term": "pathway", "FDR q-val": "fdr_q", "NOM p-val": "nominal_p"})
    table[["NES", "fdr_q", "nominal_p"]] = table[["NES", "fdr_q", "nominal_p"]].astype(float)
    return table[["pathway", "NES", "nominal_p", "fdr_q"]].sort_values("NES", ascending=False)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    counts, metadata, qc = load_inputs()
    metadata = metadata.copy()
    metadata["seq_run"] = assign_run(metadata)
    metadata["assigned_reads"] = qc["assigned_to_genes"]
    pd.crosstab(metadata["group"], metadata["seq_run"]).to_csv(OUT_DIR / "group_by_run.csv")
    depth = metadata.groupby("seq_run")["assigned_reads"].agg(["count", "min", "median", "max"])
    depth.to_csv(OUT_DIR / "depth_by_run.csv")
    u = mannwhitneyu(metadata.loc[metadata.seq_run == "earlier", "assigned_reads"],
                     metadata.loc[metadata.seq_run == "later", "assigned_reads"])

    # PCA of the 13 Q1 donors on log2(CPM + 1), all genes with reads: which axis follows the run?
    q1 = metadata[metadata["group"].isin([CASE, CONTROL])]
    sub = counts[q1.index]
    sub = sub[sub.sum(axis=1) > 0]
    log_cpm = np.log2(sub / sub.sum(axis=0) * 1e6 + 1)
    centered = (log_cpm.T - log_cpm.T.mean()).to_numpy()
    u_, s_, _ = np.linalg.svd(centered, full_matrices=False)
    pcs = pd.DataFrame(u_[:, :4] * s_[:4], index=q1.index, columns=["PC1", "PC2", "PC3", "PC4"])
    explained = (s_ ** 2 / (s_ ** 2).sum())[:4]
    is_later = (q1["seq_run"] == "later").astype(float)
    is_case = (q1["group"] == CASE).astype(float)
    pca_assoc = pd.DataFrame({
        "share_of_variation": explained,
        "abs_r_with_run": [abs(np.corrcoef(pcs[c], is_later)[0, 1]) for c in pcs],
        "abs_r_with_group": [abs(np.corrcoef(pcs[c], is_case)[0, 1]) for c in pcs],
        "abs_r_with_log_depth": [abs(np.corrcoef(pcs[c], np.log(q1["assigned_reads"]))[0, 1]) for c in pcs],
    }, index=pcs.columns).round(3)
    pca_assoc.to_csv(OUT_DIR / "pca_axes_vs_run.csv")
    pcs.join(q1[["group", "seq_run", "assigned_reads"]]).to_csv(OUT_DIR / "pca_scores_q1.csv")

    fit_a = pd.read_csv(RESULTS_DIR / "deseq2" / "fit_A_full.csv", index_col=0)
    a_hits = fit_a.index[fit_a["passes"]]
    robust = pd.read_csv(RESULTS_DIR / "deseq2" / "fit_A_hits_across_fits.csv", index_col=0)
    robust = robust.index[robust["n_fits_passing"] == 4]
    mapping = pa.symbol_map(pd.read_csv(DATA_DIR / "track01_gene_annotation.csv"))

    eope = q1[q1["group"] == CASE].copy()
    eope["group"] = eope["seq_run"]                      # reuse the contrast machinery: later vs. earlier
    plan = {
        "K_same_run_4v5": (q1[q1["seq_run"] == "earlier"], "~group", (CASE, CONTROL)),
        "L_run_adjusted": (q1, "~seq_run + group", (CASE, CONTROL)),
        "M_run_within_EOPE": (eope, "~group", ("later", "earlier")),
    }
    fits, summaries, pathways = {}, [], []
    for name, (cohort, design, contrast) in plan.items():
        print(f"Fitting {name}: {design}, n={len(cohort)}, contrast {contrast[0]} vs. {contrast[1]} ...", flush=True)
        results, _, _ = run_fit(counts, cohort, design, 1, contrast=contrast)
        results.to_csv(OUT_DIR / f"fit_{name}.csv", float_format="%.6g")
        fits[name] = results
        both = fit_a[["log2FoldChange", "baseMean"]].join(results[["log2FoldChange", "padj", "passes"]],
                                                          lsuffix="_A", rsuffix="_new")
        hits = both.loc[a_hits].dropna(subset=["log2FoldChange_new"])
        well = both[both["baseMean"] >= WELL_MEASURED_BASEMEAN].dropna(subset=["log2FoldChange_A", "log2FoldChange_new"])
        slope = np.polyfit(hits["log2FoldChange_A"], hits["log2FoldChange_new"], 1)[0]
        summaries.append({
            "fit": name, "design": design, "contrast": f"{contrast[0]} vs. {contrast[1]}", "n": len(cohort),
            "genes_with_padj": int(results["padj"].notna().sum()), "genes_passing": int(results["passes"].sum()),
            "A_hits_same_direction": int((np.sign(hits["log2FoldChange_new"]) == np.sign(hits["log2FoldChange_A"])).sum()),
            "A_hits_compared": len(hits),
            "A_hits_still_passing": int(results["passes"].reindex(a_hits).fillna(False).astype(bool).sum()),
            "robust_207_still_passing": int(results["passes"].reindex(robust).fillna(False).astype(bool).sum()),
            "A_hits_median_abs_log2fc": round(hits["log2FoldChange_new"].abs().median(), 3),
            "slope_new_on_A_over_A_hits": round(slope, 3),
            "pearson_r_with_A_well_measured": round(well["log2FoldChange_A"].corr(well["log2FoldChange_new"]), 3),
        })
        print(f"  GSEA ({GSEA_PERMUTATIONS} permutations) ...", flush=True)
        table = gsea(results, mapping)
        table.to_csv(OUT_DIR / f"hallmark_{name}.csv", index=False, float_format="%.5g")
        pathways.append(table.assign(fit=name))
    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT_DIR / "summary.csv", index=False)

    bench = pd.concat({"A_saved": fit_a.loc[BENCHMARK_GENES], **{n: f.loc[BENCHMARK_GENES] for n, f in fits.items()}},
                      names=["fit", "gene"]).reset_index()[["fit", "gene", "log2FoldChange", "ci95_low", "ci95_high", "padj", "passes"]]
    bench.to_csv(OUT_DIR / "benchmark_genes.csv", index=False, float_format="%.5g")

    saved_a = pd.read_csv(RESULTS_DIR / "pathways" / "hallmark_A_full.csv")[["pathway", "NES", "nominal_p", "fdr_q"]]
    allp = pd.concat([saved_a.assign(fit="A_saved")] + pathways, ignore_index=True)
    allp.to_csv(OUT_DIR / "hallmark_all_fits_long.csv", index=False, float_format="%.5g")
    nes = allp.pivot(index="pathway", columns="fit", values="NES")
    q = allp.pivot(index="pathway", columns="fit", values="fdr_q")
    focus = (nes.loc[SETS_OF_INTEREST].round(2).astype(str) + " (q=" + q.loc[SETS_OF_INTEREST].round(3).astype(str) + ")")
    focus.to_csv(OUT_DIR / "hallmark_focus.csv")
    run_sets = allp[(allp.fit == "M_run_within_EOPE") & (allp.fdr_q < 0.05)].sort_values("NES")

    # Bulk Ackerman ratios and scores per donor, with size factors fitted on all 30 samples.
    expr = np.log2(pd.read_csv(RESULTS_DIR / "pathways" / "log2_normalized_pathway_genes_all_30.csv", index_col=0)
                   .pipe(lambda d: 2 ** d - 1) + 1)
    donors = metadata[["group", "seq_run", "assigned_reads", "ga_delivery_weeks"]].copy()
    need = ["LDHA", "LDHB", "PDK1", "GOT1"]
    have = [g for g in need if g in expr.columns]
    if len(have) < 4:                                    # LDHB/GOT1 may be absent from the saved gene subset
        from deseq2_fits import donor_expression
        extra = donor_expression(counts, metadata, need + ["FLT1"], n_cpus=1)
        for g in need + ["FLT1"]:
            expr[g] = np.log2(extra[g] + 1)
    donors["log2_LDHA_LDHB"] = expr["LDHA"] - expr["LDHB"]
    donors["log2_PDK1_GOT1"] = expr["PDK1"] - expr["GOT1"]
    donors["log2_FLT1"] = expr["FLT1"]
    scores = pd.read_csv(RESULTS_DIR / "pathways" / "donor_scores_all_30.csv", index_col=0)
    for col in ("Hypoxia (full set)", "Glycolysis (full set)", "Oxidative phosphorylation (full set)"):
        donors[col] = scores[col]
    donors.to_csv(OUT_DIR / "bulk_ratios_and_scores_by_run.csv", float_format="%.5g")
    measures = ["log2_FLT1", "log2_LDHA_LDHB", "log2_PDK1_GOT1", "Hypoxia (full set)", "Glycolysis (full set)",
                "Oxidative phosphorylation (full set)"]
    e = donors[donors.group == CASE]
    c = donors[donors.group == CONTROL]
    by_run = pd.DataFrame({
        "EOPE_earlier_mean": e[e.seq_run == "earlier"][measures].mean(),
        "EOPE_later_mean": e[e.seq_run == "later"][measures].mean(),
        "iPTB_mean (all earlier)": c[measures].mean(),
        "run_gap_within_EOPE (later - earlier)": e[e.seq_run == "later"][measures].mean() - e[e.seq_run == "earlier"][measures].mean(),
        "same_run_EOPE_minus_iPTB": e[e.seq_run == "earlier"][measures].mean() - c[measures].mean(),
        "all_EOPE_minus_iPTB": e[measures].mean() - c[measures].mean(),
    }).round(3)
    by_run.to_csv(OUT_DIR / "bulk_measures_by_run_summary.csv")

    plot_by_run(donors, measures, OUT_DIR / "bulk_measures_by_group_and_run.png")

    info = {"run_inferred_from": "Ackerman et al. 2023 Methods; not in the released metadata",
            "earlier_run": sorted(metadata.index[metadata.seq_run == "earlier"]),
            "later_run": sorted(metadata.index[metadata.seq_run == "later"]),
            "depth_mann_whitney_p": float(u.pvalue), "gsea_permutations": GSEA_PERMUTATIONS}
    (OUT_DIR / "run_info.json").write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")

    pd.set_option("display.width", 240)
    print(pd.crosstab(metadata["group"], metadata["seq_run"]).to_string())
    print(depth.to_string(), "\nMann-Whitney p for depth by run:", round(u.pvalue, 6))
    print(pca_assoc.to_string())
    print(summary.T.to_string())
    print(bench.pivot(index="gene", columns="fit", values="log2FoldChange").round(2).to_string())
    print(bench.pivot(index="gene", columns="fit", values="passes").to_string())
    print(focus.to_string())
    print("Sets with q < 0.05 for run within EOPE (later vs. earlier):")
    print(run_sets[["pathway", "NES", "fdr_q"]].round(3).to_string(index=False))
    print(by_run.to_string())


def plot_by_run(donors, measures, output_file):
    titles = {"log2_FLT1": "FLT1 (log2 normalized count)", "log2_LDHA_LDHB": "log2 LDHA / LDHB",
              "log2_PDK1_GOT1": "log2 PDK1 / GOT1"}
    fig, axes = plt.subplots(2, 3, figsize=(14, 7.5))
    for ax, measure in zip(axes.ravel(), measures):
        for position, group in enumerate(GROUP_ORDER):
            sub = donors[donors.group == group].sort_values(measure)
            x = position + np.linspace(-0.22, 0.22, len(sub))
            later = (sub.seq_run == "later").to_numpy()
            ax.scatter(x[~later], sub[measure][~later], s=38, color="tab:blue", edgecolors="white", zorder=3,
                       label="earlier run" if position == 0 else None)
            ax.scatter(x[later], sub[measure][later], s=38, color="tab:orange", marker="s", edgecolors="white",
                       zorder=3, label="later run" if position == 0 else None)
            if group in (CASE, CONTROL):
                for xs, (sample, row) in zip(x, sub.iterrows()):
                    ax.annotate(sample.replace("_VT", ""), (xs, row[measure]), fontsize=5, alpha=0.7,
                                xytext=(2, 2), textcoords="offset points")
        ax.set_xticks(range(len(GROUP_ORDER)), GROUP_ORDER, fontsize=8)
        ax.set_title(titles.get(measure, measure.replace(" (full set)", " score")), fontsize=10)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", lw=0.4, alpha=0.4)
    handles = [plt.Line2D([], [], marker="o", ls="", color="tab:blue", label="earlier sequencing run"),
               plt.Line2D([], [], marker="s", ls="", color="tab:orange", label="later sequencing run")]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False)
    fig.suptitle("Bulk villous measures by clinical group and sequencing run (run inferred from Ackerman et al. Methods).\n"
                 "All iPTB and TB controls are in the earlier run; all FGR and EOPE+FGR are in the later run.", fontsize=10)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

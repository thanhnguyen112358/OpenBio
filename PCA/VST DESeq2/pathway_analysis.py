"""Hallmark pathway analysis for Q1 (EOPE without FGR vs. iPTB), from saved DESeq2 results.

Prespecified in the Q1 plan (step 3): rank tested genes by the signed Wald statistic and
check Hallmark Hypoxia, Glycolysis and Oxidative Phosphorylation. The other 47 Hallmark
sets are reported as exploratory context.

Nothing is refitted. The ranking comes from the saved fit tables:
  A  results/deseq2/fit_A_full.csv                     (primary)
  B, C, D  results/deseq2/fit_*.csv                    (saved delivery-age and S22 checks)
  E, F  results/additional_checks/fit_*.csv            (exploratory S24 and fetal-sex checks)

Donor-level plots use depth-normalized counts with size factors fitted on all 30 villous
samples (no model is fitted), the same scale as Figure 12 of the README.

Run from the OpenBio root:
  .\\.venv\\Scripts\\python.exe "PCA/VST DESeq2/pathway_analysis.py"
"""

import hashlib
import json
import platform
from importlib.metadata import version

import gseapy
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from common import CASE, CONTROL, DATA_DIR, GROUP_COLORS, GROUP_ORDER, RESULTS_DIR, load_inputs, q1_cohorts
from deseq2_fits import BENCHMARK_GENES, donor_expression

# Fixed before the pathway results were seen.
GMT = DATA_DIR / "track01_gene_sets.gmt"            # MSigDB Hallmark v2024.1.Hs, from get_gene_sets.py
PRESPECIFIED = ["HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS", "HALLMARK_OXIDATIVE_PHOSPHORYLATION"]
PERMUTATIONS = 10_000
SEED = 20261114
MIN_SIZE, MAX_SIZE = 15, 500                       # GSEA defaults
FDR_CUTOFF = 0.05
ANGIOGENIC_GENES = ["FLT1", "ENG"]                 # the two antiangiogenic benchmark genes

OUT_DIR = RESULTS_DIR / "pathways"
FITS = {
    "A_full": RESULTS_DIR / "deseq2" / "fit_A_full.csv",
    "B_age_adjusted": RESULTS_DIR / "deseq2" / "fit_B_age_adjusted.csv",
    "C_shared_age": RESULTS_DIR / "deseq2" / "fit_C_shared_age.csv",
    "D_without_S22": RESULTS_DIR / "deseq2" / "fit_D_without_S22.csv",
    "E_without_S24": RESULTS_DIR / "additional_checks" / "fit_E_without_S24.csv",
    "F_sex_adjusted": RESULTS_DIR / "additional_checks" / "fit_F_sex_adjusted.csv",
}
SHORT = {"HALLMARK_HYPOXIA": "Hypoxia", "HALLMARK_GLYCOLYSIS": "Glycolysis",
         "HALLMARK_OXIDATIVE_PHOSPHORYLATION": "Oxidative phosphorylation"}


def read_gmt(path):
    text = path.read_text(encoding="utf-8")
    if "<html" in text.lower():
        raise ValueError(f"{path} is an HTML page, not a GMT file. Download it again.")
    sets = {}
    for line in text.splitlines():
        fields = line.rstrip("\n").split("\t")
        if len(fields) > 2:
            sets[fields[0]] = [gene for gene in fields[2:] if gene]
    return sets


def symbol_map(annotation):
    """Source symbol -> current HGNC symbol, keeping only unambiguous one-to-one mappings."""
    usable = annotation[annotation["mapping_status"].isin(
        ["approved_symbol", "previous_symbol", "alias_symbol"])].dropna(subset=["hgnc_symbol"])
    return usable.set_index("gene")["hgnc_symbol"]


def ranking(results, mapping):
    """Signed Wald statistic for tested genes (non-missing padj), indexed by HGNC symbol.

    When two source rows map to one HGNC symbol, the row whose source symbol is already the
    approved symbol is kept; otherwise the row with the larger baseMean.
    """
    tested = results[results["padj"].notna() & results["stat"].notna()].copy()
    tested["hgnc"] = mapping.reindex(tested.index)
    log = {"tested_genes": len(tested), "dropped_unmapped_or_ambiguous": int(tested["hgnc"].isna().sum())}
    tested = tested.dropna(subset=["hgnc"])
    tested["is_approved"] = tested.index == tested["hgnc"]
    tested = tested.sort_values(["is_approved", "baseMean"], ascending=False)
    log["dropped_duplicate_hgnc"] = int(tested["hgnc"].duplicated().sum())
    tested = tested.drop_duplicates("hgnc")
    log["ranked_genes"] = len(tested)
    rank = tested.reset_index().set_index("hgnc")[["gene", "stat"]]
    return rank.sort_values("stat", ascending=False), log


def run_gsea(rank):
    result = gseapy.prerank(
        rnk=rank["stat"], gene_sets=str(GMT), permutation_num=PERMUTATIONS, min_size=MIN_SIZE,
        max_size=MAX_SIZE, weight=1, seed=SEED, threads=4, outdir=None, no_plot=True, verbose=False)
    table = result.res2d.rename(columns={"Term": "pathway", "NOM p-val": "nominal_p",
                                         "FDR q-val": "fdr_q", "FWER p-val": "fwer_p",
                                         "Lead_genes": "leading_edge"})
    for column in ("ES", "NES", "nominal_p", "fdr_q", "fwer_p"):
        table[column] = table[column].astype(float)
    table["set_size_tested"] = table["Tag %"].str.split("/").str[1].astype(int)
    table["leading_edge_size"] = table["Tag %"].str.split("/").str[0].astype(int)
    table["direction"] = np.where(table["NES"] > 0, f"higher in {CASE}", f"lower in {CASE}")
    table["prespecified"] = table["pathway"].isin(PRESPECIFIED)
    table = table.drop(columns=["Name", "Tag %", "Gene %"], errors="ignore")
    return table.sort_values("NES", ascending=False).reset_index(drop=True), result


def overlap_check(rank, sets):
    """Hypoxia and Glycolysis share 65 genes. Test each without the other, and the shared part."""
    hypoxia, glycolysis = set(sets["HALLMARK_HYPOXIA"]), set(sets["HALLMARK_GLYCOLYSIS"])
    custom = {"HYPOXIA_without_GLYCOLYSIS_genes": sorted(hypoxia - glycolysis),
              "GLYCOLYSIS_without_HYPOXIA_genes": sorted(glycolysis - hypoxia),
              "HYPOXIA_GLYCOLYSIS_shared_genes": sorted(hypoxia & glycolysis)}
    result = gseapy.prerank(
        rnk=rank["stat"], gene_sets=custom, permutation_num=PERMUTATIONS, min_size=MIN_SIZE,
        max_size=MAX_SIZE, weight=1, seed=SEED, threads=4, outdir=None, no_plot=True, verbose=False)
    table = result.res2d.rename(columns={"Term": "gene_set", "NOM p-val": "nominal_p", "Lead_genes": "leading_edge"})
    table["set_size_tested"] = table["Tag %"].str.split("/").str[1].astype(int)
    table["NES"] = table["NES"].astype(float)
    table["nominal_p"] = table["nominal_p"].astype(float)
    return table[["gene_set", "NES", "nominal_p", "set_size_tested", "leading_edge"]]


def marker_check(fit_a, annotation):
    """Organizer workflow step 7a: fit A effects for the released cell-class marker genes."""
    markers = pd.read_csv(DATA_DIR / "interpretation" / "track01_marker_matrix.csv")
    usable = annotation[annotation["mapping_status"].isin(["approved_symbol", "previous_symbol", "alias_symbol"])]
    hgnc_to_source = usable.dropna(subset=["hgnc_symbol"]).drop_duplicates("hgnc_symbol").set_index("hgnc_symbol")["gene"]
    markers["gene"] = [s if s in fit_a.index else hgnc_to_source.get(s) for s in markers["gene_symbol"]]
    return markers.join(fit_a[["baseMean", "log2FoldChange", "ci95_low", "ci95_high", "padj", "passes"]], on="gene")


def leading_edge_table(gsea_a, rank_a, fits, hits_across):
    rows = []
    for pathway in PRESPECIFIED:
        genes = gsea_a.set_index("pathway").loc[pathway, "leading_edge"].split(";")
        for hgnc in genes:
            source = rank_a.loc[hgnc, "gene"]
            row = {"pathway": pathway, "hgnc_symbol": hgnc, "gene": source}
            for name, results in fits.items():
                letter = name[0]
                hit = results.reindex([source]).iloc[0]
                row[f"log2FoldChange_{letter}"] = hit["log2FoldChange"]
                row[f"padj_{letter}"] = hit["padj"]
                row[f"passes_{letter}"] = bool(hit["passes"]) if pd.notna(hit["passes"]) else False
            row["stat_A"] = rank_a.loc[hgnc, "stat"]
            row["in_207_passing_A_to_D"] = bool(source in hits_across.index
                                                and hits_across.loc[source, "n_fits_passing"] == 4)
            rows.append(row)
    table = pd.DataFrame(rows)
    table["n_fits_passing_of_6"] = table.filter(like="passes_").sum(axis=1)
    table["same_direction_all_6"] = table.filter(like="log2FoldChange_").apply(
        lambda r: np.sign(r).nunique() == 1, axis=1)
    return table


def donor_scores(log_expr, metadata, gene_lists):
    """Mean z-score per donor. Each gene is z-scored across all 30 villous samples."""
    z = log_expr.sub(log_expr.mean()).div(log_expr.std(ddof=1))
    scores = metadata[["group", "ga_delivery_weeks", "fetal_sex"]].copy()
    for label, genes in gene_lists.items():
        present = [g for g in genes if g in z.columns]
        scores[label] = z[present].mean(axis=1)
        scores.attrs[label] = len(present)
    return scores


# ---------- plots ----------

def plot_prespecified_across_fits(across, output_file):
    fits = list(FITS)
    fig, ax = plt.subplots(figsize=(8.5, 4.2))
    width = 0.13
    colors = ["black", "tab:blue", "tab:green", "crimson", "tab:purple", "tab:brown"]
    for i, name in enumerate(fits):
        sub = across[across["fit"] == name].set_index("pathway").reindex(PRESPECIFIED)
        x = np.arange(len(PRESPECIFIED)) + (i - (len(fits) - 1) / 2) * width
        filled = sub["fdr_q"] < FDR_CUTOFF
        ax.scatter(x[filled], sub["NES"][filled], s=55, color=colors[i], zorder=3)
        ax.scatter(x[~filled], sub["NES"][~filled], s=55, facecolors="white", edgecolors=colors[i],
                   linewidths=1.5, zorder=3)
        ax.scatter([], [], s=40, color=colors[i], label=name.replace("_", " ", 1).replace("_", " "))
    ax.axhline(0, color="black", lw=0.8)
    ax.set_xticks(range(len(PRESPECIFIED)), [SHORT[p] for p in PRESPECIFIED])
    ax.set_ylabel(f"Normalized enrichment score (NES)\n+ = genes ranked toward higher in {CASE}")
    ax.set_title(f"Prespecified Hallmark sets, {CASE} vs. {CONTROL}, across six fits\n"
                 f"Filled = FDR q < {FDR_CUTOFF}; open = not", fontsize=10)
    ax.legend(fontsize=7, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.1), frameon=False)
    ax.grid(axis="y", lw=0.4, alpha=0.4)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_all_hallmarks(gsea_a, output_file):
    table = gsea_a.sort_values("NES")
    fig, ax = plt.subplots(figsize=(7.5, 11))
    y = np.arange(len(table))
    significant = table["fdr_q"] < FDR_CUTOFF
    color = np.where(table["NES"] > 0, "tab:orange", "tab:blue")
    ax.hlines(y, 0, table["NES"], color="lightgrey", lw=1.2)
    ax.scatter(table["NES"][significant], y[significant], color=color[significant], s=40, zorder=3,
               label=f"FDR q < {FDR_CUTOFF}")
    ax.scatter(table["NES"][~significant], y[~significant], facecolors="white",
               edgecolors=color[~significant], s=40, zorder=3, label="Not significant")
    labels = [p.replace("HALLMARK_", "").replace("_", " ").title() for p in table["pathway"]]
    ax.set_yticks(y, labels, fontsize=7)
    for tick, pathway in zip(ax.get_yticklabels(), table["pathway"]):
        if pathway in PRESPECIFIED:
            tick.set_fontweight("bold")
    ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel(f"NES (+ = toward higher in {CASE}, − = toward lower in {CASE})")
    ax.set_title(f"All 50 Hallmark sets, fit A ({CASE} 8 vs. {CONTROL} 5)\n"
                 "Prespecified sets in bold; the other 47 are exploratory", fontsize=10)
    ax.legend(loc="lower right", fontsize=8)
    ax.set_ylim(-1, len(table))
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_leading_edge_heatmap(edge, log_expr, metadata, output_file):
    primary = metadata[metadata["group"].isin([CASE, CONTROL])].sort_values(["group", "ga_delivery_weeks"])
    sub = log_expr.loc[primary.index]
    z = sub.sub(sub.mean()).div(sub.std(ddof=1))           # z across the 13 compared donors
    fig, axes = plt.subplots(1, len(PRESPECIFIED), figsize=(15, 9),
                             gridspec_kw={"width_ratios": [max(1, (edge["pathway"] == p).sum()) for p in PRESPECIFIED]})
    for ax, pathway in zip(axes, PRESPECIFIED):
        genes = edge[edge["pathway"] == pathway].sort_values("log2FoldChange_A", ascending=False)["gene"]
        matrix = z[genes].to_numpy()
        image = ax.imshow(matrix, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5, interpolation="nearest")
        ax.set_xticks(range(len(genes)), genes, rotation=90, fontsize=5)
        ax.set_yticks(range(len(primary)),
                      [f"{s.replace('_VT', '')} {g}" for s, g in zip(primary.index, primary["group"])], fontsize=7)
        boundary = (primary["group"] == CASE).sum() - 0.5
        ax.axhline(boundary, color="black", lw=1.5)
        ax.set_title(f"{SHORT[pathway]}\nleading edge, {len(genes)} genes", fontsize=9)
    cbar = fig.colorbar(image, ax=axes, shrink=0.5, pad=0.01)
    cbar.set_label("log2 normalized count, z-score across the 13 donors")
    fig.suptitle("Leading-edge genes of the prespecified sets, donor by donor (fit A ranking). "
                 "Genes ordered by fit A log2 fold change.\nThese genes were selected from the same "
                 "comparison, so group separation here is expected by construction.", fontsize=10)
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_donor_scores(scores, labels, output_file):
    fig, axes = plt.subplots(1, len(labels), figsize=(4.0 * len(labels), 4.4), sharey=True)
    for ax, label in zip(axes, labels):
        for position, group in enumerate(GROUP_ORDER):
            values = scores.loc[scores["group"].eq(group), label].sort_index()
            x = position + np.linspace(-0.22, 0.22, len(values))
            ax.scatter(x, values, s=24, color=GROUP_COLORS.get(group, "grey"), zorder=3)
            ax.hlines(values.median(), position - 0.33, position + 0.33, color="black", lw=1.2, zorder=4)
            if group in (CASE, CONTROL):
                for sample, xs in zip(values.index, x):
                    ax.annotate(sample.replace("_VT", ""), (xs, values[sample]), fontsize=5, alpha=0.7,
                                xytext=(2, 2), textcoords="offset points")
        ax.axhline(0, color="black", lw=0.5)
        ax.set_xticks(range(len(GROUP_ORDER)), GROUP_ORDER, fontsize=8)
        ax.set_title(f"{label}\n({scores.attrs[label]} genes)", fontsize=9)
    axes[0].set_ylabel("Donor score: mean z-score of log2 normalized counts\n(z across all 30 villous samples)")
    fig.suptitle("Donor-level scores in all five groups. Scores use every tested gene in each Hallmark set, "
                 "not only the leading edge. Black lines are group medians.", fontsize=10)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_same_donors(scores, pathway_labels, output_file):
    fig, axes = plt.subplots(1, len(pathway_labels), figsize=(4.6 * len(pathway_labels), 4.4))
    for ax, label in zip(axes, pathway_labels):
        for group in GROUP_ORDER[::-1]:
            sub = scores[scores["group"] == group]
            ax.scatter(sub["Angiogenic (FLT1, ENG)"], sub[label], s=30 if group in (CASE, CONTROL) else 16,
                       color=GROUP_COLORS.get(group, "lightgrey"), label=group, zorder=3 if group in (CASE, CONTROL) else 2,
                       edgecolors="white", linewidths=0.5)
        eope = scores[scores["group"] == CASE]
        rho, _ = spearmanr(eope["Angiogenic (FLT1, ENG)"], eope[label])
        for sample, row in scores[scores["group"].isin([CASE, CONTROL])].iterrows():
            ax.annotate(sample.replace("_VT", ""), (row["Angiogenic (FLT1, ENG)"], row[label]), fontsize=5,
                        alpha=0.7, xytext=(2, 2), textcoords="offset points")
        ax.axhline(0, color="black", lw=0.4)
        ax.axvline(0, color="black", lw=0.4)
        ax.set_xlabel("Angiogenic score (FLT1, ENG)")
        ax.set_ylabel(f"{label} score")
        ax.set_title(f"Within the 8 {CASE} donors: Spearman ρ = {rho:.2f}", fontsize=9)
    axes[0].legend(fontsize=7, loc="best")
    fig.suptitle("Do the angiogenic and metabolic patterns occur in the same donors? "
                 "(descriptive; n = 8 EOPE)", fontsize=10)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sets = read_gmt(GMT)
    missing = [p for p in PRESPECIFIED if p not in sets]
    if len(sets) != 50 or missing:
        raise ValueError(f"Expected the 50 Hallmark sets including {PRESPECIFIED}; got {len(sets)}, missing {missing}")
    annotation = pd.read_csv(DATA_DIR / "track01_gene_annotation.csv")
    mapping = symbol_map(annotation)

    fits = {name: pd.read_csv(path, index_col=0) for name, path in FITS.items()}
    ranks, rank_logs, gsea, across = {}, {}, {}, []
    for name, results in fits.items():
        print(f"GSEA preranked on {name} ...", flush=True)
        ranks[name], rank_logs[name] = ranking(results, mapping)
        gsea[name], _ = run_gsea(ranks[name])
        gsea[name].to_csv(OUT_DIR / f"hallmark_{name}.csv", index=False, float_format="%.6g")
        sub = gsea[name][gsea[name]["prespecified"]].copy()
        sub.insert(0, "fit", name)
        across.append(sub)
    across = pd.concat(across, ignore_index=True)
    across.drop(columns="leading_edge").to_csv(OUT_DIR / "prespecified_across_fits.csv", index=False,
                                               float_format="%.6g")

    # Rank stability of all 50 NES values across the six fits.
    nes = pd.DataFrame({name: table.set_index("pathway")["NES"] for name, table in gsea.items()})
    nes.round(3).to_csv(OUT_DIR / "nes_all_fits.csv")
    significant_a = gsea["A_full"].loc[gsea["A_full"]["fdr_q"] < FDR_CUTOFF, "pathway"]
    stability = pd.DataFrame({
        "NES_A": nes.loc[significant_a, "A_full"],
        "fits_with_same_sign": (np.sign(nes.loc[significant_a]).eq(np.sign(nes.loc[significant_a, "A_full"]), axis=0)).sum(axis=1),
        "fits_with_fdr_below_cutoff": pd.DataFrame({
            n: t.set_index("pathway").reindex(significant_a)["fdr_q"] < FDR_CUTOFF for n, t in gsea.items()}).sum(axis=1),
    }).sort_values("NES_A")
    stability.to_csv(OUT_DIR / "fit_A_significant_sets_across_fits.csv")

    overlap = overlap_check(ranks["A_full"], sets)
    overlap.to_csv(OUT_DIR / "hypoxia_glycolysis_overlap_fit_A.csv", index=False, float_format="%.6g")
    marker_check(fits["A_full"], annotation).to_csv(OUT_DIR / "cell_class_markers_fit_A.csv", index=False,
                                                    float_format="%.6g")

    hits_across = pd.read_csv(RESULTS_DIR / "deseq2" / "fit_A_hits_across_fits.csv", index_col=0)
    edge = leading_edge_table(gsea["A_full"], ranks["A_full"], fits, hits_across)
    edge.to_csv(OUT_DIR / "leading_edge_fit_A.csv", index=False, float_format="%.6g")

    # Donor level: size factors on all 30 samples (no model fitted), then log2(x + 1).
    counts, metadata, _ = load_inputs()
    a_tested = fits["A_full"].index[fits["A_full"]["padj"].notna()]
    set_genes = {}
    for pathway in PRESPECIFIED:
        members = set(sets[pathway])
        set_genes[pathway] = [g for g in ranks["A_full"]["gene"] if mapping.get(g) in members and g in a_tested]
    needed = sorted(set(edge["gene"]) | set(BENCHMARK_GENES) | {g for v in set_genes.values() for g in v})
    normalized = donor_expression(counts, metadata, needed, n_cpus=1)
    log_expr = np.log2(normalized[needed] + 1)
    log_expr.round(4).to_csv(OUT_DIR / "log2_normalized_pathway_genes_all_30.csv")

    gene_lists = {"Angiogenic (FLT1, ENG)": ANGIOGENIC_GENES}
    gene_lists.update({f"{SHORT[p]} (full set)": set_genes[p] for p in PRESPECIFIED})
    gene_lists.update({f"{SHORT[p]} (leading edge)": list(edge.loc[edge["pathway"] == p, "gene"])
                       for p in PRESPECIFIED})
    scores = donor_scores(log_expr, metadata, gene_lists)
    scores.round(4).to_csv(OUT_DIR / "donor_scores_all_30.csv")

    eope = scores[scores["group"] == CASE]
    co = []
    for label in [k for k in gene_lists if k != "Angiogenic (FLT1, ENG)"]:
        rho, p = spearmanr(eope["Angiogenic (FLT1, ENG)"], eope[label])
        co.append({"score": label, "n_EOPE": len(eope), "spearman_rho_with_angiogenic": round(rho, 3),
                   "nominal_p": round(p, 3)})
    pd.DataFrame(co).to_csv(OUT_DIR / "angiogenic_vs_metabolic_within_EOPE.csv", index=False)

    full_labels = [f"{SHORT[p]} (full set)" for p in PRESPECIFIED]
    plot_prespecified_across_fits(across, OUT_DIR / "prespecified_across_fits.png")
    plot_all_hallmarks(gsea["A_full"], OUT_DIR / "hallmark_all_fit_A.png")
    plot_leading_edge_heatmap(edge, log_expr, metadata, OUT_DIR / "leading_edge_heatmap_fit_A.png")
    plot_donor_scores(scores, ["Angiogenic (FLT1, ENG)"] + full_labels, OUT_DIR / "donor_scores_all_groups.png")
    plot_same_donors(scores, full_labels, OUT_DIR / "angiogenic_vs_metabolic_donors.png")

    run_info = {
        "exploratory_vs_prespecified": {"prespecified": PRESPECIFIED,
                                        "exploratory": "the other 47 Hallmark sets; fits E and F; donor scores"},
        "ranking": "signed Wald statistic (stat) of genes with non-missing padj, mapped to HGNC symbols",
        "gene_sets": {"file": str(GMT.relative_to(DATA_DIR.parent.parent.parent)),
                      "release": "MSigDB Hallmark v2024.1.Hs",
                      "sha256": hashlib.sha256(GMT.read_bytes()).hexdigest(), "n_sets": len(sets)},
        "gsea": {"method": "gseapy.prerank (gene-permutation)", "permutations": PERMUTATIONS, "seed": SEED,
                 "weight": 1, "min_size": MIN_SIZE, "max_size": MAX_SIZE},
        "ranking_logs": rank_logs,
        "baseline_fits": {name: str(path.relative_to(RESULTS_DIR)) for name, path in FITS.items()},
        "python": platform.python_version(),
        "packages": {p: version(p) for p in ("gseapy", "pydeseq2", "pandas", "numpy", "scipy")},
    }
    (OUT_DIR / "run_info.json").write_text(json.dumps(run_info, indent=2) + "\n", encoding="utf-8")

    print(across[["fit", "pathway", "NES", "nominal_p", "fdr_q", "leading_edge_size", "set_size_tested"]]
          .to_string(index=False))
    print(gsea["A_full"][["pathway", "NES", "fdr_q", "set_size_tested", "leading_edge_size"]].to_string(index=False))
    print(stability.to_string())
    print(overlap.drop(columns="leading_edge").to_string(index=False))
    print(pd.DataFrame(co).to_string(index=False))
    print(json.dumps(rank_logs, indent=1))


if __name__ == "__main__":
    main()

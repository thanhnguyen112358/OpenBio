"""Is the glycolytic-shift index a placenta-wide donor state or a cell-type-specific one?

Exploratory. In the Admati donor x cell-class pseudobulk, ask whether donors keep the same rank
for Ackerman's LDHA/LDHB (and PDK1/GOT1) log-ratio across cell classes. The statistic is the mean
pairwise Spearman correlation between cell classes, computed across donors ("concordance").

  Part 1  14 organizer cell classes. Observed concordance is compared with a null of 2,000 random
          pairs of genes that are expressed in every class (median CPM >= 20), so that donor-wide
          technical effects shared by all genes are accounted for. Run for early PE donors alone
          (10, of whom 9 delivered by cesarean) and for all 26 donors after removing group means.
          A Hallmark Hypoxia score without the four ratio genes is compared with size-matched
          random gene sets. FLT1 is shown as a cell-type-restricted contrast.
  Part 2  Teammate's 46-label ratio table (Thanh/admati-glycolytic-ratios.csv): the same
          concordance at fine resolution, using their rule (>= 20 cells, >= 1 UMI per gene).

Limits: sample handling (ischemic time), ambient RNA and delivery mode can all create
donor-wide shifts. Expression ratios are not metabolic flux.

Run from the OpenBio root:
  .\\.venv\\Scripts\\python.exe "PCA/VST DESeq2/placenta_wide_concordance.py"
"""

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata

from celltype_glycolysis import CLASS_LABELS, MIN_CELLS, PB
from common import HERE, RESULTS_DIR
from pathway_analysis import GMT, read_gmt

OUT_DIR = RESULTS_DIR / "placenta_wide"
CLASSES = ["syncytiotrophoblast", "villous_cytotrophoblast", "hofbauer_cell", "macrophage", "monocyte",
           "endothelial_cell", "stromal_cell", "pericyte", "smooth_muscle_cell"]   # >= 9 of 10 early PE donors
RATIOS = {"LDHA/LDHB": ("LDHA", "LDHB"), "PDK1/GOT1": ("PDK1", "GOT1")}
N_NULL, SEED, POOL_MIN_CPM = 2000, 20261115, 20
TEAMMATE_CSV = HERE.parent.parent / "Thanh" / "admati-glycolytic-ratios.csv"


def spearman_matrix(matrix):
    """Pairwise Spearman between rows (cell classes) over columns (donors), ignoring missing donors."""
    n = len(matrix)
    out = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(i + 1, n):
            ok = ~np.isnan(matrix[i]) & ~np.isnan(matrix[j])
            if ok.sum() >= 5:
                out[i, j] = out[j, i] = np.corrcoef(rankdata(matrix[i][ok]), rankdata(matrix[j][ok]))[0, 1]
    return out


def concordance(matrix):
    return float(np.nanmean(spearman_matrix(matrix)[np.triu_indices(len(matrix), 1)]))


def load_log_cpm():
    counts = pd.read_csv(PB / "track01_pseudobulk_counts.csv", index_col=0)
    meta = pd.read_csv(PB / "track01_pseudobulk_metadata.csv").set_index("pseudobulk_id").loc[counts.columns]
    keep = meta[(meta["n_cells"] >= MIN_CELLS) & meta["cell_class"].isin(CLASSES)]
    donors = meta.drop_duplicates("donor_id").set_index("donor_id")[
        ["clinical_group", "delivery_mode", "ga_delivery_weeks", "fgr_status"]].sort_index()
    log_cpm, median_cpm = {}, {}
    for cell_class in CLASSES:
        m = keep[keep["cell_class"] == cell_class]
        cpm = counts[m.index] / counts[m.index].sum(axis=0) * 1e6
        cpm.columns = m["donor_id"].to_numpy()
        median_cpm[cell_class] = cpm.median(axis=1)
        log_cpm[cell_class] = np.log2(cpm.reindex(columns=donors.index) + 1)
    return log_cpm, pd.DataFrame(median_cpm), donors


def class_by_donor(log_cpm, donors_used, values_for_class, center_by=None):
    """Matrix (cell class x donor) of a per-profile value; optionally remove clinical-group means."""
    matrix = np.array([values_for_class(log_cpm[c]).reindex(donors_used).to_numpy(float) for c in CLASSES])
    if center_by is not None:
        groups = center_by.reindex(donors_used).to_numpy()
        for g in np.unique(groups):
            cols = groups == g
            matrix[:, cols] -= np.nanmean(matrix[:, cols], axis=1, keepdims=True)
    return matrix


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    log_cpm, median_cpm, donors = load_log_cpm()
    pool = median_cpm.index[(median_cpm >= POOL_MIN_CPM).all(axis=1)].to_numpy()
    hypoxia = [g for g in read_gmt(GMT)["HALLMARK_HYPOXIA"]
               if g in set(pool) and g not in {"LDHA", "LDHB", "PDK1", "GOT1"}]

    def ratio(a, b):
        return lambda x: x.loc[a] - x.loc[b]

    def gene(a):
        return lambda x: x.loc[a]

    def score(genes):
        def f(x):
            sub = x.loc[genes]
            z = sub.sub(sub.mean(axis=1), axis=0).div(sub.std(axis=1, ddof=1), axis=0)
            return z.mean(axis=0)
        return f

    cohorts = {
        "early_PE_only": (donors.index[donors["clinical_group"] == "PE_early"], None),
        "all_26_group_centered": (donors.index, donors["clinical_group"]),
    }
    null_pairs = [tuple(rng.choice(pool, 2, replace=False)) for _ in range(N_NULL)]
    null_genes = rng.choice(pool, 1000, replace=False)
    null_sets = [list(rng.choice(pool, len(hypoxia), replace=False)) for _ in range(300)]

    rows, nulls, matrices = [], {}, {}
    for cohort, (used, center) in cohorts.items():
        build = lambda f: class_by_donor(log_cpm, used, f, center)  # noqa: E731
        pair_null = np.array([concordance(build(ratio(a, b))) for a, b in null_pairs])
        gene_null = np.array([concordance(build(gene(g))) for g in null_genes])
        set_null = np.array([concordance(build(score(s))) for s in null_sets])
        nulls[cohort] = {"pair": pair_null, "gene": gene_null, "set": set_null}
        tests = [(name, "ratio", ratio(*genes), pair_null) for name, genes in RATIOS.items()]
        tests += [(g, "single gene", gene(g), gene_null) for g in ("FLT1", "LDHA", "LDHB", "HK2", "NDRG1")]
        tests += [(f"Hypoxia score ({len(hypoxia)} genes, ratio genes removed)", "gene set", score(hypoxia), set_null)]
        for name, kind, f, null in tests:
            matrix = build(f)
            matrices[(cohort, name)] = matrix
            observed = concordance(matrix)
            rows.append({"cohort": cohort, "n_donors": len(used), "measure": name, "kind": kind,
                         "mean_pairwise_spearman": observed, "null_median": np.median(null),
                         "null_95th_percentile": np.percentile(null, 95),
                         "fraction_of_null_at_or_above": float((null >= observed).mean()),
                         "n_null_draws": len(null)})
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT_DIR / "concordance_summary.csv", index=False, float_format="%.4g")

    early = cohorts["early_PE_only"][0]
    ldh = pd.DataFrame(matrices[("early_PE_only", "LDHA/LDHB")], index=CLASSES, columns=early).T
    ldh.join(donors).to_csv(OUT_DIR / "ldha_ldhb_early_PE_by_class.csv", float_format="%.4g")
    pd.DataFrame(spearman_matrix(ldh.T.to_numpy()), index=CLASSES, columns=CLASSES).to_csv(
        OUT_DIR / "ldha_ldhb_class_by_class_spearman_early_PE.csv", float_format="%.3g")

    # Leave one donor out: does one extreme placenta carry the concordance?
    loo = {d: concordance(ldh.drop(index=d).T.to_numpy()) for d in ldh.index}
    pd.Series(loo, name="concordance_without_donor").to_csv(OUT_DIR / "ldha_ldhb_leave_one_donor_out.csv")

    fine = fine_resolution()
    fine["table"].to_csv(OUT_DIR / "fine_labels_ldha_ldhb_spearman.csv", float_format="%.3g")

    plot_headline(ldh, nulls["early_PE_only"]["pair"], summary, OUT_DIR / "placenta_wide_ldha_ldhb.png")

    info = {"classes": CLASSES, "pool_genes": int(len(pool)), "seed": SEED, "n_null_pairs": N_NULL,
            "leave_one_donor_out_range": [min(loo.values()), max(loo.values())],
            "fine_resolution": {k: v for k, v in fine.items() if k != "table"}}
    (OUT_DIR / "run_info.json").write_text(json.dumps(info, indent=2) + "\n", encoding="utf-8")

    pd.set_option("display.width", 220)
    print(summary.round(3).to_string(index=False))
    print("leave-one-donor-out concordance (LDHA/LDHB, early PE):", {k: round(v, 2) for k, v in loo.items()})
    print(json.dumps(info["fine_resolution"], indent=1))


def fine_resolution():
    """Concordance of log2 LDHA/LDHB across the 46 original labels (teammate's table), early PE donors."""
    t = pd.read_csv(TEAMMATE_CSV)
    t = t[(t["assay"] == "scRNA") & (t["annotation_level"] == "admati_cell_type") & (t["group"] == "PE_early")]
    t = t[(t["n_cells_or_nuclei"] >= MIN_CELLS) & (t["LDHA_umi"] >= 1) & (t["LDHB_umi"] >= 1)]
    wide = t.pivot(index="cell_type", columns="donor", values="log2_LDHA_LDHB")
    wide = wide[wide.notna().sum(axis=1) >= 8]                  # labels seen in at least 8 of 10 donors
    corr = pd.DataFrame(spearman_matrix(wide.to_numpy()), index=wide.index, columns=wide.index)
    lineage = wide.index.str.split("_").str[0]
    upper = np.triu_indices(len(wide), 1)
    same = (lineage.to_numpy()[:, None] == lineage.to_numpy()[None, :])[upper]
    values = corr.to_numpy()[upper]
    return {"table": corr, "n_labels": int(len(wide)), "labels": list(wide.index),
            "mean_pairwise_spearman": float(np.nanmean(values)),
            "fraction_of_pairs_positive": float(np.nanmean(values > 0)),
            "mean_within_lineage": float(np.nanmean(values[same])),
            "mean_between_lineages": float(np.nanmean(values[~same]))}


def plot_headline(ldh, pair_null, summary, output_file):
    """Left: each early-PE placenta as a line across cell classes. Right: observed vs. random gene pairs."""
    z = (ldh - ldh.mean()) / ldh.std(ddof=1)                    # within-class z across the 10 donors
    order = z.mean(axis=1).sort_values().index
    cmap = plt.get_cmap("viridis")
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw={"width_ratios": [2.4, 1]})
    ax = axes[0]
    for rank, donor in enumerate(order):
        ax.plot(range(len(CLASSES)), z.loc[donor], marker="o", ms=4, lw=1.6, color=cmap(rank / (len(order) - 1)),
                label=donor.replace("PE_early_", "Placenta "))
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(range(len(CLASSES)), [CLASS_LABELS[c] for c in CLASSES], fontsize=8)
    ax.set_ylabel("LDHA / LDHB ratio, relative to the other\nearly-preeclampsia placentas (z-score per cell type)")
    ax.set_title("Each line is one early-preeclampsia placenta (n = 10).\n"
                 "Placentas tend to keep their rank across cell types (highest stays highest, lowest stays low).", fontsize=10)
    ax.legend(fontsize=7, ncol=5, loc="upper center", bbox_to_anchor=(0.5, -0.16), frameon=False,
              title="ordered from lowest to highest average", title_fontsize=7)
    ax.spines[["top", "right"]].set_visible(False)

    ax = axes[1]
    row = summary[(summary.cohort == "early_PE_only")].set_index("measure")
    ax.hist(pair_null, bins=40, color="lightgrey", edgecolor="white")
    for name, color in (("LDHA/LDHB", "crimson"), ("PDK1/GOT1", "tab:purple")):
        v = row.loc[name, "mean_pairwise_spearman"]
        ax.axvline(v, color=color, lw=2)
        ax.annotate(f"{name}\n{v:.2f}", (v, ax.get_ylim()[1] * (0.92 if name == "LDHA/LDHB" else 0.7)),
                    ha="right" if name == "LDHA/LDHB" else "left", fontsize=8, color="black",
                    xytext=(-4 if name == "LDHA/LDHB" else 4, 0), textcoords="offset points")
    ax.set_xlabel("Agreement between cell types\n(mean pairwise Spearman across donors)")
    ax.set_ylabel(f"Random gene pairs (n = {len(pair_null):,})")
    ax.set_title(f"LDHA/LDHB agrees across cell types more than\n"
                 f"{100 * (1 - row.loc['LDHA/LDHB', 'fraction_of_null_at_or_above']):.1f}% of random gene pairs",
                 fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    fig.text(0.01, -0.03, "Data: Admati et al. 2023 single-cell pseudobulk (released Track 1 extension). Exploratory. "
             "RNA ratio, not metabolic flux; sample handling could also shift all cell types together.",
             fontsize=7.5, color="#555555")
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

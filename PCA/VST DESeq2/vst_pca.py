"""Delivery age, sequencing QC, and a VST PCA for the Q1 cohorts (EOPE without FGR vs. iPTB).

Plots delivery age in all five clinical groups, then refits the PCA from
PCA/organizer_pca.py with DESeq2's variance-stabilizing transformation (VST)
in place of log2(CPM + 1) and compares the two.
"""

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from sklearn.decomposition import PCA

from common import (CASE, CONTROL, FLAGGED_SAMPLE, GROUP_COLORS, GROUP_ORDER, HTSEQ_COUNTERS,
                    RESULTS_DIR, load_inputs, q1_cohorts)


N_COMPONENTS = 5
COHORT_DIR = RESULTS_DIR / "cohort"
QC_DIR = RESULTS_DIR / "qc"
PCA_DIR = RESULTS_DIR / "vst_pca"
COHORT_LABELS = {
    "eope_vs_iptb": "EOPE without FGR vs. iPTB",
    "shared_age": "EOPE vs. iPTB, shared delivery-age window",
}
# Colors for the categorical PCA colorings, keyed by the coloring's title.
CATEGORY_COLORS = {
    "group": GROUP_COLORS,
    "fetal sex": {"female": "tab:purple", "male": "tab:green"},
}


# ---------- Part 1: delivery age and sequencing QC ----------

def plot_delivery_age(metadata, window, output_file):
    """Delivery age of every pregnancy in all five groups, with the iPTB range shaded."""
    low, high = window
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.axhspan(low, high, color=GROUP_COLORS[CONTROL], alpha=0.12,
               label=f"Delivery-age range of the iPTB controls ({low:.1f}–{high:.1f} weeks)")
    labels = []
    for position, group in enumerate(GROUP_ORDER):
        ages = metadata.loc[metadata["group"].eq(group), "ga_delivery_weeks"].sort_values()
        # Spread pregnancies with similar delivery ages sideways so the dots do not overlap.
        offsets = pd.Series(0.0, index=ages.index)
        for _, close in ages.groupby(ages.round()):
            offsets[close.index] = (np.arange(len(close)) - (len(close) - 1) / 2) * 0.2
        ax.scatter(position + offsets, ages, color=GROUP_COLORS.get(group, "grey"), zorder=3)
        for sample, age in ages.items():
            ax.annotate(sample.replace("_VT", ""), (position + offsets[sample], age), fontsize=6,
                        alpha=0.7, xytext=(3, 3), textcoords="offset points")
        in_range = int(ages.between(low, high).sum())
        labels.append(f"{group}\n{len(ages)} pregnancies\n{in_range} in shaded range")
    ax.set_xticks(range(len(GROUP_ORDER)), labels)
    ax.set_xlim(-0.6, len(GROUP_ORDER) - 0.4)
    ax.set_xlabel("Clinical Group")
    ax.set_ylabel("Gestational Age at Delivery (Weeks)")
    ax.set_title(f"Delivery Age of All {len(metadata)} Pregnancies by Clinical Group")
    ax.legend(loc="upper left", frameon=False)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def sequencing_qc(counts, metadata, qc, cohort):
    """Per-sample HTSeq read accounting. Percentages are of the reads in the HTSeq summary."""
    table = qc.copy()
    total = table[HTSEQ_COUNTERS].sum(axis=1)
    table.insert(0, "group", metadata["group"])
    table.insert(1, "in_q1_comparison", np.where(table.index.isin(cohort.index), "yes", "no"))
    table["reads_in_htseq_summary"] = total
    for counter in HTSEQ_COUNTERS:
        table[f"pct_{counter}"] = (table[counter] / total * 100).round(2)
    # Gene detection: a gene is "detected" in a sample when it has at least one read.
    detected = counts > 0
    table["genes_detected"] = detected.sum(axis=0)
    in_cohort = detected.loc[:, cohort.index]
    only_here = in_cohort & in_cohort.sum(axis=1).eq(1).to_numpy()[:, None]
    table["genes_detected_in_no_other_q1_sample"] = only_here.sum(axis=0).reindex(table.index).astype("Int64")
    return table


def plot_sequencing_qc(table, output_file):
    shown = table[table["in_q1_comparison"].eq("yes")].sort_values("assigned_to_genes")
    labels = [f"{sample} ({group})" for sample, group in zip(shown.index, shown["group"])]
    fig, (left, middle, right) = plt.subplots(1, 3, figsize=(17, 5), sharey=True)

    colors = ["crimson" if sample == FLAGGED_SAMPLE else "grey" for sample in shown.index]
    left.barh(labels, shown["assigned_to_genes"] / 1e6, color=colors)
    left.set_xlabel("Reads Assigned to Genes (Millions)")
    left.set_title("Reads That Counted Toward Gene Expression")

    middle.barh(labels, shown["genes_detected"], color=colors)
    middle.set_xlim(shown["genes_detected"].min() * 0.9, shown["genes_detected"].max() * 1.03)
    middle.set_xlabel("Genes with at Least One Read")
    middle.set_title("How Many Genes Were Detected")

    parts = [("pct_assigned_to_genes", "Assigned to a gene", "tab:green"),
             ("pct_no_feature", "No annotated gene at that position", "tab:red"),
             ("pct_alignment_not_unique", "Matched more than one place", "tab:purple"),
             ("pct_ambiguous", "Overlapped more than one gene", "tab:brown")]
    start = np.zeros(len(shown))
    for column, label, color in parts:
        right.barh(labels, shown[column], left=start, label=label, color=color)
        start += shown[column].to_numpy()
    right.set_xlabel("Share of Reads in the HTSeq Summary (%)")
    right.set_title("Where Each Sample's Reads Went")
    right.legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, frameon=False)
    fig.suptitle("Sequencing QC for the 13 Pregnancies in the Q1 Comparison")
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


# ---------- Part 2: transformations and PCA ----------

def log_cpm(counts):
    # The organizer's Section 3 transformation, as used in PCA/organizer_pca.py.
    return np.log2(counts.div(counts.sum(axis=0), axis=1) * 1e6 + 1).T


def vst(counts, cohort, n_cpus):
    """DESeq2 VST fitted blind to clinical group (intercept-only dispersion trend)."""
    dds = DeseqDataSet(
        counts=counts.T.astype(int),          # PyDESeq2 wants samples as rows
        metadata=cohort[["group"]].astype(str),
        design="~group",                      # required by the constructor; unused below
        inference=DefaultInference(n_cpus=n_cpus),
        quiet=True,
    )
    dds.vst(use_design=False)
    transformed = pd.DataFrame(dds.layers["vst_counts"], index=dds.obs_names, columns=dds.var_names)
    if not np.isfinite(transformed.to_numpy()).all():
        raise ValueError("VST produced non-finite values")
    return transformed


def fit_pca(matrix):
    pca = PCA(n_components=N_COMPONENTS, random_state=0)
    scores = pd.DataFrame(pca.fit_transform(matrix), index=matrix.index,
                          columns=[f"PC{i + 1}" for i in range(N_COMPONENTS)])
    return scores, pca.explained_variance_ratio_


def abs_r(a, b):
    return abs(np.corrcoef(np.asarray(a, dtype=float), np.asarray(b, dtype=float))[0, 1])


def component_table(scores, variance, cohort, log_depth):
    """For each PC: variance explained and |Pearson r| with depth, delivery age, group, and sex."""
    is_case = cohort["group"].eq(CASE).astype(float)
    is_female = cohort["fetal_sex"].eq("female").astype(float)
    rows = []
    for i, pc in enumerate(scores.columns):
        others = scores[pc].drop(FLAGGED_SAMPLE)
        rows.append({
            "component": pc,
            "pct_variance": round(variance[i] * 100, 1),
            "abs_r_log_assigned_depth": round(abs_r(scores[pc], log_depth), 2),
            # Same scores with S22 left out of the correlation (the PCA is not refitted).
            "abs_r_log_assigned_depth_excl_s22": round(
                abs_r(others, log_depth.drop(FLAGGED_SAMPLE)), 2),
            "abs_r_delivery_age": round(abs_r(scores[pc], cohort["ga_delivery_weeks"]), 2),
            "abs_r_group": round(abs_r(scores[pc], is_case), 2),
            "abs_r_fetal_sex": round(abs_r(scores[pc], is_female), 2),
            # How many SDs S22 sits from the mean of the other samples on this axis.
            "s22_z_vs_others": round((scores.loc[FLAGGED_SAMPLE, pc] - others.mean()) / others.std(), 2),
        })
    return pd.DataFrame(rows)


def isolation_table(matrix):
    """Each sample's average distance to the others, relative to the typical sample.

    Uses the full transformed expression matrix, so it does not depend on which
    components PCA happens to pick. A ratio of 1.0 is a typical sample.
    """
    values = matrix.to_numpy()
    distances = np.sqrt(((values[:, None, :] - values[None, :, :]) ** 2).sum(axis=2))
    mean_distance = pd.Series(distances.sum(axis=1) / (len(values) - 1), index=matrix.index)
    ratio = mean_distance / mean_distance.median()
    return pd.DataFrame({"isolation_ratio": ratio.round(3),
                         "isolation_rank": ratio.rank(ascending=False).astype(int)})


def controls_outside_case_range(scores, cohort):
    """Count iPTB controls lying outside the EOPE range on PC1 (the axis sign is arbitrary)."""
    cases = scores.loc[cohort["group"].eq(CASE), "PC1"]
    controls = scores.loc[cohort["group"].eq(CONTROL), "PC1"]
    return int(max((controls < cases.min()).sum(), (controls > cases.max()).sum()))


def draw_panel(ax, scores, variance, color_by, values, fig):
    x, y = scores["PC1"], scores["PC2"]
    if color_by in CATEGORY_COLORS:
        for level, color in CATEGORY_COLORS[color_by].items():
            mask = values.eq(level)
            ax.scatter(x[mask], y[mask], label=level, color=color)
        ax.legend(title=color_by.capitalize(), fontsize=7, title_fontsize=7)
    else:
        points = ax.scatter(x, y, c=values, cmap="viridis")
        fig.colorbar(points, ax=ax, label=color_by)
    ax.scatter(x[FLAGGED_SAMPLE], y[FLAGGED_SAMPLE], s=160, facecolors="none", edgecolors="crimson")
    for sample in scores.index:
        ax.annotate(sample.replace("_VT", ""), (x[sample], y[sample]), fontsize=7, alpha=0.7)
    ax.set_xlabel(f"PC1 ({variance[0]:.0%} of variation)")
    ax.set_ylabel(f"PC2 ({variance[1]:.0%} of variation)")


def plot_views(fits, cohort, colorings, label, folder):
    # One standalone VST plot per coloring, in the style of PCA/organizer_pca.py.
    scores, variance = fits["VST"]
    for name, (title, values) in colorings.items():
        fig, ax = plt.subplots(figsize=(7, 5))
        draw_panel(ax, scores, variance, title, values, fig)
        ax.set_title(f"{label}, VST (n={len(cohort)})")
        fig.tight_layout()
        fig.savefig(folder / f"vst_{name}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)

    # Side-by-side comparison: one row per transformation, one column per coloring.
    fig, axes = plt.subplots(2, len(colorings), figsize=(6 * len(colorings), 9))
    for row, (transform, (scores, variance)) in enumerate(fits.items()):
        for col, (title, values) in enumerate(colorings.values()):
            draw_panel(axes[row, col], scores, variance, title, values, fig)
            axes[row, col].set_title(f"{transform}: colored by {title}")
    fig.suptitle(f"{label} (n={len(cohort)}). S22 is circled in red. Each row is a separate PCA fit.")
    fig.tight_layout()
    fig.savefig(folder / "logcpm_vs_vst.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-cpus", type=int, default=1, help="Worker processes for PyDESeq2")
    args = parser.parse_args()

    counts, metadata, qc = load_inputs()
    cohorts, (low, high) = q1_cohorts(metadata)
    COHORT_DIR.mkdir(parents=True, exist_ok=True)
    QC_DIR.mkdir(parents=True, exist_ok=True)

    if metadata["fetal_sex"].isna().any() or set(metadata["group"]) != set(GROUP_ORDER):
        raise ValueError("Unexpected clinical groups or missing fetal sex in the metadata")
    plot_delivery_age(metadata, (low, high), COHORT_DIR / "delivery_age_by_group.png")
    ages = metadata.groupby("group")["ga_delivery_weeks"]
    print(pd.DataFrame({"n": ages.size(), "youngest": ages.min(), "oldest": ages.max(),
                        "in_iptb_range": ages.apply(lambda a: int(a.between(low, high).sum()))
                        }).loc[GROUP_ORDER].to_string())
    print(pd.crosstab(metadata["group"], metadata["fetal_sex"]).loc[GROUP_ORDER].to_string())

    qc_table = sequencing_qc(counts, metadata, qc, cohorts["eope_vs_iptb"])
    qc_table.to_csv(QC_DIR / "sequencing_qc.csv")
    plot_sequencing_qc(qc_table, QC_DIR / "sequencing_qc.png")
    flagged = qc_table.loc[FLAGGED_SAMPLE]
    print(f"{FLAGGED_SAMPLE}: {flagged['assigned_to_genes'] / 1e6:.2f} million reads assigned to genes "
          f"({flagged['pct_assigned_to_genes']:.1f}%), {flagged['pct_no_feature']:.1f}% no feature, "
          f"{flagged['genes_detected']:,} genes detected, "
          f"{flagged['genes_detected_in_no_other_q1_sample']:,} detected in no other Q1 sample")

    log_depth_all = np.log(qc["assigned_to_genes"])
    component_rows, isolation_rows, summary_rows = [], [], []
    for name, cohort in cohorts.items():
        folder = PCA_DIR / name
        folder.mkdir(parents=True, exist_ok=True)
        # Genes with no reads in this cohort are constant under both transformations.
        cohort_counts = counts.loc[:, cohort.index]
        cohort_counts = cohort_counts[cohort_counts.sum(axis=1) > 0]
        log_depth = log_depth_all.loc[cohort.index]
        matrices = {"log2(CPM + 1)": log_cpm(cohort_counts),
                    "VST": vst(cohort_counts, cohort, args.n_cpus)}

        fits = {}
        for transform, matrix in matrices.items():
            scores, variance = fit_pca(matrix)
            fits[transform] = (scores, variance)
            scores.to_csv(folder / f"scores_{'vst' if transform == 'VST' else 'logcpm'}.csv")
            components = component_table(scores, variance, cohort, log_depth)
            components.insert(0, "transformation", transform)
            components.insert(0, "cohort", name)
            component_rows.append(components)
            isolation = isolation_table(matrix)
            isolation.insert(0, "transformation", transform)
            isolation.insert(0, "cohort", name)
            isolation_rows.append(isolation.reset_index(names="sample_id"))
            summary_rows.append({
                "cohort": name, "transformation": transform, "n_samples": len(cohort),
                "n_genes": matrix.shape[1],
                "controls_outside_case_range_on_pc1": controls_outside_case_range(scores, cohort),
                "n_controls": int(cohort["group"].eq(CONTROL).sum()),
                "s22_isolation_ratio": isolation.loc[FLAGGED_SAMPLE, "isolation_ratio"],
                "s22_isolation_rank": isolation.loc[FLAGGED_SAMPLE, "isolation_rank"],
            })

        colorings = {
            "group": ("group", cohort["group"]),
            "ga_delivery_weeks": ("delivery age (weeks)", cohort["ga_delivery_weeks"]),
            "fetal_sex": ("fetal sex", cohort["fetal_sex"]),
            "assigned_reads": ("reads assigned to genes (millions)",
                               qc.loc[cohort.index, "assigned_to_genes"] / 1e6),
        }
        label = COHORT_LABELS[name]
        if name == "shared_age":
            label += f" ({low:.3f}–{high:.3f} weeks)"
        plot_views(fits, cohort, colorings, label, folder)

    components = pd.concat(component_rows, ignore_index=True)
    components.to_csv(PCA_DIR / "component_associations.csv", index=False)
    pd.concat(isolation_rows, ignore_index=True).to_csv(PCA_DIR / "sample_isolation.csv", index=False)
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(PCA_DIR / "summary.csv", index=False)
    print(components.to_string(index=False))
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()

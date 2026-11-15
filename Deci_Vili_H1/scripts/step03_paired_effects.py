"""Estimate donor-level paired effects and their cohort and normalization sensitivities."""

from importlib.metadata import version

import matplotlib
import numpy as np
import pandas as pd
from rnanorm import TMM
from scipy.stats import gmean, ttest_ind

from step01_audit_cohort import DB_COUNTS, PANEL, PROJECT, SOURCE, VT_COUNTS, verify_source
from step02_deseq2_scores import AUDIT, COHORTS, OUTPUT as STEP2, score_panel

matplotlib.use("Agg")
import matplotlib.pyplot as plt


OUTPUT = PROJECT / "results" / "step03_paired_effects"
PRIMARY_COHORT = "inferred_same_round"
PAIRED = "db_minus_vt_c1q_trem2"
SELECTIVITY = "db_c1q_trem2_minus_broad_myeloid"
SCORE_METRICS = (
    PAIRED,
    "vt_c1q_trem2",
    "db_c1q_trem2",
    "db_broad_myeloid",
    SELECTIVITY,
    "db_minus_vt_c1q",
    "db_minus_vt_trem2",
)
MARKER_METRICS = tuple(f"{tissue}_gene_{gene}" for tissue in ("vt", "db") for gene in PANEL)
COHORT_LABELS = ("Inferred same round", "Planned age overlap", "Strict age overlap", "Full cohort")
METHOD_COLORS = {"PyDESeq2": "#28567A", "TMM": "#AF562C"}
GROUP_COLORS = {"EOPE": "#994438", "FGR": "#235D70"}


def tmm_panel(donors, deseq, tissue, count_file):
    sample_ids = donors[f"{tissue}_sample_id"]
    raw = pd.read_csv(SOURCE / count_file, index_col="gene").loc[:, sample_ids].T
    raw.index = donors.index
    library_size = raw.sum(axis=1)
    matrix = raw.to_numpy(copy=True)
    factor = TMM().fit(matrix).get_norm_factors(matrix)
    effective = library_size * np.asarray(factor)
    # Keep the +1 pseudocount on the same count scale as in Step 2.
    deseq_size_factor = deseq[f"{tissue}_size_factor"]
    scale = effective / gmean(effective) * gmean(deseq_size_factor)
    logged = np.log2(raw.loc[:, list(PANEL)].div(scale, axis=0) + 1)
    scaling = pd.DataFrame(
        {
            "sample_id": sample_ids,
            "tissue": tissue.upper(),
            "library_size": library_size,
            "tmm_factor": factor,
            "tmm_effective_library_size": effective,
            "tmm_size_factor_for_scoring": scale,
            "pydeseq2_size_factor": deseq_size_factor,
        },
        index=donors.index,
    )
    return logged, scaling


def load_scores():
    verify_source()
    deseq = pd.read_csv(STEP2 / "donor_scores.csv", index_col="donor_id")
    audit = pd.read_csv(AUDIT, index_col="donor_id")
    donors = audit.loc[audit.in_all_eope_fgr_pairs.eq(1)]
    membership = [f"in_{cohort}" for cohort in COHORTS]
    metadata = ["group", "ga_delivery_weeks", "sequencing_round", "qc_flags", *membership]
    pd.testing.assert_frame_equal(deseq[metadata], donors[metadata])
    saved_panel = pd.read_csv(STEP2 / "normalized_panel.csv")
    pydeseq_scores, tmm_scores = donors[metadata].copy(), donors[metadata].copy()
    scaling = []
    for tissue, count_file in (("vt", VT_COUNTS), ("db", DB_COUNTS)):
        panel = saved_panel.loc[saved_panel.tissue.eq(tissue.upper())]
        logged = panel.pivot(index="donor_id", columns="gene", values="log2_normalized_plus1")
        estimates = score_panel(logged).add_prefix(f"{tissue}_").loc[donors.index]
        np.testing.assert_allclose(estimates, deseq.loc[donors.index, estimates.columns])
        pydeseq_scores = pydeseq_scores.join(estimates).join(logged.add_prefix(f"{tissue}_gene_"))

        logged_tmm, library_scaling = tmm_panel(donors, deseq, tissue, count_file)
        tmm_scores = tmm_scores.join(score_panel(logged_tmm).add_prefix(f"{tissue}_"))
        tmm_scores = tmm_scores.join(logged_tmm.add_prefix(f"{tissue}_gene_"))
        scaling.append(library_scaling)

    by_method = {"PyDESeq2": pydeseq_scores, "TMM": tmm_scores}
    for rows in by_method.values():
        for component in ("c1q_trem2", "c1q", "trem2"):
            rows[f"db_minus_vt_{component}"] = rows[f"db_{component}"] - rows[f"vt_{component}"]
        rows[SELECTIVITY] = rows["db_c1q_trem2"] - rows["db_broad_myeloid"]
        rows["vt_trem2_raw_count"] = donors["vt_TREM2_raw_count"]
    np.testing.assert_allclose(pydeseq_scores[[PAIRED, SELECTIVITY]], deseq[[PAIRED, SELECTIVITY]])
    return by_method, pd.concat(scaling)


def estimate(rows, metrics):
    eope = rows.loc[rows.group.eq("EOPE"), list(metrics)]
    fgr = rows.loc[rows.group.eq("FGR"), list(metrics)]
    test = ttest_ind(eope.to_numpy(), fgr.to_numpy(), axis=0, equal_var=False, nan_policy="raise")
    interval = test.confidence_interval(confidence_level=0.95)
    return pd.DataFrame(
        {
            "n_eope": len(eope),
            "n_fgr": len(fgr),
            "mean_eope": eope.mean(),
            "mean_fgr": fgr.mean(),
            "estimate_eope_minus_fgr": eope.mean() - fgr.mean(),
            "ci95_low": interval.low,
            "ci95_high": interval.high,
            "welch_df": test.df,
            "p_value_two_sided": test.pvalue,
        }
    ).rename_axis("metric").reset_index()


def contrasts(by_method):
    results, leave_out = [], []
    for method, rows in by_method.items():
        for cohort in COHORTS:
            selected = rows.loc[rows[f"in_{cohort}"].eq(1)]
            effects = estimate(selected, SCORE_METRICS + MARKER_METRICS)
            results.append(effects.assign(normalization=method, cohort=cohort))

        primary = rows.loc[rows[f"in_{PRIMARY_COHORT}"].eq(1)]
        baseline = estimate(primary, SCORE_METRICS).set_index("metric").estimate_eope_minus_fgr
        for donor, record in primary.iterrows():
            change = estimate(primary.drop(donor), SCORE_METRICS).assign(
                normalization=method,
                dropped_donor=donor,
                dropped_group=record.group,
                dropped_donor_qc_flags=record.qc_flags,
                dropped_vt_trem2_raw_count=record.vt_trem2_raw_count,
            )
            change["change_from_primary"] = change.estimate_eope_minus_fgr - change.metric.map(baseline)
            change["sign_reversal"] = change.estimate_eope_minus_fgr * change.metric.map(baseline) < 0
            leave_out.append(change)
    return pd.concat(results, ignore_index=True), pd.concat(leave_out, ignore_index=True)


def save_figure(fig, name):
    fig.savefig(OUTPUT / "figures" / name, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_donors(scores):
    fig, ax = plt.subplots(figsize=(7.6, 5.2), layout="constrained")
    ax.axhline(0, color="#B8BEC3", linewidth=1)
    for position, (group, color) in enumerate(GROUP_COLORS.items()):
        rows = scores.loc[scores.group.eq(group)].sort_values(PAIRED)
        for x, (donor, row) in zip(position + np.linspace(-0.32, 0.32, len(rows)), rows.iterrows()):
            primary = row[f"in_{PRIMARY_COHORT}"]
            ax.scatter(x, row[PAIRED], color=color if primary else "white", edgecolor=color, s=54, zorder=3)
            ax.annotate(donor, (x, row[PAIRED]), xytext=(0, 7), textcoords="offset points",
                        ha="center", fontsize=9, color="#27323A")
        primary_rows = rows.loc[rows[f"in_{PRIMARY_COHORT}"].eq(1), PAIRED]
        ax.hlines(primary_rows.mean(), position - 0.36, position + 0.36, color=color, linewidth=2)
    ax.set(xticks=[0, 1], xticklabels=["EOPE without FGR", "Isolated FGR"], xlim=(-0.55, 1.55),
           ylabel="Decidua minus villi (log₂ score units)")
    ax.set_title("Paired C1Q–TREM2 scores by pregnancy", loc="left", weight="bold", pad=22)
    ax.text(0, 1.02, "Filled points and mean lines show the inferred same-round cohort",
            transform=ax.transAxes, fontsize=9)
    ax.margins(y=0.22)
    save_figure(fig, "paired_donor_scores.png")


def error_bars(ax, rows, y, color, label=None):
    x = rows.estimate_eope_minus_fgr
    ax.errorbar(x, y, xerr=[x - rows.ci95_low, rows.ci95_high - x], color=color,
                label=label, fmt="o", markersize=5, capsize=3)


def plot_components(effects):
    labels = {
        "db_c1q_trem2": "Decidua · C1Q–TREM2",
        "vt_c1q_trem2": "Villi · C1Q–TREM2",
        "db_broad_myeloid": "Decidua · broad myeloid",
        SELECTIVITY: "Decidua · C1Q–TREM2 minus myeloid",
        "db_gene_CD68": "Decidua · CD68",
        "db_gene_LST1": "Decidua · LST1",
        "db_gene_CSF1R": "Decidua · CSF1R",
    }
    fig, ax = plt.subplots(figsize=(8.5, 4.6), layout="constrained")
    primary = effects.loc[effects.cohort.eq(PRIMARY_COHORT) & effects.normalization.eq("PyDESeq2")]
    y = np.arange(len(labels))
    error_bars(ax, primary.set_index("metric").loc[list(labels)], y, METHOD_COLORS["PyDESeq2"])
    ax.axvline(0, color="#B8BEC3", linewidth=1)
    ax.set(yticks=y, yticklabels=list(labels.values()), xlabel="EOPE minus isolated FGR (log₂ units)")
    ax.invert_yaxis()
    ax.set_title("Tissue and marker effects in the inferred same-round cohort", loc="left", weight="bold")
    save_figure(fig, "tissue_and_marker_effects.png")


def plot_sensitivity(effects, leave_out):
    fig, (cohort_ax, omit_ax) = plt.subplots(1, 2, figsize=(12.3, 5.0), layout="constrained")
    paired = effects.loc[effects.metric.eq(PAIRED)]
    omitted = leave_out.loc[leave_out.metric.eq(PAIRED)]
    for offset, (method, color) in zip((-0.13, 0.13), METHOD_COLORS.items()):
        subset = paired.loc[paired.normalization.eq(method)].set_index("cohort").loc[list(COHORTS)]
        labels = [f"{label} ({row.n_eope} vs. {row.n_fgr})"
                  for label, row in zip(COHORT_LABELS, subset.itertuples())]
        error_bars(cohort_ax, subset, np.arange(len(labels)) + offset, color, method)
        rows = omitted.loc[omitted.normalization.eq(method)]
        donor_labels = [
            row.dropped_donor + ("*" if row.dropped_donor_qc_flags != "none" else "")
            + ("†" if row.dropped_vt_trem2_raw_count < 10 else "")
            for row in rows.itertuples()
        ]
        omit_ax.scatter(rows.estimate_eope_minus_fgr, np.arange(len(rows)) + offset, color=color, s=28)
    for axis in (cohort_ax, omit_ax):
        axis.axvline(0, color="#B8BEC3", linewidth=1)
        axis.set_xlabel("EOPE minus FGR in paired score (log₂ units)")
        axis.invert_yaxis()
    cohort_ax.set(yticks=np.arange(len(labels)), yticklabels=labels)
    cohort_ax.set_title("Cohort and normalization · 95% intervals", loc="left", weight="bold")
    cohort_ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=2, fontsize=9)
    omit_ax.set(yticks=np.arange(len(donor_labels)), yticklabels=donor_labels)
    omit_ax.set_title("Primary cohort · one donor omitted", loc="left", weight="bold", pad=22)
    omit_ax.text(0, 1.02, "* BP flag    † Villous TREM2 below 10 reads",
                 transform=omit_ax.transAxes, fontsize=9)
    save_figure(fig, "paired_effect_sensitivity.png")


def main():
    by_method, scaling = load_scores()
    effects, leave_out = contrasts(by_method)
    (OUTPUT / "figures").mkdir(parents=True, exist_ok=True)
    pd.concat(by_method, names=["normalization"]).to_csv(OUTPUT / "donor_scores.csv")
    scaling.to_csv(OUTPUT / "library_scaling.csv")
    effects.to_csv(OUTPUT / "effects.csv", index=False)
    leave_out.to_csv(OUTPUT / "leave_one_out.csv", index=False)
    plot_donors(by_method["PyDESeq2"])
    plot_components(effects)
    plot_sensitivity(effects, leave_out)
    print(f"RNAnorm {version('rnanorm')}; SciPy {version('scipy')}")
    print(effects.loc[effects.metric.eq(PAIRED), [
        "normalization", "cohort", "n_eope", "n_fgr", "estimate_eope_minus_fgr", "ci95_low", "ci95_high"
    ]].to_string(index=False))


if __name__ == "__main__":
    main()

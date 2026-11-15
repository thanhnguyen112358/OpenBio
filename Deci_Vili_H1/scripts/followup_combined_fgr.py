"""Compare the fixed RNA scores in EOPE+FGR, isolated FGR, and EOPE-only bulk samples."""

import matplotlib
import numpy as np
import pandas as pd
from pydeseq2.preprocessing import deseq2_norm_fit, deseq2_norm_transform
from scipy.stats import ttest_ind

from step01_audit_cohort import DB_COUNTS, PANEL, PROJECT, SOURCE, VT_COUNTS, verify_source
from step02_deseq2_scores import AUDIT, OUTPUT as STEP2, score_panel

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


OUTPUT = PROJECT / "results" / "followup_fgr_context"
COMBINED = "EOPE+FGR"
ISOLATED = "FGR"
EOPE_ONLY = "EOPE"
VT_SCORE = "vt_c1q_trem2"
DB_SCORE = "db_c1q_trem2"
PAIRED = "db_minus_vt_c1q_trem2"
SELECTIVITY = "db_c1q_trem2_minus_broad_myeloid"
METRICS = [VT_SCORE, DB_SCORE, PAIRED, "db_broad_myeloid", SELECTIVITY]
STYLES = {
    EOPE_ONLY: ("EOPE without FGR", "#994438", "o"),
    ISOLATED: ("Isolated FGR", "#235D70", "s"),
    COMBINED: ("EOPE with FGR", "#695686", "^"),
}


def project_tissue(donors, reference, tissue, count_file):
    sample_column = f"{tissue}_sample_id"
    available = donors.loc[donors[sample_column].notna(), sample_column]
    counts = pd.read_csv(
        SOURCE / count_file, index_col="gene", usecols=["gene", *available.tolist()]
    ).loc[:, available].T
    counts.index = available.index
    logmeans, mask = deseq2_norm_fit(counts.loc[reference.index].to_numpy())
    normalized, factors = deseq2_norm_transform(counts.to_numpy(), logmeans, mask)
    if not np.isfinite(factors).all() or (factors <= 0).any():
        raise ValueError(f"Invalid {tissue} size factors on the Step 2 reference")
    normalized = pd.DataFrame(normalized, index=counts.index, columns=counts.columns)
    scores = score_panel(np.log2(normalized.loc[:, list(PANEL)] + 1)).add_prefix(f"{tissue}_")
    scores[f"{tissue}_size_factor"] = factors
    np.testing.assert_allclose(scores.loc[reference.index], reference.loc[:, scores.columns])
    return scores


def compare(numerator, reference, metrics, label, age_window=None):
    a, b = numerator[metrics], reference[metrics]
    test = ttest_ind(a.to_numpy(), b.to_numpy(), axis=0, equal_var=False, nan_policy="raise")
    ci = test.confidence_interval()
    age_min, age_max = age_window or (np.nan, np.nan)
    return pd.DataFrame({
        "comparison": label,
        "numerator": numerator.group.iloc[0],
        "reference": reference.group.iloc[0],
        "n_numerator": len(a),
        "n_reference": len(b),
        "numerator_donor_ids": " ".join(numerator.index),
        "reference_donor_ids": " ".join(reference.index),
        "age_filter_min": age_min,
        "age_filter_max": age_max,
        "numerator_mean": a.mean(),
        "reference_mean": b.mean(),
        "estimate": a.mean() - b.mean(),
        "ci95_low": ci.low,
        "ci95_high": ci.high,
    }).rename_axis("metric").reset_index()


def age_restricted(numerator, reference):
    a, b = numerator.ga_delivery_weeks, reference.ga_delivery_weeks
    low = max(a.min(), b.min())
    high = min(a.max(), b.max())
    return numerator.loc[a.between(low, high)], reference.loc[b.between(low, high)], (low, high)


def plot_scores(same_round):
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.2), sharex=True, layout="constrained")
    for ax, field, title in zip(axes, [VT_SCORE, DB_SCORE], ["Villi", "Decidua"]):
        for group, (_, color, marker) in STYLES.items():
            selected = same_round.loc[same_round.group.eq(group) & same_round[field].notna()]
            for donor, row in selected.iterrows():
                ax.scatter(row.ga_delivery_weeks, row[field], marker=marker,
                           facecolor=color if row.complete_pair else "white", edgecolor=color, s=57)
                if donor == "S15":
                    ax.annotate(donor, (row.ga_delivery_weeks, row[field]), xytext=(8, 3),
                                textcoords="offset points", fontsize=8)
        ax.set_title(title, loc="left", weight="bold")
        ax.set_xlabel("Delivery age (weeks)")
        ax.spines[["top", "right"]].set_visible(False)
        ax.margins(x=0.12, y=0.18)
    axes[0].set_ylabel("C1Q–TREM2 score (log₂ units)")
    axes[1].set_ylabel("C1Q–TREM2 score (log₂ units)")
    fig.suptitle("C1Q–TREM2 by clinical group in one inferred sequencing round", weight="bold")
    handles = [
        Line2D([], [], marker=marker, color=color, linestyle="none", label=label)
        for label, color, marker in STYLES.values()
    ]
    handles.append(Line2D([], [], marker="^", color=STYLES[COMBINED][1], markerfacecolor="white",
                          linestyle="none", label="No paired decidua"))
    fig.legend(handles=handles, loc="outside lower center", ncols=4, frameon=False, fontsize=8.8)
    fig.savefig(OUTPUT / "scores_by_delivery_age.png", dpi=180, facecolor="white")
    plt.close(fig)


def main():
    verify_source()
    audit = pd.read_csv(AUDIT, index_col="donor_id")
    reference = pd.read_csv(STEP2 / "donor_scores.csv", index_col="donor_id")
    donors = audit.loc[audit.group.isin(STYLES)]
    scores = donors[["group", "ga_delivery_weeks", "complete_pair", "sequencing_round", "qc_flags"]].copy()
    for tissue, count_file in (("vt", VT_COUNTS), ("db", DB_COUNTS)):
        scores = scores.join(project_tissue(donors, reference, tissue, count_file))
    scores[PAIRED] = scores[DB_SCORE] - scores[VT_SCORE]
    scores[SELECTIVITY] = scores[DB_SCORE] - scores.db_broad_myeloid
    np.testing.assert_allclose(scores.loc[reference.index, METRICS], reference[METRICS])
    same_round = scores.loc[scores.sequencing_round.eq("other_round")]
    paired = same_round.loc[same_round.complete_pair.eq(1)]
    combined = paired.loc[paired.group.eq(COMBINED)]
    fgr = paired.loc[paired.group.eq(ISOLATED)]
    eope = paired.loc[paired.group.eq(EOPE_ONLY)]
    a_fgr, b_fgr, fgr_window = age_restricted(combined, fgr)
    a_eope, b_eope, eope_window = age_restricted(combined, eope)
    summary = pd.concat([
        compare(combined, fgr, METRICS, "paired_same_round"),
        compare(a_fgr, b_fgr, METRICS, "paired_shared_age_range", fgr_window),
        compare(combined, fgr.drop("S15"), METRICS, "paired_without_S15"),
        compare(same_round.loc[same_round.group.eq(COMBINED)], fgr, [VT_SCORE], "all_villi_same_round"),
        compare(combined, eope, METRICS, "paired_same_round"),
        compare(a_eope, b_eope, METRICS, "paired_shared_age_range", eope_window),
    ], ignore_index=True)
    OUTPUT.mkdir(parents=True, exist_ok=True)
    scores.to_csv(OUTPUT / "bulk_donor_scores.csv")
    summary.to_csv(OUTPUT / "bulk_contrasts.csv", index=False)
    plot_scores(same_round)
    print(summary[["comparison", "reference", "n_numerator", "n_reference",
                   "metric", "estimate", "ci95_low", "ci95_high"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()

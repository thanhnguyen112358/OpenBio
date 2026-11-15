"""Describe the RNA panel in early PE without FGR and early reference controls."""

import numpy as np
import pandas as pd

from step01_audit_cohort import PANEL, PROJECT, SOURCE, verify_source
from step02_deseq2_scores import score_panel


PB_METADATA = "extension/pseudobulk/track01_pseudobulk_metadata.csv"
PB_COUNTS = "extension/pseudobulk/track01_pseudobulk_counts.csv"
OUTPUT = PROJECT / "results" / "followup_fgr_context" / "pseudobulk"
CLASSES = ("macrophage", "hofbauer_cell")
MIN_CELLS = 10
METRICS = ("c1q", "trem2", "c1q_trem2", "broad_myeloid", "c1q_trem2_minus_myeloid")


def load_donor_scores():
    verify_source(PB_METADATA, PB_COUNTS)
    metadata = pd.read_csv(SOURCE / PB_METADATA, index_col="pseudobulk_id")
    within_scope = (
        metadata.clinical_group.isin(["PE_early", "control_early"])
        & metadata.fgr_status.eq("no_FGR")
        & metadata.cell_class.isin(CLASSES)
        & metadata.n_cells.ge(MIN_CELLS)
    )
    donors = metadata.loc[within_scope, [
        "donor_id", "cell_class", "clinical_group", "fgr_status", "ga_delivery_weeks",
        "delivery_mode", "n_cells", "total_counts", "qc_flags",
    ]].copy()
    if donors.duplicated(["donor_id", "cell_class"]).any():
        raise ValueError("Expected one eligible library per donor and cell class")
    if donors.ga_delivery_weeks.isna().any():
        raise ValueError("Delivery age is required for the observed range check")

    raw = pd.read_csv(SOURCE / PB_COUNTS, index_col="gene", usecols=["gene", *donors.index])
    sums = raw.sum(axis=0).reindex(donors.index)
    np.testing.assert_array_equal(sums, donors.total_counts)
    if sums.le(0).any():
        raise ValueError("Pseudobulk libraries must have positive total counts")
    logged = np.log2(1 + raw.loc[list(PANEL), donors.index].T.div(sums, axis=0) * 1e6)
    panel = score_panel(logged)
    panel["c1q_trem2_minus_myeloid"] = panel.c1q_trem2 - panel.broad_myeloid
    return donors.join(logged.add_prefix("log2_cpm_plus1_")).join(panel)


def summarize_comparisons(donors):
    donors = donors.copy()
    donors["in_initial_delivery_age_intersection"] = False
    summaries = []
    for cell_class, by_class in donors.groupby("cell_class"):
        eope = by_class.loc[by_class.clinical_group.eq("PE_early")]
        controls = by_class.loc[by_class.clinical_group.eq("control_early")]
        if eope.empty or controls.empty:
            raise ValueError(f"Both phenotypes are required for {cell_class}")
        lower = max(eope.ga_delivery_weeks.min(), controls.ga_delivery_weeks.min())
        upper = min(eope.ga_delivery_weeks.max(), controls.ga_delivery_weeks.max())
        overlap = by_class.ga_delivery_weeks.between(lower, upper)
        donors.loc[by_class.index, "in_initial_delivery_age_intersection"] = overlap
        for subset, rows in (("all_eligible", by_class), ("initial_age_intersection", by_class.loc[overlap])):
            case_rows = rows.loc[rows.clinical_group.eq("PE_early")]
            control_rows = rows.loc[rows.clinical_group.eq("control_early")]
            effect = case_rows[list(METRICS)].mean() - control_rows[list(METRICS)].mean()
            summaries.append({
                "cell_class": cell_class,
                "cohort": subset,
                "initial_overlap_start_weeks": lower,
                "initial_overlap_end_weeks": upper,
                "eope_only_n": len(case_rows),
                "early_control_n": len(control_rows),
                "eope_only_ages_weeks": " ".join(f"{age:.3f}" for age in sorted(case_rows.ga_delivery_weeks)),
                "early_control_ages_weeks": " ".join(
                    f"{age:.3f}" for age in sorted(control_rows.ga_delivery_weeks)
                ),
                **effect.add_prefix("eope_minus_control_").to_dict(),
            })
    return donors, pd.DataFrame(summaries)


def main():
    donors, comparisons = summarize_comparisons(load_donor_scores())
    OUTPUT.mkdir(parents=True, exist_ok=True)
    donors.sort_values(["cell_class", "clinical_group", "ga_delivery_weeks"]).to_csv(
        OUTPUT / "donor_scores.csv"
    )
    comparisons.to_csv(OUTPUT / "descriptive_contrasts.csv", index=False)
    print(comparisons.to_string(index=False))


if __name__ == "__main__":
    main()

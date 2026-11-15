"""Validate the exported Python analysis and optionally compare it with the R archive."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import false_discovery_control


CORE = ("openbio_primary", "openbio_overlap", "awamleh_primary")
FIGURES = ("age_support", "enrichment", "headline", "atlas", "pca", "volcano")
TARGETS = {"openbio_primary": 0.273, "openbio_overlap": 0.264, "awamleh_primary": 0.145,
           "awamleh_eope_vs_fgr": 0.168, "awamleh_fgr_vs_control": 0.146}
BMP = "REACTOME_SIGNALING_BY_BMP"


def verify(out, reference=None):
    files = [
        "run_manifest.tsv", "bmp_donor_scores.csv", "bmp_effects.csv", "bmp_leave_one_out.csv",
        "bmp_gene_manifest.csv", "bmp_gene_effects.tsv.gz", "cohort_manifest.csv", "cohort_support.csv",
        "comparison_membership.csv", "rank_fit_membership.csv", "rank_universe.csv", "rank_plot.csv",
        "branch_enrichment.csv", "shared_leading_edge.csv", "pca_scores.csv", "atlas_donors.csv",
        "atlas_summary.csv", "bulk_marker_effects.csv", "browser_lookups.csv",
        "paired_gene_coverage.csv", "paired_donor_scores.csv", "paired_tissue_effects.csv",
        *(f"figures/{name}.png" for name in FIGURES),
    ]
    missing = [name for name in files if not (out / name).is_file()]
    if missing:
        raise ValueError(f"Missing workflow outputs: {', '.join(missing)}")
    scores = pd.read_csv(out / "bmp_donor_scores.csv")
    membership = pd.read_csv(out / "comparison_membership.csv")
    genes = pd.read_csv(out / "bmp_gene_manifest.csv")
    effects = pd.read_csv(out / "bmp_effects.csv")
    ranks = pd.read_csv(out / "rank_universe.csv")
    enrichment = pd.read_csv(out / "branch_enrichment.csv")
    assert not scores.duplicated(["cohort", "donor_id"]).any()
    assert scores.groupby("cohort").size().to_dict() == {"Awamleh": 79, "OpenBio": 30}
    assert len(genes) == 31
    assert set(genes.loc[~genes.in_score, "source_symbol"]) == {"AMHR2", "BMP10", "CER1", "GDF2"}
    assert not ranks.symbol.duplicated().any()
    assert (enrichment.groupby("fit").size() == 3).all() and len(enrichment) == 33
    assert (enrichment.rank_size == len(ranks)).all()
    assert set(genes.loc[genes.in_rank_universe, "source_symbol"]) == set(ranks.loc[ranks.in_bmp, "symbol"])
    for _, branch in enrichment.groupby("fit"):
        np.testing.assert_allclose(branch.fdr, false_discovery_control(branch.pvalue, method="bh"))
    primary = effects[effects.score == "BMP27"].set_index("comparison")
    for fit, expected in TARGETS.items():
        assert round(primary.at[fit, "effect"], 3) == expected
        assert primary.at[fit, "loo_min"] > 0
    age = membership[membership.comparison == "openbio_overlap"]
    assert set(age.donor_id) == {"S04", "S08", "S11", "S12", "S14", "S15", "S17"}
    assert (age.ga_delivery_weeks.min(), age.ga_delivery_weeks.max()) == (26.714, 30.429)
    paired = pd.read_csv(out / "paired_tissue_effects.csv").set_index("comparison")
    np.testing.assert_allclose(paired.at["decidua_minus_villi", "effect"],
                               paired.at["decidua", "effect"] - paired.at["villi", "effect"])
    if (out / "piekos_effects.csv").is_file():
        piekos = pd.read_csv(out / "piekos_effects.csv")
        row = piekos[(piekos.score == "BMP27") &
                     (piekos.comparison == "overlap_age_adjusted_exact_ID_only")].iloc[0]
        assert (row.n_combined, row.n_fgr, round(row.effect, 3)) == (15, 16, 0.141)
    if reference is not None:
        compare(out, reference)
    print(f"Verified 109 donor scores, 11 rank fits and six figures under {out}.", flush=True)


def compare(out, reference):
    for file, keys, columns in (
        ("bmp_donor_scores.csv", ["cohort", "donor_id"], ["size_factor", "bmp27", "bmp25_without_inhibins"]),
        ("bmp_effects.csv", ["comparison", "score"], ["effect", "ci_low", "ci_high", "loo_min", "loo_max"]),
        ("paired_tissue_effects.csv", ["comparison"], ["effect", "ci_low", "ci_high"]),
        ("atlas_summary.csv", ["cell_class", "symbol"], ["median_cpm", "fraction_ge10"]),
        ("piekos_effects.csv", ["score", "comparison"], ["effect", "ci_low", "ci_high"]),
    ):
        if file == "piekos_effects.csv" and not (out / file).exists():
            continue
        py = pd.read_csv(out / file).set_index(keys).sort_index()
        r = pd.read_csv(reference / file).set_index(keys).sort_index()
        # Both implementations honor the same candidate atlas genes, plus any
        # additional genes appearing in that implementation's shared leading edge.
        if file == "atlas_summary.csv":
            common = py.index.intersection(r.index)
            assert {"BMPR2", "BMPR1B", "GREM2", "SMAD9"} <= set(common.get_level_values("symbol"))
            py, r = py.loc[common], r.loc[common]
        assert py.index.equals(r.index), file
        np.testing.assert_allclose(py[columns], r[columns], atol=1e-10, rtol=1e-9, err_msg=file)
    py = pd.read_csv(out / "branch_enrichment.csv")
    r = pd.read_csv(reference / "branch_enrichment.csv")
    py = py[py.fit.isin(CORE) & (py.pathway == BMP)]
    r = r[r.fit.isin(CORE) & (r.pathway == BMP)]
    columns = ["fit", "rank_size", "size", "NES", "fdr", "leading_edge"]
    py[columns].merge(r[columns], on="fit", suffixes=("_python", "_r"),
                      validate="one_to_one").to_csv(out / "python_vs_r_enrichment.csv", index=False)
    print("Fixed donor scores, intervals, paired estimates and matching atlas summaries agree with the R archive.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--reference", type=Path)
    args = parser.parse_args()
    verify(args.output, args.reference)

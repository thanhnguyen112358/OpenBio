"""Score the paired RNA panel and fit gene-level EOPE versus FGR contrasts."""

import argparse
import json
import platform
from importlib.metadata import version

import numpy as np
import pandas as pd
from joblib import parallel_config
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

from step01_audit_cohort import DB_COUNTS, PANEL, PROJECT, SOURCE, VT_COUNTS, verify_source


AUDIT = PROJECT / "results" / "step01_cohort_audit" / "donor_audit.csv"
OUTPUT = PROJECT / "results" / "step02_deseq2"
NGUYEN_REVISION = "7826b8d4a146a8d888269b51fd7fbf7039de502f"
COHORTS = (
    "inferred_same_round",
    "same_round_in_full_cohort_age_overlap",
    "same_round_in_within_round_age_overlap",
    "all_eope_fgr_pairs",
)
C1Q = ["C1QA", "C1QB", "C1QC"]
MYELOID = ["CD68", "LST1", "CSF1R"]


def dataset(counts, metadata, design, inference):
    return DeseqDataSet(
        counts=counts,
        metadata=metadata,
        design=design,
        refit_cooks=False,
        inference=inference,
        quiet=True,
    )


def normalize_tissue(count_file, donors, tissue, inference):
    sample_ids = donors[f"{tissue.lower()}_sample_id"]
    raw = pd.read_csv(SOURCE / count_file, index_col="gene").loc[:, sample_ids].T
    raw.index.name = "sample_id"
    metadata = donors.reset_index().set_index(sample_ids.name)[["donor_id", "group"]]
    # Adapted from Nguyen's donor_expression(). Keep full precision for scoring.
    dds = dataset(raw, metadata, "~ 1", inference)
    dds.fit_size_factors()
    factors = dds.obs["size_factors"].rename("size_factor")
    normalized = pd.DataFrame(
        dds.layers["normed_counts"], index=raw.index, columns=raw.columns
    ).loc[:, list(PANEL)]
    logged = np.log2(normalized + 1)
    score = pd.DataFrame(index=raw.index)
    score["c1q"] = logged[C1Q].mean(axis=1)
    score["trem2"] = logged["TREM2"]
    score["c1q_trem2"] = (score["c1q"] + score["trem2"]) / 2
    score["broad_myeloid"] = logged[MYELOID].mean(axis=1)
    score["size_factor"] = factors
    score.index = donors.index

    panel = pd.concat(
        {
            "raw_count": raw.loc[:, list(PANEL)].stack(),
            "normalized_count": normalized.stack(),
            "log2_normalized_plus1": logged.stack(),
        },
        axis=1,
    ).rename_axis(index=["sample_id", "gene"]).reset_index()
    panel.insert(1, "donor_id", panel["sample_id"].map(metadata["donor_id"]))
    panel.insert(2, "tissue", tissue)
    panel["group"] = panel["sample_id"].map(metadata["group"])
    return raw, factors, score, panel


def fit_gene_panel(counts, metadata, design, factors, contrast, inference):
    dds = dataset(counts, metadata, design, inference)
    design_matrix = dds.obsm["design_matrix"]
    if np.linalg.matrix_rank(design_matrix) != design_matrix.shape[1]:
        raise ValueError(f"Rank-deficient DESeq2 design: {design}")
    if design_matrix.shape[0] <= design_matrix.shape[1]:
        raise ValueError(f"No residual degrees of freedom: {design}")

    # Run the public stages of dds.deseq2() after supplying fixed size factors.
    # Calling the wrapper here would re-estimate them jointly or within a subset.
    dds.obs["size_factors"] = factors.loc[counts.index]
    dds.layers["normed_counts"] = counts.div(dds.obs["size_factors"], axis=0).to_numpy()
    dds.var["_normed_means"] = dds.layers["normed_counts"].mean(axis=0)
    dds.fit_genewise_dispersions()
    dds.fit_dispersion_trend()
    dds.fit_dispersion_prior()
    dds.fit_MAP_dispersions()
    dds.fit_LFC()
    dds.calculate_cooks()
    dds.cooks_outlier()

    if isinstance(contrast, str):
        contrast = (design_matrix.columns == contrast).astype(float)
        if contrast.sum() != 1:
            raise ValueError("Paired interaction coefficient missing from fitted design")
    stats = DeseqStats(dds, contrast=contrast, inference=inference, quiet=True)
    stats.summary()
    result = stats.results_df.loc[list(PANEL)].copy().rename_axis("gene")
    result["ci95_low"] = result["log2FoldChange"] - 1.96 * result["lfcSE"]
    result["ci95_high"] = result["log2FoldChange"] + 1.96 * result["lfcSE"]
    result["lfc_converged"] = dds.var.loc[result.index, "_LFC_converged"].astype(bool)
    result["cooks_outlier"] = dds.var.loc[result.index, "_pvalue_cooks_outlier"].astype(bool)
    result["dispersion_trend"] = dds.uns["disp_function_type"]
    return result.reset_index()


def main(n_cpus):
    verify_source()
    audit = pd.read_csv(AUDIT, index_col="donor_id")
    donors = audit.loc[audit["in_all_eope_fgr_pairs"].eq(1)].copy()
    if not donors.index.is_unique or not donors["complete_pair"].eq(1).all():
        raise ValueError("Step 1 audit must contain one complete pair per target donor")
    if donors["group"].value_counts().to_dict() != {"EOPE": 8, "FGR": 5}:
        raise ValueError("Target cohort differs from the prespecified 8 EOPE and 5 FGR")
    inference = DefaultInference(n_cpus=n_cpus)

    vt, vt_factors, vt_score, vt_panel = normalize_tissue(VT_COUNTS, donors, "VT", inference)
    db, db_factors, db_score, db_panel = normalize_tissue(DB_COUNTS, donors, "DB", inference)
    membership = [f"in_{cohort}" for cohort in COHORTS]
    scores = donors[["group", "ga_delivery_weeks", "sequencing_round", "qc_flags", *membership]]
    scores = scores.join(vt_score.add_prefix("vt_")).join(db_score.add_prefix("db_"))
    scores["db_minus_vt_c1q_trem2"] = scores["db_c1q_trem2"] - scores["vt_c1q_trem2"]
    scores["db_c1q_trem2_minus_broad_myeloid"] = (
        scores["db_c1q_trem2"] - scores["db_broad_myeloid"]
    )
    scores.to_csv(OUTPUT / "donor_scores.csv")
    pd.concat([vt_panel, db_panel], ignore_index=True).to_csv(
        OUTPUT / "normalized_panel.csv", index=False
    )

    effects = []
    fit_info = {}
    for cohort in COHORTS:
        selected = donors.loc[donors[f"in_{cohort}"].eq(1)]
        group_sizes = selected["group"].value_counts().to_dict()
        print(f"{cohort}: {group_sizes['EOPE']} EOPE vs {group_sizes['FGR']} FGR", flush=True)
        by_tissue = []
        for tissue, raw, factors in (("VT", vt, vt_factors), ("DB", db, db_factors)):
            print(f"  {tissue}: fitting EOPE minus FGR", flush=True)
            samples = selected[f"{tissue.lower()}_sample_id"]
            metadata = selected.reset_index().set_index(samples.name)[["donor_id", "group"]]
            metadata["is_decidua"] = int(tissue == "DB")
            counts = raw.loc[samples]
            result = fit_gene_panel(
                counts, metadata, "~ group", factors, ["group", "EOPE", "FGR"], inference
            )
            effects.append(result.assign(cohort=cohort, contrast=f"{tissue}_EOPE_minus_FGR"))
            by_tissue.append((counts, metadata))

        paired_counts = pd.concat([counts for counts, _ in by_tissue])
        paired_meta = pd.concat([metadata for _, metadata in by_tissue])
        paired_meta["eope_in_decidua"] = (
            paired_meta["is_decidua"] * paired_meta["group"].eq("EOPE")
        )
        paired_design = "~ donor_id + is_decidua + eope_in_decidua"
        print("  Fitting the donor-blocked tissue interaction", flush=True)
        interaction = fit_gene_panel(
            paired_counts,
            paired_meta,
            paired_design,
            pd.concat([vt_factors, db_factors]),
            "eope_in_decidua",
            inference,
        )
        effects.append(interaction.assign(cohort=cohort, contrast="paired_tissue_interaction"))
        fit_info[cohort] = {"group_sizes": group_sizes, "donor_ids": selected.index.tolist()}

    columns = ["cohort", "contrast", "gene"]
    pd.concat(effects, ignore_index=True).set_index(columns).to_csv(OUTPUT / "gene_effects.csv")
    run_info = {
        "nguyen_revision": NGUYEN_REVISION,
        "python": platform.python_version(),
        "packages": {
            name: version(name)
            for name in ("pydeseq2", "pandas", "numpy", "scipy", "scikit-learn", "anndata")
        },
        "size_factors": "Median-of-ratios, 13 complete pairs per tissue, fixed across fits",
        "normalization_donors": donors.index.tolist(),
        "n_cpus": n_cpus,
        "tissue_design": "~ group",
        "paired_design": paired_design,
        "cohorts": fit_info,
    }
    (OUTPUT / "run_info.json").write_text(json.dumps(run_info, indent=2) + "\n", encoding="utf-8")
    print(f"Saved scores for {len(scores)} donors and {sum(len(table) for table in effects)} gene effects.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--n-cpus", type=int, default=6, help="PyDESeq2 worker processes")
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    with parallel_config(temp_folder=str(OUTPUT)):
        main(args.n_cpus)

"""Compute paired-tissue and Piekos sensitivities from the released counts."""

import csv
from collections import Counter
from io import TextIOWrapper
from pathlib import Path
from zipfile import ZipFile

import numpy as np
import pandas as pd
from scipy import stats


def median_ratio_size_factors(counts: pd.DataFrame) -> pd.Series:
    """Use the DESeq2 median of log ratios over genes positive in every sample."""
    values = counts.to_numpy(dtype=float)
    if not (np.isfinite(values).all() and (values >= 0).all() and
            np.equal(values, np.floor(values)).all()):
        raise ValueError("Counts must be finite nonnegative integers")
    positive = (values > 0).all(axis=1)
    if not positive.any():
        raise ValueError("No genes are positive in every sample for median-ratio normalization")
    logs = np.log(values[positive])
    factors = np.exp(np.median(logs - logs.mean(axis=1, keepdims=True), axis=0))
    return pd.Series(factors, index=counts.columns, name="size_factor")


def welch_effect(case, reference) -> dict[str, float]:
    """Difference in means and a two-sided, approximate Welch 95% interval."""
    case, reference = np.asarray(case), np.asarray(reference)
    interval = stats.ttest_ind(case, reference, equal_var=False).confidence_interval()
    return {
        "effect": case.mean() - reference.mean(),
        "ci_low": interval.low,
        "ci_high": interval.high,
    }


def paired(
    vt_counts: pd.DataFrame,
    db_counts: pd.DataFrame,
    villi: pd.DataFrame,
    decidua: pd.DataFrame,
    fixed_genes: pd.DataFrame,
    out: Path,
) -> None:
    """Normalize each tissue across its full release before selecting paired donors."""
    genes = fixed_genes["local_gene"].tolist()
    coverage = db_counts.reindex(genes).ge(10).sum(axis=1).where(pd.Index(genes).isin(db_counts.index))
    covered = coverage.ge(5)
    paired_genes = coverage.index[covered]
    if paired_genes.empty or not paired_genes.is_unique:
        raise ValueError("Expected a nonempty, unique panel of adequately covered paired genes")
    gene_coverage = pd.DataFrame({
        "symbol": fixed_genes["source_symbol"].to_numpy(),
        "local_gene": genes,
        "decidua_libraries_ge10": coverage.to_numpy(),
        "in_paired_score": covered.to_numpy(),
    })

    villi = villi.set_index("sample_id", verify_integrity=True).loc[vt_counts.columns]
    decidua = decidua.set_index("donor_id", verify_integrity=True)
    rows = villi.loc[
        villi["group"].isin(["EOPE+FGR", "FGR"]) & villi["donor_id"].isin(decidua.index),
        ["donor_id", "group", "ga_delivery_weeks"],
    ].reset_index(names="vt_sample_id")
    db_rows = decidua.loc[rows["donor_id"]]
    if not rows["donor_id"].is_unique or not np.array_equal(rows["group"], db_rows["group"]):
        raise ValueError("Paired donors must be unique and have matching tissue diagnoses")
    rows["db_sample_id"] = db_rows["sample_id"].to_numpy()

    for tissue, library, counts in (
        ("villi", "vt_sample_id", vt_counts), ("decidua", "db_sample_id", db_counts),
    ):
        factors = median_ratio_size_factors(counts)
        normalized = counts.loc[paired_genes, rows[library]].divide(factors.loc[rows[library]], axis=1)
        rows[tissue] = np.log2(normalized + 1).mean(axis=0).to_numpy()
    rows["decidua_minus_villi"] = rows["decidua"] - rows["villi"]
    rows = rows[["donor_id", "group", "ga_delivery_weeks", "vt_sample_id", "db_sample_id",
                 "villi", "decidua", "decidua_minus_villi"]]

    cases = rows["group"].eq("EOPE+FGR")
    n_a, n_b = int(cases.sum()), int((~cases).sum())
    if (n_a, n_b) != (4, 5):
        raise ValueError(f"Expected 4 combined and 5 FGR paired donors, found {n_a} and {n_b}")
    effects = []
    for comparison in ("villi", "decidua", "decidua_minus_villi"):
        effects.append({
            "comparison": comparison, "n_genes": len(paired_genes), "n_a": n_a, "n_b": n_b,
            **welch_effect(rows.loc[cases, comparison], rows.loc[~cases, comparison]),
        })
    gene_coverage.to_csv(out / "paired_gene_coverage.csv", index=False)
    rows.to_csv(out / "paired_donor_scores.csv", index=False)
    pd.DataFrame(effects).to_csv(out / "paired_tissue_effects.csv", index=False)


def _read_piekos_counts(count_archive: Path, genes: list[str]) -> pd.DataFrame:
    with ZipFile(count_archive) as archive:
        with archive.open("transcripts-final.csv") as raw_header:
            columns = next(csv.reader(TextIOWrapper(raw_header, encoding="utf-8")))
        names = Counter(columns)
        if any(names[name] != 1 for name in ["Patient-ID", *genes]):
            raise ValueError("Expected a unique Piekos transcript identifier and each BMP gene")
        with archive.open("transcripts-final.csv") as transcript:
            table = pd.read_csv(transcript)
    patient_ids = table["Patient-ID"]
    if patient_ids.isna().any() or not patient_ids.is_unique:
        raise ValueError("Piekos transcript sample IDs must be unique and nonmissing")
    return table.iloc[:, [name not in {"", "Patient-ID"} for name in columns]].set_axis(
        patient_ids, axis=0
    ).T


def _hc3(cohort: pd.DataFrame, score: str, comparison: str, adjust: bool) -> dict:
    cases = cohort["combined"].to_numpy()
    variables = [np.ones(len(cohort)), cases]
    if adjust:
        variables.append(cohort["WksGest"].to_numpy())
    design = np.column_stack(variables)
    outcome = cohort[score].to_numpy()
    orthogonal, triangular = np.linalg.qr(design)
    beta = np.linalg.solve(triangular, orthogonal.T @ outcome)
    leverage = np.square(orthogonal).sum(axis=1)
    residuals = (outcome - design @ beta) / (1 - leverage)
    linear_weights = np.linalg.solve(triangular, orthogonal.T)
    se = np.linalg.norm(linear_weights[1] * residuals)
    halfwidth = stats.t.ppf(0.975, len(outcome) - design.shape[1]) * se
    return {
        "score": score, "comparison": comparison,
        "n_combined": int(cases.sum()), "n_fgr": int(len(cases) - cases.sum()),
        "effect": beta[1], "ci_low": beta[1] - halfwidth, "ci_high": beta[1] + halfwidth,
    }


def piekos(
    count_archive: Path,
    clinical_path: Path,
    gene_manifest: pd.DataFrame,
    out: Path,
) -> None:
    """Recreate severe-PE sensitivity, reporting the exact-linkage analysis too."""
    source_genes = gene_manifest["source_symbol"].tolist()
    counts = _read_piekos_counts(count_archive, source_genes)
    normalized = counts.loc[source_genes].divide(median_ratio_size_factors(counts), axis=1)

    identifier = "Study_ID/ID"
    clinical = pd.read_csv(clinical_path, usecols=[identifier, "Condition", "WksGest"])
    clinical = clinical.dropna(subset=[identifier]).set_index(identifier, verify_integrity=True)
    clinical["WksGest"] = pd.to_numeric(clinical["WksGest"], errors="coerce")
    rna_ids = counts.columns.to_numpy()
    if set(rna_ids).difference(clinical.index) != {"MJ-0323"} or "MJ-0323-1" not in clinical.index:
        raise ValueError("Piekos clinical and transcript sample IDs differ from the audited release")
    # Preserve the one documented suffix match and also report exact-only estimates.
    clinical_ids = pd.Index(rna_ids).where(rna_ids != "MJ-0323", "MJ-0323-1")
    if not clinical_ids.is_unique:
        raise ValueError("The Piekos suffix match would reuse a clinical record")
    cohort = clinical.loc[clinical_ids].reset_index(names=identifier)
    cohort["rna_id"] = rna_ids
    cohort["link"] = np.where(
        cohort[identifier].eq(cohort["rna_id"]), "exact", "unique clinical ID with terminal -1"
    )
    score_genes = gene_manifest.loc[gene_manifest["in_score"], "source_symbol"]
    cohort["BMP27"] = np.log2(normalized.loc[score_genes] + 1).mean(axis=0).to_numpy()
    cohort["BMP31_pilot"] = np.log2(normalized + 0.5).mean(axis=0).to_numpy()
    cohort = cohort.loc[cohort["Condition"].isin(["FGR+PE Severe", "FGR"])].copy()
    cohort["combined"] = cohort["Condition"].eq("FGR+PE Severe").astype(int)

    ranges = cohort.groupby("Condition")["WksGest"].agg(["min", "max"])
    cohort["in_age_support"] = cohort["WksGest"].between(ranges["min"].max(), ranges["max"].min())
    overlap = cohort.loc[cohort["in_age_support"]]
    exact = cohort.loc[cohort["link"].eq("exact")]
    overlap_exact = overlap.loc[overlap["link"].eq("exact")]
    if cohort["combined"].value_counts().to_dict() != {0: 36, 1: 25} or (
        overlap_exact["combined"].value_counts().to_dict() != {0: 16, 1: 15}
    ):
        raise ValueError("Piekos sensitivity donor counts differ from the audited release")
    comparisons = (
        ("full_unadjusted", cohort, False),
        ("overlap_unadjusted", overlap, False),
        ("overlap_age_adjusted", overlap, True),
        ("full_unadjusted_exact_ID_only", exact, False),
        ("overlap_age_adjusted_exact_ID_only", overlap_exact, True),
    )
    estimates = [
        _hc3(group, score, name, adjust)
        for score in ("BMP27", "BMP31_pilot")
        for name, group, adjust in comparisons
    ]
    cohort.to_csv(out / "piekos_donors.csv", index=False)
    pd.DataFrame(estimates).to_csv(out / "piekos_effects.csv", index=False)

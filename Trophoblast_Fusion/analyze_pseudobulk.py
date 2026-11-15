"""Fit candidate host transcript contrasts in donor-level placental pseudobulk."""

import argparse
from pathlib import Path
from hashlib import sha256

from importlib.metadata import version

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import t

FILES = {
    "track01_pseudobulk_counts.csv": "81817e4977e1088c8061c6fb3b13e77257b72493daa1d6c5e7ee51febeae7c15",
    "track01_pseudobulk_metadata.csv": "7ada903628e8f7271b0c3204af747b8219d2b323309f9c1751e55c077b2a0330",
}
CLASSES = ("villous_cytotrophoblast", "syncytiotrophoblast", "macrophage", "hofbauer_cell")
PANEL = ("IFITM2", "IFITM3", "GBP5")

FGR = "PEplusFGR"
PE = "PEonly"


def load(source):
    location = source / "extension" / "pseudobulk"
    for name, digest in FILES.items():
        if sha256((location / name).read_bytes()).hexdigest() != digest:
            raise ValueError(f"Digest differs for {name}")
    meta = pd.read_csv(location / "track01_pseudobulk_metadata.csv",
                       index_col="pseudobulk_id")
    meta = meta.loc[meta["clinical_group"].eq("PE_early") & meta["cell_class"].isin(CLASSES)].copy()
    meta["phenotype"] = meta["fgr_status"].map({"FGR": FGR, "no_FGR": PE})
    if not meta.index.is_unique or meta["phenotype"].isna().any():
        raise ValueError("Invalid metadata")
    counts = pd.read_csv(location / "track01_pseudobulk_counts.csv", index_col="gene",
                         usecols=["gene", *meta.index]).T.loc[meta.index]
    if not counts.columns.is_unique or not (counts.sum(axis=1) == meta["total_counts"]).all():
        raise ValueError("Invalid counts or mismatch with metadata")
    return counts, meta


def cohort(meta, cell_class, minimum_cells=20, overlap=False, cesareans=False):
    sample = meta.loc[meta["cell_class"].eq(cell_class) & meta["n_cells"].ge(minimum_cells)].copy()
    if cesareans:
        sample = sample.loc[sample["delivery_mode"].eq("c_section")].copy()
    if sample["donor_id"].duplicated().any():
        raise ValueError("More than one profile per donor in the cohort")
    if overlap:
        ages = sample.groupby("phenotype")["ga_delivery_weeks"]
        sample = sample.loc[sample["ga_delivery_weeks"].between(ages.min().max(), ages.max().min())]
    return sample


def replicate_counts(sample):
    return {"n_FGR": int(sample["phenotype"].eq(FGR).sum()), "n_PE": int(sample["phenotype"].eq(PE).sum())}


def fit(counts, sample, design):
    sample = sample.assign(ga_centered=sample["ga_delivery_weeks"] - sample["ga_delivery_weeks"].mean())
    selected = counts.loc[sample.index]
    selected = selected.loc[:, selected.ge(10).sum(axis=0).ge(3)]
    model = DeseqDataSet(counts=selected, metadata=sample, design=design,
                         size_factors_fit_type="poscounts", n_cpus=1,
                         refit_cooks=False, quiet=True)
    if np.linalg.matrix_rank(model.obsm["design_matrix"]) != model.obsm["design_matrix"].shape[1]:
        raise ValueError("Rank deficient")
    model.deseq2()
    stats = DeseqStats(model, contrast=["phenotype", FGR, PE], quiet=True, n_cpus=1,
                       independent_filter=False, cooks_filter=False)
    stats.summary()
    results = stats.results_df.reindex(PANEL).reset_index(names="gene")
    results["ci_low"] = results["log2FoldChange"] - 1.96 * results["lfcSE"]
    results["ci_high"] = results["log2FoldChange"] + 1.96 * results["lfcSE"]
    results["dispersion_trend"] = model.uns["disp_function_type"]
    cooks = pd.DataFrame(model.layers["cooks"], index=model.obs_names, columns=model.var_names)
    results["max_cooks"] = results["gene"].map(cooks.max(axis=0))
    results["lfc_converged"] = results["gene"].map(model.var["_LFC_converged"])
    size_factors = model.obs["size_factors"]
    return results, size_factors


def ols(values, metadata, age=False):
    phenotype = metadata["phenotype"].eq(FGR).astype(float).to_numpy()
    columns = [np.ones(len(metadata)), phenotype]
    if age:
        columns.append(metadata["ga_delivery_weeks"].to_numpy() - metadata["ga_delivery_weeks"].mean())
    design = np.column_stack(columns)
    return linear_interval(values.loc[metadata.index].to_numpy(), design)


def linear_interval(values, design):
    """Return the coefficient of the first non-intercept term and its t interval."""
    degrees_freedom = len(values) - design.shape[1]
    if degrees_freedom <= 0 or np.linalg.matrix_rank(design) < design.shape[1]:
        return np.nan, np.nan, np.nan
    beta = np.linalg.lstsq(design, values, rcond=None)[0]
    residual = values - design @ beta
    se = np.sqrt((residual @ residual / degrees_freedom) * np.linalg.inv(design.T @ design)[1, 1])
    half_width = t.ppf(0.975, degrees_freedom) * se
    return beta[1], beta[1] - half_width, beta[1] + half_width


def vct_maturation(counts, sample, matched, factors):
    """Check GCM1 as an exploratory donor-level differentiation proxy in VCT."""
    targets = ["IFITM2", "GCM1"]
    logged = np.log2(counts.loc[sample.index, targets].div(factors, axis=0) + 1)
    frame = sample.join(logged)
    frame["combined"] = frame["phenotype"].eq(FGR).astype(float)
    for gene in targets:
        frame[f"{gene}_z"] = (frame[gene] - frame[gene].mean()) / frame[gene].std(ddof=1)
    base = ["combined", "ga_delivery_weeks"]
    scenarios = [
        ("full", frame, "GCM1", "IFITM2", base),
        ("full", frame, "GCM1_z", "IFITM2_z", base),
        ("full_group_only", frame, "GCM1_z", "IFITM2_z", ["combined"]),
        ("overlap", frame.loc[matched.index], "GCM1_z", "IFITM2_z", base),
    ]
    scenarios.extend((f"omit_{donor}", frame.drop(donor), "GCM1_z", "IFITM2_z", base)
                     for donor in frame.index)
    rows = []
    for label, data, outcome, predictor, adjustment in scenarios:
        design = np.column_stack([np.ones(len(data)), data[[predictor, *adjustment]].to_numpy()])
        estimate, low, high = linear_interval(data[outcome].to_numpy(), design)
        rows.append({"scenario": label, "outcome": outcome, "predictor": predictor,
                     "terms": " + ".join(adjustment), "beta": estimate, "ci_low": low, "ci_high": high,
                     "n": len(data), "residual_df": len(data) - design.shape[1]})
    return pd.DataFrame(rows)


def report_counts(counts, sample, cell_class):
    rows = []
    for gene in PANEL:
        for phenotype in (FGR, PE):
            donors = sample.index[sample["phenotype"].eq(phenotype)]
            data = counts.loc[donors, gene]
            rows.append({"cell_class": cell_class, "phenotype": phenotype, "gene": gene,
                         "n": len(data), "n_zero": int(data.eq(0).sum()),
                         "n_ge10": int(data.ge(10).sum()), "min": int(data.min()),
                         "median": float(data.median()), "max": int(data.max())})
    return rows


def normalized_sensitivity(counts, sample, size_factors, cell_class, tag):
    values = np.log2(counts.loc[sample.index, "IFITM2"] / size_factors.loc[sample.index] + 1)
    rows = []
    for age in (False, True):
        center, low, high = ols(values, sample, age=age)
        row = {"cell_class": cell_class, "scenario": tag, "age_adjusted": age, "LFC_log2_norm_plus1": center,
               "ci_low": low, "ci_high": high, **replicate_counts(sample)}
        if tag == "full_floor20":
            leave_one_out = [(identifier, ols(values, sample.drop(identifier), age=age)[0]) for identifier in sample.index]
            min_donor = min(leave_one_out, key=lambda p: p[1])
            max_donor = max(leave_one_out, key=lambda p: p[1])
            row |= {"loo_min_id": min_donor[0], "loo_min": min_donor[1],
                    "loo_max_id": max_donor[0], "loo_max": max_donor[1]}
        rows.append(row)
    return rows


def main(source, output):
    counts, metadata = load(source)
    output.mkdir(parents=True, exist_ok=True)
    count_fits = []
    coverage = []
    sensitivities = []
    donor_values = []
    metadata_columns = ["donor_id", "phenotype", "cell_class", "n_cells", "total_counts",
                        "ga_delivery_weeks", "delivery_mode", "labor_induction", "qc_flags"]
    for cell_class in CLASSES:
        sample = cohort(metadata, cell_class)
        coverage.extend(report_counts(counts, sample, cell_class))
        scenarios = [
            ("full_floor20", sample, "~ phenotype"),
            ("age_adjusted_floor20", sample, "~ ga_centered + phenotype"),
            ("overlap_floor20", cohort(metadata, cell_class, overlap=True), "~ phenotype"),
        ]
        floor10 = cohort(metadata, cell_class, minimum_cells=10)
        if not sample.index.equals(floor10.index):
            scenarios.append(("full_floor10", floor10, "~ phenotype"))
        if cell_class != CLASSES[0]:
            scenarios = scenarios[:1]
        for tag, cohort_meta, design in scenarios:
            print(cell_class, tag, replicate_counts(cohort_meta), flush=True)
            results, factors = fit(counts, cohort_meta, design)
            results = results.assign(cell_class=cell_class, scenario=tag,
                                     **replicate_counts(cohort_meta))
            count_fits.append(results)
            print(results[["gene", "log2FoldChange", "ci_low", "ci_high", "lfc_converged"]].to_string(index=False), flush=True)
            if tag == "full_floor20":
                sensitivities.extend(normalized_sensitivity(counts, sample, factors,
                                                            cell_class, "full_floor20"))
                matched = cohort(metadata, cell_class, overlap=True)
                if cell_class == CLASSES[0]:
                    maturation = vct_maturation(counts, sample, matched, factors)
                    maturation.to_csv(output / "pseudobulk_vct_ifitm2_gcm1.csv", index=False)
                sensitivities.extend(normalized_sensitivity(counts, matched, factors,
                                                            cell_class, "overlap_floor20_fixed_factors"))
                c_sections = cohort(metadata, cell_class, cesareans=True)
                sensitivities.extend(normalized_sensitivity(counts, c_sections, factors,
                                                            cell_class, "cesareans_floor20_fixed_factors"))
                donor_frame = sample.loc[:, metadata_columns].copy()
                donor_frame["size_factor_poscounts"] = factors
                donor_frame = donor_frame.join(counts.loc[sample.index, list(PANEL)])
                donor_frame["IFITM2_log2_normalized_plus1"] = np.log2(donor_frame["IFITM2"] / factors + 1)
                donor_values.append(donor_frame)
    pd.concat(count_fits, ignore_index=True).to_csv(output / "pseudobulk_count_models.csv", index=False)
    pd.DataFrame(coverage).to_csv(output / "pseudobulk_coverage.csv", index=False)
    pd.DataFrame(sensitivities).to_csv(output / "pseudobulk_log_scale_sensitivity.csv", index=False)
    pd.concat(donor_values).to_csv(output / "pseudobulk_donors.csv", index=True)
    print("PyDESeq2 version", version("pydeseq2"))
    print("Outputs", output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path("openbio-datathon-2026/track-01/data"))
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    main(args.source, args.output)




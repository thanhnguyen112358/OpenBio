"""Cohort audit, fixed donor scores, and placental cell-class expression."""

import csv
import gzip
import numpy as np
import pandas as pd

from extensions import median_ratio_size_factors, welch_effect


CASE = "EOPE+FGR"
FGR = "FGR"
EXTERNAL_LABELS = {"combined": CASE, "pe": "EOPE", "fgr": FGR, "control": "Preterm control"}
BRANCHES = (
    "REACTOME_SIGNALING_BY_BMP",
    "REACTOME_SIGNALING_BY_ACTIVIN",
    "REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX",
)
SCORE_COLUMNS = {"BMP27": "bmp27", "BMP25_without_INHA_INHBA": "bmp25_without_inhibins"}
DONOR_FIELDS = ["cohort", "sample_id", "donor_id", "group", "ga_delivery_weeks"]


def read_counts(path, sep=","):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
        header = next(csv.reader(handle, delimiter=sep))
    if len(set(header[1:])) != len(header[1:]) or not all(header[1:]):
        raise ValueError(f"Duplicated or empty sample IDs in {path}")
    counts = pd.read_csv(path, sep=sep, index_col=0, low_memory=False)
    if not counts.index.is_unique or counts.index.isna().any():
        raise ValueError(f"Duplicated or missing count feature IDs in {path}")
    values = counts.to_numpy()
    if not np.issubdtype(values.dtype, np.number) or not (
        np.isfinite(values).all() and (values >= 0).all() and (values == np.floor(values)).all()
    ):
        raise ValueError(f"Invalid raw counts in {path}")
    return counts.astype(np.int64, copy=False)


def indexed(path, key, sep=","):
    table = pd.read_csv(path, sep=sep)
    if table[key].isna().any() or table[key].duplicated().any():
        raise ValueError(f"Missing or duplicate {key} in {path}")
    return table.set_index(key, drop=False)


def align(counts, samples):
    count_only = counts.columns.difference(samples.index)
    metadata_only = samples.index.difference(counts.columns)
    if len(count_only) or len(metadata_only):
        raise ValueError(f"Count columns and metadata differ. Count-only: {list(count_only)}; "
                         f"metadata-only: {list(metadata_only)}")
    return samples.loc[counts.columns].copy()


def check_groups(frame, expected):
    actual = frame.group.value_counts(dropna=False).to_dict()
    if actual != expected:
        raise ValueError(f"Unexpected clinical groups: {actual}; expected {expected}")


def read_pathways(path):
    pathways = {}
    with path.open(encoding="utf-8", newline="") as handle:
        for key, _, *genes in csv.reader(handle, delimiter="\t"):
            if key in BRANCHES:
                pathways[key] = genes
    if set(pathways) != set(BRANCHES) or len(pathways[BRANCHES[0]]) != 31:
        raise ValueError("Unexpected Reactome branch definitions")
    return {key: pathways[key] for key in BRANCHES}


def audit_tissue(counts, metadata, audit, qc, tissue):
    matched = audit.loc[metadata.donor_id].set_axis(metadata.index)
    pairs = {"group": "group", "sample_id": f"{tissue}_sample_id", "ga_delivery_weeks": "ga_delivery_weeks",
             "pe_status": "pe_status", "pe_onset": "pe_onset", "fgr_status": "fgr_status"}
    if metadata[["donor_id", "group", "ga_delivery_weeks", "fetal_sex"]].isna().any().any():
        raise ValueError(f"Missing required metadata for {tissue}")
    if metadata.donor_id.duplicated().any():
        raise ValueError(f"Duplicated donors within {tissue}")
    for original, audited in pairs.items():
        pd.testing.assert_series_equal(metadata[original], matched[audited], check_names=False, check_dtype=False)
    depth = qc.loc[metadata.sample_id]
    htseq = depth.drop(columns="sample_id").sum(axis=1)
    if not (htseq > 0).all() or not np.array_equal(depth.assigned_to_genes, counts.sum(axis=0)):
        raise ValueError(f"HTSeq library totals do not agree with gene counts in {tissue}")
    if matched.sequencing_round.isna().any():
        raise ValueError(f"Missing inferred sequencing round in {tissue}")
    return metadata.assign(cohort="OpenBio", tissue=tissue, inferred_round=matched.sequencing_round,
                           round_evidence=matched.sequencing_round_evidence,
                           assigned_gene_counts=depth.assigned_to_genes, htseq_total=htseq,
                           assignment_fraction=depth.assigned_to_genes / htseq,
                           labor=np.nan, delivery_mode=np.nan)


def cohorts(data, audit_path, awamleh, out):
    db = data / "extension/decidua"
    counts = {"vt": read_counts(data / "track01_counts_vt.csv"),
              "db": read_counts(db / "track01_counts_db.csv")}
    metadata = {"vt": align(counts["vt"], indexed(data / "track01_metadata_vt.csv", "sample_id")),
                "db": align(counts["db"], indexed(db / "track01_metadata_db.csv", "sample_id"))}
    audit = indexed(audit_path, "donor_id")
    qc = indexed(data / "qc/track01_htseq_summary.csv", "sample_id")
    all_samples = metadata["vt"].index.append(metadata["db"].index)
    all_donors = set(metadata["vt"].donor_id) | set(metadata["db"].donor_id)
    if not all_samples.is_unique or set(all_samples) != set(qc.index) or all_donors != set(audit.index):
        raise ValueError("Organizer counts, donor audit and HTSeq sample IDs do not agree")
    villi = audit_tissue(counts["vt"], metadata["vt"], audit, qc, "vt")
    decidua = audit_tissue(counts["db"], metadata["db"], audit, qc, "db")
    check_groups(villi, {"EOPE": 8, CASE: 7, FGR: 5, "iPTB": 5, "TB": 5})
    check_groups(decidua, {"EOPE": 8, CASE: 4, FGR: 5, "iPTB": 5, "TB": 5})
    external_counts = read_counts(awamleh / "gse114691_gene_counts.tsv.gz", "\t")
    geo = align(external_counts, indexed(awamleh / "gse114691_matched_metadata.tsv", "accession", "\t"))
    check_groups(geo, {"combined": 20, "pe": 20, "fgr": 18, "control": 21})
    external = pd.DataFrame({
        "cohort": "Awamleh", "tissue": "vt", "sample_id": geo.accession, "donor_id": geo.accession,
        "geo_accession": geo.accession, "group": geo.group.map(EXTERNAL_LABELS),
        "geo_diagnosis": geo.geo_diagnosis, "source_column": geo.source_column,
        "ga_delivery_weeks": np.nan, "assigned_gene_counts": external_counts.sum(axis=0),
    })
    villi["size_factor"] = median_ratio_size_factors(counts["vt"])
    external["size_factor"] = median_ratio_size_factors(external_counts)
    pd.concat([villi, decidua, external], ignore_index=True).to_csv(out / "cohort_manifest.csv", index=False)
    return counts["vt"], counts["db"], external_counts, villi, decidua, external


def audit_primary(villi):
    primary = villi[villi.group.isin([CASE, FGR])]
    if primary.inferred_round.nunique() != 1:
        raise ValueError("Primary OpenBio arms are not in the same inferred sequencing round")
    ages = primary.groupby("group").ga_delivery_weeks.agg(["min", "max"])
    lower, upper = ages["min"].max(), ages["max"].min()
    overlap = primary[primary.ga_delivery_weeks.between(lower, upper)]
    check_groups(overlap, {CASE: 2, FGR: 5})
    if (lower, upper) != (26.714, 30.429):
        raise ValueError(f"Unexpected shared age interval: {(lower, upper)}")
    flagged = primary[(primary.group == FGR) &
                      primary.qc_flags.str.contains("no_pe_with_bp_ge_140_90", na=False, regex=False)]
    if set(flagged.donor_id) != {"S04", "S15"}:
        raise ValueError("Unexpected OpenBio isolated-FGR blood pressure flags")
    return overlap.sample_id.tolist(), sorted(flagged.donor_id)


def score_manifest(annotation, pathways, local_counts, external_counts):
    fields = ["gene", "hgnc_id", "hgnc_symbol", "ensembl_gene_id"]
    mapping = annotation[fields].dropna()
    unique = ~mapping[fields[1:]].apply(lambda col: col.duplicated(keep=False)).any(axis=1)
    mapping = mapping.loc[unique].rename(columns={"gene": "local_gene", "hgnc_symbol": "source_symbol"})
    genes = pd.DataFrame({"pathway_id": "R-HSA-201451", "source_symbol": pathways[BRANCHES[0]]})
    genes = genes.merge(mapping, on="source_symbol", how="left", validate="one_to_one")
    genes["openbio_libraries_ge10"] = genes.local_gene.map(local_counts.ge(10).sum(axis=1)).astype("Int64")
    genes["awamleh_libraries_ge10"] = genes.ensembl_gene_id.map(
        external_counts.ge(10).sum(axis=1)).astype("Int64")
    coverage = genes[["openbio_libraries_ge10", "awamleh_libraries_ge10"]]
    genes["in_score"] = coverage.ge(5).fillna(False).all(axis=1).astype(bool)
    genes["score_exclusion"] = np.select(
        [genes.in_score, genes.hgnc_id.isna(), coverage.isna().any(axis=1)],
        ["", "no unique HGNC/Ensembl mapping", "absent from one or both count matrices"],
        default="fewer than five libraries with >=10 reads in one or both cohorts",
    )
    excluded = set(genes.loc[~genes.in_score, "source_symbol"])
    if excluded != {"AMHR2", "BMP10", "CER1", "GDF2"}:
        raise ValueError(f"Unexpected fixed score gene exclusions: {excluded}")
    return genes


def donor_scores(meta, counts, fixed, gene_key):
    ids = fixed[gene_key]
    log_counts = np.log2(counts.loc[ids].div(meta.size_factor, axis=1) + 1)
    without = fixed.loc[~fixed.source_symbol.isin(["INHA", "INHBA"]), gene_key]
    return meta[DONOR_FIELDS + ["size_factor"]].assign(
        bmp27=log_counts.mean(axis=0),
        bmp25_without_inhibins=log_counts.loc[without].mean(axis=0),
    )


def score_comparisons(local, external, overlap_ids, flagged, out):
    comparisons = {
        "openbio_primary": (local, CASE, FGR),
        "openbio_overlap": (local.loc[overlap_ids], CASE, FGR),
        "awamleh_primary": (external, CASE, FGR),
        "awamleh_eope_vs_fgr": (external, "EOPE", FGR),
        "awamleh_fgr_vs_control": (external, FGR, "Preterm control"),
    }
    comparisons.update({
        "openbio_without_" + "_".join(omit): (local.loc[~local.donor_id.isin(omit)], CASE, FGR)
        for omit in ([flagged[0]], [flagged[1]], flagged)
    })
    effects, omissions, membership = [], [], []
    for name, (scores, numerator, denominator) in comparisons.items():
        selected = scores[scores.group.isin([numerator, denominator])]
        a = selected.group == numerator
        membership.append(selected[DONOR_FIELDS].assign(
            comparison=name, arm=np.where(a, "numerator", "denominator")))
        for score, column in SCORE_COLUMNS.items():
            case, control = selected.loc[a, column], selected.loc[~a, column]
            omission_rows = [
                (sample, case.drop(sample).mean() - control.mean()) for sample in case.index
            ] + [(sample, case.mean() - control.drop(sample).mean()) for sample in control.index]
            for sample, effect in omission_rows:
                omissions.append({"comparison": name, "score": score,
                                  "omitted_donor_id": selected.at[sample, "donor_id"],
                                  "omitted_group": selected.at[sample, "group"], "effect": effect})
            changes = [value for _, value in omission_rows]
            effects.append({"cohort": selected.cohort.iloc[0], "comparison": name, "score": score,
                            "numerator": numerator, "denominator": denominator,
                            "n_a": len(case), "n_b": len(control),
                            **welch_effect(case.to_numpy(), control.to_numpy()),
                            "loo_min": min(changes), "loo_max": max(changes)})
    members = pd.concat(membership, ignore_index=True)
    members[["comparison", *DONOR_FIELDS, "arm"]].to_csv(out / "comparison_membership.csv", index=False)
    pd.DataFrame(effects).to_csv(out / "bmp_effects.csv", index=False)
    pd.DataFrame(omissions).to_csv(out / "bmp_leave_one_out.csv", index=False)
    support = pd.concat([members, local.assign(comparison="openbio_all_villous"),
                         external.assign(comparison="awamleh_all_villous")], ignore_index=True)
    summary = support.groupby(["cohort", "comparison", "group"]).ga_delivery_weeks.agg(
        n="size", n_age_available="count", min_weeks="min", median_weeks="median", max_weeks="max")
    summary.to_csv(out / "cohort_support.csv")


def atlas_context(data, marker_path, shared_edge, out):
    folder = data / "extension/pseudobulk"
    counts = read_counts(folder / "track01_pseudobulk_counts.csv")
    meta = align(counts, indexed(folder / "track01_pseudobulk_metadata.csv", "pseudobulk_id"))
    donors = indexed(folder / "track01_pseudobulk_donors.csv", "donor_id")
    class_map = pd.read_csv(folder / "track01_cell_class_map.csv")
    markers = pd.read_csv(marker_path)
    if (meta.duplicated(["donor_id", "cell_class"]).any() or not meta.donor_id.isin(donors.index).all()
            or not meta.total_counts.eq(counts.sum(axis=0)).all() or (meta.total_counts <= 0).any()):
        raise ValueError("Invalid donor-class identities or pseudobulk library totals")
    classes = sorted(set(class_map.cell_class) | set(meta.cell_class) | set(markers.cell_class))
    symbols = sorted(set(shared_edge) | {"BMPR2", "BMPR1B", "GREM2", "SMAD9"})
    raw = counts.reindex(symbols).T.rename_axis(index="pseudobulk_id", columns="symbol").reset_index()
    values = raw.melt(id_vars="pseudobulk_id", value_name="raw_count").join(
        meta[["donor_id", "cell_class", "clinical_group", "pe_status", "n_cells", "total_counts"]],
        on="pseudobulk_id", validate="many_to_one")
    values["cpm"] = 1e6 * values.raw_count / values.total_counts
    values["eligible_donor_class"] = values.n_cells >= 20
    values.to_csv(out / "atlas_donors.csv", index=False)
    eligible = values[values.eligible_donor_class]
    groups = eligible.groupby(["cell_class", "symbol"])
    summaries = groups.cpm.median().to_frame("median_cpm")
    summaries["fraction_ge10"] = groups.raw_count.agg(lambda x: (x >= 10).mean())
    index = pd.MultiIndex.from_product([classes, symbols], names=["cell_class", "symbol"])
    summary = summaries.reindex(index).reset_index()
    summary["donors_released"] = summary.cell_class.map(meta.cell_class.value_counts()).fillna(0).astype(int)
    n_eligible = meta.loc[meta.n_cells >= 20, "cell_class"].value_counts()
    summary["donors_eligible"] = summary.cell_class.map(n_eligible).fillna(0).astype(int)
    summary["status"] = np.select(
        [summary.donors_released == 0, summary.donors_eligible < 3, ~summary.symbol.isin(counts.index)],
        ["class absent", "fewer than 3 eligible donors", "gene absent"], default="reported")
    summary.loc[summary.status != "reported", ["median_cpm", "fraction_ge10"]] = np.nan
    summary.to_csv(out / "atlas_summary.csv", index=False)
    return len(donors), len(meta)

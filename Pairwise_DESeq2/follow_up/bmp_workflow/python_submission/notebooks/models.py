"""Count-model contrasts, variance-stabilized PCA and Reactome enrichment."""

from itertools import chain
from pathlib import Path

import gseapy
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import false_discovery_control, norm
from sklearn.decomposition import PCA


BMP = "REACTOME_SIGNALING_BY_BMP"
ACTIVIN = "REACTOME_SIGNALING_BY_ACTIVIN"
TGF = "REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX"
PRIMARY = ("openbio_primary", "openbio_overlap", "awamleh_primary")
CONTRASTS = {
    "primary": (("EOPE+FGR", "FGR"), ("combined", "fgr")),
    "eope_vs_fgr": (("EOPE", "FGR"), ("pe", "fgr")),
    "combined_vs_control": (("EOPE+FGR", "iPTB"), ("combined", "control")),
    "eope_vs_control": (("EOPE", "iPTB"), ("pe", "control")),
    "fgr_vs_control": (("FGR", "iPTB"), ("fgr", "control")),
}
EXTERNAL_LABELS = {"combined": "EOPE+FGR", "fgr": "FGR", "pe": "EOPE", "control": "Preterm control"}
STATS = ["baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "low95", "high95"]
MEMBERSHIP = ["sample_id", "donor_id", "group", "ga_delivery_weeks"]
MARKER_CLASSES = {
    "villous_cytotrophoblast", "syncytiotrophoblast", "extravillous_trophoblast",
    "endothelial_cell", "perivascular_cell", "decidual_stromal_cell", "decidual_macrophage",
}


def _local_model(counts, cohort, cpus):
    metadata = cohort.set_index("sample_id")[["fetal_sex", "group"]].copy()
    metadata["group"] = metadata["group"].replace({"EOPE+FGR": "EOPE_FGR"})
    dds = DeseqDataSet(
        counts=counts.loc[:, metadata.index].T, metadata=metadata,
        design="~ fetal_sex + group", fit_type="parametric", size_factors_fit_type="ratio",
        refit_cooks=True, n_cpus=cpus, quiet=True,
    )
    design = dds.obsm["design_matrix"].to_numpy()
    if np.linalg.matrix_rank(design) != design.shape[1] or design.shape[0] <= design.shape[1]:
        raise ValueError("OpenBio count model requires a full-rank design with residual degrees of freedom.")
    dds.deseq2()
    return dds


def _local_test(dds, counts, cohort, arms, cpus):
    labels = [arm.replace("+", "_") for arm in arms]
    stats = DeseqStats(
        dds, contrast=["group", *labels], alpha=0.05, cooks_filter=True,
        independent_filter=True, n_cpus=cpus, quiet=True,
    )
    selected = cohort.loc[cohort["group"].isin(arms), "sample_id"]
    no_arm_reads = counts.loc[dds.var_names, selected].sum(axis=1).eq(0)
    # A gene expressed only in other groups supplies a null contrast. Apply it
    # before summary() performs Cook filtering and multiple-testing adjustment.
    stats.run_wald_test()
    null_contrast = no_arm_reads & stats.p_values.notna()
    stats.p_values.loc[null_contrast] = 1.0
    stats.statistics.loc[null_contrast] = 0.0
    stats.summary()
    result = stats.results_df.loc[:, STATS[:6]].copy()
    result["both_groups_zero"] = no_arm_reads
    result.loc[null_contrast, "log2FoldChange"] = 0.0
    half_width = norm.ppf(0.975) * result["lfcSE"].mask(no_arm_reads)
    result["low95"] = result["log2FoldChange"] - half_width
    result["high95"] = result["log2FoldChange"] + half_width
    return result


def _save_pca(dds, villi, out):
    dds.vst(use_design=False)
    values = dds.layers["vst_counts"]
    variance = values.var(axis=0)
    top = np.argsort(-variance, kind="stable")[:500]
    pca = PCA(svd_solver="full")
    coordinates = pca.fit_transform(values[:, top])
    result = villi[["sample_id", "donor_id", "group", "fetal_sex", "inferred_round"]].copy()
    result["PC1"], result["PC2"] = coordinates[:, 0], coordinates[:, 1]
    result["pc1_percent"], result["pc2_percent"] = 100 * pca.explained_variance_ratio_[:2]
    result.to_csv(out / "pca_scores.csv", index=False)


def _effect_rows(results, fit, cohort, arms, mapping, key):
    if results.empty or not results.index.is_unique:
        raise ValueError(f"Empty or duplicate count-model results in {fit}.")
    frame = results[STATS + ["both_groups_zero"]].copy()
    symbols = mapping.set_index(key)
    original = symbols.reindex(frame.index)
    frame["symbol"] = original["hgnc_symbol"]
    frame["ensembl_gene_id"] = original["ensembl_gene_id"] if key == "gene" else frame.index
    frame["cohort"] = cohort
    frame["fit"] = fit
    frame["numerator"], frame["denominator"] = arms
    return frame.rename_axis("input_id").reset_index()


def _eligible(rows, annotation, key):
    # Ambiguity is resolved among genes that survive this fitted contrast.
    candidates = rows.loc[
        ~rows["both_groups_zero"]
        & np.isfinite(rows["stat"]) & np.isfinite(rows["pvalue"]) & rows["baseMean"].ge(10),
        ["input_id", "stat"],
    ]
    pairs = annotation[[key, "hgnc_symbol"]].dropna().drop_duplicates().rename(
        columns={key: "input_id", "hgnc_symbol": "symbol"}
    )
    joined = candidates.merge(pairs, on="input_id", how="inner", validate="one_to_many")
    keep = ~joined["input_id"].duplicated(keep=False) & ~joined["symbol"].duplicated(keep=False)
    return joined.loc[keep]


def _membership(fit, name, cohort, arms):
    frame = cohort[MEMBERSHIP].copy()
    frame["role"] = frame["group"].map({arms[0]: "numerator", arms[1]: "denominator"}).fillna("model_only")
    frame.insert(0, "cohort", name)
    frame.insert(0, "fit", fit)
    return frame


def _save_enrichment(ranks, pathways, out):
    branches = []
    for fit, ranks_in_fit in ranks.items():
        result = gseapy.prerank(
            rnk=ranks_in_fit, gene_sets=pathways, min_size=10, max_size=500, threads=1,
            seed=1049, ascending=None, no_plot=True, method="multilevel", sample_size=101,
            permutation_num=1000, eps=1e-10,
        ).res2d.set_index("Term")
        if set(result.index) != set(pathways):
            raise ValueError(f"Expected all three Reactome branches in {fit}.")
        result = result.loc[list(pathways)]
        pvalues = result["NOM p-val"].to_numpy(dtype=float)
        if not np.isfinite(pvalues).all() or not result["NES"].astype(float).map(np.isfinite).all():
            raise ValueError(f"Nonfinite Reactome enrichment in {fit}.")
        branches.append(pd.DataFrame({
            "fit": fit, "rank_size": len(ranks_in_fit), "pathway": result.index.to_list(),
            "pvalue": pvalues, "fdr": false_discovery_control(pvalues, method="bh"),
            "log2err": result["log2err"].to_numpy(), "ES": result["ES"].to_numpy(),
            "NES": result["NES"].to_numpy(),
            "size": [len(set(pathways[key]).intersection(ranks_in_fit.index)) for key in result.index],
            "leading_edge": result["Lead_genes"].to_numpy(),
        }))
    branch = pd.concat(branches, ignore_index=True)
    branch.to_csv(out / "branch_enrichment.csv", index=False)

    primary = branch.loc[branch["pathway"].eq(BMP) & branch["fit"].isin(PRIMARY)]
    edge_sets = [set(edge.split(";")) - {""} if nes > 0 else set()
                 for edge, nes in zip(primary["leading_edge"], primary["NES"], strict=True)]
    shared_edge = set.intersection(*edge_sets)
    pd.DataFrame({"symbol": sorted(shared_edge)}).to_csv(out / "shared_leading_edge.csv", index=False)
    plot = []
    for fit in PRIMARY:
        score = ranks[fit]
        hit = score.index.isin(pathways[BMP])
        hits = score.abs().to_numpy() * hit
        plot.append(pd.DataFrame({
            "fit": fit, "rank": np.arange(1, len(score) + 1), "symbol": score.index,
            "wald": score.to_numpy(), "in_bmp": hit, "shared_edge": score.index.isin(shared_edge),
            "running_score": np.cumsum(hits / hits.sum() - (~hit) / (~hit).sum()),
        }))
    pd.concat(plot, ignore_index=True).to_csv(out / "rank_plot.csv", index=False)


def _save_markers(tests, marker_path, out):
    markers = pd.read_csv(marker_path).drop_duplicates()
    markers = markers.loc[markers["cell_class"].isin(MARKER_CLASSES) & markers["gene_symbol"].ne("MKI67")]
    keys = ["symbol", "input_id", "log2FoldChange", "low95", "high95", "pvalue", "padj"]
    frames = []
    for fit in PRIMARY:
        genes = tests[fit][keys].dropna(subset=["symbol"])
        joined = markers.merge(genes, left_on="gene_symbol", right_on="symbol", how="left",
                               validate="many_to_one").drop(columns="symbol")
        joined.insert(0, "fit", fit)
        frames.append(joined)
    pd.concat(frames, ignore_index=True).to_csv(out / "bulk_marker_effects.csv", index=False)


def fit_and_rank(vt_counts, villi, overlap_ids, annotation, pathways, awamleh, external,
                 marker_path, out, cpus):
    """Export model, PCA, enrichment and marker tables and return the rank universe."""
    out, awamleh = Path(out), Path(awamleh)
    if set(pathways) != {BMP, ACTIVIN, TGF}:
        raise ValueError("Expected BMP, activin and TGF-beta receptor Reactome branches.")
    overlap_ids = list(overlap_ids)
    overlap = villi.set_index("sample_id").loc[overlap_ids].reset_index()
    mapping = annotation.dropna(subset=["hgnc_id", "hgnc_symbol", "ensembl_gene_id"])
    unambiguous = ~np.logical_or.reduce([
        mapping[column].duplicated(keep=False) for column in ("hgnc_id", "hgnc_symbol", "ensembl_gene_id")
    ])
    mapping = mapping.loc[unambiguous]
    tests, members = {}, []
    print("Fitting Python OpenBio count models", flush=True)
    full_model = _local_model(vt_counts, villi, cpus)
    overlap_model = _local_model(vt_counts, overlap, cpus)
    local = chain(((f"openbio_{suffix}", full_model, villi, local_arms)
                   for suffix, (local_arms, _) in CONTRASTS.items()),
                  [("openbio_overlap", overlap_model, overlap, CONTRASTS["primary"][0])])
    for fit, model, cohort, arms in local:
        result = _local_test(model, vt_counts, cohort, arms, cpus)
        tests[fit] = _effect_rows(result, fit, "OpenBio", arms, mapping, "gene")
        members.append(_membership(fit, "OpenBio", cohort, arms))
    for suffix, (_, coded_arms) in CONTRASTS.items():
        fit = f"awamleh_{suffix}"
        numerator, denominator = coded_arms
        file = awamleh / f"gse114691_{numerator}_vs_{denominator}.tsv.gz"
        result = pd.read_csv(file, sep="\t", index_col="ensembl_gene_id")
        # The original Awamleh fit uses a global feature filter and no separate both-arm filter.
        result["both_groups_zero"] = False
        arms = tuple(EXTERNAL_LABELS[label] for label in coded_arms)
        tests[fit] = _effect_rows(result, fit, "Awamleh", arms, mapping, "ensembl_gene_id")
        members.append(_membership(fit, "Awamleh", external, arms))
    pd.concat(members, ignore_index=True).to_csv(out / "rank_fit_membership.csv", index=False)
    _save_pca(full_model, villi, out)

    eligible = {fit: _eligible(table, annotation, "gene" if fit.startswith("openbio_") else "ensembl_gene_id")
                for fit, table in tests.items()}
    common = sorted(set.intersection(*(set(table["symbol"]) for table in eligible.values())))
    rank_ids = eligible["awamleh_primary"].set_index("symbol")["input_id"].reindex(common)
    universe = pd.DataFrame({"symbol": common, "ensembl_gene_id": rank_ids.to_numpy()})
    for key, branch in {"in_bmp": BMP, "in_activin": ACTIVIN, "in_tgf_beta_receptor": TGF}.items():
        universe[key] = universe["symbol"].isin(pathways[branch])
    universe.to_csv(out / "rank_universe.csv", index=False)
    ranks, exports = {}, []
    for fit, table in tests.items():
        ranked = eligible[fit].set_index("symbol")["stat"].reindex(common)
        ranks[fit] = ranked.sort_values(ascending=False, kind="stable")
        exported = table.copy()
        eligible_symbol = exported["input_id"].map(eligible[fit].set_index("input_id")["symbol"])
        exported["symbol"] = eligible_symbol.fillna(exported["symbol"])
        matched = eligible_symbol.isin(common)
        exported["in_rank_universe"] = matched
        exported.loc[matched, "ensembl_gene_id"] = exported.loc[matched, "symbol"].map(rank_ids)
        exported["in_bmp"] = exported["symbol"].isin(pathways[BMP])
        exports.append(exported)
    columns = ["fit", "cohort", "numerator", "denominator", "input_id", "symbol", "ensembl_gene_id",
               *STATS, "both_groups_zero", "in_rank_universe", "in_bmp"]
    pd.concat(exports, ignore_index=True)[columns].to_csv(
        out / "bmp_gene_effects.tsv.gz", sep="\t", index=False, compression={"method": "gzip", "mtime": 0}
    )
    print(f"Testing Reactome branches across {len(common):,} shared genes", flush=True)
    _save_enrichment(ranks, pathways, out)
    _save_markers(tests, marker_path, out)
    return set(common)

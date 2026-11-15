"""Run the targeted OpenBio bulk and paired-tissue RNA pilot from the repo root."""

import argparse
import os
from hashlib import sha256
from importlib.metadata import version
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import matplotlib
import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy import stats

matplotlib.use("Agg")
from matplotlib import pyplot as plt


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "openbio-datathon-2026/track-01/data"
ROUND_AUDIT = ROOT / "Deci_Vili_H1/results/step01_cohort_audit/donor_audit.csv"
PANEL = ("IFITM2", "IFITM3", "GBP5", "ERVW-1", "ERVFRD-1", "FURIN", "PCSK7", "GCM1")
SCORES = {
    "stb_markers": ("CGA", "PSG1", "SDC1", "CSH1"),
    "ifn": ("ISG15", "MX1", "OAS1", "OAS2", "IFIT1", "IFIT2", "IFI44L", "RSAD2",
            "CXCL9", "CXCL10", "IDO1", "IRF1", "PSMB9", "TAP1", "STAT1"),
    "hypoxia": ("LDHA", "PDK1", "SLC2A1", "EGLN3", "ALDOA", "BNIP3"),
    "myeloid": ("CD68", "LST1", "CSF1R"),
}


def verify_inputs():
    paths = ("track01_counts_vt.csv", "track01_metadata_vt.csv", "qc/track01_htseq_summary.csv",
             "extension/decidua/track01_counts_db.csv", "extension/decidua/track01_metadata_db.csv")
    manifest = dict(line.split()[::-1] for line in (DATA / "provenance/checksums.sha256").read_text().splitlines())
    for path in paths:
        if sha256((DATA / path).read_bytes()).hexdigest() != manifest[path]:
            raise ValueError(f"Count release checksum differs for {path}")


def load_tissue(tissue):
    folder = DATA if tissue == "vt" else DATA / "extension/decidua"
    raw = pd.read_csv(folder / f"track01_counts_{tissue}.csv", index_col="gene").T
    meta = pd.read_csv(folder / f"track01_metadata_{tissue}.csv", index_col="sample_id")
    assert raw.index.is_unique and raw.columns.is_unique and meta.index.is_unique
    assert meta["donor_id"].is_unique and set(raw.index) == set(meta.index)
    assert raw.ge(0).all().all() and all(pd.api.types.is_integer_dtype(d) for d in raw.dtypes)
    qc = pd.read_csv(DATA / "qc/track01_htseq_summary.csv", index_col="sample_id")
    assert raw.sum(axis=1).eq(qc.loc[raw.index, "assigned_to_genes"]).all()
    meta = meta.loc[raw.index].copy()
    meta["combined"] = meta["group"].eq("EOPE+FGR").astype(int)
    meta["age"] = meta["ga_delivery_weeks"]
    return raw, meta


def normalize(raw):
    positive = raw.loc[:, raw.gt(0).all()]
    log_positive = np.log(positive)
    size_factors = np.exp(log_positive.sub(log_positive.mean(axis=0), axis=1).median(axis=1))
    return np.log2(raw.div(size_factors, axis=0) + 1)


def overlap(meta):
    ages = meta.groupby("combined")["age"].agg(["min", "max"])
    return meta.loc[meta["age"].between(ages["min"].max(), ages["max"].min())]


def zscore(values):
    return (values - values.mean()) / values.std(ddof=0)


def ols(data, response, terms, effect, label):
    y = data[response].to_numpy(dtype=float)
    x = np.column_stack((np.ones(len(data)), data[terms].to_numpy(dtype=float)))
    assert len(data) > x.shape[1] and np.linalg.matrix_rank(x) == x.shape[1], label
    beta = np.linalg.lstsq(x, y, rcond=None)[0]
    residual = y - x @ beta
    inv = np.linalg.inv(x.T @ x)
    leverage = np.einsum("ij,jk,ik->i", x, inv, x)
    scaled_x = x * (residual / (1 - leverage))[:, None]
    covariance = inv @ scaled_x.T @ scaled_x @ inv
    j = terms.index(effect) + 1
    se = np.sqrt(covariance[j, j])
    df = len(y) - x.shape[1]
    width = stats.t.ppf(0.975, df) * se
    return {"analysis": label, "endpoint": response, "term": effect,
            "estimate": beta[j], "low95": beta[j] - width, "high95": beta[j] + width,
            "n_combined": int(data["combined"].sum()), "n_eope_only": int((1 - data["combined"]).sum())}


def differential_expression(counts, meta):
    selected = counts.columns[counts.ge(10).sum(axis=0).ge(3)]
    model_meta = meta.assign(age_centered=meta["age"] - meta["age"].mean())
    model = DeseqDataSet(counts=counts.loc[:, selected], metadata=model_meta,
                         design="~ age_centered + group", refit_cooks=False, n_cpus=1, quiet=True)
    matrix = model.obsm["design_matrix"]
    assert np.linalg.matrix_rank(matrix) == matrix.shape[1]
    model.deseq2()
    result = DeseqStats(model, contrast=["group", "EOPE+FGR", "EOPE"],
                        independent_filter=False, n_cpus=1, quiet=True)
    result.summary()
    panel = result.results_df.loc[list(PANEL)].copy()
    panel["low95"] = panel["log2FoldChange"] - 1.96 * panel["lfcSE"]
    panel["high95"] = panel["log2FoldChange"] + 1.96 * panel["lfcSE"]
    panel.loc[list(PANEL[1:]), "secondary_panel_q"] = stats.false_discovery_control(
        panel.loc[list(PANEL[1:]), "pvalue"].fillna(1))
    panel["cooks_outlier"] = model.var.loc[panel.index, "_pvalue_cooks_outlier"]
    panel["lfc_converged"] = model.var.loc[panel.index, "_LFC_converged"]
    return panel


def write_figure(panel, samples, path):
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.8), layout="constrained")
    ax, scatter = axes
    items = [(gene, row.log2FoldChange, row.low95, row.high95) for gene, row in panel.iterrows()]
    y = np.arange(len(items))
    points = np.array([item[1:] for item in items])
    ax.errorbar(points[:, 0], y, xerr=[points[:, 0] - points[:, 1], points[:, 2] - points[:, 0]],
                fmt="o", color="#235466", capsize=3)
    ax.axvline(0, color="#737780", linewidth=1)
    ax.set_yticks(y, [item[0] for item in items])
    ax.invert_yaxis()
    ax.set(xlabel="Log₂ fold change (95% CI)", title="Villous RNA, adjusted for delivery age")
    for group, color, marker in [("EOPE", "#235466", "o"), ("EOPE+FGR", "#A13D28", "D")]:
        subset = samples[samples["group"].eq(group)]
        scatter.scatter(subset["IFITM2"], subset["stb_markers"], label=group, color=color, marker=marker)
    scatter.axhline(0, color="#DDDFE2", linewidth=1, zorder=0)
    scatter.set(xlabel="IFITM2, log₂(normalized count + 1)", ylabel="Syncytiotrophoblast marker score",
                title="One point per pregnancy")
    scatter.legend(frameon=False)
    for axis in axes:
        axis.spines[["right", "top"]].set_visible(False)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main(output):
    verify_inputs()
    raw, meta = load_tissue("vt")
    db_raw, db_meta = load_tissue("db")
    audit = pd.read_csv(ROUND_AUDIT).set_index("donor_id")
    assert set(meta["donor_id"]) == set(audit.index)
    assert meta["group"].eq(meta["donor_id"].map(audit["group"])).all()
    assert meta["age"].eq(meta["donor_id"].map(audit["ga_delivery_weeks"])).all()
    meta["round"] = meta["donor_id"].map(audit["sequencing_round"])
    logged = normalize(raw)
    logged_db = normalize(db_raw)
    cases = meta[meta["group"].isin(("EOPE", "EOPE+FGR"))].copy()
    assert cases["group"].value_counts().to_dict() == {"EOPE": 8, "EOPE+FGR": 7}
    panel = differential_expression(raw.loc[cases.index], cases)

    cases = cases.join(logged.loc[:, list(PANEL)])
    for name, genes in SCORES.items():
        cases[name] = logged.loc[cases.index, list(genes)].apply(zscore).mean(axis=1)
    cases["ifitm2_z"] = zscore(cases["IFITM2"])
    cases["stb_z"] = zscore(cases["stb_markers"])
    cases["gcm1_z"] = zscore(cases["GCM1"])
    same_round = cases[cases["round"].eq("other_round")]
    checks = [ols(cases, "IFITM2", ["combined", "age"], "combined", "full_age_adjusted_ols")]
    for label, subset in (("full_unadjusted", cases), ("delivery_age_overlap", overlap(cases)),
                          ("inferred_same_round", same_round), ("same_round_age_overlap", overlap(same_round))):
        checks.append(ols(subset, "IFITM2", ["combined"], "combined", label))
    for donor in cases.index:
        checks.append(ols(cases.drop(donor), "IFITM2", ["combined", "age"], "combined", f"omit_{donor}"))
    for score in SCORES:
        checks.append(ols(cases, score, ["combined", "age"], "combined", "score_group_difference"))
        checks.append(ols(cases, "ifitm2_z", [score, "combined", "age"], score, f"IFITM2_vs_{score}"))

    for response in ("stb_z", "gcm1_z"):
        terms = ["ifitm2_z", "combined", "age"]
        for covariate in (None, "ifn", "hypoxia", "myeloid"):
            extra = [] if covariate is None else [covariate]
            checks.append(ols(cases, response, terms + extra, "ifitm2_z", f"conditional_{covariate or 'base'}"))

    db = db_meta.join(logged_db.loc[:, ["IFITM2", "IFITM3", "GBP5"]]).set_index("donor_id")
    db_fields = db[["group", "age", "IFITM2", "IFITM3", "GBP5"]].add_suffix("_db")
    pairs = cases.set_index("donor_id").join(db_fields, how="inner")
    assert pairs["group"].eq(pairs["group_db"]).all() and pairs["age"].eq(pairs["age_db"]).all()
    assert pairs["group"].value_counts().to_dict() == {"EOPE": 8, "EOPE+FGR": 4}
    for gene in ("IFITM2", "IFITM3", "GBP5"):
        pairs[f"{gene}_vt_minus_db"] = pairs[gene] - pairs[f"{gene}_db"]
        checks.append(ols(pairs, f"{gene}_vt_minus_db", ["combined", "age"], "combined", "paired_tissue"))
        for endpoint in (gene, f"{gene}_db"):
            checks.append(ols(pairs, endpoint, ["combined", "age"], "combined", "complete_pair_compartment"))

    output.mkdir(parents=True, exist_ok=True)
    panel.to_csv(output / "villous_deseq2_panel.csv")
    pd.DataFrame(checks).to_csv(output / "delivery_rna_sensitivity.csv", index=False)
    donor_values = meta[["donor_id", "group", "age", "round"]].join(logged.loc[:, list(PANEL)])
    donor_values.join(cases[list(SCORES)]).to_csv(output / "villous_donor_values.csv")
    write_figure(panel, cases, output / "delivery_rna_pilot.png")
    print(f"PyDESeq2 {version('pydeseq2')}, SciPy {version('scipy')}, pandas {version('pandas')}")
    print(f"Wrote four outputs to {output}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "results")
    main(parser.parse_args().output)

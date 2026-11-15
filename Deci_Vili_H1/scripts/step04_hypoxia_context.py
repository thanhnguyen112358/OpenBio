"""Score villous Hallmark Hypoxia and summarize placental cell-class context."""

import hashlib
import json
from importlib.metadata import version
from urllib.request import urlopen

import matplotlib
import numpy as np
import pandas as pd

from step01_audit_cohort import PROJECT, SOURCE, VT_COUNTS, verify_source
from step02_deseq2_scores import AUDIT, COHORTS, OUTPUT as STEP2

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


OUTPUT = PROJECT / "results" / "step04_hypoxia_context"
HALLMARK_URL = (
    "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/"
    "2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt"
)
HALLMARK = "HALLMARK_HYPOXIA"
ANNOTATION = "track01_gene_annotation.csv"
PB_METADATA = "extension/pseudobulk/track01_pseudobulk_metadata.csv"
PB_COUNTS = "extension/pseudobulk/track01_pseudobulk_counts.csv"
BROWSER_GENES = ("TREM2", "C1QA", "C1QC")
MIN_CELLS = 10
GROUP_STYLES = {"EOPE": ("#994438", "o"), "FGR": ("#235D70", "s")}


def load_hallmark():
    # Use the organizer's pinned release from the public Broad archive.
    with urlopen(HALLMARK_URL, timeout=20) as response:
        payload = response.read()
    matches = [
        line.split("\t")[2:]
        for line in payload.decode("utf-8-sig").splitlines()
        if line.startswith(f"{HALLMARK}\t")
    ]
    if len(matches) != 1 or not matches[0] or len(set(matches[0])) != len(matches[0]):
        raise ValueError(f"Expected one nonempty {HALLMARK} record with unique identifiers")
    return matches[0], hashlib.sha256(payload).hexdigest()


def score_hypoxia(donors, factors, genes):
    samples = donors.vt_sample_id.tolist()
    raw = pd.read_csv(SOURCE / VT_COUNTS, index_col="gene", usecols=["gene", *samples])
    annotation = pd.read_csv(SOURCE / ANNOTATION)
    mapping, aliases, unmatched = {}, {}, []
    for symbol in genes:
        if symbol in raw.index:
            mapping[symbol] = symbol
            continue
        candidates = annotation.loc[annotation.hgnc_symbol.eq(symbol), "gene"]
        if len(candidates) == 1 and candidates.iloc[0] in raw.index:
            mapping[symbol] = candidates.iloc[0]
            aliases[symbol] = candidates.iloc[0]
        else:
            unmatched.append(symbol)
    if not mapping or len(set(mapping.values())) != len(mapping):
        raise ValueError("Hypoxia gene mapping is empty or reuses a bulk count row")

    selected = raw.loc[list(mapping.values()), samples].T
    selected.index = donors.index
    logged = np.log2(selected.div(factors, axis=0) + 1)
    hypoxia = logged.sub(logged.mean(axis=0), axis=1).mean(axis=1).rename("vt_hallmark_hypoxia")
    return hypoxia, aliases, unmatched


def summarize_overlap(scores):
    rows = []
    for cohort in COHORTS:
        selected = scores.loc[scores[f"in_{cohort}"].eq(1)]
        eope = selected.loc[selected.group.eq("EOPE"), "vt_hallmark_hypoxia"]
        fgr = selected.loc[selected.group.eq("FGR"), "vt_hallmark_hypoxia"]
        low, high = max(eope.min(), fgr.min()), min(eope.max(), fgr.max())
        overlap = low <= high
        eope_shared = eope.loc[eope.between(low, high)].index
        fgr_shared = fgr.loc[fgr.between(low, high)].index
        rows.append({
            "cohort": cohort,
            "n_eope": len(eope),
            "n_fgr": len(fgr),
            "eope_min": eope.min(),
            "eope_max": eope.max(),
            "fgr_min": fgr.min(),
            "fgr_max": fgr.max(),
            "shared_min": low if overlap else np.nan,
            "shared_max": high if overlap else np.nan,
            "eope_n_in_shared_range": len(eope_shared),
            "fgr_n_in_shared_range": len(fgr_shared),
            "eope_in_shared_range": " ".join(eope_shared),
            "fgr_in_shared_range": " ".join(fgr_shared),
        })
    return pd.DataFrame(rows)


def summarize_controls():
    metadata = pd.read_csv(SOURCE / PB_METADATA, index_col="pseudobulk_id")
    controls = metadata.loc[metadata.pe_status.eq("control") & metadata.n_cells.ge(MIN_CELLS)]
    if controls.duplicated(["donor_id", "cell_class"]).any():
        raise ValueError("A reference donor appears more than once in a cell class")
    raw = pd.read_csv(SOURCE / PB_COUNTS, index_col="gene", usecols=["gene", *controls.index])
    totals = raw.sum(axis=0).reindex(controls.index)
    np.testing.assert_array_equal(totals, controls.total_counts)
    counts = raw.loc[list(BROWSER_GENES), controls.index].T
    cpm = counts.div(totals, axis=0).mul(1_000_000)
    values = pd.concat({"cpm": cpm.stack(), "detected": counts.gt(0).stack()}, axis=1)
    values = values.rename_axis(index=["pseudobulk_id", "gene"]).reset_index()
    values = values.join(controls[["donor_id", "cell_class", "onset_group"]], on="pseudobulk_id")
    values["early_control"] = values.onset_group.eq("early")
    values["late_control"] = values.onset_group.eq("late")
    summary = values.groupby(["gene", "cell_class"]).agg(
        donor_count=("donor_id", "nunique"),
        early_control_donors=("early_control", "sum"),
        late_control_donors=("late_control", "sum"),
        detected_donors=("detected", "sum"),
        median_cpm=("cpm", "median"),
    )
    return summary.reset_index().sort_values(["gene", "median_cpm"], ascending=[True, False])


def plot_hypoxia(scores, overlap):
    fig, ax = plt.subplots(figsize=(8.5, 5.5), layout="constrained")
    low, high = overlap.loc[
        overlap.cohort.eq("inferred_same_round"), ["shared_min", "shared_max"]
    ].iloc[0]
    if pd.notna(low):
        ax.axvspan(low, high, color="#E2E9ED", zorder=0)
    for donor, row in scores.iterrows():
        color, marker = GROUP_STYLES[row.group]
        primary = row.in_inferred_same_round
        x, y = row.vt_hallmark_hypoxia, row.db_c1q_trem2
        ax.scatter(x, y, marker=marker, facecolor=color if primary else "white",
                   edgecolor=color, linewidth=1.3, s=58, zorder=3)
        offset = -7 if x > 0.39 else 7
        ax.annotate(donor, (x, y), xytext=(offset, 0), textcoords="offset points",
                    ha="right" if offset < 0 else "left", va="center", fontsize=9)
    ax.set_title("Villous hypoxia RNA score and decidual C1Q–TREM2", loc="left", weight="bold")
    ax.set_xlabel("Villous Hallmark Hypoxia score (centered log₂ units)")
    ax.set_ylabel("Decidual C1Q–TREM2 score (log₂ units)")
    ax.margins(x=0.13, y=0.12)
    ax.spines[["right", "top"]].set_visible(False)
    handles = [
        Line2D([], [], marker=marker, color=color, linestyle="none", label=label)
        for label, (color, marker) in zip(("EOPE without FGR", "Isolated FGR"), GROUP_STYLES.values())
    ]
    handles.extend([
        Line2D([], [], marker="o", color="#555555", markerfacecolor="white", linestyle="none",
               label="Outside the primary cohort"),
        Patch(facecolor="#E2E9ED", label="Shared score range in the primary cohort"),
    ])
    ax.legend(handles=handles, frameon=False, loc="upper center",
              bbox_to_anchor=(0.5, -0.15), ncol=2, fontsize=9)
    fig.savefig(OUTPUT / "figures" / "hypoxia_vs_decidual_score.png", dpi=200, facecolor="white")
    plt.close(fig)


def main():
    verify_source(ANNOTATION, PB_METADATA, PB_COUNTS)
    audit = pd.read_csv(AUDIT, index_col="donor_id")
    donors = audit.loc[audit.in_all_eope_fgr_pairs.eq(1)]
    scores = pd.read_csv(STEP2 / "donor_scores.csv", index_col="donor_id")
    membership = [f"in_{cohort}" for cohort in COHORTS]
    metadata = ["group", "ga_delivery_weeks", "sequencing_round", *membership]
    pd.testing.assert_frame_equal(scores[metadata], donors[metadata])
    genes, gmt_sha256 = load_hallmark()
    hypoxia, aliases, unmatched = score_hypoxia(donors, scores.vt_size_factor, genes)
    plot_data = scores[metadata + ["db_c1q_trem2"]].join(hypoxia)
    overlap = summarize_overlap(plot_data)
    cells = summarize_controls()

    (OUTPUT / "figures").mkdir(parents=True, exist_ok=True)
    plot_data.to_csv(OUTPUT / "hypoxia_by_donor.csv")
    overlap.to_csv(OUTPUT / "hypoxia_overlap.csv", index=False)
    cells.to_csv(OUTPUT / "cell_class_expression.csv", index=False)
    run_info = {
        "gene_set": HALLMARK,
        "hallmark_archive_url": HALLMARK_URL,
        "hallmark_archive_sha256": gmt_sha256,
        "n_hallmark_genes": len(genes),
        "n_bulk_genes_used": len(genes) - len(unmatched),
        "alias_mapping_to_bulk_row": aliases,
        "unmatched_gene_symbols": unmatched,
        "normalization_donor_ids": donors.index.tolist(),
        "minimum_cells_per_control_pseudobulk": MIN_CELLS,
        "packages": {name: version(name) for name in ("pandas", "numpy", "matplotlib")},
    }
    (OUTPUT / "run_info.json").write_text(json.dumps(run_info, indent=2) + "\n", encoding="utf-8")
    plot_hypoxia(plot_data, overlap)
    print(overlap.to_string(index=False))
    print(cells.groupby("gene").head(3).to_string(index=False))


if __name__ == "__main__":
    main()

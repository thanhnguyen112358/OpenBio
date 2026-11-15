"""Shared inputs and cohort definitions for the Q1 VST PCA and DESeq2 scripts."""

from pathlib import Path

import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
DATA_DIR = HERE.parent.parent / "openbio-datathon-2026" / "track-01" / "data"
RESULTS_DIR = HERE / "results"

CASE, CONTROL = "EOPE", "iPTB"
GROUP_ORDER = [CASE, CONTROL, "EOPE+FGR", "FGR", "TB"]      # the Q1 comparison first
GROUP_COLORS = {CONTROL: "tab:blue", CASE: "tab:orange"}   # other groups are drawn in grey
FLAGGED_SAMPLE = "S22_VT"
HTSEQ_COUNTERS = ["assigned_to_genes", "no_feature", "ambiguous", "too_low_aQual",
                  "not_aligned", "alignment_not_unique"]


def load_table(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty organizer data file: {path}")
    return pd.read_csv(path, index_col=0)


def load_inputs():
    """Return villous counts (genes x samples), metadata, and the HTSeq QC table."""
    counts = load_table(DATA_DIR / "track01_counts_vt.csv")
    metadata = load_table(DATA_DIR / "track01_metadata_vt.csv")
    qc = load_table(DATA_DIR / "qc" / "track01_htseq_summary.csv")
    if not counts.columns.is_unique or not metadata.index.is_unique:
        raise ValueError("Duplicate sample IDs in the counts or metadata")
    if set(counts.columns) != set(metadata.index):
        raise ValueError("Villous counts and metadata have different sample IDs")
    if not counts.index.is_unique:
        raise ValueError("Duplicate genes in the villous count table")
    if not metadata["donor_id"].is_unique:
        raise ValueError("A donor appears more than once in the villous metadata")
    raw = counts.to_numpy()
    if not (np.isfinite(raw) & (raw >= 0) & (raw % 1 == 0)).all():
        raise ValueError("Villous counts must be finite, nonnegative whole numbers")
    missing_qc = set(counts.columns) - set(qc.index)
    if missing_qc:
        raise ValueError(f"Samples missing from the HTSeq QC table: {sorted(missing_qc)}")
    return counts, metadata.loc[counts.columns], qc.loc[counts.columns, HTSEQ_COUNTERS]


def q1_cohorts(metadata):
    """The Q1 cohorts: all EOPE-without-FGR vs. iPTB donors, and the shared-age subset.

    The shared window runs from the youngest to the oldest iPTB delivery age,
    as in PCA/organizer_pca.py.
    """
    primary = metadata[metadata["group"].isin([CASE, CONTROL])]
    control_ages = primary.loc[primary["group"].eq(CONTROL), "ga_delivery_weeks"]
    low, high = control_ages.min(), control_ages.max()
    shared = primary[primary["ga_delivery_weeks"].between(low, high)]
    return {"eope_vs_iptb": primary, "shared_age": shared}, (low, high)

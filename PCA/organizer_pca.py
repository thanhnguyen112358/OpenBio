"""Run the organizer's Track 1 log-CPM PCA on the local villous data."""

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA


PCA_DIR = Path(__file__).resolve().parent
DATA_DIR = PCA_DIR.parent / "openbio-datathon-2026" / "track-01" / "data"
RESULTS_DIR = PCA_DIR / "results"
OVERVIEW_FIELDS = (
    "group", "ga_delivery_weeks", "control_type", "fgr_status", "fetal_sex",
    "maternal_age_years",
)


def load_table(path):
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty organizer data file: {path}")
    return pd.read_csv(path, index_col=0)


def load_inputs():
    counts = load_table(DATA_DIR / "track01_counts_vt.csv")
    metadata = load_table(DATA_DIR / "track01_metadata_vt.csv")
    if not counts.columns.is_unique or not metadata.index.is_unique:
        raise ValueError("Duplicate sample IDs in the counts or metadata")
    if set(counts.columns) != set(metadata.index):
        raise ValueError("Villous counts and metadata have different sample IDs")
    if not counts.index.is_unique:
        raise ValueError("Duplicate genes in the villous count table")
    raw_counts = counts.to_numpy()
    if not (np.isfinite(raw_counts) & (raw_counts >= 0) & (raw_counts % 1 == 0)).all():
        raise ValueError("Villous counts must be finite, nonnegative whole numbers")
    if (counts.sum(axis=0) == 0).any():
        raise ValueError("Villous samples must have positive library sizes for CPM")
    return counts, metadata.loc[counts.columns]


def fit_pca(counts):
    # The organizer's Section 3 transformation and PCA. Seed the solver for reruns.
    log_cpm = np.log2(counts.div(counts.sum(axis=0), axis=1) * 1e6 + 1)
    pca = PCA(n_components=2, random_state=0)
    scores = pca.fit_transform(log_cpm.T)
    return scores, pca.explained_variance_ratio_


def plot_pca(scores, variance, metadata, color_by, title, output_file):
    fig, ax = plt.subplots(figsize=(7, 5))
    values = metadata[color_by]
    if pd.api.types.is_numeric_dtype(values):
        observed = values.notna().to_numpy()
        points = ax.scatter(scores[observed, 0], scores[observed, 1],
                            c=values[observed], cmap="viridis")
        fig.colorbar(points, ax=ax, label=color_by)
        if not observed.all():
            ax.scatter(scores[~observed, 0], scores[~observed, 1], c="grey", label="NA")
            ax.legend()
    else:
        categories = values.fillna("NA")
        for level in categories.unique():
            mask = categories.eq(level).to_numpy()
            ax.scatter(scores[mask, 0], scores[mask, 1], label=str(level))
        ax.legend(title=color_by, bbox_to_anchor=(1.02, 1), loc="upper left")
    for sample, score in zip(metadata.index, scores):
        ax.annotate(sample, score, fontsize=7, alpha=0.6)
    ax.set_xlabel(f"PC1 ({variance[0]:.0%} of variation)")
    ax.set_ylabel(f"PC2 ({variance[1]:.0%} of variation)")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def save_views(counts, metadata, fields, folder, label):
    scores, variance = fit_pca(counts)
    destination = RESULTS_DIR / folder
    destination.mkdir(parents=True, exist_ok=True)
    for field in fields:
        plot_pca(scores, variance, metadata, field, f"{label} (n={len(metadata)})",
                 destination / f"{field}.png")
    print(f"{label}: {len(metadata)} samples, {len(counts):,} genes, "
          f"PC1 {variance[0]:.1%}, PC2 {variance[1]:.1%}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--subsets", action="store_true",
                        help="Also refit PCA in the preterm and EOPE/iPTB subsets")
    args = parser.parse_args()
    counts, metadata = load_inputs()
    save_views(counts, metadata, OVERVIEW_FIELDS, "all_villous", "All villous samples")

    if args.subsets:
        # OpenBio extension. Select pregnancies from the same released count table.
        primary = metadata[metadata["group"].isin(["EOPE", "iPTB"])]
        control_ages = primary.loc[primary["group"].eq("iPTB"), "ga_delivery_weeks"]
        low, high = control_ages.min(), control_ages.max()
        cohorts = (
            ("preterm", "Preterm samples", metadata[metadata["preterm_delivery"].eq("yes")]),
            ("eope_vs_iptb", "EOPE without FGR vs. iPTB", primary),
            ("shared_age", f"EOPE vs. iPTB at {low:.3f}–{high:.3f} weeks",
             primary[primary["ga_delivery_weeks"].between(low, high)]),
        )
        for folder, label, cohort in cohorts:
            save_views(counts.loc[:, cohort.index], cohort, ("group", "ga_delivery_weeks"), folder, label)


if __name__ == "__main__":
    main()

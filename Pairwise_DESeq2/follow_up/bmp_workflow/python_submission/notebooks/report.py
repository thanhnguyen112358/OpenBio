"""Render the BMP analysis figures from exported result tables."""

import argparse
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, MaxNLocator


INK = "#182D3D"
MUTED = "#4D6271"
TEAL = "#005E68"
RUST = "#A14325"
GRID = "#E2EBEE"
NEUTRAL = "#8398A4"
GROUPS = {"FGR": ("Isolated FGR", TEAL), "EOPE+FGR": ("EOPE + FGR", RUST)}
FITS = {
    "openbio_primary": "OpenBio",
    "openbio_overlap": "OpenBio (shared age)",
    "awamleh_primary": "Awamleh",
}
STYLE = {
    "font.family": "DejaVu Sans",
    "font.size": 11,
    "axes.labelcolor": INK,
    "axes.titlecolor": INK,
    "axes.titleweight": "bold",
    "xtick.color": INK,
    "ytick.color": INK,
    "text.color": INK,
    "axes.edgecolor": GRID,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "axes.axisbelow": True,
    "grid.color": GRID,
    "grid.linewidth": 0.65,
    "legend.frameon": False,
}


def read(output, name):
    return pd.read_csv(output / f"{name}.csv")


def flag(column):
    return column.isin([True, "True", "TRUE"])


def heading(fig, title, subtitle, caption):
    fig.text(.065, .965, title, fontsize=18, va="top", weight="bold")
    fig.text(.065, .902, subtitle, fontsize=11, va="top", color=MUTED)
    fig.text(.065, .038, caption, fontsize=9, va="bottom", color=MUTED, linespacing=1.5)


def save(fig, output, name):
    fig.savefig(output / "figures" / f"{name}.png", dpi=170, facecolor="white")
    plt.close(fig)


def age_support(output, membership):
    cohort = read(output, "cohort_manifest")
    full = membership.loc[membership["comparison"].eq("openbio_primary")].merge(
        cohort[["cohort", "sample_id", "qc_flags"]],
        on=["cohort", "sample_id"], validate="one_to_one",
    )
    overlap = membership.loc[membership["comparison"].eq("openbio_overlap")]
    lo, hi = overlap["ga_delivery_weeks"].agg(["min", "max"])
    overlap_counts = overlap["group"].value_counts()

    fig, ax = plt.subplots(figsize=(11.4, 5.1))
    fig.subplots_adjust(left=.165, right=.955, top=.81, bottom=.235)
    ax.axvspan(lo, hi, color="#DBEEEE", zorder=0)
    ax.text((lo + hi) / 2, 1.44, "Shared delivery range", color=TEAL, ha="center")
    labels = []
    for i, (group, (name, color)) in enumerate(GROUPS.items()):
        donors = full.loc[full["group"].eq(group)].sort_values(["ga_delivery_weeks", "donor_id"])
        labels.append(f"{name} ({len(donors)})")
        bucket = donors["ga_delivery_weeks"].diff().gt(.38).cumsum()
        offsets = np.array([0., .23, -.23, -.37])
        heights = i + offsets[donors.groupby(bucket).cumcount().to_numpy() % offsets.size]
        for donor, y in zip(donors.itertuples(), heights):
            bp_flag = has_blood_pressure_flag(donor.qc_flags)
            ax.scatter(donor.ga_delivery_weeks, y, marker="^" if bp_flag else "o",
                       color=color, s=64, zorder=3)
            ax.annotate(donor.donor_id, (donor.ga_delivery_weeks, y), xytext=(0, 9),
                        textcoords="offset points", ha="center", fontsize=9, color=color)
    ax.set(ylim=(-.43, 1.66), xlabel="Gestational age at delivery (weeks)", yticks=[0, 1])
    ax.set_yticklabels(labels)
    ax.yaxis.grid(False)
    ax.tick_params(axis="y", length=0)
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    legend = [
        Line2D([], [], linestyle="none", marker="o", color=MUTED, label="Released diagnosis"),
        Line2D([], [], linestyle="none", marker="^", color=MUTED, label="Blood pressure flag"),
    ]
    ax.legend(handles=legend, loc="upper center", bbox_to_anchor=(.5, -.20), ncol=2, fontsize=9)
    heading(fig, "Delivery-age overlap in the OpenBio comparison",
            f"OpenBio villi  |  {lo:.3f} to {hi:.3f} weeks  |  "
            f"{overlap_counts['EOPE+FGR']} combined cases and {overlap_counts['FGR']} isolated FGR",
            "Each point is one pregnancy. Triangles mark blood pressure flags on the released diagnosis.")
    save(fig, output, "age_support")


def has_blood_pressure_flag(value):
    text = str(value).lower()
    return any(term in text for term in ("bp", "sbp", "dbp", "blood pressure"))


def headline(output, membership):
    effects = read(output, "bmp_effects")
    core = effects.loc[effects["score"].eq("BMP27")].set_index("comparison").loc[list(FITS)]
    scores = read(output, "bmp_donor_scores")
    lead = "Placental BMP-related RNA was higher with early preeclampsia"
    if (core.loc[["openbio_primary", "awamleh_primary"], "effect"] <= 0).any():
        lead = "Placental BMP-related RNA in pregnancies with growth restriction"
    fig = plt.figure(figsize=(12.4, 8.4))
    gs = fig.add_gridspec(2, 2, left=.20, right=.96, bottom=.135, top=.845,
                          hspace=.64, wspace=.28, height_ratios=[1.04, 1])
    rng = np.random.default_rng(1049)
    for col, fit in enumerate(("openbio_primary", "awamleh_primary")):
        ax = fig.add_subplot(gs[0, col])
        donors = membership.loc[membership["comparison"].eq(fit)].merge(
            scores[["cohort", "donor_id", "bmp27"]],
            on=["cohort", "donor_id"], validate="one_to_one",
        )
        labels = []
        for i, (group, (name, color)) in enumerate(GROUPS.items()):
            values = donors.loc[donors["group"].eq(group), "bmp27"]
            n = len(values)
            labels.append(f"{name}\nn = {n}")
            ax.scatter(i + rng.uniform(-.12, .12, n), values, s=35, color=color,
                       alpha=.78, edgecolors="white", linewidths=.3, zorder=3)
            ax.hlines(values.mean(), i - .18, i + .18, color=INK, linewidth=2.4, zorder=4)
        ax.set_title(FITS[fit], loc="left", pad=10)
        ax.set(xlim=(-.45, 1.45), xticks=np.arange(len(labels)),
               ylabel="Average BMP RNA score (log₂ units)")
        ax.set_xticklabels(labels)
        ax.xaxis.grid(False)
        ax.tick_params(axis="x", length=0)

    ax = fig.add_subplot(gs[1, :])
    y = np.arange(len(core))[::-1]
    ax.axvline(0., color=NEUTRAL, linestyle=(0, (4, 4)), linewidth=1.2)
    ax.hlines(y, core["ci_low"], core["ci_high"], color=TEAL, linewidth=2)
    ax.scatter(core["effect"], y, s=55, color=TEAL, zorder=3)
    for value, row in zip(core["effect"], y):
        ax.annotate(f"{value:+.3f}", (value, row), xytext=(0, 9),
                    textcoords="offset points", ha="center", fontsize=10, weight="bold")
    ax.set_title("Difference between group means", loc="left", pad=10)
    ax.set(yticks=y, ylim=(-.47, len(core) - .42),
           xlabel="EOPE + FGR minus isolated FGR (log₂ units)")
    ax.set_yticklabels(FITS.values())
    ax.margins(x=.08)
    ax.yaxis.grid(False)
    ax.tick_params(axis="y", length=0)
    heading(fig, lead, "Two cohorts of pregnancies with fetal growth restriction",
            "Points above are pregnancies and bars mark group means. Cohorts use separate score scales.\n"
            "Lines below show approximate 95% confidence intervals.")
    save(fig, output, "headline")


def enrichment(output):
    ranks = read(output, "rank_plot")
    summary = read(output, "branch_enrichment")
    bmp = summary.loc[summary["pathway"].eq("REACTOME_SIGNALING_BY_BMP")].set_index("fit")
    fig, axes = plt.subplots(1, len(FITS), figsize=(12.4, 5.55), sharex=True, sharey=True)
    fig.subplots_adjust(left=.085, right=.966, top=.745, bottom=.28, wspace=.13)
    for ax, (fit, label) in zip(axes, FITS.items()):
        table = ranks.loc[ranks["fit"].eq(fit)].sort_values("rank")
        hits = table.loc[flag(table["in_bmp"])]
        ax.axhline(0, color=NEUTRAL, linewidth=.8)
        ax.plot(table["rank"], table["running_score"], color=TEAL, linewidth=1.8)
        ax.vlines(hits["rank"], 0, .085, transform=ax.get_xaxis_transform(),
                  color=np.where(flag(hits["shared_edge"]), RUST, NEUTRAL), linewidth=1.)
        for receptor in table.loc[table["symbol"].eq("BMPR2")].itertuples():
            ax.scatter(receptor.rank, receptor.running_score, color=RUST, s=22, zorder=3)
            ax.annotate("BMPR2", (receptor.rank, receptor.running_score), xytext=(5, 5),
                        textcoords="offset points", color=RUST, fontsize=9, weight="bold")
        line = bmp.loc[fit]
        ax.set_title(f"{label}\nNES {line['NES']:.2f}  |  FDR {line['fdr']:.3g}", fontsize=11, pad=10)
        ax.xaxis.set_major_locator(MaxNLocator(nbins=4, integer=True))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:,.0f}"))
        ax.margins(y=.18)
    axes[0].set_ylabel("Running enrichment score")
    fig.text(.51, .188, "Gene rank by signed Wald statistic (higher with EOPE + FGR at left)",
             ha="center", fontsize=10)
    fig.legend(handles=[Line2D([], [], color=RUST, label="Shared leading edge"),
                        Line2D([], [], color=NEUTRAL, label="Other BMP member")],
               loc="lower center", bbox_to_anchor=(.51, .097), ncol=2, fontsize=9)
    first = bmp.loc["openbio_primary"]
    heading(fig, "Where BMP genes fall in the differential-expression ranks",
            f"{int(first['rank_size']):,} genes in the common rank universe, "
            f"including {int(first['size'])} Reactome BMP genes",
            "Ticks show Reactome BMP genes. NES is normalized enrichment. "
            "FDR adjusts for testing the three Reactome branches per fit.")
    save(fig, output, "enrichment")


def atlas(output):
    table = read(output, "atlas_summary")
    classes = table["cell_class"].drop_duplicates().tolist()
    preferred = ["BMPR2", "BMPR1B", "SMAD9", "SMAD4", "GREM2", "INHA", "INHBA"]
    observed = set(table["symbol"])
    genes = [g for g in preferred if g in observed] + sorted(observed.difference(preferred))
    table["x"] = table["symbol"].map({g: i for i, g in enumerate(genes)})
    table["y"] = table["cell_class"].map({c: i for i, c in enumerate(classes)})
    status = table["status"].eq("reported")
    available = table.loc[status]
    unavailable = table.loc[~status]
    labels = table.drop_duplicates("cell_class").set_index("cell_class")["donors_eligible"]
    fig, ax = plt.subplots(figsize=(12, 8.6))
    fig.subplots_adjust(left=.235, right=.73, bottom=.14, top=.832)
    ax.scatter(unavailable["x"], unavailable["y"], marker="x", color=NEUTRAL, s=22, linewidths=.9)
    cmap = LinearSegmentedColormap.from_list("atlas", ["#D1DDE1", "#00535F"])
    palette = Normalize(0, max(1., np.log10(available["median_cpm"].max() + 1)))
    points = ax.scatter(available["x"], available["y"], c=np.log10(available["median_cpm"] + 1),
                        s=12 + 116 * available["fraction_ge10"], cmap=cmap, norm=palette,
                        linewidths=0, zorder=3)
    ax.set(xticks=range(len(genes)), yticks=range(len(classes)), xlabel="Reactome BMP gene")
    ax.set_xticklabels(genes, weight="bold")
    donor_labels = [f"{pretty_class(c)}  ({int(labels[c])})" for c in classes]
    ax.set_yticklabels(donor_labels, fontsize=10)
    ax.set_xlim(-.5, len(genes) - .5)
    ax.set_ylim(len(classes) - .45, -.55)
    ax.xaxis.grid(False)
    ax.tick_params(axis="both", length=0)
    scale_ax = fig.add_axes([.769, .535, .018, .248])
    scale = fig.colorbar(points, cax=scale_ax)
    marks = [v for v in (0, 1, 10, 100, 1000) if np.log10(v + 1) <= palette.vmax]
    scale.set_ticks(np.log10(np.array(marks) + 1), labels=[f"{v:,}" for v in marks])
    scale.outline.set_visible(False)
    fig.text(.769, .821, "Median RNA\n(counts per million)", va="center", fontsize=10)
    handles = [Line2D([], [], linestyle="none", marker="o", color=TEAL,
                      markersize=np.sqrt(12 + 116 * val), label=f"{val:.0%}")
               for val in (0., .5, 1.)]
    fig.legend(handles=handles, title="Donors with at least 10 reads", loc="upper left",
               bbox_to_anchor=(.757, .465), fontsize=10, title_fontsize=10, labelspacing=1.1)
    heading(fig, "BMP-related transcripts occur in several placental cell classes",
            "Admati donor-by-cell-class reference  |  Parentheses show the number of eligible donors",
            "Each contributing donor-class has at least 20 cells. Crosses mark unavailable summaries.\n"
            "The reference includes preeclampsia and controls and has no isolated-FGR comparison arm.")
    save(fig, output, "atlas")


def pretty_class(value):
    text = value.replace("_", " ").capitalize()
    return text.replace("Nk cell", "NK cell").replace(" nk cell", " NK cell")


def pca(output):
    table = read(output, "pca_scores")
    pca_keys = {"group": "Diagnosis", "fetal_sex": "Fetal sex", "inferred_round": "Inferred sequencing round"}
    colors = {"EOPE+FGR": RUST, "FGR": TEAL, "EOPE": "#7666A4", "iPTB": "#5E7895", "TB": "#805A3C"}
    fig, axes = plt.subplots(1, len(pca_keys), figsize=(14.8, 5.6), sharex=True, sharey=True)
    fig.subplots_adjust(left=.08, right=.964, top=.825, bottom=.235, wspace=.13)
    axis_labels = []
    for ax, (key, title) in zip(axes, pca_keys.items()):
        positioned = []
        for i, (label, points) in enumerate(table.groupby(key, sort=True)):
            color = colors.get(label, [TEAL, RUST, "#7666A4"][i % 3])
            ax.scatter(points["PC1"], points["PC2"], s=38, color=color,
                       label=str(label).replace("_", " "), alpha=.88)
            positioned.extend((point.donor_id, point.PC1, point.PC2, color) for point in points.itertuples())
        axis_labels.append(positioned)
        ax.set_title(title, loc="left", pad=8, fontsize=12)
        ax.set_xlabel(f"PC1 ({table['pc1_percent'].iloc[0]:.1f}%)")
        ax.legend(loc="upper center", bbox_to_anchor=(.5, -.18), ncol=5 if key == "group" else 2,
                  fontsize=8, columnspacing=.6, handletextpad=.2)
        ax.margins(.12)
    axes[0].set_ylabel(f"PC2 ({table['pc2_percent'].iloc[0]:.1f}%)")
    fig.canvas.draw()
    for ax, positioned in zip(axes, axis_labels):
        annotate_pca(ax, positioned)
    heading(fig, "OpenBio villous gene-expression profiles",
            "PCA on the 500 most variable genes after variance stabilization",
            "Each point is one pregnancy. Inferred sequencing round comes from the original study methods.")
    save(fig, output, "pca")


def annotate_pca(ax, positioned):
    # Place each short donor ID on the first clear side of its point.
    candidates = [(0, 9), (0, -9), (14, 0), (-14, 0)]
    labels = []
    coords = [(ax.transData.transform((x, y)) * 72 / ax.figure.dpi) for _, x, y, _ in positioned]
    for (name, x, y, color), (xp, yp) in zip(positioned, coords):
        for dx, dy in candidates:
            anchor = (xp + dx, yp + dy)
            if any(abs(anchor[0] - px) < 10 and abs(anchor[1] - py) < 6 for px, py in coords):
                continue
            if any(abs(anchor[0] - px) < 16 and abs(anchor[1] - py) < 8 for px, py in labels):
                continue
            ax.annotate(name, (x, y), xytext=(dx, dy), textcoords="offset points",
                        fontsize=7, ha="center", va="center", color=color)
            labels.append(anchor)
            break


def volcano(output):
    columns = ["fit", "pvalue", "symbol", "log2FoldChange", "in_bmp"]
    source = pd.read_csv(output / "bmp_gene_effects.tsv.gz", sep="\t", usecols=columns)
    finite = np.isfinite(source["pvalue"]) & np.isfinite(source["log2FoldChange"])
    table = source.loc[source["fit"].eq("openbio_primary") & finite].copy()
    table["neglog"] = -np.log10(table["pvalue"].clip(lower=np.finfo(float).tiny))
    hits = table.loc[flag(table["in_bmp"])]
    fig, ax = plt.subplots(figsize=(10.8, 5.2))
    fig.subplots_adjust(left=.105, right=.955, top=.81, bottom=.235)
    ax.scatter(table["log2FoldChange"], table["neglog"], color=NEUTRAL, alpha=.22, s=9, linewidths=0)
    ax.scatter(hits["log2FoldChange"], hits["neglog"], color=TEAL, s=30, linewidths=0)
    for value in hits.loc[hits["symbol"].eq("BMPR2")].itertuples():
        ax.annotate("BMPR2", (value.log2FoldChange, value.neglog), xytext=(7, 5),
                    textcoords="offset points", color=TEAL, weight="bold", fontsize=9)
    ax.set(xlabel="Estimated log₂ fold change (positive is higher with EOPE + FGR)",
           ylabel="−log₁₀ of Wald p-value")
    heading(fig, "Gene-level differences in OpenBio villi",
            "EOPE + FGR versus isolated FGR  |  Sex-adjusted PyDESeq2",
            "Teal highlights Reactome BMP genes. The vertical scale uses unadjusted p-values.\n"
            "The exported table contains gene-wide false discovery rates and Wald intervals.")
    save(fig, output, "volcano")


def render(output: Path) -> None:
    """Write the six analysis PNGs beneath the result directory."""
    output = Path(output)
    (output / "figures").mkdir(exist_ok=True)
    membership = read(output, "comparison_membership")
    with plt.rc_context(STYLE):
        age_support(output, membership)
        headline(output, membership)
        enrichment(output)
        atlas(output)
        pca(output)
        volcano(output)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, help="workflow result directory")
    render(parser.parse_args().output)

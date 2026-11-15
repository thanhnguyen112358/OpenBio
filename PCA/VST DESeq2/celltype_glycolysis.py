"""Which placental cell classes carry the glycolytic-shift signature in early-onset preeclampsia?

Exploratory. Ackerman et al. (2023) described a bulk "glycolytic shift" in a subset of EOPE
placentas, indexed by the LDHA/LDHB and PDK1/GOT1 expression ratios, and could not assign it to
a cell type. This script computes the same two ratios, plus Hallmark Hypoxia and Glycolysis
scores, in the released Admati et al. (2023) donor x cell-class pseudobulk profiles.

Design choices fixed before looking at results:
  - a profile needs >= 20 cells; a class is analyzed only if every clinical group keeps >= 3 donors
  - counts -> CPM within each profile; ratios are log2((CPM_a + 1) / (CPM_b + 1))
  - pathway score = mean over genes of the within-class z-score of log2(CPM + 1); genes need
    mean CPM >= 10 in that class
  - effect = mean(PE) - mean(control), early and late onset separately; p-value from the exact
    permutation of donor labels (early: 286 or fewer relabelings, so p cannot go below ~0.0035)
  - primary contrast: early PE vs. early control. Late onset is a specificity check.

Caveat built into the data: all 3 early controls delivered vaginally and 9 of 10 early PE by
cesarean, so labor is confounded with early PE. Expression is not metabolic flux.

Run from the OpenBio root:
  .\\.venv\\Scripts\\python.exe "PCA/VST DESeq2/celltype_glycolysis.py"
"""

from itertools import combinations

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import DATA_DIR, RESULTS_DIR
from pathway_analysis import GMT, read_gmt

PB = DATA_DIR / "extension" / "pseudobulk"
OUT_DIR = RESULTS_DIR / "celltype_glycolysis"
MIN_CELLS, MIN_DONORS, MIN_CPM = 20, 3, 10
RATIOS = {"LDHA/LDHB": ("LDHA", "LDHB"), "PDK1/GOT1": ("PDK1", "GOT1")}
SETS = {"Hypoxia score": "HALLMARK_HYPOXIA", "Glycolysis score": "HALLMARK_GLYCOLYSIS"}
METRICS = [*RATIOS, *SETS]
GROUPS = ["control_early", "PE_early", "control_late", "PE_late"]
COLORS = {"control_early": "tab:blue", "PE_early": "tab:orange", "control_late": "#9ecae1", "PE_late": "#fdd0a2"}
CLASS_LABELS = {
    "syncytiotrophoblast": "Syncytio-\ntrophoblast", "villous_cytotrophoblast": "Cyto-\ntrophoblast",
    "hofbauer_cell": "Hofbauer\n(fetal macrophage)", "macrophage": "Macrophage\n(maternal)",
    "monocyte": "Monocyte", "endothelial_cell": "Endothelial", "stromal_cell": "Stromal",
    "pericyte": "Pericyte", "smooth_muscle_cell": "Smooth\nmuscle", "nk_cell": "NK cell", "t_cell": "T cell"}


def exact_permutation(values, is_case):
    """Mean difference (case - control) and two-sided exact permutation p-value."""
    values, is_case = np.asarray(values, float), np.asarray(is_case, bool)
    n, k = len(values), int(is_case.sum())
    observed = values[is_case].mean() - values[~is_case].mean()
    total, count, extreme = values.sum(), 0, 0
    for idx in combinations(range(n), k):
        s = values[list(idx)].sum()
        diff = s / k - (total - s) / (n - k)
        count += 1
        extreme += abs(diff) >= abs(observed) - 1e-12
    return observed, extreme / count, count


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    counts = pd.read_csv(PB / "track01_pseudobulk_counts.csv", index_col=0)
    meta = pd.read_csv(PB / "track01_pseudobulk_metadata.csv").set_index("pseudobulk_id")
    meta = meta.loc[counts.columns]
    keep = meta[meta["n_cells"] >= MIN_CELLS]
    donors_per_group = keep.groupby(["cell_class", "clinical_group"])["donor_id"].nunique().unstack(fill_value=0)
    classes = [c for c in CLASS_LABELS if c in donors_per_group.index
               and (donors_per_group.loc[c, GROUPS] >= MIN_DONORS).all()]
    donors_per_group.to_csv(OUT_DIR / "donors_per_class_and_group.csv")
    sets = read_gmt(GMT)

    rows = []
    for cell_class in classes:
        m = keep[keep["cell_class"] == cell_class]
        cpm = counts[m.index] / counts[m.index].sum(axis=0) * 1e6
        log_cpm = np.log2(cpm + 1)
        table = m[["donor_id", "clinical_group", "n_cells", "delivery_mode", "ga_delivery_weeks", "fgr_status"]].copy()
        table["cell_class"] = cell_class
        for name, (a, b) in RATIOS.items():
            table[name] = np.log2((cpm.loc[a] + 1) / (cpm.loc[b] + 1))
            for gene in (a, b):
                table[f"cpm_{gene}"] = cpm.loc[gene]
        for name, set_name in SETS.items():
            genes = [g for g in sets[set_name] if g in cpm.index and cpm.loc[g].mean() >= MIN_CPM]
            z = log_cpm.loc[genes].sub(log_cpm.loc[genes].mean(axis=1), axis=0)
            z = z.div(log_cpm.loc[genes].std(axis=1, ddof=1), axis=0)
            table[name] = z.mean(axis=0)
            table[f"n_genes_{name}"] = len(genes)
        rows.append(table)
    profiles = pd.concat(rows)
    profiles.to_csv(OUT_DIR / "profiles.csv", float_format="%.5g")

    effects = []
    for cell_class in classes:
        sub = profiles[profiles["cell_class"] == cell_class]
        for onset in ("early", "late"):
            pair = sub[sub["clinical_group"].isin([f"control_{onset}", f"PE_{onset}"])]
            is_case = pair["clinical_group"].eq(f"PE_{onset}").to_numpy()
            for metric in METRICS:
                diff, p, n_perm = exact_permutation(pair[metric], is_case)
                controls = pair.loc[~is_case, metric].to_numpy()
                cases = pair.loc[is_case, metric].to_numpy()
                # Leave one control out: the range of the effect when each control is dropped.
                loo = [cases.mean() - np.delete(controls, i).mean() for i in range(len(controls))]
                effects.append({
                    "cell_class": cell_class, "onset": onset, "metric": metric,
                    "n_PE": int(is_case.sum()), "n_control": int((~is_case).sum()),
                    "mean_PE": cases.mean(), "mean_control": controls.mean(), "difference": diff,
                    "exact_permutation_p": p, "n_relabelings": n_perm,
                    "min_possible_p": 2 / n_perm if n_perm > 1 else 1.0,
                    "leave_one_control_out_min": min(loo), "leave_one_control_out_max": max(loo),
                    "PE_donors_above_every_control": int((cases > controls.max()).sum()),
                    "PE_donors_below_every_control": int((cases < controls.min()).sum()),
                })
    effects = pd.DataFrame(effects)
    effects.to_csv(OUT_DIR / "effects_by_class.csv", index=False, float_format="%.5g")

    delivery = (profiles.drop_duplicates("donor_id").groupby(["clinical_group", "delivery_mode"]).size()
                .unstack(fill_value=0))
    delivery.to_csv(OUT_DIR / "delivery_mode_by_group.csv")

    plot_metric(profiles, effects, classes, "LDHA/LDHB", "log2 ratio LDHA / LDHB\n(higher = toward lactate production)",
                OUT_DIR / "ldha_ldhb_by_class.png")
    plot_metric(profiles, effects, classes, "PDK1/GOT1", "log2 ratio PDK1 / GOT1", OUT_DIR / "pdk1_got1_by_class.png")
    plot_metric(profiles, effects, classes, "Glycolysis score", "Hallmark Glycolysis score\n(mean z within class)",
                OUT_DIR / "glycolysis_score_by_class.png")
    plot_metric(profiles, effects, classes, "Hypoxia score", "Hallmark Hypoxia score\n(mean z within class)",
                OUT_DIR / "hypoxia_score_by_class.png")
    plot_effects(effects, classes, OUT_DIR / "effects_summary.png")

    pd.set_option("display.width", 250)
    print("classes analyzed:", classes)
    print(delivery.to_string())
    show = effects.copy()
    show["cell"] = (show["difference"].round(2).astype(str) + " (p=" + show["exact_permutation_p"].round(3).astype(str)
                    + "; LOO " + show["leave_one_control_out_min"].round(2).astype(str) + ".."
                    + show["leave_one_control_out_max"].round(2).astype(str) + ")")
    for onset in ("early", "late"):
        print(f"\n== {onset} onset: PE minus control ==")
        print(show[show["onset"] == onset].pivot(index="cell_class", columns="metric", values="cell")
              .loc[classes, METRICS].to_string())
    print("\ndetection (median CPM across profiles):")
    print(profiles.groupby("cell_class")[[c for c in profiles if c.startswith("cpm_")]].median().round(1).loc[classes].to_string())


def plot_metric(profiles, effects, classes, metric, ylabel, output_file):
    fig, ax = plt.subplots(figsize=(1.25 * len(classes) + 1.5, 4.8))
    offsets = dict(zip(GROUPS, [-0.3, -0.1, 0.1, 0.3]))
    for i, cell_class in enumerate(classes):
        sub = profiles[profiles["cell_class"] == cell_class]
        for group in GROUPS:
            values = sub.loc[sub["clinical_group"] == group, metric]
            x = i + offsets[group] + np.linspace(-0.05, 0.05, len(values))
            ax.scatter(x, values, s=22, color=COLORS[group], edgecolors="white", linewidths=0.4, zorder=3,
                       label=group.replace("_", " ") if i == 0 else None)
            ax.hlines(values.mean(), i + offsets[group] - 0.08, i + offsets[group] + 0.08, color="black", lw=1.2, zorder=4)
        row = effects[(effects.cell_class == cell_class) & (effects.onset == "early") & (effects.metric == metric)].iloc[0]
        ax.annotate(f"Δ={row['difference']:+.2f}\np={row['exact_permutation_p']:.3f}", (i - 0.2, 1.0),
                    xycoords=("data", "axes fraction"), ha="center", va="bottom", fontsize=7)
        if i:
            ax.axvline(i - 0.5, color="lightgrey", lw=0.6)
    ax.set_xticks(range(len(classes)), [CLASS_LABELS[c] for c in classes], fontsize=8)
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=7, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.14), frameon=False)
    ax.set_title(f"{metric} by placental cell class (Admati pseudobulk; each dot = one donor; bars = means)\n"
                 "Δ and exact permutation p are for early PE vs. early control (n = 10 vs. 3 where complete)\n",
                 fontsize=9, pad=22)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_effects(effects, classes, output_file):
    fig, axes = plt.subplots(1, len(METRICS), figsize=(3.6 * len(METRICS), 5), sharey=True)
    y = np.arange(len(classes))[::-1]
    for ax, metric in zip(axes, METRICS):
        for onset, dy, color in (("early", 0.15, "tab:orange"), ("late", -0.15, "#bdbdbd")):
            sub = effects[(effects.metric == metric) & (effects.onset == onset)].set_index("cell_class").loc[classes]
            ax.hlines(y + dy, sub["leave_one_control_out_min"], sub["leave_one_control_out_max"], color=color, lw=2)
            ax.scatter(sub["difference"], y + dy, color=color, s=30, zorder=3, label=f"{onset} onset")
        ax.axvline(0, color="black", lw=0.8)
        ax.set_title(metric, fontsize=10)
        ax.set_xlabel("PE minus control")
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_yticks(y, [CLASS_LABELS[c].replace("\n", " ") for c in classes], fontsize=8)
    axes[0].legend(fontsize=7, loc="lower right")
    fig.suptitle("Preeclampsia minus control, by cell class. Dot = difference in means; "
                 "line = range when each control donor is left out.", fontsize=10)
    fig.tight_layout()
    fig.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

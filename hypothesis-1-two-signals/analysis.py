"""Hypothesis 1. Two separable placental RNA patterns at delivery.

FLT1 and Ackerman's ratios (LDHA/LDHB, PDK1/GOT1) track maternal
preeclampsia. A 10-step glycolysis-enzyme score is a broader stress
signal and is present in both villi and decidua.

Predictions, fixed before reading the contrast output of this script:
  1. EOPE and EOPE+FGR sit above iPTB for FLT1 and both ratios.
  2. Isolated FGR does not. FGR and iPTB barely share a delivery-age
     window, so that contrast is directional only.
  3. The enzyme score rises in both preeclampsia groups and is not
     confined to them.
  4. The enzyme rise appears in decidua as well as villi from the
     same pregnancies.

Run from anywhere:
    python hypothesis-1-two-signals/analysis.py
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent / "openbio-datathon-2026" / "track-01" / "data"
FIG = ROOT / "figures"
RES = ROOT / "results"

STEPS = {
    "HK": ["HK1", "HK2", "HK3"],
    "GPI": ["GPI"],
    "PFK": ["PFKL", "PFKM", "PFKP"],
    "ALDO": ["ALDOA", "ALDOB", "ALDOC"],
    "TPI": ["TPI1"],
    "GAPDH": ["GAPDH"],
    "PGK": ["PGK1"],
    "PGAM": ["PGAM1", "PGAM2"],
    "ENO": ["ENO1", "ENO2", "ENO3"],
    "PK": ["PKM", "PKLR"],
}
SIGNALS = ["flt1", "pgf", "ldh_ratio", "pdk_ratio", "enzyme_score"]
GROUP_ORDER = ["iPTB", "FGR", "EOPE", "EOPE+FGR", "TB"]
COLORS = {
    "iPTB": "#4C78A8",
    "FGR": "#F58518",
    "EOPE": "#E45756",
    "EOPE+FGR": "#72B7B2",
    "TB": "#54A24B",
}
LABELS = {
    "flt1": "FLT1 activity",
    "pgf": "PGF activity",
    "ldh_ratio": "LDHA / LDHB",
    "pdk_ratio": "PDK1 / GOT1",
    "enzyme_score": "Glycolysis enzymes",
}


def load_tissue(counts_name, meta_name):
    counts = pd.read_csv(DATA / counts_name, index_col=0)
    meta = pd.read_csv(DATA / meta_name, index_col=0).loc[counts.columns].copy()
    missing = [g for genes in STEPS.values() for g in genes if g not in counts.index]
    missing += [g for g in ["FLT1", "PGF", "LDHA", "LDHB", "PDK1", "GOT1"] if g not in counts.index]
    if missing:
        raise SystemExit(f"Missing genes in {counts_name}: {missing}")
    lib = counts.sum(axis=0)
    log_cpm = np.log2(counts.div(lib, axis=1) * 1e6 + 1)
    cpm = counts.div(lib, axis=1) * 1e6
    meta["flt1"] = log_cpm.loc["FLT1"]
    meta["pgf"] = log_cpm.loc["PGF"]
    meta["ldh_ratio"] = np.log2(cpm.loc["LDHA"] / cpm.loc["LDHB"])
    meta["pdk_ratio"] = np.log2(cpm.loc["PDK1"] / cpm.loc["GOT1"])
    step_scores = []
    for genes in STEPS.values():
        step_scores.append(log_cpm.loc[genes].mean(axis=0))
    meta["enzyme_score"] = pd.concat(step_scores, axis=1).mean(axis=1)
    return meta


def welch(case, control):
    case = np.asarray(case, float)
    control = np.asarray(control, float)
    diff = case.mean() - control.mean()
    na, nb = len(case), len(control)
    va, vb = case.var(ddof=1), control.var(ddof=1)
    se = np.sqrt(va / na + vb / nb)
    df = (va / na + vb / nb) ** 2 / ((va / na) ** 2 / (na - 1) + (vb / nb) ** 2 / (nb - 1))
    tcrit = stats.t.ppf(0.975, df)
    p = stats.ttest_ind(case, control, equal_var=False).pvalue
    return {
        "n_case": na,
        "n_control": nb,
        "difference": diff,
        "ci_low": diff - tcrit * se,
        "ci_high": diff + tcrit * se,
        "p_welch": p,
    }


def contrasts(meta, tissue):
    rows = []
    iPTB = meta[meta.group.eq("iPTB")]
    lo, hi = iPTB.ga_delivery_weeks.min(), iPTB.ga_delivery_weeks.max()
    overlap = meta[meta.ga_delivery_weeks.between(lo, hi)]
    comparisons = {
        "EOPE_vs_iPTB": (meta[meta.group.eq("EOPE")], iPTB),
        "EOPE_vs_iPTB_shared_weeks": (
            overlap[overlap.group.eq("EOPE")],
            overlap[overlap.group.eq("iPTB")],
        ),
        "EOPE+FGR_vs_iPTB": (meta[meta.group.eq("EOPE+FGR")], iPTB),
        "FGR_vs_iPTB": (meta[meta.group.eq("FGR")], iPTB),
        "EOPE+FGR_vs_EOPE": (meta[meta.group.eq("EOPE+FGR")], meta[meta.group.eq("EOPE")]),
    }
    for name, (case, control) in comparisons.items():
        for signal in SIGNALS:
            if len(case) < 2 or len(control) < 2:
                continue
            row = welch(case[signal], control[signal])
            row.update(tissue=tissue, comparison=name, signal=signal)
            rows.append(row)
    return pd.DataFrame(rows), (lo, hi)


def stripplot(ax, meta, signal):
    rng = np.random.default_rng(0)
    for i, group in enumerate(GROUP_ORDER):
        values = meta.loc[meta.group.eq(group), signal]
        x = i + rng.uniform(-0.12, 0.12, len(values))
        ax.scatter(x, values, s=36, color=COLORS[group], zorder=3)
        ax.hlines(values.mean(), i - 0.22, i + 0.22, color="black", lw=1.4, zorder=4)
    ax.set_xticks(range(len(GROUP_ORDER)))
    ax.set_xticklabels(GROUP_ORDER, rotation=30, ha="right")
    ax.set_ylabel(LABELS[signal] + "\n(log2 scale)")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def main():
    FIG.mkdir(exist_ok=True)
    RES.mkdir(exist_ok=True)
    villi = load_tissue("track01_counts_vt.csv", "track01_metadata_vt.csv")
    decidua = load_tissue(
        "extension/decidua/track01_counts_db.csv",
        "extension/decidua/track01_metadata_db.csv",
    )
    villi_stats, window = contrasts(villi, "villous")
    decidua_stats, _ = contrasts(decidua, "decidua")
    stats_table = pd.concat([villi_stats, decidua_stats], ignore_index=True)
    stats_table.to_csv(RES / "contrasts.csv", index=False)

    keep = ["donor_id", "group", "ga_delivery_weeks", "fgr_status", "qc_flags", *SIGNALS]
    villi[keep].to_csv(RES / "villous_donor_scores.csv")
    decidua[keep].to_csv(RES / "decidua_donor_scores.csv")

    paired = villi.reset_index().merge(
        decidua.reset_index(), on="donor_id", suffixes=("_villous", "_decidua")
    )
    paired.to_csv(RES / "paired_donors.csv", index=False)

    fig, axes = plt.subplots(1, 3, figsize=(11, 4.2))
    for ax, signal in zip(axes, ["flt1", "ldh_ratio", "enzyme_score"]):
        stripplot(ax, villi, signal)
    axes[0].set_title("Tracks preeclampsia")
    axes[1].set_title("Ackerman ratio")
    axes[2].set_title("Broader enzyme signal")
    fig.suptitle("Villous placenta at delivery, one dot per pregnancy", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "villous_three_signals.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.2))
    for ax, signal, title in zip(
        axes,
        ["enzyme_score", "flt1"],
        ["Glycolysis enzymes", "FLT1"],
    ):
        for group in ["iPTB", "EOPE"]:
            sub = paired[paired.group_villous.eq(group)]
            for _, row in sub.iterrows():
                ax.plot(
                    [0, 1],
                    [row[f"{signal}_villous"], row[f"{signal}_decidua"]],
                    color=COLORS[group],
                    alpha=0.45,
                    lw=1,
                )
            ax.scatter([0] * len(sub), sub[f"{signal}_villous"], color=COLORS[group], s=28, label=group, zorder=3)
            ax.scatter([1] * len(sub), sub[f"{signal}_decidua"], color=COLORS[group], s=28, zorder=3)
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["Villi\n(fetal side)", "Decidua\n(maternal side)"])
        ax.set_title(title)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
    axes[0].legend(frameon=False, title="Same pregnancies")
    fig.suptitle("Each line is one pregnancy", y=1.02)
    fig.tight_layout()
    fig.savefig(FIG / "paired_villi_decidua.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    for group in GROUP_ORDER:
        sub = villi[villi.group.eq(group)]
        ax.scatter(sub.ga_delivery_weeks, sub.enzyme_score, color=COLORS[group], s=36, label=group)
    ax.axvspan(window[0], window[1], color="#4C78A8", alpha=0.08, label="iPTB age window")
    ax.set_xlabel("Weeks pregnant at delivery")
    ax.set_ylabel("Glycolysis enzyme score")
    ax.legend(frameon=False, ncol=2)
    ax.set_title("Isolated growth restriction is delivered earlier")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG / "enzyme_score_by_delivery_age.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    focus = stats_table[stats_table.comparison.isin(["EOPE_vs_iPTB", "EOPE_vs_iPTB_shared_weeks", "FGR_vs_iPTB"])]
    print(f"iPTB delivery-age window: {window[0]:.3f}–{window[1]:.3f} weeks")
    print(focus.round(3).to_string(index=False))
    print(f"\nWrote figures to {FIG}")


if __name__ == "__main__":
    main()

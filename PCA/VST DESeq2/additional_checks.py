"""Exploratory S24 and fetal-sex sensitivity checks against saved fit A.

Run from any directory with the project Python. Original fits are not rerun.
"""
import json
import platform
from importlib.metadata import version

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from common import RESULTS_DIR, load_inputs, q1_cohorts
from deseq2_fits import run_fit, fit_summary, BENCHMARK_GENES


def main():
    out = RESULTS_DIR / "additional_checks"
    out.mkdir(exist_ok=True)
    counts, metadata, _ = load_inputs()
    cohort = q1_cohorts(metadata)[0]["eope_vs_iptb"]
    pd.crosstab(cohort.group, cohort.fetal_sex).to_csv(out / "sex_balance.csv")
    original = pd.read_csv(RESULTS_DIR / "deseq2/fit_A_full.csv", index_col=0)
    original_hits = original.index[original.passes]
    robust = set(original_hits)
    for name in ["B_age_adjusted", "C_shared_age", "D_without_S22"]:
        old = pd.read_csv(RESULTS_DIR / f"deseq2/fit_{name}.csv", index_col=0)
        robust &= set(old.index[old.passes])
    specs = [
        ("E_without_S24", cohort.drop(index="S24_VT"), "~ group"),
        ("F_sex_adjusted", cohort, "~ fetal_sex + group"),
    ]
    info = {"python": platform.python_version(), "packages": {
        p: version(p) for p in ["pydeseq2", "numpy", "pandas", "scipy"]},
        "baseline": "Saved fit A, not rerun; comparisons use its rounded CSV",
        "exploratory": True, "fits": {}}
    summaries, benchmarks = [], []
    fig, axes = plt.subplots(1, 2, figsize=(11, 5))
    for ax, (name, samples, design) in zip(axes, specs):
        print(f"Running {name}: {len(samples)} donors", flush=True)
        result, cooks, cutoff = run_fit(counts, samples, design, 1)
        result.to_csv(out / f"fit_{name}.csv")
        cooks.to_csv(out / f"cooks_{name}.csv")
        info["fits"][name] = {"design": design, "samples": list(samples.index), "cooks_cutoff": cutoff}
        compare = original[["baseMean", "log2FoldChange", "padj", "passes"]].join(
            result[["log2FoldChange", "padj", "passes"]], lsuffix="_A", rsuffix="_new")
        compare["delta_lfc"] = compare.log2FoldChange_new - compare.log2FoldChange_A
        compare["same_direction"] = np.sign(compare.log2FoldChange_A) == np.sign(compare.log2FoldChange_new)
        compare.to_csv(out / f"comparison_{name}.csv")
        hits = compare.loc[original_hits]
        row = fit_summary(name, samples, design, result)
        row.update({"A_hits_retained": int(hits.passes_new.sum()),
                    "A_hits_same_direction": int(hits.same_direction.sum()),
                    "A_hits_median_abs_delta_lfc": float(hits.delta_lfc.abs().median()),
                    "A_hits_max_abs_delta_lfc": float(hits.delta_lfc.abs().max()),
                    "previous_207_retained": int(result.loc[list(robust), "passes"].sum())})
        summaries.append(row)
        for gene in BENCHMARK_GENES:
            benchmarks.append({"fit": name, "gene": gene, **result.loc[gene].to_dict()})
        measured = compare.loc[compare.baseMean.ge(10)].dropna(subset=["log2FoldChange_A", "log2FoldChange_new"])
        ax.scatter(measured.log2FoldChange_A, measured.log2FoldChange_new, s=3, alpha=.15, color="grey")
        ax.scatter(hits.log2FoldChange_A, hits.log2FoldChange_new, s=8, alpha=.6, color="tab:orange")
        lim = max(ax.get_xlim()[1], ax.get_ylim()[1], abs(ax.get_xlim()[0]), abs(ax.get_ylim()[0]))
        ax.plot([-lim, lim], [-lim, lim], "k--", lw=1)
        ax.set(xlabel="Saved fit A log2 fold change", ylabel="New fit log2 fold change", title=name)
        print(row, flush=True)
    fig.tight_layout()
    fig.savefig(out / "effect_comparison.png", dpi=180)
    pd.DataFrame(summaries).to_csv(out / "summary.csv", index=False)
    pd.DataFrame(benchmarks).to_csv(out / "benchmark_genes.csv", index=False)
    (out / "run_info.json").write_text(json.dumps(info, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

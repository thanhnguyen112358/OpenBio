# S24 and fetal-sex sensitivity checks

Exploratory follow-up to the saved A–D fits. Run with `python additional_checks.py` using the project environment. These checks reuse the saved fit A as baseline; A–D were not rerun. Full results, comparisons, package versions, sample IDs and a scatterplot are in `results/additional_checks/`.

Both new fits use raw counts, the original PyDESeq2 defaults and EOPE-versus-iPTB contrast. The selection rule remains padj < 0.05 and absolute log2 fold change >= 1. This is a zero-effect significance test followed by an estimated-effect cutoff, not a formal test that the true effect exceeds twofold.

| Check | Cohort | Total passing | Original 352 still passing | Original 352 same direction | Original 207 still passing |
|---|---|---:|---:|---:|---:|
| E: omit S24, `~ group` | 8 EOPE / 4 iPTB | 457 | 304 | 352 | 197 |
| F: `~ fetal_sex + group` | 8 EOPE / 5 iPTB | 294 | 260 | 352 | 191 |

The original 207 are the intersection of the four saved fits. These comparisons share donors and are sensitivity checks, not independent validation.

## S24

Omitting S24 preserves the direction of every original hit. Median absolute change among those hits is 0.129 log2 units, but the largest is 3.242 for LEP. Its log2 effect rises from 5.557 to 8.799 (approximately 47-fold to 445-fold). Thus S24 does not create the positive LEP association, but strongly affects its estimated magnitude. All four benchmarks still pass. Keep S24 in the primary analysis; its expression resemblance to cases alone is not grounds for exclusion.

Among the 48 original hits no longer passing, seven have missing padj, 30 have padj >= 0.05, and 15 have absolute estimated effect below 1; these categories overlap. Loss of threshold status alone does not identify its cause.

## Fetal sex

EOPE contains six female and two male pregnancies; iPTB contains two female and three male pregnancies. Both sexes occur in both groups, so the additive group and sex effects can be estimated. The sex-adjusted model preserves all original hit directions, with median absolute change 0.078 log2 units, but retains only 260 of 352 at both thresholds. All four benchmarks remain higher and pass:

| Gene | Sex-adjusted log2 fold change | Adjusted p-value |
|---|---:|---:|
| FLT1 | 3.241 | 1.08e-4 |
| ENG | 2.301 | 7.82e-7 |
| LEP | 6.879 | 5.24e-7 |
| INHBA | 2.449 | 2.94e-3 |

Four genes did not converge (C11orf86, KPRP, LOC389033, MIR93); all have missing padj and none passes. Treat their estimates as unreliable. These flags are retained in the exported tables.

## Interpretation limits

The two checks support persistence of the benchmark directions, while highlighting LEP magnitude sensitivity. They do not rule out confounding. F adjusts for sex alone, not sex and delivery age jointly. Small strata limit precision and this model assumes a common group effect across sexes.

Automatic outlier handling also changes with design: E replaced counts for 86 genes and filtered p-values for 25; F replaced no counts and filtered 21. A originally replaced 69 and filtered 18. Differences therefore include normalization, dispersion, filtering and outlier-rule changes as well as the intended sensitivity. Size factors cannot correct unmeasured library quality confounded with diagnosis.

The baseline CSV is rounded, and package versions may differ from the original run; small numerical differences should not be overinterpreted. Consult `run_info.json` for the new environment. Further biology interpretation should use effect estimates and donor plots alongside threshold counts.

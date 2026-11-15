# Step 3 paired score effects

The inferred same-round comparison included four EOPE-only and five isolated-FGR pregnancies. One paired score per pregnancy was calculated as decidua minus villi. The primary estimate then compared that score between EOPE and FGR. A negative contrast means the EOPE–FGR difference is more negative in decidua.

## Primary estimates

| Outcome | EOPE minus FGR | Approximate 95% CI |
| --- | ---: | ---: |
| Paired C1Q–TREM2, decidua minus villi | −0.25 | −1.70 to +1.20 |
| C1Q–TREM2 in decidua | −1.22 | −2.63 to +0.20 |
| C1Q–TREM2 in villi | −0.96 | −1.70 to −0.23 |
| Broad myeloid score in decidua | −0.58 | −1.38 to +0.21 |
| C1Q–TREM2 minus the broad myeloid score in decidua | −0.63 | −1.53 to +0.26 |

All values are differences in log₂ score units. The point estimates show lower C1Q–TREM2 in both tissues. The paired estimate for a larger decidual decrease is imprecise (Welch two-sided `p = 0.669`). For this paired estimate, the C1Q component was negative (−0.84) and the `TREM2` component was positive (+0.34). The three broad myeloid markers in decidua also varied in direction or magnitude (`CD68` −0.77, `LST1` −1.16, `CSF1R` +0.19). The selectivity contrast is a comparison of bulk RNA scores and cannot distinguish changes in cell composition from within-cell expression.

## Sensitivity to cohort, normalization, and individual donors

| Cohort | PyDESeq2 paired effect (95% CI) | TMM paired effect (95% CI) |
| --- | ---: | ---: |
| Inferred same round (4 vs. 5) | −0.25 (−1.70 to +1.20) | −0.21 (−1.63 to +1.20) |
| Planned age overlap (3 vs. 5) | −0.38 (−1.84 to +1.08) | −0.34 (−1.77 to +1.09) |
| Strict age overlap (3 vs. 4) | +0.08 (−0.96 to +1.12) | +0.11 (−0.90 to +1.13) |
| Full cohort (8 vs. 5) | +0.19 (−1.26 to +1.65) | +0.20 (−1.21 to +1.62) |

- **Leave one pregnancy out.** The PyDESeq2 primary estimate ranged from −0.52 to +0.21 across nine single-donor omissions. TMM ranged from −0.46 to +0.24. Under both methods, omitting `S15` was the sole sign reversal. `S15` has an organizer BP flag and falls below the strict delivery-age window. Omitting the other BP-flagged donor, `S04`, gave −0.52 with PyDESeq2 and −0.46 with TMM.
- **Low villous TREM2 coverage.** `S06` was the sole primary donor below 10 villous `TREM2` reads (five reads). Omitting that pair changed the primary estimate from −0.25 to −0.24 with PyDESeq2 and from −0.21 to −0.19 with TMM.
- **Within-decidua selectivity.** The C1Q–TREM2-minus-myeloid contrast stayed negative across all four cohorts, both normalizations, and every single-donor omission. Its cohort and leave-one-out intervals all spanned zero. The stronger claim of a selective decidual decrease remains unresolved in this cohort.

## Methods and reproduction

The analysis uses the audited cohort memberships and Step 2 PyDESeq2 values. Welch 95% intervals and two-sided p-values are calculated across pregnancies after forming each donor's paired or within-decidua difference. Intervals are approximate, and the secondary intervals and p-values are unadjusted. Sequencing rounds are inferred from the source paper. [Cohort audit](../step01_cohort_audit/README.md) · [Step 2 scores](../step02_deseq2/README.md) · [SciPy Welch interval documentation](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.ttest_ind.html)

TMM uses RNAnorm 2.3.0 and the organizer's genome-wide counts for the same 13 libraries per tissue. The effective library size is the raw library sum times the TMM factor. To preserve the scale of the `+1` pseudocount, the effective sizes are rescaled within tissue to the geometric mean of the Step 2 PyDESeq2 size factors. Counts divided by these rescaled factors enter the same `log2(count + 1)` score functions used in Step 2. Factors are kept fixed through cohort restrictions and donor omissions for both methods. [TMM method](https://link.springer.com/article/10.1186/gb-2010-11-3-r25) · [RNAnorm TMM documentation](https://rnanorm.readthedocs.io/en/stable/generated/rnanorm.TMM.html)

From the workspace root, using the Step 2 environment:

```powershell
& 'Deci_Vili_H1/.venv/Scripts/python.exe' -m pip install rnanorm==2.3.0
& 'Deci_Vili_H1/.venv/Scripts/python.exe' Deci_Vili_H1/scripts/step03_paired_effects.py
```

## Outputs

- [effects.csv](effects.csv) has the score and individual gene contrasts, group means, Welch intervals, and p-values for each cohort and normalization. The `*_gene_*` rows are contrasts of log-transformed normalized counts, alongside the raw-count gene models from Step 2.
- [leave_one_out.csv](leave_one_out.csv) records the omitted donor, coverage and BP flag, revised estimate, and sign change for each score contrast.
- [donor_scores.csv](donor_scores.csv) and [library_scaling.csv](library_scaling.csv) preserve the scored inputs and library scaling factors.
- [Paired donor scores](figures/paired_donor_scores.png) shows every audited EOPE-only and FGR pair. [Tissue and marker effects](figures/tissue_and_marker_effects.png) shows the primary component estimates and their Welch intervals. [Paired effect sensitivity](figures/paired_effect_sensitivity.png) compares cohort definitions, normalization methods, and donor omissions.

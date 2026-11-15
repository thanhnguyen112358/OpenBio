# Exploratory follow-up in the organizer data

Adding the combined EOPE+FGR group lets us compare pregnancies with and without early-onset preeclampsia when both have fetal growth restriction. This comparison reuses the same bulk study, so it provides a check across clinical groups within one cohort. The separate single-cell reference asks how the same genes behave within broad cell classes in EOPE without FGR versus early controls. [Bulk study](https://pmc.ncbi.nlm.nih.gov/articles/PMC10192030/) · [Single-cell follow-up](pseudobulk/README.md)

## Bulk score in EOPE+FGR versus isolated FGR

The combined group supplied seven villous samples and four complete decidua–villus pairs. The four paired cases and five isolated-FGR pregnancies are all assigned to the same inferred sequencing round. The table uses contrasts of normalized log₂ scores, with exploratory Welch 95% intervals available in the results file. [Donor scores](bulk_donor_scores.csv) · [Signed estimates and intervals](bulk_contrasts.csv) · [Cohort audit](../step01_cohort_audit/README.md)

| EOPE+FGR minus isolated FGR | All paired, 4 vs. 5 | Within shared delivery-age range, 2 vs. 4 | Paired without isolated-FGR donor `S15`, 4 vs. 4 |
| --- | ---: | ---: | ---: |
| C1Q–TREM2 in villi | −0.76 | −1.07 | −0.88 |
| C1Q–TREM2 in decidua | −0.79 | +0.34 | −0.44 |
| Paired decidua-minus-villi difference | −0.03 | +1.41 | +0.43 |
| C1Q–TREM2 minus broad myeloid in decidua | −0.19 | +0.38 | −0.03 |

The combined group had lower scores than isolated FGR in both tissues among complete pairs. Its mean was between the EOPE-only and isolated-FGR means in both tissues in that inferred sequencing round. The combined-versus-isolated FGR paired contrast was near zero (95% CI −2.00 to +1.94), and the decidual score relative to the broad markers was imprecise (95% CI −1.03 to +0.66) and changed sign in the smaller age subset. This follow-up gives little support for a preferentially decidual difference shared across EOPE phenotypes.

Including all seven combined villous samples moved the villous estimate from −0.76 to **−0.29** (95% CI −1.11 to +0.54). The three combined cases without decidua therefore matter to the estimated villous magnitude. The paired shared-age filter covers 27.000–30.429 weeks. The two retained combined pregnancies delivered at 27.000 and 27.714 weeks, earlier than the four retained isolated-FGR pregnancies at 29.000–30.429 weeks. The age restriction still leaves age differences between the groups. [Scores by delivery age](scores_by_delivery_age.png)

## Cell-class context in the separate single-cell cohort

In the organizer's donor-level pseudobulk, EOPE without FGR averaged lower C1Q–TREM2 scores than early controls in broad macrophages (−0.52, 3 vs. 3 donors) and higher scores in Hofbauer cells (+0.67, 4 vs. 3 donors). The signs held in the small delivery-age subsets. All eligible EOPE-only donors delivered by caesarean and all early controls vaginally. The reference has no isolated-FGR arm, so these estimates describe how the panel can behave across broad cell classes in a different comparison. They leave disease effects and subtype composition unresolved. [Cell-class analysis and exact donor counts](pseudobulk/README.md)

## Reproduce the bulk comparison

The script fits the PyDESeq2 median-ratio reference from the same 13 Step 2 pregnancies separately for villi and decidua. It projects each newly included library onto that reference and checks that the original factors and scores reproduce. The gene weights and log₂ transform are those used in Step 2. Decidua-based estimates use observed complete pairs, with the villous-only analysis used to assess the effect of including the remaining samples. The reported cohorts reuse donors and are exploratory. [Step 2 normalization](../step02_deseq2/README.md) · [Bulk analysis script](../../scripts/followup_combined_fgr.py) · [PyDESeq2 normalization reference](https://pydeseq2.readthedocs.io/en/latest/generated/pydeseq2.preprocessing.deseq2_norm_fit.html)

```powershell
& 'Deci_Vili_H1/.venv/Scripts/python.exe' Deci_Vili_H1/scripts/followup_combined_fgr.py
```

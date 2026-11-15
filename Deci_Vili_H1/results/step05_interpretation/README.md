# Lower C1Q–TREM2 RNA scores in both tissues at delivery in EOPE versus isolated FGR

We asked whether the difference in C1Q–TREM2 RNA scores between early-onset preeclampsia without fetal growth restriction (EOPE-only) and isolated fetal growth restriction (FGR) was greater in decidua than in villi. The group mean was lower in EOPE-only in both tissues across all four cohort definitions and both normalization methods. These comparisons reuse pregnancies from the same cohort. The predicted extra decrease in decidua was sensitive to who was included. Within decidua, the difference in the candidate score was consistently more negative than the difference in the broad myeloid score, giving us a lead to test in an independent cohort. [Score effects across cohorts and normalization methods](../step03_paired_effects/effects.csv)

| Finding in the primary cohort of 4 EOPE-only and 5 isolated-FGR pregnancies | EOPE minus FGR (log₂ score units) | Approximate 95% CI |
| --- | ---: | ---: |
| C1Q–TREM2 score in decidua | −1.22 | −2.63 to +0.20 |
| C1Q–TREM2 score in villi | −0.96 | −1.70 to −0.23 |
| Paired contrast for the extra decrease in decidua | −0.25 | −1.70 to +1.20 |
| Decidual C1Q–TREM2 minus broad myeloid score | −0.63 | −1.53 to +0.26 |

The table uses PyDESeq2 values from the inferred same sequencing round. The paired estimate moved to **+0.08** in the strict delivery-age subset and **+0.21** when the FGR donor `S15` alone was omitted from the primary cohort. TMM gave the same pattern of sign reversals. The paired C1Q-only component pointed toward an extra decidual decrease (−0.84), while the `TREM2` component pointed in the other direction (+0.34). The paired raw-count gene models in the primary cohort showed the same split in direction. [Paired effect and donor sensitivity](../step03_paired_effects/figures/paired_effect_sensitivity.png) · [Raw-count gene estimates](../step02_deseq2/gene_effects.csv)

## What the biology can explain

The within-decidua contrast stayed negative across cohort definitions, normalizations, and single-donor omissions. Its intervals all spanned zero. Broad myeloid markers also varied individually: `CD68` and `LST1` had negative estimates, and `CSF1R` had a small positive estimate. The contrast therefore compares two bulk RNA scores. Changes in cell abundance, macrophage subtype mix, and transcription per cell could each contribute. [Tissue and marker effects](../step03_paired_effects/figures/tissue_and_marker_effects.png) · [Leave-one-out estimates](../step03_paired_effects/leave_one_out.csv)

The reference atlas shows prominent `TREM2` expression in maternal-origin macrophages and prominent `C1QA` and `C1QC` expression in fetal Hofbauer cells as well as other macrophages. The browser and pseudobulk reflect the same atlas, which has no isolated-FGR group. They provide cell context for the genes while leaving the cellular source of the bulk difference unresolved. In the primary cohort, just two pregnancies per group fell in the shared villous hypoxia-score range, and their decidual scores interleaved. This describes partial overlap on one RNA signature at delivery. Injury burden, trophoblast death, and clearance activity would each require separate measurements. Sequencing rounds are inferred, and labor and delivery mode were unavailable in the released bulk metadata. [Cell localization and hypoxia overlap](../step04_hypoxia_context/README.md) · [Admati and Skarbianskis et al.](https://doi.org/10.1016/j.medj.2023.07.005) · [Cohort audit](../step01_cohort_audit/README.md)

## What would distinguish the explanations

In a new EOPE-only and isolated-FGR cohort matched on delivery age and labor, collect paired decidua and villi and balance both groups across sequencing batches. Quantify maternal decidual macrophage density, subtype proportions, and C1Q–TREM2 RNA per identified cell. Measure apoptotic trophoblast uptake per macrophage alongside those measurements. A lower bulk score accompanied by a shift in cell numbers or subtypes and similar uptake per cell would favor a change in cell composition. Lower uptake per cell at comparable cell abundance would support a functional difference. A preprint reporting lower uptake in EOPE than in healthy pregnancies motivates testing uptake directly in this disease comparison. [Belchamber et al. preprint](https://doi.org/10.64898/2026.02.26.708195)

## Figures for the datathon

Lead with the tissue and marker effects to show the two tissue estimates and the uncertainty around the decidual marker comparison. Use the paired donor scores, cohort and donor sensitivity plot, and villous hypoxia overlap plot to show how individual pregnancies and the available score overlap shape the interpretation. [Tissue and marker effects](../step03_paired_effects/figures/tissue_and_marker_effects.png) · [Paired scores by pregnancy](../step03_paired_effects/figures/paired_donor_scores.png) · [Cohort and donor sensitivity](../step03_paired_effects/figures/paired_effect_sensitivity.png) · [Villous hypoxia overlap](../step04_hypoxia_context/figures/hypoxia_vs_decidual_score.png)

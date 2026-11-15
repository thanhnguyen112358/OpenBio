# Hypothesis Card

## Title

Higher villous BMP-associated RNA when early preeclampsia accompanies fetal growth restriction

## Clinical/Biological Problem

Early-onset preeclampsia (EOPE) and fetal growth restriction (FGR) can occur in the same pregnancy. We ask whether chorionic villi sampled at delivery show a higher bone morphogenetic protein (BMP) RNA score when EOPE accompanies FGR, and whether the difference reflects RNA expression within placental cells, the proportions of those cells, or both.

## Dataset(s) and Relevant Cells

The OpenBio release of Ackerman and the Awamleh cohort provide bulk villous RNA for the EOPE+FGR versus isolated-FGR comparison. The Piekos cohort provides a related severe-preeclampsia comparison. We use the Admati single-cell reference, which sampled preeclampsia and control pregnancies, to locate transcripts by cell type. BMPR2 and SMAD9 RNA were highest in endothelial cells and GREM2 RNA was highest in syncytiotrophoblast. [Data sources](../README.md#data-used) · [Cell reference figure](../figures/atlas.png)

## Observed Result

The 27-gene score was higher in EOPE+FGR than in isolated FGR by **0.273 log₂ score units** in OpenBio (7 vs. 5 pregnancies, approximate 95% confidence interval 0.100 to 0.446) and **0.145** in Awamleh (20 vs. 18, 0.037 to 0.253). Both differences stayed positive when we removed any one pregnancy. In Piekos, the unadjusted severe-preeclampsia difference was +0.174 (25 vs. 36, 0.066 to 0.283). The age-adjusted estimate within the shared delivery weeks was +0.110 (−0.069 to 0.289). [Article and figures](article.md)

In OpenBio, adjusting for endothelial or syncytiotrophoblast marker RNA, or both, yielded estimates of +0.294 to +0.355 for the 27-gene score. All three intervals included zero. Endothelial marker adjustment reduced the SMAD9 estimate from +0.668 to +0.424 (−0.265 to 1.112), while the BMPR2 estimate went from +0.610 to +0.843 (0.015 to 1.672). GREM2 averaged +0.838 with an interval that included zero. [Score estimates](../results/composition_effects.csv) · [Gene estimates](../results/composition_gene_effects.csv)

## Proposed Mechanism

We hypothesize that higher BMPR2 expression within villous endothelial cells contributes to the higher bulk BMP RNA score in EOPE+FGR. We predict that BMPR2 RNA per endothelial cell will be higher in EOPE+FGR and that a positive bulk score difference will remain after accounting for measured endothelial and syncytiotrophoblast proportions. The follow-up uses nuclear phosphorylated SMAD1/5/9 to assess canonical BMP signaling, leaving the direction of any difference open.

## Why the Result Is Biologically Plausible

BMPR2 encodes a BMP receptor, and its RNA is abundant in placental endothelial cells in the single-cell reference. Yi and colleagues described altered placental BMP2 and receptor expression in early-onset preeclampsia. Deng and colleagues found that BMP2 can act through SMAD2/3 in trophoblast models. This supports measuring both phosphorylation readouts in the follow-up. [Admati et al.](https://doi.org/10.1016/j.medj.2023.07.005) · [Yi et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC8712873/) · [Deng et al.](https://www.sciencedirect.com/science/article/pii/S2352396423002293)

## Alternative Explanations/Confounders

Differences in cell proportions, delivery week, labor, delivery mode or FGR severity could contribute to the bulk signal. Marker RNA can vary with both cell abundance and expression per cell. The full OpenBio model with a linear delivery-week term gave +0.286 (0.068 to 0.503), and the shared-week subset, with only two EOPE+FGR pregnancies, gave +0.264 (−0.752 to 1.281). Removing the activin-related genes INHA and INHBA reduced the Awamleh estimate to +0.085 (−0.022 to 0.192). The OpenBio paired-tissue analysis had a larger point estimate in decidua than in villi, with an interval for the tissue difference that included zero. [Age, gene-set, and tissue analyses](article.md)

## Proposed Validation Experiment or Follow-up Analysis

In an independent cohort, compare EOPE+FGR with isolated FGR at overlapping delivery weeks. Record labor, delivery mode and FGR severity, and sample comparable villous regions. For each pregnancy, measure the bulk 27-gene score and quantify endothelial and syncytiotrophoblast proportions on adjacent sections. Measure BMPR2 RNA within endothelial cells, examine SMAD9 and GREM2 as secondary targets, and assess nuclear phosphorylated SMAD1/5/9 and SMAD2/3 in the same cell classes. Compare the score difference before and after accounting for the measured cell proportions. A precise estimate near zero for the adjusted bulk difference would favor a cell-composition explanation. A precise estimate near zero for within-endothelial BMPR2 would weaken the gene-specific prediction.

## Potential Intervention Relevance

BMPR2 and the phosphorylation measurements are candidates for studying placental vascular responses in EOPE+FGR. Diagnostic relevance would depend on whether an accessible measurement taken earlier in pregnancy predicts the outcome. Therapeutic relevance would depend on whether altering the pathway improves outcomes.

## Confidence Level and Limitations

**Medium confidence in the bulk RNA association and low confidence in the proposed cellular mechanism.** Two cohorts agree on direction, with a similar unadjusted pattern in the Piekos severe-preeclampsia comparison. Limited delivery-age overlap, imprecise marker-adjusted estimates and the contribution of activin-related RNA limit interpretation. Bulk RNA from delivery samples cannot establish earlier changes, cell-specific effects or pathway activity.

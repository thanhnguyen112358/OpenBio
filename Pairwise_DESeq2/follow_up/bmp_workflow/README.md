# BMP-associated RNA in villous tissue at delivery

This exploratory follow-up estimates the additional RNA difference associated with early-onset preeclampsia (EOPE) among pregnancies with fetal growth restriction (FGR). Both cohorts informed the choice of question during the pilot. The [article](python_submission/docs/article.md) gives the written interpretation, the [hypothesis card](python_submission/docs/hypothesis-card.md) states the prediction, and the [submission README](submission/README.md) gives the runnable procedure. The [3:44 narrated draft](presentation-draft.mp4) uses the article’s headline visual. It is a local review copy with synthetic narration. The team recording and unlisted upload remain pending.

## What this hypothesis says

The question is narrow. When a baby is already growing too slowly, does the placenta look different if the mother also has early preeclampsia? The answer these data can support is yes, in one measured way. A set of genes used for BMP signaling has more RNA in those placentas. The marker check says this is not simply more of the cells that usually use those genes, and one gene, BMPR2, carries the increase more cleanly than the rest of the list.

A placenta passes oxygen and food from the mother to the baby. The two pieces measured here are the villi, the finger-like tissue where that exchange happens, and the decidua, the mother’s tissue where the placenta attaches. Fetal growth restriction means the baby grew less than expected. Early-onset preeclampsia means the mother developed high blood pressure and related illness before 34 weeks. The comparison is pregnancies with both against pregnancies with a small baby alone, so a difference is more likely to be about the mother’s disease than about the baby’s size.

Cells send protein messages. BMP is one family of those messages. A sending cell releases a BMP protein. A receiving cell catches it with a receptor on its surface, including BMPR2. Inside the cell, SMAD proteins then change which genes are turned on. SMAD9 is one of those. Other proteins, including GREM2, can block the message before it arrives. This analysis measures RNA, the working copy a cell makes when a gene is in use. More RNA usually means the gene is more active. It does not show that the BMP signal fired. Twenty-seven genes from the pathway are averaged into one number per placenta, the BMP score. A difference of 1 on this scale is about twice as much RNA. A difference of 0.27, the main result, is about 20% more.

The hypothesis, after the marker check, is that placentas from the combined group have a higher BMP RNA score, that the score stays in that direction after delivery week and after marker RNA from blood-vessel cells and the outer placental layer, and that BMPR2 fits a within-cell change better than SMAD9. SMAD9 travels with the blood-vessel markers. The prediction for a new, age-overlapping cohort is that counted cell numbers will leave the score higher, and that BMPR2 RNA will be higher inside blood-vessel cells. Nuclear phosphorylated SMAD1/5/9 inside those cells is the protein test. If the cell counts explain the RNA score, the within-cell claim fails.

Each placenta is one dot. The reported difference is the average of the combined-group dots minus the average of the growth-restriction-only dots. The 95% interval is the range the data can still support. If it stays above zero, the combined group is higher after allowing for chance. The OpenBio comparison is 7 versus 5, so two other public studies are used as checks. A second question uses the same dots. If the combined pregnancies simply contain more of the cells that normally carry BMPR2 and GREM2, marker genes for those cells should rise with the BMP score. A marker score is built the same way, and the model asks whether including it removes the BMP difference.

In OpenBio the combined group scored 0.273 higher (interval 0.100 to 0.446). Removing any one pregnancy kept the difference positive. In Awamleh it was 0.145 higher (20 versus 18; 0.037 to 0.253). In Piekos, severe preeclampsia with growth restriction versus growth restriction alone, it was 0.174 higher (25 versus 36; 0.066 to 0.283). After a delivery-week adjustment inside Piekos’s shared weeks, that interval includes zero.

Inside the 12 OpenBio pregnancies, delivery week and the BMP score barely move together (r = 0.12). A straight-line adjustment for week leaves the difference at 0.286 (0.068 to 0.503). Keeping only the shared weeks leaves two combined pregnancies, delivered at 27.0 and 27.7 weeks, against five isolated-FGR pregnancies. That estimate is 0.264, with an interval from −0.752 to 1.281.

Blood-vessel marker RNA is higher by 0.242 (0.032 to 0.452). The outer-layer marker score is 0.319 higher, and its interval includes zero (−0.279 to 0.917). The BMP score correlates with those marker scores at 0.28 and 0.10. After they are included, the BMP difference stays around 0.29 to 0.35. Those adjusted intervals include zero because 12 pregnancies cannot support that many terms at once. The estimate does not fall toward zero.

The genes inside the average do not move together. BMPR2 is higher by 0.610 (0.054 to 1.167) and stays higher after the blood-vessel markers are included (0.843; 0.015 to 1.672). SMAD9 is higher by 0.668 (0.148 to 1.189), correlates with the blood-vessel score at 0.64, and falls to 0.424 (−0.265 to 1.112) once that score is included. GREM2 averages 0.838 higher (−0.136 to 1.812) and correlates with the outer-layer score at 0.65. SMAD1 and SMAD5 are not higher. On the nine pregnancies with both tissues, the decidual difference is 0.415 (0.058 to 0.772) and the villous difference is 0.284 (−0.024 to 0.592). Donors S01, S02 and S03, all in the combined group, have no decidual library.

The groups were not born at the same time: the combined median is 33 weeks and the isolated-FGR median is 29. Removing INHA and INHBA, two genes the list shares with activin, reduces Awamleh to 0.085 and that interval includes zero. Marker RNA is not a count of cells. Labor, delivery mode and the severity of growth restriction are not in these models. The RNA was collected at delivery, after the disease had already happened. These files cannot show that the BMP signal turned on inside a cell.

## Pregnancy-level results

Effects are the difference in mean RNA score in log₂ units. The score is the mean of log₂(normalized count + 1) across 27 genes. Intervals are approximate 95% Welch intervals. Each cohort has its own normalization.

| Comparison | Numerator vs denominator | Difference (95% interval) | Leave-one-out range |
| --- | ---: | ---: | ---: |
| OpenBio, EOPE+FGR vs FGR | 7 vs 5 | +0.273 (+0.100 to +0.446) | +0.230 to +0.304 |
| OpenBio, shared delivery-age range | 2 vs 5 | +0.264 (−0.752 to +1.281) | +0.138 to +0.391 |
| Awamleh, EOPE+FGR vs FGR | 20 vs 18 | +0.145 (+0.037 to +0.253) | +0.125 to +0.165 |
| Awamleh, EOPE only vs FGR | 20 vs 18 | +0.168 (+0.060 to +0.276) | +0.155 to +0.188 |
| Awamleh, FGR vs preterm controls | 18 vs 21 | +0.146 (+0.047 to +0.244) | +0.126 to +0.165 |

Removing the two flagged isolated-FGR pregnancies individually gave +0.301 (S04 omitted) and +0.274 (S15 omitted). Omitting both gave +0.312 (+0.079 to +0.545). Combined-case donor S12 had the lowest HTSeq assignment fraction in villi (17%) and a clinical BP flag. Leaving S12 out gave +0.296 in the full comparison and +0.391 within shared age, where only one combined case remained. These checks preserve direction and leave the shared-age analysis imprecise. The combined group had a median delivery age of 33 weeks versus 29 for isolated FGR. Within 26.714–30.429 weeks, both combined cases (27.000 and 27.714) delivered earlier than four of five isolated-FGR cases. Both arms were in the same inferred sequencing round. [Donor scores](bmp_donor_scores.csv) · [Effects](bmp_effects.csv) · [Cohort manifest](cohort_manifest.csv) · [Donors in each score comparison](comparison_membership.csv) · [Samples in each count fit](rank_fit_membership.csv)

The Awamleh context comparisons are consistent with elevation in both EOPE strata and in FGR alone. The data leave a combined-phenotype interaction unresolved. All 79 released Awamleh samples were used. The authors used additional covariates and removed outliers, leaving 16 combined and 15 isolated-FGR pregnancies for their mRNA analysis. Their supplementary PCA points are unlabeled, so the retained GSM identifiers could not be verified. The GEO matrix exposed diagnoses without per-sample ages, fetal sex or delivery mode. [Awamleh et al. and supplement](https://pmc.ncbi.nlm.nih.gov/articles/PMC6598374/#Sup2) · [GEO GSE114691](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114691)

## Rank analysis and pathway specificity

| Branch | OpenBio NES (FDR) | OpenBio shared-age NES (FDR) | Awamleh NES (FDR) |
| --- | ---: | ---: | ---: |
| BMP | +1.96 (0.0018) | +1.61 (0.0380) | +1.79 (0.0077) |
| Activin | +1.81 (0.0072) | +0.94 (0.5363) | +1.45 (0.0757) |
| TGF-β receptor complex | +1.59 (0.0072) | +1.03 (0.5363) | +1.47 (0.0229) |

The rank results reproduce the pilot’s 13,725-gene universe and 25 BMP members. FDR adjusts the three branches within each fit. A positive enrichment statistic indicates that BMP genes tended to rank toward higher expression in the combined group. It summarizes the same pregnancies as the donor scores. The seven shared leading-edge genes were BMPR1B, BMPR2, GREM2, INHA, INHBA, SMAD4 and SMAD9. BMPR2’s unshrunken log₂ fold changes were +0.553, +0.390 and +0.288, respectively. The age-restricted Wald interval for BMPR2 was −0.216 to +0.996. BMPR2 met genome-wide FDR 0.05 in Awamleh (0.0105), with FDR 0.124 in full OpenBio. [Rank universe](rank_universe.csv) · [Branch results](branch_enrichment.csv) · [Gene tests and Wald intervals](bmp_gene_effects.tsv.gz)

Removing INHA and INHBA reduced the score difference to +0.171 (+0.002 to +0.339) in OpenBio, +0.198 (−1.243 to +1.639) in shared-age OpenBio and +0.085 (−0.022 to +0.192) in Awamleh. Awamleh isolated FGR versus preterm controls also changed from +0.146 to −0.013 after the removal. The effect sizes and related branch enrichments support considering a broader TGF-β superfamily response. The gene set measures RNA from receptors, mediators and antagonists with different roles. Protein activity remains unresolved. [Effects](bmp_effects.csv) · [Reactome BMP signaling](https://reactome.org/content/detail/R-HSA-201451)

## Candidate source cells

In the organizer’s Admati reference, the highest median RNA levels for BMPR2 and SMAD9 were 427.9 and 137.7 counts per million in endothelial cells, respectively (26 donors for each, all with at least 10 reads). GREM2 was highest in syncytiotrophoblast (285.4 counts per million, 26 donors, all with at least 10 reads). BMPR1B was sparse, with its highest median in extravillous trophoblast (11.3 counts per million, 15 donors, 5 with at least 10 reads). BMPR1B’s official browser violin and t-SNE assets returned 404. BMPR2, GREM2 and SMAD9 displayed successfully. The four classes shown as absent use exact organizer marker-class names that are absent from the released pseudobulk taxonomy. Related broader atlas classes are shown by their released names. [Donor-class summary](atlas_summary.csv) · [Browser notes](browser_lookups.csv) · [Admati et al.](https://doi.org/10.1016/j.medj.2023.07.005)

PECAM1, VWF, CDH5 and CLDN5 had positive fold changes in full OpenBio, shared-age OpenBio and Awamleh. The mapped syncytiotrophoblast markers CGA, ERVW-1 and PSG1 also pointed up in all three. The follow-up below asks whether those marker RNAs account for the 27-gene score. [Bulk marker effects](bulk_marker_effects.csv) · [Hypothesis card](python_submission/docs/hypothesis-card.md)

## Marker-RNA adjustment

On the OpenBio 7 versus 5 comparison, the same size factors and the same log₂(normalized count + 1) mean give an endothelial score (PECAM1, VWF, CDH5, CLDN5) of +0.242 (+0.032 to +0.452) and a syncytiotrophoblast score (CGA, ERVW-1, PSG1) of +0.319 (−0.279 to +0.917). The 27-gene score correlates with those scores at r = 0.28 and r = 0.10, and with delivery week at r = 0.12. HC3 group coefficients for the 27-gene score:

| Adjustment | Difference (HC3 interval) |
| --- | ---: |
| None | +0.273 (+0.084 to +0.463) |
| Endothelial markers | +0.332 (−0.011 to +0.674) |
| Syncytiotrophoblast markers | +0.294 (−0.008 to +0.597) |
| Both marker scores | +0.355 (−0.087 to +0.796) |
| Delivery week, linear | +0.286 (+0.068 to +0.503) |
| Both markers and week | +0.374 (−0.203 to +0.952) |

The point estimates stay near the unadjusted difference. The marker-adjusted intervals include zero because the estimates are less precise. On the same scale, BMPR2 is +0.610 (+0.054 to +1.167) and remains +0.843 (+0.015 to +1.672) after the endothelial score. SMAD9 is +0.668 (+0.148 to +1.189), correlates with the endothelial score at r = 0.64, and falls to +0.424 (−0.265 to +1.112) after that adjustment. GREM2 averages +0.838 (−0.136 to +1.812) and correlates with the syncytiotrophoblast score at r = 0.65. BMPR1A is −0.142 (−0.410 to +0.127). SMAD1 and SMAD5 are not higher. Marker RNA can rise because a cell type is more abundant or because each cell makes more of the marker. [Donor scores](python_submission/results/composition_donor_scores.csv) · [Effects](python_submission/results/composition_effects.csv) · [Figure](python_submission/figures/composition.png)

## Severe preeclampsia context

In Piekos, the unadjusted 27-gene estimate was +0.174 (+0.066 to +0.283) for the released labels `FGR+PE Severe` versus `FGR` (25 vs 36). One released RNA ID, MJ-0323, had a unique clinical counterpart with the terminal suffix -1. With exact-ID links only, the age-adjusted estimate within 32.71–38.29 weeks was +0.141 (−0.021 to +0.303, 15 vs 16). Including the one suffix match gave +0.110 (−0.069 to +0.289, 15 vs 17). The earlier 31-gene definition with a +0.5 pseudocount is exported separately. Only one isolated-FGR pregnancy delivered before 34 weeks. Onset age was unavailable. [Piekos effects](piekos_effects.csv) · [Piekos et al.](https://www.nature.com/articles/s42003-025-08631-6)

## Paired villi and decidua

Four combined and five isolated-FGR pregnancies contributed both tissues. All 27 score genes passed the same count threshold in decidua. The group difference was +0.284 (−0.024 to +0.592) in villi and +0.415 (+0.058 to +0.772) in decidua. The estimate of how much the group difference changes between tissues was +0.131 (−0.208 to +0.471). The paired interval leaves the difference between tissues uncertain. [Paired tissue effects](paired_tissue_effects.csv) · [Donor scores](paired_donor_scores.csv) · [Decidua coverage](paired_gene_coverage.csv)

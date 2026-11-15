# BMP-associated RNA in villous tissue at delivery

This exploratory follow-up estimates the additional RNA difference associated with early-onset preeclampsia (EOPE) among pregnancies with fetal growth restriction (FGR). Both cohorts informed the choice of question during the pilot. The [article](submission/docs/article.md) gives the general-audience interpretation and the [submission README](submission/README.md) gives the runnable procedure. The [3:44 narrated draft](presentation-draft.mp4) uses the article’s headline visual. It is a local review copy with synthetic narration. The team recording and unlisted upload remain pending.

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

PECAM1, VWF, CDH5 and CLDN5 had positive fold changes in full OpenBio, shared-age OpenBio and Awamleh. The mapped syncytiotrophoblast markers CGA, ERVW-1 and PSG1 also pointed up in all three. The markers are consistent with composition contributing to the bulk shift, with within-cell regulation also plausible. At delivery the atlas provides RNA localization, and it lacks an isolated-FGR arm. The proposed independent follow-up quantifies cell proportions and nuclear pSMAD1/5/9 by cell type in the two FGR groups, with pSMAD2/3 as a complementary readout. [Bulk marker effects](bulk_marker_effects.csv) · [Hypothesis card](submission/docs/hypothesis-card.md)

## Severe preeclampsia context

In Piekos, the unadjusted 27-gene estimate was +0.174 (+0.066 to +0.283) for the released labels `FGR+PE Severe` versus `FGR` (25 vs 36). One released RNA ID, MJ-0323, had a unique clinical counterpart with the terminal suffix -1. With exact-ID links only, the age-adjusted estimate within 32.71–38.29 weeks was +0.141 (−0.021 to +0.303, 15 vs 16). Including the one suffix match gave +0.110 (−0.069 to +0.289, 15 vs 17). The earlier 31-gene definition with a +0.5 pseudocount is exported separately. Only one isolated-FGR pregnancy delivered before 34 weeks. Onset age was unavailable. [Piekos effects](piekos_effects.csv) · [Piekos et al.](https://www.nature.com/articles/s42003-025-08631-6)

## Paired villi and decidua

Four combined and five isolated-FGR pregnancies contributed both tissues. All 27 score genes passed the same count threshold in decidua. The group difference was +0.284 (−0.024 to +0.592) in villi and +0.415 (+0.058 to +0.772) in decidua. The estimate of how much the group difference changes between tissues was +0.131 (−0.208 to +0.471). The paired interval leaves the difference between tissues uncertain. [Paired tissue effects](paired_tissue_effects.csv) · [Donor scores](paired_donor_scores.csv) · [Decidua coverage](paired_gene_coverage.csv)

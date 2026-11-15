# H2 evidence for the OpenBio video

Harvard Undergraduate OpenBio Datathon 2026, Track 1. This sheet uses the archived Python results. The team estimated the contrasts below from the released datasets. [H2 article](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/docs/article.md) · [Analysis plan](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/docs/analysis-plan.md) · [Python results](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/README.md) · [Code and data](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/README.md)

## 1. Clinical question

Among pregnancies with fetal growth restriction (FGR), does early-onset preeclampsia (EOPE) accompany a higher BMP-associated RNA score in chorionic villi at delivery? EOPE refers to preeclampsia with onset before 34 weeks. BMP stands for bone morphogenetic protein. The primary comparison uses the released clinical labels and subtracts the isolated-FGR mean from the EOPE+FGR mean. Ackerman and colleagues previously studied these clinical groups. The H2 reanalysis estimates the size and consistency of the RNA association. [H2 article](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/docs/article.md) · [Ackerman et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10192030/)

## 2. Delivery age

In OpenBio, median delivery was at 33 weeks in EOPE+FGR and 29 weeks in isolated FGR. The shaded shared range is 26.714–30.429 weeks. It includes 2 of 7 combined pregnancies and all 5 isolated-FGR pregnancies. The two combined deliveries were at 27.000 and 27.714 weeks, earlier than four of the five isolated-FGR deliveries. The primary groups share the same inferred sequencing round. In the age plot each point represents one pregnancy and triangles indicate blood-pressure flags. [Age figure](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/figures/age_support.png) · [Age support table](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/cohort_support.csv) · [Cohort audit](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/cohort_manifest.csv)

Using all 12 OpenBio pregnancies and a straight-line term for delivery week gives +0.286 log₂ score units, with an approximate 95% HC3 interval of +0.068 to +0.503. Restriction to the shared range gives +0.264, with an approximate 95% Welch interval of −0.752 to +1.281. The linear model uses observations outside the shared range and assumes a common linear age relationship. The shared-range estimate is imprecise and retains age imbalance. The public Awamleh release has diagnosis labels, and individual delivery weeks are unavailable. [OpenBio sensitivity table](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/results/composition_effects.csv) · [Pregnancy score table](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/bmp_effects.csv) · [Awamleh et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6598374/)

## 3. Headline pregnancy-level score

Each pregnancy contributes one BMP27 score, the mean of log₂(normalized count + 1) across 27 genes. Counts use median-ratio size factors within each cohort. The gene set starts with 31 genes from MSigDB Reactome 2024.1.Hs. Twenty-seven have a unique mapping and at least 10 raw reads in five or more libraries in each full villous cohort. The score differences use log₂ score units and have approximate two-sided 95% Welch intervals. Compare groups within cohorts because normalization gives each cohort its own score scale. [Gene manifest](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/bmp_gene_manifest.csv) · [Analysis plan](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/docs/analysis-plan.md) · [Reactome BMP signaling](https://reactome.org/content/detail/R-HSA-201451)

| Cohort | EOPE+FGR vs isolated FGR | Difference | Approximate 95% Welch CI | Leave-one-out estimates |
| --- | --- | ---: | --- | --- |
| OpenBio | 7 vs 5 | +0.273 | +0.100 to +0.446 | +0.230 to +0.304 |
| OpenBio, shared delivery ages | 2 vs 5 | +0.264 | −0.752 to +1.281 | +0.138 to +0.391 |
| Awamleh | 20 vs 18 | +0.145 | +0.037 to +0.253 | +0.125 to +0.165 |

In the headline figure, each dot in the upper panels represents a pregnancy and the bars mark group means. In the lower panel the points mark differences, the horizontal lines show intervals, and the vertical line marks zero. Use the same headline figure as the article. [Headline figure](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/figures/headline.png) · [Pregnancy score effects](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/bmp_effects.csv) · [Leave-one-out estimates](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/bmp_leave_one_out.csv)

## 4. Enrichment and gene-set sensitivity

The Python count fits rank genes by signed Wald statistic. OpenBio full and shared-age fits include fetal sex and clinical group. Awamleh includes clinical group. GSEApy tests three Reactome branches per fit and adjusts their p-values using Benjamini–Hochberg. The Python rank universe has 13,726 genes, including 25 Reactome BMP genes. NES is the normalized enrichment score. A positive NES places BMP genes toward the higher-expression end for EOPE+FGR. The rank analysis reuses the pregnancies in the score comparison. [Enrichment figure](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/figures/enrichment.png) · [Branch enrichment](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/branch_enrichment.csv) · [Analysis plan](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/docs/analysis-plan.md)

| Python fit | BMP NES | BMP FDR |
| --- | ---: | ---: |
| OpenBio full | +1.97 | 0.000822 |
| OpenBio shared delivery ages | +1.57 | 0.0351 |
| Awamleh | +1.80 | 0.00537 |

Seven BMP genes contribute to the shared leading edge across all three fits. They are BMPR1B, BMPR2, GREM2, INHA, INHBA, SMAD4 and SMAD9. A leading-edge gene helps drive the enrichment statistic. Among these seven genes, BMPR2 and SMAD4 pass a 5% genome-wide gene-level FDR threshold in Awamleh (0.0105 for each when rounded to four decimals). The pathway FDR in the table covers the three branch tests within a fit. [Shared leading edge](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/shared_leading_edge.csv) · [Gene-level count effects](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/bmp_gene_effects.tsv.gz)

Removing INHA and INHBA, two genes shared with activin signaling, reduces the OpenBio estimate to +0.171 (+0.002 to +0.339) and the Awamleh estimate to +0.085 (−0.022 to +0.192). These are approximate 95% Welch intervals for the 25-gene score sensitivity. This 25-gene sensitivity starts from BMP27 and removes two genes. The 25 members in the rank analysis come from the separate genome-wide rank filter. Other TGF-β-family branches also shifted in the full cohorts, consistent with a broader family-level RNA response. [Score sensitivity](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/bmp_effects.csv) · [Branch enrichment](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/branch_enrichment.csv)

## 5. Cell atlas localization

In the Admati donor-by-cell-class reference, BMPR2 and SMAD9 had their highest median expression in endothelial cells (427.9 and 137.7 counts per million, respectively). GREM2 had its highest median in syncytiotrophoblast (285.4 counts per million). Each of these summaries included 26 eligible donors with at least 10 reads for the gene. BMPR1B was sparse, with its highest median in extravillous trophoblast (11.3 counts per million, 5 of 15 eligible donors with at least 10 reads). [Atlas source table](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/atlas_summary.csv) · [Admati et al.](https://doi.org/10.1016/j.medj.2023.07.005)

Darker dots indicate higher median counts per million. Larger dots indicate a higher fraction of donors with at least 10 reads for that gene. Parentheses show eligible donors per cell class, and crosses mark unavailable summaries. Every contributing donor-class has at least 20 cells. The reference includes preeclampsia and controls and has no isolated-FGR arm. It identifies candidate cell sources, and within-cell expression differences in the H2 contrast remain open. [Atlas figure](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/figures/atlas.png) · [Atlas source table](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/atlas_summary.csv)

The atlas browser inspection also found strongest visible BMPR2 and SMAD9 expression in vascular populations, including arterial, capillary and venous endothelium. GREM2 was concentrated in syncytiotrophoblast, labeled SCT in the browser. The BMPR1B browser panel was unavailable, and its localization comes from the released donor-class counts. [Browser lookups](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/browser_lookups.csv) · [Official atlas browser](https://zeisellab.org/preeclampsia/)

## 6. Marker-RNA sensitivity

In the same OpenBio 7 vs 5 pregnancy comparison, the endothelium marker panel is PECAM1, VWF, CDH5 and CLDN5. The syncytiotrophoblast panel is CGA, ERVW-1 and PSG1. They use the same size factors and log₂ mean as BMP27. The endothelial score difference is +0.242 (+0.032 to +0.452) and the syncytiotrophoblast score difference is +0.319 (−0.279 to +0.917), with approximate 95% Welch intervals. The two marker scores correlate with BMP27 at Pearson r = 0.28 and r = 0.10, respectively. [Marker effect table](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/results/composition_effects.csv) · [Marker donor scores](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/results/composition_donor_scores.csv)

Each row in the marker figure gives the group coefficient for the same 27-gene score from an ordinary least squares model. All rows in this figure, including the unadjusted row, use approximate two-sided 95% HC3 intervals. HC3 is a heteroskedasticity-robust variance estimate. All marker-adjusted point estimates stay positive and their intervals cross zero. Marker RNA can reflect cell abundance, marker expression per cell, or both. True cell fractions remain unmeasured. [Marker figure](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/figures/composition.png) · [Model outputs](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/results/composition_effects.csv) · [Model code](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/notebooks/composition.py)

| Adjustment | Difference | Approximate 95% HC3 CI |
| --- | ---: | --- |
| Unadjusted | +0.273 | +0.083 to +0.463 |
| Endothelial markers | +0.332 | −0.011 to +0.674 |
| Syncytiotrophoblast markers | +0.294 | −0.008 to +0.597 |
| Both marker scores | +0.355 | −0.087 to +0.796 |
| Delivery week, linear | +0.286 | +0.068 to +0.503 |
| Both marker scores and delivery week | +0.374 | −0.203 to +0.951 |
| Both marker scores and perivascular markers | +0.331 | −0.234 to +0.895 |

In exploratory individual-gene comparisons, bulk BMPR2 is +0.610 (+0.054 to +1.167, Welch) and +0.843 (+0.015 to +1.672, HC3) after endothelial marker adjustment. SMAD9 is +0.668 (+0.148 to +1.189, Welch) and +0.424 (−0.265 to +1.112, HC3) after the same adjustment. BMPR2 and SMAD9 correlate with the endothelial score at r = 0.15 and r = 0.64, respectively. GREM2 correlates with the syncytiotrophoblast score at r = 0.65 and has intervals crossing zero before and after adjustment. These gene effects are log₂ normalized-expression differences with nominal intervals. [Marker-adjusted gene effects](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/results/composition_gene_effects.csv)

## 7. Falsifiable follow-up

Recruit an independent set of EOPE+FGR and isolated-FGR pregnancies with overlapping delivery ages, use matched villous sampling regions, and record labor, delivery mode and FGR severity. The directional predictions are a higher BMP27 bulk score after adjustment for counted endothelial and syncytiotrophoblast fractions and higher BMPR2 RNA within villous endothelial cells. Treat SMAD9 and GREM2 as secondary RNA targets. If measured cell fractions account for the bulk association, that weakens the within-cell explanation for the score. A within-endothelial BMPR2 comparison tests the gene-specific prediction directly. [H2 hypothesis card](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/docs/hypothesis-card.md)

Nuclear phosphorylated SMAD1/5/9 within the relevant cell classes is the proposed signaling readout. The predicted direction of any phosphorylation difference is open. The hypothesis card also includes pSMAD2/3 as a complementary readout because activin-related transcripts contribute to BMP27. Summarize measurements per pregnancy. [H2 hypothesis card](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/docs/hypothesis-card.md) · [Yi et al.](https://doi.org/10.3389/fcell.2021.768669) · [Deng et al.](https://doi.org/10.1016/j.ebiom.2023.104664)

## 8. Credits and data sources

- OpenBio Track 1 supplied the curated Ackerman bulk counts, clinical metadata and Admati donor-class reference. [Pinned OpenBio data release](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/tree/fb751e91fb1eca4fb26a9658ad9d5042c382c4ba/track-01/data) · [Ackerman et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10192030/) · [Ackerman GEO GSE203507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203507) · [Admati et al.](https://doi.org/10.1016/j.medj.2023.07.005)
- Awamleh and colleagues released transcript counts and diagnoses. The team aggregated transcript IDs using Ensembl GRCh37 release 75. [Awamleh et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6598374/) · [Awamleh GEO GSE114691](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114691) · [Ensembl release 75](https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/)
- Reactome supplies the pathway definitions, using the MSigDB Reactome 2024.1.Hs gene sets. [Reactome BMP signaling](https://reactome.org/content/detail/R-HSA-201451) · [MSigDB release](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/)
- The optional severe-PE context uses the released Piekos RNA and clinical data. [Piekos et al.](https://www.nature.com/articles/s42003-025-08631-6) · [Piekos data and code](https://github.com/spiekos/Hood-Lab-Placental-Multiomics)
- Thank the study participants, contributing study teams and OpenBio organizers. Add the verified team byline to the final video. [H2 data and code](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/README.md) · [Organizer submission instructions](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/fb751e91fb1eca4fb26a9658ad9d5042c382c4ba/submissions/README.md)

## Evidence appendix

### Optional OpenBio diagnostic figures

The PCA uses the 500 most variable genes after variance stabilization. Each point represents one villous pregnancy, shown in three panels colored by diagnosis, fetal sex and inferred sequencing round. The volcano plot uses the full sex-adjusted OpenBio PyDESeq2 contrast, with estimated gene-level log₂ fold change on the horizontal axis and −log₁₀ of the nominal Wald p-value on the vertical axis. The exported table contains the genome-wide adjusted values. Forty-five genes passed a genome-wide FDR threshold of 0.05 in that fit, with zero Reactome BMP members passing. [PCA figure](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/figures/pca.png) · [PCA scores](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/pca_scores.csv) · [Volcano figure](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/figures/volcano.png) · [Gene-level count effects](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/bmp_gene_effects.tsv.gz)

### Awamleh direction check and phenotype specificity

Both OpenBio and Awamleh informed selection of the comparison in the pilot. Awamleh provides a direction check in a separate cohort within this exploratory follow-up. Prospective confirmation would require an independent cohort unseen during selection. The group-only Awamleh count fit uses all 79 public libraries. The original authors applied QC filtering, and the public records did not identify the retained sample IDs needed to reconstruct that subset. [Analysis plan](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/docs/analysis-plan.md) · [H2 workflow notes](../Pairwise_DESeq2/follow_up/bmp_workflow/README.md) · [Awamleh et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6598374/)

In Awamleh, EOPE alone vs isolated FGR was +0.168 (+0.060 to +0.276, 20 vs 18), and isolated FGR vs preterm controls was +0.146 (+0.047 to +0.244, 18 vs 21). These are BMP27 mean differences with approximate 95% Welch intervals. Elevation occurs across additional clinical comparisons, and an interaction specific to EOPE+FGR remains unresolved. [Awamleh score effects](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/bmp_effects.csv)

### RNA expression and signaling activity

The 27-gene average pools transcripts from receptors, mediators, antagonists and components shared with activin signaling. OpenBio BMPR2 and SMAD9 point up, GREM2 is positive and imprecise, and SMAD1 and SMAD5 have negative point estimates with intervals crossing zero. RNA measured at delivery estimates an association at that time. Timing of the change, protein abundance, phosphorylation and causal direction remain open. Nuclear pSMAD1/5/9 is the proposed assay for canonical BMP signaling activity. [Gene effects](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/results/composition_gene_effects.csv) · [H2 article](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/docs/article.md) · [Reactome BMP signaling](https://reactome.org/content/detail/R-HSA-201451)

### Piekos severe-PE context

The released `FGR+PE Severe` versus `FGR` labels give an unadjusted BMP27 difference of +0.174 (+0.066 to +0.283, 25 vs 36). In the shared delivery-age range, the age-adjusted estimate is +0.110 (−0.069 to +0.289, 15 vs 17). Both intervals are approximate 95% HC3. Only one isolated-FGR pregnancy delivered before 34 weeks and onset age was unavailable. This comparison supplies supporting severe-PE context. [Piekos score effects](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/piekos_effects.csv) · [H2 workflow notes](../Pairwise_DESeq2/follow_up/bmp_workflow/README.md) · [Piekos et al.](https://www.nature.com/articles/s42003-025-08631-6)

### Paired tissue context

Nine OpenBio pregnancies contributed both villi and decidua (4 combined and 5 isolated FGR). Three combined pregnancies lacked a decidual library. The analysis normalized tissues separately and compared the same 27-gene score in complete pairs. The group difference was +0.284 (−0.024 to +0.592) in villi and +0.415 (+0.058 to +0.772) in decidua. Comparing the per-pregnancy decidua-minus-villi scores between groups gives a difference of +0.131 (−0.208 to +0.471). All intervals are approximate 95% Welch. The decidual point estimate is larger in these paired pregnancies, and the tissue-difference interval leaves tissue specificity unresolved. [Paired tissue effects](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/paired_tissue_effects.csv) · [Paired donor scores](../Pairwise_DESeq2/follow_up/bmp_workflow/python_results/paired_donor_scores.csv)

# Leads from the clinical contrast grid

The [Track 1 gene browser notes](browser_lookups.csv) record TMEM45A, NTRK2, and GPIHBP1. These three show different relationships between early-onset preeclampsia (EOPE) and fetal growth restriction (FGR) in bulk. The browser provides qualitative cell-type context and has no isolated FGR arm.

**TMEM45A is the strongest preeclampsia versus growth-restriction lead in this comparison panel.** The villous log₂ fold changes for EOPE without FGR versus isolated FGR were +1.65 in OpenBio and +0.34 in Awamleh. For EOPE+FGR versus isolated FGR they were +1.67 and +0.35. The OpenBio combined comparison has contrast-wide FDR 0.059; the Awamleh FDR is 0.0053. For isolated FGR versus preterm controls, the estimates were +0.84 in OpenBio (95% Wald CI −0.10 to +1.79) and +0.07 in Awamleh (−0.10 to +0.24). For EOPE without FGR versus iPTB, OpenBio estimates range from +2.45 to +2.62 over the three restrictions; the strictest fit has five pregnancies.

This is a known PE-associated gene. [Horii et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10081686/) found an increase in PE even when both groups had maternal vascular malperfusion. [Sánchez-Corrales et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC13398505/) report TMEM45A expression changes in a hypoxia-associated syncytiotrophoblast state in severe PE. In the organizer's Admati reference, among donor profiles with at least 20 cells per class, TMEM45A had the highest median CPM in stroma (90.2) and extravillous trophoblast (64.2); it was also detected across other classes. This is an expression map and does not localize the PE–FGR bulk difference. The most informative test is a donor-level comparison of that trophoblast state and stromal expression in appropriately comparable EOPE and isolated FGR tissue. [Admati et al.](https://doi.org/10.1016/j.medj.2023.07.005)

**GPIHBP1 is a shared placental disease candidate.** In the local villous data its ratio to the geometric mean of PECAM1, KDR, CDH5, VWF, and EGFL7 was +1.87 log₂ units (unadjusted Welch CI +1.04 to +2.70; 7 vs 5) in EOPE+FGR versus FGR. The per-donor ratio uses log₂(raw count + 0.5) for the target and marker genes. This exploratory check assumes stable marker expression per endothelial cell. In Awamleh, GPIHBP1 was up in all three disease groups versus the preterm control group (+1.78 to +2.17 log₂); direct comparisons to isolated FGR were smaller and their Wald intervals included zero. The local villous Hallmark cholesterol homeostasis set did not show clear enrichment. Neither the ratio nor the bulk data measures lipid transport or separates endothelial subtype composition from a within-cell change.

**Positive control:** the villous hypoxia and glycolysis Hallmarks were positively enriched across all four local EOPE versus iPTB fits and in Awamleh. This is consistent with [Ackerman et al.](https://doi.org/10.1161/HYPERTENSIONAHA.122.20807), who also used GSE114691 in their original publication. The local sensitivity fits reuse individuals, and Awamleh supplies a separate villous cohort for the current reanalysis. No cross-tissue interaction was tested.

`candidate_effects.csv` contains unshrunk estimates, unadjusted 95% Wald intervals, and FDR within each contrast. `sensitivity/` contains all gene estimates and cohort support for the three restricted local fits: common observed delivery-age range (5 vs 5), shared inferred sequencing round (4 vs 5), and their intersection recalculated inside the shared round (3 vs 2). The models use `~ fetal_sex + group`; the age restriction does not match donors or adjust residual differences in age, labor, or delivery mode. The restricted fits omit the FGR groups. Several comparisons reuse the same pregnancies. Genes in the panel were selected after looking at the data.

The [Awamleh cohort](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114691) has 20 EOPE-only, 20 EOPE+FGR, 18 isolated FGR, and 21 preterm labor control samples. `prepare_external.py` joins the submitted libraries to case labels, sums submitted transcript-feature counts using Ensembl 75, keeps genes with at least 10 counts in five or more libraries, and fits `~ group` with PyDESeq2. The public sample metadata used here do not support adjusting for clinical or sequencing covariates. This reanalysis does not recreate the [original Awamleh analysis](https://link.springer.com/article/10.1186/s12920-019-0548-x). The candidate panel links historical GPIHBP1 `ENSG00000182851` to the current organizer `ENSG00000277494` by the source gene annotation. PECAM1 has no row under its current ID in the resulting external panel.

`hallmark_eope_vs_preterm.csv` uses `fgseaMultilevel` (min 15, max 500, `eps=1e-10`) with signed Wald statistics for genes with a finite test and base mean ≥10. The same gene universe is used across sensitivity fits within a tissue, with a separately mapped universe for Awamleh. These are competitive gene set tests. They do not measure pathway activity or add independent pregnancies. The [MSigDB v2024.1.Hs Hallmark GMT](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt) has SHA-256 `ee2463540042078bfa3f67828e1e223bb354446d9fbb4d22845866835ba5c772`.

Run from the repository root with R DESeq2, fgsea, BiocParallel, and Python with pandas, NumPy, and PyDESeq2. After running the parent `run_pairwise.R`:

```powershell
$data = 'openbio-datathon-2026/track-01/data'
$audit = 'Deci_Vili_H1/results/step01_cohort_audit/donor_audit.csv'
$sens = 'Pairwise_DESeq2/follow_up/sensitivity'
$ext = Join-Path ([IO.Path]::GetTempPath()) 'openbio-pairwise-awamleh-gse114691'
Rscript Pairwise_DESeq2/check_eope.R $data $audit vt $sens
Rscript Pairwise_DESeq2/check_eope.R $data $audit db $sens
python Pairwise_DESeq2/prepare_external.py --output $ext
python Pairwise_DESeq2/summarize_follow_up.py Pairwise_DESeq2 $ext
Rscript Pairwise_DESeq2/screen_pathways.R Pairwise_DESeq2 "$data/track01_gene_annotation.csv" h.all.v2024.1.Hs.symbols.gmt "$ext/gse114691_pe_vs_control.tsv.gz"
```

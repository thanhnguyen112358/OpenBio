# Placental BMP expression in EOPE and FGR implementation plan

**Run status:** [Verified analysis and interpretations](../../../Pairwise_DESeq2/follow_up/bmp_workflow/README.md) and [staged submission package](../../../Pairwise_DESeq2/follow_up/bmp_workflow/submission/README.md). The [local 3:44 draft](../../../Pairwise_DESeq2/follow_up/bmp_workflow/presentation-draft.mp4) uses synthetic narration. The team byline, team recording and unlisted upload remain pending. The BMPR1B browser panel was unavailable, and the released Admati pseudobulk supplied the localization. The Awamleh QC-filtered donor IDs could not be verified.

> **For agentic workers:** Use `superpowers:subagent-driven-development` or `superpowers:executing-plans` when implementing this workflow. The checkboxes track analysis and submission deliverables.

**Goal:** Estimate how much villous BMP pathway gene expression differs between pregnancies with early-onset preeclampsia and fetal growth restriction (EOPE+FGR) and pregnancies with isolated FGR.

**Architecture:** Start with the clinical comparison, fit DESeq2 on raw counts, and test the curated BMP gene set. Use one score per pregnancy to show effect size and variability. Check the same contrast in Awamleh before using individual genes and the cell atlas to interpret the pattern.

**Tools:** R 4.6.1 with Bioconductor DESeq2 1.52.0, fgsea and BiocParallel. The existing Awamleh data preparation uses Python and PyDESeq2. Record the precise package versions and source checksums in the submission run.

**Basis:** [Exploratory results](../../../Pairwise_DESeq2/follow_up/clinical_contrasts/README.md) · [Track 1 workflow](../../../openbio-datathon-2026/track-01/3-workflow.md) · [Submission requirements](../../../openbio-datathon-2026/track-01/4-submit.md)

## What we are locking in

The main tissue is **chorionic villi**. The primary contrast is **EOPE+FGR minus isolated FGR**, with the original clinical labels kept intact. The primary estimate is the difference in mean BMP RNA score across pregnancies. A positive difference means the score is higher in EOPE+FGR. The first hypothesis is that the score rises further when EOPE co-occurs with FGR. The mechanism to test later is whether that RNA pattern accompanies altered BMP signaling within specific villous cells. The pilot informed this choice, including the Awamleh results, so this is a defined exploratory follow-up. An external cohort unseen at the time the question was chosen would be needed for a prospective confirmation.

| Cohort | Comparison and role | Interpretation limit |
| --- | --- | --- |
| OpenBio, primary | 7 EOPE+FGR vs 5 isolated FGR villous biopsies. The prior audit places both groups in the same inferred sequencing round. | Delivery age varies by group. Only 2 vs 5 fall within the shared observed range. |
| Awamleh, independent cohort checked during selection | 20 EOPE+FGR vs 18 isolated FGR in the released counts. Also use 20 EOPE-only and 21 preterm controls for context. | The public rerun uses group labels. The original paper used surrogate variables and removed outliers, leaving 16 vs 15 in its mRNA analysis. |
| Piekos, optional context | 25 severe PE+FGR vs 36 isolated FGR using the exact labels. The earlier age overlap check retained 15 vs 16. | Only one isolated-FGR pregnancy delivered before 34 weeks. PE onset age is unavailable. |
| Admati, interpretation | Placental browser and organizer donor-by-cell-class reference. | The reference has no isolated FGR arm. It can locate expression, while a cell-specific diagnosis comparison needs a cohort with both arms. |

Restrict claims to an association at delivery. FGR severity, gestational age, labor, delivery mode, placental sampling and cell composition can influence expression. Keep RNA and functional activity distinct. The BMP gene set includes receptors, mediators, antagonists and components shared with activin signaling.

## Workflow files and interfaces

Keep the existing DESeq2 and Awamleh scripts as the inputs to the final analysis. The filenames below describe **planned** work. The existing pilot tables are available now. Stage new results under `Pairwise_DESeq2/follow_up/bmp_workflow/` and copy only the materials required for submission into the actual team folder.

| File | Responsibility |
| --- | --- |
| [`Pairwise_DESeq2/run_pairwise.R`](../../../Pairwise_DESeq2/run_pairwise.R) | Existing OpenBio count fits and genome-wide contrasts. |
| [`Pairwise_DESeq2/prepare_external.py`](../../../Pairwise_DESeq2/prepare_external.py) | Existing Awamleh data preparation and group-only contrasts. Use `--gene-counts` for pathway scoring. |
| `Pairwise_DESeq2/bmp/analyze.R` | Planned command line analysis for cohort audit, fixed scores, DESeq2 sensitivity fits, enrichment and numerical exports. |
| `Pairwise_DESeq2/bmp/report.R` | Planned figure generation from the exported data only. |
| `Pairwise_DESeq2/follow_up/bmp_workflow/` | Planned results directory for manifests, score and gene results, browser notes, package versions and figures. |

Prepare the existing inputs from the workspace root. Choose a cache directory inside the project and keep large raw downloads out of the submission folder.

```powershell
$data = 'openbio-datathon-2026/track-01/data'
$cache = 'Pairwise_DESeq2/cache/bmp'
Rscript Pairwise_DESeq2/run_pairwise.R $data "$cache/openbio"
python Pairwise_DESeq2/prepare_external.py --annotation "$data/track01_gene_annotation.csv" --output "$cache/awamleh" --gene-counts
Invoke-WebRequest 'https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c2.cp.reactome.v2024.1.Hs.symbols.gmt' -OutFile "$cache/reactome.gmt"
```

The planned analysis interface takes `<track01-data-dir> <donor-audit-csv> <openbio-de-dir> <awamleh-dir> <reactome-gmt> <output-dir>` as positional paths. The planned figure interface takes `<output-dir>` and writes `figures/` beneath it. These two entry points need implementation. The scripts above already exist.

## 1. Verify the clinical comparison and freeze the inputs

- [x] Join raw count columns to metadata by `sample_id`. Join to the [cohort audit](../../../Deci_Vili_H1/results/step01_cohort_audit/donor_audit.csv) by `donor_id`. Check exact sample ID agreement, nonnegative integer counts, unique sample and donor IDs within each tissue, and agreement of diagnosis across paired tissues. Record source checksums, sample counts and normalization settings in `run_manifest.tsv`.
- [x] Use `ga_delivery_weeks`, the organizer's converted week-plus-day field. Record fetal sex, inferred round, available labor or mode fields, read depth, the proportion assigned to genes when available, and all clinical or RNA quality flags. Keep flagged pregnancies in the primary analysis under their released labels. Record results excluding the flagged FGR donors as a sensitivity check. Never infer age at diagnosis from age at delivery.
- [x] Write `cohort_manifest.csv` and `cohort_support.csv` with one row per sample in the manifest and one row per cohort and comparison group in the summary. Include sample count and delivery age minimum, median and maximum. Verify the primary 7 vs 5 and 2 vs 5 counts and the shared delivery range **26.714–30.429 weeks**. Show each donor's delivery age because group ages can still differ within the overlap.
- [x] Pin `REACTOME_SIGNALING_BY_BMP` (R-HSA-201451) from [MSigDB Reactome v2024.1.Hs](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/c2.cp.reactome.v2024.1.Hs.symbols.gmt), SHA-256 `9ea1b5e656597daf423e41c5ebcaa9892bfedf3292fff768605d7b0d5e5e9703`. Start from its 31 genes. Use the organizer's annotation and unique HGNC/Ensembl mappings for score calculation. Require at least 10 raw reads in five or more libraries in **each** full villous cohort before looking at effect directions. This yields **27 score genes**. Record source ID, mapped ID, coverage and exclusion reason for all 31 in `bmp_gene_manifest.csv`. The four score exclusions are AMHR2, BMP10, CER1 and GDF2.

**Check:** Fail on a sample mismatch, duplicate join or unexpected group count. If files differ from the audited release, show the difference and resolve it before interpreting any estimate. Plot the age distribution and variance-stabilized PCA by diagnosis, fetal sex and inferred round to identify design limitations and sample anomalies.

## 2. Run DESeq2 and quantify the pathway

- [x] Reuse the full OpenBio villous DESeq2 model `~ fetal_sex + group` on raw counts and extract EOPE+FGR versus isolated FGR. All five released groups can inform the common normalization and dispersion estimates. Use the same model for the age-overlap refit if its design is full rank. Export signed Wald statistics, unshrunken log₂ fold changes, approximate Wald intervals and contrast-wide gene FDR. [DESeq2 documentation](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)
- [x] Compute the **primary donor score**. Estimate DESeq2 median-ratio size factors from the full count matrix separately within each cohort. For each pregnancy take the unweighted mean of `log2(normalized_count + 1)` across the 27 fixed members. Save `bmp_donor_scores.csv` with cohort, donor, clinical group, score and age where available. For each cohort, compute the case mean minus the FGR mean with an approximate 95% Welch interval. Show every donor. The unadjusted score estimates describe pregnancy-level differences, and the sex-adjusted count model answers the complementary gene-level question.
- [x] Reuse the same score genes and size factors in all local subsets and single-donor omissions. Report the full comparison, the 2 vs 5 age-overlap comparison, the leave-one-pregnancy-out range and the sensitivity excluding the two FGR donors with BP flags individually and jointly. Report the age-restricted interval even if it is wide. Do not treat a comparison with two cases as a well-powered test or fit a multivariable donor-score model in that subset.
- [x] Rank genes by signed Wald statistic, require finite statistics and p values and `baseMean ≥ 10`, drop ambiguous symbol mappings, and use the intersection of eligible symbols across the fixed rank fits. The pilot universe intersects five full contrasts in each cohort (each EOPE group versus isolated FGR, each EOPE group versus preterm controls, and FGR versus preterm controls) plus the local age-restricted primary contrast. Use `fgseaMultilevel` with seed 1049, serial execution, `minSize=10`, `maxSize=500` and `eps=1e-10` for the three Reactome branches. Record the exact rank universe and the **25 BMP members in the pilot rank universe**. This differs from the 27 score genes because the rank filter is stricter. Report normalized enrichment, leading-edge genes and BH correction across the three branch tests within each contrast. Activin and TGF‑β receptor signaling provide context for how specifically the RNA pattern points to BMP.
- [x] Recalculate the score without INHA and INHBA as a sensitivity check for shared BMP/activin annotation. Keep the full 27-member score as the primary estimate. Show whether this changes the direction or magnitude. Do not substitute a score restricted to the seven genes selected from the enrichment results.

**Check:** Export `bmp_effects.csv`, `bmp_gene_effects.tsv.gz` and `branch_enrichment.csv`. Positive effects must use the same group order in the score and count results. Reproduce the pilot donor score differences to rounding in the full and age-overlap cohorts (**+0.273** and **+0.264**). Keep enrichment distinct from replication, since both summaries use the same pregnancies.

## 3. Compare with independent data and decide what to claim

- [x] Reproduce Awamleh from the authors' public transcript-feature counts using the existing Ensembl 75 aggregation and audited GEO metadata. Confirm 20 vs 18 and the expected four group counts. Score the same 27 genes and run the group-only Awamleh count contrast. Document the use of all released samples. Repeat on the authors' QC-filtered subset only if the retained GSM IDs can be verified. [Awamleh paper](https://doi.org/10.1186/s12920-019-0548-x) · [GEO release](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114691)
- [x] Report the Awamleh primary comparison and two context comparisons, EOPE-only versus isolated FGR and isolated FGR versus preterm controls. Keep cohort scales separate and compare **within-cohort differences**. The pilot differences were +0.145, +0.168 and +0.146, respectively. The last contrast establishes that elevation can accompany FGR by itself. Keep the full OpenBio EOPE-only and term comparisons out of the primary inference because they add differences in sequencing round or delivery age.
- [x] If time permits, reproduce the 27-gene score with `+1` in Piekos using the severe PE+FGR and isolated FGR labels. Report the unadjusted comparison and the age-adjusted comparison within observed support. The same-score check gave an age-adjusted +0.141 (95% CI −0.021 to +0.303 with HC3 standard errors) in 15 vs 16 pregnancies at 32.71–38.29 weeks. Present it as severe PE context rather than an age-matched EOPE replication. Retain the earlier 31-gene estimate as a separate sensitivity. [Piekos paper](https://doi.org/10.1038/s42003-025-08631-6) · [Public data and code](https://github.com/spiekos/Hood-Lab-Placental-Multiomics)

| Observed pattern | Report language |
| --- | --- |
| Both main cohorts point up and leave-one-out checks preserve direction. | Higher BMP-associated RNA in EOPE+FGR than isolated FGR in two cohorts. State the effects and the delivery-age uncertainty. |
| Age restriction or one-donor omission reverses the estimate. | Group difference is sensitive to cohort composition. Show the reversal with each donor and the age distribution. |
| Awamleh is null or changes direction after a justified sensitivity check. | The comparison gives mixed evidence across cohorts. Retain the clinical question and report the uncertainty. |
| Activin and BMP have similar changes or the BMP finding shrinks after removing shared members. | Broader TGF‑β superfamily expression is a plausible explanation. Keep a specific BMP activity claim for a follow-up study. |

The same-direction estimates in the observed pilot support following the comparison. Similar changes in both EOPE strata leave a combined-phenotype interaction unresolved. Any test of such an interaction would require adequate support across all four preterm groups within the same age and sequencing context.

## 4. Explain the biology after the pathway result

- [x] Look up **BMPR2**, **BMPR1B** and **GREM2** in the official Zeisel/Admati browser. The first two encode type II and type I receptors. GREM2 is an extracellular BMP antagonist. All three contribute to the positive leading edge across the full local, age-restricted and external fits. Use **SMAD9** as an additional candidate if it helps explain the canonical signaling branch. Save `browser_lookups.csv` with gene, highest expressing cell classes, breadth, visible disease/control pattern and interpretation. [Browser guide](../../../openbio-datathon-2026/track-01/data/interpretation/zeisel_browser_guide.md) · [Reactome BMP pathway](https://reactome.org/content/detail/R-HSA-201451)
- [x] For quantitative expression context, use only the released Admati donor-by-cell-class pseudobulk. Require at least 20 cells per donor-class and at least three eligible donors per class to summarize it. Plot median CPM and the fraction of donors with at least 10 counts. Show contributing donor counts, including classes missing from the reference. A dot plot can include all shared leading-edge genes. Independently compare a small set of published class markers in bulk to see whether shifts in cell proportions could help explain the result. The atlas provides plausible source cells and the marker check provides context. Specific within-cell diagnosis effects require like-for-like donor groups and comparable delivery context. [Admati et al.](https://doi.org/10.1016/j.medj.2023.07.005)
- [x] Write one hypothesis card. Predict that the higher villous RNA score will accompany altered **nuclear pSMAD1/5/9** in defined cell types in a new cohort with the two FGR arms and comparable delivery age. Quantify cell proportions in the same sections and use pSMAD2/3 as a complementary readout because BMP can also signal through that route in trophoblasts. A shift explained by cell proportions, with similar pSMAD1/5/9 within cell types, would weaken a canonical within-cell signaling interpretation. Cite previous BMP work on EOPE and trophoblasts. [Yi et al.](https://doi.org/10.3389/fcell.2021.768669) · [Deng et al.](https://doi.org/10.1016/j.ebiom.2023.104664)

**Optional paired analysis:** Use the four combined and five isolated FGR pregnancies with both villi and decidua. Normalize within each tissue. Score the same shared members in both tissues, with an explicit coverage audit, then compare the paired decidua-minus-villus score between diagnoses. Report both tissue effects and the difference of effects with intervals. A comparison of which tissue has the smaller p value cannot establish a tissue difference. This extension follows the primary result if there is time.

## 5. Assemble the datathon story and rerun it

- [x] Produce at least four article figures. Use the same `figures/headline.png` in the article and the video. Give each axis, group label and interval a plain-language explanation.

  1. **Clinical group and age support.** One dot per pregnancy, with the shared delivery-age region marked.
  2. **Where BMP appears in the DESeq2 output.** Show the ranked enrichment and annotate a few contributing genes such as BMPR2. Report direction and leading-edge membership.
  3. **Headline.** Show pregnancy-level BMP scores in OpenBio and Awamleh as separate panels, with sample counts. Accompany them with the within-cohort effect estimates and a visibly wide age-restricted interval.
  4. **Where the transcripts are expressed.** Show the atlas expression of receptors, mediators and regulators by donor and cell class, with a caption that states what the reference contains.

- [ ] Put the runnable analysis and setup instructions in `notebooks/` and `README.md`. Include the analysis plan and one hypothesis card. Write `docs/article.md` and prepare a video of at most five minutes, using the official unlisted-video workflow. Verify the team folder against the [submission checklist](../../../openbio-datathon-2026/track-01/4-submit.md). In the article, explain that the exact comparison was used before and that our contribution is measuring the magnitude and consistency of the additional BMP-associated RNA shift. [Ackerman et al.](https://doi.org/10.1161/HYPERTENSIONAHA.122.20807) · [Aisagbonhi et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC13098728/#S12)
- [x] Run once into an empty output directory. Verify one score per pregnancy, all figure source tables present, matching estimates and signs across figures and prose, and a report of who was included in every comparison. Check that paths are relative to the project, public input downloads have checksums and a rerun can recreate the headline without manual plotting. Keep the clinical and technical limits in both the article and the video.

The core submission path is Steps 1–5 through the main two cohorts and the browser. The optional Piekos and paired tissue extensions can wait until the core figures and hypothesis card are reviewable.

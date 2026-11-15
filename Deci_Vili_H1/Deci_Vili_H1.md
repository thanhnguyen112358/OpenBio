# Feasibility plan for the paired decidua–villus hypothesis

> At delivery, pregnancies with early-onset preeclampsia (EOPE) without fetal growth restriction (FGR) will show lower C1Q–TREM2 expression in the decidua basalis than pregnancies with isolated FGR. The difference will be larger in decidua than in villi from the same pregnancies and larger than the shift in broad myeloid markers in decidua.

The pilot tests a paired-tissue RNA association at delivery. A selective decidual signal would motivate follow-up measurements of macrophage abundance, subtype composition, and uptake of apoptotic trophoblasts per macrophage. [C1Q–TREM2 signature in uterine macrophages](https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2023.1256453/full) · [Human trophoblast uptake preprint](https://doi.org/10.64898/2026.02.26.708195)

## 1. Assemble and audit the paired cohort

1. Load the released villous tissue (`VT`) and decidua basalis (`DB`) count matrices. Check count columns against `sample_id` in the corresponding metadata and join the two tissues by `donor_id`. Confirm that each pregnancy has at most one sample per tissue and that the two records agree on phenotype and delivery age.
2. Use the released `EOPE` and `FGR` labels for the primary comparison. Keep the EOPE+FGR, idiopathic preterm birth, and term birth groups available for context. Of the 30 pregnancies, 27 have complete tissue pairs. All three missing decidua samples belong to the EOPE+FGR group.
3. Make a donor audit table before inspecting group expression effects. Include delivery age, inferred sequencing round, clinical flags, and gene-assigned read counts, assignment rate, and panel gene coverage for each tissue. The organizer flags the recorded maximum blood pressure for FGR donors `S04` and `S15`. Retain their diagnoses and assess their influence in sensitivity analyses. The released bulk metadata lacks donor-level labor and delivery-mode fields.

| Analysis set | EOPE without FGR | Isolated FGR | Role |
| --- | ---: | ---: | --- |
| Inferred same sequencing round | 4 | 5 | Primary pilot estimate |
| Inferred same round and shared delivery-age range, 26.714–30.429 weeks | 3 | 5 | Delivery-age sensitivity |
| All complete pairs in the two groups | 8 | 5 | Full-cohort context |

Ackerman et al. report two sequencing rounds and identify four EOPE-only donors sequenced with the preterm and term controls (`S16`, `S20`, `S21`, and `S22`). This implies that the remaining four EOPE-only donors and all five FGR donors share the other round. Confirm the assignment for each tissue against Table S2 if available and label the assignment as inferred until it is verified. The shared-age subset excludes all four EOPE-only donors sequenced with the controls, so changing directly from 8 vs. 5 to 3 vs. 5 changes both the age and sequencing-round composition. [Ackerman et al., Methods](https://pmc.ncbi.nlm.nih.gov/articles/PMC10192030/) · [Data package](../openbio-datathon-2026/track-01/data/README.md) · [Data dictionary](../openbio-datathon-2026/track-01/data/track01_data_dictionary.md) · [Sequencing QC](../openbio-datathon-2026/track-01/data/qc/track01_htseq_summary.csv)

## 2. Calculate two fixed RNA scores in both tissues

1. Estimate DESeq2 median-ratio size factors within each tissue, using genome-wide raw counts for the same 13 EOPE-only and isolated-FGR pregnancies. Transform each normalized count as `log2(normalized count + 1)`. Use the same genes and weights in both tissues. Keep the normalized values fixed when comparing the three donor subsets.
2. Calculate both scores for each donor and tissue using the transformed values.

| Score | Definition |
| --- | --- |
| C1Q–TREM2 | `0.5 × TREM2 + 0.5 × mean(C1QA, C1QB, C1QC)` |
| Broad myeloid | `mean(CD68, LST1, CSF1R)` |

3. Plot the `TREM2` and C1Q components separately. Giving the three C1Q subunits half the total weight keeps them from dominating the score. Check coverage again in the smaller primary subset, especially for villous `TREM2`, which reached 10 raw reads in 10 of the 13 villous libraries. The broad myeloid score reflects bulk myeloid transcription and can vary with both cell composition and expression per cell. [DESeq2 normalization and contrasts](https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

## 3. Estimate the paired effect and inspect its components

1. **Primary contrast.** For each pregnancy, subtract the villous C1Q–TREM2 score from the decidual score. Compare the mean of these paired differences in EOPE against the mean in FGR. Report the signed difference (EOPE minus FGR), an approximate Welch 95% confidence interval, and all donor points. The predicted contrast is negative.
2. **Tissue effects.** Report the EOPE-minus-FGR contrast for the C1Q–TREM2 score separately in decidua and villi, using the same donors. The hypothesis predicts a negative decidual effect that is more negative than the villous effect. A negative paired contrast can also arise from higher villous expression, so inspect both estimates.
3. **Selectivity within decidua.** For each donor, subtract the broad myeloid score from the C1Q–TREM2 score in decidua. Compare that difference between EOPE and FGR. The predicted contrast is negative. Show the individual marker effects alongside this contrast so readers can assess whether the candidate and broad marker scores move together.
4. **Sensitivity analyses.** Repeat the estimates in the shared-age and full cohorts. In the primary cohort, omit one pregnancy at a time and report the range of estimates and any sign reversals, including when each BP-flagged donor is omitted. Repeat scoring with tissue-specific trimmed mean of M-values (TMM) normalization. Flag dependence on the normalization method or on low villous `TREM2` read counts. [TMM normalization method](https://link.springer.com/article/10.1186/gb-2010-11-3-r25)

## 4. Add villous stress and cell-type context

- Plot the decidual C1Q–TREM2 score against one prespecified villous `HALLMARK_HYPOXIA` score, with one point per pregnancy. Use the plot to assess overlap in hypoxia-associated transcription. Similar scores describe overlap on this RNA signature. Placental injury burden and local trophoblast death require separate measurements. Leave phenotype-specific slopes for a larger cohort. [MSigDB Hallmark Hypoxia](https://gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_HYPOXIA.html)
- Look up `TREM2`, `C1QA`, and `C1QC` in the required Zeisel/Admati browser and record which cell classes express them. The optional placental pseudobulk can add cell-class context. It contains no isolated-FGR group, so the paired EOPE–FGR association requires a separate cohort for replication. [Browser guide](../openbio-datathon-2026/track-01/data/interpretation/zeisel_browser_guide.md) · [Pseudobulk donor metadata](../openbio-datathon-2026/track-01/data/extension/pseudobulk/track01_pseudobulk_donors.csv)

## 5. Interpret the feasibility result

| Pilot result | Interpretation |
| --- | --- |
| Lower decidual C1Q–TREM2, a smaller villous shift, and a smaller broad myeloid shift, with the direction stable in the same-round and shared-age analyses | Candidate selective decidual RNA difference for follow-up in a larger cohort |
| Candidate and broad myeloid markers shift together | Bulk myeloid difference compatible with a change in cellular composition |
| The direction reverses after omitting one donor or restricting to shared ages | Estimate sensitive to the cohort definition or an individual pregnancy |
| The interval spans substantial effects in both directions | The cohort leaves the hypothesis unresolved |

Use the paired donor contrast as the headline figure if the pattern holds. Build three supporting graphs showing cohort overlap, separate tissue and marker effects, and sensitivity estimates. A follow-up cohort could compare decidual macrophage abundance, subtype composition, and apoptotic trophoblast uptake per macrophage between EOPE-only and isolated-FGR pregnancies matched on delivery age and labor. [Track 1 submission guide](../openbio-datathon-2026/track-01/4-submit.md)

**Preflight status.** The metadata joins, paired donor counts, delivery-age overlap, and pooled panel coverage have been checked. All seven genes reach at least 10 raw reads in at least 10 of the 13 target pregnancies in each tissue. The sample-level sequencing assignments and phenotype effect estimates are pending.

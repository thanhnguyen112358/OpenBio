# Step 4 villous hypoxia and cell-type context

## Overlap in villous hypoxia-associated transcription

The 4 EOPE and 5 isolated-FGR pregnancies in the inferred same sequencing round overlap on the villous `HALLMARK_HYPOXIA` score from **−0.160 to −0.064 centered log₂ score units**. The overlap contains **2 of 4 EOPE** pregnancies (`S06` and `S19`) and **2 of 5 isolated-FGR** pregnancies (`S04` and `S08`). The same four pregnancies are inside the overlap in both age-restricted subsets. The full cohort adds one EOPE pregnancy (`S20`) in that range, for 3 of 8 EOPE and 2 of 5 FGR. [Donor scatterplot](figures/hypoxia_vs_decidual_score.png) · [Overlap by cohort](hypoxia_overlap.csv)

Decidual C1Q–TREM2 scores interleave by phenotype within the shared band. In ascending order they are 7.28 (EOPE), 7.52 (FGR), 8.32 (EOPE), and 8.76 (FGR). This is a descriptive overlap on one RNA signature at delivery. A larger cohort with more shared score coverage would be needed to estimate a phenotype-specific relationship. Placental injury burden or local trophoblast death would require separate measurements.

### Scoring

The organizer pins MSigDB release **2024.1.Hs**. The download endpoint in the organizer script redirects to sign-in, so the analysis uses the identical release filename from Broad Institute's public archive. All **200** hypoxia members map to unique villous count rows. Of those, 194 are direct gene-symbol matches and six use a unique previous symbol in the organizer's HGNC annotation. There is no gene overlap with the seven-gene C1Q–TREM2 and broad myeloid panel. [Broad archive](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/) · [Hallmark Hypoxia](https://gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_HYPOXIA.html) · [Archive checksum and six symbol mappings](run_info.json)

For each gene, divide the raw villous count by its sample's fixed Step 2 size factor, take `log2(normalized count + 1)`, and subtract the mean of that gene across the 13 target pregnancies. The hypoxia score is the unweighted mean of those 200 values. The decidual values come directly from Step 2. The same hypoxia values are used in every subset. The plot shows one point per pregnancy and the shaded band gives the intersection of the two groups' observed ranges in the primary cohort. [Step 2 normalization](../step02_deseq2/README.md) · [Scores by donor](hypoxia_by_donor.csv)

## Zeisel/Admati browser lookups

The three genes were entered and viewed in the live browser on November 15, 2026. The observations below describe the visible violin and expression-map patterns. [Zeisel/Admati gene browser](https://storage.googleapis.com/www_zeisellab/PE_web/pe_web_test_v1.html)

| Gene | Qualitative browser observation | Organizer pseudobulk median CPM in nine eligible control donors per class |
| --- | --- | --- |
| `TREM2` | Strongest signal in the `MAC-TREM2` cluster, with weaker signal in other myeloid clusters including Hofbauer cells | Macrophage 248.5, Hofbauer 8.8 |
| `C1QA` | Prominent in Hofbauer and Hofbauer-IFN1 clusters and in macrophages, including `MAC-TREM2` | Hofbauer 4,298.9, macrophage 1,334.6 |
| `C1QC` | Similar localization to `C1QA`, with prominent Hofbauer, Hofbauer-IFN1, and macrophage signal | Hofbauer 1,857.7, macrophage 664.0 |

The quantitative context comes from the organizer's released donor-by-cell-class pseudobulk. It uses six late and three early non-PE controls, with one observation per donor and broad class and at least 10 captured cells per donor–class library. CPM is the gene count divided by the library's total counts and multiplied by one million. The summary takes the median across donors and reports the number of donors with a positive count. Both stages are pooled for this localization check. Some classes have fewer eligible donors. [All cell-class summaries and donor counts](cell_class_expression.csv)

The browser and pseudobulk draw on the same single-cell study. Admati and colleagues identified maternal-origin TREM2-positive macrophages and fetal Hofbauer cells in their atlas. The different distributions of `TREM2` and the two C1Q genes leave multiple plausible cellular sources for the bulk score. The reference lacks an isolated-FGR arm and cannot assign the paired decidual difference to either population or measure trophoblast clearance. [Admati and Skarbianskis et al.](https://pubmed.ncbi.nlm.nih.gov/37572658/) · [Organizer pseudobulk source](../../../openbio-datathon-2026/track-01/data/README.md)

## Reproduce

Use the Step 2 environment from the workspace root. The script reads the organizer's counts and annotation directly, checks their release checksums, and retrieves the pinned GMT in memory. It writes the tables, figure, and provenance under this folder. The qualitative browser observations above were recorded separately. [Analysis script](../../scripts/step04_hypoxia_context.py)

```powershell
& 'Deci_Vili_H1/.venv/Scripts/python.exe' Deci_Vili_H1/scripts/step04_hypoxia_context.py
```

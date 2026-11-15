# Genome-wide clinical contrasts

The analysis extracts ten unordered clinical contrasts from one R DESeq2 model per tissue. The inputs are the organizer's Track 01 raw counts for 30 villous tissue (VT) libraries and 27 decidua basalis (DB) libraries. A positive log₂ fold change means higher expression in the `numerator` group. The models are `~ fetal_sex + group` with DESeq2's default size-factor, dispersion, and Wald test procedures. They use the full set of pregnancies available in each tissue.

Ackerman et al. show nine clinical contrasts per tissue in supplementary Figures S1 and S2. Those figures omit idiopathic preterm birth versus term birth, a comparison reported in the authors' earlier study. Their 2023 models also include two surrogate variables. This run is an exploratory baseline from the released files. Gestational age and sequencing round can influence the estimates. Age overlap in the summary describes group support and does not restrict the model fit. [Ackerman et al., supplemental figures](https://www.medrxiv.org/content/medrxiv/early/2023/03/12/2023.03.08.23286998.full.pdf#page=48) · [Earlier control study](https://pubmed.ncbi.nlm.nih.gov/27452435/)

- `results/vt_contrasts.tsv.gz` and `results/db_contrasts.tsv.gz` contain the complete gene-by-contrast grids. `low95` and `high95` are unadjusted Wald intervals for the unshrunk effect. `padj` applies Benjamini–Hochberg within the contrast after DESeq2's independent filtering. The `status` column flags filtered rows, unavailable tests, and genes with no counts in either comparison group or in the entire tissue.
- `results/contrast_summary.csv` gives sample counts, delivery-age overlap, and the number of higher and lower genes at contrast-specific FDR < 0.05.
- `results/session_info.txt` records the package versions used.

Look across effect sizes and intervals for patterns worth following. A difference against term births includes delivery-age differences. The contrasts reuse pregnancies and do not provide independent replication. A direct tissue-by-group test would be needed to estimate whether a clinical association differs between VT and DB.

See [the follow-up results](follow_up/README.md) for age and inferred sequencing-round restrictions, pathway enrichment, and the comparison with a separate placental cohort.

Run from the repository root with R and Bioconductor DESeq2 installed.

```powershell
Rscript Pairwise_DESeq2/run_pairwise.R openbio-datathon-2026/track-01/data Pairwise_DESeq2/results
```

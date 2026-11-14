# Q1 · What changes in placental gene expression in preeclampsia?

**Working hypothesis.** Early-onset preeclampsia (EOPE) without fetal growth restriction (FGR) will show an antiangiogenic and stress-related RNA signature in villous tissue at delivery. Expect the angiogenic signal to recur across more donors than the metabolic signal. Ackerman’s source study found a glycolytic expression pattern in a subset of EOPE placentas. Maynard’s experiments give a mechanistic reason to look at `FLT1`, and placental deconvolution work gives a reason to track cell mixture as a competing explanation. [Ackerman et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10192030/) · [Maynard et al.](https://www.jci.org/articles/view/17189) · [Campbell et al.](https://www.nature.com/articles/s42003-023-04623-6)

## Cohort and contrasts

Use the GSE203507 core villous counts. Keep `EOPE+FGR` (7), isolated `FGR` (5), and term births, `TB` (5), visible as distinct groups. Only 2 of the 7 combined cases fall within the isolated-FGR group’s delivery-age range. [Core counts](openbio-datathon-2026/track-01/data/track01_counts_vt.csv) · [Villous metadata](openbio-datathon-2026/track-01/data/track01_metadata_vt.csv)

### Primary comparison decision

**November 14, 2026.** Compare placental villous RNA at delivery in **8 pregnancies with early-onset preeclampsia without fetal growth restriction (`EOPE`) versus 5 pregnancies with spontaneous idiopathic preterm birth (`iPTB`)**. Use the **5-vs.-5 subset** within **29.857–33.000 weeks** (29 weeks 6 days through 33 weeks) for delivery-age sensitivity.

**Why this cohort.** Both groups delivered preterm and had no diagnosed FGR, which can carry its own placental expression changes. The iPTB delivery-age range contains five of eight EOPE pregnancies without FGR and all five iPTB controls. Only one of seven EOPE+FGR pregnancies falls in that window, so pooling both EOPE groups would change the case mix when restricting by age. The five term births delivered at 39.0–39.4 weeks, beyond the delivery age of every EOPE case. [Villous metadata](openbio-datathon-2026/track-01/data/track01_metadata_vt.csv) · [Delivery-timing audit](Q06-delivery-timing.md)

**What PCA added.** The organizer-style `log2(CPM + 1)` PCA showed a broad shift along PC1 between EOPE without FGR and iPTB after removing term births. The shift was visible again when fitting PCA on the 8-vs.-5 comparison and on the 5-vs.-5 subset. PC1 explained **25.1%** and **29.3%** of within-fit expression variation, respectively. Four of five iPTB controls fell below the EOPE range on PC1 in both fits. The fifth, S24, overlapped it. Seeing a group-level shift in both fits modestly strengthens the case for testing this comparison at the gene level. [Preterm PCA](PCA/results/preterm/README.md) · [Full comparison PCA](PCA/results/eope_vs_iptb/README.md) · [Shared-age PCA](PCA/results/shared_age/README.md)

**What to carry into the analysis.** PC2 correlated with log gene-assigned read depth in the full and shared-age fits (absolute Pearson `r = 0.87` and `0.86`). S22, the lowest-depth pregnancy in the comparison, lay at a PC2 extreme in both. Check sequencing QC and donor influence, then compare gene-level effect estimates and intervals across the unadjusted and delivery-age-adjusted 8-vs.-5 fits and the 5-vs.-5 sensitivity fit. Interpret those estimates as associations relative to spontaneous iPTB at delivery, where labor and delivery context may also contribute. [PCA results](PCA/README.md)

## Analysis

1. Verify one sample per donor, count-column alignment, and sequencing QC. Plot delivery age across all five groups. Recolor the same PCA by group, delivery age, and fetal sex. Inspect any influential sample before deciding whether it warrants a sensitivity run.
2. Fit DESeq2 on raw counts for the full 8-versus-5 comparison with `~ group`. Then fit `~ ga_delivery_weeks + group` on the same donors and rerun the unadjusted model on the 5-versus-5 overlap subset. The adjusted full model relies on a common linear age trend and includes cases outside the control range. Check fetal sex separately if estimable. Record log₂ fold changes, intervals, and false discovery rate (FDR)–adjusted p-values.
3. Prespecify `FLT1`, `ENG`, `LEP`, and `INHBA` as gene benchmarks. Rank tested genes by signed Wald statistic for Hallmark Hypoxia, Glycolysis, and Oxidative Phosphorylation. Plot normalized donor expression and the pathway’s leading-edge genes. Check whether the angiogenic and metabolic patterns occur in the same donors.
4. Compare `(EOPE+FGR) − EOPE` and `FGR − iPTB` for directional context, and show how using `TB` changes the result. A metabolic signal shared with FGR would support a common placental stress hypothesis. An angiogenic effect that persists across the primary fits and is weaker in FGR alone would support an EOPE-associated component in this cohort.
5. Look up at least three result genes in the Zeisel browser. Select an angiogenic gene, a leading-edge gene, and one that helps distinguish competing explanations. If a cell class is clearly implicated, check the direction in Admati donor pseudobulk. The early arm starts with 10 cases and 3 controls. [Organizer workflow](openbio-datathon-2026/track-01/3-workflow.md) · [Browser guide](openbio-datathon-2026/track-01/data/interpretation/zeisel_browser_guide.md)

## Figures and follow-up

Build four figures covering delivery age and PCA, effect sizes across the three primary fits, donor-level expression and pathway patterns across all groups, and cell localization of three result genes. Center the hypothesis card on the pattern with the most stable effect sizes. Agreement with Ackerman benchmarks reproduces findings from the source cohort. For an external comparison, audit tissue and control comparability in GSE114691, which Ackerman also used. A claim about soluble sFLT1 requires isoform-specific RNA or protein measurements. [GSE114691](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114691)

Planning chat: “Placental Gene Expression Analysis.” [OpenBio index](Thanh/README.md)

## November 14, 2026 morning review

Reviewed at 9:35 a.m. EST after the November 14 morning organizer repo update to [`fb751e9`](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/commit/fb751e9).

The initial plan already used the released villous data and proposed four figures. The morning review verified the **8 EOPE without FGR vs. 5 iPTB** primary contrast and **5 vs. 5** shared-age subset (29.857–33.000 weeks). The updated submission guide asks for a general-audience article, so this is the proposed main question and the figures should carry one coherent story. [Released data dictionary](openbio-datathon-2026/track-01/data/track01_data_dictionary.md) · [Updated submission guide](openbio-datathon-2026/track-01/4-submit.md)

- Report the full unadjusted fit alongside the delivery-age-adjusted fit on the same 13 donors and the unadjusted 5-vs.-5 fit. Use [Q6](Q06-delivery-timing.md) for comparator context, with the term contrast clearly labeled as having no delivery-age overlap. Prespecify the comparison and thresholds before inspecting the gene results.
- Export four distinct graphs for the article. Show delivery-age overlap, sample PCA, effect sizes and intervals across the fits, and donor-level expression for the lead finding. Pick the headline graph after seeing the results. Save it as `figures/headline.png` for the article and video.
- Use [Q9](Q09-baseline-gene-localization.md) to record the three required browser lookups and explain what they suggest about cell source. Build the hypothesis card from the observed result, a plausible mechanism, an alternative such as cell mix or delivery timing, and a follow-up that could falsify it. Cite the studies behind the mechanism.

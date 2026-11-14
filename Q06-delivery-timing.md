# Q6 · How sensitive is the placental expression signal to delivery timing?

**Working hypothesis.** Some apparent PE expression differences will attenuate when term controls are replaced by preterm controls, while a subset of angiogenic or stress signals may persist in the delivery-age overlap. Maturation, labor, and spontaneous preterm-birth biology can all contribute to the comparator effect. Delivery age can also reflect disease severity and clinical decisions. The estimand is the observed sensitivity to comparator choice and delivery-age support. [Wang et al.](https://doi.org/10.1038/s41586-026-10316-x) · [Cindrova-Davies et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC1988867/)

## Comparison sequence

Define pooled `PE` from `EOPE` and `EOPE+FGR` while preserving the original five group labels. The 15 cases averaged 30.21 weeks at delivery, compared with 31.54 for `iPTB` and 39.17 for `TB`. Six cases overlap the 29.857–33.000-week `iPTB` range, including only one with FGR. Restriction therefore changes the case mix as well as the age support. [Villous metadata](openbio-datathon-2026/track-01/data/track01_metadata_vt.csv)

| Fit | Donors | Purpose |
| --- | --- | --- |
| PE versus term birth | 15 vs. 5 | Estimate the effect with a large delivery-age gap |
| PE versus idiopathic preterm birth | 15 vs. 5 | Change the comparator while keeping the cases fixed |
| PE versus iPTB in the shared window | 6 vs. 5 | Estimate effects among donors with overlapping delivery ages |
| iPTB versus term birth | 5 vs. 5 | Identify signals that also differ between the control groups |

## Analysis

1. Check donor uniqueness, count alignment, and QC. Plot each donor’s delivery age, with overlap shaded. Recolor one PCA by age, group, FGR, and sex. Plot PC1 and PC2 against age in the full cohort and inspect the preterm subset. Keep the released value for flagged term sample `S26` and check its influence.
2. Fit the four contrasts using the organizer’s raw-count DESeq2 workflow. On the same 11 overlap donors, compare `~ condition` with `~ ga_delivery_weeks + condition` to assess residual age adjustment with a fixed donor set. Also compare the five overlapping EOPE cases without FGR with the five `iPTB` controls. Keep isolated FGR as a separate contextual group. PE and term controls have zero overlap, so adjusting that contrast for age would rely on extrapolation. [Organizer workflow](openbio-datathon-2026/track-01/3-workflow.md)
3. Make a scatterplot of log₂ fold changes against term and preterm controls using a common adequately expressed gene set. Show estimates with intervals across all three PE fits for `FLT1`, `PGF`, `ENG`, `LEP`, and prespecified hypoxia and glycolysis programs. Quantify attenuation in log₂ units. Report percent attenuation only when the original estimate is meaningfully away from zero and the sign is preserved. Flag reversals. Similar estimates with wider intervals indicate lost precision.
4. Look up at least three result genes in the Zeisel browser and inspect identity markers for the affected classes. Interpret attenuation alongside labor, indication for delivery, and FGR, which this design can only partly account for. [Browser guide](openbio-datathon-2026/track-01/data/interpretation/zeisel_browser_guide.md)

## Figures and follow-up

Present age overlap, PCA, and a side-by-side effect plot. Report which signals weaken, persist, or become imprecise under each comparison. To test whether the strongest signals recur, use an independent cohort with overlapping delivery ages and documented labor, delivery indication, and FGR. Audit those fields before using GSE114691. [GSE114691](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114691)

Planning chat: “Compare Delivery Timing Effects.” [OpenBio index](Thanh/README.md)

## November 14, 2026 morning review

Reviewed at 9:35 a.m. EST after the November 14 morning organizer repo update to [`fb751e9`](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/commit/fb751e9).

The morning review verified the original pooled comparison. PE and term controls have no delivery-age overlap. The pooled PE vs. iPTB window retains **6 vs. 5** donors, including only **one** EOPE+FGR case. Restricting the pooled analysis therefore changes case composition as well as age support. [Released villous metadata](openbio-datathon-2026/track-01/data/track01_metadata_vt.csv)

For the recommended [Q1-led analysis](Q01-placental-gene-expression.md), keep EOPE without FGR as the case definition for the main **8-vs.-5** comparison and **5-vs.-5** overlap sensitivity. Show the age-adjusted fit alongside the unadjusted fit on the same donors so readers can distinguish model adjustment from donor restriction. Use Q6's term contrast to illustrate how estimates vary with the comparator. The change can reflect maturation, delivery context, or case mix. Show effect sizes and intervals across fits so attenuation and loss of precision are distinguishable.

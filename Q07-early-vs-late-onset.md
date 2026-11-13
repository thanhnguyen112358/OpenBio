# Q7 · Do early- and late-onset preeclampsia show different placental responses?

**Working hypothesis.** Early PE may have larger angiogenic and stress effects in syncytiotrophoblast and inflammatory effects in macrophages. Late PE may show a smaller shared response or a selective myeloid signal. Admati supplies the source-cohort benchmark. Lee’s independent preprint motivates retaining monocytes as a secondary class. [Admati et al.](https://pubmed.ncbi.nlm.nih.gov/37572658/) · [Lee et al., preprint](https://pubmed.ncbi.nlm.nih.gov/41472684/)

## Cohort and contrasts

Use the Admati pseudobulk extension and its released `clinical_group` labels. `ga_delivery_weeks` is the age at sampling, and the release has no separate diagnosis-week field. The early comparison has 10 PE cases and 3 controls, and the late comparison has 7 and 6. Only 7 early and 2 late cases fall within their respective control delivery-age ranges. Nine early cases delivered by cesarean and all early controls delivered vaginally. FGR occurs in 6 early cases and 2 late cases, and all 7 late cases are recorded as induced. [Donor metadata](openbio-datathon-2026/track-01/data/extension/pseudobulk/track01_pseudobulk_donors.csv)

For each cell class, estimate `PE_early − control_early`, `PE_late − control_late`, and the direct difference `(PE_early − control_early) − (PE_late − control_late)`. The third contrast measures the difference between the two observed case–control effects. Its interval is the basis for assessing an onset-related difference.

## Analysis

1. Start with `syncytiotrophoblast`, `hofbauer_cell`, `macrophage`, and `endothelial_cell`. Use ≥20 cells per profile with ≥10 and ≥30 sensitivity checks. Most priority classes retain all donors at ≥20 cells. Macrophages retain 9 early cases and all 3 early controls. EVT is descriptive with 1 early and 2 late controls at this threshold. [Profile metadata](openbio-datathon-2026/track-01/data/extension/pseudobulk/track01_pseudobulk_metadata.csv)
2. Audit ID alignment and library depth, then plot delivery age, mode, induction, and FGR by donor. Use within-class PCA to inspect the three early controls and unusual profiles. Fit one raw-count DESeq2 model `~ clinical_group` per class. Use the same expressed gene universe for its three contrasts. Report signed log₂ effects, intervals, and the scope of FDR correction.
3. Lead with SCT `FLT1` and `PGF`, followed by SCT hypoxia and angiogenesis, Hofbauer and macrophage inflammatory or interferon programs, and endothelial stress. Treat a monocyte-specific late signal as secondary. Plot the four-group donor values and compare signed pathway effects as well as gene effects.
4. Leave out each early control in turn and plot expression against delivery age within each stratum. Use overlap subsets descriptively and add a simple age term only if estimable. Repeat the priority estimates among donors without FGR, starting from 4 versus 3 early and 5 versus 6 late before profile filtering. Show captured fractions and broad-class subtype labels for context. [Organizer workflow](openbio-datathon-2026/track-01/3-workflow.md)

## Figures and follow-up

Show cohort coverage and delivery context, paired early and late effects with their direct difference, and the SCT donor values. Record at least three browser lookups as localization checks from the Admati study. If estimates align and the interval on their difference is broad, report a possible shared response with uncertain magnitude. A stable direct difference motivates independent spatial RNA follow-up with documented diagnosis week and overlapping delivery circumstances. Maternal sFLT1 and PlGF protein measures could connect a placental angiogenic result to a circulating signal. [Browser guide](openbio-datathon-2026/track-01/data/interpretation/zeisel_browser_guide.md)

Planning chat: “Preeclampsia Pseudobulk Strategy.” [OpenBio index](Thanh/README.md)

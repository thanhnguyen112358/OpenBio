# Q5 · Which placental signals accompany fetal growth restriction in early PE?

**Working hypothesis.** Fetal growth restriction (FGR) and early-onset preeclampsia (EOPE) may share parts of a villous stress response, with substantial variation across pregnancies. Prespecify Hallmark Hypoxia as the primary pathway, Hallmark Glycolysis as support, and inspect `FLT1` and `PGF` separately. Ackerman observed a metabolic pattern in a subset of the core cohort. GSE114691 provides external placental RNA data covering both conditions. [Ackerman et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC10192030/) · [Awamleh et al., GSE114691](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114691)

## Contrasts

| Role | Villous contrast | Donors | Question |
| --- | --- | --- | --- |
| Primary | `(EOPE+FGR) − EOPE` | 7 vs. 8 | Does FGR accompany a distinct response among pregnancies with EOPE? |
| Supporting | `FGR − iPTB` | 5 vs. 5 | Does a similar response appear outside PE? |
| Supporting | `EOPE − iPTB` | 8 vs. 5 | Does the program already appear in EOPE without FGR? |
| Context | `(EOPE+FGR) − iPTB` | 7 vs. 5 | Where does the combined diagnosis sit relative to preterm controls? |

Use the released `fgr_status` and `group` labels. The core provides birthweight intervals, so birthweight alone is an inadequate substitute for the clinical label. The primary comparison’s shared delivery-age range, 26.571–33.000 weeks, retains 4 of 7 and 7 of 8 donors. The `FGR − iPTB` overlap contains 2 FGR donors and 1 iPTB donor. Use that comparison mainly for directional evidence. [Villous metadata](openbio-datathon-2026/track-01/data/track01_metadata_vt.csv) · [Data dictionary](openbio-datathon-2026/track-01/data/track01_data_dictionary.md)

## Analysis

1. Check sample IDs, sequencing QC, and delivery-age overlap. Color PCA by group, age, and fetal sex. Keep source labels and examine whether flagged donors influence the primary result.
2. Fit the within-EOPE comparison on raw counts with `~ ga_delivery_weeks + group`. Show the unadjusted fit, the overlap-restricted fit, and leave-one-donor-out checks for the lead signal. Consider a separate fetal-sex sensitivity if estimable. Interpret the adjusted estimate as conditional on age at delivery.
3. Rank genes by signed Wald statistic for enrichment. Plot a fixed hypoxia score against age for each donor. Show `LDHA/LDHB` and `PDK1/GOT1` as published benchmarks and inspect their components. Place signed effects and intervals for the same genes and pathways across contrasts. If warranted, estimate `[(EOPE+FGR) − EOPE] − [FGR − iPTB]` as an exploratory interaction.
4. Look up at least three robust genes in the Zeisel browser and review markers for the implicated classes. Admati offers an exploratory within-early-PE comparison of 6 FGR versus 4 non-FGR donors before cell-class filtering. It lacks isolated-FGR controls. [Organizer workflow](openbio-datathon-2026/track-01/3-workflow.md) · [Pseudobulk donors](openbio-datathon-2026/track-01/data/extension/pseudobulk/track01_pseudobulk_donors.csv)

## Figures and follow-up

Show age overlap and PCA, donor-level pathway scores, and cross-contrast effects with cell context. A concordant within-EOPE and isolated-FGR shift supports a candidate shared FGR response. Similar effects in both EOPE groups would favor an EOPE-associated program, with the within-EOPE interval defining how much FGR effect remains plausible. For a hypoxia result, test a fixed panel against antenatal Doppler and blinded placental pathology in another cohort. For an angiogenic result, pair isoform-specific `FLT1` with maternal sFLT1 and PlGF protein measures.

Planning chat: “Placental Analysis Plan.” [OpenBio index](Thanh/README.md)

## November 14, 2026 morning review

Reviewed at 9:35 a.m. EST after the November 14 morning organizer repo update to [`fb751e9`](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/commit/fb751e9).

The released metadata confirms the original within-EOPE comparison of **7 EOPE+FGR vs. 8 EOPE** donors. Its shared delivery-age range retains **4 vs. 7**. For isolated FGR vs. iPTB the overlap retains only **2 vs. 1**, so that fit provides directional context with substantial uncertainty. The release provides clinical FGR labels and birthweight bins for the preterm groups. Preserve the clinical labels and log relevant QC flags for sensitivity checks. [Released villous metadata](openbio-datathon-2026/track-01/data/track01_metadata_vt.csv) · [Data dictionary](openbio-datathon-2026/track-01/data/track01_data_dictionary.md)

Use Q5 as the main plan if fetal growth is the team's chosen question. Lead with the within-EOPE hypoxia result across the full and overlap fits. For the updated article, export age overlap, PCA, donor-level hypoxia values, and signed effects across contrasts as four separate graphs. If [Q1](Q01-placental-gene-expression.md) leads, use Q5 as context for interpreting a possible FGR-related pattern. [Updated submission guide](openbio-datathon-2026/track-01/4-submit.md)

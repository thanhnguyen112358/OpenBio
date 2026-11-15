# Placental cell classes in the early-PE reference cohort

In the released placental pseudobulk, EOPE without FGR had lower C1Q and TREM2 scores than early controls in the broad macrophage class. Hofbauer cells showed a different pattern, with near-equal C1Q and higher TREM2. The macrophage differences retained their direction after restricting to the initial overlap of the groups' delivery-age ranges. These are descriptive differences for the EOPE-versus-early-control comparison. [Organizer data package](../../../../openbio-datathon-2026/track-01/data/README.md) · [Admati and Skarbianskis et al.](https://doi.org/10.1016/j.medj.2023.07.005)

| Cell class | Cohort | EOPE-only vs. early controls | C1Q | TREM2 | C1Q–TREM2 | Broad myeloid | C1Q–TREM2 minus broad myeloid |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Macrophage | All eligible | 3 vs. 3 | −0.43 | −0.61 | −0.52 | +0.06 | −0.57 |
| Macrophage | Initial delivery-age intersection | 2 vs. 2 | −0.13 | −0.23 | −0.18 | +0.14 | −0.32 |
| Hofbauer | All eligible | 4 vs. 3 | +0.003 | +1.35 | +0.67 | −0.14 | +0.81 |
| Hofbauer | Initial delivery-age intersection | 3 vs. 2 | +0.15 | +2.00 | +1.07 | +0.04 | +1.03 |

Each effect is the EOPE-only donor mean minus the early-control donor mean in log₂ score units. The inputs are the organizer's donor-by-class libraries with at least 10 captured cells. `PE_early_4` contributes to the Hofbauer comparison and had five captured macrophages, below the macrophage cutoff. Library size is the sum of all released raw counts for that donor and class. The analysis applies the Step 2 panel formulas to log₂(CPM + 1), giving `TREM2` and the three-gene C1Q component equal weight. [Step 2 scoring](../../step02_deseq2/README.md) · [Analysis script](../../../scripts/followup_pseudobulk_cell_context.py)

The groups' initial age ranges intersect at 29.857–34.143 weeks in both classes. Among the retained macrophage donors, the EOPE-only deliveries occurred at 29.857 and 31.857 weeks and the early-control deliveries at 32.857 and 34.143 weeks. All eligible EOPE-only donors delivered by caesarean and all eligible early controls vaginally. The age imbalance and separation by delivery mode limit attribution to disease. [Donor scores and metadata](donor_scores.csv) · [Exact contrasts and ages](descriptive_contrasts.csv)

This atlas has no isolated-FGR arm. Its broad macrophage and Hofbauer classes provide context for plausible cellular sources of the paired bulk RNA signal. A difference within a broad class can reflect the mix of subtypes as well as expression within cells. The analysis leaves trophoblast clearance and the original comparison of EOPE with isolated FGR for direct study. [Step 4 cell localization](../../step04_hypoxia_context/README.md)

Run from the workspace root in the Step 2 environment.

```powershell
& 'Deci_Vili_H1/.venv/Scripts/python.exe' Deci_Vili_H1/scripts/followup_pseudobulk_cell_context.py
```

# All 30 villous samples

Each point is one pregnancy’s placental villous bulk RNA-seq profile. Nearby points have similar expression across the two directions shown. PC1 and PC2 capture **20.9%** and **12.5%** of the variation in the `log2(CPM + 1)` expression matrix. All six plots show the same coordinates and change only the color. [PCA method](../../README.md)

## What the plots show

| View | Interpretation |
| --- | --- |
| [Clinical group](group.png) | Most idiopathic preterm birth (iPTB) samples and all term births sit toward the lower end of PC1. Early-onset preeclampsia (EOPE) without fetal growth restriction (FGR) spans a broader range, and all seven EOPE+FGR samples sit on the positive side of PC1 in this fit. S01, an EOPE+FGR sample, is farthest to the right. |
| [Delivery age](ga_delivery_weeks.png) | The five term births stand out by age, yet they share part of the PC1 range with iPTB samples delivered weeks earlier. Among preterm samples, the colors show a less consistent age gradient. |
| [Control type](control_type.png) | Preterm and term controls overlap in part of the display. The `not_control` label includes both EOPE groups and isolated FGR, so it covers several distinct clinical groups. |
| [FGR status](fgr_status.png) | FGR-labeled samples span much of PC1 because isolated FGR and EOPE+FGR occupy different parts of the plot. The clinical-group view makes that distinction visible. |
| [Fetal sex](fetal_sex.png) and [maternal age](maternal_age_years.png) | Both fetal sexes are spread across the plot, and maternal age has no clear gradient along either axis. These views describe the first two components. |

PC2 correlates with log gene-assigned read depth (absolute Pearson `r = 0.68`, calculated from the released HTSeq QC table). This makes read depth relevant when interpreting apparent vertical clusters and individual samples.

## What this contributes to the comparison

The all-sample view shows why the five clinical groups are worth keeping visible. Pooling into `PE` and `no_PE` would put isolated FGR, iPTB, and term births on the same reference side. The metadata give the stronger reason to use iPTB as the primary comparator. Term births occurred at **39.0–39.4 weeks**, beyond the delivery age of every EOPE case. The [preterm-only refit](../preterm/README.md) checks whether an EOPE–iPTB pattern is still visible once term samples are removed. [Q1 cohort plan](../../../Q01-placental-gene-expression.md) · [Q6 delivery-timing plan](../../../Q06-delivery-timing.md)

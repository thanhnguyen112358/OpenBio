# Preterm cohort · 25 pregnancies

This PCA was refit after removing the five term births. It includes 8 EOPE without FGR, 7 EOPE+FGR, 5 isolated FGR, and 5 idiopathic preterm birth (iPTB) pregnancies. Each point is a villous RNA profile. PC1 captures **19.8%** of the variation among these 25 samples and PC2 captures **13.6%**. The two plots share coordinates. The axes were fitted separately from the [all-villous PCA](../all_villous/README.md).

## What the plots show

- [Clinical group](group.png). Four of five iPTB samples sit toward the low end of PC1, while EOPE without FGR spreads across the middle and right. The iPTB sample S24 overlaps the EOPE range on PC1. Isolated FGR samples sit low on PC2, while EOPE+FGR spans a larger area and includes S01, the rightmost sample. EOPE with and without FGR partly overlap.
- [Delivery age](ga_delivery_weeks.png). PC1 has a weak linear association with gestational age at delivery (absolute Pearson `r = 0.20`). PC2 has a larger association with log gene-assigned read depth (`|r| = 0.60`) than with delivery age (`|r| = 0.39`). The age colors help inspect maturity within the preterm cohort, and the depth association adds context for vertical position.

## What this contributes to the comparison

The EOPE–iPTB shift is visible even with term births out of the fit. Isolated FGR also occupies a different region from iPTB in this view, which supports examining those pregnancies separately. The preterm cohort still spans **26.0–33.7 weeks**, and only **one of seven** EOPE+FGR pregnancies falls within the iPTB delivery-age window. The PCA shows a reason to look more closely at the focused contrast. The donor counts and delivery-age coverage determine which groups can support it. [Full EOPE–iPTB refit](../eope_vs_iptb/README.md) · [Q1 cohort plan](../../../Q01-placental-gene-expression.md)

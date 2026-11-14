# Shared delivery-age window · 5 EOPE without FGR vs. 5 iPTB

This PCA was refit on the ten pregnancies that delivered between **29 weeks 6 days and 33 weeks** (29.857–33.000 weeks in the released metadata). It includes five of the eight EOPE cases without FGR and all five iPTB controls. PC1 captures **29.3%** of expression variation in these ten samples and PC2 captures **19.2%**. The two plots show the same coordinates. The percentages describe this subset’s own PCA fit.

## What the plots show

- [Clinical group](group.png). Four iPTB controls remain at the low end of PC1, and most EOPE cases sit farther right. S24, an iPTB control, overlaps the EOPE PC1 range. The broad shift is still visible after excluding the three EOPE pregnancies that delivered before the iPTB range.
- [Delivery age](ga_delivery_weeks.png). The two groups share a delivery-age range. For example, S07 (EOPE) and S23 (iPTB) both delivered at 33 weeks and lie far apart on PC1. PC1 has a weak linear association with age across the ten pregnancies (absolute Pearson `r = 0.10`). PC2 has a strong association with log gene-assigned read depth (`|r| = 0.86`), and the lowest-depth sample, S22, lies at an extreme of that axis.

## What this contributes to the comparison

The first two components show that the broad expression shift remains visible among pregnancies with overlapping delivery ages. That supports keeping the **5-vs.-5 comparison as a sensitivity analysis** for the full 8-vs.-5 model. The window restriction also changes which cases are included, and this cohort still combines spontaneous preterm controls with EOPE cases whose delivery context can differ. To assess robustness, compare gene-level effect sizes and uncertainty across the full, age-adjusted, and shared-window fits. [Full comparison PCA](../eope_vs_iptb/README.md) · [Q1 analysis plan](../../../Q01-placental-gene-expression.md) · [Q6 delivery-timing plan](../../../Q06-delivery-timing.md)

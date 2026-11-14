# Primary comparison · 8 EOPE without FGR vs. 5 iPTB

This PCA uses the 13 pregnancies in the planned primary comparison. Each point represents one placental villous RNA profile. PC1 captures **25.1%** of the expression variation within these 13 samples and PC2 captures **17.1%**. Both plots use the same coordinates.

## What the plots show

- [Clinical group](group.png). Four of five iPTB controls sit at the low end of PC1. EOPE without FGR spans a wider range, and the fifth control, S24, falls within that range. The main pattern is a shift along PC1 with individual overlap.
- [Delivery age](ga_delivery_weeks.png). Both earlier and later deliveries are spread across the plot, and raw PC1 scores have a weak linear association with delivery age (absolute Pearson `r = 0.10`). Three EOPE pregnancies, S18, S19, and S20, still delivered earlier than the youngest iPTB control. PC2 correlates strongly with log gene-assigned read depth (`|r| = 0.87`). S22 has the lowest count of assigned reads in the comparison, about **5.6 million**, and lies at an extreme of PC2.

## What this contributes to the comparison

The 13-sample PCA shows that a broad expression difference is worth following up within the clinically chosen comparison. Its weak PC1–age correlation describes one axis and leaves age sensitivity of individual genes open. Use all 13 pregnancies for the primary estimate and the delivery-age-adjusted fit. Use the [5-vs.-5 shared-age refit](../shared_age/README.md) and corresponding gene-level model to check what happens when the three earliest cases are left out. The depth association also gives a concrete reason to inspect sequencing QC and donor influence. [Q1 analysis plan](../../../Q01-placental-gene-expression.md)

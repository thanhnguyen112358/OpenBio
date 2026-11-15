# H2 video narration

Read the paragraphs as speaker notes and use the time ranges as rehearsal targets. Coordinate speaker changes at slide boundaries so every team member appears or speaks. In speech, say the letters in BMP and BMPR2, and use “early preeclampsia” and “growth restriction” after the full terms have been introduced.

## 1. Clinical question · 0:00–0:26

How does placental RNA differ when early-onset preeclampsia accompanies fetal growth restriction? We looked at the chorionic villi, where exchange takes place between the mother and fetus. In OpenBio and Awamleh, we compared pregnancies with both conditions to pregnancies with growth restriction alone. Our measure was a bulk RNA score for 27 genes associated with bone morphogenetic protein, or BMP, signaling.

## 2. Delivery-age overlap · 0:26–0:56

We analyzed placental RNA collected at delivery. In OpenBio, each dot shows one pregnancy. The combined group usually delivered later. The shaded area shows the weeks covered by both groups, roughly 27 to 30. Only two of the seven combined cases delivered in that range, compared with all five isolated growth-restriction cases. Even inside the shaded range, those two combined cases delivered earlier than most of the isolated cases.

## 3. Pregnancy-level result · 0:56–1:44

We averaged the normalized log expression of the 27 genes to give each pregnancy one score. In the upper panels, each dot is one pregnancy and the bars show group means. The combined group scored about 0.27 log-two units higher in OpenBio and 0.15 higher in Awamleh. In each cohort, the difference remained positive when we dropped any one pregnancy. Both cohorts informed our pilot. The difference remained near 0.29 after a linear adjustment for delivery week in OpenBio. The long middle interval below shows how uncertain the estimate is when we use only pregnancies in the shared range. Awamleh’s public metadata lack individual delivery ages.

## 4. Rank enrichment and activin sensitivity · 1:44–2:15

We also ranked genes by the strength and direction of their expression differences. BMP genes tended to cluster toward the higher-expression end in all three analyses. BMPR2, a BMP receptor, contributed in each. When we removed two score genes shared with activin signaling, the difference got smaller in both cohorts. The Awamleh interval then included zero. The pattern fits a broader response involving several related pathways.

## 5. Candidate source cells · 2:15–2:44

Next we used the Admati placental reference to ask where these transcripts appear. Darker dots show more RNA, and bigger dots show detection in more donors. BMPR2 and SMAD9 were strongest in cells lining blood vessels. The inhibitor GREM2 was strongest in syncytiotrophoblast, the outer layer in contact with maternal blood. The browser showed the same locations. The reference contains preeclampsia and control samples and lacks an isolated growth-restriction group.

## 6. Marker-RNA sensitivity · 2:44–3:24

Could different cell proportions contribute to the bulk result? We added marker-RNA scores for endothelium and syncytiotrophoblast. The estimated difference stayed positive at about 0.35 log-two units, with a wider interval that includes zero. Marker RNA reflects both how many cells are present and how much RNA each cell makes. BMPR2 still had a positive estimate after accounting for endothelial marker RNA. SMAD9 tracked those markers more closely. That makes BMPR2 the first target for a cell-level follow-up.

## 7. A test in independent pregnancies · 3:24–4:04

We predict higher BMPR2 RNA inside villous endothelial cells in a new set of pregnancies with overlapping delivery ages. We would count endothelial and syncytiotrophoblast cells. We expect the bulk score difference to persist after accounting for those fractions. If the counts account for the bulk difference, that weakens the idea that the score reflects changes within cells. We would also measure phosphorylated SMAD inside the nuclei to test whether signaling differs. Its direction is open. Together these measurements address cell abundance, RNA within cells, and signaling activity.

## 8. Team and sources · 4:04–4:15

We thank the study participants, the teams who shared these datasets, and the OpenBio organizers. Our repository links the sources and the code that recreates the figures.

## Optional diagnostics insert after slide 3

As a diagnostic, we viewed the villous samples by diagnosis, sex, and inferred sequencing round. In the full OpenBio comparison, 45 genes passed a genome-wide false-discovery threshold of five percent. The highlighted BMP genes are pathway members. The pregnancy score and the ranking let us examine how those genes behaved as a set.

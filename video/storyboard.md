# H2 slide storyboard

Build a 16:9 deck. Use the supplied PNGs at their native aspect ratio on a white background and fit each entire image within the frame. Put captions and callouts in a margin outside the plot. Where a figure uses most of the slide, use the figure’s own title and one brief supporting line. Keep axes, legends, group labels, and uncertainty visible. The cues support the narration and can be done with a pointer or a simple outline.

## 1. Clinical question · 0:00–0:26

**On screen.** “Placental RNA when early preeclampsia accompanies fetal growth restriction.” Below, show “EOPE + FGR” and “Isolated FGR” side by side with equal visual weight. Caption both groups “Chorionic villi at delivery.”

**Layout and cue.** Use a plain title card and group labels. Define early-onset preeclampsia (EOPE) and fetal growth restriction (FGR) as the labels appear. Explain that the study asks what accompanies EOPE within pregnancies that already have FGR.

## 2. Delivery-age overlap · 0:26–0:56

**On screen.** “Shared delivery range ≈ 26.7–30.4 weeks · 2 EOPE+FGR vs 5 isolated FGR.”

**Layout and cue.** Center [age_support.png](figures/age_support.png) nearly full width. Reveal the caption as the narration points to the shaded band and the two combined pregnancies inside it. Leave the group counts and the gestational-age axis visible. Each point is a pregnancy. Triangles indicate blood-pressure flags on the released diagnoses. The two groups still differ in delivery timing within the band.

## 3. Pregnancy-level result · 0:56–1:44

**On screen.** “Higher mean 27-gene RNA score in both cohorts.” If adding a callout, use “Awamleh supports the direction. Both cohorts informed the pilot.”

**Layout and cue.** Fit the complete [headline.png](figures/headline.png) using the available height. Start at the upper OpenBio and Awamleh panels, where dots are pregnancies and bars are group means. The panels have separate cohort scales. Move the pointer to the lower forest plot, call out **+0.273** and **+0.145**, then hold on the long shared-age interval crossing zero. The lower plot shows within-cohort group differences with approximate 95% Welch intervals. The linear-week estimate **+0.286** is a separate fit and belongs in the speaker notes or a margin callout.

## 4. Rank enrichment and activin sensitivity · 1:44–2:15

**On screen.** “Removing INHA and INHBA attenuated both score differences.” If space permits, add “OpenBio +0.273 → +0.171 · Awamleh +0.145 → +0.085.”

**Layout and cue.** Center [enrichment.png](figures/enrichment.png) nearly full width. Keep the three Python fit labels, both axes, and the rank legend. Trace the left side of the curves, then point to BMPR2. Pause on the activin-shared-gene callout. This score sensitivity is a separate calculation from the three enrichment curves. The 25 plotted Reactome members come from the common rank universe. The pregnancy-level score uses 27 genes.

## 5. Candidate source cells · 2:15–2:44

**On screen.** “BMPR2 and SMAD9 · endothelium” and “GREM2 · syncytiotrophoblast.” A footer can read “Admati placental reference · no isolated-FGR arm.”

**Layout and cue.** Fit [atlas.png](figures/atlas.png) at maximum usable height. Keep the gene columns, cell-class rows, and both legends visible. Point to BMPR2 and SMAD9 in the endothelial row, then GREM2 in the syncytiotrophoblast row. Dot size is the share of eligible donors with at least 10 reads, and color is median counts per million. Crosses show unavailable summaries. The atlas supplies transcript localization for the follow-up.

## 6. Marker-RNA sensitivity · 2:44–3:24

**On screen.** “Both marker scores · +0.355 (95% CI −0.087 to +0.796).” Add “Marker RNA is a proxy for cell abundance” if there is room.

**Layout and cue.** Fit [composition.png](figures/composition.png) within the slide. Use a pointer on “Unadjusted” and “Both marker scores,” then show where the second interval crosses zero. Keep the additional model rows and the zero reference visible. The perivascular sensitivity can remain in the graphic. Describe BMPR2 and SMAD9 in the narration or a margin note. The plotted rows are for the bulk 27-gene score. The marker-adjusted BMPR2 and SMAD9 estimates come from individual-gene analyses.

## 7. A test in independent pregnancies · 3:24–4:04

**On screen.** “Prediction · Higher BMPR2 RNA within villous endothelial cells.” Above the workflow, show “EOPE+FGR vs isolated FGR · overlapping delivery ages.”

**Layout and cue.** Build a large editable three-step row with simple text boxes. Use “Count endothelial and syncytiotrophoblast cells,” “Measure BMPR2 RNA within endothelium,” and “Measure nuclear pSMAD1/5/9.” End on “Signaling direction open.” Keep pSMAD2/3 as a secondary readout in the notes. If measured cell fractions account for the bulk score, that weakens the within-cell explanation for the score. The BMPR2 RNA and phosphorylation measurements test distinct predictions.

## 8. Team and sources · 4:04–4:15

**On screen.** “Harvard Undergraduate OpenBio Datathon 2026 · Track 1,” `[team name and all member names]`, “OpenBio/Ackerman · Awamleh · Admati · Reactome/MSigDB,” and `github.com/thanhnguyen112358/OpenBio/tree/H2/video`. Put full data and code credits in the notes or description. Add the Piekos source if that comparison appears.

**Layout and cue.** Use a readable credits card with the names and repository link in large type. Hold it until the recording ends.

## Optional diagnostics and backup material

For a roughly 20-second diagnostics insert, show [pca.png](figures/pca.png) and [volcano.png](figures/volcano.png) after slide 3 as two consecutive figure cards. This moves the target end time to about 4:35. Fit the PCA to the width and retain the three legends for diagnosis, sex, and inferred run. Fit the volcano to the available height and retain the axes, teal BMP membership key, and effect direction. The vertical axis shows −log₁₀ of the nominal Wald p-value. If using a “45 genes” caption, label that number “Genome-wide FDR < 0.05 in the full OpenBio fit.” See the [optional narration](narration.md#optional-diagnostics-insert-after-slide-3) for these cards.

If someone asks about other cohorts or tissues, use the numbers in [evidence notes](evidence.md). Piekos is a severe-PE comparison. The paired decidual point estimate was larger, and the interval for the tissue difference included zero.

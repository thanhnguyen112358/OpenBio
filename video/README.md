# H2 video handoff

Harvard Undergraduate OpenBio Datathon 2026, Track 1. This package gives the slide editor and recording team the material for one presentation under five minutes. The proposed cut uses eight slides and runs about 4 minutes 15 seconds before rehearsal. Add the team byline and assign speakers before recording.

1. Start with the [storyboard](storyboard.md) for the slide order, on-screen copy, timings, and exact figure placement.
2. Use the [narration](narration.md) as speaker notes. Read it aloud with the slide cues and time the full recording.
3. Check numbers, methods, and attribution in the [evidence notes](evidence.md). Keep the supporting analyses available for questions or backup slides.
4. Give Claude the [paste-ready prompt](claude-prompt.md) and access to this folder, including the images.

## Figures in the recording

| Slide | File | Purpose |
| --- | --- | --- |
| 2 | [age_support.png](figures/age_support.png) | Show how many pregnancies fall within the shared delivery-age range. |
| 3 | [headline.png](figures/headline.png) | Show the full pregnancy dots and the three group-difference intervals. This is the same headline figure as in the article. |
| 4 | [enrichment.png](figures/enrichment.png) | Show the three Python rank fits and the marked BMPR2 position. |
| 5 | [atlas.png](figures/atlas.png) | Locate the candidate transcripts in the Admati reference. |
| 6 | [composition.png](figures/composition.png) | Show the positive point estimates and the wider marker-adjusted intervals. |
| Optional insert after 3 | [pca.png](figures/pca.png) and [volcano.png](figures/volcano.png) | Cover the organizer’s diagnostic plots and the OpenBio gene-level hit count. |

The bundled PNGs are exact copies of the [Python submission figures](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/figures) in the locked H2 package at `4343694`. The [submission article](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/docs/article.md) and [hypothesis card](../Pairwise_DESeq2/follow_up/bmp_workflow/python_submission/docs/hypothesis-card.md) provide the scientific basis. The package can be used without rerunning the analysis.

## Before submission

- Add the actual names to the final card and verify that every team member appears or speaks at least once.
- Replay the recording at its exported size. Check the plot labels, zero lines, confidence intervals, and audio. Confirm the runtime is under five minutes.
- Upload the team recording to YouTube as **Unlisted**, check the link signed out, and place it in the team submission README. Keep the same headline figure in the article and video. [Organizer Track 1 guide](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/fb751e91fb1eca4fb26a9658ad9d5042c382c4ba/track-01/4-submit.md)

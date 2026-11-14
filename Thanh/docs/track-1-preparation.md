# Track 1 preparation: preeclampsia and RNA sequencing

Researched November 10, 2026. This is an initial preparation checklist, not a completed analysis. Repository references are pinned to the version reviewed: `36257c75bf43a265174576c9251004ba8fdf2704`.

## What the challenge requires

The [current Track 1 overview](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/README.md) asks teams to compare placental gene expression in preeclampsia and control pregnancies, connect differences to cell types, and propose a testable explanation.

The [data guide](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/2-the-data.md) specifies:

- **Required core:** GSE203507, bulk RNA-seq from 30 placental samples, with counts and clinical metadata.
- **Required interpretation:** look up at least three genes using the Zeisel/Admati browser; a marker table provides an offline resource.
- **Optional extension:** donor-by-cell-type pseudobulk counts from the Admati/Skarbianskis single-cell study.
- **Reference atlases:** context for interpretation, rather than disease control groups.

At inspection, the repository tree listed the core counts, metadata, dictionary, pseudobulk tables, and checksum file as zero-byte placeholders. The guide explicitly marks unreleased files as pending. Their presence in the repository does not mean the analysis data are ready.

## What to submit

The [submission guide](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/4-submit.md) requests:

- A notebook or script that runs from beginning to end.
- Two to four figures, as 300-DPI PNGs or PDFs.
- A completed analysis plan and hypothesis card.
- A one- to two-page interpretation citing the datasets used.

The [submission procedure](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/submissions/README.md) puts Track 1 work under `submissions/your-team-name/` and requests a pull request before the deadline. Recheck that procedure and the exact deadline at kickoff.

The [rubric](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/docs/rubric.md) weights the biological question at 20%, rigor and reproducibility at 25%, analytical quality at 25%, hypothesis quality at 20%, and communication at 10%.

## Preparation checklist

These are proposed preparation tasks; they are not additional organizer requirements.

**Updated learning priority, November 10:** Follow the [biology-first roadmap](biology-first-roadmap.md) and [question guide](questions-biological-context.md) before selecting a computational workflow. The user's background is cell biology/biochemistry, with placental biology new to them. Software and practice-analysis tasks below are secondary to the clinical and mechanistic preparation.

- [x] Confirm the user's track preference.
- [x] Locate the current challenge guide, submission requirements, and rubric.
- [x] Check whether the released analysis data are available.
- [x] Establish the desired scientific approach and starting biology background.
- [ ] Understand the clinical syndrome and distinguish disease initiation, placental signals, and maternal organ effects.
- [ ] Learn normal placental anatomy and maternal–fetal cell identities.
- [ ] Read the source studies and map the biological meaning of all ten starter questions.
- [ ] Establish RNA-seq experience, preferred language (Python or R), and available preparation time when analysis preparation begins.
- [ ] Confirm registration, eligibility if needed, and team arrangements.
- [ ] Read the overview, data guide, and glossary in the order below.
- [ ] Choose an analysis environment and check the starter notebook with available files. A full run depends on the data release.
- [ ] Practice reading a count matrix and sample metadata using a small separate teaching dataset.
- [ ] Practice explaining PCA, effect sizes, adjusted p-values, and confounding.
- [ ] Shortlist a focused question, then finalize the comparison after reviewing the released metadata.
- [ ] Prepare to explain one result, an alternative explanation, and an experiment that could challenge the hypothesis.

## Analysis priorities

The organizer [workflow](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/3-workflow.md) proceeds through data checks, sample exploration, comparison selection, differential expression with DESeq2/PyDESeq2, and cell-type interpretation.

Preparation should emphasize matching sample IDs, distinguishing counts from normalized expression, and justifying the comparison before inspecting its gene results. Check gestational age and control type early. Treat donors as independent replicates, document exclusions, and distinguish changes in cell proportions from changes within cells. Browser expression patterns support interpretation but do not establish signaling or causation.

Two useful question directions are delivery-timing effects and the relationship between preeclampsia and fetal growth restriction. The [starter questions](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/starter-questions.md) identify both as core-data options. No question or model has been selected yet.

## Reading and setup links

1. [Track overview](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/README.md)
2. [Data guide](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/2-the-data.md)
3. [Organizer glossary](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/glossary.md)
4. [Analysis workflow](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/3-workflow.md)
5. [Setup options: Colab, local Python, or R](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/1-setup.md)
6. [Dataset and tool citation list](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/36257c75bf43a265174576c9251004ba8fdf2704/track-01/reference/citations.md)

Software setup, practice data, and analysis execution remain future work.

## November 14, 2026 morning review

Reviewed at 9:35 a.m. EST after the November 14 morning organizer repo update to [`fb751e9`](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/commit/fb751e9).

The November 10 checklist recorded placeholders and the then-current submission format. Track 1 v1.0.0 landed in `e5d6c0f` on November 13 at 1:22 a.m. EST. By the November 14 morning review, `main` was at `fb751e9`, with new submission requirements and hypothesis citation guidance. The starter questions and workflow notebooks match the data-release commit. [Data release](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/commit/e5d6c0f) · [Organizer comparison](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/compare/e5d6c0f...fb751e9)

### Confirmed data and analysis requirements

- The required villous core covers 30 pregnancies. Use the released raw counts, clinical labels, decimal `ga_delivery_weeks`, HTSeq QC, and HGNC annotation. Retain missing values as missing, report flags, and record any sensitivity exclusions. All 21 listed checksum entries matched the local files on review. [Data package](../../openbio-datathon-2026/track-01/data/README.md) · [Data dictionary](../../openbio-datathon-2026/track-01/data/track01_data_dictionary.md)
- Record at least three Zeisel/Admati gene-browser lookups. Optional additions include 27 donor-paired decidua samples, 362 donor-by-class profiles across 26 pseudobulk donors, and healthy reference atlases. The early pseudobulk comparison starts at 10 cases vs. 3 controls and is exploratory. [Track 1 data guide](../../openbio-datathon-2026/track-01/2-the-data.md)

### Current submission checklist

Submit by **Sunday, November 15 at 6:00 p.m. EST** via Devpost or a pull request to `submissions/your-team-name/` in the organizer repo. Make the team repo accessible to judges for the Devpost route. The updated Track 1 templates include a team README, article, and presentation outline. [Track 1 submission guide](../../openbio-datathon-2026/track-01/4-submit.md) · [Submission instructions](../../openbio-datathon-2026/submissions/README.md)

1. Include a team `README.md` with the question, team, video link, and run instructions. Put the analysis under `notebooks/`. After copying the starter, check its `../data` assumption and test the documented setup from a clean checkout or fresh session.
2. Write `docs/article.md` for general readers using the HUDP guide. Include **at least four graphs** with titles, axis labels, and captions. Support claims with data, cite outside sources, and describe limitations.
3. Save `figures/headline.png` and use it in the article and video. It counts toward the four graphs.
4. Record a presentation of **five minutes or less** with every team member appearing or speaking. Upload it to YouTube as unlisted and link it in the README.
5. Include `docs/analysis-plan.md` and one `docs/hypothesis-card.md`. In the card, connect a specific result to a plausible mechanism, at least one competing explanation, and a concrete test. Cite relevant literature. [HUDP article guide](../../openbio-datathon-2026/docs/article-guide.md) · [Judging rubric](../../openbio-datathon-2026/docs/rubric.md)

### Proposed scope for this submission

Start with [Q1](../../Q01-placental-gene-expression.md), EOPE without FGR vs. iPTB in villi (**8 vs. 5**). Use [Q6](../../Q06-delivery-timing.md) to show the full unadjusted fit, age-adjusted fit, and shared-age **5-vs.-5** subset. Add [Q9](../../Q09-baseline-gene-localization.md) browser lookups to explain the lead result. A focused article can use delivery-age overlap, PCA, effect estimates across fits, and donor-level expression as its four graphs. Pick the headline after inspecting the results. [Q5](../../Q05-preeclampsia-and-fetal-growth.md) is a focused alternative if the team wants fetal growth to drive the question. These are proposed analysis choices. This review verified the inputs and did not fit expression models.

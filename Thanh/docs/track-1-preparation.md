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

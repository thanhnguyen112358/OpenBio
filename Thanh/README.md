# Thanh’s OpenBio Datathon workspace

Start here for Thanh’s Track 01 preparation on placental RNA expression in preeclampsia. The ten Markdown gameplans in the parent folder synthesize the corresponding planning chats and follow the organizer’s question numbering. They contain working hypotheses, planned analyses, and proposed figures. The analysis plans use the organizer checkout at commit `e5d6c0f` as their reference for released data and metadata.

PE means preeclampsia, EOPE means early-onset preeclampsia, FGR means fetal growth restriction, and iPTB means idiopathic preterm birth. In the core release, `EOPE` is the group without FGR and `EOPE+FGR` is the group with both diagnoses.

## Gameplans

| Question | Gameplan | Main evidence |
| --- | --- | --- |
| Q1 | [Placental gene expression](../Q01-placental-gene-expression.md) | Villous RNA, EOPE versus preterm controls |
| Q2 | [Between-pregnancy heterogeneity](../Q02-between-pregnancy-heterogeneity.md) | Donor-level angiogenic and metabolic profiles |
| Q3 | [Cell-class differences](../Q03-cell-class-differences.md) | Donor pseudobulk, ranked within each onset comparison |
| Q4 | [Cell mix and within-class expression](../Q04-cell-mix-vs-expression.md) | Recovered fractions, pseudobulk, and bulk markers |
| Q5 | [Preeclampsia and fetal growth](../Q05-preeclampsia-and-fetal-growth.md) | Within-EOPE FGR contrast and isolated-FGR context |
| Q6 | [Delivery timing](../Q06-delivery-timing.md) | Term, preterm, and age-overlap comparisons |
| Q7 | [Early versus late onset](../Q07-early-vs-late-onset.md) | Direct comparison of two onset-specific effects |
| Q8 | [Cell communication](../Q08-cell-communication.md) | Donor-matched trophoblast and endothelial RNA |
| Q9 | [Baseline gene localization](../Q09-baseline-gene-localization.md) | Control pseudobulk, matched genes, and reference atlases |
| Q10 | [Healthy early-pregnancy atlas](../Q10-healthy-early-pregnancy-atlas.md) | Delivery-stage villi and paired decidua with atlas context |

## Preparation notes and papers

- [Preparation overview](preparation-overview.md) is the original project README and records Thanh’s initial biology-first focus.
- [Event brief](docs/event-brief.md) covers logistics and organizer sources.
- [Track 1 preparation](docs/track-1-preparation.md) summarizes the challenge, workflow, and submission materials.
- [Biology-first roadmap](docs/biology-first-roadmap.md) covers placental biology, disease mechanisms, and the reading sequence.
- [Biological context for the ten questions](docs/questions-biological-context.md) explains each prompt and the evidence that could answer it.
- [Primary literature](primary_lit/) contains six PDFs from Ackerman, Admati, Hiby, Maynard, Vento-Tormo, and Wang.

## Organizer checkout

The sibling `openbio-datathon-2026/` directory is the organizer’s Git checkout. Use its starter questions, data dictionary, workflow, and submission guide when implementing a gameplan. Keep the shared parent layout when copying this workspace so the relative links to the gameplans and checkout continue to resolve. [Starter questions](../openbio-datathon-2026/track-01/starter-questions.md) · [Data dictionary](../openbio-datathon-2026/track-01/data/track01_data_dictionary.md) · [Workflow](../openbio-datathon-2026/track-01/3-workflow.md) · [Submission guide](../openbio-datathon-2026/track-01/4-submit.md)

## November 14, 2026 morning review

Reviewed at 9:35 a.m. EST after the November 14 morning organizer repo update to [`fb751e9`](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/commit/fb751e9).

The root plans used `e5d6c0f` for the data release. The nested official checkout now matches `fb751e9`, which added the current submission guidance. Start with the [updated Track 1 checklist](docs/track-1-preparation.md#november-14-2026-morning-review) and [event brief update](docs/event-brief.md#november-14-2026-morning-review). Each root plan now ends with a review of what the release or new submission format means for that question. [Organizer comparison](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/compare/e5d6c0f...fb751e9)

The recommended first pass is [Q1's villous comparison](../Q01-placental-gene-expression.md), [Q6's timing sensitivity](../Q06-delivery-timing.md), and [Q9's browser lookups](../Q09-baseline-gene-localization.md). One targeted optional extension can follow if the primary finding gives us a reason to run it. The team still needs to select its final question.

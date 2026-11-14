# OpenBio Datathon 2026

This is our team workspace for Track 1 of the Harvard Undergraduate OpenBio Datathon. We’re exploring how placental gene expression and communication between maternal and fetal cells differ in pregnancies with preeclampsia.

## What’s in this repository

- `Q01-*.md` through `Q10-*.md` are working analysis plans, numbered to match the organizers’ starter questions.
- [Thanh’s workspace](Thanh/README.md) indexes the ten plans and links to the preparation notes and papers.
- [Primary literature](Thanh/primary_lit/) contains reference papers for the track.

## Get the team workspace

With Git installed and access to this private repository, run:

```sh
git clone https://github.com/thanhnguyen112358/OpenBio.git
cd OpenBio
```

## Download and update the organizer repository

The [official `openbio-datathon-2026` repository](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026) contains the challenge instructions, starter questions, and data setup guide. It’s a separate Git repository. Git ignores the local copy in our workspace, so cloning `OpenBio` doesn’t download it.

From inside `OpenBio/`, clone it once. This creates the `OpenBio/openbio-datathon-2026/` folder:

```sh
git clone https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026.git
```

To get the latest team notes and organizer materials, run these from `OpenBio/`:

```sh
git pull --ff-only
git -C openbio-datathon-2026 pull --ff-only
```

For data and software setup, follow the organizers’ [Track 1 setup guide](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/main/track-01/1-setup.md).

## November 14, 2026 morning review

Reviewed at 9:35 a.m. EST after the November 14 morning organizer repo update to [`fb751e9`](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/commit/fb751e9).

The November 14 morning checkout review found that the organizer had released the Track 1 data package and revised the submission guide since the initial event notes. The working recommendation is [Q1](Q01-placental-gene-expression.md) as the main analysis, [Q6](Q06-delivery-timing.md) for delivery-timing sensitivity, and [Q9](Q09-baseline-gene-localization.md) for cell-source interpretation. Each question plan has its own dated review. [Review and submission checklist](Thanh/docs/track-1-preparation.md#november-14-2026-morning-review)

Submit by **Sunday, November 15 at 6:00 p.m. EST** via Devpost or a pull request to the official repo. The package includes runnable code and a README, an article with at least four graphs, an unlisted video of up to five minutes, a headline visual reused in the article and video, the analysis plan, and a hypothesis card. [Official submission instructions](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/fb751e9/submissions/README.md)

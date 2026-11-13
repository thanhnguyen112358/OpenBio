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

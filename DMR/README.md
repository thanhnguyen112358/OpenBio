# DMR: Placental BMP RNA when early preeclampsia accompanies growth restriction

OpenBio 2026 Datathon, Track 1: Maternal-Fetal Cell Communication in Preeclampsia

## Team

| Name | Role |
|---|---|
| Duc Thanh Nguyen | Project owner |
| Minh Dau | Developers |
| Radley Le | Developers |

## Research Question

Among pregnancies with fetal growth restriction, is villous RNA for a BMP-centered gene set higher when early-onset preeclampsia is also present, after delivery timing and cell-marker RNA are taken into account?

## Data Used

- [x] GSE203507 bulk RNA-seq (required), through the OpenBio Track 1 release of Ackerman et al.
- [x] Zeisel/Admati gene browser, for BMPR2, GREM2, and SMAD9 (required). BMPR1B had no browser panel, so the released pseudobulk table was used for that gene.
- [x] Admati/Skarbianskis single-cell pseudobulk files (optional), to locate transcripts in cell classes.
- [x] Other public data: Awamleh et al., GSE114691; Ensembl GRCh37 release 75, to map Awamleh transcript IDs to genes; MSigDB Reactome 2024.1.Hs, for the BMP, activin, and TGF-β gene lists; Piekos et al. Hood-Lab placental multiomics, as a severe-preeclampsia context.

Organizer files are the Ackerman villous and decidua counts, the Admati donor-class tables, and the marker matrix, at commit `fb751e91fb1eca4fb26a9658ad9d5042c382c4ba`. Awamleh, Ensembl, Reactome, and Piekos are downloaded when the script runs. They are not stored in this folder.

## How to Run

1. The organizer data package is in [`track-01/data`](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/tree/fb751e91fb1eca4fb26a9658ad9d5042c382c4ba/track-01/data).
2. From this folder, with Python 3.12:

```sh
python -m pip install -r notebooks/requirements.txt
git clone https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026.git
git -C openbio-datathon-2026 checkout fb751e91fb1eca4fb26a9658ad9d5042c382c4ba
python notebooks/run.py --data openbio-datathon-2026/track-01/data --cache cache/bmp --output rerun --piekos
python notebooks/composition.py --data openbio-datathon-2026/track-01/data --scores rerun/bmp_donor_scores.csv
```

`run.py` needs network access to NCBI GEO, Ensembl, and the Broad Institute. `--piekos` also downloads the Piekos files from GitHub. The output directory must be empty. `composition.py` rebuilds the marker-RNA figure from the villous counts and the scores just written.

Reproducibility details:

- Notebook or script files: `notebooks/run.py` and `notebooks/composition.py`, with the modules they import.
- Where we ran it: our own computer.
- Random seeds: GSEApy adaptive multilevel sampling uses seed 1049. PyDESeq2 size factors are ratio estimates and do not use a seed.

## Files

| Path | Contents |
|---|---|
| `notebooks/` | The analysis scripts |
| `figures/` | Exported graphs, including `headline.png` |
| `docs/article.md` | Science communication article, following the [HUDP article guide](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/fb751e91fb1eca4fb26a9658ad9d5042c382c4ba/docs/article-guide.md) |
| `docs/analysis-plan.md` | The comparison, model, and thresholds we chose |
| `docs/hypothesis-card.md` | Our hypothesis and how to test it |

## Key Findings

In OpenBio, the 27-gene BMP score was 0.273 log₂ higher in early preeclampsia with growth restriction than in growth restriction alone (7 versus 5; interval 0.100 to 0.446). Awamleh showed the same direction (0.145; 20 versus 18). A straight-line adjustment for delivery week left the OpenBio difference at 0.286. Blood-vessel and syncytiotrophoblast marker RNA did not pull the point estimate toward zero. BMPR2 stayed higher after the blood-vessel markers were included. SMAD9 tracked those markers. The shared-week OpenBio comparison keeps only two combined pregnancies, and that interval is wide.

## Limitations

The primary comparison has 12 pregnancies, and the groups were not delivered at the same time. The age-matched slice has two combined pregnancies. The gene list shares members with activin signaling. Marker RNA is not a count of cells, and the decidual difference is at least as large as the villous one. These measurements are RNA at delivery. They do not show that a BMP signal fired inside a cell.

## Presentation

Video link: https://www.youtube.com/watch?v=YHjLy_rtnPQ

The recording must be five minutes or less. Every team member appears or speaks. Use `figures/headline.png` as the headline visual in the video and in the article.

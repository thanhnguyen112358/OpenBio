# DMR: Placental BMP RNA in early preeclampsia and fetal growth restriction

Harvard Undergraduate OpenBio Datathon 2026, Track 1 — Maternal-Fetal Cell Communication in Preeclampsia.

## Team

Duc-Thanh Nguyen, Minh Dau, Radley Le

## Research question

Among pregnancies with fetal growth restriction (FGR), is a 27-gene score for bone morphogenetic protein (BMP) pathway RNA in chorionic villi at delivery higher when early-onset preeclampsia (EOPE) also occurs, and how sensitive is the difference to delivery week and cell composition?

## Data used

- OpenBio Track 1 provides the restructured GSE203507 bulk RNA-seq tables for chorionic villi and decidua. [Organizer data package](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/tree/fb751e91fb1eca4fb26a9658ad9d5042c382c4ba/track-01/data) · [GEO GSE203507](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE203507) · [Ackerman et al.](https://doi.org/10.1161/HYPERTENSIONAHA.122.20807)
- The Zeisel lab browser and the organizer's Admati donor-by-cell-class pseudobulk and marker tables provide cell context. Browser observations for BMPR2, GREM2 and SMAD9, plus the unavailable BMPR1B panel, are bundled with the code. [Gene browser](https://storage.googleapis.com/www_zeisellab/PE_web/pe_web_test_v1.html) · [Browser observations](notebooks/browser_lookups.csv) · [Admati et al.](https://doi.org/10.1016/j.medj.2023.07.005) · [Single-cell data](https://doi.org/10.6084/m9.figshare.23264102.v1)
- Awamleh and colleagues released transcript-feature counts and diagnoses. Transcript IDs were aggregated with Ensembl GRCh37 release 75. [Awamleh et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6598374/) · [GEO GSE114691](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114691) · [Ensembl release 75](https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/)
- The pathway gene sets are from MSigDB Reactome 2024.1.Hs. [MSigDB release](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/) · [Reactome BMP signaling](https://reactome.org/content/detail/R-HSA-201451)
- The severe preeclampsia comparison uses the Piekos transcript and clinical releases. [Piekos et al.](https://www.nature.com/articles/s42003-025-08631-6) · [Data and code](https://github.com/spiekos/Hood-Lab-Placental-Multiomics)

## How to run

The archived run used Python 3.12.5 on a local computer. From this folder, use a Python 3.12 virtual environment and run the commands in order.

```sh
python -m pip install -r notebooks/requirements.txt
git clone https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026.git
git -C openbio-datathon-2026 checkout fb751e91fb1eca4fb26a9658ad9d5042c382c4ba
python notebooks/run.py --data openbio-datathon-2026/track-01/data --cache cache/bmp --output rerun --piekos
python notebooks/composition.py --data openbio-datathon-2026/track-01/data --scores rerun/bmp_donor_scores.csv
```

If you already have the organizer data, use its Track 1 data directory for both `--data` arguments. The `--output` directory must be new or empty. The workflow uses checksums for the downloaded GEO, Ensembl, Reactome and Piekos inputs. `run.py` writes the cohort tables and figures to `rerun/` and runs the built-in result checks. `composition.py` rebuilds the marker RNA estimates in `results/` and writes `figures/composition.png` using the same OpenBio normalization and BMP score. [Versions and input hashes from the archived run](run_manifest.tsv)

The score uses median-ratio normalization within each cohort. PyDESeq2 adjusts the OpenBio count analyses for fetal sex and fits the Awamleh comparison using group alone. GSEApy uses seed `1049` and adjusts for tests of three Reactome branches within each fit. To compare matching Python and R results, pass `--reference PATH_TO_R_RESULTS` to `run.py`. [Analysis plan](docs/analysis-plan.md)

## Files

| Path | Contents |
| --- | --- |
| [notebooks/](notebooks/) | Python scripts, dependencies, donor audit, and gene browser observations |
| [figures/](figures/) | Article figures and the headline visual |
| [results/](results/) | Estimates from the cell marker RNA sensitivity analysis |
| [docs/article.md](docs/article.md) | Science communication article |
| [docs/analysis-plan.md](docs/analysis-plan.md) | Comparisons, models, and thresholds |
| [docs/hypothesis-card.md](docs/hypothesis-card.md) | Interpretation and a testable follow-up |

## Key findings

The mean score was higher by 0.273 log₂ units in OpenBio (7 EOPE+FGR vs. 5 isolated FGR, 95% interval 0.100 to 0.446) and by 0.145 in Awamleh (20 vs. 18, 95% interval 0.037 to 0.253). In the full OpenBio count model, 45 genes had a genome-wide adjusted p-value below 0.05, including zero of the 25 BMP genes in the common ranking universe. The gene browser and donor-level reference showed the strongest BMPR2 and SMAD9 RNA signal in vascular endothelial cells and the strongest GREM2 RNA signal in syncytiotrophoblast. [Figures and interpretation](docs/article.md) · [Browser observations](notebooks/browser_lookups.csv)

## Limitations

Only two EOPE+FGR pregnancies fell in the shared delivery-week range. Their comparison with the five isolated-FGR pregnancies had a wide interval. Delivery conditions, growth restriction severity, and actual cell proportions could contribute to differences in bulk RNA. Marker RNA is an indirect measure of cell abundance, and pathway-associated RNA cannot determine signaling activity. [Hypothesis card](docs/hypothesis-card.md)

## Presentation

[Watch the DMR Track 1 presentation on YouTube](https://youtu.be/YHjLy_rtnPQ) · [Headline visual](figures/headline.png) · [Presentation script](docs/presentation-script.md)

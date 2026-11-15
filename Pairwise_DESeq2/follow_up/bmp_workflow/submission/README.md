# Placental BMP RNA in early preeclampsia and fetal growth restriction

Harvard Undergraduate OpenBio Datathon 2026, Track 1. The team name, member byline and unlisted presentation link are pending team handoff.

Does placental BMP-associated RNA differ when early-onset preeclampsia accompanies fetal growth restriction? See the [article](docs/article.md), [analysis plan](docs/analysis-plan.md), [hypothesis card](docs/hypothesis-card.md) and [presentation script](docs/presentation-script.md).

## Run the analysis

Install R 4.6.1, Bioconductor 3.23 and Python 3.12. The reference run used DESeq2 1.52.0, fgsea 1.38.0, BiocParallel 1.46.0 and PyDESeq2 0.5.4. The [run manifest](run_manifest.tsv) lists model settings, package versions and source SHA-256 checksums.

```r
install.packages(c("BiocManager", "digest", "ggplot2", "patchwork", "ragg", "matrixStats"))
BiocManager::install(c("DESeq2", "fgsea", "BiocParallel"), version = "3.23")
```

```sh
python -m pip install -r notebooks/requirements.txt
git clone https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026.git
git -C openbio-datathon-2026 checkout fb751e91fb1eca4fb26a9658ad9d5042c382c4ba
python notebooks/run.py --data openbio-datathon-2026/track-01/data --cache cache/bmp --output rerun --piekos
```

Run the commands from this folder. If the organizer repository is already available, point `--data` to its Track 1 data folder. Provide `--rscript` with the executable path if Rscript is outside PATH. The output directory must be empty. The script downloads the Awamleh and Reactome inputs and optionally downloads the Piekos inputs. The organizer count matrices and the included donor audit are fixed inputs. The command rebuilds the plots and checks the outputs. The cached intermediate files can remain outside the submission. The text notes from the browser were recorded by inspection and are included in `notebooks/browser_lookups.csv`.

The run writes donor membership, cohort support, score sensitivity estimates, the genome-wide count effects, the rank universe, enrichment, Admati donor-class summaries and figure source tables. `figures/headline.png` is the same visual used in the article and presentation.

## Data and credits

- OpenBio Track 1 provided the restructured Ackerman bulk data, the HTSeq summary, gene annotation, and Admati donor-class pseudobulk and marker releases. [Track 1 data](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/tree/fb751e91fb1eca4fb26a9658ad9d5042c382c4ba/track-01/data) · [Ackerman et al.](https://doi.org/10.1161/HYPERTENSIONAHA.122.20807) · [Admati et al.](https://doi.org/10.1016/j.medj.2023.07.005)
- Awamleh and colleagues released transcript-feature counts and diagnosis annotations. Transcript IDs were aggregated with Ensembl GRCh37 release 75. [Awamleh et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6598374/) · [GEO GSE114691](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114691) · [Ensembl release 75](https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/)
- The gene sets are from the MSigDB Reactome 2024.1.Hs collection. [MSigDB release](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/) · [Reactome BMP signaling](https://reactome.org/content/detail/R-HSA-201451)
- The optional severe-PE comparison uses Piekos and colleagues’ transcript and clinical releases. [Piekos et al.](https://www.nature.com/articles/s42003-025-08631-6) · [Data and code](https://github.com/spiekos/Hood-Lab-Placental-Multiomics)

Every team member must appear or speak in the final recording. Upload that recording as **Unlisted**, check that its link opens when signed out, and add the link to this README. [Organizer submission instructions](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/fb751e91fb1eca4fb26a9658ad9d5042c382c4ba/submissions/README.md)

# Placental BMP RNA in early preeclampsia and fetal growth restriction

Harvard Undergraduate OpenBio Datathon 2026, Track 1. The team byline and unlisted team presentation link are pending.

This Python workflow estimates the 27-gene BMP RNA score difference between pregnancies with early-onset preeclampsia and fetal growth restriction (EOPE+FGR) and pregnancies with isolated FGR. See the [article](docs/article.md), [analysis plan](docs/analysis-plan.md), [hypothesis card](docs/hypothesis-card.md) and [presentation script](docs/presentation-script.md).

## Reproduce the results

From this folder, with Python 3.12 in a virtual environment:

```sh
python -m pip install -r notebooks/requirements.txt
git clone https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026.git
git -C openbio-datathon-2026 checkout fb751e91fb1eca4fb26a9658ad9d5042c382c4ba
python notebooks/run.py --data openbio-datathon-2026/track-01/data --cache cache/bmp --output rerun --piekos
```

If the organizer data are already available, point `--data` to that Track 1 data directory. The output directory must be fresh or empty. The workflow retrieves checksum-pinned GEO, Ensembl 75, and MSigDB Reactome inputs. With `--piekos`, it also retrieves the pinned Piekos count and metadata releases. The audited organizer donor table and the manual browser observations are bundled in `notebooks/`. All figures are drawn from analysis tables. Paths in the archived run manifest are relative to the OpenBio repository root. [Versions and source hashes from the executed run](run_manifest.tsv)

The fixed pregnancy score uses the median of log count ratios for normalization. The count models use PyDESeq2, including the sex-adjusted OpenBio full and age-overlap fits and the group-only Awamleh fit. GSEApy runs adaptive multilevel enrichment with Benjamini–Hochberg adjustment across the three pathways per fit. The two implementations agree on the pregnancy scores and intervals and can give slightly different rank backgrounds, enrichment, and gene-level p-values. If you have the original R results, add `--reference PATH_TO_R_RESULTS` to compare the outputs. The Python workflow requires Python packages only.

The article and local presentation draft use the same Python-rendered headline figure. Every team member must appear or speak in the final recording. Upload the team recording as **Unlisted** and add the signed-out-accessible link here. [Organizer submission instructions](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/fb751e91fb1eca4fb26a9658ad9d5042c382c4ba/submissions/README.md)

On Windows, the synthetic review draft can be regenerated with the system speech voice and FFmpeg.

```sh
python -m pip install imageio-ffmpeg==0.8.1
python notebooks/render_video.py --script docs/presentation-script.md --figures figures --cache cache/video --output cache/presentation-draft.mp4
```

## Data and credits

- OpenBio Track 1 provided the restructured Ackerman bulk tables and Admati donor-class pseudobulk and marker tables. [Track 1 data](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/tree/fb751e91fb1eca4fb26a9658ad9d5042c382c4ba/track-01/data) · [Ackerman et al.](https://doi.org/10.1161/HYPERTENSIONAHA.122.20807) · [Admati et al.](https://doi.org/10.1016/j.medj.2023.07.005)
- Awamleh and colleagues released transcript-feature counts and diagnoses. Transcript IDs were aggregated with Ensembl GRCh37 release 75. [Awamleh et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC6598374/) · [GEO GSE114691](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE114691) · [Ensembl release 75](https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/)
- The gene sets are from MSigDB Reactome 2024.1.Hs. [MSigDB release](https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/) · [Reactome BMP signaling](https://reactome.org/content/detail/R-HSA-201451)
- The optional severe-PE analysis uses the Piekos transcript and clinical releases. [Piekos et al.](https://www.nature.com/articles/s42003-025-08631-6) · [Data and code](https://github.com/spiekos/Hood-Lab-Placental-Multiomics)

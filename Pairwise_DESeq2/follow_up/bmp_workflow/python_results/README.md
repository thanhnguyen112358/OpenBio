# Python analysis results

The [Python workflow](../python_submission/README.md) reruns the cohort audit, donor scores, PyDESeq2 contrasts, GSEApy enrichment, Admati donor-class summaries, paired and Piekos sensitivity analyses, and all six figures. The [original R report and tables](../README.md) remain available for comparison.

## BMP score at delivery

The difference is the mean score for EOPE+FGR minus the mean score for isolated FGR, in log₂ units.

| Fit | Pregnancies, combined vs isolated FGR | Difference (approximate 95% Welch CI) |
| --- | ---: | ---: |
| OpenBio | 7 vs 5 | +0.273 (+0.100 to +0.446) |
| OpenBio, shared delivery age | 2 vs 5 | +0.264 (−0.752 to +1.281) |
| Awamleh | 20 vs 18 | +0.145 (+0.037 to +0.253) |

The 109 pregnancy scores, the fixed score and sensitivity estimates, the optional paired and Piekos estimates, and the shared atlas summaries agree with the R run within numerical tolerance. Removing either or both flagged FGR donors preserved the OpenBio difference. Removing INHA and INHBA reduced the full-cohort estimates to +0.171 and +0.085. Delivery age remains poorly resolved in the shared-age subset. [Donor-level estimates](bmp_effects.csv) · [Piekos](piekos_effects.csv) · [Paired tissue](paired_tissue_effects.csv)

## Rank analysis

| Branch | OpenBio NES (FDR) | OpenBio shared-age NES (FDR) | Awamleh NES (FDR) |
| --- | ---: | ---: | ---: |
| BMP | +1.97 (0.0008) | +1.57 (0.0351) | +1.80 (0.0054) |
| Activin | +1.81 (0.0023) | +0.88 (0.2448) | +1.44 (0.0389) |
| TGF-β receptor complex | +1.59 (0.0023) | +1.00 (0.1997) | +1.49 (0.0188) |

PyDESeq2 admitted PPP2R2C to the common universe, for 13,726 genes and the same 25 BMP members. Both implementations identified the same seven shared BMP leading-edge genes. The Awamleh activin result is sensitive to the enrichment implementation at FDR 0.05. The BMP enrichment direction and interpretation hold across both runs. BMPR2 and SMAD9 RNA were highest in the endothelial samples in the Admati reference, and GREM2 RNA was highest in syncytiotrophoblast. Cell composition and signaling activity remain unresolved. [Python–R BMP comparison](python_vs_r_enrichment.csv) · [All enrichment fits](branch_enrichment.csv) · [Cell-class results](atlas_summary.csv)

The [Python article](../python_submission/docs/article.md) and [local video draft](presentation-draft.mp4) use the [headline figure](../python_submission/figures/headline.png) from this run. The draft uses synthetic narration and the team recording remains pending. [Run manifest](run_manifest.tsv)

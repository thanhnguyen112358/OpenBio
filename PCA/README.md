# Track 1 villous RNA-seq PCA

This script adapts Section 3 of the organizer’s Python starter notebook at [`fb751e9`](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/blob/fb751e9/track-01/notebooks/track01_workflow.ipynb). It reads the count matrix and metadata from the sibling `openbio-datathon-2026/track-01/data/` checkout and saves PNGs under `PCA/results/`. The input paths resolve from the script’s location.

From the OpenBio root, with the organizer checkout in place, run:

```sh
python -m pip install numpy pandas matplotlib scikit-learn
python PCA/organizer_pca.py
python PCA/organizer_pca.py --subsets
```

The baseline uses all released gene rows and the organizer’s `log2(CPM + 1)` transformation. It fits PCA once on all 30 samples and recolors the same coordinates by clinical group, delivery age, control type, FGR, fetal sex, and maternal age. The script seeds the PCA solver for reproducible reruns in the same environment.

`--subsets` adds three refits using the same method. The shared window is calculated from the minimum and maximum iPTB delivery ages. Each refit has its own axes.

| Results folder | Pregnancies | PC1 | PC2 |
| --- | ---: | ---: | ---: |
| [`all_villous`](results/all_villous/README.md) | 30 | 20.9% | 12.5% |
| [`preterm`](results/preterm/README.md) | 25 | 19.8% | 13.6% |
| [`eope_vs_iptb`](results/eope_vs_iptb/README.md) | 8 EOPE without FGR, 5 iPTB | 25.1% | 17.1% |
| [`shared_age`](results/shared_age/README.md) | 5 EOPE without FGR, 5 iPTB | 29.3% | 19.2% |

See the [primary comparison decision](../Q01-placental-gene-expression.md#primary-comparison-decision) for the cohort rationale, the PCA evidence, and the checks carried into the gene-level analysis.

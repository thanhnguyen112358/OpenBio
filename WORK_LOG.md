# Work log: Q1 placental villous RNA-seq (Track 1, Harvard OpenBio Datathon 2026)

**Written for:** teammates and their coding agents picking up this work.
**Last updated:** 2026-11-15, early morning (local time). Deadline is **Sunday 2026-11-15, 6:00 PM EST**.
**Status:** analysis in progress. Nothing below has been committed or pushed.

Read this file first. Then read, in order:
1. [`Q01-placental-gene-expression.md`](Q01-placental-gene-expression.md): research question and prespecified plan.
2. [`PCA/VST DESeq2/README.md`](PCA/VST%20DESeq2/README.md): teammate's original QC, PCA and DESeq2 fits A–D.
3. [`PCA/VST DESeq2/ADDITIONAL_CHECKS.md`](PCA/VST%20DESeq2/ADDITIONAL_CHECKS.md): sensitivity fits E and F.

---

## 1. Ground rules for anyone continuing

- **Do not rerun or overwrite fits A–D.** The user asked that the saved results in `PCA/VST DESeq2/results/deseq2/` stay as the baseline. New work reads those CSVs.
- **Keep prespecified and exploratory work separate.** Prespecified: fit A as primary, fits B–D as sensitivity checks, benchmark genes `FLT1`, `ENG`, `LEP`, `INHBA`, and three Hallmark sets (Hypoxia, Glycolysis, Oxidative Phosphorylation) ranked by the signed Wald statistic. Everything else is exploratory, including fits E–J, the other 47 Hallmark sets and the donor scores.
- **Gene selection rule:** `padj < 0.05` and `|log2FoldChange| >= 1`. This combines a test against zero effect with a cutoff on the estimated effect. It is not a formal test that the true change exceeds twofold.
- **Language for claims:** these are associations at delivery, relative to spontaneous preterm birth. Bulk RNA cannot separate a change in cell mix from a change within cells. Gene-level `FLT1` RNA does not measure the soluble sFLT1 isoform or circulating protein. Do not repeat the teammate README's line "delivery age does not explain the result." Write "results were relatively stable across the delivery-age sensitivity analyses."
- The user is learning RNA-seq and prefers explanations in Vietnamese, with concrete examples. Write repo documents in English.

## 2. Environment

| Item | Value |
|---|---|
| Workspace | `D:\datathon\OpenBio` |
| Remote / branch | `https://github.com/thanhnguyen112358/OpenBio.git`, branch `Nguyen` (tracks `origin/Nguyen`), HEAD `7d58af0` |
| Organizer repo | `openbio-datathon-2026/` at commit `fb751e9`. It is gitignored in this repo. |
| Python | `.\.venv\Scripts\python.exe`, run from the OpenBio root |
| Key packages | `pydeseq2 0.5.4`, `gseapy 1.1.4` (newly installed, **not** in the organizer `requirements.txt`) |
| Data checks | All 21 organizer checksums pass. Core counts: 26,485 genes × 30 samples. |

### Hallmark gene sets (manual step for a fresh checkout)

The organizer's `openbio-datathon-2026/track-01/data/get_gene_sets.py` downloads an **HTML login page** instead of the GMT file. Use the public Broad mirror for the same release:

```sh
curl -L -o openbio-datathon-2026/track-01/data/track01_gene_sets.gmt \
  "https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt"
```

Expect 50 lines (gene sets) and SHA-256 `ee2463540042078bfa3f67828e1e223bb354446d9fbb4d22845866835ba5c772`. `pathway_analysis.py` refuses to run on an HTML file.

## 3. Cohort and prior results (summary)

Primary comparison: **8 EOPE without FGR vs. 5 iPTB** (Q1). The shared delivery-age window is 29.857–33.000 weeks and contains 5 EOPE and 5 iPTB.

| Fit | Model / cohort | Genes passing | Source |
|---|---|---:|---|
| A (primary) | `~ group`, 8 vs. 5 | 352 | teammate, saved |
| B | centered delivery age + group, 8 vs. 5 | 383 | teammate, saved |
| C | `~ group`, shared-age 5 vs. 5 | 321 | teammate, saved |
| D | `~ group`, without S22, 7 vs. 5 | 326 | teammate, saved |
| E (exploratory) | `~ group`, without S24, 8 vs. 4 | 457 | `additional_checks.py` |
| F (exploratory) | `~ fetal_sex + group`, 8 vs. 5 | 294 | `additional_checks.py` |

- 207 genes pass in all of A–D. Of those, 197 still pass in E and 191 in F. All 352 fit A hits keep their direction in E and F.
- All four benchmarks pass in A–F and are higher in EOPE.
- Leaving out S24 has little effect on most genes (median absolute shift 0.129 log2), but `LEP` rises from 5.56 to 8.80 log2FC. S24 affects how large the LEP effect looks, not its direction. S22 and S24 both stay in the primary analysis.
- Fit F adjusts for fetal sex alone, not sex and delivery age jointly. Four genes did not converge in F (`C11orf86`, `KPRP`, `LOC389033`, `MIR93`); none passes.
- Outlier handling differs between designs (count replacement: A 69 genes, E 86, F 0). Comparisons between fits mix the intended change with filtering differences.
- Read depth is associated with diagnosis. Size factors cannot remove unmeasured library-quality differences that track diagnosis.

## 4. Work done in this session

### 4.1 Code changes

| File | Change |
|---|---|
| `PCA/VST DESeq2/deseq2_fits.py` | `run_fit()` now (a) keeps `fetal_sex` in its metadata and (b) takes an optional `contrast=(numerator, denominator)` defaulting to `("EOPE", "iPTB")`. A–D formulas and contrasts are unchanged. |
| `PCA/VST DESeq2/additional_checks.py` | **New.** Fits E and F, compared against the saved fit A. |
| `PCA/VST DESeq2/pathway_analysis.py` | **New.** Hallmark GSEA from saved fits, leading edges, donor scores, overlap and marker checks. |
| `PCA/VST DESeq2/context_contrasts.py` | **New.** Q1 step 4 context fits G–J plus GSEA. |
| `PCA/VST DESeq2/ADDITIONAL_CHECKS.md` | **New.** Report for E and F. |
| `.gitignore` (root) | Adds `.venv/` and `__pycache__/`. |
| `README.md` (root) | Local environment instructions. |
| 12 PNGs under `PCA/results/` | Regenerated during setup verification. They show as modified but their content should be equivalent. Decide whether to commit or `git checkout` them. |

### 4.2 Hallmark pathway analysis (`pathway_analysis.py`), Q1 step 3

Method, fixed before results were seen:
- **Ranking:** signed Wald `stat` of genes with non-missing `padj` (fit A: 15,453 genes). Source symbols are mapped to current HGNC symbols with `track01_gene_annotation.csv`, keeping only `approved`, `previous` and `alias` mappings. When two rows map to the same symbol, the row already using the approved symbol is kept, otherwise the one with higher `baseMean`.
- **Test:** `gseapy.prerank`, 10,000 gene permutations, seed `20261114`, weight 1, set size 15–500.
- **Fits ranked:** A (primary), plus B, C, D, E, F as sensitivity.

**Results so far** (NES > 0 means the set leans toward higher expression in EOPE; a reported p or q of 0 means below 1/10,000):

| Set | Fit A NES (FDR q) | Fit B | Fit C | Reading |
|---|---|---|---|---|
| Hypoxia | **+2.39** (< 0.001), top of all 50 sets | +2.52 (< 0.001) | +2.58 (< 0.001) | Strong and stable so far |
| Glycolysis | **+2.05** (< 0.001) | +2.06 (< 0.001) | +2.00 (< 0.001) | Stable, but shares 65 of 200 genes with Hypoxia |
| Oxidative phosphorylation | +1.28 (0.13) | +0.96 (0.68) | +0.59 (1.00) | **Not supported** |

- Other fit A sets with FDR q < 0.05 (exploratory): mTORC1 signaling, MYC targets V1, interferon-alpha response, protein secretion, E2F targets, unfolded protein response, inflammatory response, estrogen response late, and TNFα signaling via NF-κB. No set is significantly enriched toward lower expression in EOPE; the lowest q is Mitotic Spindle at 0.09.
- **Hypoxia leading edge (62 genes, 23 pass in fit A, 20 pass in all of A–D):** classic HIF targets, including `HK2`, `LDHA`, `BHLHE40`, `NDRG1`, `P4HA1`, `ERO1A`, `TMEM45A`, `STC1`, `AK4`, `PLIN2`, `MXI1`, `GAPDH`, `TPI1`, `PAM` and `CP`. Most also pass in E and F.
- The benchmarks are **not** in the Hypoxia set (`FLT1` is in no Hallmark set). The Hypoxia result is therefore not just the benchmark genes again. `INHA` and `PGF` are Hypoxia members.
- **Cell-mix check (organizer workflow step 7a), fit A:** trophoblast markers barely change (`KRT7` +0.06, `CDH1` −0.07, `TP63` −0.12, `ERVW-1` −0.12, `PSG1` −0.05, `CGB` +0.38, `CGA` +0.88 with padj 0.16). On the organizer's decision table, this favors a change within trophoblast over more trophoblast tissue, but the evidence is weak because there are few markers and the data are bulk.

### 4.3 Context comparisons (`context_contrasts.py`), Q1 step 4 (exploratory)

Each is `~ group` on the two groups only, with defaults and selection rule as in A–D. The groups differ in delivery age, and TB has no overlap with the preterm groups.

| Fit | Contrast | Question |
|---|---|---|
| G | EOPE+FGR vs. EOPE | Does adding growth restriction change the EOPE pattern? |
| H | FGR vs. iPTB | Is the pattern present with growth restriction alone? |
| I | EOPE vs. TB | How does a term-birth control change the answer? |
| J | EOPE+FGR vs. iPTB | Is the primary pattern also in EOPE with FGR? |

Status: running when this log was written. Fit G's table is saved; H–J and their GSEA were still pending.

## 5. Output locations

All outputs are under `PCA/VST DESeq2/results/`:

| Folder | Contents |
|---|---|
| `deseq2/` | Saved fits A–D (**do not overwrite**) |
| `additional_checks/` | Fits E, F, comparisons with A, Cook's summaries, benchmarks, `run_info.json` |
| `pathways/` | `hallmark_<fit>.csv` (all 50 sets per fit), `prespecified_across_fits.csv/.png`, `nes_all_fits.csv`, `fit_A_significant_sets_across_fits.csv`, `leading_edge_fit_A.csv`, `hypoxia_glycolysis_overlap_fit_A.csv`, `cell_class_markers_fit_A.csv`, `donor_scores_all_30.csv`, `angiogenic_vs_metabolic_within_EOPE.csv`, figures, `run_info.json` |
| `context_contrasts/` | Fits G–J, Hallmark tables, `summary.csv`, `benchmark_genes.csv`, `prespecified_pathways.csv`, `context_benchmarks_and_pathways.png`, `run_info.json` |

**Caveat about `pathways/`:** the first run started before two checks were added to the script: the Hypoxia/Glycolysis overlap split and the marker table. That run will not write `hypoxia_glycolysis_overlap_fit_A.csv` or `cell_class_markers_fit_A.csv`. **Rerun `pathway_analysis.py` once after it finishes** to get the complete, reproducible output set. Allow about 25 minutes, since each fit's GSEA takes about 4 minutes. With the fixed seed the NES values should reproduce.

### Commands (from the OpenBio root)

```powershell
.\.venv\Scripts\python.exe "PCA/VST DESeq2/additional_checks.py"   # E, F (~minutes)
.\.venv\Scripts\python.exe "PCA/VST DESeq2/pathway_analysis.py"    # ~25 min
.\.venv\Scripts\python.exe "PCA/VST DESeq2/context_contrasts.py"   # needs pathways/hallmark_A_full.csv
```

## 6. Remaining work (priority order for the deadline)

1. **Finish and verify the runs.** Check fits D–F GSEA, the overlap split (does Hypoxia hold without the 65 shared Glycolysis genes?), the donor scores and the context fits G–J. Look at every PNG for label collisions before using it.
2. **Do the angiogenic and metabolic patterns appear in the same donors?** Use `angiogenic_vs_metabolic_within_EOPE.csv` and its plot. With n = 8 this is descriptive only. Ackerman et al. reported a glycolytic pattern in a *subset* of EOPE placentas, so check for one.
3. **Three Zeisel browser lookups (required).** Proposed genes: `FLT1` (angiogenic benchmark), `NDRG1` or `HK2` (robust Hypoxia leading edge, high `baseMean`), and one EOPE-lower gene or a gene that tests the cell-mix explanation. For each, record the peak cell class, specificity, any PE vs. control difference in the browser, and notes, in `results/browser_lookups.csv` per the organizer workflow 7b. The browser is supporting context, not validation. **Ask the user before driving Chrome.**
4. **Write up.** Add the pathway and context sections to `PCA/VST DESeq2/README.md` or a new `PATHWAYS.md`, labeling exploratory work.
5. **Choose the lead finding and headline figure.** The current best candidate is a hypoxia/HIF-target signature (with FLT1) higher in EOPE than iPTB, stable across the delivery-age, S22, S24 and fetal-sex checks. OXPHOS is not supported. Save the figure to `figures/headline.png`. It must show every donor, use plain-language axes and be readable on video.
6. **Deliverables** (`openbio-datathon-2026/track-01/4-submit.md`): team README, reproducible code, an article with at least four graphs, the hypothesis card, and a video of five minutes or less. The hypothesis card needs a mechanism (for example, placental hypoxia/HIF activation driving `FLT1` and glycolytic genes), alternatives (cell mix, delivery and labor context, read depth, fetal sex), and a falsifying follow-up such as within-class pseudobulk in Admati, isoform-specific sFLT1, or an external cohort (GSE114691).
7. **Git.** Nothing is committed. Before committing, decide on the 12 regenerated PCA PNGs. Add `gseapy` to a project requirements note. Do not commit the organizer checkout or the GMT file (MSigDB license).

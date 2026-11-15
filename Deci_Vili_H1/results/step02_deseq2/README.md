# Step 2 DESeq2 normalization and panel scores

The analysis adapts `donor_expression()` and the gene-fitting pattern in Nguyen's `deseq2_fits.py` at commit `7826b8d`. Nguyen's pathway workflow also applies `log2(normalized count + 1)` before scoring. Here the inputs are the organizer's paired villous and decidua matrices. The donor groups and subsets come from the Step 1 audit. [Nguyen pipeline](https://github.com/thanhnguyen112358/OpenBio/blob/7826b8d4a146a8d888269b51fd7fbf7039de502f/PCA/VST%20DESeq2/deseq2_fits.py) · [Step 1 audit](../step01_cohort_audit/README.md)

## Normalized counts and donor scores

PyDESeq2 calculates one size factor per library using median ratios across genes. The factor estimates the relative scaling needed for that library. Dividing each observed count by its factor gives the normalized count. Factors were fitted separately for the 13 villous libraries and their 13 matched decidua libraries using the genome-wide raw count tables. The same factors are carried into every cohort subset. The factor calculation uses counts without phenotype terms. [DESeq2 documentation](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)

The score uses the full-precision library-level values `log2(normalized count + 1)`. The C1Q–TREM2 score gives half the weight to `TREM2` and half to the mean of `C1QA`, `C1QB`, and `C1QC`. The broad myeloid score averages `CD68`, `LST1`, and `CSF1R`. The donor table contains the two scores in each tissue and the within-donor differences specified in the [workflow](../../Deci_Vili_H1.md).

The descriptive EOPE-minus-FGR score means are shown below. A negative value in the last column means a larger decrease in decidua. The strict age subset reverses that contrast.

| Cohort | Villous difference | Decidua difference | Difference in the paired tissue score |
| --- | ---: | ---: | ---: |
| Inferred same round (4 vs. 5) | −0.964 | −1.215 | −0.251 |
| Same round within the full-cohort age overlap (3 vs. 5) | −1.024 | −1.403 | −0.379 |
| Same round within its own age overlap (3 vs. 4) | −1.138 | −1.056 | +0.081 |
| Full cohort (8 vs. 5) | −1.026 | −0.832 | +0.194 |

## Gene-level models

Each fit uses integer counts, the fixed size factors, and genome-wide dispersion estimation. The output reports the seven prespecified genes. The fitted log2 fold changes are group-level evidence alongside the donor-level scores.

| Model | Estimated effect |
| --- | --- |
| Villous `~ group` | EOPE minus FGR in villi |
| Decidua `~ group` | EOPE minus FGR in decidua |
| Paired `~ donor_id + is_decidua + eope_in_decidua` | EOPE minus FGR in the within-donor decidua-versus-villus effect |

All models use the four Step 1 cohorts. The inferred same-round subset is the primary comparison (4 EOPE vs. 5 FGR). The age sensitivities contain 3 vs. 5 and 3 vs. 4 donors. The full 8 vs. 5 set includes both inferred rounds and provides context. The paired model blocks on donor, and the interaction coefficient is negative when the EOPE–FGR fold change is more negative in decidua than in villi. A separate decidua fit shows the direction of the decidua effect itself. Wald intervals are approximate, and adjusted p-values are genome-wide within each fit after independent filtering. Count outlier replacement is disabled. [PyDESeq2 workflow](https://pydeseq2.readthedocs.io/en/v0.5.4/auto_examples/plot_step_by_step.html)

PyDESeq2 used its mean-based dispersion trend in the three paired fits restricted to the inferred same round. The other nine fits used the parametric trend. All 84 reported panel coefficients converged. The paired design gives each sample a distinct design row, which limits the replicate-based Cook's outlier filter. Donor influence is evaluated separately in Step 3.

## Outputs and reproduction

- [normalized_panel.csv](normalized_panel.csv) records raw, normalized, and log-transformed values by sample and gene.
- [donor_scores.csv](donor_scores.csv) records the tissue scores, size factors, within-donor differences, and cohort membership.
- [gene_effects.csv](gene_effects.csv) records the EOPE-versus-FGR fold changes, intervals, p-values, and fit diagnostics.
- [run_info.json](run_info.json) records software versions, formulas, and exact donor membership.

Run the [analysis script](../../scripts/step02_deseq2_scores.py) from the workspace root after the Step 1 audit. The local Windows environment and outputs stay under `Deci_Vili_H1/`.

```powershell
py -3.12 -m venv Deci_Vili_H1/.venv
& 'Deci_Vili_H1/.venv/Scripts/python.exe' -m pip install pydeseq2==0.5.4
& 'Deci_Vili_H1/.venv/Scripts/python.exe' Deci_Vili_H1/scripts/step02_deseq2_scores.py
```

# Hypothesis 1. Two signals, not one

At delivery, early-onset preeclampsia shows two separable placental RNA patterns.

- `FLT1` and Ackerman's ratios, `LDHA/LDHB` and `PDK1/GOT1`, are high when the mother has preeclampsia and are not clearly high when the only diagnosis is a small baby.
- A score of the ten glycolysis-enzyme steps is higher in preeclampsia and in isolated growth restriction, and the preeclampsia rise is present in both the villi and the decidua of the same pregnancies.

Adding growth restriction on top of preeclampsia changes none of these signals. The preeclampsia pattern is already there without it.

## How to rerun

From the `OpenBio` folder:

```bash
python hypothesis-1-two-signals/analysis.py
```

The script reads the Track 1 villous and decidua tables. Each pregnancy counts once. Scores are log2 gene activity. Group differences use a Welch test and a 95% confidence interval. This is not a DESeq2 genome-wide test. The gene list was fixed in `analysis.py` before the contrasts were written out.

## What the villous comparison showed

Difference versus spontaneous preterm birth (`iPTB`), log2 scale. Positive means higher in the named group.

| Signal | Preeclampsia, normal-size baby (8 vs 5) | Same, shared delivery weeks (5 vs 5) | Small baby, no preeclampsia (5 vs 5) |
|---|---:|---:|---:|
| `FLT1` | +2.66 (1.54 to 3.78) | +2.92 (1.41 to 4.42) | +0.83 (−0.45 to 2.11) |
| `LDHA/LDHB` | +1.80 (0.52 to 3.08) | +2.34 (0.45 to 4.23) | +0.17 (−0.65 to 0.99) |
| `PDK1/GOT1` | +1.37 (0.27 to 2.46) | +1.39 (−0.19 to 2.97) | +0.60 (−0.24 to 1.44) |
| Glycolysis enzymes | +0.47 (0.20 to 0.74) | +0.50 (0.09 to 0.90) | +0.39 (0.19 to 0.60) |

`PGF` is lower in preeclampsia, as the anti-angiogenic account predicts, but the villous interval includes zero (−0.52, −1.21 to 0.17). The decidua comparison is clearer (−0.81, −1.46 to −0.17).

Preeclampsia with and without a small baby do not differ: `FLT1` +0.05, enzyme score −0.03, both intervals wide and centered on zero.

## Both tissues

In the decidua of the same pregnancies, preeclampsia versus preterm birth gives `FLT1` +2.01 and an enzyme score +0.40. The shared-weeks check still shows both (`FLT1` +1.87, enzymes +0.48). The enzyme rise is not a villous-only pattern.

## Limits

Isolated growth restriction is delivered earlier than the preterm controls, and almost none of those pregnancies share an age window. The enzyme rise in that group is real in the numbers and still age-confounded. Do not treat it as a clean growth-restriction effect.

`FLT1` here is total gene RNA, not the soluble sFlt1 protein. The enzyme score is gene activity averaged over ten reaction steps, not a measurement of glucose converted to lactate. Groups are small, and one preeclampsia pregnancy sits with the preterm controls on every panel.

# Step 1 paired cohort audit

The audit links 30 villous samples and 27 decidua samples into **27 complete pregnancy pairs**. The three villous-only donors (`S01`, `S02`, and `S03`) all have EOPE+FGR. The clinical labels and delivery ages agree across both tissues for every complete pair.

## Inputs and checks

The script reads the two count matrices, two metadata tables, and HTSeq summary directly from the organizer's Track 1 checkout. That checkout was clean at commit `fb751e91fb1eca4fb26a9658ad9d5042c382c4ba`, with data package v1.0.0. All five SHA-256 checksums match the release manifest. The script checked unique donor and sample IDs, sample-column membership, the 26,485 shared genes, and nonnegative integer counts. The summed gene counts match `assigned_to_genes` for all 57 libraries. [Organizer release files](https://github.com/Harvard-Undergraduate-OpenBio/openbio-datathon-2026/tree/fb751e91fb1eca4fb26a9658ad9d5042c382c4ba/track-01/data)

## Cohort membership

| Cohort | EOPE without FGR | Isolated FGR | Delivery-age filter |
| --- | ---: | ---: | --- |
| All complete pairs in the two groups | 8 | 5 | Full ranges |
| Inferred same sequencing round | 4 | 5 | Full ranges |
| Inferred same round within the full-cohort age overlap | 3 | 5 | 26.714–30.429 weeks |
| Inferred same round within its own age overlap | 3 | 4 | 27.714–30.429 weeks |

Ackerman et al. report two sequencing rounds and name four EOPE-only donors (`S16`, `S20`, `S21`, and `S22`) sequenced alongside the preterm and term controls. The other donor assignments are inferred from that description. The Table S2 link was located, but the workbook could not be retrieved for sample-level confirmation. The cohort summary preserves the planned 3 vs. 5 subset and adds the stricter 3 vs. 4 subset. The stricter subset leaves out FGR donor `S15` because that pregnancy delivered before the earliest same-round EOPE-only donor. [Ackerman et al., Methods and supplementary materials](https://pmc.ncbi.nlm.nih.gov/articles/PMC10192030/)

## Read coverage and flags

- In the inferred same-round cohort, 8 of 9 villous libraries and all 9 decidua libraries have at least 10 reads for `TREM2`. All six other panel genes meet that count in all nine donors in both tissues. In the strict age-overlap subset, villous `TREM2` reaches 10 reads in 6 of 7 libraries. Treat the 10-read count as a coverage summary and carry those libraries into scoring.
- Both `S22_VT` and `S21_DB` stand out in the full two-phenotype cohort, with 5.59 million and 3.27 million assigned reads and assignment fractions of 22.1% and 15.8%, respectively. They remain in the full cohort and belong to the reported round with controls.
- The organizer flagged the recorded maximum blood pressures of `S04` and `S15`. Both retain their released FGR labels in the audit. Labor and delivery mode are unavailable in the released bulk metadata. [Organizer data dictionary](../../../openbio-datathon-2026/track-01/data/track01_data_dictionary.md)

## Files and reproduction

- [donor_audit.csv](donor_audit.csv) has one row per pregnancy, with source labels, ages, clinical flags, inferred round and its basis, sample IDs, tissue-specific QC and raw panel counts, and cohort membership. A missing decidua library has blank decidua fields. The assignment fraction is `assigned_to_genes` divided by the sum of that count and the five other released HTSeq counters.
- [cohort_summary.csv](cohort_summary.csv) records the four cohort definitions, group counts, age filters, and exact donor IDs.
- [step01_audit_cohort.py](../../scripts/step01_audit_cohort.py) rebuilds both tables with the Python standard library. From the workspace root, run:

  ```powershell
  python Deci_Vili_H1/scripts/step01_audit_cohort.py
  ```

These outputs establish the usable pairs, cohort membership, and library coverage ahead of scoring. C1Q–TREM2 and myeloid scores are calculated in Step 2.

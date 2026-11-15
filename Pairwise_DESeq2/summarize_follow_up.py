"""Extract the comparison panel from the full count model outputs."""

import argparse
import csv
import gzip
from pathlib import Path


GENES = set("ABCA7 AMD1 CDH5 EGFL7 FLT1 GALNT2 GPIHBP1 HTRA4 KDR LDHA LEPROTL1 LPL NTRK2 PDK1 PDE2A PECAM1 SCARB1 SLC2A1 TMEM45A VWF".split())
METRICS = ("baseMean", "log2FoldChange", "low95", "high95", "pvalue", "padj")
COMPARISONS = {("EOPE", "iPTB"), ("EOPE+FGR", "iPTB"), ("FGR", "iPTB"),
               ("EOPE", "FGR"), ("EOPE+FGR", "FGR")}
EXTERNAL_GROUPS = {"pe": "EOPE", "combined": "EOPE+FGR", "fgr": "FGR", "control": "preterm_control"}


def rows(path):
    with gzip.open(path, "rt", newline="") as source:
        yield from csv.DictReader(source, delimiter="\t")


def panel(root, external):
    for tissue in ("vt", "db"):
        for row in rows(root / "results" / f"{tissue}_contrasts.tsv.gz"):
            pair = row["numerator"], row["denominator"]
            if row["gene"] in GENES and pair in COMPARISONS:
                yield "OpenBio", tissue, "full", *pair, row["gene"], *(row[k] for k in METRICS)
        for fit in ("age_overlap", "shared_round", "round_and_age"):
            for row in rows(root / "follow_up" / "sensitivity" / f"{tissue}_{fit}.tsv.gz"):
                if row["gene"] in GENES:
                    yield "OpenBio", tissue, fit, "EOPE", "iPTB", row["gene"], *(row[k] for k in METRICS)
    paths = sorted(external.glob("gse114691_*_vs_*.tsv.gz"))
    if len(paths) != 6:
        raise ValueError("Expected six GSE114691 contrast tables")
    for path in paths:
        pair = tuple(EXTERNAL_GROUPS[g] for g in path.name.removeprefix("gse114691_").removesuffix(".tsv.gz").split("_vs_"))
        for row in rows(path):
            name = row["hgnc_symbol"]
            # Same HGNC:24945 locus; the Ensembl stable ID changed between builds.
            if row["ensembl_gene_id"] == "ENSG00000182851" and row["ensembl75_symbol"] == "GPIHBP1":
                name = "GPIHBP1"
            if name in GENES:
                yield "GSE114691", "vt", "full", *pair, name, *(row[k] for k in METRICS)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", type=Path)
    parser.add_argument("external", type=Path)
    args = parser.parse_args()
    result = sorted(panel(args.root, args.external))
    if len({row[:6] for row in result}) != len(result):
        raise ValueError("A cohort/tissue/fit/contrast/gene combination was duplicated")
    with open(args.root / "follow_up" / "candidate_effects.csv", "w", newline="", encoding="utf-8") as out:
        writer = csv.writer(out)
        writer.writerow(("cohort", "tissue", "fit", "numerator", "denominator", "gene", *METRICS))
        writer.writerows(result)

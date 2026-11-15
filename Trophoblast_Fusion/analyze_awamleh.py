"""Replicate the prespecified bulk contrast in the public GSE114691 release."""

import argparse
import csv
import gzip
import re
from hashlib import sha256
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from scipy.stats import ttest_ind


GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114691/"
SOURCES = {
    "PEONLY.txt.gz": (GEO + "suppl/GSE114691_MasterCount_PEONLY.txt.gz",
                      "db15986421b41d4e534f49438aa87223bfac91286f9317051b35096aea937a88"),
    "PEIUGR.txt.gz": (GEO + "suppl/GSE114691_MasterCount_PEIUGR.txt.gz",
                      "7df2e09365011d3757756620d8a8f7b8a28c261a2a6b8309311c464166acd3f3"),
    "series.txt.gz": (GEO + "matrix/GSE114691_series_matrix.txt.gz",
                      "2a7c799b027f0a5820d0acc469af60ac8f626fd4a1b18149ec44ccfd5b2dd7d9"),
    "ensembl75.gtf.gz": ("https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz",
                        "770f48d9aecd9bdecdf147fd141c23670b215bccbdbfa459048f45b02581c2ed"),
}
CASES = {
    "Preeclampsia only": ("PEONLY.txt.gz", "eoPE"),
    "Preeclampsia and Intrauterine Growth Restriction": ("PEIUGR.txt.gz", "eoPE_FGR"),
}
GENES = {"ENSG00000185201": "IFITM2", "ENSG00000142089": "IFITM3",
         "ENSG00000154451": "GBP5", "ENSG00000137270": "GCM1"}


def verify_sources(directory, download):
    directory.mkdir(parents=True, exist_ok=True)
    for name, (url, expected) in SOURCES.items():
        path = directory / name
        if not path.exists() and download:
            with urlopen(url, timeout=60) as stream:
                payload = stream.read()
            if sha256(payload).hexdigest() != expected:
                raise ValueError(f"Downloaded source checksum differs for {name}")
            path.write_bytes(payload)
        with path.open("rb") as stream:
            if sha256(stream.read()).hexdigest() != expected:
                raise ValueError(f"Source checksum differs for {name}")


def normalized_id(label):
    return re.sub(r"\W+", "", label.upper().removesuffix("_A")).removeprefix("X")


def transcript_counts(directory):
    with gzip.open(directory / "series.txt.gz", "rt") as stream:
        series = [row for row in csv.reader(stream, delimiter="\t") if row and row[0].startswith("!Sample_")]
    accessions = next(row[1:] for row in series if row[0] == "!Sample_geo_accession")
    descriptions = next(row[1:] for row in series if row[0] == "!Sample_description")
    clinical = next(row[1:] for row in series if row[0] == "!Sample_characteristics_ch1"
                    and row[1].startswith("disease group:"))
    records = {}
    for accession, description, field in zip(accessions, descriptions, clinical, strict=True):
        key = normalized_id(description)
        if key in records:
            raise ValueError(f"Ambiguous donor ID for {accession}")
        records[key] = accession, description, field.split(": ", 1)[1]
    frames, sample_rows = [], []
    for phenotype, (name, group) in CASES.items():
        table = pd.read_csv(directory / name, sep="\t", index_col=0)
        if not table.index.is_unique or np.any(table.to_numpy() < 0) or table.shape[1] != 20:
            raise ValueError(f"Invalid dimensions or counts for {name}")
        if not all(pd.api.types.is_integer_dtype(dtype) for dtype in table.dtypes):
            raise ValueError(f"Expected integer counts in {name}")
        for column in table.columns:
            accession, description, recorded_phenotype = records[normalized_id(column)]
            if recorded_phenotype != phenotype:
                raise ValueError(f"Mismatched phenotype for {accession}")
            sample_rows.append((accession, description, column, group))
        frames.append(table)
    if not frames[0].index.equals(frames[1].index):
        raise ValueError("Case group transcript IDs differ")
    meta = pd.DataFrame(sample_rows, columns=["accession", "description", "source_column", "group"])
    meta = meta.set_index("accession")
    if not meta.index.is_unique:
        raise ValueError("Duplicate GEO sample in the case count columns")
    counts = pd.concat(frames, axis=1).T.set_axis(meta.index)
    return counts, meta


def read_mapping(path, transcripts):
    mapping = {}
    with gzip.open(path, "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.split("\t", 8)
            if fields[2] != "transcript":
                continue
            transcript = re.search(r'transcript_id "([^".]+)', fields[8]).group(1)
            if transcript not in transcripts:
                continue
            gene = re.search(r'gene_id "([^".]+)', fields[8]).group(1)
            if mapping.setdefault(transcript, gene) != gene:
                raise ValueError(f"Transcript {transcript} maps to two genes")
    if len(mapping) != len(transcripts):
        raise ValueError(f"Mapped {len(mapping)} of {len(transcripts)} count features")
    return mapping


def sensitivity(log_counts, first):
    a, b = log_counts[first], log_counts[~first]
    test = ttest_ind(a, b, equal_var=False)
    ci = test.confidence_interval()
    return {"estimate": a.mean() - b.mean(), "low95": ci.low, "high95": ci.high, "pvalue": test.pvalue}


def main(directory, output, download):
    verify_sources(directory, download)
    counts, meta = transcript_counts(directory)
    mapping = read_mapping(directory / "ensembl75.gtf.gz", set(counts.columns))
    gene_counts = counts.T.groupby([mapping[transcript] for transcript in counts]).sum().T
    print(f"Matched {len(meta)} samples and mapped {len(mapping)} transcripts to {gene_counts.shape[1]} genes")
    gene_counts = gene_counts.loc[:, gene_counts.ge(10).sum(axis=0).ge(5)]
    model = DeseqDataSet(counts=gene_counts, metadata=meta, design="~ group", refit_cooks=False,
                         n_cpus=2, quiet=True)
    model.deseq2()
    result = DeseqStats(model, contrast=["group", "eoPE_FGR", "eoPE"], cooks_filter=False, n_cpus=2, quiet=True)
    result.summary()
    panel = result.results_df.loc[list(GENES)].rename_axis("ensembl_id").copy()
    panel.insert(0, "gene", [GENES[key] for key in panel.index])
    panel["low95"] = panel["log2FoldChange"] - 1.959963984540054 * panel["lfcSE"]
    panel["high95"] = panel["log2FoldChange"] + 1.959963984540054 * panel["lfcSE"]
    panel["lfc_converged"] = model.var.loc[panel.index, "_LFC_converged"]

    factors = model.obs["size_factors"]
    first = meta["group"].eq("eoPE_FGR")
    gene_log = np.log2(gene_counts.loc[:, list(GENES)].div(factors, axis=0) + 0.5)
    checks = [{"endpoint": gene, "analysis": "gene_welch", **sensitivity(gene_log[key], first)}
              for key, gene in GENES.items()]
    ifitm2 = gene_log["ENSG00000185201"]
    checks += [{"endpoint": "IFITM2", "analysis": f"omit_{sample}",
                **sensitivity(ifitm2.drop(sample), first.drop(sample))} for sample in meta.index]
    for transcript, gene in mapping.items():
        if gene == "ENSG00000185201":
            log_counts = np.log2(counts[transcript].div(factors) + 0.5)
            checks.append({"endpoint": transcript, "analysis": "IFITM2_transcript_welch",
                           **sensitivity(log_counts, first)})
    output.mkdir(parents=True, exist_ok=True)
    panel.to_csv(output / "awamleh_panel.csv")
    pd.DataFrame(checks).to_csv(output / "awamleh_sensitivity.csv", index=False)
    meta.to_csv(output / "awamleh_sample_join.csv")
    print(panel[["gene", "log2FoldChange", "low95", "high95", "lfc_converged"]].to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=Path(__file__).resolve().parent / "results")
    parser.add_argument("--download", action="store_true")
    args = parser.parse_args()
    main(args.input, args.output, args.download)

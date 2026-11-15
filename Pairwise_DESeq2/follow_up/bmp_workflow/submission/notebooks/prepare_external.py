"""Fit donor-level contrasts in the public Awamleh placental RNA-seq cohort."""

import argparse
import csv
import gzip
import io
import re
from concurrent.futures import ThreadPoolExecutor
from hashlib import sha256
from pathlib import Path
from urllib.request import urlopen

import numpy as np
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


GEO = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE114nnn/GSE114691/"
GTF = "https://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz"
COUNT_FILES = {
    "control": "ControlONLY",
    "pe": "PEONLY",
    "fgr": "IUGRONLY",
    "combined": "PEIUGR",
}
EXPECTED_DIAGNOSES = {
    "control": "Control",
    "pe": "Preeclampsia only",
    "fgr": "Intrauterine Growth Restriction only",
    "combined": "Preeclampsia and Intrauterine Growth Restriction",
}
SOURCE_HASHES = {
    "control": "3d429c4f8c4324cc5585f1a3dd85c36862b91eb6a6af1b89991f47b274ad020b",
    "pe": "db15986421b41d4e534f49438aa87223bfac91286f9317051b35096aea937a88",
    "fgr": "64f97cc1b88733f72777af4a0c778e06f822afd0a48a5f437d5ae0f492bed9c9",
    "combined": "7df2e09365011d3757756620d8a8f7b8a28c261a2a6b8309311c464166acd3f3",
    "series": "2a7c799b027f0a5820d0acc469af60ac8f626fd4a1b18149ec44ccfd5b2dd7d9",
    "gtf": "770f48d9aecd9bdecdf147fd141c23670b215bccbdbfa459048f45b02581c2ed",
}
CONTRASTS = [("pe", "control"), ("fgr", "control"), ("combined", "control"),
             ("pe", "fgr"), ("combined", "pe"), ("combined", "fgr")]


def download(item):
    role, url = item
    with urlopen(url, timeout=90) as response:
        payload = response.read()
    digest = sha256(payload).hexdigest()
    if digest != SOURCE_HASHES[role]:
        raise ValueError(f"Unexpected checksum for {url}: {digest}")
    print(f"{url.rsplit('/', 1)[-1]}: {len(payload):,} bytes, SHA256 {digest}", flush=True)
    return role, payload


def donor_key(value):
    return re.sub(r"\W+", "", value.upper().removesuffix("_A")).removeprefix("X")


def sample_records(payload):
    metadata = {}
    for row in csv.reader(io.StringIO(gzip.decompress(payload).decode()), delimiter="\t"):
        if row and row[0].startswith("!Sample_"):
            metadata.setdefault(row[0], []).append(row[1:])
    diagnoses = next(row for row in metadata["!Sample_characteristics_ch1"]
                     if row[0].startswith("disease group:"))
    records = {}
    for accession, description, diagnosis in zip(metadata["!Sample_geo_accession"][0],
                                                  metadata["!Sample_description"][0], diagnoses, strict=True):
        key = donor_key(description)
        if key in records:
            raise ValueError(f"Duplicate GEO identifier for {accession}")
        records[key] = accession, diagnosis.split(": ", 1)[1]
    return records


def transcript_counts(payloads, records):
    counts, sample_rows = [], []
    for group, short_name in COUNT_FILES.items():
        matrix = pd.read_csv(io.BytesIO(payloads[group]), sep="\t", compression="gzip", index_col=0)
        if not matrix.index.is_unique or not np.issubdtype(matrix.values.dtype, np.integer) or (matrix.values < 0).any():
            raise ValueError(f"Invalid transcript count table {short_name}")
        if counts and not matrix.index.equals(counts[0].index):
            raise ValueError("Transcript IDs differ between count tables")
        accession, diagnosis = zip(*(records[donor_key(column)] for column in matrix.columns), strict=True)
        if set(diagnosis) != {EXPECTED_DIAGNOSES[group]}:
            raise ValueError(f"Unexpected GEO disease label for {short_name}")
        counts.append(matrix.set_axis(accession, axis=1))
        sample_rows.extend(zip(accession, [group] * len(accession), diagnosis, matrix.columns, strict=True))
    metadata = pd.DataFrame(sample_rows, columns=["accession", "group", "geo_diagnosis", "source_column"])
    metadata = metadata.set_index("accession")
    if not metadata.index.is_unique:
        raise ValueError("Duplicate GEO sample in the joined count matrices")
    if metadata.groupby("group").size().to_dict() != {"control": 21, "pe": 20, "fgr": 18, "combined": 20}:
        raise ValueError("Unexpected clinical group counts")
    print(metadata.groupby(["group", "geo_diagnosis"]).size().to_string(), flush=True)
    return pd.concat(counts, axis=1), metadata


def annotation(payload, transcript_ids):
    transcript_to_gene, source_symbol = {}, {}
    for line in gzip.GzipFile(fileobj=io.BytesIO(payload)):
        fields = line.decode().split("\t", 8)
        if len(fields) != 9 or fields[2] != "transcript":
            continue
        attributes = fields[8]
        transcript = re.search(r'transcript_id "([^".]+)', attributes).group(1)
        if transcript not in transcript_ids:
            continue
        gene = re.search(r'gene_id "([^".]+)', attributes).group(1)
        gene_name = re.search(r'gene_name "([^"]+)', attributes)
        if transcript_to_gene.setdefault(transcript, gene) != gene:
            raise ValueError(f"Inconsistent gene assignment for {transcript}")
        source_symbol.setdefault(gene, gene_name.group(1) if gene_name else None)
    if len(transcript_to_gene) != len(transcript_ids):
        raise ValueError(f"Mapped {len(transcript_to_gene)} of {len(transcript_ids)} transcripts")
    return transcript_to_gene, pd.Series(source_symbol, name="ensembl75_symbol")


def current_symbol(path):
    table = pd.read_csv(path, usecols=["ensembl_gene_id", "hgnc_symbol"]).dropna().drop_duplicates()
    table = table.loc[~table["ensembl_gene_id"].duplicated(keep=False)]
    return table.set_index("ensembl_gene_id")["hgnc_symbol"]


def main(args):
    urls = {g: f"{GEO}suppl/GSE114691_MasterCount_{suffix}.txt.gz" for g, suffix in COUNT_FILES.items()}
    urls.update(series=f"{GEO}matrix/GSE114691_series_matrix.txt.gz", gtf=GTF)
    with ThreadPoolExecutor(max_workers=len(urls)) as executor:
        payloads = dict(executor.map(download, urls.items()))

    transcripts, metadata = transcript_counts(payloads, sample_records(payloads["series"]))
    groups, historical_symbols = annotation(payloads["gtf"], set(transcripts.index))
    genes = transcripts.groupby([groups[name] for name in transcripts.index]).sum().T
    eligible = genes.ge(10).sum(axis=0).ge(5)
    print(f"Mapped {len(transcripts):,} transcript rows to {genes.shape[1]:,} genes. "
          f"{eligible.sum():,} genes pass the count filter.", flush=True)
    symbol = current_symbol(args.annotation)

    args.output.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(args.output / "gse114691_matched_metadata.tsv", sep="\t")
    if args.gene_counts:
        genes.T.to_csv(args.output / "gse114691_gene_counts.tsv.gz", sep="\t", index_label="ensembl_gene_id")

    model = DeseqDataSet(counts=genes.loc[:, eligible], metadata=metadata, design="~ group",
                         refit_cooks=False, n_cpus=args.cpus, quiet=True)
    model.deseq2()
    for numerator, denominator in CONTRASTS:
        stats = DeseqStats(model, contrast=["group", numerator, denominator],
                           cooks_filter=False, n_cpus=args.cpus, quiet=True)
        stats.summary()
        effect = stats.results_df.rename_axis("ensembl_gene_id").join(historical_symbols).join(symbol)
        effect["low95"] = effect["log2FoldChange"] - 1.959963984540054 * effect["lfcSE"]
        effect["high95"] = effect["log2FoldChange"] + 1.959963984540054 * effect["lfcSE"]
        effect.to_csv(args.output / f"gse114691_{numerator}_vs_{denominator}.tsv.gz", sep="\t")
        print(f"Saved {len(effect):,} tests for {numerator} vs {denominator}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--annotation", type=Path, default=Path(
        "openbio-datathon-2026/track-01/data/track01_gene_annotation.csv"))
    parser.add_argument("--gene-counts", action="store_true")
    parser.add_argument("--cpus", type=int, default=2)
    main(parser.parse_args())

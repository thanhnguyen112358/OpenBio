"""Run the BMP cohort audit, models, scores, sensitivity analyses and figures in Python."""

import argparse
import csv
import hashlib
import importlib.metadata
import os
import platform
import shutil
from pathlib import Path
from types import SimpleNamespace
from urllib.request import urlopen

import pandas as pd
from threadpoolctl import threadpool_limits

from analyze import atlas_context, audit_primary, cohorts, donor_scores, indexed, read_pathways
from analyze import score_comparisons, score_manifest
from extensions import paired, piekos
from models import fit_and_rank
from prepare_external import COUNT_FILES, GEO, GTF, SOURCE_HASHES, main as prepare_awamleh
from report import render
from verify import verify


HERE = Path(__file__).resolve().parent
REACTOME = ("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/"
            "c2.cp.reactome.v2024.1.Hs.symbols.gmt")
REACTOME_SHA = "9ea1b5e656597daf423e41c5ebcaa9892bfedf3292fff768605d7b0d5e5e9703"
PIEKOS_BASE = ("https://raw.githubusercontent.com/spiekos/Hood-Lab-Placental-Multiomics/"
               "f57c4a74372a8b487c3c06a3d4ce8c4a7319323a/data/")
PIEKOS_SOURCES = {
    "normalized_cleaned/transcripts-final.csv.zip":
        "c6c7f85c396f18edc06ad31101edfd414ec51f02101e768e0123a40dcf4855f5",
    "raw/Final_Cohort.csv": "e2e79788586e80dac8501283432f77ec4aac8b95ef864b012bb01237bd9a7d25",
}
ORGANIZER_INPUTS = [
    "track01_counts_vt.csv", "track01_metadata_vt.csv", "track01_gene_annotation.csv",
    "qc/track01_htseq_summary.csv", "extension/decidua/track01_counts_db.csv",
    "extension/decidua/track01_metadata_db.csv", "extension/pseudobulk/track01_pseudobulk_counts.csv",
    "extension/pseudobulk/track01_pseudobulk_metadata.csv",
    "extension/pseudobulk/track01_pseudobulk_donors.csv", "extension/pseudobulk/track01_cell_class_map.csv",
    "interpretation/track01_marker_matrix.csv",
]


def sha256(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def fetch(url, checksum, path):
    if path.is_file() and sha256(path) == checksum:
        return path
    with urlopen(url, timeout=90) as response:
        content = response.read()
    if hashlib.sha256(content).hexdigest() != checksum:
        raise ValueError(f"Unexpected SHA-256 for {url}")
    path.write_bytes(content)
    return path


def portable_path(path):
    path = path.resolve()
    return Path(os.path.relpath(path)).as_posix() if path.drive == Path.cwd().drive else path.as_posix()


def manifest(out, files, sources, parameters):
    versions = {"Python": platform.python_version()}
    for package in ("pydeseq2", "numpy", "pandas", "scipy", "anndata", "scikit-learn",
                    "formulaic", "formulaic-contrasts", "gseapy", "matplotlib", "threadpoolctl"):
        versions[package] = importlib.metadata.version(package)
    with (out / "run_manifest.tsv").open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(("kind", "name", "value"))
        writer.writerows(("file_sha256", portable_path(path), sha256(path)) for path in files)
        writer.writerows(("source_sha256", url, value) for url, value in sources.items())
        writer.writerows(("setting", name, value) for name, value in {**versions, **parameters}.items())


def run(args):
    out, data, cache = args.output, args.data, args.cache
    if out.exists() and (not out.is_dir() or any(out.iterdir())):
        raise ValueError("--output must be a fresh or empty directory")
    cache.mkdir(parents=True, exist_ok=True)
    out.mkdir(parents=True, exist_ok=True)
    annotation_path = data / "track01_gene_annotation.csv"
    marker_path = data / "interpretation/track01_marker_matrix.csv"
    awamleh = cache / "awamleh"
    gmt = fetch(REACTOME, REACTOME_SHA, cache / "reactome.gmt")
    prepare_awamleh(SimpleNamespace(annotation=annotation_path, output=awamleh,
                                    gene_counts=True, cpus=args.cpus))
    urls = {name: f"{GEO}suppl/GSE114691_MasterCount_{suffix}.txt.gz" for name, suffix in COUNT_FILES.items()}
    urls.update(series=f"{GEO}matrix/GSE114691_series_matrix.txt.gz", gtf=GTF)
    sources = {url: SOURCE_HASHES[name] for name, url in urls.items()}
    sources[REACTOME] = REACTOME_SHA

    vt, db, ext_counts, villi, decidua, external = cohorts(data, args.audit, awamleh, out)
    overlap_ids, flagged = audit_primary(villi)
    annotation = indexed(annotation_path, "gene").reset_index(drop=True)
    pathways = read_pathways(gmt)
    genes = score_manifest(annotation, pathways, vt, ext_counts)
    fixed = genes[genes.in_score]
    local_scores = donor_scores(villi, vt, fixed, "local_gene")
    external_scores = donor_scores(external, ext_counts, fixed, "ensembl_gene_id")
    pd.concat([local_scores, external_scores], ignore_index=True).to_csv(out / "bmp_donor_scores.csv", index=False)
    score_comparisons(local_scores, external_scores, overlap_ids, flagged, out)
    paired(vt, db, villi, decidua, fixed, out)

    universe = fit_and_rank(vt, villi, overlap_ids, annotation, pathways, awamleh, external,
                            marker_path, out, args.cpus)
    genes["in_rank_universe"] = genes.source_symbol.isin(universe)
    genes.to_csv(out / "bmp_gene_manifest.csv", index=False)
    shared = pd.read_csv(out / "shared_leading_edge.csv").symbol
    n_donors, n_classes = atlas_context(data, marker_path, shared, out)
    browser = HERE / "browser_lookups.csv"
    shutil.copyfile(browser, out / browser.name)

    files = [*(data / path for path in ORGANIZER_INPUTS), args.audit, gmt, browser,
             awamleh / "gse114691_gene_counts.tsv.gz", awamleh / "gse114691_matched_metadata.tsv"]
    files += sorted(awamleh.glob("gse114691_*_vs_*.tsv.gz"))
    if args.piekos:
        piekos_inputs = []
        for name, checksum in PIEKOS_SOURCES.items():
            url = PIEKOS_BASE + name
            piekos_inputs.append(fetch(url, checksum, cache / Path(name).name))
            sources[url] = checksum
        piekos(*piekos_inputs, genes, out)
        files.extend(piekos_inputs)
    manifest(out, files, sources, {
        "normalization": "DESeq2 median of log ratios per full count matrix, separately by cohort/tissue",
        "score": "unweighted mean log2(normalized count + 1) over 27 fixed genes",
        "openbio_fit": "PyDESeq2 ~ fetal_sex + group; Wald test; Cook and independent filtering",
        "age_fit": "refit OpenBio count model within common support; score retains full size factors",
        "awamleh_fit": "PyDESeq2 ~ group; >=10 reads in >=5 libraries; Cook filtering/refit off",
        "enrichment": "GSEApy multilevel, seed 1049, threads 1, minSize 10, maxSize 500, eps 1e-10; BH within fit",
        "overlap_weeks": "26.714-30.429",
        "rank_universe_size": len(universe),
        "n_openbio_vt": len(villi), "n_openbio_db": len(decidua), "n_awamleh_vt": len(external),
        "n_atlas_donors": n_donors, "n_atlas_donor_classes": n_classes,
    })
    render(out)
    verify(out, args.reference)


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--data", type=Path, default=Path("openbio-datathon-2026/track-01/data"))
    cli.add_argument("--audit", type=Path, default=HERE / "donor_audit.csv")
    cli.add_argument("--cache", type=Path, default=Path("Pairwise_DESeq2/cache/bmp/python"))
    cli.add_argument("--output", required=True, type=Path)
    cli.add_argument("--cpus", type=int, default=2)
    cli.add_argument("--piekos", action="store_true")
    cli.add_argument("--reference", type=Path, help="Compare matching exports with the archived R run")
    options = cli.parse_args()
    if options.cpus < 1:
        cli.error("--cpus must be positive")
    with threadpool_limits(limits=1, user_api="blas"):
        run(options)

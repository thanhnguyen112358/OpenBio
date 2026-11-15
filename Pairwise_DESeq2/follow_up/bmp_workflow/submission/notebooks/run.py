"""Rebuild the OpenBio and Awamleh BMP analysis from the released counts."""

import argparse
import csv
import hashlib
import importlib.metadata
import os
from pathlib import Path
import platform
import shutil
import subprocess
import sys
from urllib.request import urlopen


HERE = Path(__file__).resolve().parent
REACTOME_URL = ("https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/"
                "c2.cp.reactome.v2024.1.Hs.symbols.gmt")
REACTOME_SHA = "9ea1b5e656597daf423e41c5ebcaa9892bfedf3292fff768605d7b0d5e5e9703"
DEFAULT_AUDIT = "Deci_Vili_H1/results/step01_cohort_audit/donor_audit.csv"
PIEKOS_BASE = ("https://raw.githubusercontent.com/spiekos/Hood-Lab-Placental-Multiomics/"
               "f57c4a74372a8b487c3c06a3d4ce8c4a7319323a/data/")
PIEKOS_INPUTS = {
    "normalized_cleaned/transcripts-final.csv.zip":
        "c6c7f85c396f18edc06ad31101edfd414ec51f02101e768e0123a40dcf4855f5",
    "raw/Final_Cohort.csv": "e2e79788586e80dac8501283432f77ec4aac8b95ef864b012bb01237bd9a7d25",
}


def run(args):
    subprocess.run([part.as_posix() if isinstance(part, Path) else str(part) for part in args], check=True)


def verified_download(url, sha256, path):
    if path.exists() and hashlib.sha256(path.read_bytes()).hexdigest() == sha256:
        return
    with urlopen(url, timeout=90) as response:
        content = response.read()
    if hashlib.sha256(content).hexdigest() != sha256:
        raise ValueError(f"Unexpected SHA-256 for {url}")
    path.write_bytes(content)


def write_provenance(awamleh):
    # The upstream preparation validates each of these source hashes before writing its outputs.
    sys.dont_write_bytecode = True
    from prepare_external import COUNT_FILES, GEO, GTF, SOURCE_HASHES

    sources = {role: f"{GEO}suppl/GSE114691_MasterCount_{suffix}.txt.gz"
               for role, suffix in COUNT_FILES.items()}
    sources.update(series=f"{GEO}matrix/GSE114691_series_matrix.txt.gz", gtf=GTF)
    rows = [("source_sha256", url, SOURCE_HASHES[role]) for role, url in sources.items()]
    rows.append(("source_sha256", REACTOME_URL, REACTOME_SHA))
    rows.append(("setting", "Python", platform.python_version()))
    for package in ("pydeseq2", "numpy", "pandas", "scipy", "anndata", "scikit-learn",
                    "formulaic", "formulaic-contrasts"):
        rows.append(("setting", package, importlib.metadata.version(package)))
    with (awamleh / "provenance.tsv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("kind", "name", "value"))
        writer.writerows(rows)


def run_piekos(options):
    files = []
    rows = []
    for name, checksum in PIEKOS_INPUTS.items():
        url = PIEKOS_BASE + name
        path = options.cache / Path(name).name
        verified_download(url, checksum, path)
        files.append(path)
        rows.append(("source_sha256", url, checksum))
    run([options.rscript, HERE / "piekos.R", *files, options.output])
    with (options.output / "run_manifest.tsv").open("a", newline="", encoding="utf-8") as handle:
        csv.writer(handle, delimiter="\t").writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("openbio-datathon-2026/track-01/data"))
    local_audit = HERE / "donor_audit.csv"
    if local_audit.drive == Path.cwd().drive:
        local_audit = Path(os.path.relpath(local_audit))
    parser.add_argument("--audit", type=Path, default=local_audit if local_audit.exists() else Path(DEFAULT_AUDIT))
    parser.add_argument("--cache", type=Path, default=Path("Pairwise_DESeq2/cache/bmp"))
    parser.add_argument("--output", type=Path, required=True, help="Fresh or empty results directory")
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--piekos", action="store_true", help="Include the separate severe-PE sensitivity")
    options = parser.parse_args()
    if options.output.exists() and any(options.output.iterdir()):
        parser.error("--output must be an empty directory")
    options.cache.mkdir(parents=True, exist_ok=True)
    options.output.mkdir(parents=True, exist_ok=True)
    reactome = options.cache / "reactome.gmt"
    openbio = options.cache / "openbio"
    awamleh = options.cache / "awamleh"
    upstream = HERE if (HERE / "run_pairwise.R").exists() else HERE.parent
    verified_download(REACTOME_URL, REACTOME_SHA, reactome)
    run([options.rscript, upstream / "run_pairwise.R", options.data, openbio])
    external = upstream / "prepare_external.py"
    run([sys.executable, "-u", external, "--annotation",
         options.data / "track01_gene_annotation.csv", "--output", awamleh, "--gene-counts"])
    sys.path.insert(0, str(upstream))
    write_provenance(awamleh)
    run([options.rscript, HERE / "analyze.R", options.data, options.audit, openbio,
         awamleh, reactome, options.output])
    shutil.copyfile(HERE / "browser_lookups.csv", options.output / "browser_lookups.csv")
    run([options.rscript, HERE / "report.R", options.output])
    if options.piekos:
        run_piekos(options)
    run([options.rscript, HERE / "verify.R", options.output])


if __name__ == "__main__":
    main()

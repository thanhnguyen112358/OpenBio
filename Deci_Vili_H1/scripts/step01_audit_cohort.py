"""Build the paired cohort audit from the organizer's Track 1 release."""

import csv
import hashlib
from decimal import Decimal
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
SOURCE = PROJECT.parent / "openbio-datathon-2026" / "track-01" / "data"
OUTPUT = PROJECT / "results" / "step01_cohort_audit"

VT_META = "track01_metadata_vt.csv"
VT_COUNTS = "track01_counts_vt.csv"
DB_META = "extension/decidua/track01_metadata_db.csv"
DB_COUNTS = "extension/decidua/track01_counts_db.csv"
HTSEQ = "qc/track01_htseq_summary.csv"

PANEL = ("TREM2", "C1QA", "C1QB", "C1QC", "CD68", "LST1", "CSF1R")
HTSEQ_FIELDS = (
    "assigned_to_genes",
    "no_feature",
    "ambiguous",
    "too_low_aQual",
    "not_aligned",
    "alignment_not_unique",
)
DONOR_FIELDS = (
    "group",
    "pe_status",
    "pe_onset",
    "fgr_status",
    "ga_delivery_source",
    "ga_delivery_weeks",
    "sbp_max_mmhg",
    "dbp_max_mmhg",
    "qc_flags",
)

# Ackerman et al., doi:10.1161/HYPERTENSIONAHA.122.20807, Methods.
CONTROL_GROUPS = {"iPTB", "TB"}
DONORS_SEQUENCED_WITH_CONTROLS = {"S16", "S20", "S21", "S22"}


def verify_source():
    manifest = SOURCE / "provenance" / "checksums.sha256"
    digests = {
        name: digest
        for digest, name in (
            line.split(maxsplit=1) for line in manifest.read_text().splitlines() if line
        )
    }
    for name in (VT_META, VT_COUNTS, DB_META, DB_COUNTS, HTSEQ):
        digest = hashlib.sha256((SOURCE / name).read_bytes()).hexdigest()
        if digest != digests[name]:
            raise ValueError(f"Organizer checksum mismatch: {name}")


def read_keyed_csv(name, key):
    rows = {}
    with (SOURCE / name).open(newline="", encoding="utf-8-sig") as handle:
        for row in csv.DictReader(handle):
            value = row[key]
            if not value or value in rows:
                raise ValueError(f"{name}: missing or repeated {key}: {value!r}")
            rows[value] = row
    return rows


def audit_counts(name, expected_samples, htseq):
    genes = set()
    panel = {}
    with (SOURCE / name).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.reader(handle)
        gene_column, *samples = next(reader)
        if gene_column != "gene" or len(samples) != len(set(samples)):
            raise ValueError(f"{name}: invalid count matrix header")
        if set(samples) != expected_samples:
            raise ValueError(f"{name}: count columns do not match tissue metadata")

        totals = [0] * len(samples)
        for gene, *counts in reader:
            if not gene or gene in genes or len(counts) != len(samples):
                raise ValueError(f"{name}: invalid count row for {gene!r}")
            values = [int(count) for count in counts]
            if min(values) < 0:
                raise ValueError(f"{name}: negative count for {gene}")
            genes.add(gene)
            totals = [total + count for total, count in zip(totals, values, strict=True)]
            if gene in PANEL:
                panel[gene] = dict(zip(samples, values, strict=True))

    if panel.keys() != set(PANEL):
        raise ValueError(f"{name}: missing panel genes {sorted(set(PANEL) - panel.keys())}")

    metrics = {}
    for sample, total in zip(samples, totals, strict=True):
        counters = {field: int(htseq[sample][field]) for field in HTSEQ_FIELDS}
        htseq_total = sum(counters.values())
        if min(counters.values()) < 0 or htseq_total <= 0:
            raise ValueError(f"{sample}: invalid HTSeq counters")
        if total != counters["assigned_to_genes"]:
            raise ValueError(f"{sample}: count total differs from HTSeq assigned reads")
        panel_counts = {f"{gene}_raw_count": panel[gene][sample] for gene in PANEL}
        metrics[sample] = {
            "assigned_to_genes": total,
            "htseq_total": htseq_total,
            "assignment_fraction": f"{total / htseq_total:.6f}",
            "panel_genes_ge10_reads": sum(count >= 10 for count in panel_counts.values()),
            **panel_counts,
        }
    return genes, metrics


def sample_ids(metadata, tissue):
    ids = set()
    for donor, row in metadata.items():
        sample = row["sample_id"]
        if sample != f"{donor}_{tissue}" or sample in ids:
            raise ValueError(f"Unexpected {tissue} sample ID for {donor}: {sample}")
        ids.add(sample)
    return ids


def age(row):
    return Decimal(row["ga_delivery_weeks"])


def age_overlap(rows):
    eope = [age(row) for row in rows if row["group"] == "EOPE"]
    fgr = [age(row) for row in rows if row["group"] == "FGR"]
    return max(min(eope), min(fgr)), min(max(eope), max(fgr))


def within_age(rows, window):
    lower, upper = window
    return [row for row in rows if lower <= age(row) <= upper]


def summarize_cohort(name, rows, window=None):
    eope = [row["donor_id"] for row in rows if row["group"] == "EOPE"]
    fgr = [row["donor_id"] for row in rows if row["group"] == "FGR"]
    return {
        "cohort": name,
        "eope_n": len(eope),
        "fgr_n": len(fgr),
        "age_filter_min_weeks": window[0] if window else "",
        "age_filter_max_weeks": window[1] if window else "",
        "eope_donor_ids": " ".join(eope),
        "fgr_donor_ids": " ".join(fgr),
    }


def write_csv(name, rows):
    with (OUTPUT / name).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    verify_source()
    vt_meta = read_keyed_csv(VT_META, "donor_id")
    db_meta = read_keyed_csv(DB_META, "donor_id")
    htseq = read_keyed_csv(HTSEQ, "sample_id")
    if db_meta.keys() - vt_meta.keys():
        raise ValueError("Decidua metadata contains donors absent from villous metadata")
    vt_ids = sample_ids(vt_meta, "VT")
    db_ids = sample_ids(db_meta, "DB")
    if set(htseq) != vt_ids | db_ids:
        raise ValueError("HTSeq sample IDs do not match released bulk metadata")

    vt_genes, vt_metrics = audit_counts(VT_COUNTS, vt_ids, htseq)
    db_genes, db_metrics = audit_counts(DB_COUNTS, db_ids, htseq)
    if vt_genes != db_genes:
        raise ValueError("Count matrices have different gene sets")

    rows = []
    empty_db = dict.fromkeys(next(iter(db_metrics.values())), "")
    for donor, vt in sorted(vt_meta.items()):
        db = db_meta.get(donor)
        if db:
            conflicts = [
                field
                for field in vt
                if field not in {"sample_id", "geo_accession"} and vt[field] != db[field]
            ]
            if conflicts:
                raise ValueError(f"{donor}: paired metadata disagrees for {conflicts}")
        with_controls = vt["group"] in CONTROL_GROUPS or donor in DONORS_SEQUENCED_WITH_CONTROLS
        vt_sample = vt["sample_id"]
        db_sample = db["sample_id"] if db else ""
        rows.append(
            {
                "donor_id": donor,
                **{field: vt[field] for field in DONOR_FIELDS},
                "complete_pair": int(bool(db)),
                "sequencing_round": "with_controls" if with_controls else "other_round",
                "sequencing_round_evidence": (
                    "reported_in_methods" if with_controls else "inferred_from_methods"
                ),
                "vt_sample_id": vt_sample,
                "db_sample_id": db_sample,
                **{f"vt_{key}": value for key, value in vt_metrics[vt_sample].items()},
                **{f"db_{key}": value for key, value in db_metrics.get(db_sample, empty_db).items()},
            }
        )

    full = [row for row in rows if row["complete_pair"] and row["group"] in {"EOPE", "FGR"}]
    same_round = [row for row in full if row["sequencing_round"] == "other_round"]
    planned_window = age_overlap(full)
    within_round_window = age_overlap(same_round)
    planned = within_age(same_round, planned_window)
    strict = within_age(same_round, within_round_window)

    cohorts = (
        ("all_eope_fgr_pairs", full, None),
        ("inferred_same_round", same_round, None),
        ("same_round_in_full_cohort_age_overlap", planned, planned_window),
        ("same_round_in_within_round_age_overlap", strict, within_round_window),
    )
    member_ids = {
        name: {row["donor_id"] for row in members} for name, members, _ in cohorts
    }
    for row in rows:
        for name, members in member_ids.items():
            row[f"in_{name}"] = int(row["donor_id"] in members)

    OUTPUT.mkdir(parents=True, exist_ok=True)
    write_csv("donor_audit.csv", rows)
    write_csv("cohort_summary.csv", [summarize_cohort(*cohort) for cohort in cohorts])
    print(f"Verified five organizer files and HTSeq totals for {len(htseq)} libraries.")
    print(f"Audited {len(rows)} pregnancies, {len(db_meta)} pairs, and {len(vt_genes):,} shared genes.")
    for name, members, window in cohorts:
        summary = summarize_cohort(name, members, window)
        print(f"{name}: {summary['eope_n']} EOPE vs {summary['fgr_n']} FGR")


if __name__ == "__main__":
    main()

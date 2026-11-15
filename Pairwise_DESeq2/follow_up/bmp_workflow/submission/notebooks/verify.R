args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 1L)
out <- args[[1L]]
expected <- c("run_manifest.tsv", "cohort_manifest.csv", "cohort_support.csv",
              "comparison_membership.csv", "bmp_gene_manifest.csv", "bmp_donor_scores.csv",
              "bmp_effects.csv", "bmp_gene_effects.tsv.gz", "branch_enrichment.csv",
              "rank_universe.csv", "rank_fit_membership.csv", "pca_scores.csv", "rank_plot.csv", "atlas_summary.csv",
              "atlas_donors.csv", "bulk_marker_effects.csv", "browser_lookups.csv",
              "paired_gene_coverage.csv", "paired_donor_scores.csv", "paired_tissue_effects.csv",
              file.path("figures", c("age_support.png", "enrichment.png", "headline.png",
                                    "atlas.png", "pca.png")))
missing <- expected[!file.exists(file.path(out, expected))]
if (length(missing)) stop("Missing workflow deliverables: ", paste(missing, collapse = ", "))
read_output <- function(filename) read.csv(file.path(out, filename))
scores <- read_output("bmp_donor_scores.csv")
genes <- read_output("bmp_gene_manifest.csv")
effects <- read_output("bmp_effects.csv")
membership <- read_output("comparison_membership.csv")
ranks <- read_output("rank_universe.csv")
enrich <- read_output("branch_enrichment.csv")
cohort <- read_output("cohort_manifest.csv")
stopifnot(!anyDuplicated(scores[c("cohort", "donor_id")]),
          sum(scores$cohort == "OpenBio") == 30, sum(scores$cohort == "Awamleh") == 79,
          !anyDuplicated(cohort[c("cohort", "sample_id")]),
          nrow(genes) == 31, sum(genes$in_score) == 27, sum(genes$in_rank_universe) == 25,
          setequal(genes$source_symbol[!genes$in_score], c("AMHR2", "BMP10", "CER1", "GDF2")),
          !anyDuplicated(ranks$symbol), nrow(ranks) == 13725, nrow(enrich) == 33,
          all(enrich$rank_size == nrow(ranks)))
target <- c(openbio_primary = 0.273, openbio_overlap = 0.264, awamleh_primary = 0.145,
            awamleh_eope_vs_fgr = 0.168, awamleh_fgr_vs_control = 0.146)
primary <- effects[effects$score == "BMP27", ]
stopifnot(all(round(primary$effect[match(names(target), primary$comparison)], 3) == target),
          sum(primary$comparison == "openbio_overlap") == 1L)
age <- membership[membership$comparison == "openbio_overlap", ]
stopifnot(setequal(age$donor_id, c("S04", "S08", "S11", "S12", "S14", "S15", "S17")),
          all(range(age$ga_delivery_weeks) == c(26.714, 30.429)),
          all(primary$ci_low <= primary$effect & primary$effect <= primary$ci_high))
for (i in seq_len(nrow(primary))) {
    row <- primary[i, ]
    members <- membership[membership$comparison == row$comparison, ]
    members$value <- scores$bmp27[match(paste(members$cohort, members$donor_id),
                                        paste(scores$cohort, scores$donor_id))]
    x <- members$value[members$arm == "numerator"]
    y <- members$value[members$arm == "denominator"]
    stopifnot(length(x) == row$n_a, length(y) == row$n_b,
              isTRUE(all.equal(mean(x) - mean(y), row$effect)))
    if (length(x) >= 2 && length(y) >= 2) {
        stopifnot(isTRUE(all.equal(as.numeric(t.test(x, y)$conf.int),
                                  c(row$ci_low, row$ci_high))))
    }
}
for (fit in unique(enrich$fit)) {
    result <- enrich[enrich$fit == fit, ]
    stopifnot(nrow(result) == 3, isTRUE(all.equal(p.adjust(result$pvalue, "BH"), result$fdr)))
}
paired <- read_output("paired_donor_scores.csv")
paired_effect <- read_output("paired_tissue_effects.csv")
stopifnot(nrow(paired) == 9L, !anyDuplicated(paired$donor_id),
          isTRUE(all.equal(paired_effect$effect[3], paired_effect$effect[2] - paired_effect$effect[1])))
if (file.exists(file.path(out, "piekos_effects.csv"))) {
    piekos <- read_output("piekos_effects.csv")
    adjusted <- subset(piekos, score == "BMP27" & comparison == "overlap_age_adjusted_exact_ID_only")
    stopifnot(nrow(adjusted) == 1L, adjusted$n_combined == 15L, adjusted$n_fgr == 16L,
              round(adjusted$effect, 3) == .141)
}
cat("Verified donor membership, Welch intervals, pilot score targets, gene coverage, rank universe, and figure inputs.\n")

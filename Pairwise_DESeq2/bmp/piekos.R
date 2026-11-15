#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) stop("Usage: piekos.R TRANSCRIPT_COUNTS_ZIP FINAL_COHORT_CSV OUTPUT_DIR")
out <- args[[3]]
gene_manifest <- read.csv(file.path(out, "bmp_gene_manifest.csv"))
source_genes <- gene_manifest$source_symbol
count_table <- read.csv(unz(args[[1]], "transcripts-final.csv"), check.names = FALSE)
stopifnot(!anyDuplicated(count_table[["Patient-ID"]]), !anyDuplicated(names(count_table)[names(count_table) %in% source_genes]),
          all(source_genes %in% names(count_table)))
count_matrix <- t(as.matrix(count_table[, !names(count_table) %in% c("", "Patient-ID")]))
stopifnot(all(is.finite(count_matrix) & count_matrix >= 0 & count_matrix == floor(count_matrix)))
factors <- estimateSizeFactorsForMatrix(count_matrix)
normalized_bmp <- sweep(count_matrix[source_genes, ], 2, factors, "/")

clinical <- read.csv(args[[2]], check.names = FALSE)[c("Study_ID/ID", "Condition", "WksGest")]
clinical <- clinical[!is.na(clinical[["Study_ID/ID"]]) & nzchar(clinical[["Study_ID/ID"]]), ]
stopifnot(!anyDuplicated(clinical[["Study_ID/ID"]]))
clinical$WksGest <- as.numeric(clinical$WksGest)
rna_id <- count_table[["Patient-ID"]]
stopifnot(identical(setdiff(rna_id, clinical[["Study_ID/ID"]]), "MJ-0323"))
# This published RNA ID differs from the unique clinical ID by a terminal replicate suffix.
# Export an exact-ID-only sensitivity so this linkage is visible and reversible.
clinical_id <- sub("^MJ-0323$", "MJ-0323-1", rna_id)
cohort <- clinical[match(clinical_id, clinical[["Study_ID/ID"]]), ]
cohort$rna_id <- rna_id
cohort$link <- ifelse(rna_id == clinical_id, "exact", "unique clinical ID with terminal -1")
cohort$BMP27 <- colMeans(log2(normalized_bmp[gene_manifest$in_score, ] + 1))
cohort$BMP31_pilot <- colMeans(log2(normalized_bmp + .5))
a <- "FGR+PE Severe"
b <- "FGR"
cohort <- cohort[cohort$Condition %in% c(a, b), ]
cohort$combined <- as.integer(cohort$Condition == a)
stopifnot(sum(cohort$combined) == 25L, sum(!cohort$combined) == 36L)
observed_support <- function(x) {
    ages <- lapply(split(x$WksGest, x$Condition), function(age) range(age, na.rm = TRUE))
    c(max(vapply(ages, `[`, numeric(1), 1L)), min(vapply(ages, `[`, numeric(1), 2L)))
}
support <- observed_support(cohort)
cohort$in_age_support <- !is.na(cohort$WksGest) & cohort$WksGest >= support[1] & cohort$WksGest <= support[2]
stopifnot(sum(cohort$combined[cohort$in_age_support]) == 15L,
          sum(!cohort$combined[cohort$in_age_support & cohort$link == "exact"]) == 16L)
hc3 <- function(y, x, adjust, name) {
    x$value <- x[[y]]
    fit <- lm(if (adjust) value ~ combined + WksGest else value ~ combined, x)
    design <- model.matrix(fit)
    bread <- solve(crossprod(design))
    weighted <- design * (residuals(fit) / (1 - hatvalues(fit)))
    covariance <- bread %*% crossprod(weighted) %*% bread
    estimate <- unname(coef(fit)["combined"])
    halfwidth <- qt(.975, df.residual(fit)) * sqrt(covariance["combined", "combined"])
    data.frame(score = y, comparison = name, n_combined = sum(x$combined),
               n_fgr = sum(!x$combined), effect = estimate, ci_low = estimate - halfwidth,
               ci_high = estimate + halfwidth)
}
results <- lapply(c("BMP27", "BMP31_pilot"), function(y) {
    overlap <- cohort[cohort$in_age_support, ]
    exact <- cohort[cohort$link == "exact", ]
    rbind(hc3(y, cohort, FALSE, "full_unadjusted"),
          hc3(y, overlap, FALSE, "overlap_unadjusted"),
          hc3(y, overlap, TRUE, "overlap_age_adjusted"),
          hc3(y, exact, FALSE, "full_unadjusted_exact_ID_only"),
          hc3(y, overlap[overlap$link == "exact", ], TRUE, "overlap_age_adjusted_exact_ID_only"))
})
write.csv(cohort, file.path(out, "piekos_donors.csv"), row.names = FALSE, na = "")
write.csv(do.call(rbind, results), file.path(out, "piekos_effects.csv"), row.names = FALSE)

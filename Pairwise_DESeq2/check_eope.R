#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L || !args[[3L]] %in% c("vt", "db")) {
    stop("Usage: Rscript check_eope.R TRACK01_DATA_DIR DONOR_AUDIT_CSV vt|db OUTPUT_DIR")
}
data_dir <- args[[1L]]
audit_path <- args[[2L]]
tissue <- args[[3L]]
out_dir <- args[[4L]]
if (tissue == "db") data_dir <- file.path(data_dir, "extension", "decidua")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

counts <- as.matrix(read.csv(file.path(data_dir, paste0("track01_counts_", tissue, ".csv")),
                             row.names = 1L, check.names = FALSE))
meta <- read.csv(file.path(data_dir, paste0("track01_metadata_", tissue, ".csv")),
                 row.names = "sample_id")
audit <- read.csv(audit_path, row.names = "donor_id")
stopifnot(setequal(rownames(meta), colnames(counts)),
          !anyDuplicated(meta$donor_id), all(meta$donor_id %in% rownames(audit)))
meta <- meta[colnames(counts), , drop = FALSE]
meta <- meta[meta$group %in% c("EOPE", "iPTB"), , drop = FALSE]
meta$group <- factor(meta$group, levels = c("iPTB", "EOPE"))
meta$round <- audit[meta$donor_id, "sequencing_round"]
stopifnot(!anyNA(meta[c("round", "ga_delivery_weeks", "fetal_sex")]),
          all(meta$round[meta$group == "iPTB"] == "with_controls"))

age_support <- function(frame) {
    ranges <- lapply(split(frame$ga_delivery_weeks, frame$group), range)
    c(max(vapply(ranges, `[[`, numeric(1L), 1L)), min(vapply(ranges, `[[`, numeric(1L), 2L)))
}

subset_by_age <- function(frame) {
    limits <- age_support(frame)
    frame[frame$ga_delivery_weeks >= limits[1L] & frame$ga_delivery_weeks <= limits[2L], , drop = FALSE]
}

round_data <- meta[meta$round == "with_controls", , drop = FALSE]
fits <- list(age_overlap = subset_by_age(meta), shared_round = round_data,
             round_and_age = subset_by_age(round_data))

cohorts <- lapply(names(fits), function(label) {
    cohort <- fits[[label]]
    cohort$fetal_sex <- factor(cohort$fetal_sex)
    design <- ~ fetal_sex + group
    matrix <- model.matrix(design, cohort)
    stopifnot(all(table(cohort$group) >= 2L), qr(matrix)$rank == ncol(matrix),
              nrow(matrix) > ncol(matrix))
    dds <- DESeqDataSetFromMatrix(counts[, rownames(cohort), drop = FALSE], cohort, design = design)
    dds <- DESeq(dds, quiet = TRUE)
    res <- as.data.frame(results(dds, contrast = c("group", "EOPE", "iPTB"), alpha = 0.05))
    interval <- qnorm(0.975) * res$lfcSE
    estimates <- data.frame(gene = rownames(res), res, low95 = res$log2FoldChange - interval,
                            high95 = res$log2FoldChange + interval, row.names = NULL)
    connection <- gzfile(file.path(out_dir, paste0(tissue, "_", label, ".tsv.gz")), "wt")
    on.exit(close(connection))
    write.table(estimates, connection, quote = FALSE, row.names = FALSE, sep = "\t")
    message(tissue, " ", label, ": ", sum(cohort$group == "EOPE"), " vs ", sum(cohort$group == "iPTB"))
    data.frame(fit = label, cohort[c("group", "round", "ga_delivery_weeks")],
               tissue = tissue, row.names = NULL)
})
write.csv(do.call(rbind, cohorts), file.path(out_dir, paste0(tissue, "_cohort_support.csv")), row.names = FALSE)

#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
    stop("Usage: Rscript Pairwise_DESeq2/run_pairwise.R TRACK01_DATA_DIR OUTPUT_DIR")
}
input_dir <- args[[1L]]
output_dir <- args[[2L]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

group_labels <- c("EOPE", "EOPE+FGR", "FGR", "iPTB", "TB")
group_levels <- c("EOPE", "EOPE_FGR", "FGR", "iPTB", "TB")
comparisons <- combn(seq_along(group_labels), 2L, simplify = FALSE)

write_tsv_gz <- function(table, path) {
    connection <- gzfile(path, "wt")
    on.exit(close(connection))
    write.table(table, connection, sep = "\t", quote = FALSE, row.names = FALSE)
}

run_tissue <- function(tissue, data_dir) {
    counts <- as.matrix(read.csv(
        file.path(data_dir, paste0("track01_counts_", tissue, ".csv")),
        row.names = 1L, check.names = FALSE
    ))
    meta <- read.csv(
        file.path(data_dir, paste0("track01_metadata_", tissue, ".csv")),
        row.names = "sample_id"
    )
    stopifnot(!anyDuplicated(colnames(counts)), !anyDuplicated(meta$donor_id),
              setequal(colnames(counts), rownames(meta)), all(counts >= 0),
              all(is.finite(counts)), all(counts == floor(counts)))
    meta <- meta[colnames(counts), , drop = FALSE]
    meta$group <- factor(meta$group, levels = group_labels, labels = group_levels)
    meta$fetal_sex <- factor(meta$fetal_sex)
    model <- ~ fetal_sex + group
    design_matrix <- model.matrix(model, meta)
    stopifnot(!anyNA(meta[c("group", "fetal_sex", "ga_delivery_weeks")]),
              all(table(meta$group) > 0L), qr(design_matrix)$rank == ncol(design_matrix))

    dds <- DESeqDataSetFromMatrix(counts, meta, design = model)
    dds <- DESeq(dds, quiet = TRUE)
    all_zero <- rowSums(counts) == 0L
    results_tables <- vector("list", length(comparisons))
    summary_tables <- vector("list", length(comparisons))

    for (k in seq_along(comparisons)) {
        pair <- comparisons[[k]]
        arm_a <- meta$group == group_levels[pair[1L]]
        arm_b <- meta$group == group_levels[pair[2L]]
        contrast <- c("group", group_levels[pair])
        estimates <- as.data.frame(results(dds, contrast = contrast, alpha = 0.05))
        both_zero <- rowSums(counts[, arm_a | arm_b, drop = FALSE]) == 0L
        status <- rep("adjusted_p_available", nrow(counts))
        status[is.na(estimates$padj)] <- "independent_filter"
        status[is.na(estimates$pvalue)] <- "p_unavailable"
        status[both_zero] <- "all_zero_in_both_groups"
        status[all_zero] <- "all_zero_in_tissue"

        # Neither group supplies evidence for an effect when both have zero reads.
        ci_eligible <- is.finite(estimates$lfcSE) & !both_zero
        half_width <- ifelse(ci_eligible, qnorm(0.975) * estimates$lfcSE, NA_real_)
        results_tables[[k]] <- data.frame(
            gene = rownames(estimates), numerator = group_labels[pair[1L]],
            denominator = group_labels[pair[2L]], estimates,
            low95 = estimates$log2FoldChange - half_width,
            high95 = estimates$log2FoldChange + half_width, status = status,
            row.names = NULL
        )

        age_a <- meta$ga_delivery_weeks[arm_a]
        age_b <- meta$ga_delivery_weeks[arm_b]
        age_lo <- max(min(age_a), min(age_b))
        age_hi <- min(max(age_a), max(age_b))
        significant <- !both_zero & !is.na(estimates$padj) & estimates$padj < 0.05
        summary_tables[[k]] <- data.frame(
            tissue = tissue, numerator = group_labels[pair[1L]],
            denominator = group_labels[pair[2L]], n_numerator = sum(arm_a),
            n_denominator = sum(arm_b),
            age_overlap_start = if (age_lo <= age_hi) age_lo else NA_real_,
            age_overlap_end = if (age_lo <= age_hi) age_hi else NA_real_,
            n_numerator_in_overlap = sum(age_a >= age_lo & age_a <= age_hi),
            n_denominator_in_overlap = sum(age_b >= age_lo & age_b <= age_hi),
            n_genes = nrow(counts), n_with_adjusted_p = sum(status == "adjusted_p_available"),
            higher_fdr05 = sum(significant & estimates$log2FoldChange > 0),
            lower_fdr05 = sum(significant & estimates$log2FoldChange < 0)
        )
    }

    table <- do.call(rbind, results_tables)
    stopifnot(nrow(table) == nrow(counts) * length(comparisons),
              !anyDuplicated(table[c("gene", "numerator", "denominator")]),
              all(table$low95 <= table$high95, na.rm = TRUE))
    write_tsv_gz(table, file.path(output_dir, paste0(tissue, "_contrasts.tsv.gz")))
    message(tissue, ": ", ncol(counts), " libraries, ", nrow(counts),
            " gene rows and ", length(comparisons), " contrasts")
    do.call(rbind, summary_tables)
}

summary <- rbind(run_tissue("vt", input_dir),
                 run_tissue("db", file.path(input_dir, "extension", "decidua")))
write.csv(summary, file.path(output_dir, "contrast_summary.csv"), row.names = FALSE)
capture.output(sessionInfo(), file = file.path(output_dir, "session_info.txt"))

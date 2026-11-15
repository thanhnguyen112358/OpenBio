#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) %in% 3:4) {
    stop("Usage: Rscript screen_pathways.R PAIRWISE_DIR ANNOTATION_CSV HALLMARK_GMT [EXTERNAL_TSV_GZ]")
}
root <- args[[1L]]
out <- file.path(root, "follow_up")
annotation <- read.csv(args[[2L]], na.strings = c("", "NA"))
annotation <- annotation[!is.na(annotation$hgnc_symbol), , drop = FALSE]
pathways <- fgsea::gmtPathways(args[[3L]])
stopifnot(length(pathways) == 50L, all(startsWith(names(pathways), "HALLMARK_")))

read_stats <- function(path) read.delim(path, check.names = FALSE)

make_ranks <- function(table, key = "gene") {
    annotation_key <- if (key == "ensembl_gene_id") "ensembl_gene_id" else "gene"
    if (key == "ensembl_gene_id") {
        table[[key]] <- sub("[.][0-9]+$", "", table[[key]])
    }
    table <- table[is.finite(table$stat) & is.finite(table$pvalue) & table$baseMean >= 10, , drop = FALSE]
    matches <- merge(table[c(key, "stat")], annotation[c(annotation_key, "hgnc_symbol")],
                     by.x = key, by.y = annotation_key)
    matches <- unique(matches)
    duplicate_symbol <- duplicated(matches$hgnc_symbol) | duplicated(matches$hgnc_symbol, fromLast = TRUE)
    matches <- matches[!duplicate_symbol, , drop = FALSE]
    setNames(matches$stat, matches$hgnc_symbol)
}

run_fgsea <- function(ranks, tissue, fit) {
    set.seed(1049)
    answer <- as.data.frame(fgsea::fgseaMultilevel(pathways, ranks, minSize = 15, maxSize = 500,
                                                  eps = 1e-10, BPPARAM = BiocParallel::SerialParam()))
    answer$leadingEdge <- vapply(answer$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1L))
    data.frame(tissue = tissue, fit = fit, n_ranked_genes = length(ranks), answer)
}

variants <- c("age_overlap", "shared_round", "round_and_age")
results <- lapply(c("vt", "db"), function(tissue) {
    full <- read_stats(file.path(root, "results", paste0(tissue, "_contrasts.tsv.gz")))
    full <- full[full$numerator == "EOPE" & full$denominator == "iPTB" &
                 !startsWith(full$status, "all_zero"), , drop = FALSE]
    ranks <- c(list(full = make_ranks(full)), setNames(lapply(variants, function(fit) {
        make_ranks(read_stats(file.path(out, "sensitivity", paste0(tissue, "_", fit, ".tsv.gz"))))
    }), variants))
    shared <- Reduce(intersect, lapply(ranks, names))
    stopifnot(length(shared) >= 5000L)
    do.call(rbind, lapply(names(ranks), function(fit) run_fgsea(ranks[[fit]][shared], tissue, fit)))
})

if (length(args) == 4L) {
    external <- read_stats(args[[4L]])
    stopifnot(all(c("ensembl_gene_id", "baseMean", "stat", "pvalue") %in% names(external)))
    ranks <- make_ranks(external, "ensembl_gene_id")
    stopifnot(length(ranks) >= 5000L)
    results[[length(results) + 1L]] <- run_fgsea(ranks, "vt", "awamleh")
}

results <- do.call(rbind, results)
stopifnot(!anyDuplicated(results[c("tissue", "fit", "pathway")]),
          all(is.finite(results$NES)), all(results$padj >= 0 & results$padj <= 1))
write.csv(results, file.path(out, "hallmark_eope_vs_preterm.csv"), row.names = FALSE)

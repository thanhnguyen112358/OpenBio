#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(DESeq2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
    stop("Usage: analyze.R TRACK01_DATA_DIR DONOR_AUDIT_CSV OPENBIO_DE_DIR AWAMLEH_DIR REACTOME_GMT OUTPUT_DIR")
}
data_dir <- args[[1]]
openbio_dir <- args[[3]]
awamleh_dir <- args[[4]]
out <- args[[6]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)

save_csv <- function(x, name) write.csv(x, file.path(out, name), row.names = FALSE, na = "")
save_tsv <- function(x, name) {
    path <- file.path(out, name)
    if (endsWith(path, ".gz")) {
        path <- gzfile(path, "wt")
        on.exit(close(path))
    }
    write.table(x, path, row.names = FALSE, quote = FALSE, sep = "\t", na = "")
}
read_counts <- function(path, sep = ",") {
    x <- as.matrix(read.table(path, header = TRUE, sep = sep, row.names = 1, check.names = FALSE))
    stopifnot(is.numeric(x), all(is.finite(x)), all(x >= 0 & x == floor(x)),
              !anyDuplicated(colnames(x)), all(nzchar(colnames(x))))
    x
}
read_index <- function(path, index, sep = ",") {
    x <- read.table(path, header = TRUE, sep = sep, check.names = FALSE, na.strings = c("", "NA"))
    stopifnot(!anyNA(x[[index]]), !anyDuplicated(x[[index]]))
    rownames(x) <- x[[index]]
    x
}
align_metadata <- function(counts, meta) {
    if (!setequal(colnames(counts), rownames(meta))) {
        stop("Count columns and metadata differ. Count-only: ",
             paste(setdiff(colnames(counts), rownames(meta)), collapse = ", "),
             "; metadata-only: ", paste(setdiff(rownames(meta), colnames(counts)), collapse = ", "))
    }
    meta[colnames(counts), , drop = FALSE]
}
check_groups <- function(group, expected) {
    observed <- table(group)
    if (anyNA(group) || !setequal(names(observed), names(expected)) ||
        any(observed[names(expected)] != expected)) {
        stop("Unexpected clinical groups: ", paste(names(observed), observed, collapse = ", "))
    }
}
sha256 <- function(path) digest::digest(file = path, algo = "sha256")

wt <- function(name) file.path(data_dir, paste0("track01_", name, "_vt.csv"))
db <- function(name) file.path(data_dir, "extension/decidua", paste0("track01_", name, "_db.csv"))
sc <- function(name) file.path(data_dir, "extension/pseudobulk", paste0("track01_", name, ".csv"))
aw <- function(name) file.path(awamleh_dir, paste0("gse114691_", name))
annotation_file <- file.path(data_dir, "track01_gene_annotation.csv")
qc_file <- file.path(data_dir, "qc/track01_htseq_summary.csv")
marker_file <- file.path(data_dir, "interpretation/track01_marker_matrix.csv")
local_de <- file.path(openbio_dir, "vt_contrasts.tsv.gz")
gmt_digest <- sha256(args[[5]])
stopifnot(gmt_digest == "9ea1b5e656597daf423e41c5ebcaa9892bfedf3292fff768605d7b0d5e5e9703")

vt_counts <- read_counts(wt("counts"))
db_counts <- read_counts(db("counts"))
vt_meta <- align_metadata(vt_counts, read_index(wt("metadata"), "sample_id"))
db_meta <- align_metadata(db_counts, read_index(db("metadata"), "sample_id"))
audit <- read_index(args[[2]], "donor_id")
qc <- read_index(qc_file, "sample_id")
all_local_ids <- c(vt_meta$sample_id, db_meta$sample_id)
stopifnot(!anyDuplicated(all_local_ids), setequal(all_local_ids, rownames(qc)),
          setequal(rownames(audit), union(vt_meta$donor_id, db_meta$donor_id)))
audit_tissue <- function(meta, counts, tissue) {
    matched <- audit[meta$donor_id, ]
    stopifnot(!anyNA(meta[c("donor_id", "group", "ga_delivery_weeks", "fetal_sex")]),
              !anyDuplicated(meta$donor_id), identical(meta$group, matched$group),
              identical(meta$sample_id, matched[[paste0(tissue, "_sample_id")]]),
              all(meta$ga_delivery_weeks == matched$ga_delivery_weeks),
              identical(meta$fgr_status, matched$fgr_status), identical(meta$pe_status, matched$pe_status),
              identical(meta$pe_onset, matched$pe_onset), !anyNA(matched$sequencing_round))
    depth <- qc[meta$sample_id, ]
    htseq_total <- rowSums(depth[setdiff(names(depth), "sample_id")])
    stopifnot(all(depth$assigned_to_genes == colSums(counts)), all(htseq_total > 0))
    data.frame(cohort = "OpenBio", tissue = tissue, meta,
               inferred_round = matched$sequencing_round,
               round_evidence = matched$sequencing_round_evidence,
               assigned_gene_counts = depth$assigned_to_genes, htseq_total = htseq_total,
               assignment_fraction = depth$assigned_to_genes / htseq_total)
}
vt <- audit_tissue(vt_meta, vt_counts, "vt")
decidua <- audit_tissue(db_meta, db_counts, "db")
check_groups(vt$group, c("EOPE" = 8, "EOPE+FGR" = 7, "FGR" = 5, "iPTB" = 5, "TB" = 5))
check_groups(decidua$group, c("EOPE" = 8, "EOPE+FGR" = 4, "FGR" = 5, "iPTB" = 5, "TB" = 5))

ext_counts <- read_counts(aw("gene_counts.tsv.gz"), "\t")
ext_meta <- align_metadata(ext_counts, read_index(aw("matched_metadata.tsv"), "accession", "\t"))
ext_labels <- c(combined = "EOPE+FGR", pe = "EOPE", fgr = "FGR", control = "Preterm control")
check_groups(ext_meta$group, c(combined = 20, pe = 20, fgr = 18, control = 21))
external <- data.frame(cohort = "Awamleh", tissue = "vt", sample_id = ext_meta$accession,
                       geo_accession = ext_meta$accession, donor_id = ext_meta$accession,
                       group = unname(ext_labels[ext_meta$group]), geo_diagnosis = ext_meta$geo_diagnosis,
                       source_column = ext_meta$source_column, ga_delivery_weeks = NA_real_,
                       assigned_gene_counts = unname(colSums(ext_counts)))
manifest_columns <- unique(c(names(vt), names(external), "labor", "delivery_mode", "size_factor"))
pad_columns <- function(x) {
    x[setdiff(manifest_columns, names(x))] <- NA
    x[manifest_columns]
}

cases <- vt$group == "EOPE+FGR"
fgr <- vt$group == "FGR"
stopifnot(length(unique(vt$inferred_round[cases | fgr])) == 1L)
age_min <- max(min(vt$ga_delivery_weeks[cases]), min(vt$ga_delivery_weeks[fgr]))
age_max <- min(max(vt$ga_delivery_weeks[cases]), max(vt$ga_delivery_weeks[fgr]))
in_age <- vt$ga_delivery_weeks >= age_min & vt$ga_delivery_weeks <= age_max
stopifnot(age_min == 26.714, age_max == 30.429, sum(cases & in_age) == 2L,
          sum(fgr & in_age) == 5L)
bp_flags <- vt$donor_id[fgr & grepl("no_pe_with_bp_ge_140_90", vt$qc_flags, fixed = TRUE)]
stopifnot(setequal(bp_flags, c("S04", "S15")))

annotation <- read_index(annotation_file, "gene")
mapping <- annotation[complete.cases(annotation[c("hgnc_id", "hgnc_symbol", "ensembl_gene_id")]), ]
unique_values <- function(x) !(duplicated(x) | duplicated(x, fromLast = TRUE))
mapping <- mapping[unique_values(mapping$hgnc_id) & unique_values(mapping$hgnc_symbol) &
                   unique_values(mapping$ensembl_gene_id), ]
pathway_names <- paste0("REACTOME_SIGNALING_BY_", c("BMP", "ACTIVIN", "TGF_BETA_RECEPTOR_COMPLEX"))
pathways <- fgsea::gmtPathways(args[[5]])[pathway_names]
bmp <- pathways[[1]]
stopifnot(length(bmp) == 31L, all(lengths(pathways) >= 10))
bmp_map <- mapping[match(bmp, mapping$hgnc_symbol), ]
coverage <- function(matrix, genes) {
    present <- genes %in% rownames(matrix)
    n <- rep(NA_integer_, length(genes))
    n[present] <- rowSums(matrix[genes[present], , drop = FALSE] >= 10)
    n
}
gene_manifest <- data.frame(pathway_id = "R-HSA-201451", source_symbol = bmp,
                            local_gene = bmp_map$gene, hgnc_id = bmp_map$hgnc_id,
                            ensembl_gene_id = bmp_map$ensembl_gene_id,
                            openbio_libraries_ge10 = coverage(vt_counts, bmp_map$gene),
                            awamleh_libraries_ge10 = coverage(ext_counts, bmp_map$ensembl_gene_id))
gene_manifest$in_score <- with(gene_manifest, !is.na(openbio_libraries_ge10) &
    !is.na(awamleh_libraries_ge10) & openbio_libraries_ge10 >= 5 & awamleh_libraries_ge10 >= 5)
gene_manifest$score_exclusion <- with(gene_manifest, ifelse(in_score, "",
    ifelse(is.na(hgnc_id), "no unique HGNC/Ensembl mapping",
    ifelse(is.na(openbio_libraries_ge10) | is.na(awamleh_libraries_ge10),
           "absent from one or both count matrices", "fewer than five libraries with >=10 reads in one or both cohorts"))))
stopifnot(setequal(bmp[!gene_manifest$in_score], c("AMHR2", "BMP10", "CER1", "GDF2")))
fixed_genes <- bmp_map[gene_manifest$in_score, ]

vt$size_factor <- unname(estimateSizeFactorsForMatrix(vt_counts))
external$size_factor <- unname(estimateSizeFactorsForMatrix(ext_counts))
save_csv(do.call(rbind, lapply(list(vt, decidua, external), pad_columns)), "cohort_manifest.csv")
score_cohort <- function(meta, counts, ids) {
    log_counts <- log2(sweep(counts[ids, , drop = FALSE], 2, meta$size_factor, "/") + 1)
    data.frame(meta[c("cohort", "sample_id", "donor_id", "group", "ga_delivery_weeks", "size_factor")],
               bmp27 = colMeans(log_counts),
               bmp25_without_inhibins = colMeans(log_counts[!fixed_genes$hgnc_symbol %in%
                                                             c("INHA", "INHBA"), , drop = FALSE]),
               row.names = NULL)
}
local_scores <- score_cohort(vt, vt_counts, fixed_genes$gene)
external_scores <- score_cohort(external, ext_counts, fixed_genes$ensembl_gene_id)
save_csv(rbind(local_scores, external_scores), "bmp_donor_scores.csv")
score_fits <- list(openbio_primary = local_scores, openbio_overlap = local_scores[in_age, ],
                   awamleh_primary = external_scores, awamleh_eope_vs_fgr = external_scores,
                   awamleh_fgr_vs_control = external_scores)
for (omitted in c(as.list(bp_flags), list(bp_flags))) {
    name <- paste(c("openbio_without", omitted), collapse = "_")
    score_fits[[name]] <- local_scores[!local_scores$donor_id %in% omitted, ]
}
score_groups <- function(name) {
    switch(name, awamleh_eope_vs_fgr = c("EOPE", "FGR"),
           awamleh_fgr_vs_control = c("FGR", "Preterm control"), c("EOPE+FGR", "FGR"))
}
membership <- lapply(names(score_fits), function(name) {
    arms <- score_groups(name)
    subset <- score_fits[[name]][score_fits[[name]]$group %in% arms, ]
    data.frame(comparison = name, subset,
               arm = ifelse(subset$group == arms[1], "numerator", "denominator"))
})
membership <- do.call(rbind, membership)
save_csv(membership[c("comparison", "cohort", "sample_id", "donor_id", "group", "arm",
                      "ga_delivery_weeks")], "comparison_membership.csv")
age_summary <- function(x) {
    age <- x$ga_delivery_weeks[is.finite(x$ga_delivery_weeks)]
    data.frame(cohort = x$cohort[1], comparison = x$comparison[1], group = x$group[1],
               n = nrow(x), n_age_available = length(age),
               min_weeks = if (length(age)) min(age) else NA_real_,
               median_weeks = if (length(age)) median(age) else NA_real_,
               max_weeks = if (length(age)) max(age) else NA_real_)
}
all_groups <- rbind(data.frame(comparison = "openbio_all_villous", local_scores),
                    data.frame(comparison = "awamleh_all_villous", external_scores))
support <- rbind(membership[names(all_groups)], all_groups)
save_csv(do.call(rbind, lapply(split(support, list(support$comparison, support$group), drop = TRUE),
                              age_summary)), "cohort_support.csv")
welch <- function(x, y) {
    difference <- mean(x) - mean(y)
    interval <- if (length(x) >= 2 && length(y) >= 2) t.test(x, y)$conf.int else c(NA, NA)
    c(effect = difference, ci_low = interval[[1]], ci_high = interval[[2]])
}
score_sets <- c(BMP27 = "bmp27", BMP25_without_INHA_INHBA = "bmp25_without_inhibins")
loo <- list()
effects <- lapply(names(score_fits), function(name) {
    members <- membership[membership$comparison == name, ]
    a <- members$arm == "numerator"
    do.call(rbind, lapply(names(score_sets), function(score) {
        values <- members[[score_sets[[score]]]]
        estimate <- welch(values[a], values[!a])
        leave_one_out <- vapply(seq_along(values), function(i) {
            mean(values[-i][a[-i]]) - mean(values[-i][!a[-i]])
        }, numeric(1))
        loo[[paste(name, score)]] <<- data.frame(
            comparison = name, score = score, omitted_donor_id = members$donor_id,
            omitted_group = members$group, effect = leave_one_out)
        data.frame(cohort = members$cohort[1], comparison = name, score = score,
                   numerator = score_groups(name)[1], denominator = score_groups(name)[2],
                   n_a = sum(a), n_b = sum(!a), as.list(estimate),
                   loo_min = min(leave_one_out), loo_max = max(leave_one_out))
    }))
})
save_csv(do.call(rbind, effects), "bmp_effects.csv")
save_csv(do.call(rbind, loo), "bmp_leave_one_out.csv")

# Optional tissue comparison uses the same shared RNA panel in both tissues.
paired_coverage <- coverage(db_counts, fixed_genes$gene)
paired_genes <- !is.na(paired_coverage) & paired_coverage >= 5
save_csv(data.frame(symbol = fixed_genes$hgnc_symbol, local_gene = fixed_genes$gene,
                    decidua_libraries_ge10 = paired_coverage, in_paired_score = paired_genes),
         "paired_gene_coverage.csv")
both_tissues <- vt[(cases | fgr) & vt$donor_id %in% decidua$donor_id, ]
db_samples <- decidua$sample_id[match(both_tissues$donor_id, decidua$donor_id)]
paired_ids <- fixed_genes$gene[paired_genes]
db_factors <- setNames(estimateSizeFactorsForMatrix(db_counts), colnames(db_counts))
paired_score <- function(counts, samples, factors) {
    colMeans(log2(sweep(counts[paired_ids, samples, drop = FALSE], 2, factors, "/") + 1))
}
paired <- data.frame(donor_id = both_tissues$donor_id, group = both_tissues$group,
    ga_delivery_weeks = both_tissues$ga_delivery_weeks, vt_sample_id = both_tissues$sample_id,
    db_sample_id = db_samples, villi = unname(paired_score(vt_counts, both_tissues$sample_id, both_tissues$size_factor)),
    decidua = unname(paired_score(db_counts, db_samples, db_factors[db_samples])))
paired$decidua_minus_villi <- paired$decidua - paired$villi
paired_cases <- paired$group == "EOPE+FGR"
stopifnot(sum(paired_cases) == 4L, sum(!paired_cases) == 5L)
save_csv(paired, "paired_donor_scores.csv")
paired_effects <- lapply(c("villi", "decidua", "decidua_minus_villi"), function(name) {
    estimate <- welch(paired[[name]][paired_cases], paired[[name]][!paired_cases])
    data.frame(comparison = name, n_genes = length(paired_ids), n_a = sum(paired_cases),
               n_b = sum(!paired_cases), as.list(estimate))
})
save_csv(do.call(rbind, paired_effects), "paired_tissue_effects.csv")

full_model_meta <- data.frame(fetal_sex = factor(vt$fetal_sex),
    group = factor(vt$group, levels = c("EOPE", "EOPE+FGR", "FGR", "iPTB", "TB"),
                   labels = c("EOPE", "EOPE_FGR", "FGR", "iPTB", "TB")), row.names = vt$sample_id)
design <- ~ fetal_sex + group
overlap_meta <- droplevels(full_model_meta[(cases | fgr) & in_age, ])
design_matrix <- model.matrix(design, overlap_meta)
stopifnot(qr(design_matrix)$rank == ncol(design_matrix), nrow(design_matrix) > ncol(design_matrix))
overlap_dds <- DESeqDataSetFromMatrix(vt_counts[, rownames(overlap_meta)], overlap_meta, design)
overlap_dds <- DESeq(overlap_dds, quiet = TRUE)
overlap_de <- as.data.frame(results(overlap_dds, contrast = c("group", "EOPE_FGR", "FGR"), alpha = .05))
overlap_de$gene <- rownames(overlap_de)
overlap_de$low95 <- overlap_de$log2FoldChange - qnorm(.975) * overlap_de$lfcSE
overlap_de$high95 <- overlap_de$log2FoldChange + qnorm(.975) * overlap_de$lfcSE
both_zero <- rowSums(vt_counts[overlap_de$gene, rownames(overlap_meta)]) == 0
overlap_de[both_zero, c("low95", "high95")] <- NA_real_

pca_dds <- DESeqDataSetFromMatrix(vt_counts, full_model_meta, design)
sizeFactors(pca_dds) <- vt$size_factor
transformed <- assay(varianceStabilizingTransformation(pca_dds, blind = TRUE))
top <- head(order(matrixStats::rowVars(transformed), decreasing = TRUE), 500)
pca <- prcomp(t(transformed[top, ]))
pca_percent <- 100 * pca$sdev^2 / sum(pca$sdev^2)
save_csv(data.frame(vt[c("sample_id", "donor_id", "group", "fetal_sex", "inferred_round")],
                    pca$x[, 1:2], pc1_percent = pca_percent[1], pc2_percent = pca_percent[2]), "pca_scores.csv")

full_de <- read.delim(local_de, check.names = FALSE)
local_pairs <- rbind(c("EOPE+FGR", "FGR"), c("EOPE", "FGR"), c("EOPE+FGR", "iPTB"),
                     c("EOPE", "iPTB"), c("FGR", "iPTB"))
ext_pairs <- rbind(c("combined", "fgr"), c("pe", "fgr"), c("combined", "control"),
                   c("pe", "control"), c("fgr", "control"))
suffix <- c("primary", "eope_vs_fgr", "combined_vs_control", "eope_vs_control", "fgr_vs_control")
ext_de <- aw(paste0(ext_pairs[, 1], "_vs_", ext_pairs[, 2], ".tsv.gz"))
effect_rows <- function(x, fit, cohort, arm_labels, ids, id_column) {
    stopifnot(!anyDuplicated(x[[id_column]]), nrow(x) > 0)
    mapped <- mapping[match(x[[id_column]], ids), ]
    data.frame(fit = fit, cohort = cohort, numerator = arm_labels[1], denominator = arm_labels[2],
               input_id = x[[id_column]], symbol = mapped$hgnc_symbol, ensembl_gene_id = mapped$ensembl_gene_id,
               x[c("baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "low95", "high95")],
               both_groups_zero = if ("status" %in% names(x)) startsWith(x$status, "all_zero") else FALSE,
               row.names = NULL)
}
rank_fits <- c(
    lapply(seq_len(nrow(local_pairs)), function(i) {
        arms <- local_pairs[i, ]
        rows <- full_de[full_de$numerator == arms[1] & full_de$denominator == arms[2], ]
        effect_rows(rows, paste0("openbio_", suffix[i]), "OpenBio", arms, mapping$gene, "gene")
    }),
    list(effect_rows(overlap_de, "openbio_overlap", "OpenBio", local_pairs[1, ], mapping$gene, "gene")),
    lapply(seq_len(nrow(ext_pairs)), function(i) {
        effect_rows(read.delim(ext_de[i], check.names = FALSE), paste0("awamleh_", suffix[i]), "Awamleh",
                    ext_labels[ext_pairs[i, ]], mapping$ensembl_gene_id, "ensembl_gene_id")
    })
)
rank_fits <- setNames(rank_fits, vapply(rank_fits, function(x) x$fit[1], character(1)))
model_membership <- lapply(names(rank_fits), function(fit) {
    model <- rank_fits[[fit]]
    cohort <- if (model$cohort[1] == "OpenBio") vt else external
    if (fit == "openbio_overlap") cohort <- cohort[cohort$sample_id %in% rownames(overlap_meta), ]
    columns <- c("cohort", "sample_id", "donor_id", "group", "ga_delivery_weeks")
    role <- ifelse(cohort$group == model$numerator[1], "numerator",
                   ifelse(cohort$group == model$denominator[1], "denominator", "model_only"))
    data.frame(fit = fit, cohort[columns], role = role)
})
save_csv(do.call(rbind, model_membership), "rank_fit_membership.csv")
eligible <- lapply(rank_fits, function(x) {
    # Resolve rank mappings among genes eligible in each fitted contrast, as in the pilot.
    key <- if (x$cohort[1] == "OpenBio") "gene" else "ensembl_gene_id"
    candidate_map <- unique(annotation[!is.na(annotation$hgnc_symbol) & !is.na(annotation[[key]]),
                                       c(key, "hgnc_symbol")])
    candidates <- subset(x, !both_groups_zero & is.finite(stat) & is.finite(pvalue) & baseMean >= 10)
    joined <- merge(candidates, candidate_map, by.x = "input_id", by.y = key)
    joined <- joined[unique_values(joined$hgnc_symbol) & unique_values(joined$input_id), ]
    joined$symbol <- joined$hgnc_symbol
    joined
})
universe <- Reduce(intersect, lapply(eligible, function(x) x$symbol))
external_rank_ids <- eligible[["awamleh_primary"]]
universe_ids <- external_rank_ids$input_id[match(universe, external_rank_ids$symbol)]
gene_manifest$in_rank_universe <- bmp %in% universe
stopifnot(sum(gene_manifest$in_rank_universe) == 25)
save_csv(gene_manifest, "bmp_gene_manifest.csv")
save_csv(data.frame(symbol = universe,
                    ensembl_gene_id = universe_ids,
                    in_bmp = universe %in% bmp,
                    in_activin = universe %in% pathways[[2]],
                    in_tgf_beta_receptor = universe %in% pathways[[3]]), "rank_universe.csv")
tests <- do.call(rbind, lapply(names(rank_fits), function(fit) {
    x <- rank_fits[[fit]]
    z <- eligible[[fit]]
    matched <- match(x$input_id, z$input_id)
    x$symbol[!is.na(matched)] <- z$symbol[matched[!is.na(matched)]]
    x$in_rank_universe <- !is.na(matched) & x$symbol %in% universe
    x$ensembl_gene_id[x$in_rank_universe] <- universe_ids[match(x$symbol[x$in_rank_universe], universe)]
    x
}))
tests$in_bmp <- tests$symbol %in% bmp
save_tsv(tests, "bmp_gene_effects.tsv.gz")

ranks <- lapply(eligible, function(x) {
    stat <- setNames(x$stat[match(universe, x$symbol)], universe)
    sort(stat, decreasing = TRUE)
})
branch <- lapply(names(ranks), function(fit) {
    set.seed(1049)
    z <- fgsea::fgseaMultilevel(pathways, ranks[[fit]], minSize = 10, maxSize = 500, eps = 1e-10,
                               BPPARAM = BiocParallel::SerialParam())
    stopifnot(nrow(z) == 3, all(is.finite(z$pval)))
    data.frame(fit = fit, rank_size = length(ranks[[fit]]), pathway = z$pathway,
               pvalue = z$pval, fdr = p.adjust(z$pval, "BH"), log2err = z$log2err,
               ES = z$ES, NES = z$NES, size = z$size,
               leading_edge = vapply(z$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1)))
})
branch <- do.call(rbind, branch)
save_csv(branch, "branch_enrichment.csv")
main_fits <- c("openbio_primary", "openbio_overlap", "awamleh_primary")
bmp_results <- branch[branch$pathway == pathway_names[1] & branch$fit %in% main_fits, ]
shared_edge <- sort(Reduce(intersect, strsplit(bmp_results$leading_edge, ";", fixed = TRUE)))
stopifnot(length(shared_edge) > 0)
save_csv(data.frame(symbol = shared_edge), "shared_leading_edge.csv")
rank_plot <- lapply(main_fits, function(fit) {
    stat <- ranks[[fit]]
    hit <- names(stat) %in% bmp
    hits <- abs(stat) * hit
    data.frame(fit = fit, rank = seq_along(stat), symbol = names(stat), wald = stat, in_bmp = hit,
               shared_edge = names(stat) %in% shared_edge,
               running_score = cumsum(hits / sum(hits) - (!hit) / sum(!hit)))
})
save_csv(do.call(rbind, rank_plot), "rank_plot.csv")

atlas_counts <- read_counts(sc("pseudobulk_counts"))
atlas_meta <- align_metadata(atlas_counts, read_index(sc("pseudobulk_metadata"), "pseudobulk_id"))
atlas_donors <- read_index(sc("pseudobulk_donors"), "donor_id")
class_map <- read.csv(sc("cell_class_map"))
markers <- read.csv(marker_file)
stopifnot(!anyDuplicated(atlas_meta[c("donor_id", "cell_class")]),
          all(atlas_meta$donor_id %in% atlas_donors$donor_id),
          all(atlas_meta$total_counts == colSums(atlas_counts)), all(atlas_meta$total_counts > 0))
classes <- sort(unique(c(class_map$cell_class, atlas_meta$cell_class, markers$cell_class)))
atlas_genes <- sort(unique(c(shared_edge, "BMPR2", "BMPR1B", "GREM2", "SMAD9")))
atlas <- expand.grid(pseudobulk_id = atlas_meta$pseudobulk_id, symbol = atlas_genes,
                     stringsAsFactors = FALSE)
atlas <- cbind(atlas, atlas_meta[atlas$pseudobulk_id,
    c("donor_id", "cell_class", "clinical_group", "pe_status", "n_cells", "total_counts")])
atlas$raw_count <- atlas_counts[cbind(match(atlas$symbol, rownames(atlas_counts)),
                                     match(atlas$pseudobulk_id, colnames(atlas_counts)))]
atlas$cpm <- 1e6 * atlas$raw_count / atlas$total_counts
atlas$eligible_donor_class <- atlas$n_cells >= 20L
save_csv(atlas, "atlas_donors.csv")
atlas_summary <- lapply(classes, function(cell_class) {
    samples <- atlas_meta[atlas_meta$cell_class == cell_class, ]
    eligible <- samples$n_cells >= 20L
    n <- sum(eligible)
    do.call(rbind, lapply(atlas_genes, function(gene) {
        has_gene <- gene %in% rownames(atlas_counts)
        values <- atlas[atlas$cell_class == cell_class & atlas$symbol == gene & atlas$eligible_donor_class, ]
        status <- if (!nrow(samples)) "class absent" else if (n < 3L) "fewer than 3 eligible donors" else
            if (!has_gene) "gene absent" else "reported"
        data.frame(cell_class = cell_class, symbol = gene, donors_released = nrow(samples),
                   donors_eligible = n, status = status,
                   median_cpm = if (status == "reported") median(values$cpm) else NA_real_,
                   fraction_ge10 = if (status == "reported") mean(values$raw_count >= 10L) else NA_real_)
    }))
})
save_csv(do.call(rbind, atlas_summary), "atlas_summary.csv")

marker_classes <- c("villous_cytotrophoblast", "syncytiotrophoblast", "extravillous_trophoblast",
                    "endothelial_cell", "perivascular_cell", "decidual_stromal_cell", "decidual_macrophage")
markers <- unique(markers[markers$cell_class %in% marker_classes & markers$gene_symbol != "MKI67", ])
marker_rows <- lapply(main_fits, function(fit) {
    x <- rank_fits[[fit]]
    joined <- cbind(markers, x[match(markers$gene_symbol, x$symbol),
                               c("input_id", "log2FoldChange", "low95", "high95", "pvalue", "padj")])
    data.frame(fit = fit, joined, row.names = NULL)
})
save_csv(do.call(rbind, marker_rows), "bulk_marker_effects.csv")

file_inputs <- c(wt("counts"), wt("metadata"), db("counts"), db("metadata"), annotation_file,
                 qc_file, args[[2]], args[[5]], sc("pseudobulk_counts"), sc("pseudobulk_metadata"),
                 sc("pseudobulk_donors"), sc("cell_class_map"), marker_file, local_de,
                 aw("gene_counts.tsv.gz"), aw("matched_metadata.tsv"), ext_de)
manifest <- data.frame(kind = "file_sha256", name = file_inputs,
                       value = vapply(file_inputs, sha256, character(1)), row.names = NULL)
settings <- c(R = as.character(getRversion()),
              vapply(c("DESeq2", "fgsea", "BiocParallel", "digest", "matrixStats", "ggplot2",
                       "patchwork", "ragg"), function(p) as.character(packageVersion(p)), character(1)),
              normalization = "DESeq2 median ratio; fit separately within each full villous cohort",
              score = "unweighted mean log2(normalized count + 1) over 27 fixed genes",
              openbio_fit = "~ fetal_sex + group; Wald test; DESeq2 independent filtering at alpha 0.05",
              age_fit = "~ fetal_sex + group; refit on common age support; score retains full size factors",
              awamleh_fit = "PyDESeq2 ~ group; >=10 raw reads in >=5 libraries; Cook filtering/refit off",
              fgsea = "seed 1049 reset per fit, SerialParam, minSize 10, maxSize 500, eps 1e-10; BH within fit",
              overlap_weeks = paste(age_min, age_max, sep = "-"), rank_universe_size = length(universe),
              n_openbio_vt = nrow(vt), n_openbio_db = nrow(decidua), n_awamleh_vt = nrow(external),
              n_atlas_donors = nrow(atlas_donors), n_atlas_donor_classes = nrow(atlas_meta))
manifest <- rbind(manifest, data.frame(kind = "setting", name = names(settings), value = unname(settings)))
provenance <- file.path(awamleh_dir, "provenance.tsv")
if (file.exists(provenance)) manifest <- rbind(manifest, read.delim(provenance))
save_tsv(manifest, "run_manifest.tsv")
capture.output(sessionInfo(), file = file.path(out, "r_session_info.txt"))

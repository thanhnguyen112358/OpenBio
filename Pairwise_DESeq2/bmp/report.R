#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(ggplot2))
library(patchwork)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: report.R OUTPUT_DIR")
out <- args[[1]]
figures <- file.path(out, "figures")
dir.create(figures, showWarnings = FALSE)
data <- function(name) read.csv(file.path(out, paste0(name, ".csv")), check.names = FALSE)
save_figure <- function(plot, name, width, height) {
    ggsave(file.path(figures, paste0(name, ".png")), plot, width = width, height = height,
           dpi = 170, bg = "#FFFFFF", device = ragg::agg_png)
}

ink <- "#182D3D"
muted <- "#4D6271"
accent <- "#005E68"
combined <- "#A14325"
group_colors <- c(FGR = accent, "EOPE+FGR" = combined)
theme_set(theme_minimal(base_size = 13) +
    theme(text = element_text(colour = ink), plot.title = element_text(face = "bold", size = 16),
          plot.subtitle = element_text(colour = muted, margin = margin(b = 12)),
          plot.caption = element_text(colour = muted, hjust = 0, lineheight = 1.1,
                                      margin = margin(t = 10)),
          axis.text = element_text(colour = ink), panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#E6EEF0", linewidth = .35),
          plot.margin = margin(12, 18, 12, 12)))
caption_style <- theme(plot.title = element_text(face = "bold", size = 21),
                       plot.subtitle = element_text(size = 14, colour = muted),
                       plot.caption = element_text(size = 11, colour = muted, hjust = 0))

membership <- data("comparison_membership")
effects <- subset(data("bmp_effects"), score == "BMP27")
full <- subset(membership, comparison == "openbio_primary")
overlap <- subset(membership, comparison == "openbio_overlap")
cohort <- data("cohort_manifest")
full$qc_flags <- cohort$qc_flags[match(full$sample_id, cohort$sample_id)]
full$flag <- ifelse(full$qc_flags == "none", "Released label", "Blood pressure flag")
full$arm_plot <- ifelse(full$group == "FGR", 1, 2)
full <- full[order(full$group, full$ga_delivery_weeks, full$donor_id), ]
full$y <- full$arm_plot + ave(full$ga_delivery_weeks, full$group, FUN = function(age) {
    cluster <- cumsum(c(TRUE, diff(age) > .38))
    ave(age, cluster, FUN = function(x) rep(c(-.07, .12, -.23), length.out = length(x)))
})
age <- ggplot(full, aes(ga_delivery_weeks, y, color = group)) +
    annotate("rect", xmin = min(overlap$ga_delivery_weeks), xmax = max(overlap$ga_delivery_weeks),
             ymin = .65, ymax = 2.45, fill = "#DBEEEE") +
    annotate("text", x = mean(range(overlap$ga_delivery_weeks)), y = 2.35,
             label = "Shared delivery range", size = 4.2, colour = accent) +
    geom_point(aes(shape = flag), size = 3.8) +
    geom_text(aes(label = donor_id), nudge_y = .095, size = 3.5, show.legend = FALSE) +
    scale_color_manual(values = group_colors) +
    scale_shape_manual(values = c("Released label" = 16, "Blood pressure flag" = 17)) +
    scale_y_continuous(breaks = 1:2, labels = c("Isolated FGR (5)", "EOPE + FGR (7)"), limits = c(.65, 2.5)) +
    scale_x_continuous(breaks = 26:34, limits = c(25.9, 34.3)) +
    guides(color = "none") +
    labs(title = "Only two combined cases share the isolated-FGR delivery range",
         subtitle = sprintf("OpenBio villi  |  %.3f to %.3f weeks  |  2 combined cases and 5 isolated FGR",
                            min(overlap$ga_delivery_weeks), max(overlap$ga_delivery_weeks)),
         x = "Gestational age at delivery (weeks)", y = NULL, shape = NULL,
         caption = "Each point is a pregnancy. Triangles mark blood pressure flags recorded with the released diagnosis.") +
    theme(panel.grid.major.y = element_blank(), legend.position = "bottom")
save_figure(age, "age_support", 11.4, 5.1)

fit_names <- c(openbio_primary = "OpenBio", openbio_overlap = "OpenBio (shared age)",
               awamleh_primary = "Awamleh")
core <- effects[match(names(fit_names), effects$comparison), ]
scores <- data("bmp_donor_scores")
donor_panel <- function(fit) {
    selected <- membership[membership$comparison == fit, ]
    selected$score <- scores$bmp27[match(paste(selected$cohort, selected$donor_id),
                                        paste(scores$cohort, scores$donor_id))]
    selected$group <- factor(selected$group, levels = c("FGR", "EOPE+FGR"))
    labels <- paste0(c("Isolated FGR", "EOPE + FGR"), "\n", "n = ", as.integer(table(selected$group)))
    ggplot(selected, aes(group, score, colour = group)) +
        geom_point(position = position_jitter(width = .15, height = 0, seed = 1049), size = 2.7, alpha = .74) +
        stat_summary(fun = mean, geom = "point", shape = 95, size = 13, colour = ink) +
        scale_colour_manual(values = group_colors) +
        scale_x_discrete(labels = labels) +
        labs(title = unname(fit_names[fit]), y = "Average BMP RNA score (log2 units)", x = NULL) +
        theme(legend.position = "none", panel.grid.major.x = element_blank(), plot.title = element_text(size = 15))
}
core$label <- factor(fit_names[core$comparison], levels = rev(fit_names))
core$estimate_label <- sprintf("%+.3f", core$effect)
forest <- ggplot(core, aes(x = effect, y = label)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "#8A9BA6") +
    geom_segment(aes(x = ci_low, xend = ci_high, yend = label), linewidth = 1, color = accent) +
    geom_point(size = 3.4, color = accent) +
    geom_text(aes(label = estimate_label), color = ink, nudge_y = .26, size = 3.6) +
    scale_x_continuous(breaks = seq(-1, 1.5, .25), limits = c(-.85, 1.36)) +
    scale_y_discrete(expand = expansion(add = .6)) +
    labs(title = "Difference between group means", x = "EOPE + FGR minus isolated FGR (log2 units)",
         y = NULL) +
    theme(panel.grid.major.y = element_blank(), plot.title = element_text(size = 15))
headline <- ((donor_panel("openbio_primary") | donor_panel("awamleh_primary")) / forest) +
    plot_layout(heights = c(1.14, 1)) +
    plot_annotation(title = "Placental BMP-related RNA was higher with early preeclampsia",
        subtitle = "Two cohorts of pregnancies with fetal growth restriction",
        caption = "Points in the upper panels are pregnancies and bars mark group means. Cohorts use separate score scales.\nLines below show approximate 95% intervals. Just two OpenBio combined cases fall in the shared age range.",
        theme = caption_style)
save_figure(headline, "headline", 12.4, 8.4)

rank <- data("rank_plot")
enrichment <- subset(data("branch_enrichment"), pathway == "REACTOME_SIGNALING_BY_BMP")
enrichment <- enrichment[match(names(fit_names), enrichment$fit), ]
panels <- sprintf("%s\nNES %.2f  |  FDR %.3f", fit_names, enrichment$NES, enrichment$fdr)
rank$panel <- factor(panels[match(rank$fit, names(fit_names))], levels = panels)
rank$hit_type <- ifelse(rank$shared_edge, "Shared leading edge", "Other BMP member")
plot_rank <- ggplot(rank, aes(rank, running_score)) +
    geom_hline(yintercept = 0, colour = "#8A9BA6", linewidth = .55) +
    geom_line(colour = accent, linewidth = .87) +
    geom_segment(data = subset(rank, in_bmp), aes(xend = rank, y = -.11, yend = -.035, colour = hit_type),
                 linewidth = .65, alpha = .85) +
    geom_point(data = subset(rank, symbol == "BMPR2"), color = combined, size = 2.5) +
    geom_text(data = subset(rank, symbol == "BMPR2"), label = "BMPR2", hjust = -.14, vjust = -.3,
              colour = combined, size = 3.5, fontface = "bold") +
    facet_wrap(~panel, nrow = 1) +
    scale_x_continuous(breaks = c(0, 4500, 9000, 13500), labels = c("0", "4,500", "9,000", "13,500")) +
    scale_color_manual(values = c("Shared leading edge" = combined, "Other BMP member" = "#8A9BA6")) +
    labs(title = "BMP genes cluster toward the higher-expression end of the rank",
         subtitle = sprintf("%s common ranked genes, including 25 Reactome BMP genes",
                            format(nrow(rank) / length(fit_names), big.mark = ",")),
         x = "Gene rank by signed Wald statistic (higher with EOPE + FGR at left)",
         y = "Running enrichment score", colour = NULL,
         caption = "Ticks show BMP genes. The shared leading edge contains genes that contribute to positive enrichment in all three fits.\nNES is the normalized enrichment score. FDR adjusts for the three Reactome branches per fit.") +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = 11),
          axis.text.x = element_text(size = 9))
save_figure(plot_rank, "enrichment", 12, 5.85)

atlas <- data("atlas_summary")
class_title <- function(x) {
    x <- gsub("_", " ", tools::toTitleCase(x))
    x <- gsub("\\bb cell\\b", "B cell", x, ignore.case = TRUE)
    x <- gsub("\\bt cell\\b", "T cell", x, ignore.case = TRUE)
    sub("\\bnk cell\\b", "NK cell", x, ignore.case = TRUE)
}
atlas$class <- paste0(class_title(atlas$cell_class), "  (", atlas$donors_eligible, ")")
atlas$class <- factor(atlas$class, levels = rev(unique(atlas$class)))
atlas$symbol <- factor(atlas$symbol, levels = c("BMPR2", "BMPR1B", "SMAD9", "SMAD4", "GREM2",
                                               "INHA", "INHBA",
                                               setdiff(unique(atlas$symbol),
                                                       c("BMPR2", "BMPR1B", "SMAD9", "SMAD4",
                                                         "GREM2", "INHA", "INHBA"))))
atlas_plot <- ggplot(atlas, aes(symbol, class)) +
    geom_point(data = subset(atlas, status != "reported"), shape = 4, size = 2, color = "#A4B3BA") +
    geom_point(data = subset(atlas, status == "reported"),
               aes(size = fraction_ge10, colour = log10(median_cpm + 1))) +
    scale_color_gradient(low = "#D1DDE1", high = "#00535F", name = "Median RNA\n(counts per million)",
                         breaks = log10(c(1, 2, 11, 101, 1001)), labels = c("0", "1", "10", "100", "1,000")) +
    scale_size_continuous(name = "Donors with at\nleast 10 reads", range = c(1.1, 7.2), limits = c(0, 1),
                          breaks = c(0, .5, 1), labels = c("0%", "50%", "100%")) +
    labs(title = "BMP-related transcripts occur in several placental cell classes",
         subtitle = "Admati donor-by-cell-class reference  |  Numbers in parentheses are eligible donors",
         x = "Reactome BMP gene", y = NULL,
         caption = "Each contributing donor-class has at least 20 cells. Crosses mark unavailable summaries.\nThe reference includes preeclampsia and controls and has no isolated-FGR comparison arm.") +
    theme(axis.text.x = element_text(face = "bold", size = 11), panel.grid.major.x = element_blank())
save_figure(atlas_plot, "atlas", 12, 8.6)

pca <- data("pca_scores")
pca_plot <- function(variable, title) {
    ggplot(pca, aes(PC1, PC2, color = .data[[variable]])) +
        geom_point(size = 2.65) +
        geom_text(aes(label = donor_id), nudge_y = 2.3, size = 2.8, show.legend = FALSE,
                  check_overlap = TRUE) +
        labs(title = title, color = NULL,
             x = sprintf("PC1 (%.1f%%)", pca$pc1_percent[1]), y = sprintf("PC2 (%.1f%%)", pca$pc2_percent[1])) +
        theme(legend.position = "bottom", legend.text = element_text(size = 8))
}
pca_all <- (pca_plot("group", "Diagnosis") | pca_plot("fetal_sex", "Fetal sex") |
            pca_plot("inferred_round", "Inferred sequencing round")) +
    plot_annotation(title = "OpenBio villous gene-expression profiles",
                    subtitle = "PCA on the 500 most variable genes after variance stabilization",
                    caption = "Each point is one pregnancy. Inferred sequencing round is based on the original study methods.",
                    theme = caption_style)
save_figure(pca_all, "pca", 14.8, 5.6)

connection <- gzfile(file.path(out, "bmp_gene_effects.tsv.gz"), "rt")
de <- read.delim(connection, check.names = FALSE)
close(connection)
de <- subset(de, fit == "openbio_primary" & is.finite(pvalue))
de$highlight <- ifelse(de$in_bmp, "Reactome BMP", "Other genes")
de$mlogp <- -log10(pmax(de$pvalue, .Machine$double.xmin))
volcano <- ggplot(de, aes(log2FoldChange, mlogp)) +
    geom_point(color = "#9EB1BA", alpha = .28, size = .8) +
    geom_point(data = subset(de, in_bmp), color = accent, size = 1.9) +
    geom_text(data = subset(de, symbol == "BMPR2"), label = "BMPR2", nudge_x = .18, hjust = 0,
              vjust = -.35, color = accent, fontface = "bold") +
    labs(title = "Gene-level differences in OpenBio villi",
         subtitle = "EOPE + FGR versus isolated FGR  |  Sex-adjusted DESeq2",
         x = "Estimated log2 fold change (positive is higher with EOPE + FGR)",
         y = "-log10 of Wald p-value",
         caption = "Teal highlights Reactome BMP genes. The vertical scale uses unadjusted p-values.\nGene-wide false discovery rates and Wald intervals are in the exported table.")
save_figure(volcano, "volcano", 9.2, 5.2)

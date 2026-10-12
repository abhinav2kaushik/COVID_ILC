### WARNING: COLORS of Network layout may differ ###

library(Seurat)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(VennDiagram)
library(RColorBrewer)
library(circlize)

setwd("/Users/abk347/Library/CloudStorage/OneDrive-HarvardUniversity/ILC/Analysis/Result_replication")

## You may have to change this libpath depending upon your library path here sclink is installed
.libPaths(c(
  "./src/",
  .libPaths()
))
library(scLink)
library(MAST)



fig2_dir <- "Figure_2"
dir.create(fig2_dir, showWarnings = FALSE, recursive = TRUE)
ilc_types <- c("ILC1", "ILC2", "ILCp")
object_file <- "covid_ilc_counts.rds"
covid_ilc <- readRDS(object_file) ## load object we created in step_1


healthy <- c("1033", "1348", "1804", "1987", "2195", "55689-0068", "55689-0071",
             "55689-0080", "WF001", "WF002", "WF006", "WF008", "WF010", "WF015",
             "WF016", "WF018", "WF019", "WF021", "WF025", "WF045", "WF026",
             "WF031", "WF033", "WF035", "WF042")
long_ter_covd <- c("55689-0007", "55689-0012", "55689-0015", "55689-0029",
                   "55689-0031", "55689-0034", "55689-0044", "55689-0050",
                   "55689-0057", "55689-0095", "55689-0165", "55689-0219",
                   "55689-0220", "OHC-001", "OHC-010", "OHC001", "OHC010")


rna <- covid_ilc
DefaultAssay(rna) <- "RNA"
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000,
                     verbose = FALSE)
rna <- ScaleData(rna, features = rownames(rna[["RNA"]]), verbose = FALSE)
rna$COVID_status <- ifelse(rna$PPID %in% healthy, "healthy", "COVID")
metadata2 <- rna[[]]
normalized_rna <- GetAssayData(rna, assay = "RNA", layer = "data")
scaled_rna <- GetAssayData(rna, assay = "RNA", layer = "scale.data")

Idents(rna) <- "handgated_cell_type2"

DE_RNA <- list()
for (cell_type in ilc_types) {
  de <- FindMarkers(rna, group.by = "COVID_status", ident.1 = "COVID", ident.2 = "healthy",
                    subset.ident = cell_type, assay = "RNA", slot = "data",
                    test.use = "MAST", latent.vars = "batch", logfc.threshold = 0,
                    min.pct = 0.1,
                    mean.fxn = function(x) log2(Matrix::rowMeans(expm1(x)) + 1),
                    verbose = FALSE)
  de$Gene <- rownames(de)
  de$ILCn <- cell_type
  DE_RNA[[cell_type]] <- de
  #write.csv(de, file.path(fig2_dir, paste0("DE_", cell_type, ".csv")), row.names = FALSE)
}
all_DE <- do.call(rbind, DE_RNA)

##################################
########## Figure 2A #############
##################################

volcano_p_cutoff <- 1e-5
volcano_fc_cutoff <- 1
volcano_labels <- list(ILC1 = c("IFITM2", "IFITM3", "SELL", "PIK3AP1", "FAM65B", "STAT5A", "CD6"),
                       ILC2 = c("IFITM2", "IFITM3", "SELL", "PIK3AP1", "FAM65B", "SPOCK2", "IL4R", "ICOS"),
                       ILCp = c("IFITM2", "IFITM3", "PIK3AP1", "SELL", "IL4R", "CCR7"))
volcano_plots <- list()
for (cell_type in ilc_types) {
  d <- DE_RNA[[cell_type]]
  d$Significance <- "Neither"
  d$Significance[abs(d$avg_log2FC) > volcano_fc_cutoff] <- "Fold change only"
  d$Significance[d$p_val_adj < volcano_p_cutoff] <- "Adjusted p only"
  d$Significance[d$p_val_adj < volcano_p_cutoff & abs(d$avg_log2FC) > volcano_fc_cutoff] <- "Both"
  volcano_plots[[cell_type]] <- ggplot(d, aes(avg_log2FC, -log10(pmax(p_val_adj, .Machine$double.xmin)))) +
    geom_point(aes(colour = Significance), size = 2, alpha = 0.7) +
    geom_vline(xintercept = c(-1, 1), linetype = 2) +
    geom_hline(yintercept = -log10(volcano_p_cutoff), linetype = 2) +
    geom_text_repel(data = d[d$Gene %in% volcano_labels[[cell_type]], ],
                    aes(label = Gene), fontface = "italic", size = 3, seed = 42, max.overlaps = Inf) +
    scale_colour_manual(values = c("Neither" = "grey40", "Fold change only" = "#4DAF4A",
                                   "Adjusted p only" = "#6685BD", "Both" = "#EE6268")) +
    labs(title = cell_type, x = expression(Log[2]~"fold change"),
          y = expression(-Log[10]~"adjusted p-value")) +
    theme_bw(base_size = 12) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"))
}
panel_A <- arrangeGrob(grobs = volcano_plots, ncol = 3, top = "A   Batch-adjusted RNA differential expression",
                       bottom = "Blue: adjusted p < 1e-5 only; red: adjusted p < 1e-5 and |log2FC| > 1; green: |log2FC| > 1 only")
ggsave(file.path(fig2_dir, "Figure_2A.png"), panel_A, width = 10, height = 4.2, dpi = 300, bg = "white")


##################################
########## Figure 2B #############
##################################

gene_sets <- lapply(DE_RNA, function(d) d$Gene[d$p_val_adj < 0.05 & abs(d$avg_log2FC) > 0.5])
venn <- venn.diagram(gene_sets, filename = NULL, category.names = ilc_types,
                      print.mode = c("raw", "percent"), fill = c("white", "white", "#E9EDE4"),
                      alpha = 0.35, col = "black", lwd = 2, cex = 1.2,
                      cat.cex = 1.3, cat.fontface = "bold", disable.logging = TRUE)


png(file.path(fig2_dir, "Figure_2B.png"), width = 1500, height = 1350, res = 300)
grid::grid.newpage(); grid::grid.draw(venn); grid::grid.text("B", x = 0.03, y = 0.98, gp = grid::gpar(fontface = "bold", fontsize = 16)); dev.off()
venn_genes <- data.frame(Gene = sort(unique(unlist(gene_sets))))
for (cell_type in ilc_types) venn_genes[[cell_type]] <- venn_genes$Gene %in% gene_sets[[cell_type]]
#write.csv(venn_genes, file.path(fig2_dir, "Venn_gene_membership.csv"), row.names = FALSE)

##################################
########## Figure 2C #############
##################################

heatmap_genes <- c("IFITM2", "SELL", "FAM65B", "PIK3AP1", "SPOCK2", "CD6", "STAT5A",
                    "TNFRSF25", "TGFB1", "ITGB7", "IL32", "HPGDS", "ICOS", "IL4R", "ICAM1", "CCR7")
heatmaps <- list()
for (cell_type in ilc_types) {
  mean_expression <- sapply(c("healthy", "COVID"), function(status) {
    cells <- rownames(metadata2)[metadata2$handgated_cell_type2 %in% cell_type & metadata2$COVID_status == status]
    rowMeans(scaled_rna[heatmap_genes, cells, drop = FALSE])
  })
  colnames(mean_expression) <- c("Healthy", "CoV")
  
  #write.csv(mean_expression, file.path(fig2_dir, paste0("Heatmap_means_", cell_type, ".csv")))
  heat_data <- as.data.frame(as.table(mean_expression))
  names(heat_data) <- c("Gene", "Status", "Mean")
  heatmaps[[cell_type]] <- ggplot(heat_data, aes(Status, Gene, fill = Mean)) +
    geom_tile(colour = "grey45", linewidth = 0.3) +
    scale_y_discrete(limits = rev(heatmap_genes)) +
    scale_fill_gradientn(colours = rev(brewer.pal(11, "RdYlBu")), name = NULL) +
    guides(fill = guide_colourbar(direction = "horizontal", barwidth = grid::unit(1.3, "in"),
                                  barheight = grid::unit(0.12, "in"))) +
    labs(title = cell_type, x = NULL, y = NULL) + theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), legend.position = "top",
          legend.text = element_text(size = 8), plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.text.y = element_text(size = 10, face = "italic", colour = "black"),
          axis.text.x = element_text(face = "bold", colour = "black"))
  if (cell_type != "ILC1") heatmaps[[cell_type]] <- heatmaps[[cell_type]] + theme(axis.text.y = element_blank())
}
panel_C <- arrangeGrob(grobs = heatmaps, ncol = 3, widths = c(1.3, 1, 1),
                       top = "C   Mean scaled gene expression")
ggsave(file.path(fig2_dir, "Figure_2C.png"), panel_C, width = 8, height = 6.5, dpi = 300, bg = "white")


##################################
########## Figure 2D #############
##################################

boxplot_genes <- c("IFITM2", "FAM65B", "SELL")
percentage_rows <- list()
for (cell_type in ilc_types) {
  samples <- unique(metadata2$Sample_id[metadata2$handgated_cell_type2 %in% cell_type])
  for (sample_id in samples) {
    md <- metadata2[metadata2$Sample_id == sample_id & metadata2$handgated_cell_type2 %in% cell_type, ]
    if (nrow(md) > 7) {
      percentage_rows[[paste(cell_type, sample_id)]] <- data.frame(
        Sample_id = sample_id, PPID = md$PPID[1], ILCn = cell_type,
        COVID_status = md$COVID_status[1], n_cells = nrow(md), Gene = boxplot_genes,
        Percent = rowMeans(normalized_rna[boxplot_genes, rownames(md), drop = FALSE] > 0) * 100)
    }
  }
}
percentages2 <- do.call(rbind, percentage_rows)
percentages2$ILCn <- factor(percentages2$ILCn, levels = ilc_types)
percentages2$COVID_status <- factor(percentages2$COVID_status, levels = c("healthy", "COVID"))
percentages2$longterm <- factor(ifelse(percentages2$PPID %in% long_ter_covd, "Long", "Other"), levels = c("Long", "Other"))
boxplots2 <- list()
for (gene in boxplot_genes) {
  d <- percentages2[percentages2$Gene == gene, ]
  labels <- all_DE[all_DE$Gene == gene, ]
  labels$ILCn <- factor(labels$ILCn, levels = ilc_types)
  labels$label <- paste0("q = ", format(labels$p_val_adj, scientific = TRUE, digits = 2))
  labels$y <- max(d$Percent) * 1.10
  boxplots2[[gene]] <- ggplot(d, aes(COVID_status, Percent)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white") +
    geom_point(aes(fill = COVID_status, shape = longterm), colour = "black", size = 2.4,
               position = position_dodge(width = 0.01), alpha = 0.85) +
    scale_fill_manual(values = c(healthy = "#008000", COVID = "#F0B548")) +
    scale_shape_manual(values = c(Long = 24, Other = 21)) +
    facet_wrap(~ILCn, nrow = 1, scales = "free") +
    geom_text(data = labels, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 3.2) +
    labs(title = paste0(gene, "+ cells"), x = NULL, y = "% of ILC subtype") +
    theme_bw(base_size = 11) + theme(legend.position = "none", axis.text.x = element_blank(),
      strip.background = element_rect(fill = "white"), plot.title = element_text(hjust = 0.5, face = "italic"))
}
panel_D <- arrangeGrob(grobs = boxplots2, ncol = 1, top = "D")
ggsave(file.path(fig2_dir, "Figure_2D.png"), panel_D, width = 5.5, height = 9.5, dpi = 300, bg = "white")

#write.csv(percentages2, file.path(fig2_dir, "Figure_2D_sample_percentages.csv"), row.names = FALSE)
#write.csv(all_DE[all_DE$Gene %in% boxplot_genes, c("ILCn", "Gene", "p_val", "p_val_adj")],
#          file.path(fig2_dir, "Figure_2D_pvalues.csv"), row.names = FALSE)




##################################
########## Figure 2E #############
##################################

raw_counts <- t(as.matrix(GetAssayData(rna, assay = "RNA", layer = "counts")))
original_status <- metadata2$COVID_status_original
correlations <- list()
for (status in c("healthy", "COVID")) {
  group_counts <- raw_counts[original_status == status, , drop = FALSE]
  group_counts <- group_counts[, colSums(group_counts > 0) > 10, drop = FALSE]
  network_expression <- sclink_norm(group_counts, scale.factor = 1e6,
                                     filter.genes = TRUE, n = ncol(group_counts))
  correlations[[status]] <- sclink_cor(network_expression, ncores = 1, nthre = 20, dthre = 0.9)
  # write.csv(correlations[[status]], file.path(fig2_dir, paste0("scLink_correlations_", status, ".csv")))
}
cor_healthy <- correlations$healthy
cor_covid <- correlations$COVID
cor_healthy[abs(cor_healthy) < 0.7] <- 0
cor_covid[abs(cor_covid) < 0.7] <- 0
common_genes <- intersect(colnames(cor_healthy), colnames(cor_covid))
healthy_only <- setdiff(colnames(cor_healthy), colnames(cor_covid))
covid_only <- setdiff(colnames(cor_covid), colnames(cor_healthy))

matrix_edges <- function(m, threshold) {
  ij <- which(m > threshold, arr.ind = TRUE)
  data.frame(from = rownames(m)[ij[, 1]], to = colnames(m)[ij[, 2]], weight = m[ij])
}
edges <- rbind(
  matrix_edges(cor_healthy[common_genes, common_genes] - cor_covid[common_genes, common_genes], 0.2),
  matrix_edges(cor_healthy[healthy_only, , drop = FALSE], 0.699),
  matrix_edges(cor_covid[covid_only, , drop = FALSE], 0.699))
focus_genes <- c("SELL", "FAM65B", "IFITM2")
edges <- edges[edges$from %in% focus_genes | edges$to %in% focus_genes, ]
network <- igraph::graph_from_data_frame(edges, directed = FALSE)
network <- igraph::simplify(network, remove.multiple = TRUE, remove.loops = TRUE,
                     edge.attr.comb = list(weight = "sum"))
edge_table <- igraph::as_data_frame(network, what = "edges")
#write.csv(edge_table, file.path(fig2_dir, "Figure_2E_network_edges.csv"), row.names = FALSE)

sector_order <- c("IGHG3-secreted", "IFITM2", "CD40", "POU2AF1", "SELL", "NINJ2", "IL12A",
                   "CR2", "CCR2", "BCL6", "FAM65B", "IL15", "FBXO22", "TCF4", "ADGRE1",
                   "TNFRSF1B", "BCL2A1", "PECAM1", "IL6R", "CXCL16", "CD33", "VNN2",
                   "RNASE6", "CHI3L2", "TNFRSF17", "OAS1")
sector_order <- intersect(sector_order, unique(c(edge_table$from, edge_table$to)))
swap <- !edge_table$from %in% focus_genes & edge_table$to %in% focus_genes
edge_table[swap, c("from", "to")] <- edge_table[swap, c("to", "from")]
network_colours <- setNames(colorRampPalette(brewer.pal(12, "Set3"))(length(sector_order)), sector_order)
network_colours[c("FAM65B", "SELL", "IFITM2")] <- c("#5E4FA2", "#8DD3C7", "#FDAE61")


png(file.path(fig2_dir, "Figure_2E.png"), width = 2400, height = 2400, res = 300)
circos.clear()
par(mar = c(1, 1, 2, 1))
chordDiagram(edge_table[, c("from", "to")], order = sector_order, grid.col = network_colours, col = network_colours[edge_table$from],
              transparency = 0, annotationTrack = "grid", preAllocateTracks = 1)
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  sector <- get.cell.meta.data("sector.index")
  xlim <- get.cell.meta.data("xlim"); ylim <- get.cell.meta.data("ylim")
  circos.text(mean(xlim), ylim[1] + 0.1, sector, facing = "clockwise", niceFacing = TRUE,
              adj = c(0, 0.5), cex = 0.65,
              col = ifelse(sector %in% focus_genes, "#B3262D", "black"))
}, bg.border = NA)
title("E   Differential co-expression network", cex.main = 1)
circos.clear(); dev.off()


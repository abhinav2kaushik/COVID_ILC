### WARNING: COLORS of Network layout may differ ###

library(Seurat)
library(ggplot2)
library(gridExtra)
library(igraph)
library(circlize)
setwd("~/Library/CloudStorage/OneDrive-HarvardUniversity/ILC/Analysis/Result_replication")
.libPaths(c(
  "./src/",
  .libPaths()
))
library(scLink)
library(MAST)


fig3_dir <- "Figure_3"
dir.create(fig3_dir, showWarnings = FALSE, recursive = TRUE)
ilc_types <- c("ILC1", "ILC2", "ILCp")
healthy <- c("1033", "1348", "1804", "1987", "2195", "55689-0068", "55689-0071",
             "55689-0080", "WF001", "WF002", "WF006", "WF008", "WF010", "WF015",
             "WF016", "WF018", "WF019", "WF021", "WF025", "WF045", "WF026",
             "WF031", "WF033", "WF035", "WF042")
long_ter_covd <- c("55689-0007", "55689-0012", "55689-0015", "55689-0029",
                   "55689-0031", "55689-0034", "55689-0044", "55689-0050",
                   "55689-0057", "55689-0095", "55689-0165", "55689-0219",
                   "55689-0220", "OHC-001", "OHC-010", "OHC001", "OHC010")

covid_ilc = readRDS("covid_ilc_counts.rds")
rna <- covid_ilc
DefaultAssay(rna) <- "RNA"
rna <- NormalizeData(rna, normalization.method = "LogNormalize", scale.factor = 10000,
                     verbose = FALSE)

rna$COVID_status <- ifelse(rna$PPID %in% healthy, "healthy", "COVID")
metadata2 <- rna[[]]
normalized_rna <- GetAssayData(rna, assay = "RNA", layer = "data")


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
  #write.csv(de, file.path(fig3_dir, paste0("DE_", cell_type, ".csv")), row.names = FALSE)
}
all_DE <- do.call(rbind, DE_RNA)

########## ########## ########## 
########## Figure 3A-C #########
########## ########## ########## 
plot_genes <- list(
  ILC1 = c("CD6", "STAT5A", "SPOCK2", "PIK3AP1", "ITGB7", "RGS1"),
  ILC2 = c("STAT3", "ICOS", "SPOCK2", "IL4R", "PIK3AP1", "HPGDS", "ITGB7", "KLRG1"),
  ILCp = c("TNFRSF25", "STAT5A", "IL2RA", "IL4R", "CCR7", "SPOCK2", "ICAM1", "PIK3AP1"))
percentage_rows <- list()
for (cell_type in ilc_types) {
  samples <- unique(metadata2$Sample_id[metadata2$handgated_cell_type2 %in% cell_type])
  for (sample_id in samples) {
    md <- metadata2[metadata2$Sample_id == sample_id & metadata2$handgated_cell_type2 %in% cell_type, ]
    if (nrow(md) > 7) {
      genes <- plot_genes[[cell_type]]
      percentage_rows[[paste(cell_type, sample_id)]] <- data.frame(
        Sample_id = sample_id, PPID = md$PPID[1], ILCn = cell_type,
        COVID_status = md$COVID_status[1], n_cells = nrow(md), Gene = genes,
        Percent = Matrix::rowMeans(normalized_rna[genes, rownames(md), drop = FALSE] > 0) * 100)
    }
  }
}
percentages3 <- do.call(rbind, percentage_rows)
percentages3$COVID_status <- factor(percentages3$COVID_status, levels = c("healthy", "COVID"))
percentages3$longterm <- factor(ifelse(percentages3$PPID %in% long_ter_covd, "Long", "Other"), levels = c("Long", "Other"))
percentages3$Group <- factor(ifelse(percentages3$COVID_status == "healthy", "Healthy",
  ifelse(percentages3$longterm == "Long", "Long term COVID-19+", "COVID-19+")),
  levels = c("Healthy", "COVID-19+", "Long term COVID-19+"))
plots3 <- list()
pvalues3 <- list()
for (i in seq_along(ilc_types)) {
  cell_type <- ilc_types[i]
  d <- percentages3[percentages3$ILCn == cell_type, ]
  d$Gene <- factor(d$Gene, levels = plot_genes[[cell_type]])
  labels <- DE_RNA[[cell_type]][plot_genes[[cell_type]], ]
  labels$Gene <- factor(labels$Gene, levels = plot_genes[[cell_type]])
  labels$y <- sapply(as.character(labels$Gene), function(g) max(d$Percent[d$Gene == g]) * 1.05)
  labels$label <- paste0("q = ", format(labels$p_val_adj, scientific = TRUE, digits = 2))
  pvalues3[[cell_type]] <- labels[, c("ILCn", "Gene", "p_val", "p_val_adj")]
  plots3[[cell_type]] <- ggplot(d, aes(COVID_status, Percent)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white") +
    geom_point(aes(fill = Group, shape = Group), colour = "black",
               alpha = 0.7, size = 2.6, position = position_dodge(width = 0.01)) +
    facet_wrap(~Gene, scales = "free", nrow = 1) +
    scale_fill_manual(values = c("Healthy" = "#008000", "COVID-19+" = "#F0B548", "Long term COVID-19+" = "#F0B548")) +
    scale_shape_manual(values = c("Healthy" = 21, "COVID-19+" = 21, "Long term COVID-19+" = 24)) +
    geom_text(data = labels, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = 3.1) +
    labs(title = paste(LETTERS[i], cell_type), x = NULL, y = "% gene-positive cells", fill = NULL, shape = NULL) +
    theme_bw(base_size = 12) +
    theme(strip.background = element_rect(fill = "white"), strip.text = element_text(face = "italic"),
          axis.text.x = element_blank(), legend.position = "bottom", plot.title = element_text(face = "bold"))
  ggsave(file.path(fig3_dir, paste0("Figure_3", LETTERS[i], ".png")), plots3[[cell_type]],
         width = length(plot_genes[[cell_type]]) * 1.55, height = 4.5, dpi = 300, bg = "white")
}
#write.csv(percentages3, file.path(fig3_dir, "Figure_3ABC_sample_percentages.csv"), row.names = FALSE)
#write.csv(do.call(rbind, pvalues3), file.path(fig3_dir, "Figure_3ABC_pvalues.csv"), row.names = FALSE)

########## ########## #########
########## Figure 3D. #########
########## ########## #########
focus_genes <- list(
  ILC1 = c("ICAM1","CD6","TGFB1","RGS1","STAT5A","PIK3AP1","ITGB7","IFITM2","FAM65B","SELL"),
  ILC2 = c("ICAM1","NAMPT","IL4R","ICOS","IL32","PIK3AP1","ITGB7","HPGDS","ITGB2","IFITM3","IFITM2","FAM65B","SELL"),
  ILCp = c("ICAM1","ICOS","IL4R","RGS1","CCR7","IL2RA","LTA","RORA","TRDC","CD300A","FCER1G","IFITM2","FAM65B","SELL"))
raw_counts <- t(as.matrix(GetAssayData(rna, assay = "RNA", layer = "counts")))
matrix_edges <- function(m, threshold) {
  ij <- which(m > threshold, arr.ind = TRUE)
  data.frame(from = rownames(m)[ij[, 1]], to = colnames(m)[ij[, 2]], weight = m[ij])
}
networks3 <- list()
for (cell_type in ilc_types) {
  correlations <- list()
  for (status in c("healthy", "COVID")) {
    message("Figure 3D: ", cell_type, " / ", status)
    use <- metadata2$handgated_cell_type2 %in% cell_type & metadata2$COVID_status == status
    group_counts <- raw_counts[use, , drop = FALSE]
    group_counts <- group_counts[, colSums(group_counts > 0) > 10, drop = FALSE]
    expression <- sclink_norm(group_counts, scale.factor = 1e6, filter.genes = TRUE, n = ncol(group_counts))
    correlations[[status]] <- sclink_cor(expression, ncores = 1, nthre = 20, dthre = 0.9)
    #write.csv(correlations[[status]], file.path(fig3_dir, paste0("Correlation_", cell_type, "_", status, ".csv")))
  }
  h <- correlations$healthy
  c <- correlations$COVID
  h[abs(h) < 0.5] <- 0
  c[abs(c) < 0.5] <- 0
  common <- intersect(colnames(h), colnames(c))
  edges <- rbind(
    matrix_edges(h[common, common, drop = FALSE] - c[common, common, drop = FALSE], 0.2),
    matrix_edges(h[setdiff(colnames(h), colnames(c)), , drop = FALSE], 0.499),
    matrix_edges(c[setdiff(colnames(c), colnames(h)), , drop = FALSE], 0.499))
  edges <- edges[edges$from %in% focus_genes[[cell_type]] | edges$to %in% focus_genes[[cell_type]], ]
  network <- graph_from_data_frame(edges, directed = FALSE)
  network <- simplify(network, remove.multiple = TRUE, remove.loops = TRUE, edge.attr.comb = list(weight = "sum"))
  networks3[[cell_type]] <- igraph::as_data_frame(network, what = "edges")
  #write.csv(networks3[[cell_type]], file.path(fig3_dir, paste0("Network_edges_", cell_type, ".csv")), row.names = FALSE)
}

network_colours <- c(PIK3AP1 = "#648D85", ICAM1 = "#F1CFEC", STAT5A = "#EF9AC5", 
FAM65B = "#6DE659", TGFB1 = "#E1BAED", ITGB7 = "#6FCB8B", TIAF1 = "#E7902E", 
DPP4 = "#BEAAA7", FYN = "#A2F2F2", C10orf54 = "#EDB595", TSPAN32 = "#8582A2", 
LIPA = "#95D8C8", TCF7 = "#EE51A5", VPS28 = "#E14071", RORA = "#66BEB0", 
ATF6B = "#F28FDF", ITK = "#E7DE42", CD247 = "#47722A", STAT3 = "#9AC2C7", 
FOXO1 = "#89BA70", CTSD = "#5852E8", RGS1 = "#A7C258", IKZF1 = "#E0F3D4", 
RPN2 = "#98EECD", CASP5 = "#EE34EC", FYB = "#3FC77B", TNFRSF25 = "#ADA268", 
BCL11B = "#BCCDAA", ITGAM = "#5635ED", ITGAE = "#4BD6F7", PTPRC = "#F1F4A8", 
CASP3 = "#EBDF88", KLRK1 = "#E9BC6A", STAT4 = "#E5F2EC", ZBTB16 = "#BA7BF0", 
IL23R = "#ABEE3C", BAX = "#DC4837", CD6 = "#B28CEA", CD3D = "#EFAFEA", 
LGALS9 = "#B62AAA", PDIA4 = "#E8E97A", GNAI2 = "#B6C1B9", CD44 = "#8EC4E8", 
`TARP-refseq` = "#776DE9", ITGA4 = "#7D4B4B", IGLC3 = "#E3F42F", 
IFITM2 = "#C17592", IL2RA = "#5542B2", XBP1 = "#734D8E", PCNA = "#9AED9D", 
PDIA6 = "#4C86D5", IL18R1 = "#ED30C0", GIMAP2 = "#BDBE69", IL4R = "#5CA583", 
APOBEC3G = "#CEB68D", KLRG1 = "#C4C8F1", FTH1 = "#C46CC7", TRBC2 = "#E159C4", 
TRAC = "#EB7A93", `IGHA1-secreted` = "#52B1E7", LTA = "#C2E8EF", 
GATA3 = "#EA67EE", HMGB2 = "#9ACE43", NAMPT = "#4C6FA7", ICOS = "#E18767", 
SELL = "#ECD0CE", HPGDS = "#B1CBE7", GAB2 = "#F8F1C3", LAP3 = "#4EEF41", 
TNFRSF1B = "#CCC7D3", TXK = "#8FB78D", IL2 = "#54F6E1", TRAT1 = "#96A5B0", 
TNFRSF4 = "#F7D598", BIRC3 = "#59BC47", ADA = "#96BDF3", SELPLG = "#298BAB", 
NCR3 = "#DFA59D", IFNGR1 = "#E0DBA4", GIMAP5 = "#852DBD", CD3E = "#E9E4EF", 
IGBP1 = "#BC7CBC", CNOT2 = "#B8D391", TRIB2 = "#87A7F3", STAT6 = "#D0F4B2", 
CD28 = "#EAB0C5", FOXO3 = "#52F3C8", CD79B = "#BB4CF2", CD300A = "#977B5F", 
LAT2 = "#731C5E", IL2RB = "#49BADC", MYC = "#A3ECB4", CCND2 = "#CDF594", 
ZAP70 = "#AA889D", CCR6 = "#5DB9C4", BCL2 = "#BFB3CF", DOCK8 = "#AAF182", 
LAMP1 = "#5DE0E9", FOSL1 = "#7CE7B5", VEGFA = "#7B73CB", FCER1G = "#B9AEE8", 
RUNX3 = "#B956D1", YBX3 = "#C1F0D9", CD63 = "#C45C97", LYN = "#4CEE89", 
NKG7 = "#71A1BD", `HLA-DPA1` = "#B191CE", KIT = "#3E76E7", IL18RAP = "#E487F2", 
KLRC1 = "#4BBC9A", CD7 = "#D2A4F2", LTB = "#D2ED61", IL7R = "#4FD0C4", 
CD69 = "#BA8B37", LAT = "#8197D3", PTGDR2 = "#6E89EC", TRDC = "#E9D5BE", 
LGALS1 = "#DFBB3F", ANXA5 = "#61EDEA", KDELR1 = "#5DEDA5", BTG1 = "#9ADDF2"
)

plot_network3 <- function(cell_type) {
  edges <- networks3[[cell_type]]
  circos.clear()
  par(mar = c(1, 1, 2, 1))
  chordDiagram(edges[, c("from", "to")], grid.col = network_colours,
               transparency = 0, annotationTrack = "grid", preAllocateTracks = 1)
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    sector <- get.cell.meta.data("sector.index")
    circos.text(mean(xlim), ylim[1] + 0.1, sector, facing = "clockwise",
                niceFacing = TRUE, adj = c(0, 0.5), cex = 0.55)
  }, bg.border = NA)
  title(cell_type)
  circos.clear()
}
for (cell_type in ilc_types) {
  png(file.path(fig3_dir, paste0("Figure_3D_", cell_type, ".png")), width = 2400, height = 2400, res = 300)
  plot_network3(cell_type)
  dev.off()
}
png(file.path(fig3_dir, "Figure_3D.png"), width = 5400, height = 1900, res = 300)
par(mfrow = c(1, 3))
for (cell_type in ilc_types) plot_network3(cell_type)
dev.off()
capture.output(sessionInfo(), file = file.path(fig3_dir, "sessionInfo.txt"))

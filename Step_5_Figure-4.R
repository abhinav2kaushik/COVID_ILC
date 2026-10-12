### WARNING: COLORS of Network layout may differ ###
library(Seurat)
library(ggplot2)
library(gridExtra)
library(pheatmap)
library(circlize)
library(dplyr)
setwd("~/Library/CloudStorage/OneDrive-HarvardUniversity/ILC/Analysis/Result_replication")
.libPaths(c(
  "./src/",
  .libPaths()
))
library(MAST)


fig4_dir <- "Figure_4"
scenic_dir <- "PyScenic_output/"
heatmap_status <- "corrected"  # "original" uses COVID_status_original, with the old OHC001 correction.
dir.create(fig4_dir, showWarnings = FALSE, recursive = TRUE)
ilc_types <- c("ILC1", "ILC2", "ILCp")
covid_ilc <- readRDS("covid_ilc_counts.rds")


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
  #write.csv(de, file.path(fig4_dir, paste0("DE_", cell_type, ".csv")), row.names = FALSE)
}
all_DE <- do.call(rbind, DE_RNA)



########## ########## ########## 
##########  Figure 4A ##########
########## ########## ########## 

binary_auc <- as.matrix(read.csv(file.path(scenic_dir, "BIN_MTX.csv"), row.names = 1, check.names = FALSE))
colnames(binary_auc) <- sub("\\(\\+\\)$", "", colnames(binary_auc))
stopifnot(all(rownames(binary_auc) %in% rownames(metadata2)), all(binary_auc %in% c(0, 1)))
ann <- metadata2[rownames(binary_auc), c("handgated_cell_type2", "COVID_status")]
names(ann) <- c("ILC", "status")
if (heatmap_status == "original") {
  ann$status <- metadata2[rownames(binary_auc), "COVID_status_original"]
  ann$status[metadata2[rownames(binary_auc), "PPID"] == "OHC001"] <- "COVID"
}
ann$ILC <- factor(ann$ILC, levels = ilc_types)
ann$status <- factor(ann$status, levels = c("healthy", "COVID"))
stopifnot(!anyNA(ann))
ann <- ann[order(ann$ILC, ann$status), , drop = FALSE]
groups <- paste(ann$ILC, ann$status)
gaps <- head(cumsum(rle(groups)$lengths), -1)
heat_matrix <- t(binary_auc[rownames(ann), , drop = FALSE])
row_tree <- as.hclust(rev(as.dendrogram(hclust(dist(heat_matrix), method = "complete"))))

pheatmap(heat_matrix,
  cluster_rows = row_tree, cluster_cols = FALSE, clustering_distance_rows = "euclidean",
  clustering_method = "complete", annotation_col = ann, gaps_col = gaps,
  color = c("white", "black"), breaks = c(-0.5, 0.5, 1.5), border_color = NA,
  show_colnames = FALSE, legend = FALSE, fontsize_row = 10,
  annotation_colors = list(ILC = c(ILC1 = "green", ILC2 = "brown", ILCp = "blue"),
                           status = c(healthy = "#008000", COVID = "#F0B548")),
  main = "A   Binarized AUCell", filename = file.path(fig4_dir, "Figure_4A.png"), width = 9, height = 8)
#write.csv(data.frame(cell_id = rownames(ann), ann), file.path(fig4_dir, "Figure_4A_cell_order.csv"), row.names = FALSE)

########## ########## ########## 
##########  Figure 4B ##########
########## ########## ########## 

DefaultAssay(covid_ilc) <- "RNA"

# Create the normalized RNA data layer
covid_ilc <- NormalizeData(
  covid_ilc,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

Idents(covid_ilc) <- "handgated_cell_type2"
covid_ilc <- SetIdent(covid_ilc, value = "handgated_cell_type2") ## by default groups will be the adt clusters, we set group to be COVID status
levels(covid_ilc) ## check identity class

DefaultAssay(covid_ilc) <- 'RNA'
ExpressMat1 = t(GetAssayData(object = covid_ilc, layer = "data")) ## rows are cell IDs and columns are genes from scRNAseq
ExpressMat1 = as.data.frame(ExpressMat1)  ## ExpressMat1 is scaled rnaseq 
mdf = covid_ilc@meta.data ## meta-data
mdf$COVID_status = ifelse(mdf$PPID %in% healthy,  "healthy","COVID")

DEgenes = list() ## singLE cell genes DE (HANDGATED)
for (clus in levels(covid_ilc)){
  DEmarkers <- FindMarkers(covid_ilc, 
                           group.by = 'COVID_status', # Regroup cells into a different identity class prior to performing differential expression
                           ident.1 = "COVID", 
                           ident.2 = 'healthy', 
                           subset.ident = clus, ## subset a particular identity class (seurat_cluster; see above) prior to regrouping.
                           latent.vars = 'batch',
                           test.use='MAST',
                           assay='RNA',
                           logfc.threshold = 0.0
  )
  DEmarkers$cluster = clus
  DEgenes[[clus]] = DEmarkers
}

d1 = DEgenes$ILC1
d2 = DEgenes$ILC2
d3 = DEgenes$ILCp

d1$genes = rownames(d1)
d2$genes = rownames(d2)
d3$genes = rownames(d3)

adjPthresh = 0.05
logFCthresh = 0.5
target_genes = c()
## ILC1 ##
covid_ilc1 <- subset(x = covid_ilc, subset = handgated_cell_type2 == "ILC1")
degenesilc1 = d1$genes[d1$p_val_adj < adjPthresh & abs(d1$avg_log2FC)> logFCthresh ]
target_genes = c(target_genes,degenesilc1) ## these genes will be plotted in heatmap 
## ILC2 ##
covid_ilc2 <- subset(x = covid_ilc, subset = handgated_cell_type2 == "ILC2")
degenesilc2 = d2$genes[d2$p_val_adj < adjPthresh & abs(d2$avg_log2FC)> logFCthresh ]
target_genes = c(target_genes,degenesilc2) ## these genes will be plotted in heatmap 
## ILC3 ##
covid_ilc3 <- subset(x = covid_ilc, subset = handgated_cell_type2 == "ILCp")
degenesilc3 = d3$genes[d3$p_val_adj < adjPthresh & abs(d3$avg_log2FC)> logFCthresh ]
target_genes = unique(c(target_genes,degenesilc3)) ## these genes will be plotted in heatmap 

TFs = c("ATF6B","BACH2","BCL6", "EGR1","EOMES","FOSB","FOSL1","FOXO1","FOXO3","FOXP1",
        "GATA3","IKZF1","IKZF2","IRF4","IRF8","JUN","JUNB","LEF1","MYC","RORC",
        "RUNX3","STAT1","STAT3","STAT6","TBX21")

targetGnes = list()
targetGnes[["ILC1"]] = TFs
targetGnes[["ILC2"]] = TFs
targetGnes[["ILCp"]] = TFs
ggdf = data.frame()
for (clus in c("ILC1","ILC2","ILCp")){
  tmpdf = DEgenes[[clus]]
  message(clus)
  zero = dplyr::filter(mdf, handgated_cell_type2==clus)
  SIDs = unique(zero$Sample_id)
  tgens = targetGnes[[clus]]
  for (m in tgens){
    for (samID in SIDs){
      ty6 = dplyr::filter(mdf, Sample_id == samID & handgated_cell_type2 == clus)
      cellIDs = ty6$cell_id ## ILC for this specific ILCn
      if (length(cellIDs) > 7) { ## if total number of ILCn is > 3 in a given sample
        tmp2 = as.data.frame(ExpressMat1[cellIDs,m])
        vals = as.numeric(tmp2[,1]) 
        vals[vals < 0] = 0
        nsum = sum(vals > 0) ## num of cells expressing this marker
        if (nsum > 0){
          total_expression = (nsum / length(cellIDs))*100 ## %age of cells per sample expressing this marker 
        } else {
          total_expression = 0
        }
        tmpodf = data.frame(value = total_expression,
                            sampleID  = samID,
                            PPID = unique(ty6$PPID),
                            Gene = m,
                            ILCn = clus,
                            COVID_status = unique(ty6$COVID_status))
        ggdf = rbind(ggdf,tmpodf)
      }
    }
  }
}
ggdf$longterm  = ifelse(ggdf$PPID %in% long_ter_covd, "Long","Other")
ggdf$COVID_status = factor(ggdf$COVID_status, levels = c("healthy","COVID"))

clrs = c("#008000","#F0B548")
plots = list()
for (clus in c("ILC1","ILC2","ILCp")){
  ggdf0 = dplyr::filter(ggdf, ILCn == clus)
  
  textdf0 = ggdf0[,c("value","Gene")] %>% dplyr::group_by(Gene) %>% dplyr::summarise_all(max)
  textdf0$max = textdf0$value
  textdf0$x = 1.5
  textdf0$Y = textdf0$max + (0.05 * textdf0$max)
  
  da1 = DEgenes[[clus]]
  da1$genes = rownames(da1)
  da1 = dplyr::filter(da1, genes %in% textdf0$Gene)
  rownames(da1) = da1$genes
  da1 = da1[textdf0$Gene,]
  
  keygenes = dplyr::filter(da1, p_val_adj < 0.05)$genes
  allgenes =  c("LEF1","RUNX3","STAT3","STAT6","IKZF1","GATA3","JUNB")
  com = intersect(keygenes,allgenes)
  other = setdiff(allgenes,keygenes)
  
  if (length(keygenes) < 7){
    diff = 7 - length(keygenes)
    keygenes = c(keygenes, other[1:diff])
  }
  
  label = paste("q =", format(da1$p_val_adj, scientific =T, digits=2))
  textdf0$label=label
  
  ggdf0$Gene = factor(ggdf0$Gene, levels = targetGnes[[clus]] )
  textdf0$Gene = factor(textdf0$Gene, levels = targetGnes[[clus]])
  
  ggdf0 = dplyr::filter(ggdf0, Gene %in% keygenes)
  textdf0 = dplyr::filter(textdf0, Gene %in% keygenes)
  
  g2 <- ggplot(ggdf0, aes(y=value, x=COVID_status)) +
    geom_boxplot(width = 0.5, alpha = 0.4,outlier.color = NA) +
    geom_point(aes(fill=COVID_status,shape = longterm), color = "black", alpha = 0.7, position = position_dodge(width = 0.01), size=2.6) + 
    theme_bw(base_size = 14) + 
    facet_wrap(~Gene, scales = "free", nrow = 1) +
    scale_fill_manual(values = clrs) +
    scale_shape_manual(values = c(24,21)) +
    theme(legend.position = "none",
          strip.background = element_rect(fill = "white"),
          strip.text = element_text(size=10, face = "italic"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
    ) + ylab(label = "% cells expression") +
    geom_text(data=textdf0, mapping=aes(x=x, y=Y, label=label), size=3.1)
  
  plots[[clus]] = g2
}


fn = paste0(fig4_dir, "/Figure4B.png")

png(fn, height =8, width = 10, units = "in", res = 300)
do.call("grid.arrange", plots)
dev.off()


########## ########## ########## 
##########  Figure 4C ##########
########## ########## ########## 
TFs <- c("ATF6B","BACH2","BCL6","EGR1","EOMES","FOSB","FOSL1","FOXO1","FOXO3","FOXP1",
         "GATA3","IKZF1","IKZF2","IRF4","IRF8","JUN","JUNB","LEF1","MYC","RORC",
         "RUNX3","STAT1","STAT3","STAT6","TBX21")
adjacencies <- read.csv(file.path(scenic_dir, "adj.csv"), check.names = FALSE)
# Saved original palette is embedded below; additional genes in fresh runs get fixed colours.
original_colours <- c(FOSB = "#E699E7", JUN = "#56E9A8", JUNB = "#A57EE4", LEF1 = "#ED40AB", 
EGR1 = "#E0E23A", RUNX3 = "#85CBE3", IRF8 = "#A9A494", DUSP1 = "#ACEFDB", 
CD37 = "#E4E6CC", `HLA-A` = "#BD53A2", DUSP2 = "#DCB3E6", BTG1 = "#69A9DE", 
C10orf54 = "#73EF84", LAMP1 = "#68C5AB", IL32 = "#E98A8A", IFITM2 = "#E389B2", 
IL4R = "#DFAD55", CD52 = "#E0623A", IFITM3 = "#86E94A", LTB = "#73BC66", 
CD69 = "#DF3F6D", VPS28 = "#9D47CA", GNAI2 = "#E1AEB6", CXCR4 = "#A5ADCD", 
CTSW = "#E0B893", TNFRSF25 = "#C3E3E6", SPOCK2 = "#939EEB", TRBC2 = "#9C7497", 
S100A10 = "#DC40E6", PIK3IP1 = "#5A65E3", ITGB2 = "#ABE9AF", 
ITGB7 = "#E3D2E3", NCR3 = "#E66FDC", RGS1 = "#E8EAA7", LCK = "#4A61A4", 
ICOS = "#7F37E6", FYB = "#DBE777", CD7 = "#558D8F", ZAP70 = "#64EBE8", 
TSPAN32 = "#B97A5A", CCR7 = "#929F56")

all_genes <- sort(unique(c(adjacencies$TF, adjacencies$target)))
network_colours <- setNames(grDevices::hcl.colors(length(all_genes), "Dynamic"), all_genes)
network_colours[names(original_colours)] <- original_colours
networks4 <- list()
for (cell_type in ilc_types) {
  de <- DE_RNA[[cell_type]]
  targets <- de$Gene[de$p_val_adj < 0.05]
  networks4[[cell_type]] <- adjacencies[adjacencies$TF %in% TFs & adjacencies$target %in% targets &
                                        adjacencies$importance > 50, c("TF", "target", "importance")]
}
plot_network4 <- function(cell_type) {
  edges <- networks4[[cell_type]]
  tfs <- unique(edges$TF)
  targets <- setdiff(unique(edges$target), tfs)
  circos.clear()
  par(mar = c(1, 1, 2, 1))
  chordDiagram(edges, order = c(tfs, targets), grid.col = network_colours,
    transparency = 0, annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.04, track.margin = c(0, 0)))
  highlight.sector(tfs, track.index = 1, col = "#A52A2A")
  highlight.sector(targets, track.index = 1, col = "#1E90FF")
  circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    circos.text(mean(xlim), ylim[1] + 1.2, get.cell.meta.data("sector.index"),
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
  }, bg.border = NA)
  title(cell_type)
  circos.clear()
}
for (cell_type in ilc_types) {
  png(file.path(fig4_dir, paste0("Figure_4C_", cell_type, ".png")), width = 2400, height = 2400, res = 300)
  plot_network4(cell_type)
  dev.off()
}
capture.output(sessionInfo(), file = file.path(fig4_dir, "sessionInfo.txt"))

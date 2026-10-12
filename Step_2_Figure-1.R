library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(Seurat)
library(harmony)
library(scales)
library(gridExtra)
library(grid)
library(MASS)
library(png)

setwd("~/Library/CloudStorage/OneDrive-HarvardUniversity/ILC/Analysis/Result_replication")

################################
########## Figure 1B ###########
################################

object_file <- "covid_ilc_counts.rds"
output_dir <- "Figure_1"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
covid_ilc <- readRDS(object_file) ## load object we created in step_1
DefaultAssay(covid_ilc) <- "Abseq"
covid_ilc <- NormalizeData(covid_ilc, normalization.method = "CLR", margin = 2)
ilc <- subset(covid_ilc, subset = handgated_cell_type2 %in% c("ILC1", "ILC2", "ILCp"))

gating <- data.frame(
  cell_id = colnames(covid_ilc),
  handgated_cell_type2 = covid_ilc$handgated_cell_type2
)
gating <- gating[!is.na(gating$handgated_cell_type2), ]

ilc <- ilc[, gating$cell_id]  # to Preserve the original UMAP input cell order.

ilc <- RunUMAP(ilc, assay = "Abseq", slot = "data",
               features = c("CD127", "CD161", "CD117", "CD45", "CD294"),
               umap.method = "uwot", n.neighbors = 30, n.components = 2,
               metric = "cosine", min.dist = 0.3, spread = 1,
               learning.rate = 1, seed.use = 42,
               reduction.name = "ilc.umap", reduction.key = "ILCUMAP_")


new_coordinates <- Embeddings(ilc, "ilc.umap")
umap <- data.frame(
  cell_id = rownames(new_coordinates),
  UMAP1 = -new_coordinates[, 1],
  UMAP2 = new_coordinates[, 2]
)
ilc <- ilc[, umap$cell_id]  # Keep the original figure's point order.
umap$UMAP1 <- umap$UMAP1
umap$ILC <- factor(ilc$handgated_cell_type2, levels = c("ILC1", "ILC2", "ILCp"))
metadata <- ilc[[]]
umap$COVID_status <- metadata$COVID_status[match(umap$cell_id, rownames(metadata))]

expression <- t(as.matrix(GetAssayData(ilc, assay = "Abseq", layer = "data")))
umap <- cbind(umap, expression[umap$cell_id, c("CD127", "CD161", "CD117", "CD294")])

x_limits <- range(umap$UMAP1)
y_limits <- range(umap$UMAP2)
plot_theme <- theme_bw(base_size = 12) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        panel.border = element_rect(colour = "grey75"),
        panel.grid = element_line(colour = "grey92"),
        plot.title = element_text(face = "bold", size = 13),
        legend.key.height = unit(0.30, "in"))

p_ilc <- ggplot(umap, aes(UMAP1, UMAP2, colour = ILC)) +
  geom_point(size = 0.8, alpha = 0.8) +
  scale_colour_manual(values = c(ILC1 = "#85BF52", ILC2 = "#FFA500", ILCp = "#4C5E9C")) +
  coord_cartesian(xlim = x_limits, ylim = y_limits) + plot_theme +
  labs(title = "B   Manually gated ILCs", x = "UMAP 1", y = "UMAP 2") +
  theme(legend.position = "bottom", legend.title = element_blank())

marker_labels <- c(CD127 = "CD127", CD161 = "CD161", CD117 = "CD117", CD294 = "CRTH2")
marker_plots <- list()
for (marker in names(marker_labels)) {
  marker_plots[[marker]] <- ggplot(umap, aes(UMAP1, UMAP2, colour = .data[[marker]])) +
    geom_point(size = 0.45, alpha = 0.8) +
    scale_colour_gradientn(name = marker_labels[[marker]],
                           colours = rev(brewer.pal(10, "Spectral"))) +
    coord_cartesian(xlim = x_limits, ylim = y_limits) + plot_theme +
    labs(title = marker_labels[[marker]], x = NULL, y = NULL) +
    theme(plot.margin = margin(4, 4, 4, 4), legend.text = element_text(size = 8),
          legend.title = element_text(size = 9))
}
marker_grid <- arrangeGrob(marker_plots$CD127, marker_plots$CD117,
                           marker_plots$CD161, marker_plots$CD294, ncol = 2)
panel_B <- arrangeGrob(p_ilc, marker_grid, ncol = 2, widths = c(1, 1.45))
ggsave(file.path(output_dir, "Figure_1B.png"), panel_B, width = 12, height = 4.8, dpi = 300)


################################
########## Figure 1c ###########
################################
kde_cells <- umap[umap$ILC %in% c("ILC1", "ILC2", "ILC3"), ]
rainbow_density <- colorRampPalette(c("#000099", "#00FEFF", "#45FE4F",
                                      "#FCFF00", "#FF9400", "#FF3100"))
contour_colours <- colorRampPalette(c("#FFFFFF", brewer.pal(9, "YlOrRd")))(10)
plot_kde <- function(status, label) {
  cells <- kde_cells[kde_cells$COVID_status == status, ]
  xy <- cells[, c("UMAP1", "UMAP2")]
  smoothScatter(xy, transformation = function(z) z^0.4,
                colramp = rainbow_density, axes = FALSE, xlab = "", ylab = "")
  density <- MASS::kde2d(xy$UMAP1, xy$UMAP2, n = 50)
  contour(density, drawlabels = FALSE, nlevels = 9,
          col = contour_colours, add = TRUE)
  box(col = "grey60")
  mtext(label, side = 4, line = 0.5, font = 2, cex = 0.9)
}

png(file.path(output_dir, "Figure_1C.png"), width = 3.4, height = 4.8, units = "in", res = 300)
par(mfrow = c(2, 1), mar = c(1.6, 1.7, 1.6, 2.7), oma = c(0, 0, 0, 0))
plot_kde("healthy", "Healthy")
mtext("C                         KDE", side = 3, line = 0.35, font = 2, cex = 1)
plot_kde("COVID", "COVID-19 +")
mtext("UMAP 1", side = 1, line = 0.4, font = 2, cex = 0.8)
mtext("UMAP 2", side = 2, line = 0.4, font = 2, cex = 0.8)
dev.off()


################################
########## Figure 1D ###########
################################
mdf = covid_ilc@meta.data
sampleIDs = unique(mdf$Sample_id)
ggdf = data.frame()
ggdf2 = data.frame()
for (sid in sampleIDs){
  tmpdf = dplyr::filter(mdf, Sample_id %in% sid)
  tmp2 = dplyr::filter(tmpdf, handgated_cell_type2 %in% c("ILC1","ILC2","ILCp"))
  status = unique(tmpdf$COVID_status)
  freqdf = as.data.frame(table(tmp2$handgated_cell_type2))
  freqdf$Prop = freqdf$Freq/dim(tmpdf)[[1]] # % of total ilcs
  freqdf$Sample_id = sid
  freqdf$COVID_status = status
  ggdf = rbind(ggdf, freqdf)
  
  freqdf = as.data.frame(sum(table(tmp2$handgated_cell_type2)))
  colnames(freqdf) = "Tot_ILCs"
  freqdf$Prop = freqdf$Tot_ILCs/dim(tmpdf)[[1]] #tmpdf includes ilc and other cells
  freqdf$Sample_id = sid
  freqdf$COVID_status = status
  ggdf2 = rbind(ggdf2, freqdf)
  
}

clrs = c("#008000","#F0B548")

ggdf$COVID_status = factor(ggdf$COVID_status, levels = c("healthy","COVID"))

textdf = data.frame()
for (ilc in c("ILC1","ILC2","ILCp")){
  tmpdf = dplyr::filter(ggdf,  Var1 %in% ilc)
  x = tmpdf[tmpdf$COVID_status == "healthy",]$Prop
  y =  tmpdf[tmpdf$COVID_status == "COVID",]$Prop
  ymax = max(x,y) + 0.1 * (max(x,y))
  pval = wilcox.test(x,y)$p.value
  
  tdf = data.frame(x = 1.5,
                   y = ymax * 100,
                   Var1 = ilc,
                   COVID_status = "UNK",
                   label = paste("p = ",sprintf("%.4f", pval), sep=""))
  textdf = rbind(textdf,tdf)
}

### calculte p-value of tottal ILCs in coivid vs healthy ###
x = ggdf2$Prop[ggdf2$COVID_status=="COVID"]
y = ggdf2$Prop[ggdf2$COVID_status=="healthy"]
pval = wilcox.test(x,y)$p.value

textdf1 = data.frame(x = 1.5,
                     y = 90,
                     Var1 = "ILCs",
                     COVID_status = "UNK",
                     label = paste("p = ",sprintf("%.4f", pval), sep=""))

ggdf2$Var1 = "ILCs"
ggdf2$COVID_status = factor(ggdf2$COVID_status, levels = c("healthy","COVID"))


ggdf0 = rbind(ggdf2[,c("Prop","Sample_id","COVID_status","Var1")],
              ggdf[,c("Prop","Sample_id","COVID_status","Var1")])

textdf0 = rbind(textdf,textdf1)

ggdf0$Var1 = factor(ggdf0$Var1, levels = c("ILCs","ILC1","ILC2","ILCp"))
textdf0$Var1 = factor(textdf0$Var1, levels = c("ILCs","ILC1","ILC2","ILCp"))

mm = match(ggdf0$Sample_id, mdf$Sample_id)
ggdf0$PPID = mdf$PPID[mm]
ggdf0$longterm  = ifelse(ggdf0$PPID %in% long_ter_covd, "Long","Other")

g1 <- ggplot(ggdf0, aes(y=Prop*100, x=COVID_status)) +
  geom_boxplot(width = 0.5, alpha = 0.5,outlier.color = NA) +
  geom_point(aes(fill=COVID_status,shape = longterm), color = "black", alpha = 1.0, position = position_dodge(width = 0.01), size=3.2) + 
  theme_bw(base_size = 14) + 
  facet_wrap(~Var1, scales = "free", nrow = 1) +
  scale_fill_manual(values = clrs) +
  scale_shape_manual(values = c(24,21)) +
  theme(legend.position = "none",
        strip.background = element_rect(fill = "white"),
        strip.text = element_text(size=18),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
  ) +
  ylab(label = "% Total live cells") +
  geom_text(data=textdf0, mapping=aes(x=x, y=y, label=label), size=5)

png(file.path(output_dir, "Figure_1D.png"), height = 3.8, width = 6.2, units = "in", res = 300)
grid.arrange(g1, nrow=1)
dev.off()


######################################
########## Figure 1E and G ###########
######################################
ilc_types <- c("ILC1", "ILC2", "ILCp")
markers_EG <- c("CD69", "CD62L")

# Scale across ALL retained cells before selecting ILCs, as in the old analysis.
efg <- covid_ilc
DefaultAssay(efg) <- "Abseq"
efg <- NormalizeData(efg, normalization.method = "CLR", margin = 2, verbose = FALSE)
efg <- ScaleData(efg, features = rownames(efg[["Abseq"]]), verbose = FALSE)

healthy <- c("1033", "1348", "1804", "1987", "2195",
             "55689-0068", "55689-0071", "55689-0080", "WF001", "WF002",
             "WF006", "WF008", "WF010", "WF015", "WF016", "WF018", "WF019",
             "WF021", "WF025", "WF045", "WF026", "WF031", "WF033", "WF035", "WF042")
long_ter_covd <- c("55689-0007", "55689-0012", "55689-0015", "55689-0029",
                   "55689-0031", "55689-0034", "55689-0044", "55689-0050",
                   "55689-0057", "55689-0095", "55689-0165", "55689-0219",
                   "55689-0220", "OHC-001", "OHC-010", "OHC001", "OHC010")
efg$COVID_status <- ifelse(efg$PPID %in% healthy, "healthy", "COVID")
mdf_EFG <- efg[[]]
scaled_protein <- GetAssayData(efg, assay = "Abseq", layer = "scale.data")

Idents(efg) <- "handgated_cell_type2"
DE_list <- list()
for (cell_type in ilc_types) {
  de <- FindMarkers(efg, group.by = "COVID_status",
                    ident.1 = "COVID", ident.2 = "healthy", subset.ident = cell_type,
                    test.use = "wilcox", assay = "Abseq", slot = "data",
                    logfc.threshold = 0, min.pct = 0.1, verbose = FALSE)
  de$Protein <- rownames(de)
  de$ILCn <- cell_type
  DE_list[[cell_type]] <- de
}
DE_proteins <- do.call(rbind, DE_list)
statistics_EG <- DE_proteins[DE_proteins$Protein %in% markers_EG,
                             c("ILCn", "Protein", "p_val", "p_val_adj")]
# Preserve the old figure's "q" label. These are Seurat's Bonferroni-adjusted
# p-values across 29 proteins, not a separately calculated FDR/q-value procedure.
statistics_EG$label <- paste0("q = ", format(statistics_EG$p_val_adj,
                                             scientific = TRUE, digits = 2))
print(statistics_EG)


# PERCENTAGE-POSITIVE CELLS PER SAMPLE AND ILC SUBTYPE
percentage_list <- list()
for (cell_type in ilc_types) {
  sample_ids <- unique(mdf_EFG$Sample_id[mdf_EFG$handgated_cell_type2 %in% cell_type])
  for (sample_id in sample_ids) {
    cell_metadata <- mdf_EFG[mdf_EFG$Sample_id == sample_id &
                               mdf_EFG$handgated_cell_type2 %in% cell_type, ]
    cell_ids <- rownames(cell_metadata)
    if (length(cell_ids) >= 3) {
      percentage_list[[paste(cell_type, sample_id, sep = "_")]] <- data.frame(
        Sample_id = sample_id, PID = cell_metadata$PPID[1], ILCn = cell_type,
        COVID_status = cell_metadata$COVID_status[1], n_cells = length(cell_ids),
        CD69 = mean(scaled_protein["CD69", cell_ids] > 0) * 100,
        CD62L = mean(scaled_protein["CD62L", cell_ids] > 0) * 100)
    }
  }
}
percentages <- do.call(rbind, percentage_list)
rownames(percentages) <- NULL
percentages$ILCn <- factor(percentages$ILCn, levels = ilc_types)
percentages$COVID_status <- factor(percentages$COVID_status, levels = c("healthy", "COVID"))
percentages$longterm <- factor(ifelse(percentages$PID %in% long_ter_covd, "Long", "Other"),
                               levels = c("Long", "Other"))
status_colours <- c(healthy = "#008000", COVID = "#F0B548")
point_shapes <- c(Long = 24, Other = 21)


# FIGURES 1E and G BOXPLOTS 
boxplots_EG <- list()
for (marker in markers_EG) {
  annotation <- statistics_EG[statistics_EG$Protein == marker, ]
  annotation$ILCn <- factor(annotation$ILCn, levels = ilc_types)
  annotation$x <- 1.5
  annotation$y <- max(percentages[[marker]]) * 1.10
  panel_title <- c(CD69 = "E   CD69+ cells", CD62L = "G   CD62L (SELL)+ cells")[[marker]]
  boxplots_EG[[marker]] <- ggplot(percentages, aes(x = COVID_status, y = .data[[marker]])) +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = "white", colour = "grey35") +
    geom_point(aes(fill = COVID_status, shape = longterm), colour = "black",
               size = 3.2, position = position_dodge(width = 0.01)) +
    facet_wrap(~ILCn, scales = "free", nrow = 1) +
    scale_fill_manual(values = status_colours, labels = c("Healthy", "COVID-19+")) +
    scale_shape_manual(values = point_shapes, labels = c("Long-term COVID-19+", "Other")) +
    geom_text(data = annotation, aes(x = x, y = y, label = label),
              inherit.aes = FALSE, size = 4) +
    labs(title = panel_title, x = NULL, y = "% of the indicated ILC subtype") +
    theme_bw(base_size = 13) +
    theme(legend.position = "none", strip.background = element_rect(fill = "white"),
          strip.text = element_text(size = 14), axis.text.x = element_blank(),
          plot.title = element_text(face = "bold"))
}

ggsave(file.path(output_dir,  "Figure_1E.png"), boxplots_EG$CD69, width = 6.4, height = 4, dpi = 300)
ggsave(file.path(output_dir,  "Figure_1G.png"), boxplots_EG$CD62L, width = 6.4, height = 4, dpi = 300)


################################
########## Figure 1F ###########
################################
correlation_list <- list()
for (cell_type in ilc_types) {
  sample_data <- percentages[percentages$ILCn == cell_type, ]
  correlation <- cor.test(sample_data$CD69, sample_data$CD62L, method = "spearman",
                          exact = FALSE)  # Asymptotic calculation: these data have ties.
  correlation_list[[cell_type]] <- data.frame(ILCn = cell_type, n = nrow(sample_data),
                                              rho = unname(correlation$estimate), p_value = correlation$p.value)
}
statistics_F <- do.call(rbind, correlation_list)
statistics_F$ILCn <- factor(statistics_F$ILCn, levels = ilc_types)
statistics_F$label <- sprintf("rho = %.2f\np = %.3f", statistics_F$rho, statistics_F$p_value)
print(statistics_F)

percentages$display_group <- factor(
  ifelse(percentages$COVID_status == "healthy", "Healthy",
         ifelse(percentages$longterm == "Long", "Long-term COVID-19+", "COVID-19+")),
  levels = c("Healthy", "COVID-19+", "Long-term COVID-19+"))

plot_F <- ggplot(percentages, aes(x = CD69, y = CD62L)) +
  geom_point(aes(fill = display_group, shape = display_group), colour = "black", size = 3.4) +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, colour = "#6078BD", fill = "grey70") +
  facet_wrap(~ILCn, scales = "free", nrow = 1) +
  scale_fill_manual(values = c("Healthy" = "#008000", "COVID-19+" = "#F0B548",
                               "Long-term COVID-19+" = "#F0B548")) +
  scale_shape_manual(values = c("Healthy" = 21, "COVID-19+" = 21, "Long-term COVID-19+" = 24)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.25))) +
  geom_text(data = statistics_F, aes(x = Inf, y = Inf, label = label),
            inherit.aes = FALSE, hjust = 1.08, vjust = 1.1, size = 3.6) +
  labs(title = "F", x = "CD69+ (% of ILC subtype)",
       y = "CD62L (SELL)+ (% of ILC subtype)", fill = NULL, shape = NULL) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 13), plot.title = element_text(face = "bold"))

ggsave(file.path(output_dir, "Figure_1F.png"), plot_F, width = 9, height = 4, dpi = 300)




# Reconstruct the paper's RNA + AbSeq count object using CSV inputs only.
# Includes all cells retained by the original sample-tag and RNA quality filters.
# Requires Seurat and Matrix. Edit the three paths below before running.
library(Seurat)
library(Matrix)

setwd("./")

########  File paths ##########

raw_dir <- "./DataFiles/processed_data/count_and_tag_files/"
gating_file <- "./DataFiles/processed_data/paper_handgated_cell_labels.csv"
output_file <- "covid_ilc_counts.rds"
st = "./DataFiles/processed_data/count_and_tag_files/"
fp = "./DataFiles/processed_data/handgatedILCs/"
fileinfo = "./DataFiles/processed_data/FileInfo.csv"

healthy <- c("1033", "1348", "1804", "1987", "2195", "55689-0068", "55689-0071",
             "55689-0080", "WF001", "WF002", "WF006", "WF008", "WF010", "WF015",
             "WF016", "WF018", "WF019", "WF021", "WF025", "WF045", "WF026",
             "WF031", "WF033", "WF035", "WF042")
long_ter_covd <- c("55689-0007", "55689-0012", "55689-0015", "55689-0029",
                   "55689-0031", "55689-0034", "55689-0044", "55689-0050",
                   "55689-0057", "55689-0095", "55689-0165", "55689-0219",
                   "55689-0220", "OHC-001", "OHC-010", "OHC001", "OHC010")

##### LOAD FILES ######

md <- read.csv(fileinfo)
gating_files <- list.files(fp,pattern = "ILC",full.names = TRUE)

gating <- do.call(rbind, lapply(gating_files, function(file) {
  filename <- basename(file)
  cell_type <- sub("_.*", "", sub("^DBEC_", "", filename))
  batch <- sub(".*_", "", sub("\\.csv$", "", filename))
  cells <- names(read.csv(file, nrows = 0, check.names = FALSE))[-1]
  cells <- sub("^X", "", cells)
  
  calls <- read.csv(paste0(st,"/Batch_", batch, "_Sample_Tag_Calls.csv"))
  
  batch_md <- md[md$batch == batch, ]
  tags <- calls$Sample_Tag[match(cells, calls$Cell_Index)]
  sample_ids <- batch_md$Sample_id[match(tags, batch_md$Sample_Tag)]
  
  data.frame(
    cell_id = paste0(batch, "_", cells[!is.na(sample_ids)]),
    handgated_cell_type2 = cell_type
  )
}))

batches <- c("131a", "131b", "131c", "131d")
count_files <- file.path(
  raw_dir, paste0("Batch_", batches, "_DBEC_MolsPerCell.csv")
)
tag_files <- file.path(
  raw_dir, paste0("Batch_", batches, "_Sample_Tag_Calls.csv")
)
counts <- lapply(count_files, read.csv, skip = 7, check.names = FALSE)
tags <- lapply(tag_files, read.csv, check.names = FALSE)

#####  PROCESS FILES ##### 

feature_names <- names(counts[[1]])[-1]  # All columns except Cell_Index.

for (i in seq_along(batches)) {
  counts[[i]]$batch <- batches[i]
  counts[[i]]$cell_id <- paste0(batches[i], "_", counts[[i]]$Cell_Index)
  tag_match <- match(counts[[i]]$Cell_Index, tags[[i]]$Cell_Index)
  counts[[i]]$Sample_Tag <- tags[[i]]$Sample_Tag[tag_match]
  counts[[i]]$Sample_Name <- tags[[i]]$Sample_Name[tag_match]
}

counts <- do.call(rbind, counts)
counts <- counts[!counts$Sample_Tag %in% c("Multiplet", "Undetermined"), ]

protein_features <- feature_names[grepl("pAbO", feature_names)]
rna_features <- setdiff(feature_names, protein_features)
rna_counts <- t(as.matrix(counts[, rna_features]))
protein_counts <- t(as.matrix(counts[, protein_features]))
colnames(rna_counts) <- colnames(protein_counts) <- counts$cell_id
rownames(rna_counts) <- sub("\\|.*", "", rna_features)
rownames(protein_counts) <- sub("[|:].*", "", protein_features)
rownames(rna_counts) <- gsub("_", "-", rownames(rna_counts))
rownames(protein_counts) <- gsub("_", "-", rownames(protein_counts))
keep <- colSums(rna_counts > 0) > 15 & colSums(rna_counts) < 10000
rna_counts <- Matrix(rna_counts[, keep], sparse = TRUE)
protein_counts <- Matrix(protein_counts[, keep], sparse = TRUE)
cell_info <- counts[keep, c("cell_id", "batch", "Cell_Index", "Sample_Tag", "Sample_Name")]
rownames(cell_info) <- cell_info$cell_id


##### CREATE SEURAT OBJECT ##### 

covid_ilc <- CreateSeuratObject(counts = rna_counts, project = "COVID_ILC")
covid_ilc[["Abseq"]] <- CreateAssayObject(counts = protein_counts)
DefaultAssay(covid_ilc) <- "RNA"

##### LOAD METADATA ##### 
sample_metadata <- read.csv(file.path(raw_dir, "BD_Sample_tag_metadata.csv"),
                            colClasses = "character", na.strings = c("", "NA", "N/A"))
sample_key <- paste(sample_metadata$batch, sample_metadata$Sample_Tag, sep = "_")
cell_key <- paste(cell_info$batch, cell_info$Sample_Tag, sep = "_")
sample_match <- match(cell_key, sample_key)
metadata_fields <- setdiff(names(sample_metadata), c("batch", "Sample_Tag"))
metadata <- cbind(cell_info, sample_metadata[sample_match, metadata_fields])
rownames(metadata) <- metadata$cell_id
metadata$Sample_id <- paste0(metadata$batch, "_", metadata$Sample_id)

healthy <- c("1033", "1348", "1804", "1987", "2195",
             "55689-0068", "55689-0071", "55689-0080", "WF001", "WF002",
             "WF006", "WF008", "WF010", "WF015", "WF016", "WF018", "WF019",
             "WF021", "WF025", "WF045", "WF026", "WF031", "WF033", "WF035", "WF042")
long_ter_covd = c("55689-0007","55689-0012","55689-0015","55689-0029",
                  "55689-0031","55689-0034","55689-0044","55689-0050",
                  "55689-0057","55689-0095","55689-0165","55689-0219",
                  "55689-0220","OHC-001","OHC-010", 'OHC001', 'OHC010')
other_covid = c("55689-0038","55689-0046","55689-0047","55689-0053",
                "55689-0091","55689-0157","55689-0218")

metadata$COVID_status_original <- metadata$COVID_status
metadata$COVID_status <- ifelse(metadata$PPID %in% healthy, "healthy", "COVID")
gating_match <- match(metadata$cell_id, gating$cell_id)

# labels <- data.frame(
#   cell_id = covid_ilc$cell_id,
#   handgated_cell_type2 =
#     gating$handgated_cell_type2[match(covid_ilc$cell_id, gating$cell_id)]
# )
# labels <- labels[!is.na(labels$handgated_cell_type2), ]

#write.csv(labels,gating_file,row.names = FALSE)
#gating <- read.csv(gating_file, colClasses = "character")

metadata$handgated_cell_type2 <- gating$handgated_cell_type2[gating_match]
covid_ilc <- AddMetaData(covid_ilc, metadata = metadata)

saveRDS(covid_ilc, file = output_file)

###########################
#### Counting ILCs ########
###########################
# Loading seurat object #
totalILCs  = sum(table(covid_ilc@meta.data$handgated_cell_type2))
#> totalILCs
#[1] 6247
#lets count ILC1 
sum(covid_ilc@meta.data$handgated_cell_type2=="ILC1", na.rm = T)
#[1] 2103
#lets count ILC2
sum(covid_ilc@meta.data$handgated_cell_type2=="ILC2", na.rm = T)
#[1] 2509
#lets count ILCp 
sum(covid_ilc@meta.data$handgated_cell_type2=="ILCp", na.rm = T)
#[1] 1635

##############################
## Now lets see CD161+ ILCs ##
DefaultAssay(covid_ilc) <- "Abseq"
ExpressMat1 = t(GetAssayData(object = covid_ilc, layer = "counts")) ## rows are cell IDs and columns are genes from scRNAseq
ExpressMat1 = as.data.frame(ExpressMat1)  ## ExpressMat1 is scaled rnaseq 
mdf = covid_ilc@meta.data ## meta-data
mdf$COVID_status = ifelse(mdf$PPID %in% healthy,  "healthy","COVID")
idx = mdf$handgated_cell_type2 %in% c(c("ILC1","ILC2","ILCp"))
ILCmatrix = ExpressMat1[idx,]
summary(ILCmatrix$CD161)
#> summary(ILCmatrix$CD161)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#10.00   32.00   69.00   98.44  135.00  882.00 
## All HANDGATED ILCs has CD161 (see supp. fig. 1 for schema) ##

######################################################
#### Lets review ICAM1 in ILCs across 4 batches ######
######################################################
DefaultAssay(covid_ilc) <- 'RNA'
ExpressMat1 = GetAssayData(object = covid_ilc, layer = "counts") ## other option is "scale.data" ## rows are cell IDs and columns are genes from scRNAseq
ExpressMat1 = as.data.frame(ExpressMat1)  ## ExpressMat1 is scaled rnaseq 
mdf = covid_ilc@meta.data ## meta-data
mdf$COVID_status = ifelse(mdf$PPID %in% healthy,  "healthy","COVID")
dmat = as.data.frame(t(ExpressMat1))
mm = match(rownames(dmat), mdf$cell_id)
dmat$covidstatus = mdf$COVID_status[mm]
dmat$ILCn = mdf$handgated_cell_type2[mm] ## cell ID to ILCn 
dmat$batch = mdf$batch[mm] ## cell ID to batchn 
idx = dmat$ILCn %in% c(c("ILC1","ILC2","ILCp"))
ILCmatrixRNA = dmat[idx,]

(batch_averages <- ILCmatrixRNA %>% group_by(batch) %>%
    summarize(average_ICAM1 = mean(ICAM1, na.rm = TRUE)))
# A tibble: 4 × 2
#batch average_GAPDH
#<chr>         <dbl>
#  1 131a          0.193
#. 2 131b          0.311
#  3 131c          0.218
#  4 131d          0.207

### Any observed batch effect wahtsoever was resolved or adjusted during modeling and p-value calculation ###

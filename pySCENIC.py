#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 22 18:12:26 2022
https://www.nature.com/articles/s41596-020-0336-2#Sec13
https://link.springer.com/protocol/10.1007/978-1-0716-1534-8_10
http://htmlpreview.github.io/?https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/SCENIC%20Protocol%20-%20Case%20study%20-%20Cancer%20data%20sets.html
https://github.com/aertslab/SCENICprotocol/blob/master/notebooks/PBMC10k_SCENIC-protocol-CLI.ipynb
@author: abhinavkaushik
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
import json
import zlib
import base64
import umap
from pyscenic.binarization import binarize
from pyscenic.export import export2loom, add_scenic_metadata
from pyscenic.utils import load_motifs
from pyscenic.transform import df2regulons
from pyscenic.aucell import aucell
from pyscenic.binarization import binarize
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_binarization, plot_rss

# set a working directory
wdir = "/Users/abk347/Library/CloudStorage/OneDrive-HarvardUniversity/ILC/Analysis/Result_replication/"
os.chdir( wdir )


f_loom_path_unfilt = "PyScenic_output/unfilt.loom" # path to anndata object, which will be updated to store Scanpy results as they are generated below
f_loom_path_scenic = "PyScenic_output/ILCs_filtered_scenic.loom"
f_tfs = "PyScenic_output/allTFs_hg38.txt" # human

f_anndata_path = "PyScenic_output/anndata.h5ad"
f_exprMat = 'PyScenic_output/exprMat.tsv' ## unfiltered count matrix
fname_meta="PyScenic_output/MataData.csv"
f_PyScenic_output_output = "PyScenic_output/pyscenic_output.loom"

adata = sc.read_text( f_exprMat, delimiter=' ', first_column_names=True )

row_attrs = {
    "Gene": np.array(adata.var.index) ,
}
col_attrs = {
    "CellID":  np.array(adata.obs.index) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}

lp.create( f_loom_path_unfilt, adata.X.transpose(), row_attrs, col_attrs )

#############################
## Initial/basic filtering ##
#############################

# read unfiltered data from a loom file
adata = sc.read_loom( f_loom_path_unfilt )
nCountsPerGene = np.sum(adata.X, axis=0)
nCellsPerGene = np.sum(adata.X>0, axis=0)

# Show info
print("Number of counts (in the dataset units) per gene:", nCountsPerGene.min(), " - " ,nCountsPerGene.max())
print("Number of cells in which each gene is detected:", nCellsPerGene.min(), " - " ,nCellsPerGene.max())
nCells=adata.X.shape[0]

# PyScenic_output thresholds
minCountsPerGene=3*.01*nCells # 3 counts in 1% of cells
print("minCountsPerGene: ", minCountsPerGene)

minSamples=.01*nCells # 1% of cells
print("minSamples: ", minSamples)

# simply compute the number of genes per cell (computers 'n_genes' column)
sc.pp.filter_cells(adata, min_genes=0)
# add the total counts per cell as observations-annotation to adata
adata.obs['n_counts'] = adata.X.sum(axis=1).A1


sc.pl.violin(adata, ['n_genes', 'n_counts'],
    jitter=0.4, multi_panel=True )

#####################################
### Carry out the filtering steps ###
#####################################
# initial cuts
sc.pp.filter_cells(adata, min_genes=10 )
sc.pp.filter_genes(adata, min_cells=3 )


## Finalize the selected filters
adata.write( f_anndata_path )
# create basic row and column attributes for the loom file:
row_attrs = {
    "Gene": np.array(adata.var_names) ,
}
col_attrs = {
    "CellID": np.array(adata.obs_names) ,
    "nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
    "nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
}
lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs)

################################################
### Further pre-processing of expression data ##
################################################
# save a copy of the raw data
adata.raw = adata
# Total-count normalize (library-size correct) to 10,000 reads/cell
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
# log transform the data.
sc.pp.log1p(adata)
# identify highly variable genes.
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pl.highly_variable_genes(adata)
# keep only highly variable genes:
adata = adata[:, adata.var['highly_variable']]
# scale each gene to unit variance, clip values exceeding SD 10.
sc.pp.scale(adata, max_value=10)
# update the anndata file:
adata.write( f_anndata_path )

##################
### PCA & UMAP ###
##################
sc.tl.pca(adata, svd_solver='arpack')
sc.pl.pca_variance_ratio(adata, log=True)
adata.write( f_anndata_path )
# neighborhood graph of cells (determine optimal number of PCs here)
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=5)
# compute UMAP
sc.tl.umap(adata)
# tSNE
tsne = TSNE( n_jobs=20 )
adata.obsm['X_tsne'] = tsne.fit_transform( adata.X )
adata.write( f_anndata_path )

################
# cluster graph #
################
# cluster the neighbourhood graph
sc.tl.louvain(adata,resolution=0.4)
sc.pl.umap(adata, color=['louvain'] )

######################
### SCENIC steps ####
## STEP 1: Gene regulatory network inference, and generation of co-expression modules
## Phase Ia: GRN inference using the GRNBoost2 algorithm
######################
## [CLI] pyscenic grn {f_loom_path_scenic} {f_tfs} -o adj.csv --num_workers 20


! pyscenic grn PyScenic_output/ILCs_filtered_scenic.loom PyScenic_output/allTFs_hg38.txt -o PyScenic_output/adj.csv --num_workers 6

adjacencies = pd.read_csv("PyScenic_output/adj.csv", index_col=False)
adjacencies.head()

#######################################################
### STEP 2-3: Regulon prediction aka cisTarget from CLI
import glob
# ranking databases
f_db_glob = "PyScenic_output/gene_based/*rankings.feather"
f_db_names = ' '.join( glob.glob(f_db_glob) )
f_motif_path = "PyScenic_output/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl"

## agin use CLI
! pyscenic ctx PyScenic_output/adj.csv {f_db_names}  \
--annotations_fname PyScenic_output/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname PyScenic_output/ILCs_filtered_scenic.loom \
--output PyScenic_output/reg.csv \
--mask_dropouts \
--num_workers 6

#######################################################
## STEP 4: Cellular enrichment (aka AUCell) from CLI
nGenesDetectedPerCell = np.sum(adata.X>0, axis=1)
nGenesDetectedPerCell = pd.DataFrame(nGenesDetectedPerCell)
percentiles = nGenesDetectedPerCell.quantile([.01, .05, .10, .50, 1])
print(percentiles)


fig, ax = plt.subplots(1, 1, figsize=(8, 5), dpi=150)
sns.distplot(nGenesDetectedPerCell, norm_hist=False, kde=False, bins='fd')
for i,x in enumerate(percentiles):
    fig.gca().axvline(x=x, ymin=0,ymax=1, color='red')
    ax.text(x=x, y=ax.get_ylim()[1], s=f'{int(x)} ({percentiles.index.values[i]*100}%)', color='red', rotation=30, size='x-small',rotation_mode='anchor' )
ax.set_xlabel('# of genes')
ax.set_ylabel('# of cells')
fig.tight_layout()

# path to PyScenic_output output
! pyscenic aucell PyScenic_output/ILCs_filtered_scenic.loom PyScenic_output/reg.csv \
    --output {f_PyScenic_output_output} --num_workers 20


###### Visualization of SCENIC's AUC matrix ########
# collect SCENIC AUCell output
lf = lp.connect( f_PyScenic_output_output, mode='r+', validate=False )
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
lf.close()

# UMAP
runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4, metric='correlation').fit_transform
dr_umap = runUmap( auc_mtx )
pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "PyScenic_output/scenic_umap.txt", sep='\t')
# tSNE
tsne = TSNE( n_jobs=20 )
dr_tsne = tsne.fit_transform( auc_mtx )
pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "PyScenic_output/scenic_tsne.txt", sep='\t')

# scenic output
lf = lp.connect( f_PyScenic_output_output, mode='r+', validate=False )
meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
#exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID)
auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
regulons = lf.ra.Regulons
dr_umap = pd.read_csv( 'PyScenic_output/scenic_umap.txt', sep='\t', header=0, index_col=0 )
dr_tsne = pd.read_csv( 'PyScenic_output/scenic_tsne.txt', sep='\t', header=0, index_col=0 )
###
#Fix regulon objects to display properly in SCope:

auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(')
regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )
# regulon thresholds
rt = meta['regulonThresholds']
for i,x in enumerate(rt):
    tmp = x.get('regulon').replace("(","_(")
    x.update( {'regulon': tmp} )

##Concatenate embeddings (tSNE, UMAP, etc.)

########### plotting and analysis ##############
metadata = pd.read_csv(fname_meta)
df_obs = metadata[['cell_id', 'Sample_id', 'handgated_cell_type2', 'batch', 'PPID', 'COVID_status','Long_Term']].set_index('cell_id').sort_index()
adata.obs = df_obs
adata.var_names_make_unique()
adata.raw = adata
adata.write_h5ad("PyScenic_output/adata.hd5")


auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
auc_mtx.to_csv("PyScenic_output/AUC_MTX.csv")

##########################################
#### Binarize the AUcell value matrix ####
bin_mtx, thresholds = binarize(auc_mtx)
bin_mtx.to_csv("PyScenic_output/BIN_MTX.csv")
thresholds.to_frame().rename(columns={0:'threshold'}).to_csv("PyScenic_output/Threshold.csv")
thresholds = thresholds.to_frame().rename(columns={0:'threshold'})

###############################################################################
#### Create cluster map; columns are cells/cell types and rows are regulon; ###
#### black means value of '1' for a given regulon in a given cell ##

N_COLORS = len(adata.obs.handgated_cell_type2.dtype.categories)
COLORS = [color['color'] for color in mpl.rcParams["axes.prop_cycle"]]
cell_type_color_lut = dict(zip(adata.obs.handgated_cell_type2.dtype.categories, COLORS))
cell_id2cell_type_lut = metadata.set_index('cell_id').handgated_cell_type2.to_dict()
bw_palette = sns.xkcd_palette(["white", "black"])

sns.set()
sns.set(font_scale=1.0)
sns.set_style("ticks", {"xtick.minor.size": 1, "ytick.minor.size": 0.1})
g = sns.clustermap(bin_mtx.T,
               col_colors=auc_mtx.index.map(cell_id2cell_type_lut).map(cell_type_color_lut),
               cmap=bw_palette, figsize=(20,20))
g.ax_heatmap.set_xticklabels([])
g.ax_heatmap.set_xticks([])
g.ax_heatmap.set_xlabel('Cells')
g.ax_heatmap.set_ylabel('Regulons')
g.ax_col_colors.set_yticks([0.5])
g.ax_col_colors.set_yticklabels(['Cell Type'])
g.cax.set_visible(False)
g.fig.savefig("PyScenic_output/ClusterMap.png", format='png', dpi=400)

# #############################
# ## tsne based on AUC values #
# sc.tl.tsne(adata, use_rep='X_aucell')
# embedding_aucell_tsne = pd.DataFrame(adata.obsm['X_tsne'], columns=[['_X', '_Y']], index=adata.obs_names)
# sc.set_figure_params(frameon=False, dpi=150, fontsize=8)
# sc.pl.tsne(adata, color=['handgated_cell_type2'], ncols=1, palette=COLORS,save='PyScenic_output/AUCell_tSNE.png')

# #############################
# ## UMAP based on AUC values #
# sc.tl.umap(adata, use_rep='X_aucell')
# embedding_aucell_tsne = pd.DataFrame(adata.obsm['X_tsne'], columns=[['_X', '_Y']], index=adata.obs_names)
# sc.set_figure_params(frameon=False, dpi=150, fontsize=8)
# sc.pl.tsne(adata, color=['handgated_cell_type2'], ncols=1, palette=COLORS,save='PyScenic_output/AUCell_tSNE.png')
# sc.pl.umap(adata, color=['handgated_cell_type2'], ncols=1, palette=COLORS,save='PyScenic_output/AUCell_umap.png')


#### CELL TYPE SPECIFIC REGULATORS - Z-SCORE ####
## compute the mean value of AUcell values for each cell type ##
## so that for each regulon we have a Z score (scaled score) ##
## those cell type where Z-score < 3, then make Z score NULL ##

## use R to do the below steps #
# df_scores = auc_mtx[signature_column_names]

# ## formula (mean - mean)/standard_deviation
# df_results = ((df_scores.groupby(by='handgated_cell_type2').mean() - df_obs[signature_column_names].mean())/
#                df_obs[signature_column_names].std()).stack().reset_index().rename(columns={'level_1': 'regulon', 0:'Z'})

# df_results['regulon'] = list(map(lambda s: s[8:-1], df_results.regulon))
# df_results[(df_results.Z >= 3.0)].sort_values('Z', ascending=False).head()


#### CELL TYPE SPECIFIC REGULATORS - RSS ###
rss = regulon_specificity_scores(auc_mtx, adata.obs.COVID_status)
rss.head()
rss.to_csv("PyScenic_output/RSS_1_MTX.csv")

rss = regulon_specificity_scores(auc_mtx, adata.obs.handgated_cell_type2)
rss.head()
rss.to_csv("PyScenic_output/RSS_2_MTX.csv")

adata.obs['CT']  = adata.obs['handgated_cell_type2'].astype(str) + "_" + adata.obs['COVID_status'].astype(str)
rss = regulon_specificity_scores(auc_mtx, adata.obs.CT)
rss.head()
rss.to_csv("PyScenic_output/RSS_3_MTX.csv")


sns.set()
sns.set(style='whitegrid', font_scale=0.8)
fig, ((ax1, ax2, ax3)) = plt.subplots(1, 3, figsize=(8, 6), dpi=100)
plot_rss(rss, 'ILC1', ax=ax1)
ax1.set_xlabel('')
plot_rss(rss, 'ILC2', ax=ax2)
ax2.set_xlabel('')
ax2.set_ylabel('')
plot_rss(rss, 'ILCp', ax=ax3)
ax3.set_xlabel('')
ax3.set_ylabel('')
plt.tight_layout()
savesvg('PyScenic_output/rss.svg', fig)






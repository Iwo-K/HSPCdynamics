# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.12.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% tags=[]
#Singularity container used:
# !echo $SINGULARITY_CONTAINER
# !hostname

# %% [markdown]
# # Integrating SS2 and 10x data into a single landscape

# %% [markdown]
# ## Setup

# %% tags=[]
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import pathlib
import re
import anndata
import pickle
from scipy.sparse import csr_matrix
import utils.proc, utils.plots
from utils.load10x import load_10xfiles
import cellproject as cp
import utils.annotator as an
from anndata import AnnData
from utils.HSCscore import calculate_HSCscore
from utils.fixes import unclog_umap_caching

import warnings
warnings.filterwarnings("ignore", category=FutureWarning) 

sc.settings.verbosity = 3
sc.settings.figdir = './figures/04script/'
base_figures = './figures/04script/'
base_procdata = './procdata/04script/'
for i in [sc.settings.figdir, base_figures, base_procdata]:
    pathlib.Path(i).mkdir(parents=True, exist_ok=True)

# %% tags=[]
#Figure settings
sc.set_figure_params(color_map='viridis', dpi_save=350)
mpl.rcParams['figure.figsize'] = (4, 4) #1:1 plot ratio
mpl.rc('xtick', labelsize=14)
mpl.rc('ytick', labelsize=14)
mpl.rc('axes', labelsize=16)
mpl.rc('axes', labelsize=16)
mpl.rcParams['pdf.fonttype'] = 42 #Ensures readable fonts in illustrator
mpl.rcParams['ps.fonttype'] = 42
plt.rc('axes', axisbelow=True) #Ensure that gridlines are behind points
mplparams = mpl.rcParams.copy()

# %% tags=[]
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# Loading genesets

# %% tags=[]
# Getting cell cycle genes per phase - for cell cycle assignment
cycle_genes = pd.read_csv('./data/genesets/macosko2015_hum_ccgenes.csv') #These are human cell cycle genes with just names swapped for mouse
Sgenes = cycle_genes.loc[cycle_genes['phase'] == 'S','gene']
G2Mgenes = cycle_genes.loc[cycle_genes['phase'] == 'G2_M','gene']

# %% tags=[]
# loading 10x data
meta10x = pd.read_csv('./data/10x/scB5Tom_10x_samplemeta.csv')
meta10x.index = meta10x.biosample_id
adata10x = sc.read('./procdata/02script/combined10x_qced_logn.h5ad')

# loading SS2 data
adataSS2 = sc.read('./procdata/01script/SS2data.h5ad')
adataSS2 = adataSS2[:, adata10x.var.gene_ids]
adataSS2.var.index = adata10x.var.index

# We will be using logn data for everything
utils.proc.lognorm(adataSS2)

# %% tags=[]
# Selecting genes which should be excluded from the highly-variable gene list
ygenes = pd.read_csv('./data/genesets/Y_genes.csv')
realY_filter = [not bool(re.search("predicted gene", i)) for i in ygenes['Gene description']]
ygenes = ygenes.loc[realY_filter,'Gene name']
toremove =  ygenes.append(pd.Series(['Xist']))

# Using method from Weinreb et al. 2020, finding genes which are correlated with
# cell cycle signature, at least 0.1 Pearson R
ccgenes = utils.proc.find_correlated_genes(adata10x,
                      data_scaled = False,
                      genes = ['Ube2c', 'Hmgb2', 'Hmgn2', 'Tuba1b', 'Ccnb1', 'Tubb5', 'Top2a', 'Tubb4b'])
toremove = toremove.append(pd.Series(ccgenes))

# %% [markdown]
# ## Looking at processed SS2 data alone

# %% tags=[]
#Assigning cell cycle phase to SS2 data
adataSS2_scaled = sc.pp.scale(adataSS2, copy = True)
Sgenes_ = Sgenes.loc[Sgenes.isin(adataSS2.var.symbol)].copy()
G2Mgenes_ = G2Mgenes.loc[G2Mgenes.isin(adataSS2.var.symbol)].copy()
utils.proc.score_genes_cell_cycle_fixed(adataSS2_scaled, Sgenes_, G2Mgenes_, use_raw = False)
adataSS2.obs = adataSS2_scaled.obs.copy()

# Quick processing of SS2 data to get overview and see if there are any strong batch effects
adataSS2_proc = utils.proc.process(adataSS2,
                  compute_hivar = True,
                  n_pcs = 50,
                  n_neighbors = 12,
                  lognorm_data = False,
                  n_variable = 5000,
                  remove_genes = toremove,
                  Sgenes = Sgenes,
                  G2Mgenes = G2Mgenes, 
                  regress_cc = False,
                  slim = True)
sc.pl.umap(adataSS2_proc, color = 'biosample_id')
sc.pl.umap(adataSS2_proc, color = ['phase', 'Procr', 'Klf1', 'Dntt', 'Elane', 'Irf8', 'Cma1', 'Pf4'])

utils.plots.umap_subgroups(adataSS2_proc, key = 'biosample_id', toplot = adataSS2.obs.biosample_id.unique())

# %% [markdown] tags=[]
# ## SS2 and 10x data integration

# %% [markdown]
# Basic processing of 10x data

# %%
#Remove the 3d outlier sample
adata10x = adata10x[adata10x.obs.biosample_id !='3d7d_SITTB7_146814.46e',:].copy()

#Removing the all 0 genes
adata10x = adata10x[:, adata10x.var.index[(adata10x.X.sum(axis=0).A1 > 0)]].copy()

# %% tags=[]
# Running basic processing to get highly variable genes, assign cell cycle etc, 
# storing only the highly variable genes .X while keeping all log-normalised values in .raw.X

adata10x = utils.proc.process(adata10x,
                  compute_hivar = True,
                  n_pcs = 50,
                  n_neighbors = 12,
                  lognorm_data = False,
                  n_variable = 5000,
                  remove_genes = toremove,
                  Sgenes = Sgenes,
                  G2Mgenes = G2Mgenes, 
                  regress_cc = False,
                  slim = True)
#For Seurat correction, we will need the logn values and not scaled, getting from .raw.X
adata10x.X = adata10x.raw[:, adata10x.var.index].X.copy()

#Checking if the landscape looks as expected
sc.pl.umap(adata10x, color = ['phase', 'leiden'], wspace = 1.2)
sc.pl.umap(adata10x, color = ['Procr', 'Klf1', 'Dntt', 'Elane', 'Irf8', 'Cma1', 'Pf4'])

# %% [markdown]
# Joining SS2 and 10x data

# %% tags=[]
adataSS2.raw = adataSS2[:,adata10x.raw.var.index]
adataSS2 = adataSS2[:, adata10x.var.index].copy()
comb = adata10x.concatenate(adataSS2, batch_key = 'data_type', batch_categories = ['10x', 'SS2'])
comb.var['highly_variable'] = comb.var['highly_variable-10x'].copy() #Using variable genes from 10x

# %% [markdown]
# Combining metadata

# %% tags=[]
#Unifyign the mouse sex information (in 10x inferred, in SS2 from metadata)
sex_comb = [comb.obs.loc[x, 'sex_adata'] if comb.obs.loc[x, 'data_type'] == '10x' else comb.obs.loc[x, 'sex']  for x in comb.obs.index]
comb.obs['sex_combined'] = sex_comb

#Creating handy metadata with each rows as a biological sample
meta = comb.obs[['biosample_id', 'data_type', 'start_age', 'sex_combined', 'timepoint_tx_days', 'tom', 'batch', 'mouse_id']]
meta = meta.drop_duplicates(ignore_index = True)
meta.index = meta.biosample_id.values
meta['ncells'] = comb.obs.biosample_id.value_counts()[meta.index]
meta['longname'] = (meta.biosample_id.astype('str') + 
    meta.start_age.astype('str') + "\n" + 
    meta.ncells.astype('str') + "cells_" +
    meta.sex_combined.astype('str') + "_" + 
    meta.timepoint_tx_days.astype('str') + "d" +
    "_Tom" + meta.tom.astype(str))

meta.to_csv(base_procdata + 'combined_sample.meta.csv')
meta

# %% [markdown]
# Running Seurat integration

# %% [markdown]
# ### SS2 - 10x batch correction

# %% tags=[]
comb_cor = cp.run_SeuratCCA(comb, debug = False)
# As the correction step takes a couple hours, saving/reading as a checkpoint
comb_cor.write(base_procdata + 'comb_10x_SS2_seurat_cor.h5ad', compression ='lzf')
comb_cor = sc.read(base_procdata + 'comb_10x_SS2_seurat_cor.h5ad')

# %% [markdown]
# ## PCA landscape

# %% [markdown]
# ### Without Harmony correction

# %% tags=[]
# Setting seurat-corrected expression as .X and saving it for later use
comb.X = comb_cor[comb.obs.index, comb.var.index].X.copy()
comb.layers['seurat'] = comb.X.copy()

sc.pp.scale(comb)
sc.tl.pca(comb, n_comps = 50)

# %% [markdown]
# ### With Harmony correction

# %% tags=[]
sc.external.pp.harmony_integrate(comb, key = 'sample_id', max_iter_harmony = 20)
sc.pp.neighbors(comb, n_pcs = 50, n_neighbors = 12, use_rep = 'X_pca_harmony')
sc.tl.leiden(comb, resolution=0.85)

#UMAP
umapref = cp.quick_umap(comb, use_rep='X_pca_harmony')
comb.obsm['X_umap_2d'] = comb.obsm['X_umap'].copy()
umapref_3d = cp.quick_umap(comb, n_components=3, use_rep='X_pca_harmony')
comb.obsm['X_umap_3d'] = comb.obsm['X_umap'].copy()
comb.obsm['X_umap'] = comb.obsm['X_umap_2d'].copy()

sc.pl.umap(comb, color = ['biosample_id', 'phase', 'data_type', 'leiden', 'doublet_scores', 'predicted_doublets'], wspace = 2.4, save='PCA_info.pdf')
sc.pl.umap(comb, color = ['Procr', 'Klf1', 'Dntt', 'Elane', 'Irf8', 'Cma1', 'Pf4', 'Ms4a2'], save='PCA_markers1.pdf')

comb.obs['longname'] = meta.loc[comb.obs.biosample_id, 'longname'].values
utils.plots.umap_subgroups(comb, key = 'longname', toplot = comb.obs.longname.unique(), file=base_figures+'umapPCA_by_biosample.pdf')

# %% [markdown]
# ### Annotation

# %% tags=[]
sc.pl.umap(comb, color = 'leiden', legend_loc = 'on data', save = 'PCA_clusters.pdf')
utils.plots.umap3d(comb, color='leiden', key='X_umap_3d', filename=base_figures+'umap3dPCA_leiden.html')
#DE to help with annotation
sc.tl.rank_genes_groups(comb, groupby='leiden')
sc.pl.rank_genes_groups(comb, sharey=False, save = '_PCA_rankgenesgroups.pdf', n_genes=40)

#Plotting markers to help with annotation
sc.pl.umap(comb, color = ['Ccr2', 'Spn', 'Cx3cr1', 'Itgam', 'Car1', 'Klf1', 'Hba-a1', 'Bcl11a'], save='PCA_addmarkers2.pdf')
sc.pl.umap(comb, color = ['Dntt', 'Ly6d', 'Cd79a', 'Cd79b', 'Cd3e', 'Bcl11b', 'Gata3', 'Id2',
                         'Cd8a', 'Gzmb'], save = 'PCA_Lymarkers2.pdf')
sc.pl.umap(comb, color = ['Ms4a2', 'Cma1', 'Gzmb', 'Prss34', 'Mcpt8', 'Prg2', 'Prg3'], save = 'PCA_baseomc_markers3.pdf')
sc.pl.umap(comb, color = ['Mpo', 'Irf8', 'Irf4', 'Sirpa', 'Icam1', 'Siglech',
                           'Adgre1', 'Bst2', 'Csf1r', 'Ly6c1', 'Itgam', 'Spn',
                           'Sell', 'Itgax', 'Fcgr2b'], save = '_markers4.pdf')
sc.pl.umap(comb, color = ['Xcr1', 'Clec9a', 'Fcer1a', 'Cd14', 'Cd8a', 'Ctsg', 'Fcgr4', 'H2-Aa', 'Ly6c1', 'Ly6c2'], save='PCA_markes5.pdf')
sc.pl.umap(comb, color = 'doublet_scores')


# %%
def plot_LKmarkers(adata, groupby='leiden', save=None):
  
    markers = {'HSC' : ['Procr', 'Mecom'],
             'Meg' : ['Pf4'],
             'Ery' : ['Klf1', 'Gata1', 'Hba-a1'],
             'Neu' : ['Elane', 'Prtn3', 'Ms4a3'],
             'Bas' : ['Prss34', 'Mcpt8'], 
             'Eos' : ['Prg2', 'Prg3', 'Epx'],
             'MC': ['Kit', 'Cma1', 'Gzmb'],
             'Mono/DC prog' : ['Csf1r', 'Mpo'],
             'Ly' : ['Il7r', 'Dntt'],
             'T cell' : ['Cd3e', 'Cd8a', 'Bcl11b'],
             'B cell' : ['Ly6d', 'Vpreb3', 'Cd79a'],
             'Ilc' : ['Gata3', 'Id2'],
             'Mono' : ['Cd14', 'Ctsg', 'Ly6c1', 'Ly6c2'],
             'pDC' : ['Irf8', 'Siglech', 'Ly6d'],
             'DCs' : ['H2-Aa', 'Xcr1', 'Itgax'],
             'Ifn-act' : ['Ifitm3', 'Irf7', 'Isg15'],
             'Myo C1' : ['C1qa', 'C1qb', 'Mpo']
             }

    genes = []
    for i in markers.values():
        genes = genes + i
    labels = list(markers.keys())

    positions = []
    p = 0
    for i in labels:
        n = len(markers[i]) - 1
        x = (p, n+p)
        p = p + n + 1
        positions.append(x)

    sc.pl.dotplot(adata,
               var_names = genes,
               groupby = groupby,
               var_group_positions = positions,
               var_group_labels = labels,
               save = save)
plot_LKmarkers(comb)

# %%
# # Detailed manual annotation - analysis excluding the 3d outlier
anno = {'0' : 'HSC', #Procr, Mecom
          '1' : 'Ery', #Klf1, hemoglobins
          '2' : 'Int prog',
          '3' : 'Int prog',
          '4' : 'Int prog',
          '5' : 'Int prog',
          '6' : 'Mono/DC prog', #Csf1r, Mpo
          '7' : 'Meg', #Pf4
          '8' : 'Int prog',
          '9' : 'Ery', #Klf1, hemoglobins
          '10' : 'Neu', #Prtn3, Elane, Ms4a3
          '11' : 'Ery', #Klf1, hemoglobins
          '12' : 'Bas/MC prog', #Ms4a2, tip has Mcpt8 (Bas marker), sugroup has Gzmb (MC marker)
          '13' : 'Mono/Doublets', #Cd14, Mpo, Ctsg, Ly6c1 and Ly6c2 but also different parts have all other lineages
          '14' : 'pDC', #Irf8, Siglech, Ly6d 
          '15' : 'T cell', #Cd3e, Gata3, Id2 
          '16' : 'Ly', #Dntt
          '17' : 'ILC', #Gata3, Id2, Il7r 
          '18' : 'HiMito', #Hi mitochondrial gene expression
          '19' : 'B cell', #Cd79a, Cd79b
          '20' : 'Ery', #hemoglobins
          '21' : 'Unknown/dublet', #Erythroid and Pf4 and Elane
          '22' : 'Myo C1', #Expressing complement genes and monocyte 
          '23' : 'Ifn-act prog.', #A series of IFN response gens: Ifitm3, Irf7, Isg15, Ifi203...
          '24' : 'DCs', #H2-Aa (MHC II),Xcr1, Specific for DCs: Itgax, 
          '25' : 'Eos', #Prg2, Prg3 
          '26' : 'Bas', #Mcpt8 and Prss34 
          '27' : 'T cell/ILC', #Cd3e, Gata3, Id2 
          '28' : 'B prog', #Dntt, Ly6d, Cd79a in the tip 
          '29' : 'ILC', #Id2, Gata3,
          '30' : 'Unknown'}

comb.obs['anno_man'] = [anno[i] if i in anno.keys() else 'unassigned' for i in comb.obs.leiden]

sc.pl.umap(comb, color = 'anno_man', save = 'PCA_anno_man.pdf')
utils.plots.umap3d(comb, color='anno_man', key='X_umap_3d', filename=base_figures+'umap3dPCA_anno_man.html')

# %% [markdown]
# ### Checkpoint

# %%
# We won't need the scaled data, removing to save space
comb.X = comb.layers['seurat'].copy()
del comb.layers['seurat']
comb.write(base_procdata + 'combined.h5ad', compression = 'lzf')

# %%
comb = sc.read(base_procdata + 'combined.h5ad')

# %% [markdown]
# ### Removing outlier clusters

# %%
tokeep=['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', 
        '14', '16', '20' ,'24', '25', '26' ,'28']
comb_filt = comb[comb.obs.index[comb.obs.leiden.isin(tokeep)],:].copy()
sc.pl.umap(comb_filt, color='anno_man')

sc.tl.paga(comb_filt, groups='leiden')
sc.pl.paga_compare(comb_filt, threshold = 0.05)
sc.pl.paga(comb_filt, threshold= 0.05)

# %% [markdown]
# ### UMAP on filtered data

# %%
#Getting the PAGA cluster position as initilisation coordinates for UMAP 
init_pos = np.zeros((comb_filt.n_obs, 2))
for n, i in enumerate(comb_filt.obs.leiden.cat.categories):
    init_pos[comb_filt.obs.leiden == i, :] = comb_filt.uns['paga']['pos'][n,:]
comb_filt.uns['umap_init_pos'] = init_pos.copy()

temp = cp.quick_umap(comb_filt, use_rep='X_pca_harmony', init_coords = init_pos)
unclog_umap_caching(temp)
umapref = cp.quick_umap(comb_filt, use_rep='X_pca_harmony', init_coords = init_pos)
comb_filt.obsm['X_umap_2d'] = comb_filt.obsm['X_umap'].copy()

#Saving umapref for future use
with open(base_procdata + 'combined_filt_umapref.pkl', 'wb') as f:
    pickle.dump(umapref, f)

temp3d = cp.quick_umap(comb_filt, n_components=3, use_rep='X_pca_harmony')
unclog_umap_caching(temp3d)
umapref_3d = cp.quick_umap(comb_filt, n_components=3, use_rep='X_pca_harmony')

#Saving umapref for future use
with open(base_procdata + 'combined_filt_umapref3d.pkl', 'wb') as f:
    pickle.dump(umapref_3d, f)

comb_filt.obsm['X_umap_3d'] = comb_filt.obsm['X_umap'].copy()
comb_filt.obsm['X_umap'] = comb_filt.obsm['X_umap_2d'].copy()

comb_filt.uns['leiden_colors'] = ['#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b',
       '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a',
       '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5',
       '#ad494a', '#8c6d31']

sc.pl.umap(comb_filt, color='leiden', legend_loc='on data', save = 'PCA_filt_leiden.pdf')
sc.pl.umap(comb_filt, color='anno_man', save = 'PCA_filt_annoman.pdf')
sc.pl.paga_compare(comb_filt, threshold = 0.05, save = 'PCA_filt.pdf')

# %%
comb_filt.uns['pagaPCA'] = comb_filt.uns['paga'].copy()

# %%
# comb_filt.write(base_procdata + 'combined_filt.h5ad', compression = 'lzf')

# %% [markdown]
# ## Diffusion pseudotime

# %% [markdown]
# ### HSC score
# First we compute the HSC score for each cell to identify the most immature cells and set the root for DPT

# %%
comb_counts = comb_filt.raw.to_adata()
comb_counts = utils.proc.delognorm(comb_counts)
comb_counts.write(base_procdata + 'combined_filt_counts.h5ad', compression='lzf')

# %%
HSCscore_genes = pd.Series(np.genfromtxt('data/HSCscore/model_molo_genes.txt', dtype='str'))
missing = HSCscore_genes[~HSCscore_genes.isin(comb_counts.var.index)]

tempX = np.concatenate((comb_counts.X.toarray(), np.zeros((comb_counts.n_obs, len(missing)))),
               axis = 1)
temp = AnnData(X = tempX,
              obs = comb_counts.obs,
              var = pd.DataFrame(index = np.concatenate((comb_counts.var.index.values, missing.values))))

comb_filt.obs['HSCscore'] = calculate_HSCscore(temp)
sc.pl.umap(comb_filt, color = 'HSCscore', save = '_HSCscore.pdf')


#To choose the root cell we need to find the cell with the highest HSC score.
#Data is quite noisy so we will use the nearest neighbor information to find cells
#with the highest HSCscore among the nearest neighbors (including the cells itself)
x = comb_filt.obsp['connectivities']
inds = x.tolil().rows #Simple way to get the indices for nearest neighbors. lil is a lis of list format.

y = comb_filt.obs.HSCscore.values
nn_HSCscores = [y[i + [n]].mean() for n, i in enumerate(inds)] #indices are missing the cell itself so adding it back
comb_filt.obs['nn_HSCscore'] = nn_HSCscores

sc.pl.umap(comb_filt, color = 'nn_HSCscore', save = '_nnHSCscore.pdf')

# %% [markdown]
# ### DPT

# %%
root = np.where(comb_filt.obs.nn_HSCscore == comb_filt.obs.nn_HSCscore.max())[0][0]
# #Manually specifying root
# root = np.where(comb_clean.obs.index == '')[0][0]

comb_filt.uns['iroot'] = root
comb_filt.obs['isroot'] = False
comb_filt.obs.isroot[root] = True
comb_filt.obs['isroot'] = pd.Categorical(comb_filt.obs.isroot)
sc.pl.umap(comb_filt, color = 'isroot')
utils.plots.umap3d(comb_filt, color = 'isroot', key = 'X_umap_3d', filename = base_figures+ 'rootcheck.html')

sc.tl.dpt(comb_filt, n_dcs = 15)

# %%
comb_filt.write(base_procdata + 'combined_filt.h5ad', compression = 'lzf')

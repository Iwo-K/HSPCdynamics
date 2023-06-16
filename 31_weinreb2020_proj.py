# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.13.8
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Cell projections Reizis

# %% [markdown]
# ## Setup

# %%
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import pathlib
import re
import anndata
import utils.proc, utils.plots
import cellproject as cp
import pickle
from utils.fixes import unclog_umap_caching

import warnings
warnings.filterwarnings("ignore", category=FutureWarning) 

sc.settings.verbosity = 3
base_figures = './figures/07script/'
sc.settings.figdir = base_figures
base_procdata = './procdata/07script/'
for i in [base_figures, base_procdata]:
    pathlib.Path(i).mkdir(parents=True, exist_ok=True)

# %%
#Figure settings
sc.set_figure_params(color_map = 'viridis', dpi_save = 350, frameon=False)
mpl.rcParams['figure.figsize'] = (4, 4) #1:1 plot ratio
mpl.rc('xtick', labelsize=14) 
mpl.rc('ytick', labelsize=14) 
mpl.rc('axes', labelsize=16) 
mpl.rc('axes', labelsize=16)
mpl.rcParams['pdf.fonttype'] = 42 #Ensures readable fonts in illustrator
mpl.rcParams['ps.fonttype'] = 42
plt.rc('axes', axisbelow=True) #Ensure that gridlines are behind points
mplparams = mpl.rcParams.copy()

# %%
# %load_ext autoreload
# %autoreload 2

# %%
hoxb5 = sc.read('procdata/04script/combined_filt.h5ad')
meta = pd.read_csv('procdata/04script/combined_sample.meta.csv', index_col = 0)
meta = meta.sort_values(by ='timepoint_tx_days', ascending = False)
hoxb5.obs['longname'] = [meta.loc[i, 'longname'] for i in hoxb5.obs.biosample_id]

# sc.pp.subsample(hoxb5, n_obs=5000)
#Getting the UMAP object
with open('procdata/04script/combined_filt_umapref.pkl', 'rb') as f:
    umapref = pickle.load(f)

sc.pl.umap(hoxb5)
hoxb5_pos = hoxb5[hoxb5.obs.tom == 'pos'].copy()
hoxb5 = hoxb5[hoxb5.obs.data_type=='10x',:].copy()

# %% [markdown]
# ## Fates

# %%
for clu in [6, 10, 11, 16, 7]:
    table = pd.read_csv( './pseudodynamics/tables/table_all_parameters_clu_' + str(clu) + '.csv',index_col=0)
    #table = table.sort_values(by = 'dpt_pseudotime')
    hoxb5_pos.obs['fate_' + str(clu)] = hoxb5_pos.obs_names.isin(table.index)
    hoxb5_pos.obs['dpt_fate' + str(clu)] = table.dpt_pseudotime_scaled

# %%
hoxb5_pos = hoxb5_pos[hoxb5_pos.obs.data_type == '10x'].copy()

# %%
hoxb5_pos.obs['fate_6'] = hoxb5_pos.obs.fate_6.astype('category')
hoxb5_pos.obs['fate_7'] = hoxb5_pos.obs.fate_7.astype('category')
hoxb5_pos.obs['fate_10'] = hoxb5_pos.obs.fate_10.astype('category')
hoxb5_pos.obs['fate_11'] = hoxb5_pos.obs.fate_11.astype('category')
hoxb5_pos.obs['fate_16'] = hoxb5_pos.obs.fate_16.astype('category')

# %% [markdown]
# # Weinreb projections and abundance

# %%
klein = sc.read('klein2020/procdata/01script/klein2020_dcounts.h5ad')

# %%
klein.obs.Time_point.value_counts()

# %%
sc.pp.filter_cells(klein,min_counts=0)
sc.pp.filter_genes(klein,min_cells=0)
klein.var['mt'] = klein.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(klein, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# %%
# klein_basic = klein.copy()

# %% [markdown]
# ## Day 2

# %%
klein2 = klein[klein.obs.Time_point == 2.0].copy()

# %%
klein2 = klein2[klein2.obs.n_counts > 1800, :].copy()# Iwo day 2 threshold

# %%
utils.proc.lognorm(klein2)

# %%
klein2_basic = klein2.copy()

# %%
ref = hoxb5.copy()
#Unifying genes
ref = ref[:, ref.var.index[ref.var.index.isin(klein2.var.index)]].copy()
ref.X = ref.raw[:,ref.var.index].X.copy()
del ref.raw
klein2 = klein2[:, ref.var.index].copy()
#Combining data
comb = ref.concatenate(klein2, batch_key='batch', batch_categories=['ref', 'klein'], index_unique=None)
comb.var['highly_variable'] = comb.var['highly_variable-ref']
#The harmony integration does not work at all, scanorama does not work either, trying Seurat

# %%
comb.obs.batch.value_counts()

# %%
comb_cor = cp.run_SeuratCCA(comb, batch_key='batch', reference='ref', k_filter=10,debug=True) crushes with 20 nodes, ok with 24

# %%
comb_cor[0].write('./procdata/comb_temp2d.h5')

# %%
comb_cor = sc.read('./procdata/comb_temp2da.h5')

# %%
#We will need PCs based on the common space
sc.pp.scale(ref)
sc.pp.pca(ref, n_comps=50)
sc.pp.neighbors(ref, n_neighbors=15)

#Setting the Seurat-corrected values for the target data (ref data is already scaled)
klein2.X = comb_cor[klein2.obs.index, klein2.var.index].X.copy()

# %%
#Projecting into the ref PC space and identifying neighbors
project_cells(klein2, ref,
                 obs_columns=['leiden','dpt_pseudotime'],
                 fit_pca=True,
                 scale_data=True)

#Swapping for the original PCs, based on cell harmony-corrected data
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()
#nn-regressing the PCs
cp.nnregress(klein2,
             ref,
             regress=['pca'],
             weighted=True)

#fitting into the original umap
project_cells(klein2, ref,
                 obs_columns=['leiden','dpt_pseudotime'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)


# %%
klein2.obs.head()

# %%
sc.pl.umap(klein2, color = ['ref_dpt_pseudotime','ref_leiden'],size = 20)

# %%
klein2_basic.obs = pd.concat([klein2_basic.obs,klein2.obs.ref_leiden],axis=1)

# %%
klein2_basic.obs = pd.concat([klein2_basic.obs,klein2.obs.ref_dpt_pseudotime],axis=1)

# %%
klein2_basic.write(base_procdata + 'Klein2_2020_hoxb5projection.h5ad', compression='lzf')

# %%
klein2_basic.obsm = klein2.obsm.copy()

# %%
klein2 = sc.read(base_procdata + 'Klein2_2020_hoxb5projection.h5ad')

# %%
ref = hoxb5.copy()

# %%
#blu dots figure 
comb_plot = ref.concatenate(klein2, batch_key='batch', batch_categories=['ref', 'klein2'], index_unique=None)

# %%
sc.pl.umap(comb_plot, color = ['ref_dpt_pseudotime','ref_leiden'],size = 20)

# %%

# %%
# Counting cells in each cluster seprately for Tom+ and Tom- cells
pop = comb_plot[comb_plot.obs.batch == 'klein2',:].copy()

pop_leiden = pd.crosstab(pop.obs.ref_leiden, pop.obs.Time_point)
pop_leiden.to_csv(base_procdata + 'klein2_cluster_cellcounts.csv')

# %%
pop_leiden.sum()

# %%
pop_leiden

# %% [markdown]
# # Project on label

# %%
klein2 = klein[klein.obs.Time_point == 2.0].copy()

# %%
klein2 = klein2[klein2.obs.n_counts > 1800, :].copy()# Iwo day 2 threshold

# %%
utils.proc.lognorm(klein2)

# %%
ref = hoxb5_pos.copy()
#Unifying genes
ref = ref[:, ref.var.index[ref.var.index.isin(klein2.var.index)]].copy()
ref.X = ref.raw[:,ref.var.index].X.copy()

del ref.raw
klein2 = klein2[:, ref.var.index].copy()

comb = ref.concatenate(klein2, batch_key='batch', batch_categories=['ref', 'klein2'], index_unique=None)
comb.var['highly_variable'] = comb.var['highly_variable-ref']

# %%
comb.obs.batch.value_counts()

# %%
#We will need PCs based on the common space
sc.pp.scale(ref)
sc.pp.pca(ref, n_comps=50)
sc.pp.neighbors(ref, n_neighbors=15)

#Setting the Seurat-corrected values for the target data (ref data is already scaled)
klein2.X = comb_cor[klein2.obs.index, klein2.var.index].X.copy()

# %%
cp.project_cells(klein2, ref,
                 obs_columns=['fate_6','fate_7','fate_10','fate_11','fate_16'],
                 fit_pca=True,
                 scale_data=True)

#Swapping for the original PCs, based on cell harmony-corrected data
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()
#nn-regressing the PCs
cp.nnregress(klein2,
             ref,
             regress=['pca'],
             weighted=True)

#fitting into the original umap
cp.project_cells(klein2, ref,
                 obs_columns=['fate_6','fate_7','fate_10','fate_11','fate_16'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)

# %%
sc.pl.umap(klein2,color=['ref_fate_6','ref_fate_7','ref_fate_10','ref_fate_11','ref_fate_16'],size = 10)

# %%
klein2_basic.obsm = klein2.obsm.copy()

# %%
klein2_basic.obs = pd.concat([klein2_basic.obs,klein2.obs.iloc[:,-5:]],axis=1)

# %%
klein2_basic.obs.head()

# %%

# %%
sc.pl.umap(klein2_basic,color=['ref_fate_6','ref_fate_7','ref_fate_10','ref_fate_11','ref_fate_16'],size = 10)

# %%
#klein2_basic.write(base_procdata + 'Klein2_2020_hoxb5projection.h5ad', compression='lzf')

# %% [markdown]
# # Project dpt fate

# %%
for clu in [6,7,10,11,16]:

    klein2_temp = klein2[klein2.obs['ref_fate_'+str(clu)] == True].copy()
    ref_temp = ref[ref.obs['fate_'+str(clu)] == True].copy()

    project_cells(klein2_temp, ref_temp,
                 obs_columns=['dpt_fate'+str(clu)],
                 fit_pca=False,
                 scale_data=False)
    
    klein2.obs['ref_dpt_fate'+str(clu)] = klein2_temp.obs['ref_dpt_fate'+str(clu)].copy()


# %%
sc.pl.umap(klein2,color=['ref_dpt_fate6','ref_dpt_fate7','ref_dpt_fate10','ref_dpt_fate11','ref_dpt_fate16'],size = 30)

# %%
ref = hoxb5.copy()

# %%
comb_plot = ref.concatenate(klein2, batch_key='batch', batch_categories=['ref', 'klein2'], index_unique=None)

# %%
sc.pl.umap(comb_plot,color=['ref_dpt_fate6','ref_dpt_fate11','ref_dpt_fate16'],size = 30)

# %%

# %%

# %%
klein2_basic.obs = pd.concat([klein2_basic.obs,klein2.obs.iloc[:,-5:]],axis=1)

# %%
klein2_basic.write(base_procdata + 'Klein2_2020_hoxb5projection.h5ad', compression='lzf')

# %%
sc.pl.umap(klein2_basic,color=['ref_dpt_fate6','ref_dpt_fate7','ref_dpt_fate10','ref_dpt_fate11','ref_dpt_fate16'],size = 30)

# %% [markdown]
# # Day 4

# %%
klein4 = klein[klein.obs.Time_point == 4.0].copy()

# %%
klein4 = klein4[klein4.obs.n_counts > 1800, :].copy()# Iwo day 2 threshold

# %%
klein4

# %%
utils.proc.lognorm(klein4)

# %%
klein4_basic = klein4.copy()

# %%
ref = hoxb5.copy()
#Unifying genes
ref = ref[:, ref.var.index[ref.var.index.isin(klein4.var.index)]].copy()
ref.X = ref.raw[:,ref.var.index].X.copy()
del ref.raw
klein4 = klein4[:, ref.var.index].copy()
#Combining data
comb = ref.concatenate(klein4, batch_key='batch', batch_categories=['ref', 'klein'], index_unique=None)
comb.var['highly_variable'] = comb.var['highly_variable-ref']

#The harmony integration does not work at all, scanorama does not work either, trying Seurat

# %%
comb.obs.batch.value_counts()

# %%
#Running the batch correction with Seurat

comb_cor = cp.run_SeuratCCA(comb, batch_key='batch', reference='ref', k_filter=10,debug=True)

# %%
comb_cor[0].write('comb_temp4d.h5')

# %%
comb_cor = sc.read('./procdata/comb_temp4d.h5')

# %%
comb_cor.obs.batch.value_counts()

# %%
#We will need PCs based on the common space
sc.pp.scale(ref)
sc.pp.pca(ref,n_comps=50)
sc.pp.neighbors(ref, n_neighbors=15)

#Setting the Seurat-corrected values for the target data (ref data is already scaled)
klein4.X = comb_cor[klein4.obs.index, klein4.var.index].X.copy()

# %%
#Projecting into the ref PC space and identifying neighbors
project_cells(klein4, ref,
                 obs_columns=['leiden','dpt_pseudotime'],
                 fit_pca=True,
                 scale_data=True)

#Swapping for the original PCs, based on cell harmony-corrected data
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()
#nn-regressing the PCs
cp.nnregress(klein4,
             ref,
             regress=['pca'],
             weighted=True)

#fitting into the original umap
project_cells(klein4, ref,
                 obs_columns=['leiden','dpt_pseudotime'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)


# %%
sc.pl.umap(klein4, color = ['ref_dpt_pseudotime','ref_leiden'],size = 20)

# %%

# %%
klein4_basic.obs = pd.concat([klein4_basic.obs,klein4.obs.ref_leiden],axis=1)

# %%
klein4_basic.obs = pd.concat([klein4_basic.obs,klein4.obs.ref_dpt_pseudotime],axis=1)

# %%
klein4_basic.write(base_procdata + 'Klein4_2020_hoxb5projection.h5ad', compression='lzf')

# %%
klein4_basic.obsm = klein4.obsm.copy()

# %%

# %%
ref = hoxb5.copy()

# %%
comb_plot = ref.concatenate(klein4, batch_key='batch', batch_categories=['ref', 'klein4'], index_unique=None)

# %%
sc.pl.umap(comb_plot, color = ['ref_dpt_pseudotime','ref_leiden'],size = 20)

# %%
# Counting cells in each cluster seprately for Tom+ and Tom- cells
pop = comb_plot[comb_plot.obs.batch == 'klein4',:].copy()

pop_leiden = pd.crosstab(pop.obs.ref_leiden, pop.obs.Time_point)
pop_leiden.to_csv(base_procdata + 'klein4_cluster_cellcounts.csv')

# %%
pop_leiden.sum()

# %%
pop_leiden

# %% [markdown]
# # Project on label

# %%
klein4 = klein[klein.obs.Time_point == 4.0].copy()

# %%
klein4 = klein4[klein4.obs.n_counts > 1800, :].copy()# Iwo day 2 threshold

# %%
utils.proc.lognorm(klein4)

# %%
klein4

# %%
ref = hoxb5_pos.copy()
#Unifying genes
ref = ref[:, ref.var.index[ref.var.index.isin(klein4.var.index)]].copy()
ref.X = ref.raw[:,ref.var.index].X.copy()

del ref.raw
klein4 = klein4[:, ref.var.index].copy()

comb = ref.concatenate(klein4, batch_key='batch', batch_categories=['ref', 'klein4'], index_unique=None)
comb.var['highly_variable'] = comb.var['highly_variable-ref']

# %%
comb.obs.batch.value_counts()

# %%
#We will need PCs based on the common space
sc.pp.scale(ref)
sc.pp.pca(ref, n_comps=50)
sc.pp.neighbors(ref, n_neighbors=15)

#Setting the Seurat-corrected values for the target data (ref data is already scaled)
klein4.X = comb_cor[klein4.obs.index, klein4.var.index].X.copy()

# %%
cp.project_cells(klein4, ref,
                 obs_columns=['fate_6','fate_7','fate_10','fate_11','fate_16'],
                 fit_pca=True,
                 scale_data=True)

#Swapping for the original PCs, based on cell harmony-corrected data
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()
#nn-regressing the PCs
cp.nnregress(klein4,
             ref,
             regress=['pca'],
             weighted=True)

#fitting into the original umap
cp.project_cells(klein4, ref,
                 obs_columns=['fate_6','fate_7','fate_10','fate_11','fate_16'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)

# %%
sc.pl.umap(klein4,color=['ref_fate_6','ref_fate_7','ref_fate_10','ref_fate_11','ref_fate_16'],size = 10)

# %%
klein4_basic.obsm = klein4.obsm.copy()

# %%
klein4_basic.obs = pd.concat([klein4_basic.obs,klein4.obs.iloc[:,-5:]],axis=1)

# %%
klein4_basic.obs.head()

# %%
sc.pl.umap(klein4,color=['ref_fate_6','ref_fate_7','ref_fate_10','ref_fate_11','ref_fate_16'],size = 10)

# %%
sc.pl.umap(klein4_basic,color=['ref_fate_6','ref_fate_7','ref_fate_10','ref_fate_11','ref_fate_16'],size = 10)

# %%
sc.pl.umap(klein4_basic,color=['ref_leiden','ref_dpt_pseudotime'],size = 10)

# %%
klein4_basic.write(base_procdata + 'Klein4_2020_hoxb5projection.h5ad', compression='lzf')

# %% [markdown]
# # Project dpt fate

# %%
for clu in [6,7,10,11,16]:

    klein4_temp = klein4[klein4.obs['ref_fate_'+str(clu)] == True].copy()
    ref_temp = ref[ref.obs['fate_'+str(clu)] == True].copy()

    project_cells(klein4_temp, ref_temp,
                 obs_columns=['dpt_fate'+str(clu)],
                 fit_pca=False,
                 scale_data=False)
    
    klein4.obs['ref_dpt_fate'+str(clu)] = klein4_temp.obs['ref_dpt_fate'+str(clu)].copy()


# %%
sc.pl.umap(klein4,color=['ref_dpt_fate6','ref_dpt_fate7','ref_dpt_fate10','ref_dpt_fate11','ref_dpt_fate16'],size = 30)

# %%
ref = hoxb5.copy()

# %%
comb_plot = ref.concatenate(klein4, batch_key='batch', batch_categories=['ref', 'klein4'], index_unique=None)

# %%
sc.pl.umap(comb_plot,color=['ref_dpt_fate6','ref_dpt_fate11','ref_dpt_fate16'],size = 30)

# %%

# %%

# %%
klein4_basic.obs = pd.concat([klein4_basic.obs,klein4.obs.iloc[:,-5:]],axis=1)

# %%
klein4_basic.write(base_procdata + 'Klein4_2020_hoxb5projection.h5ad', compression='lzf')

# %%
sc.pl.umap(klein4_basic,color=['ref_dpt_fate6','ref_dpt_fate7','ref_dpt_fate10','ref_dpt_fate11','ref_dpt_fate16'],size = 30)

# %% [markdown]
# # Day 6

# %%
klein6 = klein[klein.obs.Time_point == 6.0].copy()

# %%
klein6 = klein6[klein6.obs.n_counts > 1800, :].copy()# Iwo day 2 threshold

# %%
klein6

# %%
utils.proc.lognorm(klein6)

# %%
klein6_basic = klein6.copy()

# %%
ref = hoxb5.copy()
#Unifying genes
ref = ref[:, ref.var.index[ref.var.index.isin(klein6.var.index)]].copy()
ref.X = ref.raw[:,ref.var.index].X.copy()
del ref.raw
klein6 = klein6[:, ref.var.index].copy()
#Combining data
comb = ref.concatenate(klein6, batch_key='batch', batch_categories=['ref', 'klein'], index_unique=None)
comb.var['highly_variable'] = comb.var['highly_variable-ref']

#The harmony integration does not work at all, scanorama does not work either, trying Seurat

# %%
comb.obs.batch.value_counts()

# %%
##Running the batch correction with Seurat

comb_cor = cp.run_SeuratCCA(comb, batch_key='batch', reference='ref', k_filter=10,debug=True)

# %%
comb_cor[0].write('comb_temp6d.h5')

# %%
comb_cor = sc.read('./procdata/comb_temp6d.h5')

# %%
comb_cor.obs.batch.value_counts()

# %%
#We will need PCs based on the common space
sc.pp.scale(ref)
sc.pp.pca(ref,n_comps=50)
sc.pp.neighbors(ref, n_neighbors=15)

#Setting the Seurat-corrected values for the target data (ref data is already scaled)
klein6.X = comb_cor[klein6.obs.index, klein6.var.index].X.copy()

# %%
#Projecting into the ref PC space and identifying neighbors
project_cells(klein6, ref,
                 obs_columns=['leiden','dpt_pseudotime'],
                 fit_pca=True,
                 scale_data=True)

#Swapping for the original PCs, based on cell harmony-corrected data
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()
#nn-regressing the PCs
cp.nnregress(klein6,
             ref,
             regress=['pca'],
             weighted=True)

#fitting into the original umap
project_cells(klein6, ref,
                 obs_columns=['leiden','dpt_pseudotime'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)


# %%
sc.pl.umap(klein6, color = ['ref_dpt_pseudotime','ref_leiden'],size = 20)

# %%

# %%
klein6_basic.obs = pd.concat([klein6_basic.obs,klein6.obs.ref_leiden],axis=1)

# %%
klein6_basic.obs = pd.concat([klein6_basic.obs,klein6.obs.ref_dpt_pseudotime],axis=1)

# %%

# %%
klein6_basic.write(base_procdata + 'Klein6_2020_hoxb5projection.h5ad', compression='lzf')

# %%
klein6_basic.obsm = klein6.obsm.copy()

# %%

# %%
ref = hoxb5.copy()

# %%
comb_plot = ref.concatenate(klein6, batch_key='batch', batch_categories=['ref', 'klein6'], index_unique=None)

# %%
sc.pl.umap(comb_plot, color = ['ref_dpt_pseudotime','ref_leiden'],size = 20)

# %%
# Counting cells in each cluster seprately for Tom+ and Tom- cells
pop = comb_plot[comb_plot.obs.batch == 'klein6',:].copy()

pop_leiden = pd.crosstab(pop.obs.ref_leiden, pop.obs.Time_point)
pop_leiden.to_csv(base_procdata + 'klein6_cluster_cellcounts.csv')

# %%
pop_leiden.sum()

# %%
pop_leiden

# %% [markdown]
# # Project on label

# %%
klein6 = klein[klein.obs.Time_point == 6.0].copy()

# %%
klein6 = klein6[klein6.obs.n_counts > 1800, :].copy()# Iwo day 2 threshold

# %%
utils.proc.lognorm(klein6)

# %%
klein6

# %%
ref = hoxb5_pos.copy()
#Unifying genes
ref = ref[:, ref.var.index[ref.var.index.isin(klein6.var.index)]].copy()
ref.X = ref.raw[:,ref.var.index].X.copy()

del ref.raw
klein6 = klein6[:, ref.var.index].copy()

comb = ref.concatenate(klein6, batch_key='batch', batch_categories=['ref', 'klein6'], index_unique=None)
comb.var['highly_variable'] = comb.var['highly_variable-ref']

# %%
comb.obs.batch.value_counts()

# %%
#We will need PCs based on the common space
sc.pp.scale(ref)
sc.pp.pca(ref, n_comps=50)
sc.pp.neighbors(ref, n_neighbors=15)

#Setting the Seurat-corrected values for the target data (ref data is already scaled)
klein6.X = comb_cor[klein6.obs.index, klein6.var.index].X.copy()

# %%
cp.project_cells(klein6, ref,
                 obs_columns=['fate_6','fate_7','fate_10','fate_11','fate_16'],
                 fit_pca=True,
                 scale_data=True)

#Swapping for the original PCs, based on cell harmony-corrected data
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()
#nn-regressing the PCs
cp.nnregress(klein6,
             ref,
             regress=['pca'],
             weighted=True)

#fitting into the original umap
cp.project_cells(klein6, ref,
                 obs_columns=['fate_6','fate_7','fate_10','fate_11','fate_16'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)

# %%
sc.pl.umap(klein6,color=['ref_fate_6','ref_fate_7','ref_fate_10','ref_fate_11','ref_fate_16'],size = 10)

# %%

# %%
klein6_basic.obsm = klein6.obsm.copy()

# %%
klein6_basic.obs = pd.concat([klein6_basic.obs,klein6.obs.iloc[:,-5:]],axis=1)

# %%
klein6_basic.obs.head()

# %%
sc.pl.umap(klein6,color=['ref_fate_6','ref_fate_7','ref_fate_10','ref_fate_11','ref_fate_16'],size = 10)

# %%
sc.pl.umap(klein6_basic,color=['ref_fate_6','ref_fate_7','ref_fate_10','ref_fate_11','ref_fate_16'],size = 10)

# %%
sc.pl.umap(klein6_basic,color=['ref_leiden','ref_dpt_pseudotime'],size = 10)

# %%
klein6_basic.write(base_procdata + 'Klein6_2020_hoxb5projection.h5ad', compression='lzf')

# %% [markdown]
# # Project dpt fate

# %%
for clu in [6,7,10,11,16]:

    klein6_temp = klein6[klein6.obs['ref_fate_'+str(clu)] == True].copy()
    ref_temp = ref[ref.obs['fate_'+str(clu)] == True].copy()

    project_cells_mela(klein6_temp, ref_temp,
                 obs_columns=['dpt_fate'+str(clu)],
                 fit_pca=False,
                 scale_data=False)
    
    klein6.obs['ref_dpt_fate'+str(clu)] = klein6_temp.obs['ref_dpt_fate'+str(clu)].copy()


# %%
sc.pl.umap(klein6,color=['ref_dpt_fate6','ref_dpt_fate7','ref_dpt_fate10','ref_dpt_fate11','ref_dpt_fate16'],size = 30)

# %%

# %%
ref = hoxb5.copy()

# %%
comb_plot = ref.concatenate(klein6, batch_key='batch', batch_categories=['ref', 'klein6'], index_unique=None)

# %%
sc.pl.umap(comb_plot,color=['ref_dpt_fate6','ref_dpt_fate11','ref_dpt_fate16'],size = 30)

# %%

# %%

# %%
klein6_basic.obs = pd.concat([klein6_basic.obs,klein6.obs.iloc[:,-5:]],axis=1)

# %%
klein6_basic.write(base_procdata + 'Klein6_2020_hoxb5projection.h5ad', compression='lzf')

# %%
sc.pl.umap(klein6_basic,color=['ref_dpt_fate6','ref_dpt_fate7','ref_dpt_fate10','ref_dpt_fate11','ref_dpt_fate16'],size = 30)

# %%

# %%

# %%

# %%

# %%

# %%

# %% [markdown]
# # Create input PD

# %%
klein2 = sc.read(base_procdata + 'Klein2_2020_hoxb5projection.h5ad')
klein4 = sc.read(base_procdata + 'Klein4_2020_hoxb5projection.h5ad')
klein6 = sc.read(base_procdata + 'Klein6_2020_hoxb5projection.h5ad')

# %%
sc.pl.umap(klein2,ncols=5,color=['ref_dpt_fate6','ref_dpt_fate7','ref_dpt_fate10','ref_dpt_fate11','ref_dpt_fate16'],size = 30)
sc.pl.umap(klein4,ncols=5,color=['ref_dpt_fate6','ref_dpt_fate7','ref_dpt_fate10','ref_dpt_fate11','ref_dpt_fate16'],size = 30)
sc.pl.umap(klein6,ncols=5,color=['ref_dpt_fate6','ref_dpt_fate7','ref_dpt_fate10','ref_dpt_fate11','ref_dpt_fate16'],size = 30)

# %%

# %%
plt.hist(klein2[klein2.obs.ref_fate_6 == 'True'].obs.ref_dpt_fate6,bins = 100);
plt.hist(klein4[klein4.obs.ref_fate_6 == 'True'].obs.ref_dpt_fate6,bins = 100);
plt.hist(klein6[klein6.obs.ref_fate_6 == 'True'].obs.ref_dpt_fate6,bins = 100);


# %%
plt.hist(klein2[klein2.obs.ref_fate_7 == 'True'].obs.ref_dpt_fate7,bins = 100);
plt.hist(klein4[klein4.obs.ref_fate_7 == 'True'].obs.ref_dpt_fate7,bins = 100);
plt.hist(klein6[klein6.obs.ref_fate_7 == 'True'].obs.ref_dpt_fate7,bins = 100);

# %%
plt.hist(klein2[klein2.obs.ref_fate_10 == 'True'].obs.ref_dpt_fate10,bins = 100);
plt.hist(klein4[klein4.obs.ref_fate_10 == 'True'].obs.ref_dpt_fate10,bins = 100);
plt.hist(klein6[klein6.obs.ref_fate_10 == 'True'].obs.ref_dpt_fate10,bins = 100);

# %%
plt.hist(klein2[klein2.obs.ref_fate_11 == 'True'].obs.ref_dpt_fate11,bins = 100);
plt.hist(klein4[klein4.obs.ref_fate_11 == 'True'].obs.ref_dpt_fate11,bins = 100);
plt.hist(klein6[klein6.obs.ref_fate_11 == 'True'].obs.ref_dpt_fate11,bins = 100);

# %%
plt.hist(klein2[klein2.obs.ref_fate_16 == 'True'].obs.ref_dpt_fate16,bins = 100);
plt.hist(klein4[klein4.obs.ref_fate_16 == 'True'].obs.ref_dpt_fate16,bins = 100);
plt.hist(klein6[klein6.obs.ref_fate_16 == 'True'].obs.ref_dpt_fate16,bins = 100);

# %%
klein

# %%
klein = klein2.concatenate(klein4,klein6)

# %%
klein.obs.batch.value_counts()

# %%
df = pd.DataFrame(list(klein.obs.Time_point),columns=['stage'],index = klein.obs_names)

# %%
df['sample'] = klein.obs.Well

# %%
df

# %%
for clu in [6, 10, 11, 16, 7]:
    df_temp = df.loc[klein.obs['ref_fate_' + str(clu)] == 'True'].copy()
    df_temp['dpt_pseudotime'] = klein.obs['ref_dpt_fate' + str(clu)]
    df_temp.to_csv(base_procdata + 'klein_input_pseudo_dyn_clu_' + str(clu) + '_dpt.csv')

# %%
df_temp

# %%
plt.hist(klein2[klein2.obs.ref_fate_16 == 'True'].obs.ref_dpt_fate16,bins = 50,density=True);
plt.hist(klein4[klein4.obs.ref_fate_16 == 'True'].obs.ref_dpt_fate16,bins = 50,density=True);
plt.hist(klein6[klein6.obs.ref_fate_16 == 'True'].obs.ref_dpt_fate16,bins = 50,density=True);

# %%
plt.hist(klein2[klein2.obs.ref_fate_7 == 'True'].obs.ref_dpt_fate7,bins = 50,density=True);
plt.hist(klein4[klein4.obs.ref_fate_7 == 'True'].obs.ref_dpt_fate7,bins = 50,density=True);
plt.hist(klein6[klein6.obs.ref_fate_7 == 'True'].obs.ref_dpt_fate7,bins = 50,density=True);

# %%
plt.hist(klein2[klein2.obs.ref_fate_11 == 'True'].obs.ref_dpt_fate11,bins = 50,density=True);
plt.hist(klein4[klein4.obs.ref_fate_11 == 'True'].obs.ref_dpt_fate11,bins = 50,density=True);
plt.hist(klein6[klein6.obs.ref_fate_11 == 'True'].obs.ref_dpt_fate11,bins = 50,density=True);

# %%
plt.hist(klein2[klein2.obs.ref_fate_10 == 'True'].obs.ref_dpt_fate10,bins = 50,density=True);
plt.hist(klein4[klein4.obs.ref_fate_10 == 'True'].obs.ref_dpt_fate10,bins = 50,density=True);
plt.hist(klein6[klein6.obs.ref_fate_10 == 'True'].obs.ref_dpt_fate10,bins = 50,density=True);

# %%
plt.hist(klein2[klein2.obs.ref_fate_6 == 'True'].obs.ref_dpt_fate6,bins = 50,density=True);
plt.hist(klein4[klein4.obs.ref_fate_6 == 'True'].obs.ref_dpt_fate6,bins = 50,density=True);
plt.hist(klein6[klein6.obs.ref_fate_6 == 'True'].obs.ref_dpt_fate6,bins = 50,density=True);

# %%

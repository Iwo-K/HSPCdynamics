# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %%
#Singularity container used:
# !echo $SINGULARITY_CONTAINER
# !hostname

# %% [markdown]
# # Figure 2 - cell projections

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
hoxb5 = hoxb5[hoxb5.obs.data_type=='10x',:].copy()

# %% [markdown] tags=[]
# ## Nestorowa et al. 2016 data projection

# %%
#Loading Nestorowa et al. 2016 data
sfdata = sc.read('./data/ext_data/sfdata/sfdata_nlog.h5ad')
sfdata.raw = sfdata

# %% [markdown]
# ### Projecting data

# %%
ref = hoxb5.copy()
#Unifying genes
ref = ref[:, ref.var.index[ref.var.index.isin(sfdata.var.index)]].copy()
ref.X = ref.raw[:,ref.var.index].X.copy()
del ref.raw
sfdata = sfdata[:, ref.var.index].copy()
#Combining data
comb = ref.concatenate(sfdata, batch_key='batch', batch_categories=['ref', 'sfdata'], index_unique=None)
comb.var['highly_variable'] = comb.var['highly_variable-ref']

#Running the batch correction with Seurat
comb_cor = cp.run_SeuratCCA(comb, batch_key='batch', reference='ref')

# %%
#We will need PCs based on the common space
sc.pp.scale(ref)
sc.pp.pca(ref, n_comps=50)
sc.pp.neighbors(ref, n_neighbors=15)

#Setting the Seurat-corrected values for the target data (ref data is already scaled)
sfdata.X = comb_cor[sfdata.obs.index, sfdata.var.index].X.copy()

#Projecting into the ref PC space and identifying neighbors
cp.project_cells(sfdata, ref,
                 obs_columns=['leiden'],
                 fit_pca=True,
                 scale_data=True)

#Swapping for the original PCs, based on cell harmony-corrected data
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()
#nn-regressing the PCs
cp.nnregress(sfdata,
             ref,
             regress=['pca'],
             weighted=True)

#fitting into the original umap
cp.project_cells(sfdata, ref,
                 obs_columns=['leiden'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)

# %%
sfdata.write(base_procdata + 'sfdata_hoxb5projection.h5ad', compression='lzf')

# %% [markdown]
# ### Figures

# %% tags=[]
sc.pl.umap(ref, color='leiden')
sc.pl.umap(sfdata[~sfdata.obs.celltype.isin(('unassigned', 'CMP')),:],
           color=['celltype', 'ref_leiden', 'Procr', 'Dntt', 'Prtn3', 'Pf4', 'Klf1'],
           size = 30,
           wspace = 0.5,
          save = '_sfdata_projection.pdf')

# %% tags=[]
ref.obs['celltype'] = 'unassigned'
comb = ref.concatenate(sfdata, batch_key='batch', batch_categories=['ref', 'sfdata'], index_unique=None)
utils.plots.umap_subgroups(comb, key = 'celltype', toplot = sfdata.obs.celltype.unique(), file=base_figures + 'sfdata_projection_subgroups.pdf')

# %% tags=[]
utils.plots.umap_subgroups(comb, key = 'celltype_e', toplot = sfdata.obs.celltype_e.unique(), file=base_figures + 'sfdata_projection_subgroups_ESLAM.pdf')

# %% [markdown] tags=[]
# ## Bowling et al. 2020 data projection

# %%
bhsc = sc.read('./data/ext_data/Bowling2020/HSC_labelled_counts.h5ad')
utils.proc.lognorm(bhsc)

# %%
ref = hoxb5.copy()
#Unifying genes
ref = ref[:, ref.var.index[ref.var.index.isin(bhsc.var.index)]].copy()
ref.X = ref.raw[:,ref.var.index].X.copy()
del ref.raw
bhsc = bhsc[:, ref.var.index].copy()

ref.obs['HSC_labels'] = 'unassigned'
comb = ref.concatenate(bhsc, batch_key='batch', batch_categories=['ref', 'bhsc'], index_unique=None)

# %%
#First finding a common space for both batches
sc.pp.scale(comb)
sc.pp.pca(comb, n_comps=50)
sc.external.pp.harmony_integrate(comb, key='batch')
#Correcting for batch effects
sc.pp.neighbors(comb, n_neighbors=15, use_rep='X_pca_harmony')

#How well does the integation look like?
umap2 = cp.quick_umap(comb, use_rep='X_pca_harmony')
sc.pl.umap(comb, color=['batch', 'Procr'])
utils.plots.umap_subgroups(comb, key = 'HSC_labels',
                           toplot = bhsc.obs.HSC_labels.unique())

# %%
bhsc.obsm['X_pca'] = comb[bhsc.obs.index,:].obsm['X_pca_harmony'].copy()
ref.obsm['X_pca'] = comb[ref.obs.index,:].obsm['X_pca_harmony'].copy()

# %%
cp.project_cells(bhsc, ref,
                 obs_columns=['leiden'],
                 fit_pca=False,
                 scale_data=False,
                k=8)

#Looking at cross nearest neighbors
parent_nns = bhsc.uns['cross_nn'][bhsc.obs.HSC_labels == 'Parent',:].sum(axis=0).A1
ref.obs['parent_nns'] = parent_nns
childless_nns = bhsc.uns['cross_nn'][bhsc.obs.HSC_labels == 'Childless',:].sum(axis=0).A1
ref.obs['childless_nns'] = childless_nns
sc.pl.umap(ref, color=['parent_nns', 'childless_nns'])

#Swapping for the original PCs, based on cell harmony-corrected data
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()
#nn-regressing the PCs
cp.nnregress(bhsc,
             ref,
             regress=['pca'],
             weighted=True)

#fitting into the original umap
cp.project_cells(bhsc, ref,
                 obs_columns=['leiden'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)

#Looking at cross nearest neighbors

parent_nns = bhsc.uns['cross_nn'][bhsc.obs.HSC_labels == 'Parent',:].sum(axis=0).A1
ref.obs['parent_nns'] = parent_nns
childless_nns = bhsc.uns['cross_nn'][bhsc.obs.HSC_labels == 'Childless',:].sum(axis=0).A1
ref.obs['childless_nns'] = childless_nns
sc.pl.umap(ref, color=['parent_nns', 'childless_nns'])

# %%
bhsc.write(base_procdata + 'Bowling2020_hoxb5projection.h5ad', compression='lzf')

# %% [markdown]
# ### Figures

# %%
sc.pl.umap(ref, color='leiden')
sc.pl.umap(bhsc,
           color=['HSC_labels'],
           size = 30,
          save = '_Bowling2020_projection.pdf')

# %%
comb = ref.concatenate(bhsc, batch_key='batch', batch_categories=['ref', 'bhsc'], index_unique=None)
utils.plots.umap_subgroups(comb, key = 'HSC_labels',
                           toplot = bhsc.obs.HSC_labels.unique(), file=base_figures + 'Bowling2020_projection_subgroups.pdf')

# %% [markdown] tags=[]
# ## Projecting Weinreb et al 2020 data

# %%
d2clones = sc.read('./data/ext_data/Weinreb2020/Weinreb2020_d2withfates_logn.h5ad')
d2clones.raw = d2clones

# %%
x = d2clones.raw.to_adata()
x = utils.proc.delognorm(x)
mts = list(filter(lambda x: re.search('^mt-', x), x.var.index))

mt = x[:,mts].copy()
x.obs['mt_count'] = mt.X.todense().sum(axis = 1).A1
x.obs['mt_frac'] = x.obs['mt_count'] / x.obs['n_counts']
d2clones.obs['mt_frac'] = x.obs['mt_frac'].copy()

# sc.pl.umap(d2clones, color=['Annotation', 'n_counts', 'mt_frac'])
#There is a subset of cells with very low counts and higher mitochondrial fraction
#It has very high expression of mito genes but also ribosomal genes, Gapdh etc.
#Not sure what these cells are but some of them project aberrantly, so we remove them

# %%
d2clones = d2clones[d2clones.obs.n_counts >= 1800,:].copy()

# %%
ref = hoxb5.copy()
#Unifying genes
ref = ref[:, ref.var.index[ref.var.index.isin(d2clones.var.index)]].copy()
ref.X = ref.raw[:,ref.var.index].X.copy()
del ref.raw
d2clones = d2clones[:, ref.var.index].copy()
#Combining data
comb = ref.concatenate(d2clones, batch_key='batch', batch_categories=['ref', 'd2clones'], index_unique=None)
comb.var['highly_variable'] = comb.var['highly_variable-ref']

#The harmony integration does not work at all, scanorama does not work either, trying Seurat

# %%
#Running the batch correction with Seurat
comb_cor = cp.run_SeuratCCA(comb, batch_key='batch', reference='ref', k_filter=10)

# %%
#We will need PCs based on the common space
sc.pp.scale(ref)
sc.pp.pca(ref, n_comps=50)
sc.pp.neighbors(ref, n_neighbors=15)

#Setting the Seurat-corrected values for the target data (ref data is already scaled)
d2clones.X = comb_cor[d2clones.obs.index, d2clones.var.index].X.copy()

#Projecting into the ref PC space and identifying neighbors
cp.project_cells(d2clones, ref,
                 obs_columns=['leiden'],
                 fit_pca=True,
                 scale_data=True)

#Swapping for the original PCs, based on cell harmony-corrected data
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()
#nn-regressing the PCs
cp.nnregress(d2clones,
             ref,
             regress=['pca'],
             weighted=True)

#fitting into the original umap
cp.project_cells(d2clones, ref,
                 obs_columns=['leiden'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)

# %%
d2clones.write(base_procdata + 'Weirenb2020_hoxb5projection.h5ad', compression='lzf')

# %%
fates = d2clones.obs.fateclass.value_counts()
tokeep = fates[fates > 7].index
tokeep = tokeep[tokeep != 'nofate']
d2clones = d2clones[d2clones.obs.fateclass.isin(tokeep),:].copy()

# %% [markdown]
# ### Figures

# %%
sc.pl.umap(ref, color='leiden')
sc.pl.umap(d2clones,
           color=['fateclass', 'Procr', 'Dntt', 'Pf4', 'Prtn3', 'Klf1'],
           size = 30,
          save = '_Weinreb2020_projection.pdf')

# %% tags=[]
comb = ref.concatenate(d2clones, batch_key='batch', batch_categories=['ref', 'bhsc'], index_unique=None)
utils.plots.umap_subgroups(comb, key = 'fateclass',
                           toplot = d2clones.obs.fateclass.unique(), file=base_figures + 'Weinreb2020_projection_subgroups.pdf')

# %% [markdown]
# ## Projecting Dong et al 2020 data

# %%
#Loading Nestorowa et al. 2016 data
txdata = sc.read('./data/ext_data/Dong2020/dong2020_BM_Tx_proc_logn.h5ad')
txdata.raw = txdata

# %%
txdata.obs.Time = txdata.obs.Time.astype(str)
txdata.obs.loc[txdata.obs.Time == 'nan', 'Time'] = 'source'
sc.pl.umap(txdata, color = ['Time'])

# %% [markdown]
# ### Projecting data

# %%
ref = hoxb5.copy()
#Unifying genes
ref = ref[:, ref.var.index[ref.var.index.isin(txdata.var.index)]].copy()
ref.X = ref.raw[:,ref.var.index].X.copy()
del ref.raw
txdata = txdata[:, ref.var.index].copy()
#Combining data
comb = ref.concatenate(txdata, batch_key='batch', batch_categories=['ref', 'txdata'], index_unique=None)
comb.var['highly_variable'] = comb.var['highly_variable-ref']

#Running the batch correction with Seurat
comb_cor = cp.run_SeuratCCA(comb, batch_key='batch', reference='ref')

# %%
#We will need PCs based on the common space
sc.pp.scale(ref)
sc.pp.pca(ref, n_comps=50)
sc.pp.neighbors(ref, n_neighbors=15)

#Setting the Seurat-corrected values for the target data (ref data is already scaled)
txdata.X = comb_cor[txdata.obs.index, txdata.var.index].X.copy()

#Projecting into the ref PC space and identifying neighbors
cp.project_cells(txdata, ref,
                 obs_columns=['leiden'],
                 fit_pca=True,
                 scale_data=True)

#Swapping for the original PCs, based on cell harmony-corrected data
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()
#nn-regressing the PCs
cp.nnregress(txdata,
             ref,
             regress=['pca'],
             weighted=True)

#fitting into the original umap
cp.project_cells(txdata, ref,
                 obs_columns=['leiden'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)

# %%
# umapref2 = cp.quick_umap(hoxb5, use_rep='X_pca_harmony')
# #fitting into the original umap
# cp.project_cells(sfdata, ref,
#                  obs_columns=['leiden'],
#                  fit_pca=False,
#                  scale_data=False,
#                  umap_ref=umapref2)

# %% [markdown]
# ### Figures

# %%
txdata.write(base_procdata + 'Dong2020_hoxb5projection.h5ad', compression='lzf')

# %%
txdata.obs['ref_leiden'] = pd.Categorical(txdata.obs.ref_leiden, categories=hoxb5.obs.leiden.cat.categories)
tx_cellnos = pd.crosstab(txdata.obs.ref_leiden, [txdata.obs.Time], dropna=False)
tx_cellnos.to_csv(base_procdata + 'tx_cellnos.csv')
tx_cellnos

# %%
sc.pl.umap(ref, color='leiden')
sc.pl.umap(txdata,
           color=['anno_quick', 'ref_leiden', 'Procr', 'Dntt', 'Prtn3', 'Pf4', 'Klf1'],
           size = 30,
           wspace = 0.5,
          save = '_txdata_projection.pdf')

# %% tags=[]
ref.obs['anno_quick'] = 'unassigned'
comb = ref.concatenate(txdata, batch_key='batch', batch_categories=['ref', 'txdata'], index_unique=None)
utils.plots.umap_subgroups(comb, key = 'anno_quick', toplot = txdata.obs.anno_quick.unique(), file=base_figures + 'txdata_projection_subgroups.pdf')

# %% tags=[]
utils.plots.umap_subgroups(comb, key = 'Time', toplot = txdata.obs.Time.unique(), file=base_figures + 'txdata_projection_subgroups_Time.pdf')

# %%
tx_cellnos_permouse = pd.crosstab(txdata.obs.Mouse_ID, [txdata.obs.Time, txdata.obs.ref_leiden])
tx_cellnos_permouse.to_csv(base_procdata + 'tx_cellnos_permouse.csv')
tx_cellnos_permouse

# %% [markdown]
# ## Pei2020 data integration

# %%
plx = sc.read('./data/ext_data/Pei2020/Pei2020_small_HSCMPPclones.h5ad')
# sc.pp.subsample(plx, n_obs=2000)

# %%
a

# %%
ref = hoxb5.copy()
# sc.pp.subsample(ref, n_obs=2000)
#Unifying genes
ref = ref[:, ref.var.index[ref.var.index.isin(plx.var.index)]].copy()
ref.X = ref.raw[:,ref.var.index].X.copy()
del ref.raw
plx = plx[:, ref.var.index].copy()
#Combining data
comb = ref.concatenate(plx, batch_key='batch', batch_categories=['ref', 'plx'], index_unique=None)
comb.var['highly_variable'] = comb.var['highly_variable-ref']

#Running the batch correction with Seurat
comb_cor = cp.run_SeuratCCA(comb, batch_key='batch', reference='ref')

# %%
#We will need PCs based on the common space
sc.pp.scale(ref)
sc.pp.pca(ref, n_comps=50)
sc.pp.neighbors(ref, n_neighbors=15)

#Setting the Seurat-corrected values for the target data (ref data is already scaled)
plx.X = comb_cor[plx.obs.index, plx.var.index].X.copy()

#Projecting into the ref PC space and identifying neighbors
cp.project_cells(plx, ref,
                 obs_columns=['leiden'],
                 fit_pca=True,
                 scale_data=True)

#Swapping for the original PCs, based on cell harmony-corrected data
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()
#nn-regressing the PCs
cp.nnregress(plx,
             ref,
             regress=['pca'],
             weighted=True)

#fitting into the original umap
cp.project_cells(plx, ref,
                 obs_columns=['leiden'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)

# %%
sc.pl.umap(ref, color='leiden')
sc.pl.umap(plx,
           color=['population'],
           size = 30,
          save = '_Pei2020_projection.pdf')

# %%
plx = plx[~plx.obs.population.isin(['preT', 'proB']),:]
plx.write(base_procdata + 'Pei2020_hoxb5projection.h5ad', compression='lzf')

# %%
sc.pl.umap(plx, color=['Dntt', 'Il7r'])
sc.pl.umap(hoxb5, color=['Dntt', 'Il7r'])

# %%
comb = ref.concatenate(plx, batch_key='batch', batch_categories=['ref', 'plx'], index_unique=None)
utils.plots.umap_subgroups(comb, key = 'population',
                           toplot = plx.obs.population.unique(), file=base_figures + 'Pei2020_projection_population_subgroups.pdf')
utils.plots.umap_subgroups(comb, key = 'fate',
                           toplot = plx.obs.fate.unique(), file=base_figures + 'Pei2020_projection_fate_subgroups.pdf')

# %%

# %%

# %%

# %%

# %%
ref = hoxb5.copy()
# sc.pp.subsample(ref, n_obs=2000)
#Unifying genes
ref = ref[:, ref.var.index[ref.var.index.isin(plx.var.index)]].copy()
ref.X = ref.raw[:,ref.var.index].X.copy()
del ref.raw
plx = plx[:, ref.var.index].copy()

ref.obs['HSC_labels'] = 'unassigned'
comb = ref.concatenate(plx, batch_key='batch', batch_categories=['ref', 'plx'], index_unique=None)

# %%
#First finding a common space for both batches
sc.pp.scale(comb)
sc.pp.pca(comb, n_comps=50)
sc.external.pp.harmony_integrate(comb, key='batch', sigma=0.2) #try setting sigma >0.1 if not try seurat correction
#Correcting for batch effects
sc.pp.neighbors(comb, n_neighbors=15, use_rep='X_pca_harmony')

#How well does the integation look like?
umap2 = cp.quick_umap(comb, use_rep='X_pca_harmony')
sc.pl.umap(comb, color=['batch', 'Procr'])
utils.plots.umap_subgroups(comb, key = 'population',
                           toplot = plx.obs.population.unique())

# %%
plx.obsm['X_pca'] = comb[plx.obs.index,:].obsm['X_pca_harmony'].copy()
ref.obsm['X_pca'] = comb[ref.obs.index,:].obsm['X_pca_harmony'].copy()

# %%
cp.project_cells(plx, ref,
                 obs_columns=['leiden'],
                 fit_pca=False,
                 scale_data=False,
                k=5)

#Looking at cross nearest neighbors
# parent_nns = plx.uns['cross_nn'][plx.obs.HSC_labels == 'Parent',:].sum(axis=0).A1
# ref.obs['parent_nns'] = parent_nns
# childless_nns = plx.uns['cross_nn'][plx.obs.HSC_labels == 'Childless',:].sum(axis=0).A1
# ref.obs['childless_nns'] = childless_nns
# sc.pl.umap(ref, color=['parent_nns', 'childless_nns'])

#Swapping for the original PCs, based on cell harmony-corrected data
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()
#nn-regressing the PCs
cp.nnregress(plx,
             ref,
             regress=['pca'],
             weighted=True)

#fitting into the original umap
cp.project_cells(plx, ref,
                 obs_columns=['leiden'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)

#Looking at cross nearest neighbors

# parent_nns = plx.uns['cross_nn'][plx.obs.HSC_labels == 'Parent',:].sum(axis=0).A1
# ref.obs['parent_nns'] = parent_nns
# childless_nns = plx.uns['cross_nn'][plx.obs.HSC_labels == 'Childless',:].sum(axis=0).A1
# ref.obs['childless_nns'] = childless_nns
# sc.pl.umap(ref, color=['parent_nns', 'childless_nns'])

# %%
plx = plx[~plx.obs.population.isin(['preT', 'proB']),:]
plx.write(base_procdata + 'Pei2020_hoxb5projection.h5ad', compression='lzf')

# %% [markdown]
# ### Figures

# %%
sc.pl.umap(ref, color='leiden')
sc.pl.umap(plx,
           color=['population'],
           size = 30,
          save = '_Pei2020_projection.pdf')

# %%
comb = ref.concatenate(plx, batch_key='batch', batch_categories=['ref', 'plx'], index_unique=None)
utils.plots.umap_subgroups(comb, key = 'population',
                           toplot = plx.obs.population.unique(), file=base_figures + 'Pei2020_projection_population_subgroups.pdf')
utils.plots.umap_subgroups(comb, key = 'fate',
                           toplot = plx.obs.fate.unique(), file=base_figures + 'Pei2020_projection_fate_subgroups.pdf')

# %% [markdown]
# With default harmony
# when using a random subsample of 2000 it works fine
# but with full dataset the CLP cells get forced into the main ladnscape

# %%

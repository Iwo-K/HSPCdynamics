# -*- coding: utf-8 -*-
# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.10.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# %%
#Singularity container used:
# !echo $SINGULARITY_CONTAINER
# !hostname

# %% [markdown]
# # Basic analysis of label distribution per cluster

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
import pickle
from scipy.sparse import csr_matrix
import utils.proc
from utils.load10x import load_10xfiles
import utils.proc, utils.plots

import cellproject as cp

sc.settings.verbosity = 3
sc.settings.figdir = './figures/05script/'
base_figures = './figures/05script/'
base_procdata = './procdata/05script/'
for i in [sc.settings.figdir, base_figures, base_procdata]:
    pathlib.Path(i).mkdir(parents=True, exist_ok=True)

# %%
#Figure settings
sc.set_figure_params(color_map = 'viridis', dpi_save = 350)
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
comb = sc.read('procdata/04script/combined_filt.h5ad')
meta = pd.read_csv('procdata/04script/combined_sample.meta.csv', index_col = 0)
meta = meta.sort_values(by ='timepoint_tx_days', ascending = False)
comb.obs['longname'] = [meta.loc[i, 'longname'] for i in comb.obs.biosample_id]

# %% [markdown] tags=[]
# ## PCA landscape

# %%
sc.pl.umap(comb, color = 'leiden', legend_loc = 'on data')

# %% [markdown]
# ## Plotting cells per timepoint

# %%
comb.obs['tom_timepoint'] = comb.obs.tom.astype(str) + '_' + comb.obs.timepoint_tx_days.astype(str)

# %%
utils.plots.umap_subgroups(comb, key = 'tom_timepoint', toplot = ['pos_7', 'pos_269'],file=base_figures + './cells_timepoints.pdf')

# %% [markdown]
# ## Abundance per cluster

# %%
# Counting cells in each cluster seprately for Tom+ and Tom- cells
tompos = comb[comb.obs.tom == 'pos',:].copy()
tomneg = comb[comb.obs.tom == 'neg',:].copy()

tompos_clcounts = pd.crosstab(tompos.obs.leiden, tompos.obs.biosample_id)
tompos_clcounts.to_csv(base_procdata + 'tompos_cluster_cellcounts.csv')
# Calculating fraction of cells in each clusters
tompos_clcountsN = tompos_clcounts / tompos_clcounts.sum(axis = 0)
tompos_clcountsN.to_csv(base_procdata + 'tompos_cluster_cellfractions.csv')

tomneg_clcounts = pd.crosstab(tomneg.obs.leiden, tomneg.obs.biosample_id)
tomneg_clcounts.to_csv(base_procdata + 'tomneg_cluster_cellcounts.csv')
# Calculating fraction of cells in each clusters
tomneg_clcountsN = tomneg_clcounts / tomneg_clcounts.sum(axis = 0)
tomneg_clcountsN.to_csv(base_procdata + 'tomneg_cluster_cellfractions.csv')


# %%
# Calculating the ratios between Tom+ cluster fractions and Tom- cluster fractions, they need to match the start_age and tx
meta['start_age_tx_days'] = meta.start_age + '_' + meta.timepoint_tx_days.astype('str')
tompos_clcountsN_rel = tompos_clcountsN.copy()

for i in tompos_clcountsN.columns:
    #Finding the matching references and calculating mean
    group = meta.loc[meta.biosample_id == i,'start_age_tx_days'][0]  
    r = meta.loc[(meta.start_age_tx_days == group) & (meta.tom == 'neg'),'biosample_id']
    r = tomneg_clcountsN.loc[:,r]
    r = r.mean(axis = 1)
    #Calculating the relative number of cells
    tompos_clcountsN_rel.loc[:,i] = np.log2((tompos_clcountsN_rel.loc[:,i] / r) + 0.01)
    
tompos_clcountsN_rel.to_csv(base_procdata + 'tompos_cluster_cellfractions_relative.csv')

# %% tags=[]
meta = meta.sort_values(by = 'timepoint_tx_days', ascending = False)
toplot = []

for i in meta.biosample_id[meta.tom == 'pos']:
    fc = tompos_clcountsN_rel.loc[:,i]
    name = meta.loc[i, 'longname']
    tomneg.obs[name + '_fc'] = [fc[i] for i in tomneg.obs.leiden]
    toplot.append(name + '_fc')

from matplotlib.colors import LinearSegmentedColormap
cmap2 = utils.plots.cmap_RdBu2(values = None, vmin = -5, vmax = 5)
sc.pl.umap(tomneg, color = toplot, cmap = cmap2, vmin = -5, vmax = 5, save='_clusterfractions_all.pdf')

# %%
avs = tompos_clcountsN_rel.melt(ignore_index = False, var_name = 'biosample_id')
avs['cluster'] = avs.index
avs['timepoint_tx_days'] = meta.loc[avs.biosample_id, 'timepoint_tx_days'].values

avsum = avs.groupby(['timepoint_tx_days', 'cluster']).mean()
print(avsum)

toplot = []
for i in avsum.index.get_level_values(0).unique():
    ratio = avsum.loc[i,:].copy()
    ratio = ratio.value

    #Clipping at 5
    clip = 5
    ratio.iloc[ratio > clip] = clip
    ratio.iloc[ratio < -clip] = -clip

    tomneg.obs[str(i) + 'd_avratio'] = [ratio[i] for i in tomneg.obs.leiden]
    toplot.append(str(i) + 'd_avratio')

cmap2 = utils.plots.cmap_RdBu2(values = None, vmin = -clip, vmax = clip)
sc.pl.umap(tomneg, color = toplot, cmap = cmap2, vmin = -clip, vmax = clip, save = 'rel_abundance_percluster_average.pdf')


# %% [markdown]
# ## Abundance per cluster - adjusted with flow data
# To have a better estimate of Tom+ labelling frequency we integrate the scRNA-Seq data with observed total Tom+ cell number from flow cytometry data

# %%
flow = pd.read_csv('./data/flow/flow_summarise_data.csv')
flow = flow.loc[flow.data_type.isin(('10x', 'SS2')),:]
flow.index = flow.biosample_id.values
# flow = flow[['biosample_id', 'central_singlets_live_linneg_Lgate_Tompos__estim_total']]

# %% [markdown]
# ### Tom+ cells

# %% tags=[]
#sc_ncells_filt is the total number of cells in the filtered clusters, sc_ncells_total is the total number of cells in all clusters (ie cells after QC)
tompos_cl_flow = tompos_clcounts.T
tompos_cl_flow.columns = tompos_cl_flow.columns.astype(str)
tompos_cl_flow['sc_ncells_filt'] = tompos_cl_flow.sum(axis=1)

cols = ['central_singlets_live_linneg_Lgate_Tompos__estim_total', 'start_age', 'days_postTx']
tompos_cl_flow = tompos_cl_flow.join(flow.loc[tompos_cl_flow.index, cols])
tompos_cl_flow = tompos_cl_flow.join(meta.loc[tompos_cl_flow.index, 'ncells'])
tompos_cl_flow = tompos_cl_flow.rename(columns={'central_singlets_live_linneg_Lgate_Tompos__estim_total': 'flow_total',
                                               'days_postTx' : 'time',
                                               'ncells' : 'sc_ncells_total'})
#Saving separately numbers for aged animals
tompos_cl_flow_aged = tompos_cl_flow.loc[tompos_cl_flow.start_age == 'aged',:].copy()
tompos_cl_flow_aged.to_csv(base_procdata + 'model_input_tompos_aged.csv')

tompos_cl_flow = tompos_cl_flow.loc[tompos_cl_flow.start_age == 'young',:].copy()
tompos_cl_flow.to_csv(base_procdata + 'model_input_tompos.csv')
tompos_cl_flow

# %%
m = np.array([-11 for i in range(len(np.unique(tompos_cl_flow.time)))],dtype = 'float') 
s = np.array([-11 for i in range(len(np.unique(tompos_cl_flow.time)))],dtype = 'float') 
k = -1
for j in np.unique(tompos_cl_flow.time):
    k+=1
    m[k] = np.mean(tompos_cl_flow.loc[tompos_cl_flow.time == j].flow_total)
    s[k] = np.std(tompos_cl_flow.loc[tompos_cl_flow.time == j].flow_total)/np.sqrt(np.sum(tompos_cl_flow.time == j))
    
plt.errorbar(range(len(np.unique(tompos_cl_flow.time))),m,s)

dfsum = pd.DataFrame({'time' : np.unique(tompos_cl_flow.time), 'tot_mean' : m, 'tot_sd' : s})
dfsum.to_csv(base_procdata + 'mean_std_percluster_tompos.csv')

# %% [markdown]
# ### Tom- cells

# %%
# The flow data and metadata is missing biosample_ids fo the Tom- samples, need to assign them.

# %%
#Subsetting for Tom- samples and finding mapping mouse_id to biosample_id
mapid = meta.loc[meta.tom == 'neg', ['mouse_id', 'biosample_id']].copy()
if mapid.mouse_id.duplicated().any():
    raise Exception('There are some duplicated records!')
mapid.index = mapid.mouse_id.astype(str).values

tomneg_flow = flow.copy()
tomneg_flow.index = tomneg_flow.mouse_id.values
if tomneg_flow.mouse_id.duplicated().any():
    raise Exception('There are some duplicated records!') 
    
#Subsetting for Tom- samples and adding biosample_id for Tom- samples
tomneg_flow= tomneg_flow.loc[tomneg_flow.mouse_id.isin(mapid.mouse_id),:]
tomneg_flow['biosample_id_tomneg'] = mapid.loc[tomneg_flow.mouse_id, 'biosample_id']
tomneg_flow.index = tomneg_flow['biosample_id_tomneg'].values

# %%
#Estimating total number of Tom- cells
posno = tomneg_flow.central_singlets_live_linneg_Lgate_Tompos__estim_total
posfreq = tomneg_flow.central_singlets_live_linneg_Lgate_Tompos__freq_of_parent
negno = (posno/posfreq) - posno
tomneg_flow['central_singlets_live_linneg_Lgate_Tomneg__estim_total'] = negno

# %%
tomneg_cl_flow = tomneg_clcounts.T

tomneg_cl_flow = tomneg_clcounts.T
tomneg_cl_flow.columns = tomneg_cl_flow.columns.astype(str)
tomneg_cl_flow['sc_ncells_filt'] = tomneg_cl_flow.sum(axis=1)


cols = ['central_singlets_live_linneg_Lgate_Tomneg__estim_total', 'start_age', 'days_postTx']
tomneg_cl_flow = tomneg_cl_flow.join(tomneg_flow.loc[tomneg_cl_flow.index, cols])
tomneg_cl_flow = tomneg_cl_flow.join(meta.loc[tomneg_cl_flow.index, 'ncells'])
tomneg_cl_flow = tomneg_cl_flow.rename(columns={'central_singlets_live_linneg_Lgate_Tomneg__estim_total': 'flow_total',
                                               'days_postTx' : 'time',
                                               'ncells' : 'sc_ncells_total'})

tomneg_cl_flow_aged = tomneg_cl_flow.loc[tomneg_cl_flow.start_age == 'aged',:].copy()
tomneg_cl_flow_aged.to_csv(base_procdata + 'model_input_tomneg_aged.csv')

tomneg_cl_flow = tomneg_cl_flow.loc[tomneg_cl_flow.start_age == 'young',:].copy()
tomneg_cl_flow.to_csv(base_procdata + 'model_input_tomneg.csv')
tomneg_cl_flow

#There is a problem because in the 10x metadata there is an accidental space in one of the ids!

# %% [markdown] jp-MarkdownHeadingCollapsed=true tags=[]
# ## Cluster0 substructure

# %% tags=[]
cl0 = comb[comb.obs.leiden == '0', :].copy()
sc.pp.neighbors(cl0)
sc.tl.leiden(cl0, resolution=0.4)

# %% tags=[]
sc.pl.umap(cl0, color=['leiden', 'data_type', 'start_age', 'dpt_pseudotime', 'HSCscore', 'phase', 'S_score',
                       'Procr', 'Satb1', 'Dntt', 'Elane', 'Gata1', 'Pf4', 'Vwf', 'Apoe'], size=25)

# %%
sc.tl.umap(cl0)
sc.pl.umap(cl0, color=['leiden', 'data_type', 'start_age', 'dpt_pseudotime', 'HSCscore', 'phase', 'S_score',
                       'Procr', 'Satb1', 'Dntt', 'Elane', 'Gata1', 'Pf4', 'Vwf', 'Apoe'], size=25)

# %% tags=[]
cl0.obs.timepoint_tx_days = cl0.obs.timepoint_tx_days.astype(str)
utils.plots.umap_subgroups(cl0[(cl0.obs.tom == 'pos') & (cl0.obs.start_age == 'young') , :],
                           key='timepoint_tx_days', toplot=['3', '7', '12',
                                                            '27', '49', '76', 
                                                           '112', '161', '269'])

# %%
tompos = cl0[cl0.obs.tom == 'pos',:].copy()
tomneg = cl0[cl0.obs.tom == 'neg',:].copy()

tompos_clcounts = pd.crosstab(tompos.obs.leiden, tompos.obs.biosample_id)
tompos_clcountsN = tompos_clcounts / tompos_clcounts.sum(axis = 0)
tomneg_clcounts = pd.crosstab(tomneg.obs.leiden, tomneg.obs.biosample_id)
tomneg_clcountsN = tomneg_clcounts / tomneg_clcounts.sum(axis = 0)

# %%
#Normalising data with matching references, they need to match the start_age and tx
meta['start_age_tx_days'] = meta.start_age + '_' +meta.timepoint_tx_days.astype('str')
tompos_clcountsN_rel = tompos_clcountsN.copy()

for i in tompos_clcountsN.columns:
    #Finding the matching references and calculating mean
    group = meta.loc[meta.biosample_id == i,'start_age_tx_days'][0]  
    r = meta.loc[(meta.start_age_tx_days == group) & (meta.tom == 'neg'),'biosample_id']
    r = tomneg_clcountsN.loc[:,r]
    r = r.mean(axis = 1)
    #Calculating the relative number of cells
    tompos_clcountsN_rel.loc[:,i] = np.log2((tompos_clcountsN_rel.loc[:,i] / r) + 0.01)

# %% tags=[]
meta = meta.sort_values(by = 'timepoint_tx_days', ascending = False)
toplot = []

for i in meta.biosample_id[meta.tom == 'pos']:
    fc = tompos_clcountsN_rel.loc[:,i]
    name = meta.loc[i, 'longname']
    tomneg.obs[name + '_fc'] = [fc[i] for i in tomneg.obs.leiden]
    toplot.append(name + '_fc')

from matplotlib.colors import LinearSegmentedColormap
cmap2 = utils.plots.cmap_RdBu2(values = None, vmin = -5, vmax = 5)
sc.pl.umap(tomneg, color = toplot, cmap = cmap2, vmin = -5, vmax = 5)

# %% tags=[]
avs = tompos_clcountsN_rel.melt(ignore_index = False, var_name = 'biosample_id')
avs['cluster'] = avs.index
avs['timepoint_tx_days'] = meta.loc[avs.biosample_id, 'timepoint_tx_days'].values

avsum = avs.groupby(['timepoint_tx_days', 'cluster']).mean()
print(avsum)

toplot = []
for i in avsum.index.get_level_values(0).unique():
    ratio = avsum.loc[i,:].copy()
    ratio = ratio.value

    #Clipping at 5
    clip = 5
    ratio.iloc[ratio > clip] = clip
    ratio.iloc[ratio < -clip] = -clip

    tomneg.obs[str(i) + 'd_avratio'] = [ratio[i] for i in tomneg.obs.leiden]
    toplot.append(str(i) + 'd_avratio')

cmap2 = utils.plots.cmap_RdBu2(values = None, vmin = -clip, vmax = clip)
sc.pl.umap(tomneg, color = toplot, cmap = cmap2, vmin = -clip, vmax = clip)


# %%
# df = pd.DataFrame(index=tompos_clcounts.index)
# for i in meta.timepoint_tx_days.unique():
#     print(i)
#     negs = meta.loc[(meta.timepoint_tx_days == i) & (meta.tom == 'neg') & (meta.start_age == 'young'),'biosample_id'].values
#     negs = tomneg_clcounts[negs].mean(axis=1)
#     poss = meta.loc[(meta.timepoint_tx_days == i) & (meta.tom == 'pos') & (meta.start_age == 'young'),'biosample_id'].values
#     poss = (tompos_clcounts[poss].T / negs).T
#     df = pd.concat((df, poss), axis=1)

# %%

# %%

# %%

# %%

# %% [markdown] tags=[] jp-MarkdownHeadingCollapsed=true
# ## DM landscape

# %%
comb.obsm['X_umap'] = comb.obsm['X_umap_DM_2d'].copy()
sc.pl.umap(comb, color = 'leiden_DM', legend_loc = 'on data')

# %%
tompos = comb[comb.obs.tom == 'pos',:].copy()
tomneg = comb[comb.obs.tom == 'neg',:].copy()
tompos = tompos[tompos.obs.biosample_id !='SLX20264_SITTB7_146814.46e',:].copy()
meta = meta.loc[meta.biosample_id != 'SLX20264_SITTB7_146814.46e',:].copy()

# %%
tompos_clcounts = pd.crosstab(tompos.obs.leiden_DM, tompos.obs.biosample_id)
tompos_clcounts.to_csv(base_procdata + 'DM_tompos_cluster_cellcounts.csv')
tompos_clcountsN = tompos_clcounts / tompos_clcounts.sum(axis = 0)
tompos_clcountsN.to_csv(base_procdata + 'DM_tompos_cluster_cellfractions.csv')


tomneg_clcounts = pd.crosstab(tomneg.obs.leiden_DM, tomneg.obs.biosample_id)
tomneg_clcounts.to_csv(base_procdata + 'DM_tomneg_cluster_cellcounts.csv')
tomneg_clcountsN = tomneg_clcounts / tomneg_clcounts.sum(axis = 0)
tomneg_clcountsN.to_csv(base_procdata + 'DM_tomneg_cluster_cellfractions.csv')


# %%
#Normalising data with matching references, they need to match the start_age and tx
meta['start_age_tx_days'] = meta.start_age + '_' +meta.timepoint_tx_days.astype('str')
tompos_clcountsN_rel = tompos_clcountsN.copy()

for i in tompos_clcountsN.columns:
    #Finding the matching references and calculating mean
    group = meta.loc[meta.biosample_id == i,'start_age_tx_days'][0]  
    r = meta.loc[(meta.start_age_tx_days == group) & (meta.tom == 'neg'),'biosample_id']
    r = tomneg_clcountsN.loc[:,r]
    r = r.mean(axis = 1)
    #Calculating the relative number of cells
    tompos_clcountsN_rel.loc[:,i] = np.log2((tompos_clcountsN_rel.loc[:,i] / r) + 0.01)
    
tompos_clcountsN_rel.to_csv(base_procdata + 'DM_tompos_cluster_cellfractions_relative.csv')

# %% tags=[]
meta = meta.sort_values(by = 'timepoint_tx_days', ascending = False)
toplot = []

for i in meta.biosample_id[meta.tom == 'pos']:
    fc = tompos_clcountsN_rel.loc[:,i]
    name = meta.loc[i, 'longname']
    tomneg.obs[name + '_fc'] = [fc[i] for i in tomneg.obs.leiden_DM]
    toplot.append(name + '_fc')

from matplotlib.colors import LinearSegmentedColormap
cmap2 = utils.plots.cmap_RdBu2(values = None, vmin = -5, vmax = 5)
sc.pl.umap(tomneg, color = toplot, cmap = cmap2, vmin = -5, vmax = 5, save='DM_clusterfractions_all.pdf')

# %% tags=[]
avs = tompos_clcountsN_rel.melt(ignore_index = False, var_name = 'biosample_id')
avs['cluster'] = avs.index
avs['timepoint_tx_days'] = meta.loc[avs.biosample_id, 'timepoint_tx_days'].values

avsum = avs.groupby(['timepoint_tx_days', 'cluster']).mean()
print(avsum)

for i in avsum.index.get_level_values(0).unique():
    ratio = avsum.loc[i,:].copy()
    ratio = ratio.value

    #Clipping at 3
    clip = 5
    ratio.iloc[ratio > clip] = clip
    ratio.iloc[ratio < -clip] = -clip

    tomneg.obs[str(i) + 'd_avratio'] = [ratio[i] for i in tomneg.obs.leiden_DM]
    cmap2 = utils.plots.cmap_RdBu2(values = None, vmin = -clip, vmax = clip)
    sc.pl.umap(tomneg, color = str(i) + 'd_avratio', cmap = cmap2, vmin = -clip, vmax = clip)

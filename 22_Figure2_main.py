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
# # Plots for Figure 2

# %% [markdown]
# ## Setup

# %%
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib as mpl
import pathlib
import anndata
import utils.proc, utils.plots

import warnings
warnings.filterwarnings("ignore", category=FutureWarning) 

sc.settings.verbosity = 3
base_figures = './figures/22script/'
sc.settings.figdir = base_figures
base_procdata = './procdata/22script/'
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

# %% [markdown]
#
# ## Plotting unfiltered data - basic information for supplement

# %%
comb_unfilt = sc.read('procdata/04script/combined.h5ad')
x = np.where(comb_unfilt.obs.anno_man.cat.categories == 'Int prog')[0][0]
comb_unfilt.uns['anno_man_colors'][x] = "#C5C5C5"
sc.pl.umap(comb_unfilt, color=['anno_man', 'leiden'], wspace=0.7, save='unfiltered_annotation_leiden.pdf')
sc.pl.umap(comb_unfilt, color=['leiden'], legend_loc='on data', save='unfiltered_clusters_ondata.pdf')

# %% [markdown]
# ## Plotting filtered data - basic information

# %%
hoxb5 = sc.read('procdata/04script/combined_filt.h5ad')
meta = pd.read_csv('procdata/04script/combined_sample.meta.csv', index_col = 0)
meta = meta.sort_values(by ='timepoint_tx_days', ascending = False)
hoxb5.obs['longname'] = [meta.loc[i, 'longname'] for i in hoxb5.obs.biosample_id]
hoxb5.obsm['X_umap'] = hoxb5.obsm['X_umap_2d']

# %%
from scanpy.plotting.palettes import default_20
hoxb5.uns['anno_man_colors'] = default_20[0 : len(hoxb5.obs.anno_man.cat.categories)]
x = np.where(hoxb5.obs.anno_man.cat.categories == 'Int prog')[0][0]
hoxb5.uns['anno_man_colors'][x] = "#C5C5C5"
sc.pl.umap(hoxb5, color=['anno_man', 'leiden', 'phase'], wspace=0.4, save='_annotation_leiden.pdf')

# %%
sc.pl.umap(hoxb5, color=['leiden'], legend_loc='on data', wspace=0.4, save='_clusters_ondata.pdf')
sc.pl.umap(hoxb5, color=['nn_HSCscore'], save='_HSCscore.pdf')

# %%
sc.pl.umap(hoxb5, color=['dpt_pseudotime'], cmap='rainbow', save='_dpt.pdf')

# %%
utils.plots.umap_subgroups(hoxb5, key = 'data_type',
                           toplot = ['SS2'], 
                           file=base_figures + 'umap_datatype.pdf')

# %%
from matplotlib.colors import LinearSegmentedColormap
cmap2 = LinearSegmentedColormap.from_list('mycmap', [(0, '#e8e8e8'),
                                                     (1, '#cc0404')])

# %%
markers = ['Procr', 'Klf1', 'Pf4', 'Csf1r', 'Siglech', 'Dntt', 'Vpreb3', 'Elane', 'H2-Aa', 'Prg2', 'Mcpt8', 'Cma1'] #Itgax
sc.pl.umap(hoxb5, color=markers, cmap=cmap2, save='_markers.pdf', size=5)

# %%
markers = ['Hba-a1', 'Pf4', 'Satb1', 'Elane']
sc.pl.umap(hoxb5, color=markers, cmap=cmap2)

# %%
markers = ['Cish', 'Socs2']
sc.pl.umap(hoxb5, color=markers, cmap=cmap2, size=5, save='_Cish_Socs2.pdf')

# %% [markdown]
# ## Cells per timepoint

# %%
hoxb5.obs['tom_timepoint'] = hoxb5.obs.tom.astype(str) + '_' + hoxb5.obs.timepoint_tx_days.astype(str)
x = hoxb5.obs.tom_timepoint.value_counts()

toplot = ['pos_3', 'pos_12', 'pos_27', 'pos_269']
utils.plots.umap_subgroups(hoxb5, key = 'tom_timepoint',
                           toplot = toplot, 
                           subsample_to={i : 8000 if x[i] > 8000 else x[i].max() for i in x.index},
                           file=base_figures + 'cells_timepoints.pdf')

toplot = ['neg_269']
utils.plots.umap_subgroups(hoxb5, key = 'tom_timepoint',
                           toplot = toplot, 
                           subsample_to={'neg_269' : 8000},
                           file=base_figures + 'cells_neg269.pdf')

# %% [markdown]
# ## Plotting relative fractions of cells per timepoint

# %%
# Loading number of cells per cluster
tompos_clcounts = pd.read_csv('procdata/05script/model_input_tompos.csv', index_col=0)
tomneg_clcounts = pd.read_csv('procdata/05script/model_input_tomneg.csv', index_col=0)

# %%
#Calculating fraction of cells in each cluster
totpos = tompos_clcounts.sc_total.copy()
tompos_clcounts = tompos_clcounts.drop(labels=['flow_total', 'start_age', 'time', 'sc_total'], axis=1)
tompos_clcountsN = (tompos_clcounts.T / totpos)

totneg = tomneg_clcounts.sc_total.copy()
tomneg_clcounts = tomneg_clcounts.drop(labels=['flow_total', 'start_age', 'time', 'sc_total'], axis=1)
tomneg_clcountsN = (tomneg_clcounts.T / totneg)

# %%
# Calculating the ratios between Tom+ cluster fractions and Tom- cluster fractions
# The calculation is performed per timepoint, so they need to match the start_age and tx_timepoint
meta['start_age_tx_days'] = meta.start_age + '_' + meta.timepoint_tx_days.astype('str')
tompos_clcountsN_rel = tompos_clcountsN.copy()

for i in tompos_clcountsN.columns:
    print(i)
    #Finding the matching references and calculating mean
    group = meta.loc[meta.biosample_id == i,'start_age_tx_days'][0]  
    r = meta.loc[(meta.start_age_tx_days == group) & (meta.tom == 'neg'),'biosample_id']
    r = tomneg_clcountsN.loc[:,r]
    r = r.mean(axis = 1)
    #Calculating the relative number of cells
    tompos_clcountsN_rel.loc[:,i] = np.log2((tompos_clcountsN_rel.loc[:,i] / r) + 0.01)

# %%
#Plotting the average relative fractions per timepoint
avs = tompos_clcountsN_rel.melt(ignore_index = False, var_name = 'biosample_id')
avs['cluster'] = avs.index
avs['timepoint_tx_days'] = meta.loc[avs.biosample_id, 'timepoint_tx_days'].values

avsum = avs.groupby(['timepoint_tx_days', 'cluster']).mean()
print(avsum)

toplot = []
for i in avsum.index.get_level_values(0).unique():
    ratio = avsum.loc[i,:].copy()
    ratio = ratio.value
    hoxb5.obs[str(i) + ' days'] = [ratio[i] for i in hoxb5.obs.leiden]
    toplot.append(str(i) + ' days')

cmap2 = utils.plots.cmap_RdBu2(values=None, vmin=avsum.value.min(), vmax=avsum.value.max())
sc.pl.umap(hoxb5, color = toplot, cmap = cmap2,
           vmin=avsum.value.min(), vmax=avsum.value.max(),
           save = 'rel_abundance_percluster_average.pdf')

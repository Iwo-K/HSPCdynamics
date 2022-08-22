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

# %% tags=[]
#Singularity container used:
# !echo $SINGULARITY_CONTAINER
# !hostname

# %% [markdown]
# # Plots for figure 2 with projections of external data

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
import utils.proc, utils.plots
import cellproject as cp
import pickle
from utils.fixes import unclog_umap_caching
import plotnine as pn

import warnings
warnings.filterwarnings("ignore", category=FutureWarning) 

sc.settings.verbosity = 3
base_figures = './figures/23script/'
sc.settings.figdir = base_figures
base_procdata = './procdata/23script/'
for i in [base_figures, base_procdata]:
    pathlib.Path(i).mkdir(parents=True, exist_ok=True)

# %% tags=[]
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

# %% tags=[]
# %load_ext autoreload
# %autoreload 2

# %% [markdown] tags=[]
# ## Nestorowa et al. 2016 data projection

# %% tags=[]
hoxb5 = sc.read('procdata/04script/combined_filt.h5ad')
meta = pd.read_csv('procdata/04script/combined_sample.meta.csv', index_col = 0)
meta = meta.sort_values(by ='timepoint_tx_days', ascending = False)
hoxb5.obs['longname'] = [meta.loc[i, 'longname'] for i in hoxb5.obs.biosample_id]

# %% tags=[]
#Loading Nestorowa et al. 2016 data
sfdata = sc.read('./procdata/07script/sfdata_hoxb5projection.h5ad')

# %% [markdown] tags=[]
# ### UMAPs

# %% tags=[]
ref = hoxb5.copy()
ref.obs['celltype'] = 'unassigned'
comb = ref.concatenate(sfdata, batch_key='batch', batch_categories=['ref', 'sfdata'], index_unique=None)
toplot = ['ESLAM', 'STHSC', 'MPP1', 'MPP3', 'LMPP', 'CMP', 'GMP', 'MEP']
# toplot = sfdata.obs.celltype_e.unique()
# toplot = toplot[~toplot.isin(['unassigned', 'MPP2'])]
utils.plots.umap_subgroups(comb, key = 'celltype_e', toplot = toplot, size=48,
                           file=base_figures + 'sfdata_projection_subgroups.pdf')

# %% tags=[]
from utils.proc import subsample_tomax
from scanpy.plotting.palettes import default_20

x = sfdata.copy()
print(x.obs.celltype_e.value_counts())
x = subsample_tomax(x, 'celltype_e', n_obs_max=60)
print(x.obs.celltype_e.value_counts())

# %% tags=[]
fig, ax = plt.subplots()

fig.set_size_inches(5, 4)
sc.pl.umap(hoxb5, alpha =0.35, ax=ax, show=False)
y = x[x.obs.index[x.obs.celltype_e.isin(['ESLAM', 'STHSC', 'MPP1'])],:].copy()
y.uns['celltype_e_colors'] = default_20[0:3]
["#580053", "#698100", "#f224de"]
sc.pl.umap(y,
           color='celltype_e', ax=ax, show=False, size = 55, frameon=False)
# fig.tight_layout()
ax.set_box_aspect(1)
plt.savefig(base_figures + 'sfdata_hoxb5proj_ESLAM_ST_MPP1.pdf')

# %% tags=[]
fig, ax = plt.subplots()
fig.set_size_inches(5, 4)
sc.pl.umap(hoxb5, alpha=0.35, ax=ax, show=False)
y = x[x.obs.index[x.obs.celltype.isin(['MPP3', 'LMPP', 'GMP', 'MEP'])],:].copy()
y.uns['celltype_colors'] = default_20[4:8]

sc.pl.umap(y,
           color='celltype', ax=ax, show=False, size = 55, frameon=False)
ax.set_box_aspect(1)
plt.savefig(base_figures + 'sfdata_hoxb5proj_progenitors.pdf')

# %% [markdown] tags=[]
# ### Frequency of clusters within each immunophenotypic gate

# %%
df = pd.crosstab(sfdata.obs.ref_leiden, sfdata.obs.celltype_e, normalize='columns')
dfm = df.melt(ignore_index=False)
dfm.index.name = None
dfm['ref_leiden'] = dfm.index
dfm = dfm.loc[~dfm.celltype_e.isin(['unassigned', 'MPP2']),:]

dfm['ref_leiden'] = dfm.ref_leiden.astype(int)
print(dfm.sort_values(by='ref_leiden'))

dfm['ref_leiden'] = pd.Categorical(dfm.ref_leiden.astype(str), categories=sfdata.obs.ref_leiden.cat.categories)

# %% tags=[]
lcol = {leiden : color for leiden, color in zip(sfdata.obs.ref_leiden.cat.categories, sfdata.uns['ref_leiden_colors'])}

g = (pn.ggplot(dfm, pn.aes(x='celltype_e', y='value', fill='ref_leiden'))
 + pn.geom_col()
 + pn.scales.scale_fill_manual(lcol, name = "Clusters")
 + pn.xlab('Population')
 + pn.ylab('Fraction of cells')
 + pn.theme_bw(base_size=14)    
    
)
print(g)
g.save(base_figures + 'sfdata_clusterfractions.pdf')

# %% [markdown]
# ## Weinreb et al. 2020 projections

# %%
d2clones = sc.read('./procdata/07script/Weirenb2020_hoxb5projection.h5ad')
fates = d2clones.obs.fateclass.value_counts()

#Keeping fates with at least cells
tokeep = fates[fates > 10].index
tokeep = tokeep[tokeep != 'nofate']
d2clones = d2clones[d2clones.obs.fateclass.isin(tokeep),:].copy()

ref.obs['fateclass'] = 'unassigned'
comb = ref.concatenate(d2clones, batch_key='batch', batch_categories=['ref', 'bhsc'], index_unique=None)

# %%
fig, ax = plt.subplots()
fig.set_size_inches(5, 4)
sc.pl.umap(hoxb5, alpha=0.35, ax=ax, show=False)

y = d2clones[d2clones.obs.fateclass.isin(['Monocyte', 'Neutrophil', 'Monocyte_Neutrophil']),:].copy()
y.uns['fateclass_colors'] = ['#37bded', '#715544', '#fc95ca']
#default_20[0:3]# + [default_20[5]]

sc.pl.umap(y,
           color='fateclass', ax=ax, show=False, size = 35, frameon=False, alpha=0.7)
# fig.tight_layout()
ax.set_box_aspect(1)
plt.savefig(base_figures + 'Weinreb2020_hoxb5proj_MonoNeu.pdf')

# %%
toplot = ['Monocyte', 'Neutrophil', 'Monocyte_Neutrophil', 'Baso',
 'Mast', 'Baso_Neutrophil', 'Baso_Mast', 'Lymphoid',
  'Lymphoid_Monocyte', 'Meg', 'Erythroid_Meg', 'Erythroid']

# %% tags=[]
utils.plots.umap_subgroups(comb, key = 'fateclass', size=48,
                           toplot = toplot, file=base_figures + 'Weinreb2020_projection_subgroups.pdf')

# %% [markdown]
# ## Bowling et al. 2020 projections

# %%
bhsc = sc.read('./procdata/07script/Bowling2020_hoxb5projection.h5ad')

# %%
fig, ax = plt.subplots()
fig.set_size_inches(5, 4)
sc.pl.umap(hoxb5, alpha=0.35, ax=ax, show=False)

sc.pl.umap(bhsc,
           color='HSC_labels', ax=ax, show=False, size = 55, frameon=False, alpha=0.7)
# fig.tight_layout()
ax.set_box_aspect(1)
plt.savefig(base_figures + 'Bowling2020_hoxb5proj_HSClabels.pdf')

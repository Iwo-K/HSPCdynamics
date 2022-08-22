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

# %%
#Singularity container used:
# !echo $SINGULARITY_CONTAINER
# !hostname

# %% [markdown]
# # Combining 10x data

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
import utils.proc
import utils.plots
from utils.load10x import load_10xfiles

sc.settings.verbosity = 3
sc.settings.figdir = './figures/02script/'

base_figures = './figures/02script/'
base_procdata = './procdata/02script/'
for i in [sc.settings.figdir, base_figures, base_procdata]:
    pathlib.Path(i).mkdir(parents=True, exist_ok=True)

# %%
#Figure settings
sc.set_figure_params(color_map='viridis', dpi_save=350)
mpl.rcParams['figure.figsize'] = (4, 4) # 1:1 plot ratio
mpl.rc('xtick', labelsize=14) 
mpl.rc('ytick', labelsize=14) 
mpl.rc('axes', labelsize=16) 
mpl.rc('axes', labelsize=16)
mpl.rcParams['pdf.fonttype'] = 42 # Ensures readable fonts in illustrator
mpl.rcParams['ps.fonttype'] = 42
plt.rc('axes', axisbelow=True) # Ensure that gridlines are behind points
mplparams = mpl.rcParams.copy()

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# Loading genesets

# %%
ygenes = pd.read_csv('./data/genesets/Y_genes.csv')
realY_filter = [not bool(re.search("predicted gene", i)) for i in ygenes['Gene description']]
ygenes = ygenes.loc[realY_filter,:]
ygenes = ygenes['Gene name']

ccgenes = pd.read_csv('./data/genesets/CC_genes_forremoval.txt', header = None)
ccgenes = ccgenes[0]

toremove = ccgenes.copy()
toremove = toremove.append(pd.Series(['Xist']))
toremove = toremove.append(ygenes)

#Getting cell cycle genes per phase
cycle_genes = pd.read_csv('./data/genesets/macosko2015_hum_ccgenes.csv') #These are human cell cycle genes with just names swapped for mouse
Sgenes = cycle_genes.loc[cycle_genes['phase'] == 'S','gene']
G2Mgenes = cycle_genes.loc[cycle_genes['phase'] == 'G2_M','gene']

# %%
meta = pd.read_csv('./data/10x/scB5Tom_10x_samplemeta.csv')
meta.head()

# %% tags=[]
adatas = load_10xfiles(meta, ygenes, qc = True)

# %% tags=[]
import anndata
comb = anndata.concat(adatas, label = 'biosample_id', merge = 'same')
comb.obs_names_make_unique()
comb.obs['cellid'] = comb.obs.index

#Merging with sample metadata
x = comb.obs.merge(meta, left_on = 'biosample_id', right_on = 'biosample_id', suffixes = ['_adata', '_meta'])
x.index = x.cellid.values
if((comb.obs['cellid'] == comb.obs.index).all()):
    comb.obs = x
else: print("something wrong")   

#Removing all 0 genes
comb = comb[:, comb.var.index[(comb.X.sum(axis=0).A1 > 0)]].copy()
#Checkpoint
comb.write(base_procdata + 'combined10x_qced_logn.h5ad', compression = 'lzf')

# %% [markdown]
# ## All sample overview - checking for batch effect, doublets etc

# %% [markdown]
# ### All samples

# %% tags=[]
comb = sc.read(base_procdata + 'combined10x_qced_logn.h5ad')
comb_proc = comb.copy()
utils.proc.process(comb_proc,
                  compute_hivar = True,
                  n_pcs = 50,
                  n_neighbors = 12,
                  lognorm_data = False,
                  n_variable = 5000,
                  remove_genes = toremove,
                  Sgenes = Sgenes,
                  G2Mgenes = G2Mgenes, 
                  regress_cc = False)

# %% tags=[]
comb_proc.obs.predicted_doublets = pd.Categorical(comb_proc.obs.predicted_doublets)
sc.pl.umap(comb_proc, color = ['biosample_id', 'phase', 'leiden', 'predicted_doublets', 'doublet_scores'], save = '_allsamples_info.pdf', wspace = 1.6)
sc.pl.umap(comb_proc, color = ['Procr', 'Klf1', 'Dntt', 'Elane', 'Irf8', 'Cma1', 'Pf4'], save = '_allsamples_markers.pdf')
utils.plots.umap_subgroups(comb_proc, key = 'sample_id', toplot = comb.obs.sample_id.cat.categories, file = base_figures + 'umap_bysampleid.pdf')

# %% tags=[]
utils.plots.umap_subgroups(comb_proc, key = 'biosample_id', toplot = comb.obs.biosample_id.cat.categories, file = base_figures + 'umap_by_biosampleid.pdf')

# %%
# Removing the SLX-20620_SITTA3, which has technical issues (diagnosed as partial wetting failure)
comb = comb[comb.obs.index[comb.obs.sample_id != '7w11w_SITTA3'],:].copy()

# %% [markdown]
# ### Doublets

# %%
# Removing potential doublets
doubs = pd.crosstab(comb.obs.sample_id, [comb.obs.predicted_doublets])
doubs['ncells'] = comb.obs.sample_id.value_counts()[doubs.index]
doubs['ratio'] = doubs[True]/(doubs[True]+doubs[False])
print(doubs)

# %%
comb = comb[comb.obs.index[comb.obs.predicted_doublets != True],:].copy()

#Removing all 0 genes
sc.pp.filter_genes(comb, min_counts=1)
comb.write(base_procdata + 'combined10x_qced_logn.h5ad')

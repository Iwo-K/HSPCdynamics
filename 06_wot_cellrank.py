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
# # Trajectory analysis

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
import utils.proc, utils.plots
import scvelo as scv
from cellrank.tl.estimators import GPCCA

import cellrank as cr

sc.settings.verbosity = 3
sc.settings.figdir = './figures/06script/'
base_figures = './figures/06script/'
base_procdata = './procdata/06script/'
for i in [sc.settings.figdir, base_figures, base_procdata]:
    pathlib.Path(i).mkdir(parents=True, exist_ok=True)

# %%
#Figure settings
sc.set_figure_params(dpi=100, frameon=False, figsize=(4,4), color_map='viridis', dpi_save=350)
scv.settings.set_figure_params("scvelo")
scv.set_figure_params(dpi=100, frameon=False, figsize=(4,4), dpi_save=350)

scv.settings.set_figure_params("scvelo")
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
# ## Loading data

# %%
comb = sc.read('procdata/04script/combined_filt.h5ad')
meta = pd.read_csv('procdata/04script/combined_sample.meta.csv', index_col = 0)
meta = meta.sort_values(by ='timepoint_tx_days', ascending = False)
comb.obs['longname'] = [meta.loc[i, 'longname'] for i in comb.obs.biosample_id]

# %%
sc.pl.umap(comb, color = 'leiden', legend_loc = 'on data')

# %% tags=[]
comb = comb[comb.obs.tom == 'pos', :].copy() #using only Tom+ cells
#Recomputing neighbors
sc.pp.neighbors(comb, n_neighbors=12, use_rep='X_pca_harmony')

# %% [markdown]
# ## Pseudotime Kernel

# %%
# use scVelo's `moments` function for imputation - note that hack we're using here:
# we're copying our `.X` matrix into the layers because that's where `scv.tl.moments`
# expects to find counts for imputation
comb.layers["spliced"] = comb.X
comb.layers["unspliced"] = comb.X
scv.pp.moments(comb, n_pcs=30, n_neighbors=30, use_rep='X_pca_harmony')

# %% tags=[]
dptk = cr.tl.kernels.PseudotimeKernel(comb)
dptk.compute_transition_matrix(threshold_scheme="soft", nu=0.5)

dptk.compute_projection(basis="umap", key_added='T_fwd')
scv.pl.velocity_embedding_stream(
    comb, color="dpt_pseudotime", vkey="T_fwd", basis="umap", legend_loc="right", cmap='jet_r',
    save=base_figures+'pseudotime_kernel_stream.pdf', dpi=350
)

scv.pl.velocity_embedding_stream(
    comb, color="dpt_pseudotime", vkey="T_fwd", basis="umap", legend_loc="right", cmap='jet_r', vmax=0.55,
    save=base_figures+'pseudotime_kernel_stream_vmax055.pdf', dpi=350
)


# %% [markdown]
# ### Terminal states
#
# %% tags=[]
g = GPCCA(dptk)
g.compute_schur(n_components=20)
g.plot_spectrum(real_only=True)

# %%
g.compute_macrostates(n_states=12, cluster_key="leiden")
g.plot_macrostates(discrete=True, basis="umap", legend_loc="right")

# %%
sc.pl.umap(comb, color = 'leiden', legend_loc = 'on data')

# %% [markdown]
# We will set terminal states manually, they are largely consistent with the macrostates identified with GPCCA

# %%
def get_terminal_states(adata, clusters, cell_no=50, clusters_key='leiden',
                        pseudotime_key='dpt_pseudotime'):
    cells = pd.Series(np.nan, index=adata.obs_names)
    for i in clusters:
        cluster_adata = adata[adata.obs[clusters_key] == i, :]
        topT = cluster_adata.obs[pseudotime_key].sort_values(ascending=False)
        topT = topT.index[:cell_no].values
        cells[topT] = i
    return(pd.Series(cells, dtype="category"))

fates = ['20', '7', '12', '25', '26', '6', '24', '14', '16', '28', '11', '10']
ends = get_terminal_states(comb, clusters=fates)

# %%
g.set_terminal_states(ends)

# %%
g.compute_absorption_probabilities(solver="gmres", tol=1e-9, preconditioner='ilu')
#This sometimes works and sometimes does not (values not summing up to 1)

# %%
scv.settings.set_figure_params("scvelo")
scv.set_figure_params(dpi=100, frameon=False, figsize=(4,4), dpi_save=350)
scv.settings.set_figure_params("scvelo")

g.plot_absorption_probabilities(same_plot=False, basis="umap",
                                dpi=350,
                                wspace=0.5,
                                perc=[0, 99],
                                save=base_figures + 'pseudotime_kernel_fates.pdf')

# %%
fates = pd.DataFrame(g.absorption_probabilities.X, columns = g.absorption_probabilities.names)
fates.to_csv(base_procdata + 'pseudotime_kernel_fates.csv')

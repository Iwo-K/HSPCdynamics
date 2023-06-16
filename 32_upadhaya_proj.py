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

# %%
#Singularity container used:
# !echo $SINGULARITY_CONTAINER
# !hostname

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

# %%
sc.pl.umap(hoxb5,color = 'leiden',legend_loc='on data')

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

# %%
hoxb5_store = hoxb5.copy()
hoxb5_store.obs = pd.concat([hoxb5.obs,hoxb5_pos.obs.iloc[:,-10:]],axis=1)

# %% [markdown]
# ## Upadhaya 2018

# %%
counts1 = pd.read_csv('reizis2018/GSE120239_RNA1_counts.csv',index_col=0)

# %%
counts2 = pd.read_csv('reizis2018/GSE120239_RNA2_counts.csv',index_col=0)

# %%
counts = pd.concat([counts1,counts2],axis = 1,join = 'outer').fillna(0)

# %%
counts.index = [name.capitalize() for name in counts.index]

# %%
time = [name.split('-')[0] for name in counts.columns]
batch = [name.split('-')[1] + '-' + name.split('-')[1] for name in counts.columns]

# %%
reizis = sc.AnnData(counts.T)

# %%
reizis.obs['time_reizis'] = pd.Categorical(time)
reizis.obs['batch_reizis'] = pd.Categorical(batch)

# %%
reizis.var['mt'] = reizis.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(reizis, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

# %%
sc.pp.filter_cells(reizis, min_genes=1400)
sc.pp.filter_cells(reizis, max_genes=9000)

reizis = reizis[reizis.obs.pct_counts_mt < 8, :]

# %%
reizis.obs.batch_reizis.value_counts()

# %%
reizis.obs.time_reizis.value_counts()

# %%
utils.proc.lognorm(reizis)

# %%
reizis_basic = reizis.copy()

# %% [markdown]
# Projection

# %%
ref = hoxb5.copy()
#Unifying genes
ref = ref[:, ref.var.index[ref.var.index.isin(reizis.var.index)]].copy()
ref.X = ref.raw[:,ref.var.index].X.copy()
del ref.raw
reizis = reizis[:, ref.var.index].copy()

ref.obs['time_reizis'] = 'unassigned'
comb = ref.concatenate(reizis, batch_key='batch', batch_categories=['ref', 'reizis'], index_unique=None)

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
utils.plots.umap_subgroups(comb, key = 'time_reizis',toplot = reizis.obs.time_reizis.unique())

# %%
sc.pl.umap(comb, color=['batch', 'Procr','leiden'],s = 10)
utils.plots.umap_subgroups(comb, key = 'time_reizis',toplot = reizis.obs.time_reizis.unique())

# %%
sc.pl.umap(comb, color=['Ly6a', 'Csf1r','Pf4'],s = 10)

# %%
reizis.obsm['X_pca'] = comb[reizis.obs.index,:].obsm['X_pca_harmony'].copy()
ref.obsm['X_pca'] = comb[ref.obs.index,:].obsm['X_pca_harmony'].copy()

# %%
cp.project_cells(reizis, ref,
                 obs_columns=['leiden'],
                 fit_pca=False,
                 scale_data=False,
                k=8)


#Swapping for the original PCs, based on cell harmony-corrected data
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()
#nn-regressing the PCs
cp.nnregress(reizis,
             ref,
             regress=['pca'],
             weighted=True)

#fitting into the original umap
cp.project_cells(reizis, ref,
                 obs_columns=['leiden'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)

# %%
reizis_basic.obs = pd.concat([reizis_basic.obs,reizis.obs.ref_leiden],axis=1)

# %%
reizis_basic.write(base_procdata + 'Reizis2018_hoxb5projection.h5ad', compression='lzf')

# %% [markdown]
# Cluster aboundance

# %%
# Counting cells in each cluster seprately for Tom+ and Tom- cells
# pop = comb[comb.obs.batch == 'reizis',:].copy()

pop_leiden = pd.crosstab(reizis_basic.obs.ref_leiden, reizis_basic.obs.time_reizis)
pop_leiden.to_csv(base_procdata + 'reizis_cluster_cellcounts.csv')

# %% [markdown]
# Project dpt

# %%
reizis.obsm['X_pca'] = comb[reizis.obs.index,:].obsm['X_pca_harmony'].copy()
ref.obsm['X_pca'] = comb[ref.obs.index,:].obsm['X_pca_harmony'].copy()

# %%
project_cells(reizis, ref,
                 obs_columns=['dpt_pseudotime'],
                 fit_pca=False,
                 scale_data=False,
                k=8)

#Swapping for the original PCs, based on cell harmony-corrected data
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()
#nn-regressing the PCs
cp.nnregress(reizis,
             ref,
             regress=['pca'],
             weighted=True)

#fitting into the original umap
project_cells(reizis, ref,
                 obs_columns=['dpt_pseudotime'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)

# %%
reizis_basic.obs = pd.concat([reizis_basic.obs,reizis.obs.ref_dpt_pseudotime],axis=1)

# %%
reizis_basic.write(base_procdata + 'Reizis2018_hoxb5projection.h5ad', compression='lzf')

# %% [markdown]
# Replot

# %%
#reizis = sc.read(base_procdata + 'Reizis2018_hoxb5projection.h5ad')
comb_plot = ref.concatenate(reizis, batch_key='batch', batch_categories=['ref', 'reizis'], index_unique=None)

# %%
sc.pl.umap(ref,color='leiden')
sc.pl.umap(reizis,color=['time_reizis'],size = 30,save = '_Reizis2018_projection2.pdf')
#upload this
utils.plots.umap_subgroups(comb_plot, key = 'time_reizis',toplot = reizis.obs.time_reizis.unique(),file=base_figures + 'Reiz2018_projection_subgroups2.pdf')

# %%

# %%
sc.pl.umap(reizis, color='time_reizis')

# %%
sc.pl.umap(reizis, color=['Pf4', 'Procr',  'Elane', 'Csf1', 'Dntt'])#. recover Gata1 or Klf1

# %%
reizis_sub = reizis[reizis.obs.time_reizis == 'd03'].copy()

# %%
sc.pl.umap(reizis_sub, color=['Pf4', 'Procr',  'Elane', 'Csf1', 'Dntt'])#. recover Gata1 or Klf1

# %%
sc.pl.umap(reizis,color=['ref_leiden'],size = 30)

# %%
sc.pl.umap(reizis,color=['ref_dpt_pseudotime'],size = 30)

# %%
sc.pl.umap(ref,color=['dpt_pseudotime'],size = 30)

# %% [markdown]
# Project on label

# %%
reizis.obsm['X_pca'] = comb[reizis.obs.index,:].obsm['X_pca_harmony'].copy()
ref.obsm['X_pca'] = comb[ref.obs.index,:].obsm['X_pca_harmony'].copy()

# %%
ref = ref[hoxb5_pos.obs_names].copy()

# %%
ref.obs = hoxb5_pos.obs.copy()

# %%
cp.project_cells(reizis, ref,
                 obs_columns=['fate_6','fate_7','fate_10','fate_11','fate_16'],
                 fit_pca=False,
                 scale_data=False,
                k=8)


#Swapping for the original PCs, based on cell harmony-corrected data
ref.obsm['X_pca'] = ref.obsm['X_pca_harmony'].copy()
#nn-regressing the PCs
cp.nnregress(reizis,
             ref,
             regress=['pca'],
             weighted=True)

#fitting into the original umap
cp.project_cells(reizis, ref,
                 obs_columns=['fate_6','fate_7','fate_10','fate_11','fate_16'],
                 fit_pca=False,
                 scale_data=False,
                 umap_ref=umapref)

# %%
sc.pl.umap(reizis,color=['ref_fate_6','ref_fate_7','ref_fate_10','ref_fate_11','ref_fate_16'],size = 30)

# %%
reizis_basic.obs = pd.concat([reizis_basic.obs,reizis.obs.iloc[:,-5:]],axis=1)

# %%
reizis_basic.write(base_procdata + 'Reizis2018_hoxb5projection.h5ad', compression='lzf')

# %%


# %%
reizis.obsm['X_pca'] = comb[reizis.obs.index,:].obsm['X_pca_harmony'].copy()
ref.obsm['X_pca'] = comb[ref.obs.index,:].obsm['X_pca_harmony'].copy()

# %%

# %%
for clu in [6, 10, 11, 16, 7]:

    reizis_temp = reizis[reizis.obs['ref_fate_'+str(clu)] == True].copy()
    ref_temp = ref[ref.obs['fate_'+str(clu)] == True].copy()

    project_cells(reizis_temp, ref_temp,
                 obs_columns=['dpt_fate'+str(clu)],
                 fit_pca=False,
                 scale_data=False,
                 k=8)
    


    #Swapping for the original PCs, based on cell harmony-corrected data
    ref_temp.obsm['X_pca'] = ref_temp.obsm['X_pca_harmony'].copy()
    #nn-regressing the PCs
    cp.nnregress(reizis_temp,
                 ref_temp,
                 regress=['pca'],
                 weighted=True)

    #fitting into the original umap
    project_cells_mela(reizis_temp, ref_temp,
                     obs_columns=['dpt_fate'+str(clu)],
                     fit_pca=False,
                     scale_data=False,
                     umap_ref=umapref)
    
    reizis.obs['ref_dpt_fate'+str(clu)] = reizis_temp.obs['ref_dpt_fate'+str(clu)]

# %%
sc.pl.umap(reizis,color=['ref_dpt_fate6','ref_dpt_fate7','ref_dpt_fate10','ref_dpt_fate11','ref_dpt_fate16'],size = 30)

# %%
reizis_basic.obs = pd.concat([reizis_basic.obs,reizis.obs.iloc[:,-5:]],axis=1)

# %%
reizis_basic.write(base_procdata + 'Reizis2018_hoxb5projection.h5ad', compression='lzf')

# %%

# %% [markdown]
# Create input Pseudodynamics 

# %%
reizis = sc.read(base_procdata + 'Reizis2018_hoxb5projection.h5ad')

# %%
df = pd.DataFrame(list(reizis.obs.time_reizis),columns=['stage'],index = reizis.obs_names)

# %%
df['sample'] = reizis.obs.batch_reizis

# %%
for clu in [6, 10, 11, 16, 7]:
    df_temp = df.loc[reizis.obs['ref_fate_' + str(clu)] == 'True'].copy()
    df_temp['dpt_pseudotime'] = reizis.obs['ref_dpt_fate' + str(clu)]
    df_temp.to_csv(base_procdata + 'reizis_input_pseudo_dyn_clu_' + str(clu) + '_dpt.csv')

# %%
df6 = pd.read_csv(base_procdata + 'reizis_input_pseudo_dyn_clu_' + str(6) + '_dpt.csv')

# %%
df7 = pd.read_csv(base_procdata + 'reizis_input_pseudo_dyn_clu_' + str(7) + '_dpt.csv')

# %%
df11 = pd.read_csv(base_procdata + 'reizis_input_pseudo_dyn_clu_' + str(10) + '_dpt.csv')

# %%
plt.hist(df11[df11.stage == 'd03'].dpt_pseudotime,bins = 20,density=True);
plt.hist(df11[df11.stage == 'd10'].dpt_pseudotime,bins = 20,density=True);
plt.hist(df11[df11.stage == 'd17'].dpt_pseudotime,bins = 20,density=True);

# %%
plt.hist(df7[df7.stage == 'd03'].dpt_pseudotime,bins = 20,density=True);
plt.hist(df7[df7.stage == 'd10'].dpt_pseudotime,bins = 20,density=True);
plt.hist(df7[df7.stage == 'd17'].dpt_pseudotime,bins = 20,density=True);


# %%
plt.hist(df6[df6.stage == 'd03'].dpt_pseudotime,bins = 10,density=True);
plt.xlim([0,1])
# plt.show()
plt.hist(df6[df6.stage == 'd10'].dpt_pseudotime,bins = 10,density=True);
plt.xlim([0,1])

# plt.show()

plt.hist(df6[df6.stage == 'd17'].dpt_pseudotime,bins = 30,density=True);
plt.xlim([0,1])

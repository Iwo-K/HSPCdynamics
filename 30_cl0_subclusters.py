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
# # Subclustering cluster 0

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
adata = sc.read('procdata/04script/combined_filt.h5ad')
meta = pd.read_csv('procdata/04script/combined_sample.meta.csv', index_col = 0)
meta = meta.sort_values(by ='timepoint_tx_days', ascending = False)
adata.obs['longname'] = [meta.loc[i, 'longname'] for i in adata.obs.biosample_id]

# sc.pp.subsample(hoxb5, n_obs=5000)
#Getting the UMAP object
# with open('procdata/04script/combined_filt_umapref.pkl', 'rb') as f:
#     umapref = pickle.load(f)

sc.pl.umap(adata)
adata_pos = adata[adata.obs.tom == 'pos'].copy()

# %%
adata_10x = adata[adata.obs.data_type == '10x'].copy() 

# %% [markdown]
# ## Cluster 0: reclustering

# %%
adata0 = adata_10x[adata_10x.obs.leiden == '0'].copy()

# %%
del adata0.obsm
del adata0.varm
del adata0.obsp

# %%
adata0 = adata0.raw.to_adata()
adata0_pos = adata0[adata0.obs.tom == 'pos'].copy()

# %%
# Selecting genes which should be excluded from the highly-variable gene list
ygenes = pd.read_csv('../../human_preimplantation/data_input/Y_genes_iwo.csv')
realY_filter = [not bool(re.search("predicted gene", i)) for i in ygenes['Gene description']]
ygenes = ygenes.loc[realY_filter,'Gene name']
toremove =  ygenes.append(pd.Series(['Xist']))

# Using method from Weinreb et al. 2020, finding genes which are correlated with
# cell cycle signature, at least 0.1 Pearson R
ccgenes = utils.proc.find_correlated_genes(adata0,
                      data_scaled = False,
                      genes = ['Ube2c', 'Hmgb2', 'Hmgn2', 'Tuba1b', 'Ccnb1', 'Tubb5', 'Top2a', 'Tubb4b'])
toremove = toremove.append(pd.Series(ccgenes))

# %%
sc.pp.highly_variable_genes(adata0,min_mean=0.03,min_disp=0.3)
sc.pl.highly_variable_genes(adata0)
np.sum(adata0.var.highly_variable)

# %%
print("Removing: " + str(len(toremove)) + " genes")
adata0.var['gene_removed'] = adata0.var.index.isin(toremove)
adata0.var.loc[adata0.var.gene_removed, adata0.var.columns == 'highly_variable'] = False
print(sum(adata0.var.highly_variable))

# %%
adata0.raw = adata0
adata0 = adata0[:,adata0.var.highly_variable].copy()

# %%

# %%
sc.pp.scale(adata0)

# %%
sc.tl.pca(adata0, svd_solver='arpack')

# %% tags=[]
adata0.obs['time_cat'] = 't_' + adata0.obs.timepoint_tx_days.astype('str')

# %%
sc.pl.pca(adata0,color = ['time_cat','phase'],s = 30)


# %%
sc.external.pp.harmony_integrate(adata0, key='biosample_id')#do it for pos and neg
#Correcting for batch effects
sc.pp.neighbors(adata0, n_neighbors=15, use_rep='X_pca_harmony')

#How well does the integation look like?
umapref_t = cp.quick_umap(adata0, use_rep='X_pca_harmony')
unclog_umap_caching(umapref_t)
umapref = cp.quick_umap(adata0, use_rep='X_pca_harmony')

# with open(base_procdata + 'adata0_umapref.pkl', 'wb') as f:
#     pickle.dump(umapref, f)

sc.pl.umap(adata0, color=['batch', 'Procr'])


# %%
adata0.obs['time_cat'] = 't_' + adata0.obs.timepoint_tx_days.astype('str')

# %%
sc.pl.umap(adata0,ncols = 2,color = ['tom','biosample_id','batch','Procr','Ly6a','HSCscore'],wspace=.4, hspace = .5,save = 'info_subclusters.png')

# %%
sc.tl.leiden(adata0,resolution=1.2)
sc.pl.umap(adata0,color = ['leiden'],legend_loc = 'on data',save = 'leiden_subclusters.png')
#1.3

# %%
sc.pl.umap(adata0,color = ['leiden'],legend_loc = 'on data')

# %%
# adata0.write(base_procdata + 'adata0_for_subcusters.h5ad')

# %%
ll = adata0[adata0.obs.tom == 'neg'].obs.leiden.value_counts()/adata0[adata0.obs.tom == 'neg'].shape[0]

# %%

# %%
plt.figure(figsize =(9,9))

plt.bar(ll.index,adata0[adata0.obs.tom == 'neg'].obs.leiden.value_counts()/adata0[adata0.obs.tom == 'neg'].shape[0])
plt.hlines(0.03,-1,14,colors='r')
plt.hlines(0.023,-1,14,colors='k')
plt.legend(['bestfit','max/min bounds'])


plt.hlines(0.034,-1,14,colors='k')

plt.xlabel('subcluster')
plt.ylabel('relative size')

plt.savefig('./figures/relative_size_subclusters.png')

# %%

# %%
sc.pl.umap(adata0, color=['tom'])

# %%
sc.pl.violin(adata0,keys = ['Procr'],groupby='leiden',wspace = 10,save = 'Procr_subclusters.png')
sc.pl.violin(adata0,keys = ['Ly6a'],groupby='leiden',wspace = 10,save = 'Ly6a_subclusters.png')

# sc.pl.violin(adata0,keys = ['Vwf'],groupby='leiden',wspace = 10)
sc.pl.violin(adata0,keys = ['HSCscore'],groupby='leiden',wspace = 10,save = 'HSCS_score_subclusters.png')

# %%
adata0_pos.obsm['X_umap'] = adata0[adata0_pos.obs_names].obsm['X_umap'].copy()

# %%
adata0_pos.obs['leiden'] = adata0[adata0_pos.obs_names].obs['leiden'].copy()

# %%
sc.pl.umap(adata0,color = ['phase'],size = 30)

# %% [markdown]
# # Prepare input absolute numbers

# %%
# positive

# %%
df_pos = pd.DataFrame(list(np.unique(adata0[adata0.obs.tom == 'pos'].obs.timepoint_tx_days)),columns=['day'],index = None)

# %%
df1 = pd.read_csv('./procdata/model_input_tompos.csv',index_col=0)
df1.head()

# %%
df1 = df1.dropna()

# %%
df2 = pd.crosstab(adata0[adata0.obs.tom == 'pos'].obs.biosample_id,adata0[adata0.obs.tom == 'pos'].obs.leiden)

# %%
df2.head()

# %%
name_list = [name for name in df2.index if name in df1.index]

# %%
df2 = df2.loc[name_list]

# %%
df2.columns = df2.columns.add_categories(['actual_cells'])

# %%
df2['actual_cells'] = np.sum(df2.iloc[:,:len(np.unique(adata0.obs.leiden))],axis=1)/df1['sc_ncells_total']*df1['flow_total']

# %%
df2['tot_leiden'] = np.sum(df2.iloc[:,:len(np.unique(adata0.obs.leiden))],axis=1)

# %%
df2.head()

# %%
mean_pos = np.array([[-11 for i in range(len(np.unique(adata0.obs.leiden)))] for j in range(len(np.unique(df1.time)))],dtype = 'float') 
std_pos = np.array([[-11 for i in range(len(np.unique(adata0.obs.leiden)))]for j in range(len(np.unique(df1.time)))],dtype = 'float') 
k = -1
for j in np.unique(df1.time):
    k+=1
    df_temp = df2.loc[df1.time == j].copy()

#     for i in range(len(np.unique(adata0.obs.leiden))):
#         h+=1
    
        
    mean_pos[k,:] = np.mean((df_temp.iloc[:,:len(np.unique(adata0.obs.leiden))]).T/(df_temp.tot_leiden) * df_temp.actual_cells,axis = 1)
    std_pos[k,:] = np.std((df_temp.iloc[:,:len(np.unique(adata0.obs.leiden))]).T/(df_temp.tot_leiden) * df_temp.actual_cells,axis = 1)/np.sqrt(np.sum(df1.time == j))

# %%
# negative

# %%
df_neg = pd.DataFrame(list(np.unique(adata0[adata0.obs.tom == 'neg'].obs.timepoint_tx_days)),columns=['day'],index = None)

# %%
df1 = pd.read_csv('./procdata/model_input_tomneg.csv',index_col=0)
df1.head()

# %%
df1 = df1.dropna()

# %%
df2 = pd.crosstab(adata0[adata0.obs.tom == 'neg'].obs.biosample_id,adata0[adata0.obs.tom == 'neg'].obs.leiden)

# %%
df2.head()

# %%
name_list = [name for name in df2.index if name in df1.index]

# %%
df2 = df2.loc[name_list]

# %%
df2.columns = df2.columns.add_categories(['actual_cells'])

# %%
df2['actual_cells'] = np.sum(df2.iloc[:,:len(np.unique(adata0.obs.leiden))],axis=1)/df1['sc_ncells_total']*df1['flow_total']

# %%
df2['tot_leiden'] = np.sum(df2.iloc[:,:len(np.unique(adata0.obs.leiden))],axis=1)

# %%
df2.head()

# %%
df1.shape

# %%
df2.shape

# %%
(df2.iloc[:,:len(np.unique(adata0.obs.leiden))]).T/(df2.tot_leiden)

# %%
mean_neg = np.array([[-11 for i in range(len(np.unique(adata0.obs.leiden)))] for j in range(len(np.unique(df1.time)))],dtype = 'float') 
std_neg = np.array([[-11 for i in range(len(np.unique(adata0.obs.leiden)))]for j in range(len(np.unique(df1.time)))],dtype = 'float') 
k = -1
for j in np.unique(df1.time):
    k+=1
    df_temp = df2.loc[df1.time == j].copy()


    
    mean_neg[k,:] = np.mean((df_temp.iloc[:,:len(np.unique(adata0.obs.leiden))]).T/(df_temp.tot_leiden) * df_temp.actual_cells,axis = 1)
    std_neg[k,:] = np.std((df_temp.iloc[:,:len(np.unique(adata0.obs.leiden))]).T/(df_temp.tot_leiden) * df_temp.actual_cells,axis = 1)/np.sqrt(np.sum(df1.time == j))

# %%
plt.figure(figsize =(9,9))
err = np.sqrt( (std_pos[:,8]/mean_pos[:,8])**2 +  (std_neg[:,8]/mean_neg[:,8])**2  ) * mean_pos[:,8]/mean_neg[:,8]
plt.errorbar(np.unique(df1.time),mean_pos[:,8]/mean_neg[:,8],err)
plt.hlines(0.023,0,300,colors='r')
plt.hlines(0,0,300,colors='k')
# plt.xticks()

plt.legend(['bestfit','max/min bounds','subcluster 8'])
plt.hlines(0.049,0,300,colors='k')
plt.xlabel('time')
plt.ylabel('labelling frequency')

plt.savefig('./figures/lab_frequency_sub8.png')


# %%
plt.figure(figsize =(9,9))
err = np.sqrt( (std_pos[:,12]/mean_pos[:,12])**2 +  (std_neg[:,12]/mean_neg[:,12])**2  ) * mean_pos[:,12]/mean_neg[:,12]
plt.errorbar(np.unique(df1.time),mean_pos[:,12]/mean_neg[:,12],err)
plt.hlines(0.023,0,300,colors='r')
plt.hlines(0,0,300,colors='k')
# plt.xticks()

plt.legend(['bestfit','max/min bounds','subcluster 12'])
plt.hlines(0.049,0,300,colors='k')
plt.xlabel('time')
plt.ylabel('labelling frequency')

plt.savefig('./figures/lab_frequency_sub12.png')

# %%

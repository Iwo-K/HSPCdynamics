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
# ## Create input mKO2 

# %% tags=[]
import numpy as np
import pandas as pd
import h5py
import scanpy as sc
import anndata as ad
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

# %% [markdown]
# Setup

# %% tags=[]
dir_data = '/Users/mb2338/Desktop/hoxb5_fate_mapping/mko_2023/data/'

# %% tags=[]
cmap = LinearSegmentedColormap.from_list(name = 'gene_cmap', colors = ['lightgrey', 'thistle', 'red', 'darkred'])

# %% [markdown]
# Read files

# %% tags=[]
counts_temp = pd.read_csv(dir_data + 'scB5_star_fcounts.txt' , sep = '\t',skiprows = 1, index_col = 0)

# %% tags=[]
counts = counts_temp.iloc[:,5:]

# %% tags=[]
counts.shape

# %% tags=[]
summary = pd.read_csv(dir_data + 'scB5_star_fcounts_summary.txt' , sep = '\t', index_col = 0)

# %% tags=[]
summary.shape

# %% tags=[]
meta_data = pd.read_csv(dir_data + 'scb5_meta2.csv' , index_col = 0)

# %% tags=[]
index_new = counts.columns.str.replace('./STARout/SLX-','SLX.')

# %% tags=[]
index_new2 = index_new.str.slice(0,19)

# %% tags=[]
meta_data = meta_data.reindex(index_new2)

# %% tags=[]
meta_data

# %% tags=[]
meta_data.shape

# %% tags=[]
summary.columns = index_new2

# %% tags=[]
counts.columns = index_new2

# %% tags=[]
ensembl_gene_table = pd.read_csv(dir_data + 'ensembl_93.txt', sep=',', index_col=0)

# %% tags=[]
gene_names = ensembl_gene_table.loc[counts.index[:54232]]['Gene name']

# %% tags=[]
total = np.sum(summary, axis = 0)
mito_genes = np.where(gene_names.str.startswith('mt-'))
mito_number=  np.sum(counts.iloc[mito_genes[0],:],axis = 0)
nuclear_number = np.sum(counts.iloc[:54232,:],axis = 0)-np.sum(counts.iloc[mito_genes[0],:],axis=0)
exonic_number = np.sum(counts.iloc[:54232,:],axis = 0)
spike = np.sum(counts.iloc[54232:,:],axis = 0)

# %% [markdown] tags=[]
# Create adata

# %% tags=[]
datatr = counts.T
adata = ad.AnnData(datatr.iloc[:,:54232])

# %% tags=[]
adata

# %% tags=[]
adata_var_ids = adata.var_names
adata.var_names = gene_names
adata.var['ids'] = adata_var_ids
adata_obs_ids = datatr.index
adata.obs_names =  adata_obs_ids

# %% tags=[]
adata.var_names_make_unique()

# %% tags=[]
adata.obs['exonic_number'] = np.array(exonic_number)
adata.obs['mitotic_number'] = np.array(mito_number)
adata.obs['nuclear_number'] = np.array(nuclear_number)
adata.obs['total_number'] = np.array(total)
adata.obs['spike'] = np.array(spike)
adata.obs['well'] = np.array(meta_data['well'])
adata.obs['mouse'] = np.array(meta_data['mouse'])
adata.obs['population'] = np.array(meta_data['population'])
adata.obs['mKO'] = np.array(meta_data['mKO'])
adata.obs['plate'] = np.array(meta_data['plate'])
groups = np.array(['HSC+' for i in range(adata.shape[0])], dtype='U8')
indices = np.where((adata.obs['population'] == 'LTHSC') & (adata.obs['mKO'] == 'mKOneg'))
groups[indices] = 'HSC-'
indices = np.where((adata.obs['population'] == 'STHSC') & (adata.obs['mKO'] == 'mKOpos'))
groups[indices] = 'ST+'
indices = np.where((adata.obs['population'] == 'STHSC') & (adata.obs['mKO'] == 'mKOneg'))
groups[indices]='ST-'
adata.obs['groups'] = groups

# %% tags=[]
sc.pp.filter_cells(adata,min_genes = 0)
sc.pp.filter_genes(adata,min_cells = 0)

# %% tags=[]
mito_genes = adata.var_names.str.startswith('mt-')
adata.obs['percent_mito'] = list(np.sum(adata[:, mito_genes].X, axis = 1) / np.sum(adata.X, axis = 1))
adata.obs['n_counts'] = adata.X.sum(axis = 1)

# %% tags=[]
sc.pl.scatter(adata, x = 'n_counts', y = 'n_genes')

# %% tags=[]
plt.scatter(np.log10(adata.obs['total_number']),np.log10(adata.obs['nuclear_number']))

# %% tags=[]
plt.scatter(np.log10(adata.obs['total_number']),np.log10(adata.obs['exonic_number']))

# %% tags=[]
plt.scatter(np.log10(adata.obs['total_number']),np.log10(adata.obs['mitotic_number']))

# %% tags=[]
plt.scatter(np.log10(adata.obs['total_number']),np.log10(adata.obs['spike']))

# %% tags=[]
exonic_fraction = (adata.obs['mitotic_number'] + adata.obs['nuclear_number']) / adata.obs['total_number']

# %% tags=[]
plt.hist(np.log10(adata.obs['exonic_number']),20);

# %% tags=[]
plt.hist(adata.obs['exonic_number']/adata.obs['total_number'],40);

# %% tags=[]
plt.hist( adata.obs['spike'] / adata.obs['total_number'],20);

# %% tags=[]
plt.scatter(np.log10(adata.obs['total_number']),exonic_fraction)

# %% tags=[]
adata = adata[adata.obs['exonic_number'] > 10**5,:]
adata = adata[(adata.obs['exonic_number'] / adata.obs['total_number']) > 0.1,:]
adata = adata[adata.obs['percent_mito'] < 0.1,:]
adata = adata[adata.obs['nuclear_number'] > 10**5.45,:]
adata = adata[adata.obs['spike'] / adata.obs['total_number'] < 0.1,:]

# %% tags=[]
adata

# %%
# adata.write(dir_data + 'anndata_filtered_mKO.h5ad')

# %% [markdown]
# Write file for R pipeline

# %%
counts_filt = counts[adata.obs_names]
counts_filt.shape

# %%
a = np.array(adata.var['ids']) 
b = np.array(counts.index[54232:])
c = np.append(a,b)

counts_filt = counts_filt.T[c]
counts_filt.shape

# %% tags=[]
# counts_filt.T.to_csv(dir_data + 'input_mKO_R.csv')

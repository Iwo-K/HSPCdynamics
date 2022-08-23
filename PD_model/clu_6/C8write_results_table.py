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
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt

# %%
dir_adata = '/Users/mb2338/Documents/MATLAB/hoxb5_fate_mapping/october_2021/input_data/'
dir_fates = '/Users/mb2338/Documents/MATLAB/hoxb5_fate_mapping/october_2021/input_data/'

# %%

# %% [markdown]
# read data

# %%
adata = sc.read(dir_adata + 'combined_filt.h5ad')
adata = adata[adata.obs.tom == 'pos', :].copy() #maybe use all cells?
del adata.obsm['X_diffmap'] #For some reason presence of X_diffmap breaks things

# %%

# %%
fates = pd.read_csv(dir_fates + 'pseudotime_kernel_fates.csv',index_col=0)

# %%
adata.obs['DCfate'] = fates['6'].values
adata.obs['Neufate'] = fates['10'].values
adata.obs['NeuDCratio'] = adata.obs['Neufate']/adata.obs['DCfate']
adata.obs['DCNeuratio'] = adata.obs['DCfate']/adata.obs['Neufate']

# %%
filt = ((adata.obs.DCfate > 0.18) & (adata.obs.Neufate < 0.49) & ~adata.obs.leiden.isin(['12', '25', '16']))
adata.obs['filt'] = filt.astype(str)

# %%
adata_rid = adata[adata.obs.filt == 'True',:].copy()

# %%
vec_keep = (adata_rid.obs.start_age == 'young') & (adata_rid.obs.tom == 'pos')

adata_rid = adata_rid[vec_keep,:].copy()

# %%
adata_rid

# %%
vec_rem= ((adata_rid.obs.leiden=='8') | (adata_rid.obs.leiden=='14')|
                      (adata_rid.obs.leiden=='24') )


adata_rid = adata_rid[~vec_rem,:].copy()

# %%
vec_keep = adata_rid.obs.dpt_pseudotime<=0.11

adata_rid = adata_rid[vec_keep,:].copy()

# %%
adata_rid

# %%

# %%
df = pd.DataFrame(adata_rid.obs.dpt_pseudotime)

# %%
parameters1 = pd.read_csv('./tables/parameters_clu_6_dpt_diff_drift.txt',
                          index_col=None,header=None)

# %%
parameters1.columns = ['dpt','diffusion','drift']

# %%
parameters1.head()

# %%
plt.scatter(parameters1.dpt,parameters1.drift)

# %%
m = np.min(adata_rid.obs['dpt_pseudotime'])
M = np.max(adata_rid.obs['dpt_pseudotime'])

mn = np.min(parameters1.dpt)
Mn = np.max(parameters1.dpt)

# %%
dpt = np.array(adata_rid.obs['dpt_pseudotime'])

# %%
dpt_s = (Mn - mn)/(M - m) * (dpt-m) + mn

# %%
df['dpt_pseudotime_scaled'] = dpt_s

# %%
vec_diffusion = np.array([-1.1 for i in range(df.shape[0])],dtype = 'float')
vec_drift = np.array([-1.1 for i in range(df.shape[0])],dtype = 'float')


for i in range(df.shape[0]):
    am = np.argmin(np.abs(df.dpt_pseudotime_scaled[i] - parameters1['dpt']))
    vec_diffusion[i] = parameters1['diffusion'][am]
    vec_drift[i] = parameters1['drift'][am]

# %%
df['diffusion'] = vec_diffusion
df['drift'] = vec_drift

# %%
df.head()

# %%
parameters2 = pd.read_csv('./tables/parameters_clu_6_dpt_growth.txt',
                          index_col=None,header=None)

# %%
parameters2.columns = ['dpt','growth']

# %%
parameters2.head()

# %%
mn = np.min(parameters2.dpt)
Mn = np.max(parameters2.dpt)

# %%
dpt_s = (Mn - mn)/(M - m) * (dpt-m) + mn

# %%
df['dpt_pseudotime_scaled2'] = dpt_s

# %%
vec_growth = np.array([-1.1 for i in range(df.shape[0])],dtype = 'float')


for i in range(df.shape[0]):
    am = np.argmin(np.abs(df.dpt_pseudotime_scaled2[i] - parameters2['dpt']))
    vec_growth[i] = parameters2['growth'][am]

# %%
df['growth'] = vec_growth

# %%
df.head()

# %%
df.to_csv('./tables/table_all_parameters_clu_6.csv')

# %%
df

# %%
plt.scatter(df.dpt_pseudotime_scaled,df.drift)

# %%

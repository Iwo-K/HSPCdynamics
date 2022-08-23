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
sc.settings.vector_friendly = False
sc.set_figure_params( dpi=300, dpi_save = 300, frameon=False, figsize = (7,7), format='jpg',fontsize=30)

# %%
dir_flow = '/Users/mb2338/Documents/MATLAB/hoxb5_fate_mapping/december_2021/input_data/'
dir_adata = '/Users/mb2338/Documents/MATLAB/hoxb5_fate_mapping/october_2021/input_data/'
dir_fates = '/Users/mb2338/Documents/MATLAB/hoxb5_fate_mapping/october_2021/input_data/'

# %%

# %% [markdown]
# read data

# %%
adata = sc.read(dir_adata + 'combined_filt.h5ad')
adata = adata[(adata.obs.tom == 'pos'), :].copy() #maybe use all cells?
del adata.obsm['X_diffmap'] #For some reason presence of X_diffmap breaks things

# %%
adata

# %%
np.unique(adata.obs.sample_id)

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
sc.pl.scatter(adata,basis='umap_2d',color=['filt'],
              size = 30,palette=['k','y'],save='_select_cells_clu_6',title = 'DC')

# %%
adata_rid = adata[adata.obs.filt == 'True',:].copy()

# %%
vec_keep = (adata_rid.obs.start_age == 'young') & (adata_rid.obs.tom == 'pos')

adata_rid = adata_rid[vec_keep,:].copy()

# %%
adata_rid

# %%

# %% [markdown]
# prepare input dpt

# %%
df = pd.DataFrame(list(adata_rid.obs.timepoint_tx_days),columns=['stage'],index = adata_rid.obs_names)

# %%
df['sample'] = adata_rid.obs.biosample_id

# %%
df['dpt_pseudotime'] = adata_rid.obs.dpt_pseudotime

# %%
vec_rem= ((adata_rid.obs.leiden=='8') | (adata_rid.obs.leiden=='14')|
                      (adata_rid.obs.leiden=='24') )

df = df.loc[~vec_rem].copy()

adata_rid = adata_rid[~vec_rem,:].copy()

# %%
plt.hist(adata_rid.obs.dpt_pseudotime,bins=50);
plt.vlines(0.11,0,1000)

# %%
adata_rid

# %%
vec_keep = adata_rid.obs.dpt_pseudotime<=0.11

df = df.loc[vec_keep].copy()

df['stage'] = df['stage']-3
df['dpt_pseudotime'] = df['dpt_pseudotime'] + 0.000001

adata_rid = adata_rid[vec_keep,:].copy()

# %%
adata_rid

# %%
plt.hist(adata_rid.obs.dpt_pseudotime,bins=50);
plt.vlines(0.11,0,1000)

# %%
sc.pl.scatter(adata_rid,basis = 'umap_2d',color = ['leiden'],
              legend_loc='on data',size=  30,title = '',save='_selected_cells_clu_6')

# %%
sc.pl.scatter(adata_rid,basis = 'umap_2d',color = ['dpt_pseudotime'],
              legend_loc='on data',save='_dpt_clu_6',size = 30)

# %%
adata_rid

# %%

# %%

# %%
m = np.min(adata_rid.obs['dpt_pseudotime'])
M = np.max(adata_rid.obs['dpt_pseudotime'])

dpt = np.array(adata_rid.obs['dpt_pseudotime'])

mn = 0
Mn = 1

dpt_s = (Mn - mn)/(M - m) * (dpt-m) + mn
adata_rid.obs['dpt_pseudotime_scaled'] = dpt_s

# %%

# %%
adata_rid.obs['dpt_pseudotime_scaled'] = adata_rid.obs['dpt_pseudotime'] / np.max(adata_rid.obs['dpt_pseudotime'])

# %%
sc.pl.scatter(adata_rid,basis = 'umap_2d',color = ['dpt_pseudotime_scaled'],
              legend_loc='on data',size=  30,save='_scaled_dpt_clu_6')

# %%
adata_rid

# %%
df.shape

# %%
df.to_csv('./tables/input_pseudo_dyn_clu_6_dpt.csv',index=None)

# %%

# %% [markdown]
# prepare input abs num

# %%
df = pd.DataFrame(list(np.unique(adata_rid.obs.timepoint_tx_days)),columns=['day'],index = None)

# %%
df = pd.DataFrame(list(np.unique(adata_rid.obs.timepoint_tx_days)),columns=['day'],index = None)

# %%
df1 = pd.read_csv(dir_flow + 'model_input_tompos.csv',index_col=0)
df1.head()

# %%
df1 = df1.dropna()

# %%
df1.shape

# %%
df2 = pd.crosstab(adata_rid.obs.biosample_id,adata_rid.obs.leiden)

# %%
df2.head()

# %%
df2 = df2.loc[df1.index]

# %%
df2.columns = df2.columns.add_categories(['traj_DC'])

# %%
df2['traj_DC'] = list(np.sum(df2,axis=1)/df1['sc_ncells_total']*df1['flow_total'])

# %%
df2.head()

# %%
mean = np.array([-11 for i in range(len(np.unique(df1.time)))],dtype = 'float') 
std = np.array([-11 for i in range(len(np.unique(df1.time)))],dtype = 'float') 
k = -1
for j in np.unique(df1.time):
    k+=1
 
    mean[k] = np.mean(df2.loc[df1.time == j].traj_DC)
    std[k] = np.std(df2.loc[df1.time == j].traj_DC)/np.sqrt(np.sum(df1.time == j))

# %%
plt.errorbar(range(len(np.unique(df1.time))),mean,std)

# %%
df = df -3

# %%
df['mean'] = mean

# %%
df['sd'] = std

# %%
df['replicates'] = list(df1.groupby('time').count().sc_ncells_total)

# %%
df

# %%
df.to_csv('./tables/input_pseudo_dyn_clu_6_size.csv',index=None)

# %%
adata_rid

# %%

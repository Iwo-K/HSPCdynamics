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
# # Figure PD

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
from scipy.sparse import csr_matrix
import utils.proc, utils.plots
import utils.paga_curved as pagac
from utils.pretty_arrows import plot_arrows, plot_widthbar
from mpl_toolkits.axes_grid1 import make_axes_locatable
import plotnine as pn
from scipy.stats import gaussian_kde

sc.settings.verbosity = 3
base_figures = './figures/26script/'
base_procdata = './procdata/26script/'
sc.settings.figdir  = base_figures
for i in [sc.settings.figdir, base_figures, base_procdata]:
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

def cm2inch(x):
    return x/2.54


# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown]
# ## TITLE

# %%
hoxb5 = sc.read('procdata/04script/combined_filt.h5ad')
meta = pd.read_csv('procdata/04script/combined_sample.meta.csv', index_col = 0)
meta = meta.sort_values(by ='timepoint_tx_days', ascending = False)
hoxb5.obs['longname'] = [meta.loc[i, 'longname'] for i in hoxb5.obs.biosample_id]

# %%
hoxb5.uns['paga'] = hoxb5.uns['pagaPCA']
sc.pl.paga_compare(hoxb5, threshold=0.05, save='_clusters.pdf')

# %% [markdown]
# ## MK

# %%
mkpar = pd.read_csv('./PD_model/clu_7/tables/table_all_parameters_clu_7.csv', index_col=0)
mkpar = mkpar.sort_values(by='dpt_pseudotime')

mk = hoxb5[mkpar.index,:].copy()
mk.obs['drift'] = mkpar['drift']
mk.obs['growth'] = mkpar['growth']
sc.pl.umap(mk, color=['leiden', 'dpt_pseudotime', 'drift', 'growth', 'Pf4', 'Vwf', 'Apoe'], save='_mk_PD.pdf')

# %%
sc.set_figure_params(dpi=150, dpi_save = 400, frameon=False, figsize = (cm2inch(6),cm2inch(6)), format='png',fontsize=6.5)


fig, axes = plt.subplots(nrows = 2)
axes[0].plot(mkpar.dpt_pseudotime, mkpar.drift, lw=2)
axes[1].plot(mkpar.dpt_pseudotime, mkpar.growth, lw=2)
plt.xlabel("Pseudotime")
axes[0].set_ylabel("Drift", fontsize=7.5)
axes[1].set_ylabel("Growth", fontsize=7.5)
fig.tight_layout()
fig.savefig(base_figures + 'mk_drift_growth.pdf')
fig.show()

# %%
sc.pl.umap(mk, color=['S_score', 'G2M_score', 'phase'], save='_mk_cc_scores.pdf')

# %% [markdown]
# Plotting cell densities over pseudotime

# %% [markdown]
# ### Calculating mean density along pseudotime

# %%
# Loading cells with assigned real time and pseudotime
mkcells = pd.read_csv('./PD_model/clu_7/tables/input_pseudo_dyn_clu_7_dpt.csv')
mkcells['stage'] = mkcells.stage + 3 #Time was offset to 0
# mkcells = mkcells.loc[mkcells.stage <=49,:]
points = np.linspace(0, mkcells.dpt_pseudotime.max(), num=100)

ds = pd.DataFrame(columns=['d', 'points', 'sample', 'stage'])
for i in mkcells['sample'].unique():
    x = gaussian_kde(mkcells.loc[mkcells['sample'] == i, 'dpt_pseudotime'], bw_method=0.4)(points)
    df = pd.DataFrame(dict(d=x, points=points, sample=i, stage=mkcells.loc[mkcells['sample'] == i, 'stage'].iloc[0]))
    ds = pd.concat((ds, df))
    
dsg = ds.groupby(['stage', 'points']).mean()
dsg = dsg.reset_index()
dsg['stage'] = pd.Categorical(dsg.stage)

dsgsub = dsg.loc[dsg.stage.isin([3,7,12,27,49]),:]
dsgsub.stage = dsgsub.stage.cat.remove_unused_categories()
                                
cols = {3 : '#0069B4', 7 : '#009FE3', 12 : '#8FC89A', 27 : '#EAD303', 49 : '#F39200'} # , 76 : 'grey', 112 : 'grey', 169 : 'grey'}
g1 = (pn.ggplot(dsgsub, pn.aes(x='points', y='d', group='stage', color='stage')) + pn.geom_line(size=2, alpha=0.5) +
 pn.scale_color_manual(values=cols) +
pn.theme_bw())
g1.draw()
             
dsgsub = dsg.loc[dsg.stage.isin([3,7,12,27]),:]
dsgsub.stage = dsgsub.stage.cat.remove_unused_categories()

dsgbak = dsg.loc[dsg.stage.isin([49, 112, 161, 269]),:]
dsgbak = dsgbak.groupby(['points']).mean()
dsgbak = dsgbak.reset_index()
dsgbak['stage'] = 49 #column needed otherwise plotnine complains

g3 = (pn.ggplot(dsgsub, pn.aes(x='points', y='d', group='stage', color='stage')) + 
      pn.geom_area(data=dsgbak, mapping=pn.aes(x='points', y='d'), color=None, size=2, alpha=0.15) +
      pn.geom_line(size=1, alpha=0.7) + 
      pn.scale_color_manual(values=cols) +
      pn.theme_classic(base_size=7.5))
g3.save(base_figures + 'mk_densities_pseudotime.pdf', width=72/25.4, height=38/25.4)
g3

# %%
mksizes = pd.read_csv('./PD_model/clu_7/tables/input_pseudo_dyn_clu_7_size.csv')
sizes = []
for i in dsg.stage:
    sizes.append(mksizes.loc[mksizes.day == (i-3), 'mean'].iloc[0])
dsg['ncells'] = sizes
dsg['dsize'] = dsg.d * dsg.ncells
dsgsub = dsg.loc[dsg.stage.isin([3,7,12,27]),:]
dsgsub.stage = dsgsub.stage.cat.remove_unused_categories()

g2 = (pn.ggplot(dsgsub, pn.aes(x='points', y='dsize', group='stage', color='stage')) + pn.geom_line(size=2, alpha=0.5) +
      pn.scale_color_manual(values=cols) +
      pn.theme_bw(base_size=7.5) + 
      pn.scale_y_log10())
g2.save(base_figures + 'mk_dsizes_pseudotime_log.pdf', width=50/25.4, height=34/25.4)
g2.draw()

mkcells = pd.read_csv('./PD_model/clu_7/tables/input_pseudo_dyn_clu_7_dpt.csv')
mkcells['stage'] = mkcells.stage + 3
mkcells = mkcells.loc[mkcells.stage <=49,:]
mkcells['stage'] = pd.Categorical(mkcells.stage)
(pn.ggplot(mkcells, pn.aes(x="dpt_pseudotime", group="sample", color="stage")) +
 pn.geom_density() + 
 pn.theme_bw() + 
 pn.scales.scale_color_discrete() + 
pn.facet_wrap("~stage"))

# %% [markdown]
# ## Ery

# %%
erypar = pd.read_csv('./PD_model/clu_11/tables/table_all_parameters_clu_11.csv', index_col=0)
ery = hoxb5[erypar.index,:].copy()
ery.obs['drift'] = erypar['drift']
ery.obs['growth'] = erypar['growth']
sc.pl.umap(ery, color=['leiden', 'dpt_pseudotime', 'drift', 'growth', 'Klf1', 'Gata1', 'Hba-a1'], save='_ery_PD.pdf')

# %% [markdown]
# ## DC

# %%
DCpar = pd.read_csv('./PD_model/clu_6/tables/table_all_parameters_clu_6.csv', index_col=0)
DC = hoxb5[DCpar.index,:].copy()
DC.obs['drift'] = DCpar['drift']
DC.obs['growth'] = DCpar['growth']
sc.pl.umap(DC, color=['leiden', 'dpt_pseudotime', 'drift', 'growth', 'Csf1r', 'Mpo'], save='_DC_PD.pdf')

# %% [markdown]
# ## Neu

# %%
neupar = pd.read_csv('./PD_model/clu_10/tables/table_all_parameters_clu_10.csv', index_col=0)
neu = hoxb5[neupar.index,:].copy()
neu.obs['drift'] = neupar['drift']
neu.obs['growth'] = neupar['growth']
sc.pl.umap(neu, color=['leiden', 'dpt_pseudotime', 'drift', 'growth', 'Elane', 'Mpo'], save='_neu_PD.pdf')

# %% tags=[]
#Flt3 coexpressed genes
g = 'Cpa3','Srgn','Ms4a3','Rnase6','Prss57','Lat2','Rab37','Irf8','Chst13','Cpxm1','Prom1','Cd34','Lair1','Nucb2'
sc.pl.umap(neu, color=g)
sc.pl.umap(neu, color='Flt3')

# %% tags=[]
#Il2/Stat5 pathway
sc.pl.umap(neu, color=['Selp', 'Ifngr1', 'Irf8', 'Bcl2', 'Ikzf2', 'Cst7', 'Ndrg1', 'Dhrs3'])
#EMT genes
sc.pl.umap(neu, color=['Vcam1', 'Serpine2', 'Tpm4', 'Basp1', 'Itga2', 'Bgn', 'Tgfbi', 'Sat1'])
#KRAS signaling up
sc.pl.umap(neu, color=['Hsd11b1', 'Lat2', 'Gprc5b', 'Itga2', 'Epb41l3', 'Irf8', 'G0s2'])

# %%
from matplotlib.colors import LinearSegmentedColormap
cmap2 = LinearSegmentedColormap.from_list('mycmap', [(0, '#e8e8e8'),
                                                     (1, '#cc0404')])
sc.pl.umap(hoxb5, color=['Gfi1', 'Flt3', 'Irf8'], cmap=cmap2, save='_Gfi1_Flt3_Irf8.pdf')

# %% [markdown]
# ## Ly

# %%
lypar = pd.read_csv('./PD_model/clu_16/tables/table_all_parameters_clu_16.csv', index_col=0)
ly = hoxb5[lypar.index,:].copy()
ly.obs['drift'] = lypar['drift']
ly.obs['growth'] = lypar['growth']
sc.pl.umap(ly, color=['leiden', 'dpt_pseudotime', 'drift', 'growth', 'Dntt', 'Il7r'], save='_ly_PD.pdf')


# %% [markdown]
# ## Waiting times

# %%
def get_colormap(adata, obs_column):
    if pd.api.types.is_categorical_dtype(adata.obs[obs_column]):
        cats = adata.obs[obs_column].cat.categories.values
        cols = adata.uns[obs_column + '_colors']
        if len(cats) != len(cols):
            raise Exception('The number of categories and colors do not match, are you sure they are valid?')
        m = {cats[i] : cols[i] for i in range(len(cols))}
        return m
    else:
        raise Exception('This column is not categorical')
        
from scanpy.plotting.palettes import default_20
hoxb5.uns['anno_man_colors'] = default_20[0 : len(hoxb5.obs.anno_man.cat.categories)]
x = np.where(hoxb5.obs.anno_man.cat.categories == 'Int prog')[0][0]
hoxb5.uns['anno_man_colors'][x] = "#C5C5C5"

leiden_colmap = get_colormap(hoxb5, 'leiden')
anno_colmap = get_colormap(hoxb5, 'anno_man')
leiden_anno_map = hoxb5.obs[['anno_man', 'leiden']]
leiden_anno_map = leiden_anno_map.drop_duplicates()
leiden_anno_map = dict(zip(leiden_anno_map.leiden, leiden_anno_map.anno_man))

# %%
for clu in [6, 10, 11, 16, 7]:
    table = pd.read_csv(f"PD_model/clu_{str(clu)}/tables/table_all_parameters_clu_{str(clu)}.csv",index_col=0)
    wt = pd.read_csv(f"PD_model/clu_{str(clu)}/tables/waiting_times_clu_{str(clu)}.txt",header = None)
    wt[1] = wt[1]-1
    table = table.sort_values(by = 'dpt_pseudotime')
    wt.index = table.index
    table['waiting_time'] = wt[1]
    hoxb5.obs[f"wt{str(clu)}_log"] = np.log10(wt[1]+100)
    hoxb5.obs[f"wt{str(clu)}"] = wt[1]

# %% [markdown]
# Plotting thresholded trajectories

# %%
toplot = []
for clu in [6, 10, 11, 16, 7]:
    colnam = f'Trajectory {str(clu)}'
    hoxb5.obs[colnam] = hoxb5.obs[f"wt{str(clu)}"]
    hoxb5.obs.loc[~hoxb5.obs[colnam].isna(),[colnam]] = 'YES'
    toplot.append(colnam)
sc.pl.umap(hoxb5, color=toplot, na_color='lightgrey', legend_loc=None, save='_thresholded_trajectories.pdf')


# %% [markdown]
# Plotting waiting times

# %% tags=[]
def plot_cellfront(adata, cluster, color, tpoints, save=None):
    fig, axes = plt.subplots(ncols=len(tpoints))
    fig.set_size_inches(1.45*len(tpoints), 1.7)
    for i,ax in zip(tpoints, axes):
        sc.pl.umap(adata, ax=ax, show=False, alpha=0.5)
        sc.pl.umap(adata[adata.obs['wt' + cluster] == i,:],
                   na_color=color,
                   ax=ax,
                   title=f"Day {str(i)}",
                   size=16,
                   show=False)
        fig.subplots_adjust(top=0.80)
        plt.suptitle(f'Cluster {cluster}', fontproperties = {'size' : 10}, color=color)
        
    plt.tight_layout()
    if save is not None:
        plt.savefig(sc.settings.figdir / (save))
    fig.show()


# %%
sc.set_figure_params(fontsize=7, frameon=False)
col = anno_colmap[leiden_anno_map['6']]
plot_cellfront(hoxb5, cluster='6', color=col, tpoints=[0, 4, 16, 24, 32], save='cl6_wt.pdf')

col = anno_colmap[leiden_anno_map['7']]
# plot_cellfront(hoxb5, cluster='7', color=col, tpoints=[0, 1, 2, 3, 4, 5, 6, 8, 9])
plot_cellfront(hoxb5, cluster='7', color=col, tpoints=[0, 1, 3, 5, 10], save='cl7_wt.pdf')

col = anno_colmap[leiden_anno_map['10']]
plot_cellfront(hoxb5, cluster='10', color=col, tpoints=[0, 4, 24, 30, 42], save='cl10_wt.pdf')

col = anno_colmap[leiden_anno_map['11']]
plot_cellfront(hoxb5, cluster='11', color=col, tpoints=[0, 4, 8, 12, 18], save='cl11_wt.pdf')

col = anno_colmap[leiden_anno_map['16']]
plot_cellfront(hoxb5, cluster='16', color=col, tpoints=[0, 4, 12, 30, 47], save='cl16_wt.pdf')

# for i in ['6', '7', '10', '11', '16']:
#     col = anno_colmap[leiden_anno_map[i]]
#     plot_cellfront(hoxb5, cluster=i, color=col, tpoints=[0, 4, 24, 39, 60], save=f'cl{i}_wt.pdf')

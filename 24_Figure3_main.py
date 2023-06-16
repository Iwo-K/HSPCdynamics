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
# # Plots for Figure 3

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
from sklearn.metrics import pairwise_distances
import networkx as nx
from scipy.stats import linregress

sc.settings.verbosity = 3
sc.settings.figdir = './figures/24script/'
base_figures = './figures/24script/'
base_procdata = './procdata/24script/'
for i in [sc.settings.figdir, base_figures, base_procdata]:
    pathlib.Path(i).mkdir(parents=True, exist_ok=True)

# %%
#Figure settings
sc.set_figure_params(color_map = 'viridis', dpi_save = 350)
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
# ## Loading data

# %%
hoxb5 = sc.read('procdata/04script/combined_filt.h5ad')
meta = pd.read_csv('procdata/04script/combined_sample.meta.csv', index_col = 0)
meta = meta.sort_values(by ='timepoint_tx_days', ascending = False)
hoxb5.obs['longname'] = [meta.loc[i, 'longname'] for i in hoxb5.obs.biosample_id]

# %%
hoxb5.uns['paga'] = hoxb5.uns['pagaPCA']
sc.pl.paga_compare(hoxb5, threshold=0.05, save='_clusters.pdf')

# %% [markdown]
# ## Pltting clusters and PAGA with annotation

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


# %%
#Swapping Int prog to grey color
from scanpy.plotting.palettes import default_20
hoxb5.uns['anno_man_colors'] = default_20[0 : len(hoxb5.obs.anno_man.cat.categories)]
x = np.where(hoxb5.obs.anno_man.cat.categories == 'Int prog')[0][0]
hoxb5.uns['anno_man_colors'][x] = "#C5C5C5"
clustmap = get_colormap(hoxb5, 'leiden')
annomap = get_colormap(hoxb5, 'anno_man')

#Map from clusters to annotation
m = hoxb5.obs[['leiden', 'anno_man']]
m = m.loc[~m.duplicated(),:]
m.index = m.leiden
m = m['anno_man'].to_dict()

#Substituting the leiden colors for anno_man colors (to avoid using the color argument in paga which plots piecharts)
hoxb5.uns['leiden_colors_bak'] = hoxb5.uns['leiden_colors']
hoxb5.uns['leiden_colors'] = [annomap[m[i]] for i in hoxb5.obs.leiden.cat.categories]

# %%
sc.set_figure_params( dpi=150, dpi_save = 400, frameon=False, figsize = (5,4), format='png',fontsize=14)

sc.pl.paga(hoxb5,
           pos = hoxb5.uns['pagaPCA']['pos'],
           text_kwds={'color' : 'black'},
           node_size_scale = 2,
           save='_anno_man.pdf')

# %% [markdown]
# ## Plotting cell numbers on PAGA over time

# %%
# Getting cell numbers
cellno = pd.read_csv('./procdata/05script/model_input_tompos.csv', index_col=0)
cellno['cellno'] = cellno['sc_ncells_filt'] / cellno['sc_ncells_total'] * cellno['flow_total']
cellno = cellno.loc[~cellno.cellno.isna(),:]

cellnoN = cellno.loc[:, cellno.columns.isin(hoxb5.obs.leiden.cat.categories)]
cellnoN = (cellnoN.T / cellnoN.sum(axis=1)).T
cellnoN['time'] = cellno.time
cellnoN = cellnoN.groupby('time').mean()

avtotal = cellno[['cellno', 'time']]
avtotal = avtotal.groupby('time').mean()
cellnoN = (cellnoN.T * avtotal['cellno']).T

# %%
# Constructing graph
x = hoxb5.uns['paga']['connectivities'].todense()
x[x <= 0.01 ] = 0
g = nx.Graph(x)
cats = hoxb5.obs.leiden.cat.categories

# Plotting
ncols = len(cellnoN.index)
fig, axes = plt.subplots(nrows=1, ncols=ncols)
fig.set_size_inches(ncols*4, 4, forward=True)
divider = make_axes_locatable(axes[-1])
cax = divider.append_axes('right', size='5%', pad=0.05)
cax = [cax]

for i,ax in zip(cellnoN.index,axes):
    nx.draw(g,
            pos=hoxb5.uns['pagaPCA']['pos'],
            node_size=np.sqrt(hoxb5.uns['leiden_sizes'])*2,
            node_color=np.log10(cellnoN.loc[i,:]+1),
            vmin=0,
            vmax=3.8,
            with_labels=False,
            cmap='Reds',
            edge_color='lightgrey',
            width=0.95,
            style='dashed',
           ax = ax)
    labels = {n:cats[n] for n in range(len(cats))}
    nx.draw_networkx_labels(g, pos=hoxb5.uns['pagaPCA']['pos'], labels=labels)
    cmap=plt.cm.Reds
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin = 0, vmax=3.8))
    ax.set_title(str(i))
    # sm = plt.cm.ScalarMappable(cmap=cmap, norm=mpl.colors.LogNorm(vmin = 0, vmax=3.8))
    
    sm._A = []
    plt.colorbar(sm, cax=cax[0])

fig.show()
fig.tight_layout()
fig.savefig(base_figures + 'paga_timecourse.pdf')

# %% [markdown]
# ## Plotting Discrete model

# %% [markdown]
# ### Getting data from discrete model - Self-renewal rates (k), net proliferation, cluster sizes etc

# %%
k = pd.read_csv('./discrete_model/output/self_renewal.txt', header=None)[0].values
#Cluster sizes, using initial conditions
leiden_sizes =pd.read_csv('./discrete_model/output/best.txt', header=None).T[0:22][0].values
leiden_sizes = leiden_sizes / leiden_sizes[0]
net_prolif = pd.read_csv('./discrete_model/output/net_proliferation.txt', header=None)[0].values
df = pd.DataFrame({'leiden' : np.concatenate((hoxb5.obs.leiden.cat.categories.values, ['30'], ['40'])),
                   'leiden_sizes' : leiden_sizes,
                  'net_prolif' : net_prolif,
                  'k' : k})

# Loading matlab output
df.leiden = df.leiden.astype(str)
#Transforming k rates with a pseudovalue of 0.01
df['logk'] = np.log(df.k + 0.01)
df['SR'] = 1/df.k
df.loc[df.SR < 0,'SR'] = np.inf
#Switching the infinity to 1.2* maximum value
df.loc[df.SR == np.inf,'SR'] = df.SR[df.SR != np.inf].max()*1.2
df['log10SR'] = np.log10(df.SR)
df.loc[df.SR < 0,'SR'] = np.inf

# Threshold for considering self renewal
print(f'Threshold for considering self-renewal on log scale: {np.log(0.01 + 0.01)}')

# %%
# Hacking AnnData/PAGA to think there are two additional clusters
hoxb5.obs['leiden_orig'] = hoxb5.obs.leiden.copy()
vec_leiden = list(hoxb5.obs.leiden)
vec_leiden[0] = '30'
vec_leiden[1] = '40'
hoxb5.obs.leiden = pd.Categorical(vec_leiden, categories = df.leiden)

# Adding logk rates to each cell
vec = np.zeros(hoxb5.n_obs, dtype = 'float')
for clu in df.leiden:
    vec[hoxb5.obs.leiden == clu] = df.loc[df.leiden == clu, 'logk'] 
hoxb5.obs['logk'] = vec

# Adding net proliferation rates to each cell
vec = np.zeros(hoxb5.n_obs, dtype = 'float')
for clu in df.leiden:
    vec[hoxb5.obs.leiden == clu] = df.loc[df.leiden == clu, 'net_prolif'] 
hoxb5.obs['net_prolif'] = vec

# Adding self-renewal/residence time to each cell
vec = np.zeros(hoxb5.n_obs, dtype = 'float')
for clu in df.leiden:
    vec[hoxb5.obs.leiden == clu] = df.loc[df.leiden == clu, 'log10SR'] 
hoxb5.obs['log10SR'] = vec

# Adding cluster sies
hoxb5.uns['leiden_sizes'] = np.ceil(df.leiden_sizes.values*10000)

# %% [markdown]
# ### Plotting differentiation rates on PAGA-like graphs
# The plots are composed of several layers (with different types of arrows)

# %% [markdown]
# #### Preparing for plotting

# %%
# Loading matrix with differentiation rates (best fit)
drates = pd.read_csv('discrete_model/output/differentiation_matrix.txt',
                delimiter='\t', header=None)
drates.columns = df.leiden
drates.index = df.leiden
#transposing as paga/networkx read columns to rows
drates = drates.T

# Adding matrix with differentiation rates to AnnData 
logdrates = np.log1p(drates.values)
logdrates = csr_matrix(logdrates)
hoxb5.uns['paga']['transitions'] = logdrates

# %%
# Adding positions for cluster 30 and 40
pos = hoxb5.uns['pagaPCA']['pos'].copy()
pos2 = np.concatenate((pos, np.array([[8.25,15]]), np.array([[5.25,15]])), axis=0)
hoxb5.uns['pagaPCA']['pos'] = pos2

#Renaming clusters 0, 30 and 40 to 0c, 0b, 0a
hoxb5.obs['leiden'] = hoxb5.obs.leiden.cat.rename_categories({'0' : '0c', '30' : '0a', '40' : '0b'})

sc.set_figure_params( dpi=300, dpi_save = 400, frameon=False, figsize = (6,4), fontsize=12)

# %% [markdown]
# #### net proliferation and differentiation rates plot

# %%
#Plotting differentiation rates (arrows) + net proliferation rates (colors)
#All rates below 0.01 plotted as dashed lines
fig, ax = plt.subplots(1,1, figsize = (6.5,5))
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
cax = [cax]

#Plotting parameters
nodesizes = np.power(leiden_sizes, 0.5)  # This is to match node sizes in paga
trm = logdrates.toarray()
trm1 = trm.copy()
trmin = 0.01

#Plotting arrows above the threshold
trm1[trm1 < trmin] = 0
plot_arrows(hoxb5,
            transitions=trm1,
            connectionstyle='arc3, rad=0.2',
            edge_width_scale=5,
            node_size=nodesizes*50,
            arrowsize=7,
            ax=ax)

#Plotting arrows below the threshold (dashed)
trm2 = trm.copy()
trm2[trm2 >= trmin] = 0
plot_arrows(hoxb5,
            transitions=trm2,
            connectionstyle='arc3, rad=0.2',
            width= 0.3,
            node_size=nodesizes*40,
            arrowsize=7,
            style=(0, (1,6)),
            edge_color='grey',
            ax=ax)

pagac.paga(hoxb5,
           transitions='transitions',
           pos = hoxb5.uns['pagaPCA']['pos'],
           text_kwds={'color' : '#00cc52'}, # #00cc52 or #90bfa3 or #05d65e or springgreen - #00FF7F
           threshold = 10^16, 
           color = 'net_prolif',
           node_size_scale = 1.8, cmap='plasma',
           ax=ax,
           cax=cax,
           show=False)
plt.savefig(base_figures + 'paga_netprolif_drates.pdf')
plt.show()

# %% [markdown]
# #### Plotting residence time and cell flux

# %%
# Getting the differentiation rates (row to columns direction)
xflux = drates.values
# Adjusting for source cluster sie
xflux = np.einsum('ij,j->ij', xflux, hoxb5.uns['leiden_sizes'])
np.sort(xflux[xflux.nonzero()])

# %% tags=[]
# Figure setup
fig, ax = plt.subplots(1,1, figsize = (6.5,5))
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
cax = [cax]
cax2 = divider.append_axes('bottom', size='20%', pad=0.05)

# Plotting parameters
trm = xflux #  transitions matrix
nodesizes = np.power(leiden_sizes, 0.5)  # This is to match node sizes in paga
trmin = 500  # transitions below thi value will be dashed
trmax=trm.max() # transitions above this value will be clipped to this value
trmax=7.5e+05 # transitions above this value will be clipped to this value
# edge_width_scale = 0.000009
edge_width_scale = 0.000013

#Plotting only thick arrows (values are clipped at trmax)
trm0 = trm.copy()
trm0[trm0 >= trmax] = trmax
trm0[trm0 < trmax] = 0

plot_arrows(hoxb5,
            transitions=trm0,
            connectionstyle='arc3, rad=0.2',
            edge_width_scale=edge_width_scale,
            node_size=nodesizes*30, # Slightly different to prettify the figure
            arrowsize=6,
            ax=ax)

#Plotting mid-sized arrows
trm1 = trm.copy()
trm1[np.logical_or(trm1 < trmin, trm1 >= trmax)] = 0

plot_arrows(hoxb5,
            transitions=trm1,
            connectionstyle='arc3, rad=0.2',
            edge_width_scale=edge_width_scale,
            node_size=nodesizes*50,
            arrowsize=6,
            ax=ax)

# Plotting dashed arrows (values below the threshold)
trm2 = trm.copy()
trm2[np.logical_and(trm2 >= trmin, trm2 != 0)] = 0
trm2[np.logical_and(trm2 < trmin, trm2 != 0)] = 1/trmin
plot_arrows(hoxb5,
            transitions=trm2,
            connectionstyle='arc3, rad=0.2',
            width= 0.3,
            node_size=nodesizes*50,
            arrowsize=6,
            style=(0, (1,6)),
            edge_color='grey',
            ax=ax)

# Plotting nodes
pagac.paga(hoxb5,
           transitions='transitions',
           pos = hoxb5.uns['pagaPCA']['pos'],
           text_kwds={'color' : '#00cc52'}, # #00cc52 or #90bfa3 or #05d65e or springgreen - #00FF7F
           threshold = 10^6, 
           color = 'log10SR',
           node_size_scale = 1.8,
           cmap='plasma',
           ax=ax,
           cax=cax,
           show=False)

plot_widthbar(vmin=trm1[np.nonzero(trm1)].min(),
              vmax=trmax,
              edge_width_scale=edge_width_scale,
              ndigits=-3,
              ax=cax2)
cax2.margins(0.2)

plt.savefig(base_figures + 'paga_log10SR_flux.pdf')
plt.show()

# %% [markdown] tags=[]
# ## Plotting waiting times

# %%
leiden_anno_map = hoxb5.obs[['anno_man', 'leiden']]
leiden_anno_map = leiden_anno_map.drop_duplicates()
leiden_anno_map = dict(zip(leiden_anno_map.leiden, leiden_anno_map.anno_man))

# %%
wtimesD = pd.read_csv('./discrete_model/output/waiting_times.txt', header=None)
wtimesD.columns = ['cluster', 'wtime']
wtimesD['cluster'] = wtimesD.cluster.astype(int).astype(str)
wtimesD.index = wtimesD.cluster.values

wtimesD['anno'] = [leiden_anno_map[i] for i in wtimesD.cluster]
wtimesD['anno_color'] = [annomap[i] for i in wtimesD.anno]

wtimesD['dest_type'] = 'unassigned'
terminal = ['20', '10', '24', '14', '7', '16']
small = ['25', '26', '28']
for i in wtimesD.index:
    if i in terminal:
        wtimesD.loc[i,'dest_type'] = 'terminal'
    elif i in small:
        wtimesD.loc[i,'dest_type'] = 'small'
    else:
        wtimesD.loc[i,'dest_type'] = 'intermediate'

# %%
g2 = (pn.ggplot(wtimesD, pn.aes('cluster', 'wtime', fill='anno'))
 + pn.scale_fill_manual(values = annomap)
 + pn.geom_bar(stat = 'identity')
 + pn.facet_wrap('~dest_type', scales='free')
 + pn.theme_bw(base_size=7.5)
 + pn.theme(subplots_adjust={'wspace' : 0.25}))
g2.save(base_figures + './wtimes_deterministic.pdf', width = cm2inch(12), height = cm2inch(5))
g2

# %%
hoxb5.obs['wtimeD'] = [wtimesD.wtime[i] if i != '0' else None for i in hoxb5.obs.leiden_orig]
# hoxb5.obs['wtimeD'] = hoxb5.obs.wtime
sc.set_figure_params(figsize = (5,4), frameon=False)
sc.pl.umap(hoxb5, color='wtimeD', norm=mpl.colors.LogNorm(vmin=10), save='_wtimes_deterministic.pdf')

# %% [markdown] jp-MarkdownHeadingCollapsed=true tags=[]
# ## Switch off (cluster 0)

# %%
swoff = pd.read_csv('./discrete_model/output/switchoff_deterministic1  21  22.txt',
                    index_col=0, header=None)
t = np.array(range(1,271))
swoff.columns = t
swoff = swoff.T

# %%
swoff

# %%
swoff.columns = hoxb5.obs.leiden.cat.categories.values
swoff['time'] = swoff.index.values
swoff = swoff[['time', '10', '7', '16', '11']]
swoff = swoff.melt(id_vars = 'time', var_name='cluster', value_name='relative size')
swoff['anno'] = [leiden_anno_map[i] for i in swoff.cluster]
swoff['anno_color'] = [annomap[i] for i in swoff.anno]
swoff

# %%
g2 = (pn.ggplot(swoff, pn.aes(x='time', y='relative size', color='anno'))
+ pn.geom_line(size = 1) 
+ pn.theme_bw(base_size=7.5)
+ pn.scale_color_manual(values=annomap)
+ pn.scale_x_continuous(breaks=[0, 10, 30, 50, 100, 200, 300])
+ pn.xlab('time (days)')
+ pn.ylab('rel. cluster size'))      
g2.save(base_figures + './swoff_0_30_40.pdf', width = cm2inch(6), height = cm2inch(6))
print(g2)

# %% [markdown] tags=[] jp-MarkdownHeadingCollapsed=true tags=[]
# ## Comparison of differentiation rates with connectivities or changes in pseudotime

# %%
sc.set_figure_params( dpi=150, dpi_save = 400, frameon=False, figsize = (4,4), fontsize=10)

# %%
keep1 = np.array(['0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])
mask1 = hoxb5.obs.leiden_orig.cat.categories.isin(keep1)

#Connectivities matrix for relevant clusters
cons = hoxb5.uns['pagaPCA']['connectivities'].toarray()
cons = cons[mask1,:][:,mask1]
#To 1d array
cons = cons.ravel(order='F')

#Differentiation rates matrix for relevant clusters
keep2 = np.array(['0c', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12'])
mask2 = hoxb5.obs.leiden.cat.categories.isin(keep2)
drates2 = drates.loc[mask2, mask2].values
drates2 = drates2.ravel(order='F')

#Change in pseudotime for relevant clusters
dpt = hoxb5.obs[['leiden_orig', 'dpt_pseudotime']]
dpt = dpt.groupby('leiden_orig').dpt_pseudotime.mean()
dpt = dpt[dpt.index.isin(keep1)].values
ddpt = pairwise_distances(dpt.reshape(-1, 1), metric='sqeuclidean')
ddpt = ddpt.ravel(order='F')

#Making indices for annotation
inds = [[j + '_' + i for i in keep1] for j in keep1]
inds = np.array(sum(inds, []))
# np.array(np.meshgrid(keep, keep)).T.reshape(-1, 2)

#aseembling into dataframe
df5 = pd.DataFrame({'cons' : cons, 'drates' : drates2, 'ddpt' : ddpt}, index=inds)

#Looking only at transitions with diff rates greater than 
df5 = df5.loc[df5.drates > 1e-12,:]

# %%
ds = hoxb5.obs[['leiden_orig', 'dpt_pseudotime']]
ds = ds.groupby('leiden_orig').dpt_pseudotime.mean()
ds = ds[ds.index.isin(keep1)].values
from sklearn.metrics import pairwise_distances
dpt = pairwise_distances(ds.reshape(-1, 1), metric='sqeuclidean').ravel(order='F')

# %%
sc.set_figure_params(dpi=150, dpi_save = 400, frameon=False, figsize = (cm2inch(6),cm2inch(6)), format='png',fontsize=7.5)

# %%
#Hack adapted from
import networkx as nx
def repel_labels(ax, x, y, labels, k=0.01, ex_mins=0.05, ex_maxs=0.05):
    G = nx.DiGraph()
    data_nodes = []
    init_pos = {}
    for xi, yi, label in zip(x, y, labels):
        data_str = 'data_{0}'.format(label)
        G.add_node(data_str)
        G.add_node(label)
        G.add_edge(label, data_str)
        data_nodes.append(data_str)
        init_pos[data_str] = (xi, yi)
        init_pos[label] = (xi, yi)

    pos = nx.spring_layout(G, pos=init_pos, fixed=data_nodes, k=k)

    # undo spring_layout's rescaling
    pos_after = np.vstack([pos[d] for d in data_nodes])
    pos_before = np.vstack([init_pos[d] for d in data_nodes])
    scale, shift_x = np.polyfit(pos_after[:,0], pos_before[:,0], 1)
    scale, shift_y = np.polyfit(pos_after[:,1], pos_before[:,1], 1)
    shift = np.array([shift_x, shift_y])
    for key, val in pos.items():
        pos[key] = (val*scale) + shift

    for label, data_str in G.edges():
        ax.annotate(label,
                    xy=pos[data_str], xycoords='data',
                    fontsize=4,
                    xytext=pos[label], textcoords='data',
                    arrowprops=dict(arrowstyle="->",
                                    shrinkA=0, shrinkB=0,
                                    connectionstyle="arc3", 
                                    color='grey'), )
    # expand limits
    all_pos = np.vstack(pos.values())
    x_span, y_span = np.ptp(all_pos, axis=0)
    mins = np.min(all_pos-x_span*ex_mins, 0)
    maxs = np.max(all_pos+y_span*ex_maxs, 0)
    ax.set_xlim([mins[0], maxs[0]])
    ax.set_ylim([mins[1], maxs[1]])

import seaborn as sns
slope, intercept, r, p, se = linregress(x=np.log10(df5.cons.values), y=np.log10(df5.drates.values))
print(r**2)

sns.regplot(np.log10(df5.cons.values), np.log10(df5.drates.values), scatter_kws={'s':16})
plt.xlabel("log10(connectivity)")
plt.ylabel("log10(diff. rate)")
ax = plt.gca()
repel_labels(ax, np.log10(df5.cons).values, np.log10(df5.drates).values, labels=df5.index.values,
             k=0.25, ex_mins=0.07, ex_maxs=0.01)
plt.savefig(base_figures + 'cons_drates.pdf')
plt.show()


# %% tags=[]
slope, intercept, r, p, se = linregress(x=np.log10(df5.ddpt.values), y=np.log10(df5.drates.values))
print(r**2)

# %%
sc.set_figure_params( dpi=150, dpi_save = 400, frameon=False, figsize = (cm2inch(6),cm2inch(6)), format='png',fontsize=7.5)

sns.regplot(x=np.log10(df5.ddpt.values), y=np.log10(df5.drates.values))
plt.xlabel("log10(pseudotime change)")
plt.ylabel("log10(differentiation rates)")
for i in range(len(df5.index)):
    plt.annotate(df5.index[i], (np.log10(df5.ddpt[i]), np.log10(df5.drates[i])),
                 fontsize=5,
                 xytext=(4, 1), textcoords='offset points',)

plt.savefig(base_figures + 'drate_ddpt.pdf')
plt.show()

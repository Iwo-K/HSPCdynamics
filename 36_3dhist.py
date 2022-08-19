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
# # Surface marker analysis

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
from anndata import AnnData
from scipy.sparse import csr_matrix
import networkx as nx

sc.settings.verbosity = 3
base_figures = './figures/36script/'
base_procdata = './procdata/36script/'
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

# %%
# %load_ext autoreload
# %autoreload 2

# %% [markdown] tags=[]
# ## Preparing data

# %%
hoxb5 = sc.read('procdata/04script/combined_filt.h5ad')
hoxb5 = hoxb5[(hoxb5.obs.tom == 'pos') & (hoxb5.obs.start_age == 'young'), :].copy()

# %% [markdown]
# ## Plotting bins

# %%
xc = hoxb5.obsm['X_umap'][:,0]
yc = hoxb5.obsm['X_umap'][:,1]
h = plt.hist2d(xc, yc, bins=40)

#Getting coordinates for later
xmids = np.array([(h[1][n] + h[1][n+1])/2 for n in range(len(h[1])-1)])
ymids = np.array([(h[2][n] + h[2][n+1])/2 for n in range(len(h[2])-1)])
X, Y = np.meshgrid(xmids, ymids, indexing='ij')
Z = h[0]

# %%
fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
surf = ax.plot_surface(X, Y, h[0], cmap='coolwarm',
                       linewidth=0, antialiased=False)
ax.view_init(60, -80)
fig.show()

#Smoothing with a gaussian filter (optional)
import plotly.graph_objects as go
import scipy.ndimage
sigma = [1.7, 1.7]
Zsmooth = scipy.ndimage.filters.gaussian_filter(h[0], sigma)
fig = go.Figure(data=[go.Surface(z=Zsmooth, x=X, y=Y)])
fig.update_layout(title='TITLE', autosize=False,
                  width=500, height=500,
                  margin=dict(l=65, r=50, b=65, t=90))
fig.show()

# %%
### Smoothing bins with kde

# %%
# smalladata = sc.pp.subsample(hoxb5, n_obs=1000, copy=True)
from scipy import stats
kde = stats.gaussian_kde(hoxb5.obsm['X_umap'].T, bw_method=0.2)

dens = kde(np.vstack((X.ravel(), Y.ravel())))
dens = dens.reshape((40,40), order='C')


fig = go.Figure(data=[go.Surface(z=dens, x=X, y=Y)])
fig.update_layout(title='TITLE', autosize=False,
                  width=500, height=500,
                  margin=dict(l=65, r=50, b=65, t=90))
fig.show()

# %% [markdown]
# ## Plotting smoothed cell density over time

# %%
cellnos = pd.read_csv('./procdata/05script/model_input_tompos.csv')
cellnos['cellno'] = cellnos['sc_ncells_filt'] / cellnos['sc_ncells_total'] * cellnos['flow_total']
cellnos = cellnos.loc[~cellnos.cellno.isna(),:]
cellnos = cellnos.groupby(by='time').cellno.mean()
cellnos

# %%
ni=60 #Number of intervals
pad = 6
X, Y = np.meshgrid(np.linspace(start=xc.min()-pad, stop=xc.max()+pad, num=ni), np.linspace(start=yc.min()-pad, stop=yc.max()+pad, num=ni), indexing='ij')

def get_density(adata, n=40, bw_method=0.2):
    kde = stats.gaussian_kde(adata.obsm['X_umap'].T, bw_method=bw_method)
    dens = kde(np.vstack((X.ravel(), Y.ravel())))
    dens = dens.reshape((n,n), order='C')
    return(dens)

for i in [3, 7, 27, 269]:
# for i in [3, 7
    a = hoxb5[hoxb5.obs.timepoint_tx_days == i,:]
    d = get_density(a, n=ni)
    cells = np.log10(d*cellnos[i]+2)
    cells[cells<0.1]=0.1 #this preventa a plotting artifact (a line close to 0 values))
    
    fig = go.Figure(data=[go.Surface(z=cells, x=X, y=Y,
                                     surfacecolor=d,
                                     colorscale='Turbo',
                                     showscale=False,
                                     cmax=260)])
    fig.update_layout(title=str(i), autosize=False,
                  width=500, height=500,
                  # margin=dict(l=65, r=50, b=65, t=90),
                  scene_camera=dict(up=dict(x=0, y=0, z=1),
                                    center=dict(x=0, y=0, z=0),
                                    eye=dict(x=0.3, y=-0.65, z=1.9)),
                  scene_dragmode='orbit',
                  scene = dict(
                      xaxis = dict(
                          backgroundcolor="rgba(0, 0, 0,0)",
                          gridcolor="white",
                          showbackground=True,
                          zerolinecolor="white",
                          visible=False,
                          range=[X.min()-1,X.max()+1]),
                      yaxis = dict(
                          backgroundcolor="rgba(0, 0, 0,0)",
                          gridcolor="white",
                          showbackground=True,
                          zerolinecolor="white",
                          visible=False,
                          range=[Y.min()-1,Y.max()+1]),
                      zaxis = dict(
                          backgroundcolor="rgba(0, 0, 0,0)",
                          gridcolor="white",
                          showbackground=True,
                          zerolinecolor="white",
                          visible=False,
                          range=[0,np.log10(275+2)])))
    fig.show()
    fig.write_image(base_figures + 'celldens_surface_' + str(i) + 'd.pdf', scale=2)


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

# %%
cols = get_colormap(hoxb5, obs_column='anno_man')
cols = [cols[i] for i in hoxb5.obs.anno_man]

# %%
fig = go.Figure(data=[go.Scatter3d(z=np.zeros(xc.shape), x=xc, y=yc,
                                   mode='markers',
                                   marker=dict(
                                   size=2,
                                   color=cols,                # set color to an array/list of desired values
                                   colorscale='Viridis',   # choose a colorscale
                                   opacity=1))])
fig.update_layout(autosize=False,
                  width=500, height=500,
                  # margin=dict(l=65, r=50, b=65, t=90),
                  scene_camera=dict(up=dict(x=0, y=0, z=1),
                                    center=dict(x=0, y=0, z=0),
                                    eye=dict(x=0.3, y=-0.65, z=1.9)),
                  scene_dragmode='orbit',
                  scene = dict(
                      xaxis = dict(
                          backgroundcolor="rgba(0, 0, 0,0)",
                          gridcolor="white",
                          showbackground=True,
                          zerolinecolor="white",
                          visible=False,
                          range=[X.min()-1,X.max()+1]),
                      yaxis = dict(
                          backgroundcolor="rgba(0, 0, 0,0)",
                          gridcolor="white",
                          showbackground=True,
                          zerolinecolor="white",
                          visible=False,
                          range=[Y.min()-1,Y.max()+1]),
                      zaxis = dict(
                          backgroundcolor="rgba(0, 0, 0,0)",
                          gridcolor="white",
                          showbackground=True,
                          zerolinecolor="white",
                          visible=False,
                          range=[0,np.log10(275+2)])))
fig.show()
fig.write_image(base_figures + 'annotation.pdf', scale=2)

# %%

# %%

# %%

# %%

ni=20 #Number of intervals
pad = 6
X, Y = np.meshgrid(np.linspace(start=xc.min()-pad, stop=xc.max()+pad, num=ni), np.linspace(start=yc.min()-pad, stop=yc.max()+pad, num=ni), indexing='ij')

def get_density(adata, n=40, bw_method=0.2):
    kde = stats.gaussian_kde(adata.obsm['X_umap'].T, bw_method=bw_method)
    dens = kde(np.vstack((X.ravel(), Y.ravel())))
    dens = dens.reshape((n,n), order='C')
    return(dens)

for i in [3, 7, 27, 269]:
# for i in [3, 7
    a = hoxb5[hoxb5.obs.timepoint_tx_days == i,:]
    d = get_density(a, n=ni)
    cells = np.log10(d*cellnos[i]+2)
    cells[cells<0.1]=0.1 #this preventa a plotting artifact (a line close to 0 values))
    
    fig = go.Figure(data=[go.Surface(z=d, x=X, y=Y,
#                                      surfacecolor=d,
#                                      colorscale='Turbo',
#                                      showscale=False,
                                     opacity=0.3)])
#     fig.update_layout(title=str(i), autosize=False,
#                   width=500, height=500,
#                   # margin=dict(l=65, r=50, b=65, t=90),
# #                   scene_camera=dict(up=dict(x=0, y=0, z=1),
# #                                     center=dict(x=0, y=0, z=0),
# #                                     eye=dict(x=0.3, y=-0.65, z=1.9)),
# #                   scene_dragmode='orbit',
# #                   scene = dict(
# #                       xaxis = dict(
# #                           backgroundcolor="rgba(0, 0, 0,0)",
# #                           gridcolor="white",
# #                           showbackground=True,
# #                           zerolinecolor="white",
# #                           visible=False,
# #                           range=[X.min()-1,X.max()+1]),
# #                       yaxis = dict(
# #                           backgroundcolor="rgba(0, 0, 0,0)",
# #                           gridcolor="white",
# #                           showbackground=True,
# #                           zerolinecolor="white",
# #                           visible=False,
# #                           range=[Y.min()-1,Y.max()+1]),
# #                       zaxis = dict(
# #                           backgroundcolor="rgba(0, 0, 0,0)",
# #                           gridcolor="white",
# #                           showbackground=True,
# #                           zerolinecolor="white",
# #                           visible=False,
# #                           range=[0,np.log10(275+2)]))
#                      )
    fig.show()
    fig.write_image(base_figures + 'celldens_surface_' + str(i) + 'd.pdf')

# %%
z1 = np.array([
    [8.83,8.89,8.81,8.87,8.9,8.87],
    [8.89,8.94,8.85,8.94,8.96,8.92],
    [8.84,8.9,8.82,8.92,8.93,8.91],
    [8.79,8.85,8.79,8.9,8.94,8.92],
    [8.79,8.88,8.81,8.9,8.95,8.92],
    [8.8,8.82,8.78,8.91,8.94,8.92],
    [8.75,8.78,8.77,8.91,8.95,8.92],
    [8.8,8.8,8.77,8.91,8.95,8.94],
    [8.74,8.81,8.76,8.93,8.98,8.99],
    [8.89,8.99,8.92,9.1,9.13,9.11],
    [8.97,8.97,8.91,9.09,9.11,9.11],
    [9.04,9.08,9.05,9.25,9.28,9.27],
    [9,9.01,9,9.2,9.23,9.2],
    [8.99,8.99,8.98,9.18,9.2,9.19],
    [8.93,8.97,8.97,9.18,9.2,9.18]
])

z2 = z1 + 1
z3 = z1 - 1

fig = go.Figure(data=[
    go.Surface(z=z3, surfacecolor=z3, colorscale='Turbo', showscale=False, opacity=0.999)

])

fig.show()

# %%

# %%

# %%

# %%
# #Getting coordinates for later (OLD APPROACH WITH 2D HISTOGRAM)
nbins=200 #Number of bins
h = plt.hist2d(hoxb5.obsm['X_umap'][:,0], hoxb5.obsm['X_umap'][:,1], bins=200)
xmids = np.array([(h[1][n] + h[1][n+1])/2 for n in range(len(h[1])-1)])
ymids = np.array([(h[2][n] + h[2][n+1])/2 for n in range(len(h[2])-1)])
X, Y = np.meshgrid(xmids, ymids, indexing='ij')
Z = h[0]

def get_density(adata, n=40, bw_method=0.2):
    kde = stats.gaussian_kde(adata.obsm['X_umap'].T, bw_method=bw_method)
    dens = kde(np.vstack((X.ravel(), Y.ravel())))
    dens = dens.reshape((n,n), order='C')
    return(dens)

for i in [3, 7, 27, 269]:
# for i in [3, 7
    a = hoxb5[hoxb5.obs.timepoint_tx_days == i,:]
    d = get_density(a, n=nbins)
    cells = np.log10(d*cellnos[i]+1)
    cells[cells<0.1]=0.1 #this preventa a plotting artifact (a line close to 0 values))
    
    fig = go.Figure(data=[go.Surface(z=cells, x=X, y=Y, surfacecolor=d, colorscale='Turbo', cmax=260)])
    fig.update_layout(title=str(i), autosize=False,
                  width=500, height=500,
                  # margin=dict(l=65, r=50, b=65, t=90),
                  scene_camera=dict(up=dict(x=0, y=0, z=1),
                                    center=dict(x=0, y=0, z=0),
                                    eye=dict(x=0.3, y=-0.65, z=1.9)),
                  scene_dragmode='orbit',
                  scene = dict(
                      xaxis = dict(
                          backgroundcolor="rgba(0, 0, 0,0)",
                          gridcolor="white",
                          showbackground=True,
                          zerolinecolor="white",
                          visible=False,
                          range=[X.min()-1,X.max()+1]),
                      yaxis = dict(
                          backgroundcolor="rgba(0, 0, 0,0)",
                          gridcolor="white",
                          showbackground=True,
                          zerolinecolor="white",
                          visible=False,
                          range=[Y.min()-1,Y.max()+1]),
                      zaxis = dict(
                          backgroundcolor="rgba(0, 0, 0,0)",
                          gridcolor="white",
                          showbackground=True,
                          zerolinecolor="white",
                          visible=False,
                          range=[0,np.log10(275+1)])))
    fig.show()
    fig.write_image(base_figures + 'celldens_surface_' + str(i) + 'd.pdf')

# %%

    # %%
    # fig2, ax = plt.subplots(subplot_kw={"projection": "3d"})
    # surf = ax.plot_surface(X, Y, d, cmap='coolwarm',
    #                    linewidth=10, antialiased=False)
    # ax.view_init(60, -80)
    # fig2.show()

# %%
# gx = np.array([[1,2,3], [4,5,6]])
# print(gx)
# gy = np.array([[10,20,30], [40,50,60]])
# print(gy)

# v = np.vstack((gx.ravel(), gy.ravel()))
# print(v)
# v = v.sum(axis=0)
# print(v)
# v.reshape((2,3), order='C')

# %%

# %%

# %%

# %%
# X, Y = np.mgrid[-1:1:20j, -1:1:20j]

# %%
# https://stackoverflow.com/questions/35157650/smooth-surface-plot-with-pyplot

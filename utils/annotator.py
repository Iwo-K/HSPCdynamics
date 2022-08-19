import scanpy as sc
import pandas as pd
import numpy as np
from anndata import AnnData


def check_markers_mincells(adata, min_cells=10):
    """
    Checking if there is more than 'min_cells' with >0
    expression (all genes considered)
    """

    cells_withexpr = sc.pp.filter_cells(adata, min_counts=1, inplace=False)
    cells_withexpr = cells_withexpr[0].sum()
    return(cells_withexpr >= min_cells)


def check_markers_meanexpr(adata, mean_tr=1):
    """Checking if mean expression is above threshold"""
    meanexpr = adata.X.mean()
    return meanexpr >= mean_tr


def annotate_lineages(adata,
                      group_key='leiden',
                      key_added='anno',
                      markers='LK',
                      min_cells=10,
                      mean_tr=1,
                      use_raw=True):
    """
    Basic cell population annotation. Looks at minimum number of cells expressing
    markers and whether the expression is above a pre-define threshold.
    Supply to markers either a pre-defined set of markers or a dictionary
    with the name: markers lists.
    Expect log-normalised input (uses .raw by default)
    """

    # Using sc.pp.filter_cells which prints a lot of text
    verbosity = sc.settings.verbosity
    sc.settings.verbosity = 1

    if markers == 'LK':
        # A few marker taken from Weinreb2020
        markers = dict(HSC=['Procr', 'Mllt3', 'Hoxb5', 'Mecom', 'Ly6a'],
                       Ery=['Klf1', 'Gata1', 'Hba-a1'],
                       MC=['Cma1', 'Gzmb'],
                       Bas=['Prss34', 'Mcpt8'],
                       Eos=['Prg2', 'Prg3'],
                       Mk=['Pf4', 'Ppbp', 'Vwf'],
                       Mo_DC=['Csf1r', 'Ctss', 'Ms4a6c', 'Irf8',
                              'Ccr2', 'F13a1'],
                       Neu=['Elane', 'Prtn3', 'Cebpe', 'Gfi1'],
                       Ly=['Dntt', 'Il7r', 'Jchain', 'Ighm', 'Flt3'],
                       Bcell=['Cd79a', 'Cd79b', 'Cd19'],
                       Tcell=['Cd3e', 'Bcl11b', 'Lck'])

    adata2 = adata
    if use_raw:
        adata_temp = AnnData(X=adata.raw.X, var=adata.raw.var, obs=adata.obs)

    cellgroups = adata_temp.obs[group_key].unique()
    anno = {i : [] for i in cellgroups}

    for lineage, genes in markers.items():
        x = adata_temp[:, genes].copy()
        x_scaled = sc.pp.scale(x, copy=True)
        for group in cellgroups:
            x_group = x[x.obs.index[x.obs[group_key] == group],:].copy()
            x_group_scaled = x_scaled[x_group.obs.index,:].copy()
            # Checking if sufficient number of cells express markers
            if check_markers_mincells(x_group, min_cells=min_cells):
                if check_markers_meanexpr(x_group_scaled, mean_tr=mean_tr):
                    anno[group].append(lineage)

    sc.settings.verbosity = verbosity
    anno = {k : '/'.join(v) for k, v in anno.items()}
    adata.obs[key_added] = [anno[i] for i in adata.obs[group_key]]
    adata.obs.loc[adata.obs[key_added] == '', key_added] = np.nan

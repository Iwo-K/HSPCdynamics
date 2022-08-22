import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.distance import euclidean
import scrublet as scr
from scipy.sparse import issparse, csr_matrix
from scipy.stats import pearsonr
import anndata

#NOTE: most of these function work with mouse data and needs adapting to human datasets

def find_correlated_genes(adata, genes, data_scaled = True, min_corr = 0.1):
    if not data_scaled:
        adata = adata.copy()
        sc.pp.scale(adata)

    X = adata.X
    if issparse(X):
        X = X.toarray()

    genes = pd.Series(genes)
    genes = genes[genes.isin(adata.var.index)]

    score = X[:,adata.var.index.isin(genes)]
    score = score.sum(axis=1)
 
    c = np.zeros(len(adata.var.index))
    for i in range(len(c)):
        c[i],_ = pearsonr(score, X[:,i])

    return adata.var.index[c >= min_corr].values


def do_qc(adata, min_genes = 1000, min_counts = 0, mt_frac = 0.1, species = "mouse"):
    ''' Performs basic QC of anndata, uses qc_plots to print diagnostics'''

    print("BEFORE QC")
    qc_plots(adata, species = species)

    sc.pp.filter_cells(adata, min_genes = min_genes)
    sc.pp.filter_cells(adata, min_counts = min_counts)

    print('Cells with mitochondrial counts> 10%:' + str(sum(adata.obs['mt_frac'] > 0.1)))
    #Skipping mitochondrial count filtering
    #adata = adata[~(adata.obs['mt_frac'] > mt_frac),:]

    print("AFTER QC")
    qc_plots(adata, species = species)

def qc_plots(data, species = "mouse"):
    '''Prints plots with no of counts and no of genes per cell'''
    sc.pp.filter_cells(data, min_genes = 0)
    sc.pp.filter_cells(data, min_counts = 0)

    fig, axes = plt.subplots(nrows = 1, ncols = 3)
    fig.set_size_inches(15, 4)

    axes[0].hist(np.log10(data.obs.n_counts+1), range = [0,4.75], bins = 50, facecolor='green', alpha=0.75)
    axes[0].set(title = 'Log10 No of counts per cell + 1')

    ngenes = data.obs['n_genes']
    axes[1].hist(ngenes, 50, range = [0,6000], facecolor='green', alpha=0.75)
    axes[1].set(title = 'No of genes per cell')

    #Printing fraciton of mitochondrial genes
    import re
    if species == "mouse":
        mts = list(filter(lambda x: re.search('^mt-', x), data.var.index))
    elif species == "human":
        mts = list(filter(lambda x: re.search('^MT-', x), data.var.index))
    else: print("species needs to be either mouse or human")

    mt = data[:,mts].copy()
    if issparse(mt.X):
        mt.X = mt.X.toarray()
    data.obs['mt_count'] = mt.X.sum(axis = 1)
    data.obs['mt_frac'] = data.obs['mt_count'] / data.obs['n_counts']
    axes[2].hist(data.obs['mt_frac'], 50, range=[0, 0.2], density=1, facecolor='green', alpha=0.75)
    axes[2].set(title = 'Fraction of counts in mitochondrial genes')
    plt.show()

def doubletperrun(adata, obscol):
    '''Loops through batches provided in obscol (subsetting data) and calculates doublet scores with scrublet
    Returns a dataframe for all input cells with doublet scores and predicted doublets (boolean)'''
    adata = adata.copy()
    sc.pp.filter_genes(adata, min_counts=1)

    scoreDF = pd.DataFrame(columns=['run' , 'doublet_score', 'predicted_doublets'])
    runs =  pd.unique(adata.obs[obscol])

    for i in runs:
        print(i)
        tempdata = adata[adata.obs[obscol] == i, :].copy()
        print('Captured cell no: ', + tempdata.shape[0])

        scrub = scr.Scrublet(tempdata.X.todense(), expected_doublet_rate=0.1)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
        scores = pd.DataFrame({'run' : tempdata.obs[obscol], 'doublet_score' : doublet_scores,
                               'predicted_doublets' : predicted_doublets}, index = tempdata.obs.index)

        scoreDF = pd.concat((scoreDF, scores))
        print(scoreDF.shape)

        scrub.plot_histogram()
        plt.show()
    return(scoreDF)


import anndata
def remove_doublets(adata, tr = 0.25, qtr = None, groupby = 'batch'):
    if not 'doublet_score' in adata.obs.columns:
        raise Exception('Error: no doublet score column detected in the .obs slot')
    if not groupby in adata.obs.columns:
        raise Exception('Error: groupby column not found in the .obs slot')
    if (tr is not None) and (qtr is not None):
        raise Exception('Error: cannot provide both tr and qtr arguments')
    
    if tr is not None:
        print('Filtering with common threshold: ' + str(tr))
        subs = []
        for i in adata.obs[groupby].unique():
            adata_temp = adata[adata.obs.index[adata.obs[groupby] == i],:].copy()
            tokeep = adata_temp.obs.index[adata_temp.obs.doublet_score < tr]
            print(f'Batch {i}: filtering out {adata_temp.n_obs - len(tokeep)} cells')
            adata_temp = adata_temp[tokeep, :]
            subs.append(adata_temp)
    
    if qtr is not None:
        print('Filtering with a quantile threshold: ' + str(qtr))

        subs = []
        for i in adata.obs[groupby].unique():
            adata_temp = adata[adata.obs.index[adata.obs[groupby] == i],:].copy()
            tokeep = adata_temp.obs.index[adata_temp.obs.doublet_score < adata_temp.obs.doublet_score.quantile(1-qtr)]
            trval = adata_temp.obs.doublet_score.quantile(1-qtr)
            print(f'Batch {i}: filtering out {adata_temp.n_obs - len(tokeep)} cells with a threshold of {trval}')
            adata_temp = adata_temp[tokeep, :]
            subs.append(adata_temp)

    fin_adata = anndata.concat(subs, label = 'remove_doublets_groups')
    fin_adata.var = adata.var.copy()
    return(fin_adata)

def assign_sex(data, ygenes, use_raw = True):
    '''Assigns sex based on Xist and Y-chromosome gene expression. Expects log-normalised input'''

    if use_raw == True:
        xist_expr = data.raw[:,'Xist'].X.copy()
        if issparse(xist_expr):
            xist_expr = xist_expr.toarray()
        data.obs['Xist_logn'] = xist_expr[:,0]
        
        ygenes = data.raw.var.index.isin(ygenes)
        yexpr = data.raw.X[:,ygenes]
        if issparse(yexpr): 
            yexpr = yexpr.toarray()
    else:
        xist_logn = data[:,'Xist'].X.copy()
        if issparse(xist_logn): 
            xist_logn = xist_logn.toarray()
        data.obs['xist_logn'] = xist_logn[:,0]

        ygenes = data.var.index.isin(ygenes)
        yexpr = data.X[:,ygenes]
        if issparse(yexpr): 
            yexpr = yexpr.toarray()

    yexpr = yexpr.mean(axis = 1)
    data.obs['Ygene_logn'] = yexpr
    
    plt.scatter(data.obs.xist_logn, data.obs.Ygene_logn, alpha = 0.05)
    plt.show()

    data.obs['xist_bin'] = data.obs['xist_logn'] > 0
    data.obs['Ygene_bin'] = data.obs['Ygene_logn'] > 0
    print(pd.crosstab(index = data.obs.xist_bin, columns = data.obs.Ygene_bin))

    data.obs['sex'] = 'unassigned'
    data.obs.loc[data.obs.xist_bin & ~data.obs.Ygene_bin, data.obs.columns == 'sex'] = 'female'
    data.obs.loc[data.obs.xist_bin & data.obs.Ygene_bin, data.obs.columns == 'sex'] = 'mfdoublet'
    data.obs.loc[~data.obs.xist_bin & data.obs.Ygene_bin, data.obs.columns == 'sex'] = 'male'

def lognorm(adata):
    '''log-normalise adata (counts per 10000)'''
    sc.pp.normalize_per_cell(adata, counts_per_cell_after = 10000) # normalize with total UMI count per cell
    sc.pp.log1p(adata)

def delognorm(adata, counts_per_cell_after=10000):
    """
    Returning counts from log-normalised adata
    Requires the total number of counts in the adata.obs.n_counts
    """

    if 'n_counts' not in adata.obs.columns:
        raise Exception('n_counts not found in adata.obs')

    ratio = (adata.obs.n_counts/counts_per_cell_after).values

    if issparse(adata.X):
        X = adata.X.toarray()
    else:
        X = adata.X

    X = np.asarray(np.exp(X)-1)
    X = ratio[:, np.newaxis] * X
    X = np.around(X, decimals=0)
    if issparse(adata.X):
        X = csr_matrix(X)

    adata.X = X
    return adata


def process(adata,
            n_pcs = 50,
            n_neighbors = 12,
            compute_hivar = True ,
            leiden_resolution = 1,
            n_variable = 5000,
            lognorm_data = True,
            remove_genes = None,
            Sgenes = None,
            G2Mgenes = None,
            regress_cc = False,
           batch_correct = None,
           batch_key = 'batch',
           slim = False,
           copy = False):
    """Pre-processed the data as follows: normalisatio to 10k reads (optional), log-transformation,
    (optional) selection of variable genes, removal of selected genes from variable list (optional), scaling,
    assignment of cell cycle scores (optional), regression of cell cycle effect (optional),
    pca, batch correction with bbknn or harmony (optional), nearest neighbour selection, clustering, umap calculation
    CAUTION: All operations in-place
    CAUTION: slim means that while running data is subset for highly variable genes, .X is reduced in size but .raw still contains all logn value.
    Copy is returned and not modified in place.
    """
    if slim or copy:
        adata = adata.copy()
    
    if lognorm_data:
        sc.pp.normalize_per_cell(adata, counts_per_cell_after = 10000) # normalize with total UMI count per cell
        sc.pp.log1p(adata)

    logn_bak = adata.copy()

    if compute_hivar:
        sc.pp.highly_variable_genes(adata, flavor='seurat', n_top_genes = n_variable)
        sc.pl.highly_variable_genes(adata)
    if remove_genes is not None:
        print("Removing: " + str(len(remove_genes)) + " genes")
        adata.var['gene_removed'] = adata.var.index.isin(remove_genes)
        adata.var.loc[adata.var.gene_removed, adata.var.columns == 'highly_variable'] = False
        print(sum(adata.var.highly_variable))
    
    sc.pp.scale(adata)

    if (Sgenes is not None) and (G2Mgenes is not None):
        Sgenes = Sgenes.loc[Sgenes.isin(adata.var.symbol)].copy()
        G2M = G2Mgenes.loc[G2Mgenes.isin(adata.var.symbol)].copy()
        score_genes_cell_cycle_fixed(adata, Sgenes, G2Mgenes, use_raw = False)
        
    #A memory and runtime saving options which subsets data for highly variable,
    #It in that case the aata at the end will have logn values
    if slim:
        print('slim is True, final .X will contain only variable genes')
        adata = adata[:, adata.var.index[adata.var.highly_variable]].copy()
        
    if regress_cc:
        adata.X = logn_bak[:,adata.var.index].X.copy()
        ##it does not look like there is a big difference whether the regression is done on scale or lognn data, at least on a small scale
        sc.pp.regress_out(adata, ['S_score', 'G2M_score'])
        sc.pp.scale(adata)

    adata.raw = logn_bak
        
    sc.tl.pca(adata, n_comps=n_pcs, svd_solver = 'arpack')
    
    if batch_correct  == 'bbknn':
        sc.external.pp.bbknn(adata, n_pcs = n_pcs, batch_key = batch_key, neighbors_within_batch = 5)
    elif batch_correct == 'harmony':
        #This will require installing harmonypy
        sc.external.pp.harmony_integrate(adata, key = batch_key)
        sc.pp.neighbors(adata, n_pcs = n_pcs, n_neighbors = n_neighbors, use_rep = 'X_pca_harmony')
    else:
        sc.pp.neighbors(adata, n_pcs = n_pcs, n_neighbors = n_neighbors, use_rep = 'X_pca')
        
    sc.tl.leiden(adata, resolution = leiden_resolution)
    sc.tl.umap(adata)

    if slim or copy:
        return adata
    

def subsample_adata(adata, sub):
    ''' Subsample adata to a given fraction of reads, return a copy'''
    subn  = adata.obs.n_counts*sub
    subn = subn.round()
    subn = subn.astype(int)
    subn = subn.values
    adata = sc.pp.downsample_counts(adata, counts_per_cell = subn, copy = True)
    adata.X = adata.X.astype('float32')
    return(adata)



def subsample_each(adata, group, n_obs = 1000):
    '''Subsample adata for each level in the group, group indicates a .obs column'''
    import anndata

    subs = []
    for i in adata.obs[group].unique():
        adata_temp = adata[adata.obs.index[adata.obs[group] == i],:].copy()
        sc.pp.subsample(adata_temp, n_obs = n_obs)
        subs.append(adata_temp)
    
    return(anndata.concat(subs, label = 'groups_concat'))

def subsample_tomax(adata, group, n_obs_max=8000):
    """Subsample adata for each level in the group if exceeding n_obs_max"""
    import anndata

    subs = []
    for i in adata.obs[group].unique():
        adata_temp = adata[adata.obs.index[adata.obs[group] == i],:].copy()
        if adata_temp.n_obs > n_obs_max:
            sc.pp.subsample(adata_temp, n_obs=n_obs_max)
        subs.append(adata_temp)
    
    return(anndata.concat(subs, label = 'groups_concat'))


############################ SCANPY FUNCTION FIXES ################################
#This is longer required in newer version of scanpy (used in container v4), kept here just in case
def score_genes_fixed(
        adata = None,
        gene_list = None,
        ctrl_size = 50,
        gene_pool = None,
        n_bins = 25,
        score_name = 'score',
        random_state = 0,
        copy = False,
        use_raw = False):
    """\
    Fixed version of score_gene function, which uses float64 precision.
    The original implementation from scanpy does not specify precision
    and can give rise to slightly different outputs on different machines,
    even when using same containers.

    Score a set of genes [Satija15]_.
    The score is the average expression of a set of genes subtracted with the
    average expression of a reference set of genes. The reference set is
    randomly sampled from the `gene_pool` for each binned expression value.
    This reproduces the approach in Seurat [Satija15]_ and has been implemented
    for Scanpy by Davide Cittaro.
    Parameters
    ----------
    adata
        The annotated data matrix.
    gene_list
        The list of gene names used for score calculation.
    ctrl_size
        Number of reference genes to be sampled from each bin. If `len(gene_list)` is not too
        low, you can set `ctrl_size=len(gene_list)`.
    gene_pool
        Genes for sampling the reference set. Default is all genes.
    n_bins
        Number of expression level bins for sampling.
    score_name
        Name of the field to be added in `.obs`.
    random_state
        The random seed for sampling.
    copy
        Copy `adata` or modify it inplace.
    use_raw
        Whether to use `raw` attribute of `adata`. Defaults to `True` if `.raw` is present.
        .. versionchanged:: 1.4.5
           Default value changed from `False` to `None`.
    Returns
    -------
    Depending on `copy`, returns or updates `adata` with an additional field
    `score_name`.
    Examples
    --------
    See this `notebook <https://github.com/theislab/scanpy_usage/tree/master/180209_cell_cycle>`__.
    """
    adata = adata.copy() if copy else adata

    if random_state is not None:
        np.random.seed(random_state)

    gene_list_in_var = []
    var_names = adata.raw.var_names if use_raw else adata.var_names
    genes_to_ignore = []
    for gene in gene_list:
        if gene in var_names:
            gene_list_in_var.append(gene)
        else:
            genes_to_ignore.append(gene)
    gene_list = set(gene_list_in_var[:])

    if len(gene_list) == 0:
        raise ValueError("No valid genes were passed for scoring.")

    if gene_pool is None:
        gene_pool = list(var_names)
    else:
        gene_pool = [x for x in gene_pool if x in var_names]

    # Trying here to match the Seurat approach in scoring cells.
    # Basically we need to compare genes against random genes in a matched
    # interval of expression.

    # use_raw = _check_use_raw(adata, use_raw)
    _adata = adata.raw if use_raw else adata

    _adata_subset = _adata[:, gene_pool] if len(gene_pool) < len(_adata.var_names) else _adata
    if issparse(_adata_subset.X):
        obs_avg = pd.Series(
            np.array(_sparse_nanmean(_adata_subset.X, axis=0)).flatten(), index=gene_pool)  # average expression of genes
    else:
        obs_avg = pd.Series(
            np.nanmean(_adata_subset.X, axis=0), index=gene_pool)  # average expression of genes

    obs_avg = obs_avg[np.isfinite(obs_avg)] # Sometimes (and I don't know how) missing data may be there, with nansfor

    n_items = int(np.round(len(obs_avg) / (n_bins - 1)))
    obs_cut = obs_avg.rank(method='min') // n_items
    control_genes = set()

    # now pick `ctrl_size` genes from every cut
    for cut in np.unique(obs_cut.loc[gene_list]):
        r_genes = np.array(obs_cut[obs_cut == cut].index)
        np.random.shuffle(r_genes)
        # uses full r_genes if ctrl_size > len(r_genes)
        control_genes.update(set(r_genes[:ctrl_size]))

    # To index, we need a list – indexing implies an order.
    control_genes = list(control_genes - gene_list)
    gene_list = list(gene_list)

    X_list = _adata[:, gene_list].X
    if issparse(X_list):
        X_list = np.array(_sparse_nanmean(X_list, axis=1)).flatten()
    else:
        X_list = np.nanmean(X_list, axis=1, dtype = 'float64')

    X_control = _adata[:, control_genes].X
    if issparse(X_control):
        X_control = np.array(_sparse_nanmean(X_control, axis=1)).flatten()
    else:
        X_control = np.nanmean(X_control, axis=1, dtype = 'float64')

    score = X_list - X_control

    adata.obs[score_name] = pd.Series(np.array(score).ravel(), index=adata.obs_names).astype('float32')

    return adata if copy else None


def score_genes_cell_cycle_fixed(
        adata,
        s_genes,
        g2m_genes,
        copy = False,
    **kwargs):
    """\
    Modified version of scanpy's function that uses the score_genes_fixed
    function instead of score_genes.

    Score cell cycle genes [Satija15]_.
    Given two lists of genes associated to S phase and G2M phase, calculates
    scores and assigns a cell cycle phase (G1, S or G2M). See
    :func:`~scanpy.tl.score_genes` for more explanation.
    Parameters
    ----------
    adata
        The annotated data matrix.
    s_genes
        List of genes associated with S phase.
    g2m_genes
        List of genes associated with G2M phase.
    copy
        Copy `adata` or modify it inplace.
    **kwargs
        Are passed to :func:`~scanpy.tl.score_genes`. `ctrl_size` is not
        possible, as it's set as `min(len(s_genes), len(g2m_genes))`.
    Returns
    -------
    Depending on `copy`, returns or updates `adata` with the following fields.
    **S_score** : `adata.obs`, dtype `object`
        The score for S phase for each cell.
    **G2M_score** : `adata.obs`, dtype `object`
        The score for G2M phase for each cell.
    **phase** : `adata.obs`, dtype `object`
        The cell cycle phase (`S`, `G2M` or `G1`) for each cell.
    See also
    --------
    score_genes
    Examples
    --------
    See this `notebook <https://github.com/theislab/scanpy_usage/tree/master/180209_cell_cycle>`__.
    """

    adata = adata.copy() if copy else adata
    ctrl_size = min(len(s_genes), len(g2m_genes))
    # add s-score
    score_genes_fixed(adata, gene_list=s_genes, score_name='S_score', ctrl_size=ctrl_size, **kwargs)
    # add g2m-score
    score_genes_fixed(adata, gene_list=g2m_genes, score_name='G2M_score', ctrl_size=ctrl_size, **kwargs)
    scores = adata.obs[['S_score', 'G2M_score']]

    # default phase is S
    phase = pd.Series('S', index=scores.index)

    # if G2M is higher than S, it's G2M
    phase[scores.G2M_score > scores.S_score] = 'G2M'

    # if all scores are negative, it's G1...
    phase[np.all(scores < 0, axis=1)] = 'G1'

    adata.obs['phase'] = phase
    return adata if copy else None

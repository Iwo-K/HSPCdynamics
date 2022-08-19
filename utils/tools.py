def clip(x, min = -3, max = 3):
    '''Clip values to a given min and max values in an array'''
    x = x.copy()
    x[x < min] = min
    x[x > max] = max
    return(x)

def cell_kde_pair(adata,
                      groupby = 'condition',
                      ref = 'WT',
                      coords = 'X_umap',
                      bw_method = None,
                      estimate_cell_densities = True,
                      output = 'ratio'):
    '''Takes in adata and densities for two categories in the groupby columns,
    returns densities, difference or ratio between them. Densities are
    calculated for the split data and then estimated for the whole data provided.
    if estimate_cell_densities is False it returns the density functions'''

    conds = adata.obs[groupby].unique()
    if len(conds) !=2:
        print("Error: need only two conditions (levels)")
        return()
    treat = conds[conds != ref][0]
    print(ref)
    print(treat)

    refadata = adata[adata.obs.index[adata.obs[groupby] == ref],:]
    treatadata = adata[adata.obs.index[adata.obs[groupby] == treat],:]

    refdens = stats.gaussian_kde(refadata.obsm[coords].T, bw_method = bw_method)
    treatdens = stats.gaussian_kde(treatadata.obsm[coords].T, bw_method = bw_method)

    if estimate_cell_densities == False:
        return({treat : treatdens, ref : refdens})
    else:
        treatdens = treatdens(adata.obsm[coords].T)
        refdens = refdens(adata.obsm[coords].T)

        if output == 'diff':
            return(treatdens - refdens)
        elif output == 'ratio':
            ratio = treatdens/refdens
            return(ratio)
        elif output == 'dens':
            return((treatdens, refdens))
        else: print("Error: output can be either diff or ratio")





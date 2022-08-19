import matplotlib.pyplot as plt
import plotly as plotly
import plotly.graph_objs as go
import scanpy as sc
import pandas as pd
import numpy as np
from scipy.sparse import issparse


def umap_subgroups(adata, key, toplot, alpha=0.5, size=16, subsample_to = None, file=None):

    n = len(toplot)
    nrows = -(-n // 4) # hack for ceil division
    print(nrows)
    fig, axes = plt.subplots(nrows=nrows, ncols=4)
    fig.set_size_inches(20, nrows * 5)
    axes = axes.ravel()
    
    for n,i in enumerate(toplot):
        sc.pl.umap(adata, show=False, ax=axes[n], size = 16, alpha = alpha)
        #In the old version of scanpy there is no na_color and the palette
        # argument does not work for me, As a workaround copying the object
        # and plotting in default blue
        # In newer just set the na_color to a given value
        temp = adata[adata.obs.index[adata.obs[key] == i],:].copy()
        temp.obs['temp'] = 'a'
        if subsample_to is not None:
            if not isinstance(subsample_to, dict):
                raise Exception("subsample_to is expected to be a dictionary")
            sc.pp.subsample(temp, n_obs=subsample_to[i])
        print(temp.n_obs)
        sc.pl.umap(temp, color='temp', show=False, ax=axes[n], size=size, legend_loc=None, title=i)

        axes[n].margins(0.05)

    fig.tight_layout()
    fig.show()
    if file is not None:
        plt.savefig(file)


def umap_cell_abundance(adata, key, subsample_to = None, file = None, density_plots = False):
    import anndata as an
    import seaborn as sns
    from scanpy.plotting.palettes import default_20

    catno = len(adata.obs[key].unique())
    
    #Getting the limits
    xlim = (1.1*min(adata.obsm['X_umap'][:,0]), 1.1*max(adata.obsm['X_umap'][:,0]))
    ylim = (1.1*min(adata.obsm['X_umap'][:,1]), 1.1*max(adata.obsm['X_umap'][:,1]))

    if catno <= 8:
        fig, axes = plt.subplots(nrows=1, ncols=catno)
        fig.set_size_inches(7*catno, 5)

        x = []
        for n, i in enumerate(adata.obs[key].cat.categories):
            sub = adata[adata.obs[key].isin([i])].copy()
            if subsample_to is not None:
                sc.pp.subsample(sub, n_obs = subsample_to)
                if density_plots:
                    #Note this kdeplot syntax works with seaborn 10.1 but not 11!
                    sns.kdeplot(data = sub.obsm['X_umap'][:,0], data2 = sub.obsm['X_umap'][:,1], shade = True, shade_lowest = False,
                            ax = axes[n]) 
                    # color = default_20[n]
                    axes[n].set_title(i)
                else:
                    sc.pl.umap(sub, color = [key], 
                           alpha = 0.7,
                           ax = axes[n], show = False, frameon = True)
            axes[n].set_xlim(xlim)
            axes[n].set_ylim(ylim)
            for i in axes: i.axes.set_box_aspect(1)

        
            fig.tight_layout()   
            plt.savefig(file)
    else: print('Sorry, only up to 8 categories are supported')

def get_cell_abundance(adata, key, mincells = 0, dropna = True):
    ''' Returns a normalised number of cells for each category in key, this is number of cells per 10,000 cells sampled'''
    exp_no = adata.obs[key].value_counts(sort = False)
    leiden_no = adata.obs['leiden'].value_counts(sort = False)

    op = pd.crosstab(adata.obs[key], adata.obs.leiden, dropna = dropna)
    print(op)
    op = op.loc[exp_no.index, leiden_no.index]

    op.loc[:,op.sum(axis = 0) < mincells] = np.nan
    oprel1 = op.T/op.sum(axis = 1)*10000

    return(oprel1)

def heatmap_cell_abundance(adata, key, file = None, logcolor = False, vmin = None, vmax = None):
    import seaborn as sns
    catno = len(adata.obs[key].unique())
    sns.set_style("whitegrid", {'axes.grid' : False})

    #Calculating relative cell abundance
    oprel1 = get_cell_abundance(adata = adata, key = key)

    from matplotlib.colors import LinearSegmentedColormap
    plt.figure(figsize=(19,6))
    sns.set(font_scale=1.4)

    if logcolor:
        oprel2 = np.log2(oprel1.T/oprel1.mean(axis = 1)+0.01)
        cmap2 = cmap_RdBu2(oprel2, vmin = vmin, vmax = vmax)

        sns.heatmap(oprel2, cmap = cmap2, annot = oprel1.T, fmt='.1f', annot_kws={"size": 9}, vmin = vmin, vmax = vmax)

    else:
        oprel2 = oprel1.T/oprel1.T.sum(axis=0) #Normalising all as expected fractions. So if there are 2 samples the expected is 0.5, if 3 then 0.33 etc.    
        cmap2 = LinearSegmentedColormap.from_list('mycmap', [(0, 'blue'),
                                                             (1/catno, 'white'),
                                                             (1, 'red')])
        vmin = 0
        vmax = 1
        sns.heatmap(oprel2, cmap = cmap2, annot = oprel1.T, fmt='.1f', annot_kws={"size": 9}, vmin = vmin, vmax = vmax)

    plt.savefig(file)
    plt.show()

def cmap_RdBu2(values, vmin = None, vmax = None):
    """Generates a blue/red colorscale (potentially asymmetric) with white value centered around the value 0, needs both negative and positive values.

    Parameters
    ----------
    values : 2d numpy array
        List of values to be used for creating the color map
    vmin : type
        Minimum value in the color map, if None then the min(values) is used
    vmax : type
        Maximum value in the color map, if None then the max(values) is used

    Returns
    -------
    type
        Description of returned object.

    """
    if vmin != None:
        scoremin = vmin
    else:
        scoremin = values.min().min()
    if vmax != None:
        scoremax = vmax
    else:
        scoremax = values.max().max()

    from matplotlib.colors import LinearSegmentedColormap

    cmap2 = LinearSegmentedColormap.from_list('mycmap', [(0, 'blue'),
                                                        (-scoremin/(scoremax-scoremin), 'white'),
                                                        (1, 'red')])
    return(cmap2)


######## 3d plots ########
def umap3d(adata, color, key = 'X_umap', filename = None, colorscale = 'Viridis', range_color = None,
           cmid = None,
           cmax = None,
           cmin = None):
    """Create colourcoded 3d UMAP plots with plotly.

    Outputs 3d plots (using plotly framework) of the UMAp representation. Cells
    are colourcoded according to either gene expression or a column in the adata.obs
    slot.

    Parameters
    ----------
    adata : AnnData object
        AnnData objects with gene expression values and UMAP coordinates
    color : str
        Name of a column in the adata.obs (will be treated as categorical) slot
        or name of the gene (in adata.var) to be used for colourcoding
    key : str
        Name of the dimensionality reduction coordinates to be used, defaults to 'X_umap'
    filename : str
        If filename is provided, the plot is saved to a filename.html file,
        otherwise the plot will be lanuched inside of jupyter notebook

    Returns
    -------
    Output plots, eiher to the jupyter notebook or to the indicated file.

    """

    if key == "X_diffmap":
        x, y, z = adata.obsm[key][:, 1], adata.obsm[key][:, 2], adata.obsm[key][:, 3]
    else:
        x, y, z = adata.obsm[key][:, 0], adata.obsm[key][:, 1], adata.obsm[key][:, 2]

    traces = []

    if color in adata.obs.columns:
        print("testing")
        #Consider using this code for a different colour pallete
        # color = adata.obs[color].astype('str').values
        # color = pd.factorize(color)[0]
        # pal = plt.get_cmap('tab10')
        # color = [pal(x) for x in color]

        if hasattr(adata.obs[color], 'cat'):

            for n, i in enumerate(adata.obs[color].cat.categories):

                xI = x[adata.obs[color] == i]
                yI = y[adata.obs[color] == i]
                zI = z[adata.obs[color] == i]

                trace1 = go.Scatter3d( 
                    x=xI,
                    y=yI,
                    z=zI,
                    mode='markers',
                    marker=dict(
                        size=5,
                        color = adata.uns[color + '_colors'][n]
                    ),
                    text = adata.obs.index[adata.obs[color] == i],
                    name = str(i)
            )

                traces.append(trace1)

                data = traces
                
        elif adata.obs[color].dtype == 'float32' or adata.obs[color].dtype == 'float64' or adata.obs[color].dtype == 'int':

            color = adata.obs[color]

            traces = go.Scatter3d(
                x=x,
                y=y,
                z=z,
                text = [str(x) for x in color],
                mode='markers',
                marker=dict(
                    size=5,
                    color=color,                # set color to an array/list of desired values
                    colorscale=colorscale)# choose a colorscale
            )
            #assembling the traces into data
            data = [traces]

        else: print("The column in the .obs slot needs to be either category, int or float32/float64")

    elif color in adata.var.index:

        color = adata.var.index.get_loc(color)
        color = adata.X[:,color]
        if issparse(color):
            color = color.todense()
            color = color.A1
        
        traces = go.Scatter3d(
            x=x,
            y=y,
            z=z,
            text = [str(x) for x in color],
            mode='markers',
            marker=dict(
                size=5,
                color=color,                # set color to an array/list of desired values
                colorscale=colorscale,   # choose a colorscale
                )
            )
        #assembling the traces into data
        data = [traces]

    else:
        print("Color not found in adata.obs or in the adata.var")
        return

    #Adjusting the colorscale
    if cmid != None: #setting the center of the continuous color scale (for instance with negative and positive values one often wants this at 0)
         traces.marker['cmid'] = cmid
    if cmin != None:
         traces.marker['cmin'] = cmin
    if cmax != None:
         traces.marker['cmax'] = cmax


    layout = go.Layout(
        margin=dict(
            l=0,
            r=0,
            b=0,
            t=30
        )
    )

    fig = go.Figure(data=data)#, layout=layout)

    if filename is not None:
        plotly.offline.plot(fig, filename=filename, auto_open=False)
    else:
        return(fig)


def get_cluster_coord(adata, group = 'leiden', key = 'X_umap'):
    '''Extracts mean coordinate positions for a given cluster'''
    dims = adata.obsm[key].shape[1]
    cl_coord = np.empty((0, dims))

    for i in adata.obs[group].cat.categories:
        a = adata[adata.obs[group] == i,:].copy()
        coord = a.obsm[key].mean(axis = 0)
        cl_coord = np.append(cl_coord, np.array([coord]), axis = 0)
        
    return(cl_coord)


def plot_paga3d(adata, edge_tr = 0.15, cluster_colors = None, filename = None):
    
    #Getting average cluster coordinates (needs to be 3d!)
    cluster_coords = get_cluster_coord(adata)
    #Getting cluster colors
    if cluster_colors == None:
        cluster_colors = adata.uns['leiden_colors']
    
    #Extracting PAGA connectivity (adjacency matrix) and converting into an annotated DataFrame
    con = adata.uns['paga']['connectivities'].todense()
    
    node1seq = np.repeat(list(range(0, con.shape[0])), con.shape[0])
    node2seq = list(range(0, con.shape[0]))*con.shape[0]
    df = pd.DataFrame({'node1' : node1seq, 'node2' : node2seq, 'weight' : con.A1})

    #Filtering low weight edges
    df = df.loc[df.weight > edge_tr]
    
    #Adding the cluster coordinates to the DataFrame
    df['node1_Xcoord'] = cluster_coords[df['node1'], 0]
    df['node1_Ycoord'] = cluster_coords[df['node1'], 1]
    df['node1_Zcoord'] = cluster_coords[df['node1'], 2]

    df['node2_Xcoord'] = cluster_coords[df['node2'], 0]
    df['node2_Ycoord'] = cluster_coords[df['node2'], 1]
    df['node2_Zcoord'] = cluster_coords[df['node2'], 2]
    
    #Setting up lists with X, Y, Z coordinates of edges. Plotly does not support plotting graphs in 3d
    #so this uses a little hack by providing a list of point coordinates, each pair separated by the None value
    #This wasy pairs of points are connected but nothing else.
    Xe = []
    Ye = []
    Ze = []
    for i in df.index:
        Xe+= [df.loc[i, 'node1_Xcoord'], df.loc[i, 'node2_Xcoord'], None]
        Ye+= [df.loc[i, 'node1_Ycoord'], df.loc[i, 'node2_Ycoord'], None]
        Ze+= [df.loc[i, 'node1_Zcoord'], df.loc[i, 'node2_Zcoord'], None]
        
    # Plotting the graph
    
    #Plotly does not support different line width, only one value can be provided at a time. Thus
    #I am plotting each edge as a separate trace. For this purpose I am looping through the Xe,Ye,Ze lists
    # taking paired values at the time and passing the width.
    traces = []
    for n, i in enumerate(range(0, len(Xe), 3)):
        a = go.Scatter3d(x=Xe[i:i+3],
                   y=Ye[i:i+3],
                   z=Ze[i:i+3],
                   mode='lines',
                   line=dict(color='rgb(125,125,125)', width=df.loc[df.index[n],'weight']*15),
                   hoverinfo='none'
                   )
        traces.append(a)

    #Plotting the points
    trace2=go.Scatter3d(x=cluster_coords[:,0],
                            y=cluster_coords[:,1],
                            z=cluster_coords[:,2],
                   mode='markers',
                   marker=dict(symbol='circle',
                                 size=9,
                                 color = cluster_colors,
                                 colorscale='Viridis',
                                 line=dict(color='rgb(50,50,50)', width=0.5)
                                 ),
                        text = ['cluster: ' + str(i) for i in range(0, len(cluster_colors))],
                   )
    #Joining the traces and plotting
    traces.append(trace2)
    data = traces

    fig=go.Figure(data=data)
    fig.update_layout(showlegend=False)
    
    if filename is not None:
        plotly.offline.plot(fig, filename=filename, auto_open=False)
    else:
        fig.show()


############ UNTESTED CODE!!!!!#################

import scipy
import numpy as np
import pandas as pd
import scanpy as sc
import anndata

def cell_kde_pergroup(adata,
                     groupby = 'condition',
                     coords = 'X_umap',
                     bw_method = None,
                     estimate_cell_densities = True):
    from scipy import stats
    if pd.api.types.is_categorical_dtype(adata.obs[groupby]):
        cats = adata.obs[groupby].cat.categories.values
    else:
        cats = adata.obs[groupby].unique()
    d = pd.DataFrame(index = adata.obs.index)
    
    for i in cats:
        x = adata[adata.obs.index[adata.obs[groupby] == i],:]
        dens = stats.gaussian_kde(x.obsm[coords].T, bw_method = bw_method)
        if estimate_cell_densities: dens = dens(adata.obsm[coords].T)
        d[i] = dens
    return(d)

def plot_kde_pergroup(adata, 
                      groupby = 'condition',
                      coord = 'X_umap',
                      zscore = True,
                      clip = (-3,3),
                      bw_method = None,
                      save = None):
    
    bdens = cell_kde_pergroup(adata, groupby = groupby,
                              bw_method = bw_method)
    
    if zscore:
        bdens.apply(scipy.stats.zscore, axis = 1)
        bdensZ = bdens.T.apply(scipy.stats.zscore, axis = 0).T
        bdensZ.columns = [x + '_densZ' for x in bdensZ.columns ]

    x = anndata.AnnData(obs = adata.obs.copy())
    x.obsm[coord] = adata.obsm[coord]
    x.obs = x.obs.merge(bdensZ, left_index = True, right_index = True)
    
    cmap2 = cmap_RdBu2(None, vmin = clip[0], vmax = clip[1])
    if zscore:
        sc.pl.umap(x, color = bdensZ.columns, cmap = cmap2, save = save, vmin = clip[0], vmax = clip[1])
    else:
        sc.pl.umap(x, color = bdens.columns, save = save)

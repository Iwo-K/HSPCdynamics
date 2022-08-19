#!/usr/bin/env ipython
import networkx as nx
import pandas as pd
from scipy.sparse import issparse
import matplotlib.pyplot as plt
import numpy as np


def plot_arrows(adata,
                transitions=None,
                pos=None,
                width=None,
                edge_width_scale=1,
                edge_color='k',
                style='solid',
                alpha=None,
                arrowstyle='-|>',
                arrowsize=10,
                edge_cmap=None,
                edge_vmin=None,
                edge_vmax=None,
                ax=None,
                arrows=True,
                node_size=300,
                connectionstyle='arc3'):
    """
    Plot pretty arrows to a PAGA graph, allows much more customisation than
    the original Essentially a wrapper around draw_networkx_edges, which
    exposes a series of edge customisation options/arguments

    If pos is None, node positions are extracted from adata.uns['paga']['pos']
    If transitions is None, the transitions matrix is extracted from
    adata.uns['paga']['transitions']

    Note that width is set to None by default, which uses the weights
    from the transition matrix (either supplied or adata object)
    otherwise set to
    """
    groups = adata.uns['paga']['groups']
    labels = adata.obs[groups].cat.categories

    if transitions is None:
        transitions = adata.uns['paga']['transitions']
    if issparse(transitions):
        transitions = transitions.toarray()
    transitions = pd.DataFrame(transitions)
    transitions.index = labels
    transitions.columns = labels

    if pos is None:
        pos = adata.uns['paga']['pos']

    pos = {n: [p[0], p[1]] for n, p in zip(labels, pos)}

    G = nx.DiGraph(transitions.T)

    if width is None:
        width = [x[-1]['weight']*edge_width_scale for x in G.edges(data=True)]

    nx.draw_networkx_edges(G,
                           pos=pos,
                           width=width,
                           edge_color=edge_color,
                           style=style,
                           alpha=alpha,
                           arrowstyle=arrowstyle,
                           arrowsize=arrowsize,
                           edge_cmap=edge_cmap,
                           edge_vmin=edge_vmin,
                           edge_vmax=edge_vmax,
                           ax=ax,
                           arrows=arrows,
                           node_size=node_size,
                           connectionstyle=connectionstyle)


def plot_widthbar(vmin=0,
                  vmax=4,
                  steps=4,
                  edge_width_scale=1,
                  font_size=2,
                  ndigits=0,
                  ax=None):
    """
    Plot a series of edges with thickness labels - basic legend.
    """
    steps = 4
    nodes = [[i, i+1] for i in range(0, steps*2, 2)]

    g = nx.DiGraph(nodes)
    pos = {}
    length = 2  # line length
    for n, i in enumerate(nodes):
        pos[i[0]] = (0, n)
        pos[i[1]] = (length, n)

    w = np.round(np.linspace(vmin, vmax, steps), ndigits)
    width = w * edge_width_scale
    width_labels = w

    nx.draw_networkx_edges(g,
                           pos=pos,
                           width=width,
                           ax=ax)
    labels = {i: n for i, n in zip(g.edges, width_labels)}
    nx.draw_networkx_edge_labels(g,
                                 pos=pos,
                                 edge_labels=labels,
                                 font_size=font_size,
                                 horizontalalignment='center',
                                 ax=ax)
    if ax is None:
        plt.grid(False)
    else:
        ax.grid(None)

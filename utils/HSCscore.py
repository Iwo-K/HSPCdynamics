import pickle
import pandas as pd
import numpy as np
import scanpy as sc

def total_count_normalise(count_matrix):
    """Normalise count matrix for input into hscScore model.
    Performs read depth normalisation normalising each cell so that normalised 
    counts sum to the same value.
    
    Parameters
    ----------
    count_matrix : pandas dataframe
        Gene count matrix of dimension cells x genes with column names as genes
        and index as cell names
    
    Returns
    -------
    **norm_matrix** : pandas dataframe
        Normalised count matrix of dimension cells x genes
    """
    
    # Set the value normalised counts will sum to for each cell
    wilson_molo_genes_median_counts = 18704.5
    
    # Scale rows
    count_matrix_expression = np.array(count_matrix, dtype='float')
    counts_per_cell = np.sum(count_matrix_expression, axis=1)
    counts_per_cell += (counts_per_cell == 0)
    counts_per_cell /= wilson_molo_genes_median_counts
    norm_matrix_expression =  count_matrix_expression/counts_per_cell[:, None]
    norm_matrix = pd.DataFrame(norm_matrix_expression, index=count_matrix.index,
                               columns=count_matrix.columns)
    # log + 1 transform the data
    norm_matrix = np.log(norm_matrix + 1)
    
    return norm_matrix

def calculate_HSCscore(adata, modelpath='./data/HSCscore/'):
    with open(modelpath + '/hscScore_model_0.22.2post1.pkl', 'rb') as f:
        HSCscore_model = pickle.load(f)
    model_genes = np.genfromtxt(modelpath + '/model_molo_genes.txt', dtype='str')
    
    X = adata[:,model_genes].X.copy()
    X = pd.DataFrame(X)
    X = total_count_normalise(X)
    
    HSCscore = HSCscore_model.predict(np.array(X))
    return HSCscore
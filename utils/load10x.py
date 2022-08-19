import scanpy as sc
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
from .proc import do_qc, lognorm, assign_sex
import scrublet as scr


def load_10xfiles(meta, ygenes, qc = True):
    adatas = {}
    verbosity = sc.settings.verbosity
    sc.settings.verbosity = 1
    
    files = meta.countfile.unique()

    for i in files:
        metaI = meta.loc[meta.countfile == i,:]
        #' Checks on the number of samples in combination with sex_mixed
        if len(metaI.index) > 2:
            print("Error: to many samples per countfile")
            return()
        elif (metaI.sex_mixed == 'no').all() and len(metaI.index) != 1:
            print("Error: sex_mixed set to no but sample number does not equal 1")
        elif (metaI.sex_mixed == 'yes').all() and len(metaI.index) != 2:
            print("Error: sex_mixed set to yes but sample number does not equal 2")

        print(f'''################ Processing file {i} ##################''')
        adata = sc.read_10x_h5(i)
        #Setting unique cellnames
        adata.obs.index = metaI.sample_id.values[0] + "_" + adata.obs.index
        
        #Setting unique gene names
        adata.var['symbol'] = adata.var.index
        adata.var_names_make_unique()
        if qc: do_qc(adata, min_genes = 1000, mt_frac = 0.075)

        #Estimating doublets
        scrub = scr.Scrublet(adata.X.todense(), expected_doublet_rate=0.1)
        doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
        scrub.plot_histogram()
        adata.obs['doublet_scores'] = doublet_scores.copy()
        adata.obs['predicted_doublets'] = predicted_doublets.copy()

        lognorm(adata)

        #Not all samples are a mixture of M+F samples so checking and splitting the sample by sex
        if (metaI.sex_mixed == 'yes').all():
            #Subsetting for Y chromosome genes existing in dataset to avoid NAs
            ygenes = ygenes[ygenes.isin(adata.var.symbol)]
            print("Results of sex assignment:")
            assign_sex(adata, ygenes = ygenes, use_raw = False)
            #' Filtering out the unassigned
            adata = adata[adata.obs.index[adata.obs.sex.isin(('male', 'female'))],:]

            adata_female = adata[adata.obs.index[adata.obs.sex == 'female'],:]
            adata_male = adata[adata.obs.index[adata.obs.sex == 'male'],:]

            female_id = metaI.loc[metaI.sex == 'female', 'biosample_id'].values[0]
            male_id = metaI.loc[metaI.sex == 'male', 'biosample_id'].values[0]
            adatas[female_id] = adata_female
            adatas[male_id] = adata_male

            print(f'''Returning:\n
            - female adata: {female_id} with {adata_female.n_obs} observations\n
            - male adata: {male_id} with {adata_male.n_obs} observations''')

        elif (metaI.sex_mixed == 'no').all():
            adatas[metaI.biosample_id.values[0]] = adata
            adata.obs['sex'] = metaI.sex.values[0]
            print(f'''Returning:\n
            - adata: {metaI.biosample_id.values[0]} with {adata.n_obs} observations''')

        print(f'''###################################################################''')

    sc.settings.verbosity = verbosity
    return(adatas)

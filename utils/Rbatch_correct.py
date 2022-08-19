import subprocess
import scanpy as sc
import pathlib
import os
import tempfile


TRANSFER_DIR = '.anndata_transfers_'


def temp_save_adata(adata):
    ''' Create a temporary directory and save a provided adata.
    Returns a tuple with the TemporaryDirectory, name of the saved file
    and a name for the future processed file
    '''

    x = tempfile.TemporaryDirectory(prefix = TRANSFER_DIR, dir = './')
    infile = x.name + '/source_adata.h5ad'
    outfile = x.name + '/corrected_adata.h5ad'
    adata.write(infile, compression = 'lzf')

    # Need to return the TemporaryDirectory otherwise it leaves the
    # scope and directory gets deleted
    return (x, infile, outfile)

def run_SeuratCCA(adata,
                  batch_key='data_type',
                  reference='10x',
                  anchor_dims=30,
                  integrate_dims=30,
                  k_anchor=5,
                  k_filter=200,
                  k_score=30,
                  debug=False):
    ''' Python wrapper around Routine run_SeuratCCA, which
    performs batch correction on scRNA-Seq data.
    1. Transfers the AnnData object via disk to R (temporary directory)
    2. Runs the R function run_SeuratCCA (passing the arguments)
    by creating a temporary script
    3. Reads the run_SeuratCCA output file and returns to the user

    If debug=True, prints stdout and stderr from the R subprocess
    and returns a tuple of corrected adata and the
    TemporaryDirectory object to prevent the directory cleanup.
    '''


    d, infile, outfile = temp_save_adata(adata)
    print(f'Storing data in "{d.name}" temporary directory')

    script = f'''
source("./utils/Routines.R")
run_SeuratCCA(infile = "{infile}",
              anchor.dims = 1:{anchor_dims},
              integrate.dims = 1:{integrate_dims},
              batch_key = "{batch_key}",
              reference = "{reference}",
              k.anchor = {k_anchor},
              k.filter = {k_filter},
              k.score = {k_score},
              outfile = "{outfile}")
    '''
    script_file = d.name + '/run_seuratCCA.R'
    with open(script_file, 'w') as f:
        f.write(script)

    out = subprocess.run(f'R --no-save < {script_file}',
                         capture_output=True,
                         text=True,
                         shell=True)

    corrected = sc.read(outfile)

    if debug:
        print('out:' + out.stdout)
        print('err:' + out.stderr)
        return corrected,d
    else:
        out.check_returncode()
        return corrected

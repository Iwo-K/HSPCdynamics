library(anndata)
library(batchelor)
library(scater)
library(Seurat)

anndata_to_sce = function(file, subset_to_var = TRUE){

  data = read_h5ad(file)
  sce <- SingleCellExperiment(assays = List(logcounts = as(t(data$X), "CsparseMatrix")),
                              colData=data$obs,
                              rowData=data$var,
                              reducedDims = data$obsm)
  if (subset_to_var){
    scesub = sce[rowData(sce)$highly_variable,]
  }
  return(sce)
}

anndata_to_Seurat = function(file){

  sce = anndata_to_sce(file)
  seu = as.Seurat(sce, counts = NULL, data = 'logcounts')
  return(seu)
}

run_SeuratCCA = function(infile,
                         anchor.dims = 1:20,
                         integrate.dims = 1:30,
                         batch_key = 'data_type',
                         reference = '10x',
                         k.anchor = 5,
                         k.filter = 200,
                         k.score = 30,
                         outfile = 'corrected_adata.h5ad'){
  print('Running SeuratCCA')

  seu = anndata_to_Seurat(infile)

  seu_list = SplitObject(seu, split.by = batch_key)
  reference = which(names(seu_list) == reference)
  print(paste0('Reference data: ', names(seu_list)[reference]))
  print(reference)

  anchors = FindIntegrationAnchors(object.list = seu_list,
                                   dims = anchor.dims,
                                   reference = reference,
                                   anchor.features = row.names(seu),
                                   k.anchor = k.anchor,
                                   k.filter = k.filter,
                                   k.score = k.score)
  print(anchors)

  combined = IntegrateData(anchorset = anchors,
                           dims = integrate.dims)
  print(combined)

  #Making sure cell/gene orders are the same
  combined = combined[row.names(seu), colnames(seu)]
  #Convertin to AnnData and saving
  x = combined@assays$integrated@data
  x = AnnData(X = t(x),
              obs = combined@meta.data,
              var = seu@assays$originalexp@meta.features)
  x$write_h5ad(outfile, compression = 'lzf')
}



## run_fastMNN = function(file1,
##                        file2,
##                        d = 50,
##                        ndist = 10,
##                        k = 15,
##                        cos.norm = FALSE){
##   sce1 = anndata_to_sce(file1)
  ## sce2 = anndata_to_sce(file2)

  ## set.seed(101)
  ## integrated <- fastMNN(A=sce1,
  ##                B=sce1,
  ##                d = d,
  ##                ndist = ndist,
  ##                k = k,
  ##                cos.norm = cos.norm)
  ## corPCA = reducedDim(integrated, 'corrected')
  ## corPCA = as.data.frame(corPCA)
  ## return(corPCA)
## 

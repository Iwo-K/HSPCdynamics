#' # scB5Tom experiment SS2 data - 3d7d, 3dr2, 7d and 2w batches
#' Author: Iwo Kucinski
#' Last updated: `r Sys.time()`
#+ setup, include=FALSE
library(devtools)
library(biomaRt)
library(reshape2)
library(ggrepel)
library(viridis)
library(anndata)
library(ggplot2)
theme_set(theme_bw())
source("./utils/loadSS2.R")
source("./utils/scfuns.R")

#' Output directories
procdata_dir = './procdata/01script/'
dir.create(procdata_dir, showWarnings = FALSE)
figures_dir = './figures/01script/'
dir.create(figures_dir, showWarnings = FALSE)

#' ## Loading data
meta = read.csv("./data/SS2/scB5Tom_SS2_cellmeta.csv", as.is = TRUE)

load_ss2data = function(meta){

  data = list()
  for (i in unique(meta$countfolder)){
    print(i)
    data[[i]] = get_counts(i)
  }

  cellids = lapply(data, colnames)
  cellids = unlist(cellids)
  if (sum(duplicated(cellids)) > 0) stop('WARNING some of the names are not unique!')

  names(data) = c()
  data = do.call(cbind, data)
  return(data)

}

data = load_ss2data(meta)

#' ### Ensembl genes
## genedata = makeGENEDATA(release = "jul2018.archive.ensembl.org")
## fwrite(genedata, './data/genedata_v93.csv')
genedata = fread('./data/genedata_v93.csv')
genedata = as.data.frame(genedata)
row.names(genedata) = genedata$geneid

#' ### Removing wells which there were some problems while sorting
## meta = meta[meta$touse == 'TRUE',]
data = data[, as.character(meta$cellid)]

meta$batch_plate_sorted = paste0(meta$batch, "_", meta$plate_sorted)

#' ## QC
dataQC = scQC(data,
              metadata = meta,
              splitby = "batch_plate_sorted",
              mitogenes = genedata[genedata$is_mito,"geneid"],
              id = "cellid",
              tr_mapped = 0.23,
              tr_mincounts=1e+05,
              tr_maxcounts=3e+07,
              tr_erccfrac=0.085,
              tr_higeneno = 2000,
              tr_mitofrac = 0.1,
              filebase = paste0(procdata_dir, "QCperplate"))
metaQC = meta[meta$cellid %in% colnames(dataQC), ]
dataQC = dataQC[,metaQC$cellid]

#' Selectign variable genes
vargenes = hivar(dataQC,
                 meanquantile = 0.4,
                 cv2tr = 0.1,
                 padjtr = 0.1,
                 plotfile = paste0(figures_dir, "hivar.pdf"))

#' ## Normalising data
#' Extracting only genes in the counttable
dataQC = dataQC[grepl("^ENSM*", row.names(dataQC)),]

dataN = scnormalise_v2(dataQC, norm_method = "scran")
logdata = log2(dataN + 1)
logdatavar = logdata[row.names(logdata) %in% row.names(vargenes),]

#' ### PCA - per set
pca12 = PCAproj(logdatavar, id = 'cellid', metaQC, sets = c("batch_plate_sorted", "batch", "mouse_id", "sort_method", "timepoint_tx_days"), pdffile = paste0(figures_dir, "pca_sets12.pdf"))

tsne12 = TSNEproj(logdatavar, id ='cellid', metaQC, sets = c("batch_plate_sorted", "batch", "mouse_id", "sort_method", "timepoint_tx_days"), pdffile = paste0(figures_dir, "tsne_sets12.pdf"))

pdf(paste0(figures_dir, "pca_chosengenes.pdf"))
for (gene in c("Procr", "Slamf1", "Klf1", "Pf4", "Hoxb5", "Cd48", "Elane", "Irf8", "Dntt")){
  geneid = genedata[genedata$symbol == gene,"geneid"]
  expr = unlist(logdata[geneid,])
  g = ggplot(pca12, aes(x = PC1, y = PC2, color = expr)) + geom_point() + scale_color_viridis() + ggtitle(gene)
  plot(g)
}
dev.off()

#' ## Saving an AnnData file
genedata = genedata[row.names(dataQC),]
adata = AnnData(X = t(dataQC), obs = metaQC, var = genedata)
adata$write_h5ad(paste0(procdata_dir, 'SS2data.h5ad'), compression = 'lzf')

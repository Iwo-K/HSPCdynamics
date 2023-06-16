# ---
# jupyter:
#   jupytext:
#     formats: ipynb,R:percent
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.11.4
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# %% tags=[]
library(tradeSeq)
library(anndata)
library(SingleCellExperiment)
library(data.table)
library(ggplot2)
library(scater)
library(RColorBrewer)
source('utils/DE.R')
source('utils/scfuns.R')
library(plotly)
library(pheatmap)
library(data.table)

BPPARAM <- BiocParallel::bpparam()
BPPARAM$workers <- 4

#ggplot theme
theme_paper = theme_bw(base_size = 7, base_family = "",
                       base_line_size = 0.5, base_rect_size = 0.5)
theme_set(theme_paper)

DIR_DE = './DE/'


# %% tags=[]
adata_to_sce = function(adata){
    require(SingleCellExperiment)
    sce = SingleCellExperiment(assays = List(counts = t(adata$X)),
                               colData = adata$obs,
                               rowData = adata$var,
                               reducedDims = adata$obsm)
}

# %% tags=[]
data  = read_h5ad('./procdata/04script/combined_filt_counts.h5ad')

# %% tags=[]
#Getting transcription factors
allTFs = fread('./data/genesets/TFs_Ravasi2010short.csv')
allTFs = unlist(allTFs[,'Symbol (Mouse)'])

# %% [markdown]
# ## Megakaryocyte trajectory

# %% tags=[]
mkdir = paste0(DIR_DE, 'Mk/')
dir.create(mkdir, recursive=TRUE)

mkpar = fread('./PD_model/clu_7/tables/table_all_parameters_clu_7.csv')
mk = data[mkpar$V1, ]

# %% tags=[]
mk = adata_to_sce(mk)
colData(mk) = cbind(colData(mk), data.frame(drift = mkpar$drift,
                                            growth = mkpar$growth,
                                           dpt_pseudotime = mkpar$dpt_pseudotime))

plotReducedDim(mk, dimred = "X_umap_2d", colour_by = "drift")

# %% tags=[]
plot_cellparams = function(sce, drift_range, growth_range, ratio, save='test.pdf'){
    celldens = ggplot(data.frame(colData(sce)), aes(x=dpt_pseudotime)) + 
      geom_density() 

    drift = ggplot(data.frame(colData(sce)), aes(x = dpt_pseudotime, y = drift)) +
      geom_line() +
      geom_vline(xintercept = drift_range, alpha = 0.5, color='blue')


    growth = ggplot(data.frame(colData(sce)), aes(x = dpt_pseudotime, y = growth)) +
      geom_line() +
      geom_vline(xintercept = growth_range, alpha = 0.5, color='red')

    comb = ggplot(data.frame(colData(sce)), aes(x = dpt_pseudotime)) +
      geom_line(aes(x = dpt_pseudotime, y = growth), color = 'red', size = 0.8) +
      geom_line(aes(x = dpt_pseudotime, y = drift*ratio), color = 'blue', size = 0.8) +
      scale_y_continuous(name = 'Net prolif.', sec.axis = sec_axis(~.*1/ratio, name="Diff. rate")) + 
      geom_vline(xintercept = dptrange_drift, alpha = 0.5, color='blue', size=0.5, linetype='dashed') + 
      geom_vline(xintercept = dptrange_growth, alpha = 0.5, color='red', size=0.5, linetype='dashed') +
      theme(axis.text.y.right = element_text(color = "blue"),
        axis.title.y.right = element_text(color = "blue"),
        axis.text.y.left = element_text(color = "red"),
        axis.title.y.left = element_text(color = "red"))

    pdf(save, width=8/2.54, height=12/2.54)
    multiplot(celldens, drift, growth, comb)
    dev.off()
    multiplot(celldens, drift, growth, comb)
    }


# %% tags=[]
dptrange_drift = c(0.0125, 0.0375)
dptrange_growth = c(0.01, 0.028)
plot_cellparams(mk,
                drift_range = dptrange_drift,
                growth_range = dptrange_growth,
                ratio = 10,
                save=paste0(mkdir, 'mk_dpt_plots.pdf'))

# %%
#Excluding SS2 samples to not affect the DE for early DPT
mk = mk[colData(mk)$data_type != 'SS2',]
mk = logNormCounts(mk)

# %%
tokeep = gene_filter(assays(mk)$logcounts, min_expr = 0, min_fraction = 0.025, min_meanexpr = 0.05)
genes = row.names(mk)[tokeep]

demk1 = run_transitionDE(mk,
                         dpt_col = 'dpt_pseudotime',
                         dpt_range = dptrange_drift,
#                   genes = c('Actb', 'Pf4', 'Procr', 'Gata1'),
                         genes = genes,
                         allTFs = allTFs,
                         traj_name = 'mk_drift',
                         nknots = 7,
                         dir = mkdir,
                         BPPARAM = BPPARAM)

# %%

# %%
#calculating derivative of the drift
diff_drift = function(x, dpt, drift, plot=FALSE){
    af = approxfun(dpt, drift)


    y = af(x)

    #getting the midpoints
    mids = x[1:(length(x)-1)] + (x[2:length(x)] - x[1:(length(x)-1)])/2
    #getting derivatives
    dy = diff(y)

    #calculating the derivative for each cell
    daf = approxfun(mids, dy)

    if(plot){
        g = ggplot(data.frame(), aes(x=x, y=y)) + 
        geom_point(size=3) + #black original points
        geom_point(data=data.frame(), aes(x=mids, y=dy), color='red', size=3) + #red mid points by derivative
        geom_point(data=data.frame(), aes(x=x, daf(x)), color='green', size=3) #green interpolated points in the same spots
        print(g)
    }
    return(daf(x))
    }


# Runs an associationTest to find dynamic genes and then tests their correlation with the differentiation rates
analyse_cor = function(sim, pseudotime, drift, TFs, prefix){
                      
    #First running an association test to seleting genes with varying along pseudotime with l2fc of at least1
    de = associationTest(sim, l2fc=1, contrastType='end')
    de$fdr = p.adjust(de$pvalue, method='fdr')
    de = de[complete.cases(de),]
    desig = de[de$fdr < 0.01,]
            
    #predicting smoothed expression for each de gene
    yhat = predictSmooth(sim, gene = row.names(desig), nPoints = 100, tidy = TRUE)
    yhat = data.table(yhat)[, yhat_scaled := scale(yhat), by='gene']
    yhatw = dcast(yhat, time~gene, value.var='yhat')
            
            
    ddrift = diff_drift(x = yhatw$time, pseudotime, drift, plot=TRUE)
    keept = !is.na(ddrift)
    
    #calculating correlations
    ddrift_cors = apply(yhatw[,colnames(yhatw) != "time", with=FALSE], 2, function(x) cor(x[keept], ddrift[keept]))
    ddrift_cors = ddrift_cors[order(ddrift_cors, decreasing=TRUE)]
                    
    toptf = head(ddrift_cors[names(ddrift_cors) %in% TFs], n = 10)    
    print(toptf)                    
    bottomtf = tail(ddrift_cors[names(ddrift_cors) %in% TFs], n = 10)    
    print(bottomtf)
                               
    default_20 = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61',
               '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2',
               '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31')

    gup = ggplot(yhat[gene %in% names(toptf),], aes(x=time, y=yhat_scaled, color=gene)) + 
        geom_line(size=0.8) +
        geom_line(data= data.frame(), aes(x = yhatw$time, y = ddrift*250), , color='grey', linetype='dashed', size = 0.8) +
        scale_y_continuous(name = 'Scaled expression', sec.axis = sec_axis(~.*250, name="Diff. rate derivative")) +
        theme(axis.text.y.right = element_text(color = "grey45"),
            axis.title.y.right = element_text(color = "grey45"),
            axis.text.y.left = element_text(color = "black"),
            axis.title.y.left = element_text(color = "black")) +
           scale_color_manual(values=default_20) +
           theme_paper +
           xlab('Pseudotime') +
           ggtitle('Top positively correlated TFs') 

    gdown = ggplot(yhat[gene %in% names(bottomtf),], aes(x=time, y=yhat_scaled, color=gene)) + 
        geom_line(size=0.8) +
        geom_line(data= data.frame(), aes(x = yhatw$time, y = ddrift*250), , color='grey', linetype='dashed', size = 0.8) +
        scale_y_continuous(name = 'Scaled expression', sec.axis = sec_axis(~.*250, name="Diff. rate derivative")) +
        theme(axis.text.y.right = element_text(color = "grey45"),
            axis.title.y.right = element_text(color = "grey45"),
            axis.text.y.left = element_text(color = "black"),
            axis.title.y.left = element_text(color = "black")) +
           scale_color_manual(values=default_20) +
           theme_paper +
           xlab('Pseudotime') +
           ggtitle('Top negatively correlated TFs')

    corplots = cowplot::plot_grid(plotlist = list(gup, gdown), align='h', ncol=1)
    ggsave(paste0(prefix, '_topcorplots.pdf'), corplots, height= 6, width=5)
}
                        
analyse_cor(demk1$sim, pseudotime=mkpar$dpt_pseudotime, drift=mkpar$drift, TFs=allTFs, prefix='test_')

# %%

# %%

# %% [markdown]
# ## Analysing correlation with drift (function development)

# %% tags=[]
#First running an association test to seleting genes with varying along pseudotime with l2fc of at least1
de = associationTest(demk1$sim, l2fc=1, contrastType='end')
de$fdr = p.adjust(de$pvalue, method='fdr')

de = de[complete.cases(de),]
desig = de[de$fdr < 0.01,]

# %%
#predicting smoothed expression for each de gene
yhat = predictSmooth(demk1$sim, gene = row.names(desig), nPoints = 100, tidy = TRUE)
yhat = data.table(yhat)[, yhat_scaled := scale(yhat), by='gene']
yhatw = dcast(yhat, time~gene, value.var='yhat')

# %%
#calculating derivative of the drift
diff_drift = function(x, dpt, drift, plot=FALSE){
    af = approxfun(dpt, drift)


    y = af(x)

    #getting the midpoints
    mids = x[1:(length(x)-1)] + (x[2:length(x)] - x[1:(length(x)-1)])/2
    #getting derivatives
    dy = diff(y)

    #calculating the derivative for each cell
    daf = approxfun(mids, dy)

    if(plot){
        g = ggplot(data.frame(), aes(x=x, y=y)) + 
        geom_point(size=3) + #black original points
        geom_point(data=data.frame(), aes(x=mids, y=dy), color='red', size=3) + #red mid points by derivative
        geom_point(data=data.frame(), aes(x=x, daf(x)), color='green', size=3) #green interpolated points in the same spots
        print(g)
    }
    return(daf(x))
    }

ddrift = diff_drift(x = yhatw$time, mkpar$dpt_pseudotime, mkpar$drift, plot=TRUE)
keept = !is.na(ddrift)

# %%
#calculating correlations
ddrift_cors = apply(yhatw[,colnames(yhatw) != "time", with=FALSE], 2, function(x) cor(x[keept], ddrift[keept]))
ddrift_cors = ddrift_cors[order(ddrift_cors, decreasing=TRUE)]
                    
toptf = head(ddrift_cors[names(ddrift_cors) %in% allTFs], n = 10)    
print(toptf)                    
bottomtf = tail(ddrift_cors[names(ddrift_cors) %in% allTFs], n = 10)    
print(bottomtf)

# %%
head(yhat)

# %%
default_20 = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61',
               '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2',
               '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31')

gup = ggplot(yhat[gene %in% names(toptf),], aes(x=time, y=yhat_scaled, color=gene)) + 
    geom_line(size=0.8) +
    geom_line(data= data.frame(), aes(x = yhatw$time, y = ddrift*250), , color='grey', linetype='dashed', size = 0.8) +
    scale_y_continuous(sec.axis = sec_axis(~.*250, name="Diff. rate derivative")) +
    theme(axis.text.y.right = element_text(color = "grey45"),
        axis.title.y.right = element_text(color = "grey45"),
        axis.text.y.left = element_text(color = "black"),
        axis.title.y.left = element_text(color = "black")) +
       scale_color_manual(values=default_20) +
       theme_paper

gdown = ggplot(yhat[gene %in% names(bottomtf),], aes(x=time, y=yhat_scaled, color=gene)) + 
    geom_line(size=0.8) +
    geom_line(data= data.frame(), aes(x = yhatw$time, y = ddrift*250), , color='grey', linetype='dashed', size = 0.8) +
    scale_y_continuous(sec.axis = sec_axis(~.*250, name="Diff. rate derivative")) +
    theme(axis.text.y.right = element_text(color = "grey45"),
        axis.title.y.right = element_text(color = "grey45"),
        axis.text.y.left = element_text(color = "black"),
        axis.title.y.left = element_text(color = "black")) +
       theme_paper

corplots = cowplot::plot_grid(plotlist = list(gup, gdown), align='h', ncol=1)
ggsave('test.pdf', corplots, height= 6, width=5)

# %%

# %%

# %%

# %%

# %%

# %%

# %%
# ccf(yhatw$Gata1[!is.na(ddrift)], ddrift[!is.na(ddrift)])

# %%

# %%
af = approxfun(mkpar$dpt_pseudotime, mkpar$drift)

#Calculating derivatives with 2x linear interpolation

#equally spaced points interpolated
x = yhatw$time[1:40]
y = af(x)
# plot(x, y)


#getting the midpoints
mids = x[1:(length(x)-1)] + (x[2:length(x)] - x[1:(length(x)-1)])/2
#getting derivatives
dy = diff(y)

# plot(mids, dy)

#calculating the derivative for each cell
daf = approxfun(mids, dy)

# plot(x, daf(x))

#Because of differentiation and interpolation we lose the most left and most right point

ggplot(data.frame(), aes(x=x, y=y)) + 
geom_point(size=3) + #black original points
geom_point(data=data.frame(), aes(x=mids, y=dy), color='red', size=3) + #red mid points by derivative
geom_point(data=data.frame(), aes(x=x, daf(x)), color='green', size=3) #green interpolated points in the same spots

# %%

# %%

# %%

# %%

# %%

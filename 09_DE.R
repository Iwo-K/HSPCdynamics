# ---
# jupyter:
#   jupytext:
#     formats: ipynb,R:percent
#     text_representation:
#       extension: .R
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.12.0
#   kernelspec:
#     display_name: R
#     language: R
#     name: ir
# ---

# %% [markdown]
# # Differential expression analysis for trajectories
#
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
# Some handy functions
adata_to_sce = function(adata){
    require(SingleCellExperiment)
    sce = SingleCellExperiment(assays = List(counts = t(adata$X)),
                               colData = adata$obs,
                               rowData = adata$var,
                               reducedDims = adata$obsm)
}

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

# %% tags=[]
tokeep = gene_filter(assays(mk)$logcounts, min_expr = 0, min_fraction = 0.025, min_meanexpr = 0.05)
genes = row.names(mk)[tokeep]

demk1 = run_transitionDE(mk,
                         dpt_col = 'dpt_pseudotime',
                         dpt_range = dptrange_drift,
                         genes = genes,
                         allTFs = allTFs,
                         traj_name = 'mk_drift',
                         nknots = 7,
                         dir = mkdir,
                         BPPARAM = BPPARAM)

# %% tags=[]
tokeep = gene_filter(assays(mk)$logcounts, min_expr = 0, min_fraction = 0.025, min_meanexpr = 0.05)
genes = row.names(mk)[tokeep]
demk2 = run_transitionDE(mk,
                  dpt_col = 'dpt_pseudotime',
                  dpt_range = dptrange_growth,
                  genes = genes,
                  allTFs = allTFs,
                  traj_name = 'mk_growth',
                 nknots = 7,
                  dir = mkdir,
                  BPPARAM = BPPARAM)

# %% [markdown]
# ### Additional plots for the growth DE

# %%
# Running enrichment with enrichr
gen = run_enrichr(row.names(demk2$desig),
            gene_sets=c('MSigDB_Hallmark_2020', 'Reactome_2016'),
            figure_prefix=paste0(mkdir, '/mk_growth_'))
#Extracting genes from categories of interest

e2f = gen[gen$Gene_set == 'MSigDB_Hallmark_2020' & gen$Term == 'E2F Targets', 'Genes']
e2f = str2mouse(e2f)

g2mc = gen[gen$Gene_set == 'MSigDB_Hallmark_2020' & gen$Term == 'G2-M Checkpoint', 'Genes']
g2mc = str2mouse(g2mc)

cc = gen[gen$Gene_set == 'Reactome_2016' & gen$Term == 'Cell Cycle Homo sapiens R-HSA-1640170', 'Genes']
cc = str2mouse(cc)

# %%
#Preparing data for plotting
npoints=100

yhat = predictSmooth(demk2$sim, gene = row.names(demk2$desig), nPoints = 100, tidy = TRUE)
yhatw = dcast(yhat, time~gene, value.var='yhat')
xtime = yhatw$time

yhatw = as.matrix(yhatw[,-1])
row.names(yhatw) = xtime
yhatw = t(yhatw)
yhatw <- t(scale(t(yhatw)))

dpt_anno = rep('out', npoints)
dpt_anno[(xtime > dptrange_growth[1] & xtime < dptrange_growth[2])] = 'in'
annotation_col = data.frame(region = dpt_anno, row.names=colnames(yhatw))

annotation_row = data.frame(row.names = row.names(yhatw),
                            E2F = ifelse(row.names(yhatw) %in% e2f, 'E2F targets', 'other'),
                            G2MC = ifelse(row.names(yhatw) %in% g2mc, 'G2-M checkpoint', 'other'),
                            CC = ifelse(row.names(yhatw) %in% cc, 'Cell cycle', 'other'))
ann_colors = list(
    region = c(out="#00BFC4", `in`="#F8766D"),
    E2F = c(`E2F targets` ='blue', other='lightgrey'),
    G2MC = c(`G2-M checkpoint` ='purple', other='lightgrey'),
    CC = c(`Cell cycle` ='darkorange', other='lightgrey')
    )

# %%
p1 = pheatmap(yhatw,
              cluster_cols=FALSE,
              annotation_col=annotation_col,
              annotation_row=annotation_row,
              annotation_colors = ann_colors,
              fontsize_row=4.5,
              treeheight_row=0,
              show_colnames=FALSE)
ggsave(paste0(mkdir, '/mk_growth_heatmap_annotated.pdf'), p1, height= 20, width=14)

# %%
toplot = yhat[yhat$gene %in% c('Pf4', 'Procr',
                             'Mcm5', 'Mki67', 'Hist2h2ac', 'E2f8', 'Lmnb1'),]
p2 = ggplot(toplot,
           aes(x=time, y=log2(yhat+1), color=gene)) +
           geom_line() +
           xlab('Pseudotime') + 
           geom_vline(xintercept = dptrange_growth, linetype='dashed', color = 'red')
ggsave(paste0(mkdir, '/mk_growth_example_genes.pdf'), p2, height=2, width=4)
p2

# %% [markdown]
# ## Neutrophil trajectory

# %% tags=[]
neudir = paste0(DIR_DE, 'Neu/')
dir.create(neudir, recursive=TRUE)

neupar = fread('./PD_model/clu_10/tables/table_all_parameters_clu_10.csv')
neu = data[neupar$V1, ]

# %% tags=[]
neu = adata_to_sce(neu)
colData(neu) = cbind(colData(neu), data.frame(drift = neupar$drift,
                                            growth = neupar$growth,
                                           dpt_pseudotime = neupar$dpt_pseudotime))

plotReducedDim(neu, dimred = "X_umap_2d", colour_by = "drift")

# %% tags=[]
dptrange_drift = c(0.065, 0.11)
dptrange_growth = c(0.01, 0.028)
plot_cellparams(neu,
                drift_range = dptrange_drift,
                growth_range = dptrange_growth,
                ratio = 50,
                save=paste0(neudir, 'neu_dpt_plots.pdf'))

# %%
#Excluding SS2 samples to not affect the DE for early DPT
neu = neu[colData(neu)$data_type != 'SS2',]
neu = logNormCounts(neu)

# %% [markdown]
# ### Neutrophil trajectory drift DE analysis

# %% tags=[]
tokeep = gene_filter(assays(neu)$logcounts, min_expr = 0, min_fraction = 0.025, min_meanexpr = 0.05)
genes = row.names(neu)[tokeep]

deneu1 = run_transitionDE(neu,
                 dpt_col = 'dpt_pseudotime',
                 dpt_range = dptrange_drift,
                 genes = genes,
                 allTFs = allTFs,
                 nknots = 8,
                 traj_name = 'neu_drift',
                 dir = neudir,
                 BPPARAM = BPPARAM)

# %% [markdown]
# #### Additional figures

# %% tags=[]
#Preparing data for plotting
npoints=100

#Long version of data unscaled
y = predictSmooth(deneu1$sim, gene = row.names(deneu1$desig), nPoints = 100, tidy = TRUE)
y = data.table(y)[, yhat_scaled := scale(yhat), by='gene']
yw = dcast(y, time~gene, value.var='yhat')
xtime = yw$time
yw = yw[,-1]

#Wide data scaled
yw_scaled = as.matrix(yw)
row.names(yw_scaled) = xtime
yw_scaled = t(yw_scaled)
yw_scaled <- t(scale(t(yw_scaled)))

#Annotation for heatmaps
dpt_anno = rep('out', npoints)
dpt_anno[(xtime > dptrange_drift[1] & xtime < dptrange_drift[2])] = 'in'
annotation_col = data.frame(region = dpt_anno, row.names=colnames(yw_scaled))

# %% [markdown]
# ##### Plotting transcription factors

# %% tags=[]
TFs = deneu1$desig[row.names(deneu1$desig) %in% allTFs,]
TFs = row.names(TFs)

#Plotting heatmaps
# First getting the hierarchical clustering (default settings do a good job)
tfheat = pheatmap(yw_scaled[row.names(yw_scaled) %in% TFs,],
              cluster_cols=FALSE,
              annotation_col=annotation_col,
              fontsize_row=4.5,
              show_colnames=FALSE)

cl = cutree(tfheat$tree_row, 4)
cl_str = paste0('TF group ', cl)
annotation_row = data.frame(row.names = names(cl), TF_group = cl_str)

tfheat = pheatmap(yw_scaled[row.names(yw_scaled) %in% TFs,],
              cluster_cols=FALSE,
              annotation_col=annotation_col,
              annotation_row=annotation_row,
              fontsize_row=4.5,
              show_colnames=FALSE)
ggsave(paste0(neudir, '/neu_drift_TFheatmap.pdf'), tfheat, height= 20, width=14)

# %% tags=[]
#Plotting as lines
toplot = y[y$gene %in% TFs,]
toplot$cluster = cl_str[match(toplot$gene, names(cl))]


default_20 = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61',
               '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2',
               '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31')

plotlist = list()
for (i in unique(toplot$cluster)){
    plotlist[[i]] = ggplot(toplot[toplot$cluster == i,],
       aes(x=time, y=yhat_scaled, color=gene)) +
       geom_line(data = as.data.frame(colData(neu)), aes(x = dpt_pseudotime, y=drift*200), color='grey', linetype='dashed') +
       geom_line(alpha=0.65) +
       xlab('Pseudotime') + 
       theme_paper +
       theme(legend.key.height=unit(0.35, 'cm')) +
       ggtitle(i) +
       scale_color_manual(values=default_20) +
       scale_y_continuous(name = 'Scaled expression', sec.axis = sec_axis(~.*1/200, name="Diff. rate"))
    }
tfplots = cowplot::plot_grid(plotlist = plotlist, align='v')
ggsave(paste0(neudir, '/neu_drift_TFplots.pdf'), tfplots, height= 4, width=7)
tfplots

# %% [markdown]
# ##### Genes associated with neutrophil and monocyte lineages

# %%
genes_toplot = list(Neu_genes = c('Elane', 'Cst7', 'Cebpe', 'Fcgr3', 'Prtn3', 'Gfi1', #Early neutrophil genes
             'S100a8', 'Clec4a2', 'Wfdc21', 'G0S2'), #Late neutrohil genes
             Transition_genes = c('Irf8', 'Flt3', 'Cd34', 'Ikzf2', 'Cxcr2', 'Slc22a3', 'Ctsh'), #handchosen transient genes from the heatmap
             Lineage_balance = c('Gfi1', 'Irf8', 'Flt3')) #LeeGrimes paper poses these two factors as opposing forces and a binary outcome in fate choice

# %% tags=[]
plotlist2 = list()
for (i in names(genes_toplot)){
    plotlist2[[i]] = ggplot(y[y$gene %in% genes_toplot[[i]],],
       aes(x=time, y=yhat_scaled, color=gene)) +
       geom_line(data = as.data.frame(colData(neu)), aes(x = dpt_pseudotime, y=drift*200), color='grey', linetype='dashed') +
       geom_line(alpha=0.65) +
       xlab('Pseudotime') +
       ggtitle(i) +
       theme_paper +
       theme(legend.key.height=unit(0.35, 'cm')) +
       scale_color_manual(values=default_20) +
       scale_y_continuous(name = 'Scaled expression', sec.axis = sec_axis(~.*1/200, name="Diff. rate"))
    }
neuplots = cowplot::plot_grid(plotlist = plotlist2, align='v')
neuplots
ggsave(paste0(neudir, '/neu_drift_chosengenes.pdf'), neuplots, height= 4, width=7)

# %%
# re-plotting to match dimensions in the figure
blankaxis = theme(axis.text.x = element_blank(),
                              axis.ticks.x = element_blank(),
                              axis.title.x = element_blank())
neu_combplot = cowplot::plot_grid(plotlist = list(plotlist2[[1]] + blankaxis,
                                                  plotlist[[1]] + blankaxis,
                                                  plotlist[[2]]), ncol = 1, align = 'v')
ggsave(paste0(neudir, '/neu_drift_neumarkers_TF12.pdf'), neu_combplot, height= 4.82, width=3.5)


# %% [markdown]
# ### Neutrophil trajectory growth DE analysis

# %%
deneu2 = run_transitionDE(neu,
                 dpt_col = 'dpt_pseudotime',
                 dpt_range = dptrange_growth,
                 genes = genes,
                 allTFs = allTFs,         
                 nknots = 8,
                 traj_name = 'neu_growth',
                 dir = neudir,
                 BPPARAM = BPPARAM)

# %% [markdown]
# ## Ery trajectory

# %% tags=[]
erydir = paste0(DIR_DE, 'Ery/')
dir.create(erydir, recursive=TRUE)

erypar = fread('./PD_model/clu_11/tables/table_all_parameters_clu_11.csv')
ery = data[erypar$V1, ]

# %% tags=[]
ery = adata_to_sce(ery)

colData(ery) = cbind(colData(ery), data.frame(drift = erypar$drift,
                                            growth = erypar$growth,
                                           dpt_pseudotime = erypar$dpt_pseudotime))

plotReducedDim(ery, dimred = "X_umap_2d", colour_by = "drift")
plotReducedDim(ery, dimred = "X_umap_2d", colour_by = "growth")

# %% tags=[]
dptrange_drift = c(0.09, 0.19)
dptrange_growth = c(0.22, 0.27)
plot_cellparams(ery,
                drift_range = dptrange_drift,
                growth_range = dptrange_growth,
                ratio = 10,
                save=paste0(erydir, 'ery_dpt_plots.pdf'))

# %%
#Excluding SS2 samples to not affect the DE for early DPT
ery = ery[colData(ery)$data_type != 'SS2',]
ery = logNormCounts(ery)

# %% [markdown]
# ### Ery trajectory drift analysis

# %% tags=[]
source('utils/DE.R')
tokeep = gene_filter(assays(ery)$logcounts, min_expr = 0, min_fraction = 0.025, min_meanexpr = 0.05)
genes = row.names(ery)[tokeep]

deery1 = run_transitionDE(ery,
                 dpt_col = 'dpt_pseudotime',
                 dpt_range = dptrange_drift,
                 genes = genes,
                 allTFs = allTFs,
                 nknots = 8,
                 traj_name = 'ery_drift',
                 dir = erydir,
                 BPPARAM = BPPARAM)

# %% [markdown]
# ### Additional figures

# %% tags=[]
#Preparing data for plotting
npoints=100

#Long version of data unscaled
y = predictSmooth(deery1$sim, gene = row.names(deery1$desig), nPoints = 100, tidy = TRUE)
y = data.table(y)[, yhat_scaled := scale(yhat), by='gene']
yw = dcast(y, time~gene, value.var='yhat')
xtime = yw$time
yw = yw[,-1]

#Wide data scaled
yw_scaled = as.matrix(yw)
row.names(yw_scaled) = xtime
yw_scaled = t(yw_scaled)
yw_scaled <- t(scale(t(yw_scaled)))


#Annotation for heatmaps
dpt_anno = rep('out', npoints)
dpt_anno[(xtime > dptrange_drift[1] & xtime < dptrange_drift[2])] = 'in'
annotation_col = data.frame(region = dpt_anno, row.names=colnames(yw_scaled))

# %% tags=[]
eryheat = pheatmap(yw_scaled,
              cluster_cols=FALSE,
              annotation_col=annotation_col,
              fontsize_row=4.5,
              show_colnames=FALSE,
                  clustering_method='ward')

cl = cutree(eryheat$tree_row, 6)
cl_str = paste0('group ', cl)
annotation_row = data.frame(row.names = names(cl), group = cl_str)


pal20 = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b',
  '#e377c2', '#b5bd61', '#17becf', '#aec7e8', '#ffbb78', '#98df8a',
  '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5',
  '#ad494a', '#8c6d31')
groupcol = pal20[1:length(unique(cl_str))]
names(groupcol) = unique(cl_str)
ann_colors = list(
    region = c(out="#00BFC4", `in`="#F8766D"),
    group = groupcol)


eryheat = pheatmap(yw_scaled,
              cluster_cols=FALSE,
              annotation_col=annotation_col,
              annotation_row=annotation_row,
              treeheight_row = 0,
              treeheight_col = 0,
              fontsize_row=4.5,
              show_colnames=FALSE,
              show_rownames=FALSE,
              annotation_colors = ann_colors,
                  clustering_method='ward')
eryheat
ggsave(paste0(erydir, '/ery_drift_heatmap_clusters.pdf'), eryheat, height= 16, width=8)
ggsave(paste0(erydir, '/ery_drift_heatmap_clusters.png'), eryheat, height= 16, width=8)

# %% tags=[]
#Plotting as lines
toplot = y
toplot$cluster = cl_str[match(toplot$gene, names(cl))]
default_20 = c('#1f77b4', '#ff7f0e', '#279e68', '#d62728', '#aa40fc', '#8c564b', '#e377c2', '#b5bd61',
               '#17becf', '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2',
               '#dbdb8d', '#9edae5', '#ad494a', '#8c6d31')

clnames = unique(toplot$cluster)
clnames = clnames[order(clnames)]

plotlist = list()
n = 1
for (i in clnames){
    plotlist[[i]] = ggplot(toplot[toplot$cluster == i,],
       aes(x=time, y=yhat_scaled, group=gene)) +
       geom_line(alpha=0.15, color=default_20[n]) +
       geom_line(data = as.data.frame(colData(ery)), aes(x = dpt_pseudotime, y=drift*20, group='test'), color='grey', linetype='dashed') +
       xlab('Pseudotime') + 
       theme_paper +
       theme(legend.key.height=unit(0.35, 'cm')) +
       ggtitle(i) +
       scale_y_continuous(name = 'Scaled expression', sec.axis = sec_axis(~.*1/20, name="Diff. rate")) + 
       theme(axis.text.y.right = element_text(color = "grey"),
        axis.title.y.right = element_text(color = "grey"))
    n = n + 1
    }
groupplots = cowplot::plot_grid(plotlist = plotlist)
ggsave(paste0(erydir, '/ery_drift_groups.pdf'), groupplots, height= 4, width=7)
groupplots


# %% [markdown] tags=[]
# #### Chosen erythroid regulators

# %% tags=[]
genes_toplot = list(Ery_genes = c('Klf1', 'Gata1', 'Zfpm1', 'Tal1'), #also could include Nfe2, Ldb1
                   globins = c('Hba-a1', 'Hba-a2'))

yery = predictSmooth(deery1$sim, gene = unlist(genes_toplot), nPoints = 100, tidy = TRUE)
yery = data.table(yery)[, yhat_scaled := scale(yhat), by='gene']
ywery = dcast(yery, time~gene, value.var='yhat')
xtime = ywery$time
ywery = ywery[,-1]


plotlist2 = list()
for (i in names(genes_toplot)){
    plotlist2[[i]] = ggplot(yery[yery$gene %in% genes_toplot[[i]],],
       aes(x=time, y=yhat, color=gene)) +
       geom_line(data = as.data.frame(colData(ery)), aes(x = dpt_pseudotime, y=drift*30), color='grey', linetype='dashed') +
       geom_line(alpha=0.65, size=1) +
       xlab('Pseudotime') +
       ggtitle(i) +
       theme_paper +
       theme(legend.key.height=unit(0.35, 'cm')) +
       scale_color_manual(values=default_20) +
       scale_y_continuous(name = 'Norm. expression', sec.axis = sec_axis(~.*1/30, name="Diff. rate")) + 
       theme_bw(base_size = 8.5, base_family = "",
                       base_line_size = 0.5, base_rect_size = 0.5) + 
       theme(axis.text.y.right = element_text(color = "grey"),
        axis.title.y.right = element_text(color = "grey"))
    }
eryplots = cowplot::plot_grid(plotlist = plotlist2)
eryplots
ggsave(paste0(erydir, '/ery_drift_chosengenes.pdf'), eryplots, height= 2.5, width=7)

# %% [markdown]
# ### Ery trajectory growth analysis

# %%
deery2 = run_transitionDE(ery,
                 dpt_col = 'dpt_pseudotime',
                 dpt_range = dptrange_growth,
#                   genes = c('Actb', 'Elane', 'Ms4a3', 'Mpo'),
                 genes = genes,
                 allTFs = allTFs,
                 nknots = 8,
                 traj_name = 'ery_growth',
                 dir = erydir,
                 BPPARAM = BPPARAM)

# %% [markdown]
# ## DC trajectory

# %% tags=[]
DCdir = paste0(DIR_DE, 'DC/')
dir.create(DCdir, recursive=TRUE)

DCpar = fread('./PD_model/clu_6/tables/table_all_parameters_clu_6.csv')
DC = data[DCpar$V1, ]


# %% tags=[]
DC = adata_to_sce(DC)

colData(DC) = cbind(colData(DC), data.frame(drift = DCpar$drift,
                                            growth = DCpar$growth,
                                           dpt_pseudotime = DCpar$dpt_pseudotime))

plotReducedDim(DC, dimred = "X_umap_2d", colour_by = "drift")
plotReducedDim(DC, dimred = "X_umap_2d", colour_by = "growth")

# %%
#Clipping the long tail with very few cells
DC = DC[, colData(DC)$dpt_pseudotime < 0.115]

# %% tags=[]
dptrange_drift = c(0.037, 0.085)
dptrange_growth = c(0.08, 0.105)
plot_cellparams(DC,
                drift_range = dptrange_drift,
                growth_range = dptrange_growth,
                ratio = 200,
                save=paste0(DCdir, 'DC_dpt_plots.pdf'))

# %%
#Excluding SS2 samples to not affect the DE for early DPT
DC = DC[colData(DC)$data_type != 'SS2',]
DC = logNormCounts(DC)

# %%
source('utils/DE.R')
tokeep = gene_filter(assays(DC)$logcounts, min_expr = 0, min_fraction = 0.025, min_meanexpr = 0.05)
genes = row.names(DC)[tokeep]

deDC1 = run_transitionDE(DC,
                 dpt_col = 'dpt_pseudotime',
                 dpt_range = dptrange_drift,
                 genes = genes,
                 allTFs = allTFs,
                 nknots = 8,
                 traj_name = 'DC_drift',
                 dir = DCdir,
                 BPPARAM = BPPARAM)

# %%
deDC2 = run_transitionDE(DC,
                 dpt_col = 'dpt_pseudotime',
                 dpt_range = dptrange_growth,
                 genes = genes,
                 allTFs = allTFs,        
                 nknots = 8,
                 traj_name = 'DC_growth',
                 dir = DCdir,
                 BPPARAM = BPPARAM)


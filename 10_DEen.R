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
library(data.table)
library(ggplot2)
source('utils/DE.R')
library(pheatmap)
library(data.table)

#ggplot theme
theme_paper = theme_bw(base_size = 7, base_family = "",
                         base_line_size = 0.5, base_rect_size = 0.5)
theme_set(theme_paper)

DIR_DE = './DE/'


# %% [markdown]
# ## Megakaryocyte trajectory

# %%
mkdir = paste0(DIR_DE, 'Mk/')
dir.create(mkdir, recursive=TRUE)

# %%
gde = fread('./DE/Mk/mk_growth__desig.csv')
colnames(gde)[1] = 'symbol'

# %%
# Running enrichment with enrichr
gen = run_enrichr(gde$symbol,
            gene_sets=c('MSigDB_Hallmark_2020', 'Reactome_2016'),
            figure_prefix=paste0(mkdir, '/mk_growth_'))
#Extracting genes from categories of interest
e2f = gen[gen$Gene_set == 'MSigDB_Hallmark_2020' & gen$Term == 'E2F Targets', 'Genes']
e2f = strsplit(e2f, split = ';')[[1]]
e2f = tools::toTitleCase(tolower(e2f))
g2mc = gen[gen$Gene_set == 'MSigDB_Hallmark_2020' & gen$Term == 'G2-M Checkpoint', 'Genes']
g2mc = strsplit(g2mc, split = ';')[[1]]
g2mc = tools::toTitleCase(tolower(g2mc))
cc = gen[gen$Gene_set == 'Reactome_2016' & gen$Term == 'Cell Cycle Homo sapiens R-HSA-1640170', 'Genes']
cc = strsplit(cc, split = ';')[[1]]
cc = tools::toTitleCase(tolower(cc))

# %%
#Preparing data for plotting
data = fread('./DE/Mk/mk_growth__yhatsmooth.csv')
dataw = dcast(data, time~gene, value.var='yhat')
datawM = as.matrix(dataw[,-1])
row.names(datawM) = dataw$time
datawM = t(datawM)
datawM <- t(scale(t(datawM)))

npoints=100
dptrange_growth = c(0.01, 0.028)

xtime = as.numeric(colnames(datawM))
dpt_anno = rep('out', npoints)
dpt_anno[(xtime > dptrange_growth[1] & xtime < dptrange_growth[2])] = 'in'
annotation_col = data.frame(region = dpt_anno, row.names=colnames(datawM))

annotation_row = data.frame(row.names = row.names(datawM),
                            E2F = ifelse(row.names(datawM) %in% e2f, 'E2F targets', 'other'),
                            G2MC = ifelse(row.names(datawM) %in% g2mc, 'G2-M checkpoint', 'other'),
                            CC = ifelse(row.names(datawM) %in% cc, 'Cell cycle', 'other'))
ann_colors = list(
    region = c(out="#00BFC4", `in`="#F8766D"),
    E2F = c(`E2F targets` ='blue', other='lightgrey'),
    G2MC = c(`G2-M checkpoint` ='purple', other='lightgrey'),
    CC = c(`Cell cycle` ='darkorange', other='lightgrey')
    )

# %%
p1 = pheatmap(datawM,
              cluster_cols=FALSE,
              annotation_col=annotation_col,
              annotation_row=annotation_row,
              annotation_colors = ann_colors,
              fontsize_row=4.5,
              treeheight_row=0,
              show_colnames=FALSE)
ggsave(paste0(mkdir, '/mk_growth_heatmap_annotated.pdf'), p1, height= 20, width=14)

# %%
p2 = ggplot(data[data$gene %in% c('Pf4', 'Procr',
                             'Mcm5', 'Mki67', 'Hist2h2ac', 'E2f8', 'Lmnb1'),],
       aes(x=time, y=log2(yhat+1), color=gene)) +
       geom_line() +
       xlab('Pseudotime') + 
       ylab('Log2(Expression + 1)') + 
       geom_vline(xintercept = dptrange_growth, linetype='dashed', color = 'red')
ggsave(paste0(mkdir, '/mk_growth_example_genes.pdf'), p2, height=2, width=4)

# %%

# %%

# %%
neudir = paste0(DIR_DE, 'Neu/')
dir.create(neudir, recursive=TRUE)

# %%
dde = fread('./DE/Neu/neu_growth__desig.csv')
colnames(dde)[1] = 'symbol'

# %%
#Preparing data for plotting
data = fread('./DE/Neu/neu_drift__yhatsmooth.csv')
dataw = dcast(data, time~gene, value.var='yhat')
datawM = as.matrix(dataw[,-1])
row.names(datawM) = dataw$time
datawM = t(datawM)
datawM <- t(scale(t(datawM)))

npoints=100
dptrange_growth = c(0.01, 0.028)

xtime = as.numeric(colnames(datawM))
dpt_anno = rep('out', npoints)
dpt_anno[(xtime > dptrange_growth[1] & xtime < dptrange_growth[2])] = 'in'
annotation_col = data.frame(region = dpt_anno, row.names=colnames(datawM))

# %%
p1 = pheatmap(datawM,
              cluster_cols=FALSE,
              annotation_col=annotation_col,
              fontsize_row=4.5,
              treeheight_row=0,
              show_colnames=FALSE)
# ggsave(paste0(neudir, '/neu_growth_heatmap_annotated.pdf'), p1, height= 20, width=14)

# %%
theme_paper = theme_bw(base_size = 14, base_family = "",
                         base_line_size = 0.5, base_rect_size = 0.5)
dptrange_drift = c(0.075, 0.11)
earlygenes = c('Elane', 'Cst7', 'Ms4a3', 'Cebpe')
lategenes = c('S100a8', 'Clec4a2', 'Ms4a3', 'Cebpe')
othergenes = c('Irf8', 'Cd34', 'Flt3', 'Elane', 'S100a8')

p3 = ggplot(data[data$gene %in% c(earlygenes, lategenes, othergenes),],
       aes(x=time, y=log2(yhat+1), color=gene)) +
       geom_line() +
       xlab('Pseudotime') + 
       ylab('Log2(Expression + 1)') + 
       geom_vline(xintercept = dptrange_drift, linetype='dashed', color = 'blue') + 
        theme_paper
p3
# ggsave(paste0(mkdir, '/neu_drift_example_genes.pdf'), p2, height=2, width=4)

# %%
genes = c('Dstn', 'Elane', 'Ly6c2', 'S100a8')

p3 = ggplot(data[data$gene %in% genes,],
       aes(x=time, y=log2(yhat+1), color=gene)) +
       geom_line() +
       xlab('Pseudotime') + 
       ylab('Log2(Expression + 1)') + 
       geom_vline(xintercept = dptrange_drift, linetype='dashed', color = 'blue') + 
        theme_paper
p3

# %%
genes = c('Dstn', 'Ly6c2', 'Elane', 'Ikzf2', 'Wfdc17')

p3 = ggplot(data[data$gene %in% genes,],
       aes(x=time, y=log2(yhat+1), color=gene)) +
       geom_line() +
       xlab('Pseudotime') + 
       ylab('Log2(Expression + 1)') + 
       geom_vline(xintercept = dptrange_drift, linetype='dashed', color = 'blue') + 
        theme_paper
p3

# %% [markdown]
# Highlight the different groups:
# correlated best with the differentation rate?
# gradually increasing throughout, gradually decreasing throughout
# with a bump before (these include irf8 and Ikzf2, cd34)
#

# %%

# %%

# %%

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
mk = logNormCounts(mk)
plotReducedDim(mk, dimred = "X_umap_2d", colour_by = "drift")

# %%
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


# %%
dptrange_drift = c(0.0125, 0.0375)
dptrange_growth = c(0.01, 0.028)
plot_cellparams(mk,
                drift_range = dptrange_drift,
                growth_range = dptrange_growth,
                ratio = 10,
                save=paste0(mkdir, 'mk_dpt_plots.pdf'))

# %% tags=[]
tokeep = gene_filter(assays(mk)$logcounts, min_expr = 0, min_fraction = 0.025, min_meanexpr = 0.05)
genes = row.names(mk)[tokeep]

run_transitionDE(mk,
                  dpt_col = 'dpt_pseudotime',
                  dpt_range = dptrange_drift,
#                   genes = c('Actb', 'Pf4', 'Procr', 'Gata1'),
                  genes = genes,
                  traj_name = 'mk_drift',
                 nknots = 7,
                  dir = mkdir,
                  BPPARAM = BPPARAM)

# %%
run_transitionDE(mk,
                  dpt_col = 'dpt_pseudotime',
                  dpt_range = dptrange_growth,
#                   genes = c('Actb', 'Pf4', 'Procr', 'Gata1'),
                  genes = genes,
                  traj_name = 'mk_growth',
                 nknots = 7,
                  dir = mkdir,
                  BPPARAM = BPPARAM)

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
neu = logNormCounts(neu)
plotReducedDim(neu, dimred = "X_umap_2d", colour_by = "drift")

# %%
dptrange_drift = c(0.075, 0.11)
dptrange_growth = c(0.01, 0.028)
plot_cellparams(neu,
                drift_range = dptrange_drift,
                growth_range = dptrange_growth,
                ratio = 25,
                save=paste0(neudir, 'neu_dpt_plots.pdf'))

# %%
tokeep = gene_filter(assays(neu)$logcounts, min_expr = 0, min_fraction = 0.025, min_meanexpr = 0.05)
genes = row.names(neu)[tokeep]

run_transitionDE(neu,
                 dpt_col = 'dpt_pseudotime',
                 dpt_range = dptrange_drift,
#                   genes = c('Actb', 'Elane', 'Ms4a3', 'Mpo'),
                 genes = genes,
                 nknots = 8,
                 traj_name = 'neu_drift',
                 dir = neudir,
                 BPPARAM = BPPARAM)

# %%
run_transitionDE(neu,
                 dpt_col = 'dpt_pseudotime',
                 dpt_range = dptrange_growth,
#                   genes = c('Actb', 'Elane', 'Ms4a3', 'Mpo'),
                 genes = genes,
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
ery = logNormCounts(ery)
plotReducedDim(ery, dimred = "X_umap_2d", colour_by = "drift")
plotReducedDim(ery, dimred = "X_umap_2d", colour_by = "growth")

# %%
dptrange_drift = c(0.15, 0.195)
dptrange_growth = c(0.22, 0.27)
plot_cellparams(ery,
                drift_range = dptrange_drift,
                growth_range = dptrange_growth,
                ratio = 10,
                save=paste0(erydir, 'ery_dpt_plots.pdf'))

# %%
source('utils/DE.R')
tokeep = gene_filter(assays(ery)$logcounts, min_expr = 0, min_fraction = 0.025, min_meanexpr = 0.05)
genes = row.names(ery)[tokeep]

run_transitionDE(ery,
                 dpt_col = 'dpt_pseudotime',
                 dpt_range = dptrange_drift,
#                   genes = c('Actb', 'Elane', 'Ms4a3', 'Mpo'),
                 genes = genes,
                 nknots = 8,
                 traj_name = 'ery_drift',
                 dir = erydir,
                 BPPARAM = BPPARAM)

# %%
run_transitionDE(ery,
                 dpt_col = 'dpt_pseudotime',
                 dpt_range = dptrange_growth,
#                   genes = c('Actb', 'Elane', 'Ms4a3', 'Mpo'),
                 genes = genes,
                 nknots = 8,
                 traj_name = 'ery_growth',
                 dir = erydir,
                 BPPARAM = BPPARAM)

# %%

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
DC = logNormCounts(DC)
plotReducedDim(DC, dimred = "X_umap_2d", colour_by = "drift")
plotReducedDim(DC, dimred = "X_umap_2d", colour_by = "growth")

# %%
dptrange_drift = c(0.055, 0.08)
dptrange_growth = c(0.08, 0.105)
plot_cellparams(DC,
                drift_range = dptrange_drift,
                growth_range = dptrange_growth,
                ratio = 200,
                save=paste0(DCdir, 'DC_dpt_plots.pdf'))

# %%
source('utils/DE.R')
tokeep = gene_filter(assays(DC)$logcounts, min_expr = 0, min_fraction = 0.025, min_meanexpr = 0.05)
genes = row.names(DC)[tokeep]

run_transitionDE(DC,
                 dpt_col = 'dpt_pseudotime',
                 dpt_range = dptrange_drift,
#                   genes = c('Actb', 'Elane', 'Ms4a3', 'Mpo'),
                 genes = genes,
                 nknots = 8,
                 traj_name = 'DC_drift',
                 dir = DCdir,
                 BPPARAM = BPPARAM)

# %%
run_transitionDE(DC,
                 dpt_col = 'dpt_pseudotime',
                 dpt_range = dptrange_growth,
#                   genes = c('Actb', 'Elane', 'Ms4a3', 'Mpo'),
                 genes = genes,
                 nknots = 8,
                 traj_name = 'DC_growth',
                 dir = DCdir,
                 BPPARAM = BPPARAM)

# %%

# %% [markdown]
# ## Ly trajectory

# %% tags=[]
Lydir = paste0(DIR_DE, 'Ly/')
dir.create(Lydir, recursive=TRUE)

Lypar = fread('./PD_model/clu_16/tables/table_all_parameters_clu_16.csv')
Ly = data[Lypar$V1, ]

# %% tags=[]
Ly = adata_to_sce(Ly)

colData(Ly) = cbind(colData(Ly), data.frame(drift = Lypar$drift,
                                            growth = Lypar$growth,
                                           dpt_pseudotime = Lypar$dpt_pseudotime))
Ly = logNormCounts(Ly)
plotReducedDim(Ly, dimred = "X_umap_2d", colour_by = "drift")
plotReducedDim(Ly, dimred = "X_umap_2d", colour_by = "growth")

# %%
dptrange_drift = c(0.005, 0.04)
dptrange_growth = c(0.01, 0.08)
plot_cellparams(Ly,
                drift_range = dptrange_drift,
                growth_range = dptrange_growth,
                ratio = 100,
                save=paste0(Lydir, 'Ly_dpt_plots.pdf'))

# %%
source('utils/DE.R')
tokeep = gene_filter(assays(Ly)$logcounts, min_expr = 0, min_fraction = 0.025, min_meanexpr = 0.05)
genes = row.names(Ly)[tokeep]

run_transitionDE(Ly,
                 dpt_col = 'dpt_pseudotime',
                 dpt_range = dptrange_drift,
#                   genes = c('Actb', 'Elane', 'Ms4a3', 'Mpo'),
                 genes = genes,
                 nknots = 6,
                 traj_name = 'Ly_drift',
                 dir = Lydir,
                 BPPARAM = BPPARAM)

# %%
run_transitionDE(Ly,
                 dpt_col = 'dpt_pseudotime',
                 dpt_range = dptrange_growth,
#                   genes = c('Actb', 'Elane', 'Ms4a3', 'Mpo'),
                 genes = genes,
                 nknots = 6,
                 traj_name = 'Ly_growth',
                 dir = Lydir,
                 BPPARAM = BPPARAM)

# %%
# dptrange_drift = c(0.0125, 0.0375)
# mkdens_drift = ggplot(data.frame(colData(mk)), aes(x = dpt_pseudotime)) +
#   geom_density(size = 1.5) +
#   geom_line(aes(x = dpt_pseudotime, y = drift*100), color = 'red') +
#   scale_y_continuous(name = 'cell density', sec.axis = sec_axis(~.*0.01, name="drift Axis")) + 
#   geom_vline(xintercept = dptrange_drift, alpha = 0.5)
# mkdens_drift

# dptrange_growth = c(0.01, 0.028)
# mkdens_growth = ggplot(data.frame(colData(mk)), aes(x = dpt_pseudotime)) +
#   geom_density(size = 1.5) +
#   geom_line(aes(x = dpt_pseudotime, y = growth*10), color = 'red') +
#   scale_y_continuous(name = 'cell density', sec.axis = sec_axis(~.*0.1, name="growth Axis")) + 
#   geom_vline(xintercept = dptrange_growth, alpha = 0.5)
# mkdens_growth

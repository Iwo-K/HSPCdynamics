##' Generates a boolean filter to select genes passing the expression levels criteria.
##'
##' Genes are filtered adccording to two criteria: global mean expression (meanexpr_tr)
##' and fraction of samples with gene expressed (expr_tr and min_fraction)
##' Global mean expression threshold
##' Genes with values greater than the threshold as set to TRUE.
##' @title Get filter for genes with low expression
##' @param x data.frame or data.tablet to be filtered
##' @param min_expr threshold for genes to be considered expressed in the min_fraction argument
##' @param min_fraction fraction of samples in which genes should be above the expr_tr
##' @param min_meanexpr threshold for the gene global mean expression level
##' @return boolean vector for genes passing the thresholds, can be used for filtering
##' @author Iwo Kucinski
gene_filter = function(x, min_expr = 0, min_fraction = 0, min_meanexpr = 0){
  n = ncol(x)
  above_tr = apply(x, 2, function(x) x > min_expr)
  expr_no = rowSums(above_tr)

  keep1 = expr_no/n > min_fraction
  print(paste0("Genes expressed (min_expr >", min_expr, ") in >", min_fraction, " of cells: ", sum(keep1)))
# Filtering by mean expr
  keep2 = rowMeans(x) > min_meanexpr
  print(paste0("Genes with mean expression > ", min_meanexpr, ": ", sum(keep2)))

  return(keep1 & keep2)
}


##' Find genes DE within a chosen pseudotime region.
##'
##' Uses the tradeSeq framework (startVsEndTest speicifically) to
##' test if a gene changes expresion level within the specified
##' pseudotime interval
##' @title
##' @param sce SingleCellExperiment object with counts in assays slot
##' @param dpt_col  name of the column in the colData(sce) with pseudotime
##' @param dpt_range 2-element vector specifying the psueoditme interval
##' @param genes vector with gene names (need to match the gene names in sce)
##' @param fdr_tr float, FDR threshold used to call genes significanta
##' @param nknots int, number of knots for the GAM bmodel
##' @param l2fc_tr float, log2FoldChange threshold for DE call
##' @param traj_name str, name of the trajectory (used for saving)
##' @param dir str, path where output files should be saved
##' @param BPPARAM BiocParallel:bpparam object with specified workers
##' @return nothing, just saves files to the indicated folder
##' @author John Doe
run_transitionDE = function(sce,
                 dpt_col = 'dpt_pseudotime',
                 dpt_range = NULL,
                 genes = NULL,
                 allTFs = NULL,
                 fdr_tr = 0.1,
                 nknots = 5,
                 l2fc_tr = 1,
                 traj_name = 'mk',
                 dir = './',
                 BPPARAM = NULL) {

  require(tradeSeq)
  require(pheatmap)
  require(ggplot2)

  # Path + name for the output files
  filename = paste0(dir, traj_name, '_')

  # Running GAMs
  cellWeights = data.frame(lin = rep(1, dim(sce)[2]))
  pseudotime = data.frame(lin = colData(sce)[[dpt_col]])
  sim = fitGAM(as.matrix(assays(sce)$counts),
               pseudotime = pseudotime,
               cellWeights = cellWeights,
               genes = genes,
               nknots = nknots,
               verbose = TRUE,
               parallel = TRUE,
               BPPARAM = BPPARAM
               )

  knots = metadata(sim)$tradeSeq$knots
  print(knots)

  # DE test
  de = startVsEndTest(sim, pseudotimeValues = dpt_range, l2fc = l2fc_tr)
  de$padj <- p.adjust(de$pvalue, "fdr")
  de$padj2 = stats::pchisq(de$waldStat, df = de$df, lower.tail = FALSE, log.p = TRUE)
  de = de[order(de$padj2, decreasing = FALSE),]

  write.csv(de, paste0(filename, 'de.csv'))
  desig = de[de$padj < fdr_tr,]
  write.csv(desig, paste0(filename, '_desig.csv'))
  print(head(desig))
  print(paste0('No of DE genes: ', nrow(desig)))

  # Plotting top genes
  if (nrow(desig) > 0){
    pdf(paste0(filename, 'topgenes.pdf'))
    if (nrow(desig) > 30) n = 30 else n = nrow(desig)
    for (i in row.names(desig)[1:n]){
      g = plotSmoothers(sim, assays(sim)$counts, gene = i) +
        ggtitle(i) + geom_vline(xintercept=dpt_range, alpha = 0.5)
      print(g)
    }
dev.off()
}

  if (nrow(desig) >=2){
    yhatSmooth <- predictSmooth(sim, gene = row.names(desig), nPoints = 100, tidy = TRUE)
    write.csv(yhatSmooth, paste0(filename, '_yhatsmooth.csv'))     
      
    # Plotting heatmaps
    npoints = 100
    xtime = predictSmooth(sim, gene = row.names(de)[1], nPoints = npoints, tidy = TRUE)$time
    dpt_anno = rep('out', npoints)
    dpt_anno[(xtime > dpt_range[1] & xtime < dpt_range[2])] = 'in'
      
    #Heatmap for all significant genes
    n = nrow(desig)
    y <- predictSmooth(sim, gene = row.names(desig)[1:n], nPoints = npoints, tidy = FALSE)
    yscaled <- t(scale(t(y)))
    annotation_col = data.frame(region = dpt_anno, row.names=colnames(yscaled))

    #Only top 400 genes are drawn if number of DE is >400
    if (nrow(desig) > 400) nheat = 400 else nheat = nrow(desig)
    heat <- pheatmap(yscaled[1:nheat,],
                    cluster_cols = FALSE,
                    cluster_rows = TRUE,
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    fontsize_row=4.5,
                    annotation_col = annotation_col,
                    legend = TRUE,
                    silent = TRUE,
                    treeheight_row = 0,
                    treeheight_col = 0,
                    border_color = NA)
    ggsave(paste0(filename, 'heatmap.pdf'), heat, height = 2 + 0.06 * nrow(yscaled[1:nheat,]))
    print(heat)
    
    #Heatmap for TFs  
    if (!is.null(allTFs)){
        yscaled_tf = yscaled[row.names(yscaled) %in% allTFs,]
        heatTF <- pheatmap(yscaled_tf,
                    cluster_cols = FALSE,
                    cluster_rows = TRUE,
                    show_rownames = TRUE,
                    show_colnames = FALSE,
                    fontsize_row=6,
                    annotation_col = annotation_col,
                    legend = TRUE,
                    silent = TRUE,
                    treeheight_row = 0,
                    treeheight_col = 0,
                    border_color = NA)
        ggsave(paste0(filename, 'heatmapTF.pdf'), heatTF, width = 4.5, height = 2 + 0.06 * nrow(yscaled_tf))
#     print(heatTF)
        }

    #Log-transformed heatmap for all significant genes

#     ylog  = apply(y, 2, function(x) log2(x/mean(x)))
#     heatlog <- pheatmap(ylog,
#                     cluster_cols = FALSE,
#                     cluster_rows = TRUE,
#                     show_rownames = TRUE,
#                     show_colnames = FALSE,
#                     fontsize_row=4.5,
#                     annotation_col = annotation_col,
#                     legend = TRUE,
#                     silent = TRUE,
#                     treeheight_row = 0,
#                     treeheight_col = 0,
#                     border_color = NA)
#     ggsave(paste0(filename, 'heatmap_log.pdf'), heatlog, height = 2 + 0.06 * nrow(y))
                  
    return(list(sim = sim, desig = desig))
#     print(heatlog)


#     #Heatmap for top  50 significant genes
#     if (nrow(desig) > 50) n = 50 else n = nrow(desig)
#     y <- predictSmooth(sim, gene = row.names(desig)[1:n], nPoints = npoints, tidy = FALSE)
#     yscaled <- t(scale(t(y)))
#     annotation_col = data.frame(region = dpt_anno, row.names=colnames(yscaled))

#     heat <- pheatmap(yscaled,
#                     cluster_cols = FALSE,
#                     cluster_rows = TRUE,
#                     show_rownames = TRUE,
#                     show_colnames = FALSE,
#                     annotation_col = annotation_col,
#                     legend = TRUE,
#                     silent = TRUE,
#                     treeheight_row = 0,
#                     treeheight_col = 0,
#                     border_color = NA)
#     ggsave(paste0(filename, 'heatmap_top50.pdf'), heat, height = 1 + 0.12 * nrow(y))
#     print(heat)
  }

}

##' Run enrichR enrichment analysis
##'
##' Uses gseapy python package (via reticulate) to perform gene enrichment
##' analysis. Requires R/reticulate + python/matplotlib/gseapy environment setup.
##' Figure plotting depends on scanpy for figure settings
##' @title
##' @param genes character vector with gene symbols
##' @param gene_sets character vector with categories to be used, see enrichR for details
##' @param organism string, name of the organism
##' @param cutoff float, cutoff used for enrichment
##' @param output_dir string, path where gseapy output is saved
##' @param figure_prefx, string, if provided figures are saved for each gene_set
##' @param csv_prefx, string, if provided the output table is saved as csv
##' @return data.frame with enrichment results
##' @author Iwo Kucinski
run_enrichr = function(genes,
                       gene_sets = c('MSigDB_Hallmark_2020'),
                       organism = 'Mouse',
                       cutoff=0.05,
                       output_dir = 'gseapy_output/test',
                       figure_prefix = NULL,
                       csv_prefix = NULL){
  require(reticulate)
  gseapy <- import("gseapy")
  sc <- import("scanpy")
  mpl <- import('matplotlib')
  sc$set_figure_params(dpi=100,
                       figsize=c(4,4),
                       color_map = 'viridis',
                       dpi_save = 350)
  mpl$rcParams['pdf.fonttype'] = 42 #Ensures readable fonts in illustrator
  mpl$rcParams['ps.fonttype'] = 42

  plt <- import("matplotlib.pyplot")

  enr = gseapy$enrichr(gene_list = genes,
               gene_sets = gene_sets,
               organism='Mouse',
               outdir='gseapy_out/enrichr',
               cutoff=cutoff)
  res = enr$results

  if (!is.null(figure_prefix)){
    for (i in gene_sets){
      ressub = res[res$Gene_set == i,]
      gseapy$plot$dotplot(ressub, title=i,
                          cmap='viridis_r',
                          ofname=paste0(figure_prefix, '_', i,  'enrichr.pdf'))
    }
  }

  if (!is.null(csv_prefix)){
    write.csv(res, paste0(csv_prefix, 'enrichr.csv'))
  }
  return(res)

}

##' Little convenient function which takes a ; separated gene names (from enrichhr) and convert to Title case (for mouse)
str2mouse = function(x){
    x = strsplit(x, split = ';')[[1]]
    x = tools::toTitleCase(tolower(x))
    }


                   
                   
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
analyse_cor = function(sim, pseudotime, drift, pseudotime_lim=NULL, TFs, prefix, scale_ratio=250){
                      
    #First running an association test to seleting genes with varying along pseudotime with l2fc of at least1
    de = associationTest(sim, l2fc=1, contrastType='end')
    de$fdr = p.adjust(de$pvalue, method='fdr')
    de = de[complete.cases(de),]
    desig = de[de$fdr < 0.01,]
            
    #predicting smoothed expression for each de gene
    yhat = predictSmooth(sim, gene = row.names(desig), nPoints = 100, tidy = TRUE)
    yhat = data.table(yhat)[, yhat_scaled := scale(yhat), by='gene']
    
    if(!is.null(pseudotime_lim)){
        print(paste0("Clipping pseudotime to: ", pseudotime_lim))
    yhat = yhat[(time > pseudotime_lim[1]) & (time < pseudotime_lim[2]),]
    }
    
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
                        
    ratio=1/max(ddrift)*10

    gup = ggplot(yhat[gene %in% names(toptf),], aes(x=time, y=yhat_scaled, color=gene)) + 
        geom_line(size=0.8) +
        geom_line(data= data.frame(), aes(x = yhatw$time, y = ddrift*scale_ratio), , color='grey', linetype='dashed', size = 0.8) +
        scale_y_continuous(name = 'Scaled expression', sec.axis = sec_axis(~.*1/scale_ratio, name="Diff. rate derivative")) +
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
        geom_line(data= data.frame(), aes(x = yhatw$time, y = ddrift*scale_ratio), , color='grey', linetype='dashed', size = 0.8) +
        scale_y_continuous(name = 'Scaled expression', sec.axis = sec_axis(~.*1/scale_ratio, name="Diff. rate derivative")) +
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
                        
    fwrite(data.table(gene = names(ddrift_cors), cor = ddrift_cors), paste0(prefix, '_ddrift_cors.csv'))
}
                        

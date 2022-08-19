library("scran")
library(ggplot2)

#' Normalise a data.frame with counts
#'
#' Function takes raw counts from and RNAseq experiments and normalises them using a chosen method. Updated to work with newer version of scran (bioconductor >=3.11 )
#' @param counts a dataframe with rows as genes and samples/cells as columns
#' @param norm_method, DESeq2 or scran
#' @param whichgenes, specifies which genes are supposed to be used for normalisation, choice between mouse (all mouse genes "ENSM.."), "ERCC" (the ERCC spikeins) and custom (need to specify the genes in the genelist parameter)
#' @param genelist a vector with gene names used for normalisation
#' @param doclustering, for scran option only, logical whether clustering should be performed prior to normalisation (this helps when there are wildly different cell populations e.g. with majority of genes differentially expressed
#' @param minclusterize, for scran option only, numeric, a minimal size of clusters to be used for scran quick clustering
#' @param maxpoolsize, for scran option only, numeric a maximum pool size that scran uses for normalisation, needs to be <= minclustersize
#' @return A data.frame with normalised expression values
#' @export
scnormalise_v2 = function(counts,
                       norm_method = "DESeq2",
                       whichgenes = "mouse",
                       genelist = NULL,
                       doclustering = FALSE,
                       minclustersize = 200,
                       maxpoolsize = 100){

  if(whichgenes == "mouse") countsbase = counts[grep("^EN.*", row.names(counts)),]
  else if(whichgenes == "ERCC") countsbase = counts[grep("^ERCC.*", row.names(counts)),]
  else if(whichgenes == "custom") countsbase = counts[genelist,]

  if (norm_method == "DESeq2"){
    sf = DESeq2::estimateSizeFactorsForMatrix(as.matrix(countsbase))
    ncounts = t(t(counts)/sf)

    histogram = ggplot(data.frame(), aes(x=sf)) +
      geom_histogram(binwidth=.5, colour="black", fill="white") +
      theme_bw() +
      xlab("Size factors per cell")

    sumcorrelation = ggplot(data.frame(), aes(x = sf, y = colSums(countsbase))) +
      geom_point(alpha = 0.35) +
      theme_bw() +
      xlab("Size factors per cell") +
      ylab("Total mapped reads per cell")

    multiplot(histogram, sumcorrelation, cols = 2)
  }
  if (norm_method == "scran"){
    sfscran = getSCRANSF_v2(as.matrix(countsbase), doclustering = doclustering, minclustersize = minclustersize, maxpoolsize = maxpoolsize)
    ncounts = t(t(counts) / sfscran)

    histogram = ggplot(data.frame(), aes(x=sfscran)) +
      geom_histogram(binwidth=.5, colour="black", fill="white") +
      theme_bw() +
      xlab("Size factors per cell")

    sumcorrelation = ggplot(data.frame(), aes(x = sfscran, y = colSums(countsbase))) +
      geom_point(alpha = 0.35) +
      theme_bw() +
      xlab("Size factors per cell") +
      ylab("Total mapped reads per cell")

    multiplot(histogram, sumcorrelation, cols = 2)
  }
  return(ncounts)  
}

#' Compute sum factors using the scran package ComputeSumFactors normalisation
#'
#' A simple wrapper of the ComputeSumFactors functin froms cran (combines clustering)
#' @param counts dataframe with rows as genes and samples/cells as columns
#' @param doclustering, logical whether clustering should be performed prior to normalisation (this helps when there are wildly different cell populations e.g. with majority of genes differentially expressed
#' @param minclusterize, numeric, a minimal size of clusters to be used for scran quick clustering
#' @param maxpoolsize, numeric a maximum pool size that scran uses for normalisation, needs to be <= minclustersize
#' @return A data.frame with normalised expression values
#' @export
getSCRANSF_v2 = function(counts, doclustering = FALSE, minclustersize = 200, maxpoolsize = 100){
  
  if(doclustering){
    clusters <- scran::quickCluster(as.matrix(counts), min.size = minclustersize)
      #The new scran version changed behaviour and computeSumFactors works only on SingleCell Experiment objects  
    #sf = scran::computeSumFactors(as.matrix(counts), clusters = clusters, size = seq(20, maxpoolsize, 5))
    sf = scran::calculateSumFactors(as.matrix(counts), clusters = clusters, size = seq(20, maxpoolsize, 5))

  }
  else {
      sf = scran::calculateSumFactors(as.matrix(counts))
      #sf = scran::computeSumFactors(as.matrix(counts))
    }
  return(sf)
}


#' Function running DESeq2 differential expression
#' @export
#' @param counts dataframe with rows as genes and samples/cells as columns (HTSeq output)
#' @param metadata dataframe with samples/cells as row and columns specifying metadata for each samples
#' @param genedata dataframe with rows as genes and gene information in columns, this is the output of the function makeGENEDATA
#' @param contrast a vector with three element, the first element specifies the column in metadata which will be used for pairwise comparison (source of the factor), the second elemetn specifies the factor level of a treated sample (numerator in the fold change) and the third element the factor level of the control sample (the denominator in the fold change)
#' @param blocking a name of the column in the metadata to be used as a blocking factor, if diferent then NULL this will initiate DESeq2 call as: ~condition + blocking. The contrast used for pairwise comparison stays the same.
#' @param shrinkage logical, whether shrinkage option of the DESeq2 package should be used. This runs the DESeq2 with betaPrior=TRUE, affects differential expression and replaces the fold changes with shrunken fold changes (using the DESeq2 normal type of shrinkage)
#' @param calcshrinkage logical, if shrinkage=FALSE, this calculates shrunken fold changes and adds them to the output data frame in addition to the normal fold changes and generates respective MA plots (as labelled). There is no harm in calculating the shrinkages each time but this significantly increases computation time, so this function allows to switch it off
#' @param independentFiltering  logical, whether DESeq2 should perform independent filtering
#' @param exprfilter, numberic, sets the minimum expression value (minimum rowMean of normalised values) for filtering, all the genes below the number are removed from the dataset. make sure that the independent filtering is switched to FALSE for this to work, usually this is not necessary anyway.
#' @param padj_threshold numeric, specifies the padj cutoff
#' @param lfc_threshold numeric, specifies the log2FoldChange cutoff (for Deseq2 testing) - this is not applied to the shrinkage calculation (lfcShrink function)
#' @param pdfname string, specifies the base name for the pdf report generated
#' @param output string, either "sig" or "all", to output the list of all genes or only the significant genes respectively
#' @param subsetcomparison logical, specifies it the data should be first subset only for the factor levels compared and rejecting samples unused in the comparison. This affects the dispersion estimation, and can make the method more sensitive for samples with little DE and more conservative for samples with a lot of DE genes. Note however, that if there are too few samples left, the estimation of the dispersion will innaccurate.
#' @param normmethod, DESeq2, scran (no clustering version) or other. If "sf" is selected, then the a vector of sizefactors needs to be provided as the sf argument
#' @param sf a vector of size factors, used if normmethod is "sf"
runDESeq2 = function(counts, metadata, genedata,
                     contrast = c("condition", "treated", "control"),
                     blocking = NULL,
                     shrinkage = FALSE,
                     calcshrinkage = TRUE,
                     independentFiltering = TRUE,
                     exprfilter = NULL,
                     padj_threshold = 0.1,
                     lfc_threshold = 0,
                     filename = "DEseq2_report",
                     output = "sig",
                     subsetcomparison = FALSE,
                     normmethod = "DESeq2",
                     sf,
                     parallel = FALSE,
                     save = TRUE){

  #' Libraries and settings
  require("DESeq2", quietly = TRUE)
  require("BiocParallel", quietly = TRUE)
  require("scran", quietly = TRUE)
  theme_set(theme_bw())
  #if(parallel) register(MulticoreParam(4))

  #' Mode of analysis (is subsetcomparison is TRUE, then the data is first subset for the factor levels of intersts
  if(subsetcomparison){
    sub = metadata[,contrast[1]] %in% c(contrast[2], contrast[3])
    metadata = metadata[sub,]
    counts = counts[,sub]
  }

  if(all(colnames(counts) == metadata$cellid)) print("Ids are matching!")
  else stop("Please match the cell IDs in data columns and the cellid column of metadata.")

  design = data.frame(row.names = colnames(counts))
  design[,contrast[1]] = metadata[,contrast[1]]
  design[,contrast[1]] = as.factor(design[,contrast[1]])
  design[,contrast[1]] = relevel(design[,contrast[1]], contrast[3]) #Here changed contrast[2] to contrast[3] to make the correct reference level

  dds = DESeqDataSetFromMatrix(countData = counts, colData = design, design = as.formula(paste0("~", contrast[1])))

  if(!is.null(blocking)){
    design[,blocking] = metadata[,blocking]
    dds = DESeqDataSetFromMatrix(countData = counts, colData = design, design = as.formula(paste0("~", contrast[1], "+", blocking)))
  }

  #' Normalising the data
  print(paste("Estimating size factors with: ", normmethod))
  if(normmethod == "DESeq2") sizeFactors(dds) = estimateSizeFactorsForMatrix(counts)
  else if(normmethod == "scran") sizeFactors(dds) = scran::computeSumFactors(as.matrix(counts))
  else if(normmethod == "sf") sizeFactors(dds) = sf[sub]

  #' Filtering out low expressed genes
    if(!is.null(exprfilter)){
        mnc <- rowMeans(counts(dds, normalized=TRUE))
        dds <- dds[mnc > exprfilter,] 
        print(paste0("Genes left after filtering: ", sum(mnc > exprfilter)))
    }

  #' Running DEseq2
    dds = DESeq(dds, minReplicatesForReplace = Inf, betaPrior = shrinkage, parallel = parallel)

    DEres = results(dds,
                    contrast = contrast,
                    independentFiltering = independentFiltering,
                    alpha = padj_threshold,
                    parallel = parallel,
                    lfcThreshold = lfc_threshold)

    print(paste0("Coeficients:", resultsNames(dds)))
    print(summary(DEres))

    DEresCLEAN = DEres[!is.na(DEres$padj),]
    DEresCLEAN$isDE = DEresCLEAN$padj < padj_threshold

  #' Plotting the MA plot (to screen)
    DESeq2::plotMA(DEres, ylim = c(-2.5,2.5), alpha = padj_threshold)

  #' Calculating the shrunk fold changes (if shrinkage is set to FALSE
    if(!shrinkage & calcshrinkage){
        resshrink_normal <- lfcShrink(dds, contrast = contrast)

        resshrink_ashr <- lfcShrink(dds, contrast = contrast, type = "ashr")
        DESeq2::plotMA(resshrink_ashr, ylim=c(-2.5,2.5), main = "shrinkage_ashr")
    }

    #' Preparing the data.frame to export
    countsN = counts(dds, normalized=TRUE)
    DFout = as.data.frame(DEres)

    #' Adding the shrunken fold changes if shrinkage is set to FALSE
    if(!shrinkage & calcshrinkage){
        DFout$log2FoldChange_normshrink = resshrink_normal$log2FoldChange
        DFout$log2FoldChange_ashrshrink = resshrink_ashr$log2FoldChange
    }

    DFout$dispersion = dispersions(dds)
    DFout[,paste0(contrast[2], "Mean")] = rowMeans(countsN[,metadata[,contrast[1]] == contrast[2]])
    DFout[,paste0(contrast[3], "Mean")] = rowMeans(countsN[,metadata[,contrast[1]] == contrast[3]])
    DFout$symbol = genedata[match(row.names(DFout), genedata$geneid), "symbol"]
    DFout = data.frame(id = row.names(DFout), DFout, stringsAsFactors = FALSE)

    DFoutSIG = DFout[!is.na(DFout$padj),]
    DFoutSIG = DFoutSIG[DFoutSIG$padj < padj_threshold,]

  #' Plotting the report pdfs
    if(save){

    if(is.null(blocking)){
      pdf(paste0(filename, "_", contrast[1], "_", contrast[2], "_", contrast[3], ".pdf"))
    }

    else if(!is.null(blocking)){
      pdf(paste0(filename, "_", contrast[1], "_", contrast[2], "_", contrast[3], "_blocking_", blocking, ".pdf"))
    }

    plotDispEsts(dds)
    DESeq2::plotMA(DEres, ylim = c(-2.5,2.5), alpha = padj_threshold)

    if(!shrinkage & calcshrinkage){
      DESeq2::plotMA(resshrink_ashr, ylim=c(-2.5,2.5), main = "shrinkage_ashr")
    }

    hist(DEres$pvalue)
    
    volcano = ggplot(as.data.frame(DEresCLEAN), aes(x = log2FoldChange, y = -log10(pvalue), colour = isDE)) +
      geom_point(shape = 16)
    print(volcano)

    if(independentFiltering){
      plot(metadata(DEres)$filterNumRej, 
           type="b", ylab="number of rejections",
           xlab="quantiles of filter")
      lines(metadata(DEres)$lo.fit, col="red")
      abline(v=metadata(DEres)$filterTheta)
    }

    dev.off()
    write.csv(DFout, paste0(filename, "_", contrast[1], "_", contrast[2], "_", contrast[3], "_", normmethod, "_allgenes.csv"))
    write.csv(DFoutSIG, paste0(filename, "_", contrast[1], "_", contrast[2], "_", contrast[3], "_", normmethod, "_sigDE.csv"))

    #' Plotting violin plots of TOP genes
    if(nrow(DEresCLEAN) != 0) {
      DEresSORT = DEresCLEAN[order(DEresCLEAN$padj),]
      TOP = plotGeneExpr(log2(as.matrix(counts(dds, normalized = TRUE)+1)), metadata = metadata, genedata = genedata, genes = row.names(DEresSORT[1:15,]), xaxis = contrast[1], wrapby = "geneid", ncol = 3, freescale = TRUE)
      if (is.null(blocking)) {
        ggsave(paste0(filename, "_TOPgenes", "_", contrast[1], "_", contrast[2],
                      "_", contrast[3], ".pdf"), TOP, height = 25, width = 25)
      }
      else if (!is.null(blocking)) {
        ggsave(paste0(filename, "_TOPgenes", "_", contrast[1], "_", contrast[2],
                      "_", contrast[3], "_blocking_", blocking, ".pdf"), TOP, height = 25, width = 25)
      }
    }
  }
  
  if(output == "sig") return(DFoutSIG)
  else if (output == "all") return(DFout)
  else print("Please provide the output argument with 'sig' or 'all'")
}

#' Run quality control for single cell RNA-Seq data
#'
#' The function takes raw counts from an RNA-Seq experiments and performs a number of quality control checks, prints a full report and returns the counts table excluding low quality cells (the exact threshold can be adjusted accordingly)
#' @param counts dataframe with rows as genes and samples/cells as columns
#' @param metadata dataframe with rows as samples/cells and meta information stored in columns
#' @param mitogenes a vector of gene names indicating the mitochondrial genes
#' @param filebase basename for the report files generated, default to QC
#' @param tr_mapped treshold for filtering based on fraction of the reads mapped to the genome
#' @param tr_mincounts treshold for filtering based on the minimum number of reads per cell
#' @param tr_maxcounts treshold for filtering based on the maximum number of reads per cell
#' @param tr_erccfrac treshold for filtering based on the fraction of the reads mapped to the ERCCs
#' @param tr_higeneno treshold for filtering based on the number of genes exceeding 50 cpm
#' @param tr_mitofrac treshold for filtering based on the fraction of the reads mapped to the mitochondrial genes
#' @return A data.frame with counts after filtering out law quality cells
#' @export
scQC <- function(counts, metadata = NULL, mitogenes, filebase = "QCperplate", id = "SLX_index", splitby = "Plate",
                tr_mapped = 0.15,
                tr_mincounts = 100000,
                tr_maxcounts = 4000000,
                tr_erccfrac = 0.005,
                tr_higeneno = 2500,
                tr_mitofrac = 0.2){

  require(dplyr)
  theme_set(theme_bw())
  # Splitting the genes into categories (mouse genome, ERCCs, notaligned etc)
  muscounts = counts[grep("^EN.*", row.names(counts)),]
  ERCCcounts = counts[grep("^ERCC.*", row.names(counts)),]
  notcounts = counts[grep("^__.*", row.names(counts)),]

  # Calculating basic statistics
  statdf = data.frame(cellid = colnames(counts),
                      muscounts_frac = colSums(muscounts)/colSums(counts),
                      ercccounts_frac = colSums(ERCCcounts) / colSums(counts))

  statdf$mus_ERCC_prop = statdf$muscounts_frac / statdf$ercccounts_frac
  statdf$mus_countno = colSums(muscounts)

  #Handy function, which returns the number of genes above a certain expression threshold (in cpms)
  hiexpr  = function(x, tr = 10){
    x = x/(sum(x) / 1000000) #The treshold is in cpm
    return(sum(x>tr))
  }
  higeneno = apply(muscounts, 2, FUN = hiexpr, tr = 50)
  
  statdf$higeneno = higeneno
  statdf$mitoreads = colSums(muscounts[row.names(muscounts) %in% mitogenes,])
  statdf$mitofrac = statdf$mitoreads / statdf$mus_countno

  print(paste("Number of cells with less than 10% read mapped to mouse genome: ", sum(statdf$muscounts_frac < 0.1)))
  print(paste("Number of cells with less than 100 000 reads: ", sum(colSums(muscounts) < 100000)))

  pdf(paste(filebase, ".pdf", sep = ""))

  musfrac_hist= ggplot(statdf, aes(muscounts_frac)) +
    geom_histogram(colour = "black", fill = "white") +
    geom_vline(xintercept = tr_mapped, linetype="dashed", color = "red") +
    ggtitle("Fraction of reads mapped to mouse genome")
  print(musfrac_hist)

  counthist = ggplot(statdf, aes(log10(mus_countno))) +
    geom_histogram(colour = "black", fill = "white") +
    geom_vline(xintercept = log10(tr_mincounts), linetype="dashed", color = "red") +
    geom_vline(xintercept = log10(tr_maxcounts), linetype="dashed", color = "red") +
    ggtitle("No of counts mapped to the mouse genome")
  print(counthist)

  erccfrac_hist= ggplot(statdf, aes(ercccounts_frac)) +
    geom_histogram(colour = "black", fill = "white") +
    geom_vline(xintercept = tr_erccfrac, linetype="dashed", color = "red") +
    ggtitle("Fraction of reads mapped to ERCC")
  print(erccfrac_hist)

  higeneno= ggplot(statdf, aes(x = mus_countno, y = higeneno, colour = higeneno < tr_higeneno)) +
    geom_point(alpha = 0.8) +
    ggtitle("Highly expressed genes vs total counts") +
    theme(legend.position="none")
  print(higeneno)

  mito_hist= ggplot(statdf, aes(mitofrac)) +
    geom_histogram(colour = "black", fill = "white") +
    geom_vline(xintercept = tr_mitofrac, linetype="dashed", color = "red") +
    ggtitle("Reads mapped to mitochondrial genes/mapped genes")
  print(mito_hist)


  if(!is.null(metadata)){
    statdf$splitset = metadata[match(statdf$cellid, metadata[,id]), splitby]
    print(head(statdf))

    musfrac_hist2= ggplot(statdf, aes(x = muscounts_frac, fill = splitset)) +
      geom_histogram(alpha = 0.75, position="identity") +
      facet_wrap(~splitset) + 
      geom_vline(xintercept = tr_mapped, linetype="dashed", color = "red") +
      ggtitle("Fraction of reads mapped to mouse genome")
    print(musfrac_hist2)

    counthist2 = ggplot(statdf, aes(log10(mus_countno), fill = splitset)) +
      geom_histogram(alpha = 0.75, position = "identity") +
      facet_wrap(~splitset) +
      geom_vline(xintercept = log10(tr_mincounts), linetype="dashed", color = "red") +
      geom_vline(xintercept = log10(tr_maxcounts), linetype="dashed", color = "red") +
      ggtitle("No of counts mapped to the mouse genome")
    print(counthist2)
    
    erccfrac_hist2= ggplot(statdf, aes(ercccounts_frac, fill = splitset)) +
      geom_histogram(alpha = 0.75, position = "identity") +
      facet_wrap(~splitset) +
      geom_vline(xintercept = tr_erccfrac, linetype="dashed", color = "red") +
      ggtitle("Fraction of reads mapped to ERCC")
    print(erccfrac_hist2)

    mitofrac_hist2= ggplot(statdf, aes(mitofrac, fill = splitset)) +
      geom_histogram(alpha = 0.75, position = "identity") +
      facet_wrap(~splitset) +
      geom_vline(xintercept = tr_mitofrac, linetype="dashed", color = "red") +
      ggtitle("Reads mapped to mitochondrial genes/mapped genes")
    print(mitofrac_hist2)

  higeneno2= ggplot(statdf, aes(x = mus_countno, y = higeneno, fill = splitset, colour = higeneno < tr_higeneno)) +
    geom_point(alpha = 0.8) +
    facet_wrap(~splitset, scales = "free") +
    ggtitle("Highly expressed genes vs total counts") +
    theme(legend.position="none")
    print(higeneno2)

    higeneno2S= ggplot(statdf, aes(x = mus_countno, y = higeneno, fill = splitset, colour = higeneno < tr_higeneno)) +
    geom_point(alpha = 0.8) +
    facet_wrap(~splitset) +
    ggtitle("Highly expressed genes vs total counts (same scale)") +
    theme(legend.position="none")
    print(higeneno2S)
  }

  dev.off()
  
  countsQC = counts[,(statdf$muscounts_frac > tr_mapped) & (statdf$mus_countno > tr_mincounts) & (statdf$mus_countno < tr_maxcounts) & (statdf$ercccounts_frac < tr_erccfrac) & (statdf$higeneno > tr_higeneno) & (statdf$mitofrac < tr_mitofrac)]
  
  sumdf = data.frame(cellname = statdf$cellid,
                     splitset = statdf$splitset,
                     TotalReads = colSums(counts),
                     MappedReads = statdf$mus_countno,
                     MappedFraction = statdf$muscounts_frac,
                     ERCCfraction = statdf$ercccounts_frac,
                     higenes = statdf$higeneno,
                     MitoFrac = statdf$mitofrac,
                     QCpass = statdf$cellid %in% colnames(countsQC))
   
  sumdf_bysplitset = sumdf %>%
    group_by(splitset) %>%
    summarise(TotalReads = sum(TotalReads),
              MappedReads = sum(MappedReads),
              MappedFractionAV = mean(MappedFraction),
              ERCCfractionAV = mean(ERCCfraction),
              higenesAV = mean(higenes),
              MitoFrac = mean(MitoFrac),
              QCpass = sum(QCpass))

  print(as.data.frame(sumdf_bysplitset))

  write.csv(sumdf, paste(filebase, "percell.csv", sep = ""))

  write.csv(sumdf_bysplitset, paste(filebase, "perset.csv", sep = ""))

  return(countsQC)

}
#' Plot multiple ggplots function
#'
#' Function which takes multiple ggpltos and plots them in a single window. Ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects).
#' #' If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
#' then plot 1 will go in the upper left, 2 will go in the upper right, and
#' 3 will go all the way across the bottom.
#' Source: Cookbook R website: http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#' @param cols:   Number of columns in layout
#' @param layout: A matrix specifying the layout. If present, 'cols' is ignored.
#' @export
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#' Function which adds the linear model equation to the geom_smooth(method = lm) call in ggplot2
#' The function was created by KDauria, soure https://gist.github.com/kdauria/524eade46135f6348140
#' @param mapping 
#' @param data 
#' @param geom 
#' @param position 
#' @param ... 
#' @param method 
#' @param formula 
#' @param se 
#' @param n 
#' @param span 
#' @param fullrange 
#' @param level 
#' @param method.args 
#' @param na.rm 
#' @param show.legend 
#' @param inherit.aes 
#' @param xpos 
#' @param ypos 
#' @export
stat_smooth_func <- function(mapping = NULL, data = NULL,
                             geom = "smooth", position = "identity",
                             ...,
                             method = "auto",
                             formula = y ~ x,
                             se = TRUE,
                             n = 80,
                             span = 0.75,
                             fullrange = FALSE,
                             level = 0.95,
                             method.args = list(),
                             na.rm = FALSE,
                             show.legend = NA,
                             inherit.aes = TRUE,
                             xpos = NULL,
                             ypos = NULL) {
  layer(
    data = data,
    mapping = mapping,
    stat = StatSmoothFunc,
    geom = geom,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      method = method,
      formula = formula,
      se = se,
      n = n,
      fullrange = fullrange,
      level = level,
      na.rm = na.rm,
      method.args = method.args,
      span = span,
      xpos = xpos,
      ypos = ypos,
      ...
    )
  )
}
StatSmoothFunc <- ggproto("StatSmooth", Stat,
                          
                          setup_params = function(data, params) {
                                        # Figure out what type of smoothing to do: loess for small datasets,
                                        # gam with a cubic regression basis for large data
                                        # This is based on the size of the _largest_ group.
                            if (identical(params$method, "auto")) {
                              max_group <- max(table(data$group))
                              
                              if (max_group < 1000) {
                                params$method <- "loess"
                              } else {
                                params$method <- "gam"
                                params$formula <- y ~ s(x, bs = "cs")
                              }
                            }
                            if (identical(params$method, "gam")) {
                              params$method <- mgcv::gam
                            }
                            
                            params
                          },
                          
                          compute_group = function(data, scales, method = "auto", formula = y~x,
                                                   se = TRUE, n = 80, span = 0.75, fullrange = FALSE,
                                                   xseq = NULL, level = 0.95, method.args = list(),
                                                   na.rm = FALSE, xpos=NULL, ypos=NULL) {
                            if (length(unique(data$x)) < 2) {
                                        # Not enough data to perform fit
                              return(data.frame())
                            }
                            
                            if (is.null(data$weight)) data$weight <- 1
                            
                            if (is.null(xseq)) {
                              if (is.integer(data$x)) {
                                if (fullrange) {
                                  xseq <- scales$x$dimension()
                                } else {
                                  xseq <- sort(unique(data$x))
                                }
                              } else {
                                if (fullrange) {
                                  range <- scales$x$dimension()
                                } else {
                                  range <- range(data$x, na.rm = TRUE)
                                }
                                xseq <- seq(range[1], range[2], length.out = n)
                              }
                            }
                                        # Special case span because it's the most commonly used model argument
                            if (identical(method, "loess")) {
                              method.args$span <- span
                            }
                            
                            if (is.character(method)) method <- match.fun(method)
                            
                            base.args <- list(quote(formula), data = quote(data), weights = quote(weight))
                            model <- do.call(method, c(base.args, method.args))
                            
                            m = model
                            eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                                             list(a = format(coef(m)[1], digits = 3), 
                                                  b = format(coef(m)[2], digits = 3), 
                                                  r2 = format(summary(m)$r.squared, digits = 3)))
                            func_string = as.character(as.expression(eq))
                            
                            if(is.null(xpos)) xpos = min(data$x)*0.9
                            if(is.null(ypos)) ypos = max(data$y)*0.9
                            data.frame(x=xpos, y=ypos, label=func_string)
                            
                          },
                          
                          required_aes = c("x", "y")
                          )
#' Prepare a gene data table
#'
#' The function pulls from the Biomart/Ensembl information about genes. The version/release is specified in the argument. In current version the function works only for mouse genes.
#' NOTE: make sure to sort the gene names according to the data
#' @param release specify which release of the ensembl should be used
#' @param species either "human" or "mouse"
#' @return A data.frame rows as genes and additional information such as gene symbols etc in columns
#' @export
makeGENEDATA <- function(release = "jul2015.archive.ensembl.org", species = "mouse"){
#"jul2015.archive.ensembl.org" is Ensembl 81 version

    require(biomaRt)
    if(species == "mouse"){
        ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host=release)
        ensembl <- useDataset("mmusculus_gene_ensembl", ensembl)
        genedata <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                          mart = ensembl)
        colnames(genedata) <- c("geneid", "symbol", "chromosome")
        genedata <- genedata[!duplicated(genedata$geneid),]
        genedata$is_mito = genedata$chromosome == "MT"
        row.names(genedata) = genedata$geneid
        genedata = genedata[order(genedata$geneid, decreasing = FALSE),]
    }
    else if(species == "human"){
        ensembl <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", host=release)
        ensembl <- useDataset("hsapiens_gene_ensembl", ensembl)
        genedata <- getBM(attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name"),
                          mart = ensembl)
        colnames(genedata) <- c("geneid", "symbol", "chromosome")
        genedata <- genedata[!duplicated(genedata$geneid),]
        genedata$is_mito = genedata$chromosome == "MT"
        row.names(genedata) = genedata$geneid
        genedata = genedata[order(genedata$geneid, decreasing = FALSE),]


    }
    else print("species needs to be either human or mouse")

  return(genedata)
}

##' Plot expression for selected genes (violin + jitter plot)
##'
##' This function takes in the gene counts + accessory information in form of metadata and genedata data.frames and plots expression for selected genes. Additional graphical parameters can be passed onto the functio to specify, which factor should be put on the x axis or plotted as colours or plot wrapped.
##' @param counts a dataframe with rows as genes and samples/cells as columns
##' @param metadata metadata a data.frame containing columns with sets and subsets arguments
##' @param genedata a data.frame containing information for genes (gene ids need to be specified in the id column)
##' @param genes a vector with names of genes to plot (the values are checked in the id column or in the symbol column of the genedata)
##' @param xaxis which factor should be plotted as x axis
##' @param wrapby according to which factor the the plot should be wrapped (default checked in metadata first)
##' @param colourby according to which factor the the plot should be colourcoded (default checked in metadata first)
##' @param ncol indicates the number of columns when wrapping the plot
##' @param freescale indicate if the scales should be free when wrapping the plot
##' @return returns a ggplot2 object ready to plot
##' @author idk25
##' @export
plotGeneExpr = function(counts, metadata, genedata, genes, xaxis = "symbol", colourby = NULL, wrapby = NULL, ncol = 3, freescale = FALSE){
  require("dplyr")
  require("reshape2")
  theme_set(theme_bw())


  print("Make sure that the metadata, counts and genedata data frames have the same order of elements!")
  if(grepl("ENS", genes[1])){
    genecounts = counts[genes, , drop = FALSE]
  }
  else{
    geneids = genedata[match(genes, genedata$symbol), "geneid"]
    genecounts = counts[geneids, , drop = FALSE]
  }

  dftoplot = reshape2::melt(genecounts)
  colnames(dftoplot) = c("geneid", "cellid", "expr")

  dftoplot$symbol = genedata[match(dftoplot$geneid, genedata$geneid), "symbol"]

  if(xaxis != "symbol" & xaxis != "geneid") dftoplot[,xaxis] = metadata[match(dftoplot$cellid, metadata$cellid), xaxis]
  

  geneplot = ggplot(dftoplot, aes_string(y = "expr", x = xaxis))

  if(!is.null(colourby)){
    if(colourby %in% colnames(metadata)) dftoplot$colourfactor = metadata[match(dftoplot$cellid, metadata$cellid), colourby]
    else dftoplot$colourfactor = dftoplot[,colourby]

    geneplot = geneplot + geom_violin(aes(colour = dftoplot$colourfactor), draw_quantiles = c(0.5)) +
      geom_jitter(aes(colour = dftoplot$colourfactor), alpha = 0.6, shape = 16, width = 0.25)
  }
  else{
    geneplot = geneplot + geom_violin(draw_quantiles = c(0.5)) +
      geom_jitter(alpha = 0.6, shape = 16, width = 0.25)
  }

  if(!is.null(wrapby)){
    if(wrapby %in% colnames(metadata)) dftoplot$wrapfactor = metadata[match(dftoplot$cellid, metadata$cellid), wrapby]
    else dftoplot$wrapfactor = dftoplot[,wrapby]

    if(freescale) geneplot = geneplot + facet_wrap(~dftoplot$wrapfactor, ncol = ncol, scales = "free")
    else geneplot = geneplot + facet_wrap(~dftoplot$wrapfactor, ncol = ncol)
  }


  #print(head(dftoplot))

  return(geneplot)
}


#' Recursively merges data frame, simply supply a list of data frames, specify the column and select whether you want a union or an intersectoin (all=FALSE for intersection)
#' @param dflist a list of data.frames
#' @param by specifies the name of the column which is used as id for the merge
#' @param all logical, TRUE outputs a union, FALSE outputs an intersection
#' @return outputs the merged data.frame (a single column of by, and all the remaining columns)
#' @export

mergeDFlist = function(dflist, by = "id", all = FALSE){
    a = dflist[[1]]
    for (i in 2:length(dflist)){
        a = merge(a, dflist[[i]], by = by, all = all)
    }
    return(a)
}

#' Select highly variable genes, based on the Brennecke et al. 2013 method
#'
#' The function takes not-normalised counts (all genes included, mouse genes, ERCCs etc), normalises the samples/cells internally (DESeq2 method) and fits function to CV^2 vs expression level based on the ERCC spike in data. Following that it tests for genes which variation exceedes variation inferred from ERCC spike ins. The function returns the raw counts for the selected genes (mouse genes only) and the associated adjusted pvalue.
#' @export
#' @param counts dataframe with rows as genes and samples/cells as columns (recommended after QC)
#' @param minBiolVar parameter of minimal biological variance for the statistical tests
#' @param cv2tr minimal threshold for the CV^2 variance for ERCCs
#' @param meanquantile threshold for minimal expression level for ERCCs (excluding too low expressed)
#' @param padjtr threshold for adjusted p value
#' @param plotfile name of the file to plot the diagnostic plots
#' @return a data.frame with raw counts for highly variable genes (only the mouse genes) and associated adj p-value
hivar = function(counts,
                 minBiolVar = 0.5^2,
                 cv2tr = 0.3,
                 meanquantile = 0.8,
                 padjtr = 0.05,
                 plotfile = "highvar_plots.pdf"){
  require("DESeq2")
  require("statmod")
  require("genefilter")
  
  #' Separating the mouse and ERCC rows
  countsMmus = counts[grep("^EN.*", row.names(counts)),]
  countsERCC = counts[grep("^ERCC.*", row.names(counts)),]

  #' Getting the sizefactors from the DESeq package for each dataset
  sfMmus <- estimateSizeFactorsForMatrix( countsMmus )
  sfERCC <- estimateSizeFactorsForMatrix( countsERCC )

  #' Normalising the data
  nCountsERCC <- t( t(countsERCC) / sfERCC )
  nCountsMmus <- t( t(countsMmus) / sfMmus )

  #' Calculating the means, variances and cv2 for each dataset
  meansERCC <- rowMeans( nCountsERCC )
  varsERCC <- rowVars( nCountsERCC )
  cv2ERCC <- varsERCC / meansERCC^2
  meansMmus <- rowMeans( nCountsMmus )
  varsMmus <- rowVars( nCountsMmus )
  cv2Mmus <- varsMmus / meansMmus^2

  #' Excluding genes with too low mean, which could skew the fit
  minMeanForFitA <- unname( quantile( meansERCC[ which( cv2ERCC > cv2tr ) ], meanquantile ) )
  useForFitA <- meansERCC >= minMeanForFitA
  print("The minimum mean used for fit:")
  print(minMeanForFitA) #The minimum mean used
  print("The number of genes used for fit:)")
  print(table( useForFitA )) #How many genes are used

  #' Fitting E w ≈ (Xi + atilde)/mean + alpha0
  #' where Xi= mean( 1 / sizefactors )
  fitA <- glmgam.fit( cbind( a0 = 1, a1tilde = 1/meansERCC[useForFitA] ),
                     cv2ERCC[useForFitA] )
  fitA$coefficients
  #' The intercept is the estimator for alpha0, the coefficient for 1/mean is the estimator for a1 = Xi+atilde
  
  #' Looking at residuals
  residualA <- var( log( fitted.values(fitA) ) - log( cv2ERCC[useForFitA] ) )
  totalA <- var( log( cv2ERCC[useForFitA] ) )
  #' How much variance is explained?
  print("Proportion of variance explained:")
  print(c( A = 1 - residualA / totalA))
  
  #' Plotting the ERCC cv2 vs means
  pdf(plotfile)
  
  plot( meansERCC, cv2ERCC, log="xy", col=1+useForFitA, main="A" )
  xg <- 10^seq( -3, 5, length.out=100 )
  lines( xg, coefficients(fitA)["a0"] + coefficients(fitA)["a1tilde"]/xg )
  segments( meansERCC[useForFitA], cv2ERCC[useForFitA],
           meansERCC[useForFitA], fitA$fitted.values, col="gray" )
                                        #text(meansERCC, cv2ERCC, labels=names(meansERCC), cex= 0.7)
  
  #' Performing the statistical test
  minBiolDisp <- minBiolVar #Setting the minimum biological variance to 0.5
  
  xi <- mean( 1 / sfERCC )
  m <- ncol(countsMmus)
  
  psia1thetaA <- mean( 1 / sfERCC ) +
    ( coefficients(fitA)["a1tilde"] - xi ) * mean( sfERCC / sfMmus )
  cv2thA <- coefficients(fitA)["a0"] + minBiolDisp + coefficients(fitA)["a0"] * minBiolDisp
  
  #' Denominator for the test
  testDenomA <- ( meansMmus * psia1thetaA + meansMmus^2 * cv2thA ) / ( 1 + cv2thA/m )
  #' ChiSquare test:
  pA <- 1 - pchisq( varsMmus * (m-1) / testDenomA, m-1 )
  
  #' Adjusting the p-values with Benjamini-Hochberg
  padjA <- p.adjust( pA, "BH" )
  #' How many genes passed
  print(paste("The number of genes with padj<", padjtr ,sep = ""))
  print(table( padjA < padjtr ))
  
  #' Plotting the fit and genes
  plot( NULL, xaxt="n", yaxt="n",
       log="xy", xlim = c( 1e-1, 3e5 ), ylim = c( .005, 100 ),
       xlab = "average normalized read count", ylab = "squared coefficient of variation (CV^2)")
  axis( 1, 10^(-1:5), c( "0.1", "1", "10", "100", "1000",
                        expression(10^4), expression(10^5) ))
  axis( 2, 10^(-2:2), c( "0.01", "0.1", "1", "10" ,"100"), las=2 )
  abline( h=10^(-2:1), v=10^(-1:5), col="#D0D0D0", lwd=2 )
  #' Plot the plant genes, use a different color if they are highly variable
  points( meansMmus, cv2Mmus, pch=20, cex=.2,
         col = ifelse( padjA < padjtr, "#C0007090", "#70500040" ) )
  #' Add the technical noise fit, as before
  xg <- 10^seq( -2, 6, length.out=1000 )
  lines( xg, coefficients(fitA)["a1tilde"] / xg + coefficients(fitA)["a0"], col="#FF000080", lwd=3 )
  #' Add a curve showing the expectation for the chosen biological CV^2 thershold
  lines( xg, psia1thetaA/xg + coefficients(fitA)["a0"] + minBiolDisp,
        lty="dashed", col="#C0007090", lwd=3 )
  #' Add the normalised ERCC points
  points( meansERCC, cv2ERCC, pch=20, cex=1, col="#0060B8A0" )

  dev.off()

  dfout = cbind(countsMmus, padjA)
  dfout = dfout[!is.na(dfout$padjA),]
  dfout = dfout[dfout$padjA < padjtr,]
  return(dfout)
  
}

#' Plot t-SNE projections for scRNA-Seq data (report type)
#'
#' This function takes in a counttable (QCed and normalised), calculates the tSNE coordinates and plots the projection. Additional argument-sets plots the plot is colourcoded according to the sets factor levels. Additional argument subsets recalculates the tSNE  and plots it for an indicated subset of the data, colourcoding according to the specified arguments. The dimensions to be plotted can be specified in ndims argument.
#' @export
#' @param counts a dataframe with rows as genes and samples/cells as columns
#' @param metadata a data.frame containing columns with sets and subsets arguments
#' @param sets vector with names of columns in the metadata, which indicates the colourcoding scheme for the plots
#' @param subsets a list of vectors, where first element of each vector is a name of the column in the metadata with factors of interest, and the following arguments are levels of the factors that will be used for tSNE calculation and plotting
#' @param pdffile base file name for the plots
#' @param dims the number of dimension to be calculated by tSNE (default to 2)
#' @param ndims a 2-element vector indicating which dimensions of the tSNE should be plotted
#' @param id of the column in the metadata with cell/sample ids (default to SLX_index)
#' @param perplexity parameter for tSNE calculation, see the tSNE documentation
#' @param seed number to ensure the exact same plot is generated each time
#' @return data.frame with tSNE coordinates
TSNEproj = function(counts, metadata = NULL, sets = NULL, subsets = NULL,
                    pdffile = "tsne.pdf",
                    dims = 2,
                    ndims = c(1,2),
                    id = "SLX_index",
                    perplexity = 30,
                    seed = 123){
  require("dplyr")
  require("Rtsne")
  require("reshape2")
  theme_set(theme_bw())

  #' Settting seed to ensure reproducibility of the projections
  print("!!!! SETTING THE SEED GLOBALLY !!!!! (to remove: rm(.Random.seed, envir=.GlobalEnv)")
  set.seed(seed)

  #' Calculating the tSNE for the whole data
  tsne_out <- Rtsne(as.matrix(t(counts)), dims = dims, perplexity = perplexity) # Run TSNE
  tsneDF = data.frame(x = tsne_out$Y[,ndims[1]], y = tsne_out$Y[,ndims[2]])
  Xname = paste0("tsne", ndims[1])
  Yname = paste0("tsne", ndims[2])
  colnames(tsneDF) = c(Xname, Yname)
  tsneDF$SLX_i = colnames(counts)

  pdf(pdffile)

  #' Plotting the tSNE for the entire set
  g = ggplot(tsneDF, aes_string(x = Xname, y = Yname)) +
    geom_point()
  print(g)

  #' Plotting the versions of the above plot colourcoded for the factors indicates in sets
  if(!is.null(sets)){
    for(i in sets){
      tsneDF[,i] = metadata[match(tsneDF$SLX_i, metadata[,id]), i]
      print(head(tsneDF))
      
      g = ggplot(tsneDF, aes_string(Xname, Yname, colour = i)) +
        geom_point(size = 1.7, alpha = 0.7) +
        theme(legend.position="bottom")
      print(g)
    }
  }
#' Calculating and plotting tSNE for subsets of the data 
  if(!is.null(subsets)){
    for(i in subsets){

      chosencells = metadata[metadata[,i[1]] %in% i[-1], id]
      countsub = counts[,colnames(counts) %in% chosencells]

      if(ncol(countsub) > 120) tsne_out <- Rtsne(as.matrix(t(countsub)), dims = 2, perplexity = 30) # Run TSNE
      else tsne_out <- Rtsne(as.matrix(t(countsub)), dims = 3, perplexity = 12, theta = 0.5) # Run TSNE
      tsneDFsub = data.frame(x = tsne_out$Y[,1], y = tsne_out$Y[,2])
      colnames(tsneDFsub) = c(Xname, Yname)
      tsneDFsub$SLX_i = colnames(countsub)
      tsneDFsub[,i[1]] = metadata[match(tsneDFsub$SLX_i, metadata[,id]), i[1]]
      print(head(tsneDFsub))
      
      g = ggplot(tsneDFsub, aes_string(Xname, Yname, colour = i[1])) +
        geom_point(size = 1.7, alpha = 0.7) +
        theme(legend.position="bottom")
      print(g)
    }
  }

  dev.off()

  return(tsneDF)
}



#' Plot PCA projections for scRNA-Seq data (report type)
#'
#' This function takes in a counttable (QCed and normalised), calculates the PCA coordinates and plots the component variance and the projection. Additional argument-sets plots the PCA colourcoded according to the sets factor levels. Additional argument subsets recalculates the PCA and plots it for an indicated subset of the data, colourcoding according to the specified arguments. The dimensions to be plotted can be specified in ndims argument.
#' @export
#' @param counts a dataframe with rows as genes and samples/cells as columns
#' @param metadata a data.frame containing columns with sets and subsets argumenta
#' @param sets vector with names of columns in the metadata, which indicates the colourcoding scheme for the PCA plots
#' @param subsets a list of vectors, where first element of each vector is a name of the column in the metadata with factors of interest, and the following arguments are levels of the factors that will be used for PCA calculation and plotting
#' @param ndims a 2-element vector indicating which dimensions of the PCA should be plotted
#' @param id of the column in the metadata with cell/sample ids (default to SLX_index)
#' @param scaled a logical specifyin if the values should be scaled (by using the scale function, which calculates essentially a zscore)

#' @return data.frame with PCA cooordinates
PCAproj = function(counts,
                   metadata,
                   sets = NULL,
                   subsets = NULL,
                   pdffile = "PCA.pdf",
                   ndims = c(1,2),
                   id = "SLX_index",
                   scaled = FALSE){
  require("dplyr")
  require("reshape2")
  require("statmod")
  theme_set(theme_bw())

  #' Defining names for the PC dimensions
  PCs = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")

  pdf(pdffile)
  
  #'Calculating the PCA for the entire dataset
  pca = prcomp(t(counts), scale. = scaled)
  psum = summary(pca)
  #'Plotting the variance proportion
  Variance_proportion = psum$importance[2,1:10]
  plot(Variance_proportion, type = "b")

  #' Plotting the PCA for the entire set
  pcadata = as.data.frame(pca$x)
  g = ggplot(pcadata, aes(x = PC1, y = PC2)) +
    geom_point()
  print(g)

  pcadata$SLX_i = colnames(counts)

  #' Plotting the versions of the above plot colourcoded for the factors indicates in sets
  if(!is.null(sets)){
    for(i in sets){
      pcadata[,i] = metadata[match(pcadata$SLX_i, metadata[,id]), i]

      g = ggplot(pcadata, aes_string(PCs[ndims[1]], PCs[ndims[2]], colour = i)) +
        geom_point(size = 1.7, alpha = 0.7) +
        theme(legend.position="bottom")
      print(g)
    }
  }

  #' Calculating and plotting PCA for subsets of the data 
  if(!is.null(subsets)){
    for(i in subsets){

      chosencells = metadata[metadata[,i[1]] %in% i[-1], id]
      countsub = counts[,colnames(counts) %in% chosencells]

      pca = prcomp(t(countsub), scale. = scaled)
      pcadata = as.data.frame(pca$x)
      
      pcadata$SLX_i = colnames(countsub)
      pcadata[,i[1]] = metadata[match(pcadata$SLX_i, metadata[,id]), i[1]]
      
      g = ggplot(pcadata, aes_string(PCs[ndims[1]], PCs[ndims[2]], colour = i[1])) +
        geom_point(size = 1.7, alpha = 0.7) +
        theme(legend.position="bottom")
      print(g)
    }
  }
  dev.off()
  return(pcadata)
}



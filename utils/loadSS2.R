
clean_SLXnames = function(x){
  x = gsub("-", "_", x)
  x = gsub("SLX_", "SLX", x)
  return(x)
}

##' Load and concatenate fcounts
##'
##' Function for loading all fcount files in a given directory. Data is loaded using a very fast fread function and subsequently converted to data.frame and all data is concatenated by columns
##' @title Load and concatenate fcounts
##' @param path string, path to the fcount files
##' @return data.frame with concatenated columns of each fcount table
##' @author idk25
get_counts = function(path){
    require(data.table)

    #' Loading the counts
    x = list.files(path, ".*fcounts\\.txt$", full.names = TRUE)
    print(x)

    allcounts = list()
    for (i in x){
        counts = fread(i)
        counts2 = as.data.frame(counts, stringsAsFactors = FALSE)
        row.names(counts2) = counts2$Geneid
        counts2 = counts2[,grepl(".*bam", colnames(counts2))]
        counts2 = counts2[order(row.names(counts2)),]
        print(dim(counts2))
        ## print(counts2[1:5,1:2])

        allcounts[[i]] = counts2
    }
    allcounts = do.call(cbind, allcounts)
    allcounts = allcounts[, !grepl("lostreads", colnames(allcounts))]

    #Substituting all dashes for dots, so that they are valid column names in R and getting the sample names
    print(colnames(allcounts)[1:5])
    colnames(allcounts) = clean_SLXnames(colnames(allcounts))
    colnames(allcounts) = gsub("(.*STARout\\/)(SLX[0-9][0-9][0-9][0-9][0-9]\\.i7[0-9][0-9]_i5[0-9][0-9])(.*)", "\\2", colnames(allcounts))
    print(colnames(allcounts)[1:5])

   #' Loading the stats
    y = list.files(path, ".*fcounts\\.txt.summary", full.names = TRUE)
    print(y)

    allstats = list()
    for (i in y){
        stats = fread(i)
        stats2 = as.data.frame(stats, stringsAsFactors = FALSE)
        row.names(stats2) = stats2$Status
        stats2 = stats2[,grepl(".*bam", colnames(stats2))]

        allstats[[i]] = stats2
    }
    allstats = do.call(cbind, allstats)
    allstats = allstats[, !grepl("lostreads", colnames(allstats))]

    colnames(allstats) = clean_SLXnames(colnames(allstats))
    colnames(allstats) = gsub("(.*STARout\\/)(SLX[0-9][0-9][0-9][0-9][0-9]\\.i7[0-9][0-9]_i5[0-9][0-9])(.*)", "\\2", colnames(allstats))

    allstats = as.data.frame(t(allstats))
    print(allstats[1:5,1:5])

    htseq_stats = data.frame(row.names = row.names(allstats),
                             no_feature = allstats$Unassigned_NoFeatures,
                             ambiguous = allstats$Unassigned_Ambiguity,
                             too_low_aQual = allstats$Unassigned_MappingQuality,
                             not_aligned = allstats$Unassigned_Unmapped,
                             alignment_not_unique = allstats$Unassigned_MultiMapping,
                             stringsAsFactors = FALSE)
    colnames(htseq_stats) = paste0("__", colnames(htseq_stats))

    htseq_stats = as.data.frame(t(htseq_stats))
    htseq_stats = htseq_stats[, colnames(allcounts)]
    print(dim(htseq_stats))

    allcounts = rbind(allcounts, htseq_stats)

    return(allcounts)
}

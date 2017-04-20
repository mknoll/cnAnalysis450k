#' @title 
#' Create overview for all samples of calculated Bins
#' 
#' @description
#' Creates a aligned matrix of all calculated bins with
#' samples as columns and rows as ordered bins.
#'
#' @param data list created with createBinsFast()
#' containing CN values and p-values for the calculated
#' bins
#'
#' @param pval pvalue for which bins are filtered. 
#' If one bin has a p-value in any of the samples,
#' the bin is retained
#'
#' @return 
#' Matrix with samples as columns and CN Bin values 
#' as rows
#'
#' @import plyr
#' 
#' @export
#' 
#' @examples
createBinMatrix <- function(data, pval) {    
    ## use parallelization
    no_cores <- parallel::detectCores() - 1
    no_cores <- ifelse(no_cores == 0, 1, no_cores)
    doParallel::registerDoParallel(no_cores)
    
    out <- foreach (i=1:length(data)) %dopar% {
        ##select rows for which p <= pval
        sel <- apply(data[[i]]$p.val, 1, function(x) any(x <= pval))
        subD <- data[[i]]$median[sel,]
        rownames(subD) <- paste(data[[i]]$chr, 
                                data[[i]]$startCGs[sel],
                                rownames(data[[i]]$p.val)[sel], 
                                sep=":")
        subD
    }
    
    stopImplicitCluster()
    
    return (do.call(rbind, out))
}

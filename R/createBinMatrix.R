#' @title Create Matrix from calculated Bins
#' 
#' @import plyr
#' 
#' @export
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
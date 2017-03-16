require(plyr)

#' @title Binarize Data
#'
#' @description
#' Binarizes Data (-1, 0, +1) according to given cutoff list
#'
#' @param data inputdata
#' @param cutoffs cutoffs, calculated e.g. with \code{findCutoffs()}
#' @param effectsize minimal required effectsize as vector (loss, gain)
#'
#' @import utils 
#'
#' @export
#'
#' @return binarized values
#'
#' @examples
#' candMatrix <- data.frame(
#' smp1=c(-0.097, -1.208,-0.134, 1.732),
#' smp2=c(-0.006, 0.004, 0.004, -0.001),
#' smp3=c(0.050, 0.008, 0.008,0.046 ))
#' rownames(candMatrix) <- c(
#' "chr1:15865", "chr1:110230252",
#' "chr1:110254692", "chr3:45838226"
#' )
#' cutoffs <- data.frame(
#' rn=c(
#' "chr1:15865", "chr1:110230252",
#' "chr1:110254692", "chr3:45838226"
#' ),
#' i=c(1,2,3,1),
#' cutoffLoss=c(-0.85,NA,NA,-0.05),
#' cutoffGain=c(NA,NA,NA,NA),
#' cutoffLossWP=c(-0.80, -0.83, -0.83, NA), 
#' cutoffGainWP=c(NA,NA,NA,NA),
#' baseline=c(0.01,0.01,0.01,0.01)
#' )
#' segmentData(candMatrix, cutoffs)
segmentData <- function(data, cutoffs, effectsize = c(0, 0)) {
    print("Binarize data ... ")
    completeDATA <- data
    signCluster <- cutoffs
    
    candDataLG <-
        completeDATA[
            which(rownames(completeDATA) %in% signCluster$rn), ,drop=FALSE]    
    ind <- match(rownames(candDataLG), signCluster$rn)
    
    total <- length(candDataLG[, 1])
    for (i in 1:total) {
        tm <- Sys.time()
        ## progress
        if (i %% 10 == 0) { 
            utils::flush.console() 
            cat("\r",round(i/total*100,2),"%  ") 
        }
        
        ## gains / losses
        #idx <- which(signCluster$rn == rownames(candDataLG)[i])
        idx <- ind[i]
        
        ## Baseline
        base <- signCluster[idx, 7] #"cutoffGainWP"]
        
        ## fallback auf wp cutoff, wenn keine echten maxima identifiziert 
        #werden konnten
        cutG <- ifelse(is.na(signCluster[idx, 4]), 
                    signCluster[idx, 6],signCluster[idx, 4])
        cutL <- ifelse(is.na(signCluster[idx, 3]), 
                    signCluster[idx, 5],signCluster[idx, 3])
        
        ## minimale effectgroesse ueberschritten?
        diff <- abs(base - cutL)
        if (!is.na(diff) && diff < effectsize[1]) { cutL <- NA } 
        if (!is.na(diff) && diff < effectsize[2]) { cutG <- NA }
        
        candDataLG[i, ] <-
            ifelse(!is.na(cutG) & candDataLG[i, ] > cutG, 1, candDataLG[i, ])
        candDataLG[i, ] <-
            ifelse(!is.na(cutL) & candDataLG[i, ] < cutL,-1, candDataLG[i, ])
        tm <- c(tm, Sys.time())
    }
    cat("\n")
    candDataLG <- apply(candDataLG, 2, 
                        function(x) ifelse (x == 1 | x == -1, x, 0))
    
    return(candDataLG)
}

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
#' segmentsData(candMatrix, cutoffs)
segmentData <- function(data, cutoffs, effectsize = c(0, 0)) {
    print("Binarize data ... ")
    completeDATA <- data
    signCluster <- cutoffs
    
    candDataLG <-
        completeDATA[which(rownames(completeDATA) %in% signCluster$rn), ,drop=F]
    gainsLosses <- NULL
    
    total <- length(candDataLG[, 1])
    for (i in 1:total) {
        ## progress
        if (i %% 10 == 0) { flush.console(); cat("\r",round(i/total*100,2),"%  ") }
        
        ## gains / losses
        cutG <-
            signCluster[which(signCluster$rn == rownames(candDataLG)[i]), 
                        "cutoffGain"]
        cutL <-
            signCluster[which(signCluster$rn == rownames(candDataLG)[i]), 
                        "cutoffLoss"]
        ## Wendepunkte?
        cutGWP <-
            signCluster[which(signCluster$rn == rownames(candDataLG)[i]), 
                        "cutoffGainWP"]
        cutLWP <-
            signCluster[which(signCluster$rn == rownames(candDataLG)[i]), 
                        "cutoffLossWP"]
        ## Baseline
        base <-
            signCluster[which(signCluster$rn == rownames(candDataLG)[i]), 
                        "cutoffGainWP"]
        
        ## fallback auf wp cutoff, wenn keine echten maxima identifiziert 
        #werden konnten
        if (is.na(cutG) || length(cutG) == 0) {
            cutG <- cutGWP
        }
        if (is.na(cutL) || length(cutL) == 0) {
            cutL <- cutLWP
        }
        
        ## minimale effectgroesse ueberschritten?
        if (!is.na(cutL) &&
            !is.na(base) && abs(base - cutL) < effectsize[1]) {
            cutL <- NA
        }
        if (!is.na(cutG) &&
            !is.na(base) && abs(base - cutG) < effectsize[2]) {
            cutG <- NA
        }
        #print(paste(i, "G:", cutG, " L:", cutL))
        
        vec <- data.frame(gain = cutG, loss = cutL)
        gainsLosses <- rbind.fill(gainsLosses, vec)
        
        candDataLG[i, ] <-
            ifelse(!is.na(cutG) & candDataLG[i, ] > cutG, 1, candDataLG[i, ])
        candDataLG[i, ] <-
            ifelse(!is.na(cutL) & candDataLG[i, ] < cutL,-1, candDataLG[i, ])
        candDataLG[i, ] <-
            ifelse(candDataLG[i, ] == 1 |
                        candDataLG[i, ] == -1, candDataLG[i, ], 0)
    }
    cat("\n")
    
    return(candDataLG)
}

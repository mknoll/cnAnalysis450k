require(plyr)

#' @title Binarize Data
#'
#' @description
#' Binarizes Data (-1, 0, +1) according to given cutoff list
#'
#' @param data inputdata
#' @param cutoffs cutoffs, calculted e.g. with \code{findCutoffs()}
#' @param effectsize minimal required effectsize as vector (loss, gain)
#'
#' @export
#'
#' @return binarized values
#'
#' @examples
#' print("Please refer to the 'completeWorkflow' vignette!")
segmentData <- function(data, cutoffs, effectsize = c(0, 0)) {
    print("Binarize data ... ")
    completeDATA <- data
    signCluster <- cutoffs
    
    candDataLG <-
        completeDATA[which(rownames(completeDATA) %in% signCluster$rn), ]
    gainsLosses <- NULL
    for (i in 1:length(candDataLG[, 1])) {
        cat(".")
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
    
    return(candDataLG)
}

#' @title
#' Define gain / loss with given cutoffs
#'
#' @description
#' Use absolute cutoffs to create CNV gain / loss matrix.
#'
#' @param data segmentMatrix, e.g. acquired with createSegmentMatrix
#' @param upper Upper cutoff (+1)
#' @param lower Lower cutoff (-1)
#' @param ylim y-axis limits of cutoff plot
#' @param plot Should the data be plotted? Might take a long time.
#'
#' @return matrix, containing -1,0,+1 for losses / gains
#'
#' @export
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
#' segmentDataAbs(candMatrix, upper=0.1, lower=-0.15)
segmentDataAbs <-
    function(data,
            upper,
            lower,
            ylim = NULL,
            plot = TRUE) {
        # Upper < lower?
        if (upper < lower) {
            tmp <- upper
            upper <- lower
            lower <- tmp
        }
        
        col <- c()
        for (i in 1:length(data[1, ])) {
            co <- ifelse(i %% 2 == 0, 3, 4)
            col <- c(col, rep(co, length(data[, 1])))
        }
        if (plot) {
            plot(unlist(c(data)),
                col = col,
                ylab = "values",
                ylim = ylim)
            abline(h = upper, col = 2)
            abline(h = lower, col = 2)
        }
        
        tmpdat <- data
        tmpdat[] <- NA
        tmpdat[data > upper] <- 1
        tmpdat[data < lower] <- -1
        tmpdat[data <= upper & data >= lower] <- 0
        
        return(tmpdat)
    }

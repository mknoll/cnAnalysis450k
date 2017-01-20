#' @title
#' Apply Fishers exact test
#'
#' @description
#' Calculates p values (Fishers exact test) for gain/loss 
#' and given group assignment
#'
#' @param data previously created CN matrix
#' @param cluster group assignent, same order as data columns
#' @param bin Aggregate Fishers p.values to < 0.001, <0.01, 
#' <0.05, <0.1, >0.1
#'
#' @return vector with all p-values
#'
#' @export
#'
#' @examples
#' data <- data.frame(smp1=c(0,1,1,-1,0),
#' smp2=c(0,1,0,1,1),
#' smp3=c(1,1,0,0,0))
#' rownames(data) <- c("chr1:10000", "chr1:50000", 
#' "chr1:100000", "chr1:150000", "chr1:200000")
#' group <- c(1,2,2)
#' calcFisher(data,group,FALSE)
#' calcFisher(data,group,TRUE)
calcFisher <- function(data, cluster, bin = TRUE) {
    fVal <-
        apply(data, 1, function(x) {
            if (length(levels(factor(unlist(x)))) >= 2) {
                stats::fisher.test(table(x, cluster))$p.value
            } else {
                NA
            }
        })
    if (bin) {
        fVal <-
            ifelse(fVal < 0.001,
                    "<0.001",
                    ifelse(
                        fVal < 0.01,
                        "<0.01",
                        ifelse(fVal < 0.05, "<0.05", 
                                ifelse(fVal < 0.1, "<0.1", ">0.1"))
                    ))
    }
    return (fVal)
}

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
#' @return vector with all p values
#'
#' @export
#'
#' @examples
#' print("Please refer to the 'completeWorkflow' vignette!")
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

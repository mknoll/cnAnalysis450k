#' @title
#' Create matrix from conumee data
#'
#' @description
#' Create matrix from conumee segment / bin data.
#'
#' @param data output of runConumee(), might contain segments,
#' bins or transcripts
#'
#' @return n x m matrix, n beeing the amoutn of segments / bins;
#' m equals the amount of samples
#'
#' @import plyr
#'
#' @export
#' 
#' @examples 
#' print("Please refer to the 'completeWorkflow' vignette!")
createConumeeMatrix <- function(data) {
    #sanity check
    if (data$what == "segments") {
        return(createConumeeMatrix.segment(data$data))
    } else if (data$what == "bins") {
        return(createConumeeMatrix.bins(data$data))
    } else if (data$what == "transcripts") {
        return(createConumeeMatrix.tx(data$data))
    } else {
        stop("Wrong function call!")
    }
    ## TODO: Allow for custom call - e.g. isoforms
}

createConumeeMatrix.segment <- function(data, mod = "seg.median") {
    ## all sample names (IDs)
    ids <- levels(factor(data$ID))
    
    ## matrix of segments
    completeDATA <- NULL
    
    ## positions of Segments
    completePos <- NULL
    completeChrom <- NULL
    
    chrs <- paste(rep("chr"), 1:22, sep = "")
    for (ch in chrs) {
        cat(".")
        chrData <- data[data$chrom == ch, ]
        allPosCh <- unique(c(chrData$loc.start, chrData$loc.end))
        
        ## build segments data / chromosome
        tmpData <- matrix(ncol = length(allPosCh), nrow = length(ids), NA)
        rownames(tmpData) <- ids
        colnames(tmpData) <- allPosCh
        allPosCh <- as.numeric(allPosCh)
        
        for (id in ids) {
            idData <- chrData[which(chrData$ID == id), ]
            for (j in 1:length(tmpData[1, ])) {
                if (allPosCh[j] %in% idData$loc.start) {
                    end <- idData[which(idData$loc.start == allPosCh[j]), 
                                    "loc.end"]
                    value <-
                        idData[which(idData$loc.start == allPosCh[j]), mod]
                    k <- 0
                    while (allPosCh[j + k] <= end &&
                            ((j + k) != (length(allPosCh) + 1))) {
                        tmpData[id, j + k] <- value
                        k <- k + 1
                    }
                }
            }
        }
        
        tmp <- t(tmpData)
        tmp <- tmp[!duplicated(rownames(tmp)), ]
        completeDATA <- rbind.fill(completeDATA, data.frame(tmp))
        completeChrom <- c(completeChrom, rep(ch, length(tmp[, 1])))
        completePos <- c(completePos, rownames(tmp))
    }
    
    rownames(completeDATA) <-
        paste(completeChrom, ":", completePos, sep = "")
    
    return(completeDATA)
}

createConumeeMatrix.bins <- function(data) {
    ## matrix of bins
    completeDATA <- data[, -c(1:4)]
    rownames(completeDATA) <-
        paste(data$Chromosome, ":", data$Start, sep = "")
    
    return(completeDATA)
}

createConumeeMatrix.tx <- function(data) {
    ## matrix of transcripts
    tmpdata <- data[!duplicated(data[, 4]) & !is.na(data[, 4]), ]
    completeDATA <- tmpdata[, -c(1:7)]
    rownames(completeDATA) <- tmpdata[, 4]
    
    return(completeDATA)
}

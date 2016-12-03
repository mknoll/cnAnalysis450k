#' @title
#' Create Segment-Matrix
#'
#' @description
#' Create Segment Matrix from Segment-list (created with \code{findSegments()})
#'
#' @param data input data
#' @param p.select 0.05 selection of candidates which have a p.value
#' < p.select for at least 1 sample
#'
#' @return all aligned segments
#'
#' @export
#'
#' @examples
#' print("Please refer to the 'completeWorkflow' vignette!")
createSegmentMatrix <- function(data, p.select = 0.05) {
    print("### Create Segment Matrix ...")
    
    anno <- minfi::getAnnotation(CopyNumber450kData::RGcontrolSetEx)
    annoSorted <- anno[order(anno$chr, anno$pos),]
    
    ##add position
    patDataPos <- data
    patDataPos$POS_START <-
        paste(annoSorted[as.character(patDataPos$startCG), "chr"], ":",
            annoSorted[as.character(patDataPos$startCG), "pos"], sep =
                    "")
    patDataPos$POS_END <-
        paste(annoSorted[as.character(patDataPos$endCG), "chr"], ":",
            annoSorted[as.character(patDataPos$endCG), "pos"], sep =
                    "")
    
    ### create segment matrix
    allPos <-
        unique(c(patDataPos$POS_START[], patDataPos$POS_END[]))
    allSegments <-
        matrix(ncol = length(levels(factor(patDataPos$smp))), 
                nrow = length(allPos), NA)
    colnames(allSegments) <- levels(factor(patDataPos$smp))
    posDF <-
        data.frame(do.call(rbind, strsplit(allPos, ":")), 
                    substr(do.call(rbind, strsplit(allPos, ":"))[, 1], 
                        4, nchar(as.character(
            do.call(rbind, strsplit(allPos, ":"))[, 1]
        ))))
    rownames(posDF) <- allPos
    colnames(posDF) <- c("chr", "pos", "Chromosome")
    posDF <-
        posDF[order(as.numeric(as.character(posDF$Chromosome)), 
                    as.numeric(as.character(posDF$pos))),]
    rownames(allSegments) <- rownames(posDF)
    for (smp in levels(factor(patDataPos$smp))) {
        print(paste("Processing", smp, "..."))
        smpPos <- patDataPos[which(patDataPos$smp == smp),]
        smpPos$ST <-
            do.call(rbind, strsplit(smpPos$POS_START, ":"))[, 2]
        smpPos$EN <-
            do.call(rbind, strsplit(smpPos$POS_END, ":"))[, 2]
        vec <- rep(NA, length(allSegments[, 1]))
        for (i in 1:length(smpPos[, 1])) {
            start <- which(rownames(posDF) == smpPos[i, "POS_START"])
            end <- which(rownames(posDF) == smpPos[i, "POS_END"])
            vec[start:end] <- smpPos[i, "median"]
        }
        allSegments[, smp] <- vec
    }
    
    ## remove all rows which do not contain at least one significant sample
    p.val <- data.frame(
        p.val = data$p.val,
        start = as.numeric(as.character(do.call(
            rbind, strsplit(patDataPos$POS_START, ":")
        )[, 2])),
        end = as.numeric(as.character(do.call(
            rbind, strsplit(patDataPos$POS_END, ":")
        )[, 2])),
        chr = do.call(rbind, strsplit(patDataPos$POS_START, ":"))[, 1]
    )
    
    pos <-
        data.frame(do.call(rbind, strsplit(rownames(allSegments), ":")))
    pos[, 2] <- as.numeric(as.character(pos[, 2]))
    rownames(pos) <- rownames(allSegments)
    pos$SELECT <- NA
    for (i in 1:(length(pos[, 1]) - 1)) {
        sub <-
            p.val[which(p.val$chr == pos$X1[i] &
                            p.val$start >= pos$X2[i] &
                            p.val$end <= pos$X2[i + 1]),]
        pos$SELECT[i] <- any(sub$p.val <= p.select)
    }
    
    return(allSegments[which(pos$SELECT == TRUE),])
}

#' @title Create Matrix from calculated Bins
#' 
#' @export
createBinMatrix <- function(data, val) {
    ##create empty matrix
    m <- matrix(ncol=length(levels(factor(data$smp))), 
                nrow=length(levels(factor(data$startCG))))
    colnames(m) <- levels(factor(data$smp))
    rownames(m) <- levels(factor(data$startCG))
    
    for (p in colnames(m)) {
        subCh <- data[which(data$smp == p),val]
        names(subCh) <- data$startCG[which(data$smp == p)]
        subCh <- subCh[fastmatch::fmatch(rownames(m), names(subCh))]
        m[,p] <- subCh
    }
    
    ##add position/ chrom data
    anno <- getAnnoData(determineArrayType(data))
    anno$chr <- factor(anno$chr, levels=c(paste("chr", 1:22,sep=""), "chrX", "chrY"))
    anno <- anno[which(rownames(anno) %in% rownames(m)),]
    anno <- anno[order(anno$chr, anno$pos),]
    m <- m[fastmatch::fmatch(rownames(anno), rownames(m)),]
    rownames(m) <- paste(rownames(m), anno$chr, anno$pos, sep=":")
    
    return (m)
}
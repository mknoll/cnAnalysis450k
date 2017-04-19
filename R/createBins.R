#' @title
#' Bin CN data
#'
#' @description
#' Calculate mean, median, sd for each bin for a given binsize.
#'
#' @param data sample data
#' @param ctrl control data
#' @param ctrlAll CN data of all control samples
#' @param statistic statistics test to compare groups (controls 
#' vs smp); "t.test" or "wilcoxon" (Mann-Whitney-U Test)
#' @param binsize binsize
#' @param output "ratio" (Sample/Ctrl) or "diff" (Sample-Ctrl)
#' @param arrayType "auto","450k", "EPIC"; auto -> tries to automatically 
#' determine the array type (450k, EPIC)
#' @param noCores number of cores for parallelization: needed 
#' to pass the cran examples :/ 
#' 
#' @return bins with their corresponding value
#'
#' @import foreach
#' @import doParallel
#' @import parallel
#'
#' @export
#'
#' @examples
#' data <- minfi::getCN(minfi::preprocessRaw(minfiData::RGsetEx))
#' ctrlAll <- data[,4,drop=FALSE]
#' ctrl <- data[,5] # apply(ctrlAll, 1, "median")
#' data <- data[,1,drop=FALSE]
#' createBins(data,ctrl,ctrlAll, binsize=50000000, noCores=2)[1:3,]
createBins <-
    function(data,
            ctrl,
            ctrlAll,
            statistic = "wilcoxon",
            binsize = 200000,
            output = "diff",
            arrayType="auto",
            noCores=-1) {
        if (noCores == -1) {
            no_cores <- parallel::detectCores() - 1
            no_cores <- ifelse(no_cores == 0, 1, no_cores)
        } else {
            no_cores <- noCores
        }
        doParallel::registerDoParallel(no_cores)

        print(paste("### Bin CN Data: binsize=", binsize, " ..."))
        
        ##get annotation
        if (arrayType=="auto") {
            anno <- getAnnoData(determineArrayType(data))
        } else {
            anno <- getAnnoData(arrayType)
        }
        anno$chr <- factor(anno$chr, levels=c(paste("chr", 1:22, sep=""), "chrX", "chrY"))
        annoSorted <- anno[order(anno$chr, anno$pos), ]
        
        # Get cg Borders of Chromosomes
        chrs <- paste("chr", 1:22, sep = "")
        chBorder <- NULL
        ch <- c()
        out <- foreach (ch=chrs) %dopar% {
            subCh <- data.frame(anno[which(anno$chr == ch), ])
            subCh <- subCh[order(subCh$pos), ]
            vec <-
                data.frame(
                    chr = ch,
                    start = rownames(subCh)[1],
                    end = rownames(subCh)[length(subCh[, 1])],
                    startPos = annoSorted[which(rownames(annoSorted) 
                                                == rownames(subCh)[1]), "pos"],
                    endPos = annoSorted[which(rownames(annoSorted) 
                                              == rownames(subCh)[length(subCh[, 1])]),"pos"]
                )
            vec
        }
        chBorder <- do.call(rbind, out)
        rownames(chBorder) <- chBorder[, 1]
        
        ##Calculate Bin-Values
        patData <- NULL
        for (j in 1:length(data[1, ])) {
            #different patients
            ct <- ctrl
            ct <- ct[fastmatch::fmatch(rownames(annoSorted), names(ct))]
            ctAll <- ctrlAll
            ctAll <- ctAll[fastmatch::fmatch(rownames(annoSorted), 
                rownames(ctAll)), ,drop=FALSE]
            da <- data[, j]
            names(da) <- rownames(data)
            da <- da[fastmatch::fmatch(rownames(annoSorted), names(da))]
            smpName <- colnames(data)[j]
            cat("\n")
            print(paste("Processing", smpName, "..."))
            
            ##splitpos
            sampleBins <- 
                foreach::foreach(ch=levels(factor(chBorder$chr))) %dopar% {
                cat("\n#", ch, " ")
                subCh <-
                    data.frame(anno[which(anno$chr == ch), c("pos", "chr")])
                subCh <- subCh[order(subCh$pos), ]
                
                fromPos <- chBorder[ch, "startPos"]
                toPos <- chBorder[ch, "endPos"]
                
                ## find bin-borders
                starts <- fromPos
                ends <- c()
                z1 <- ceiling(fromPos / binsize)
                z2 <- floor(toPos / binsize)
                if (z2 > z1) { 
                    for (z in z1:z2) {
                        starts <- c(starts, z * binsize)
                        ends <- c(ends, z * binsize - 1)
                    }
                }
                ends <- c(ends, toPos)
                
                ## caluclate bin values
                startCgs <- c()
                endCgs <- c()
                median <- c()
                mean <- c()
                sd <- c()
                p.val <- c()
                ##use for fastmatch
                namesDa <- names(da)
                namesCt <- names(ct)
                namesCtAll <- rownames(ctAll)
                for (p in 1:length(starts)) {
                    cat(".")
                    daTmp <-
                        da[fastmatch::fmatch(
                            rownames(subCh)[subCh$pos >= starts[p] &
                                subCh$pos <= ends[p]], namesDa)]
                    ctTmp <-
                        ct[fastmatch::fmatch(
                            rownames(subCh)[subCh$pos >= starts[p] &
                                subCh$pos <= ends[p]], namesCt)]
                    ctAllTmp <-
                        ctAll[fastmatch::fmatch(
                            rownames(subCh)[subCh$pos >= starts[p] &
                                subCh$pos <= ends[p]], namesCtAll), ]
                    
                    if (length(daTmp) == 0 || length(ctTmp) == 0) {
                        next
                    }
                    
                    ## get cg-IDs
                    rnTmp <-
                        rownames(subCh)[subCh$pos >= starts[p] & 
                                            subCh$pos <= ends[p]]
                    startCgs <- c(startCgs, rnTmp[1])
                    endCgs <- c(endCgs, rnTmp[length(rnTmp)])
                    
                    ##ratio or difference?
                    if (output == "ratio") {
                        rat <- daTmp / ctTmp
                    } else if (output == "diff") {
                        rat <- daTmp - ctTmp
                    }
                    
                    rat <- rat[!is.na(rat) & !is.infinite(rat)]
                    median <- c(median, median(rat))
                    mean <- c(mean, mean(rat))
                    sd <- c(sd, sd(rat))
                    ##statistics
                    if (statistic == "t.test") {
                        #standard t-test
                        p.val <-
                            c(p.val, t.test(daTmp,
                                ctAllTmp)$p.value)
                    } else if (statistic == "wilcoxon") {
                        #mann-whitney-wilcoxon test
                        p.val <-
                            c(p.val,
                                wilcox.test(daTmp, 
                                ctAllTmp)$p.value)
                    } else {
                        stop("Unknown statistic!")
                    }
                }
                
                chrDF <- data.frame(
                    chr = ch,
                    startCG = startCgs,
                    #start=starts,
                    endCG = endCgs,
                    #end=ends,
                    median = median,
                    mean = mean,
                    sd = sd,
                    smp = smpName,
                    p.val = p.val,
                    statistic = statistic
                )
            }
            patData <- rbind(patData, do.call(rbind, sampleBins))
        }

        stopImplicitCluster()

        return(patData)
    }


#' @title Faster binning of CN data
#' 
#' @export
createBinsFast <- function(data,
                           ctrl,
                           ctrlAll,
                           binsize=50000) {
    ## use parallelization
    no_cores <- parallel::detectCores() - 1
    no_cores <- ifelse(no_cores == 0, 1, no_cores)
    doParallel::registerDoParallel(no_cores)
    
    ##get annotation
    anno <- cnAnalysis450k::getAnnoData(
        cnAnalysis450k::determineArrayType(data))
    anno$chr <- factor(anno$chr, levels=c(paste("chr", 1:22, sep=""), "chrX", "chrY"))
    annoSorted <- anno[order(anno$chr, anno$pos), ]
    
    ##Order controls and data
    ct <- ctrl
    ct <- ct[fastmatch::fmatch(rownames(annoSorted), names(ct))]
    ctAll <- ctrlAll
    ctAll <- ctAll[fastmatch::fmatch(rownames(annoSorted), 
                                     rownames(ctAll)), ,drop=FALSE]
    data <- data[fastmatch::fmatch(rownames(annoSorted), rownames(data)),]
    
    # Get cg Borders of Chromosomes
    chrs <- paste("chr", 1:22, sep = "")
    chBorder <- NULL
    ch <- c()
    out <- foreach (ch=chrs) %dopar% {
        subCh <- data.frame(anno[which(anno$chr == ch), ])
        subCh <- subCh[order(subCh$pos), ]
        vec <-
            data.frame(
                chr = ch,
                start = rownames(subCh)[1],
                end = rownames(subCh)[length(subCh[, 1])],
                startPos = annoSorted[which(rownames(annoSorted) 
                               == rownames(subCh)[1]), "pos"],
                endPos = annoSorted[which(rownames(annoSorted) 
                             == rownames(subCh)[length(subCh[, 1])]),"pos"]
            )
        vec
    }
    chBorder <- do.call(rbind, out)
    rownames(chBorder) <- chBorder[, 1]
    
    
    ##Calculate bins per chromosome
    patData <- foreach::foreach(ch=levels(factor(chBorder$chr))) %dopar% {
        cat("\r#", ch, "             ")
        subCh <- data.frame(anno[which(anno$chr == ch), c("pos", "chr")])
        subCh <- subCh[order(subCh$pos), ]
        
        fromPos <- chBorder[ch, "startPos"]
        toPos <- chBorder[ch, "endPos"]
        
        ## find bin-borders
        starts <- fromPos
        ends <- c()
        z1 <- ceiling(fromPos / binsize)
        z2 <- floor(toPos / binsize)
        if (z2 > z1) { 
            for (z in z1:z2) {
                starts <- c(starts, z * binsize)
                ends <- c(ends, z * binsize - 1)
            }
        }
        ends <- c(ends, toPos)
        
        ## caluclate bin values
        startCgs <- c()
        endCgs <- c()
        
        ### medians
        medians <- matrix(ncol=length(colnames(data)),
                         nrow=length(starts), NA)
        colnames(medians) <- colnames(data)
        
        ##p-values / wilcoxon
        pvals <- matrix(ncol=length(colnames(data)),
                        nrow=length(starts), NA)
        colnames(pvals) <- colnames(data)
        
        namesDa <- rownames(data)
        rmRows <- c()
        
        for (p in 1:length(starts)) {
            mtch <- fastmatch::fmatch(
                rownames(subCh)[subCh$pos >= starts[p] &
                                    subCh$pos <= ends[p]], namesDa)
            
            daTmp <- data[mtch,]
            ctTmp <- ct[mtch]
            ctAllTmp <- ctAll[mtch,]
            
            if (length(daTmp) == 0 || length(ctTmp) == 0) {
                rmRows <- c(rmRows, p)
                next
            }
            
            ## get cg-IDs
            startCgs <- c(startCgs, rownames(subCh)[which(subCh$pos >= starts[p])][1])
            
            ##difference?
            rat <- daTmp - ctTmp
            rat[is.infinite(rat)] <- NA 
            if (length(rat) > length(colnames(data))) {
                medians[p, ] <- apply(rat, 2, "median")
            } else {
                medians[p, ] <- rat
            }
            
            ##statistics
            if (length(rat) > length(colnames(data))) {
                pvals[p, ] <- apply(daTmp, 2, function(x) wilcox.test(x, ctAllTmp)$p.value)
            } else {
                pvals[p, ] <- NA
            }
        }
        medians <- medians[-rmRows,]
        pvals <- pvals[-rmRows,]
        rownames(medians) <- startCgs
        rownames(pvals) <- starts[-rmRows]
        
        chrDF <- list(
            chr=ch, 
            startCGs=startCgs,
            median=medians,
            p.val=pvals,
            statistic="wilcoxon",
            type="diff"
        )
        chrDF
    }

    stopImplicitCluster()
    
    return(patData)
}

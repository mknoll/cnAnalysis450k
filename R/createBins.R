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
#'
#' @return bins with their corresponding value
#'
#' @export
#'
#' @examples
#' print("Please refer to the 'completeWorkflow' vignette!")
createBins <-
    function(data,
            ctrl,
            ctrlAll,
            statistic = "wilcoxon",
            binsize = 200000,
            output = "diff") {
        print(paste("### Bin CN Data: binsize=", binsize, " ..."))
        
        # Get annotation
        anno <- minfi::getAnnotation(CopyNumber450kData::RGcontrolSetEx)
        annoSorted <- anno[order(anno$chr, anno$pos), ]
        
        # Get cg Borders of Chromosomes
        chrs <- paste("chr", 1:22, sep = "")
        chBorder <- NULL
        for (ch in chrs) {
            subCh <- data.frame(anno[which(anno$chr == ch), ])
            subCh <- subCh[order(subCh$pos), ]
            vec <- data.frame(
                chr = ch,
                start = rownames(subCh)[1],
                end = rownames(subCh)[length(subCh[, 1])],
                startPos = subCh$pos[1],
                endPos = subCh$pos[length(subCh[, 1])]
            )
            chBorder <- rbind.fill(chBorder, vec)
        }
        rownames(chBorder) <- chBorder[, 1]
        
        ##Calculate Bin-Values
        patData <- NULL
        for (j in 1:length(data[1, ])) {
            #different patients
            ct <- ctrl
            ct <- ct[fastmatch::fmatch(rownames(annoSorted), names(ct))]
            ctAll <- ctrlAll
            ctAll <- ctAll[match(rownames(annoSorted), rownames(ctAll)), ]
            da <- data[, j]
            da <- da[fastmatch::fmatch(rownames(annoSorted), names(da))]
            smpName <- colnames(data)[j]
            cat("\n")
            print(paste("Processing", smpName, "..."))
            
            sampleBins <- NULL
            ##splitpos
            for (ch in levels(factor(chBorder$chr))) {
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
                for (z in z1:z2) {
                    starts <- c(starts, z * binsize)
                    ends <- c(ends, z * binsize - 1)
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
                            c(p.val, t.test(daTmp[starts[p]:ends[p]],
                                ctAllTmp[starts[p]:ends[p]])$p.value)
                    } else if (statistic == "wilcoxon") {
                        #mann-whitney-wilcoxon test
                        p.val <-
                            c(p.val,
                                wilcox.test(daTmp[starts[p]:ends[p]], 
                                ctAllTmp[starts[p]:ends[p]])$p.value)
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
                sampleBins <- rbind.fill(sampleBins, chrDF)
            }
            patData <- rbind.fill(patData, sampleBins)
        }
        
        return(patData)
    }
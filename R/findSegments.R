#' @title Find Segments in the provided CN data
#'
#' @description
#' Uses data from minfis getCN() function and normalizes 
#' probe-wise against control CN data.
#' Segments are identified with changepoints cpr.var() 
#' function (BinSeg)
#'
#' @param data CN data to evaluate
#' @param ctrl CN data of controls, levels to test to (1 mean / median
#'  over all ctrl samples)
#' @param ctrlAll CN data of all control samples
#' @param statistic statistics test to compare groups (controls vs smp);
#' "t.test" or "wilcoxon" (Mann-Whitney-U Test)
#' @param plot plot changepoints
#' @param ylim if plot=T give ylim
#' @param delta area around changepoints to plot for plot=T
#' @param output "ratio" (Sample/Ctrl) or "diff" (Sample-Ctrl)
#' @param arrayType "auto","450k", "EPIC"; auto -> tries to automatically 
#' determine the array type (450k, EPIC)
#'
#' @return data containing chr, startCG, endCG, segmentmedian, -mean, 
#' SD and samplename
#'
#' @import plyr
#'
#' @export
#'
#' @examples
#' norm <- minfi::getCN(minfi::preprocessRaw(minfiData::RGsetEx))
#' ctrlAll <- norm[,5,drop=FALSE]
#' ctrl <- norm[,4] #ctrl <- apply(ctrlAll, 1, "median")
#' samples <- norm[,1,drop=FALSE]
#' findSegments(samples,ctrl, ctrlAll)[1:4,]
findSegments <-
    function(data,
            ctrl,
            ctrlAll,
            statistic = "wilcoxon",
            plot = FALSE,
            delta = 500,
            output = "diff",
            ylim = NULL,
            arrayType="auto") {
        print("### Find Segments in CN Data ...")
        
	##check if fast version can be used
	if (!plot && output == "diff" && statistic=="wilcoxon") {
	    return (findSegmentsFast(data, ctrl, ctrlAll))
	}

        ##get annotation
        if (arrayType=="auto") {
            anno <- getAnnoData(determineArrayType(data))
        } else {
            anno <- getAnnoData(arrayType)
        }
        annoSorted <- anno[order(anno$chr, anno$pos), ]
        
        # Get cg Borders of Chromosomes
        chrs <- paste("chr", 1:22, sep = "")
        chBorder <- NULL
        for (ch in chrs) {
            subCh <- data.frame(anno[which(anno$chr == ch), ])
            subCh <- subCh[order(subCh$pos), ]
            vec <-
                data.frame(
                    chr = ch,
                    start = rownames(subCh)[1],
                    end = rownames(subCh)[length(subCh[, 1])]
                )
            chBorder <- rbind.fill(chBorder, vec)
        }
        rownames(chBorder) <- chBorder[, 1]
        
        ##Find segments
        patData <- NULL
        
        ct <- ctrl
        ct <- ct[match(rownames(annoSorted), names(ct))]
        ctAll <- ctrlAll
        ctAll <- ctAll[match(rownames(annoSorted), rownames(ctAll)), ]
        
        for (j in 1:length(data[1, ])) {
            #different patients
            da <- data[, j]
            da <- da[match(rownames(annoSorted), names(da))]
            smpName <- colnames(data)[j]
            print(paste("Processing", smpName, "..."))
            
            sampleSegments <- NULL
            ##splitpos
            for (ch in levels(factor(chBorder$chr))) {
                from <- which(names(da) == chBorder[ch, "start"])
                to <- which(names(da) == chBorder[ch, "end"])
                
                ##ratio or difference?
                if (output == "ratio") {
                    rat <- da[from:to] / ct[from:to]
                } else if (output == "diff") {
                    rat <- da[from:to] - ct[from:to]
                }
                
                namesRat <- names(da)[from:to]
                sel <- !is.na(rat) & !is.infinite(rat)
                rat <- rat[sel]
                namesRat <- namesRat[sel]
                
                ##find changepoints
                res <- changepoint::cpt.var(rat, method = "BinSeg")
                
                ##calculate segment data
                if (length(res@cpts) > 1) {
                    ends <- res@cpts
                    starts <- c(0, ends[1:(length(ends) - 1)]) + 1
                    median <- c()
                    mean <- c()
                    sd <- c()
		    p.val <- c()
                    for (pos in 1:length(starts)) {
                        median <- c(median, 
                                    median(rat[starts[pos]:ends[pos]]))
                        mean <- c(mean, 
                                mean(rat[starts[pos]:ends[pos]]))
                        sd <- c(sd, sd(rat[starts[pos]:ends[pos]]))
                        ##statistics
                        if (statistic == "t.test") {
                            #standard t-test
                            p.val <- c(p.val,
                                t.test(da[starts[pos]:ends[pos]], 
                                ctAll[starts[pos]:ends[pos]])$p.value)
                        } else if (statistic == "wilcoxon") {
                            #mann-whitney-wilcoxon test
                            p.val <- c(p.val,
                                wilcox.test(da[starts[pos]:ends[pos]], 
                                ctAll[starts[pos]:ends[pos]])$p.value)
                        } else {
                            stop("Unknown statistic!")
                        }
                    }
                    segDF <- NULL
                } else {
                    starts <- 1
                    ends <- res@cpts
                    median <- median(rat[starts:ends])
                    mean <- mean(rat[starts:ends])
                    sd <- sd(rat[starts:ends])
                    ##statistics
                    if (statistic == "t.test") {
                        #standard t-test
                        p.val <-
                            t.test(da[starts:ends], 
                                    ctAll[starts:ends])$p.value
                    } else if (statistic == "wilcoxon") {
                        #mann-whitney-wilcoxon test
                        p.val <-
                            wilcox.test(da[starts:ends], 
                                    ctAll[starts:ends])$p.value
                    } else {
                        stop("Unknown statistic!")
                    }
                }
                
                ## plot segment borders
                if (plot) {
                    changepoints <- res@cpts
                    par(mfrow = c(1, 1))
                    if (!is.null(ylim)) {
                        plot(rat, main = ch, ylim = ylim)
                    } else {
                        plot(rat, main = ch)#, ylim=c(0.9,1.1))
                    }
                    abline(v = changepoints,
                            col = 2,
                            lwd = 2)
                    par(mfrow = c(3, 3))
                    for (q in 1:length(changepoints)) {
                        cg <- changepoints[q]
                        ##get borders
                        cgMin <- ifelse(cg - delta < 0,
                            0,
                            ifelse(
                                q == 1,
                                cg - delta,
                                ifelse(
                                    cg - delta < changepoints[q - 1] - delta,
                                    changepoints[q - 1],
                                    cg - delta
                                )
                            ))
                        segMedian1 <- median(rat[(cgMin):cg])
                        sd1 <- sd(rat[(cgMin):cg])
                        
                        cgMax <- ifelse(
                            cg + delta > length(rat),
                            length(rat),
                            ifelse(
                                q == length(changepoints),
                                cg + delta,
                                ifelse(
                                    cg + delta > changepoints[q + 1],
                                    changepoints[q + 1],
                                    cg + delta
                                )
                            )
                        )
                        segMedian2 <- median(rat[cg:cgMax])
                        sd2 <- sd(rat[cg:cgMax])
                        
                        ##plot
                        if (!is.null(ylim)) {
                            plot(
                                rat,
                                xlim = c(cgMin, cgMax),
                                ylim = ylim,
                                main = paste(
                                    ch,
                                    ": ",
                                    round(segMedian1, 3),
                                    ",",
                                    round(sd1, 3),
                                    " / ",
                                    round(segMedian2, 3),
                                    ",",
                                    round(sd2, 3)
                                )
                            )
                        } else {
                            plot(
                                rat,
                                xlim = c(cgMin, cgMax),
                                #, ylim=c(0.9, 1.1),
                                main = paste(
                                    ch,
                                    ": ",
                                    round(segMedian1, 3),
                                    ",",
                                    round(sd1, 3),
                                    " / ",
                                    round(segMedian2, 3),
                                    ",",
                                    round(sd2, 3)
                                )
                            )
                        }
                        abline(v = cg,
                                col = 4,
                                lwd = 3)
                        abline(h = segMedian1,
                                col = 3,
                                lwd = 2)
                        abline(h = segMedian2,
                                col = 2,
                                lwd = 2)
                    }
                    abline(v = changepoints, col = 2)
                }
                
                segDF <- data.frame(
                    chr = ch,
                    startCG = namesRat[starts],
                    #start=starts,
                    endCG = namesRat[ends],
                    #end=ends,
                    median = median,
                    mean = mean,
                    sd = sd,
                    smp = smpName,
                    p.val = p.val
                )
                
                sampleSegments <- rbind.fill(sampleSegments, segDF)
            }
            
            #plot(sampleSegments$mean, col=sampleSegments$chr)
            patData <- rbind.fill(patData, sampleSegments)
        }
        return (patData)
    }

findSegmentsFast <-
    function(data,
            ctrl,
            ctrlAll,
            arrayType="auto") {
        print("### Find Segments in CN Data ...")
        
        ##get annotation
        ## ggf cachen
        if (arrayType=="auto") {
            anno <- cnAnalysis450k::getAnnoData(
	    cnAnalysis450k::determineArrayType(data))
        } else {
            anno <- cnAnalysis450k::getAnnoData(arrayType)
        }
        annoSorted <- anno[order(anno$chr, anno$pos), ]
        
        ## parallel
        no_cores <- parallel::detectCores() - 1
        no_cores <- ifelse(no_cores == 0, 1, no_cores)
        doParallel::registerDoParallel(no_cores)
        
        ##Find segments
        patData <- NULL
        ct <- ctrl[match(rownames(annoSorted), names(ctrl))]
        ctAll <- ctrlAll[match(rownames(annoSorted), rownames(ctrlAll)), ]
        smp <- data[match(rownames(annoSorted), rownames(data)),]
        
        # Get cg Borders of Chromosomes
        chrs <- paste("chr", 1:22, sep = "")
        chBorder <- NULL
        out <- foreach (ch=chrs) %dopar% {
            subCh <- data.frame(anno[which(anno$chr == ch), ])
            subCh <- subCh[order(subCh$pos), ]
            vec <-
                data.frame(
                    chr = ch,
                    start = rownames(subCh)[1],
                    end = rownames(subCh)[length(subCh[, 1])],
                    startI = which(rownames(smp) 
                        == rownames(subCh)[1]),
                    endI = which(rownames(smp) 
                        == rownames(subCh)[length(subCh[, 1])])
                )
            vec
        }
        chBorder <- do.call(rbind, out)
        rownames(chBorder) <- chBorder[, 1]
        
        for (j in 1:length(data[1, ])) {
            smpName <- colnames(data)[j]
            print(paste("Processing", smpName, "..."))
            
            #different patients
            da <- smp[,j]
            
            sampleSegments <- NULL
            ##splitpos
            out <- foreach(ch=unique(chBorder$chr)) %dopar% {
                from <- chBorder[ch,"startI"]
                to <- chBorder[ch,"endI"]
                
                rat <- da[from:to] - ct[from:to]
                namesRat <- names(da)[from:to]
                
                sel <- !is.na(rat) & !is.infinite(rat)
                rat <- rat[sel]
                namesRat <- namesRat[sel]
                
                ##find changepoints
                res <- changepoint::cpt.var(rat, method = "BinSeg")
                
                ##calculate segment data
                if (length(res@cpts) > 1) {
                    ends <- res@cpts
                    starts <- c(0, ends[1:(length(ends) - 1)]) + 1
                    median <- c()
		    p.val <- c()
                    for (pos in 1:length(starts)) {
                        median <- c(median, 
                                    median(rat[starts[pos]:ends[pos]]))
                        ##statistics
                        p.val <- c(p.val, 
                            wilcox.test(
                                da[starts[pos]:ends[pos]], 
                                ctAll[starts[pos]:ends[pos]])$p.value)
                    }
                    segDF <- NULL
                } else {
                    starts <- 1
                    ends <- res@cpts
                    median <- median(rat[starts:ends])
                    p.val <-
                            wilcox.test(da[starts:ends], 
                                        ctAll[starts:ends])$p.value
                }
                segDF <- data.frame(
                    chr = ch,
                    startCG = namesRat[starts],
                    #start=starts,
                    endCG = namesRat[ends],
                    #end=ends,
                    median = median,
                    smp = smpName,
                    p.val = p.val
                )
                segDF
            }
            sampleSegments <- do.call(rbind, out)
            
            patData <- rbind.fill(patData, sampleSegments)
        }
        
        doParallel::stopImplicitCluster()
        
        
        return (patData)
    }

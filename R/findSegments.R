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
#' @param arrayType "auto","450k", "EPIC"; auto -> tries to automatically determine 
#' the array type (450k, EPIC)
#'
#' @return data containing chr, startCG, endCG, segmentmedian, -mean, 
#' SD and samplename
#'
#' @import plyr
#'
#' @export
#'
#' @examples
#' data <- data.frame(
#' smp1=c(-8.12, -5.225, -3.24, -3.61),
#' smp2=c(-5.0, -3.98, -4.06, -4.5),
#' smp3=c(NA, -2.48, -2.27, -2.1)
#' )
#' ctrlAll <- data.frame(
#' ctl1=c(1.0, -3.6, 0.7, -0.73),
#' ctl2=c(-0.4, -4.1, -4.2, -3.9),
#' ctl2=c(0.74, -1.12, -2.8, -1.67)
#' )
#' ctrl <- apply(ctrlAll, 1, "median")
#' findSegments(data,ctrlAll,ctrl)
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
        for (j in 1:length(data[1, ])) {
            #different patients
            ct <- ctrl
            ct <- ct[match(rownames(annoSorted), names(ct))]
            ctAll <- ctrlAll
            ctAll <- ctAll[match(rownames(annoSorted), rownames(ctAll)), ]
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
                    for (pos in 1:length(starts)) {
                        median <- c(median, 
                                    median(rat[starts[pos]:ends[pos]]))
                        mean <- c(mean, 
                                mean(rat[starts[pos]:ends[pos]]))
                        sd <- c(sd, sd(rat[starts[pos]:ends[pos]]))
                        ##statistics
                        if (statistic == "t.test") {
                            #standard t-test
                            p.val <-
                                t.test(da[starts[pos]:ends[pos]], 
                                ctAll[starts[pos]:ends[pos]])$p.value
                        } else if (statistic == "wilcoxon") {
                            #mann-whitney-wilcoxon test
                            p.val <-
                                wilcox.test(da[starts[pos]:ends[pos]], 
                                ctAll[starts[pos]:ends[pos]])$p.value
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

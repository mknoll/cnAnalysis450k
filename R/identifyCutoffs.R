#' @title
#' Finding cutoffs
#'
#' @description
#' Analyzes row-wise histograms by selecting the global maximum (0),
#' as well as one further local maximum right (+) and left (-).
#' Can make use of previous/following rows.
#'
#' @param data Inputdata
#' @param zeichne Debugging
#' @param ignoreNAs F / T: if NAs occur, proceed, otherwise fail
#' @param eps difference between x points (internal calculation)
#' @param proximity how many previous / following rows should be 
#' considered for cutoff determination? c(prev, folw)
#'
#' @import plyr
#' @import graphics
#' @import stats
#'
#' @export
#'
#' @return identified cutoffs
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
#' findCutoffs(candMatrix, proximity=c(1,0))
findCutoffs <-
    function(data,
            zeichne = FALSE,
            ignoreNAs = FALSE,
            eps = 0.01,
            proximity = c(0, 0)) {
        # Check for NAs
        if (any(is.na(data)) && !ignoreNAs) {
            stop(
                "NAs found! Set ignoreNAs to T if 
                you want to proceed anyhow. Otherwise remove NA rows!"
            )
        }
        
        if (proximity[1] == 0 && proximity[2] == 0 && zeichne == FALSE) {
            ## Keine vorgaenger / Nachfogler beruecksichtigen
            return (findCutoffsRunFast(data, ignoreNAs, eps))
        } else {
            chrR <- do.call(rbind, strsplit(rownames(data), ":"))
            outData <- NULL
            for (ch in levels(factor(chrR[, 1]))) {
                print(paste("Find cutoffs:", ch))
                chrDat <- data[which(chrR[, 1] == ch), ]
                outData <-
                    rbind.fill(outData,
                                findCutoffsRun(chrDat, zeichne, ignoreNAs,
                                            eps, proximity))
            }
            return(outData)
        }
    }



#' @title
#' Finding cutoffs
#'
#' @description
#' Analyzes row-wise histograms by selecting the global maximum (0),
#' as well as one further local maximum right (+) and left (-)
#'
#' @param data Inputdata
#' @param zeichne Debugging
#' @param ignoreNAs F / T: if NAs occure, proceed, otherwise fail
#' @param eps difference between x points (internal calculation)
#' @param proximity how many previous / following rows should be 
#' considered for cutoff determination? c(prev, folw)
#'
#' @import plyr
#' @import graphics
#' @import stats
#' @import utils
#'
#' @return identified cutoffs
findCutoffsRun <-
    function(data, zeichne, ignoreNAs, eps, proximity) {
        print("### Find Cutoffs (+/-) ... ")
        subD <- data.matrix(data)
        
        groupsDATA <- NULL
        #find differing groups: row wise
        iBak <- NULL
        total <- length(subD[, 1])
        for (i in 1:total) {
            ## progress
            if (i %% 10 == 0) { 
                utils::flush.console(); 
                cat("\r",round(i/total*100,2),"%   ") 
            }
            
            ##Determine borders for i
            iBak <- i
            i <-
                c((i - ifelse(i < proximity[1], 0, proximity[1])):(i + ifelse(
                    proximity[2] + i > length(subD[, 1]), 0, proximity[2]
                )))
            if (zeichne) {
                hist(subD[i, ,drop=FALSE], breaks = 100)
            }
            
            epsilon <- eps
            if (any(is.na(subD[i, ]))) {
                tmpNA <- subD[i, ]
                warning(paste(
                    rownames(subD)[i],
                    ", debug(i:",
                    i,
                    "); ",
                    length(tmpNA[is.na(tmpNA)]),
                    " NAs found!"
                ))
                #Enough points to proceed?
                if (length(subD[i, ]) - length(tmpNA[is.na(tmpNA)]) <= 2) {
                    warning(
                        paste(
                            "Only 2 datapoints for",
                            rownames(subD)[i],
                            " (debug: i=",
                            i,
                            ") -> skipping!"
                        )
                    )
                    next
                }
            }
            
            minV <- min(subD[i, ], na.rm = TRUE) - epsilon
            maxV <- max(subD[i, ], na.rm = TRUE) + epsilon
            ##fit function
            df <- approxfun(density(subD[i, ], na.rm = TRUE))
            
            #Sollte mit einer numerisch etwas stabieleren 
            #Loesung ersetzt werden.
            xVal <- seq(from = minV,
                        to = maxV,
                        by = eps)
            if (zeichne)  {
                plot(density(subD[i, ,drop=FALSE], na.rm = TRUE))
            }
            
            ##find
            d1 <- diff(df(xVal))
            ## Minima
            posD1Min <- c()
            for (j in 2:length(d1)) {
                pos <- which(d1[j - 1] < 0 & d1[j] > 0)
                if (length(pos) == 1) {
                    posD1Min <- c(posD1Min, j)
                }
            }
            if (zeichne) {
                points(xVal[posD1Min], df(xVal[posD1Min]), col = 2)
            }
            minV <- data.frame(key = xVal[posD1Min], 
                                val = df(xVal[posD1Min]))
            
            ## Maxima
            posD1Max <- c()
            for (j in 2:length(d1)) {
                pos <- which(d1[j - 1] > 0 & d1[j] < 0)
                if (length(pos) == 1) {
                    posD1Max <- c(posD1Max, j)
                }
            }
            if (zeichne) {
                points(xVal[posD1Max], df(xVal[posD1Max]), col = 3)
            }
            #order maxima
            maxV <- data.frame(key = xVal[posD1Max], val = df(xVal[posD1Max]))
            maxV <- maxV[rev(order(maxV$val)), ]
            if (zeichne) {
                abline(v = maxV$key[1], col = 6)
            }
            assumedBaseline <- maxV$key[1]
            
            ## wendepunkte
            wendepunkte <- c()
            d2 <- diff(df(xVal))
            if (length(d2) >= 3) {
                for (k in 2:(length(d2) - 1)) {
                    if (is.na(d2[k - 1]) || is.na(d2[k]) || is.na(d2[k + 1])) {
                        next
                    }
                    if (d2[k - 1] < d2[k] &&
                            d2[k + 1] < d2[k] || d2[k - 1] > d2[k] &&
                            d2[k + 1] > d2[k])  {
                        wendepunkte <- c(wendepunkte, k)
                    }
                }
                if (zeichne) {
                    points(xVal[wendepunkte], df(xVal[wendepunkte]), col = 4)
                }
            }
            
            ###########################
            ##potentielle cutoffs: drei wendepunkte zwischen zwei 
            #extremalstellen od extremalstelle u rand
            wp <-
                data.frame(key = xVal[wendepunkte], 
                            val = df(xVal[wendepunkte]))
            #rechts
            assumedSepGainWP <- NA
            extrR <- rbind(minV[which(minV$key > maxV$key[1]), ])
            if (dim(extrR)[[1]] == 0) {
                ##rechts kein weiteres extremum -> randwert
                wpTmp <- wp[which(wp$key > maxV$key[1]), ]
                wpTmp <- wpTmp[order(wpTmp$key), ]
                if (length(wpTmp[, 1]) >= 3) {
                    assumedSepGainWP <- mean(wpTmp[2:3, "key"])
                    if (length(assumedSepGainWP) > 1) {
                        assumedSepGainWP <- mean(assumedSepGainWP)
                    }
                    if (zeichne) {
                        abline(v = assumedSepGainWP,
                                lwd = 1,
                                col = 4)
                    }
                }
            } else {
                ## rechts weiteres extremum; sind vorher 3 wp vorhanden?
                wpTmp <-
                    wp[which(wp$key > maxV$key[1] & wp$key < extrR$key[1]), ]
                if (length(wpTmp[, 1]) >= 3) {
                    assumedSepGainWP <-  mean(wpTmp[2:3, "key"])
                    if (length(assumedSepGainWP) > 1) {
                        assumedSepGainWP <- mean(assumedSepGainWP)
                    }
                    if (zeichne) {
                        abline(v = assumedSepGainWP,
                                lwd = 1,
                                col = 4)
                    }
                }
            }
            #links
            assumedSepLossWP <- NA
            extrL <- rbind(minV[which(minV$key < maxV$key[1]), ])
            if (dim(extrL)[[1]] == 0) {
                ##links kein weiteres extremum -> randwert
                wpTmp <- wp[which(wp$key < maxV$key[1]), ]
                wpTmp <- wpTmp[order(wpTmp$key), ]
                if (length(wpTmp[, 1]) >= 3) {
                    assumedSepLossWP <-  mean(wpTmp[2:3, "key"])
                    if (length(assumedSepLossWP) > 1) {
                        assumedSepLossWP <- mean(assumedSepLossWP)
                    }
                    if (zeichne) {
                        abline(v = assumedSepLossWP,
                                lwd = 1,
                                col = 4)
                    }
                }
            } else {
                ## links weiteres extremum; sind vorher 3 wp vorhanden?
                wpTmp <-
                    wp[which(wp$key < maxV$key[1] & wp$key > extrL$key[1]), ]
                if (length(wpTmp[, 1]) >= 3) {
                    assumedSepLossWP <- mean(wpTmp[2:3, "key"])
                    if (length(assumedSepLossWP) > 1) {
                        assumedSepLossWP <- mean(assumedSepLossWP)
                    }
                    if (zeichne) {
                        abline(v = assumedSepLossWP,
                                lwd = 1,
                                col = 4)
                    }
                }
            }
            
            
            ###########################
            #### find losses
            maxVLoss <- maxV[which(maxV$key <= maxV$key[1]), ]
            #find separating value for first and second max
            minCand <-
                minV[which(minV$key <= maxVLoss$key[1] &
                                minV$key >= maxVLoss$key[2]), ]
            assumedSepLoss <- NA
            if (dim(minCand)[[1]] > 0) {
                assumedSepLoss <-
                    minCand[which(minCand$val == min(minCand$val)), "key"]
                if (length(assumedSepLoss) > 1) {
                    assumedSepLoss <- mean(assumedSepLoss)
                }
                if (zeichne) {
                    abline(v = maxVLoss$key[2], col = 6)
                    abline(v = assumedSepLoss,
                            lwd = 2,
                            col = 4)
                }
            }
            
            
            #### find gains
            maxVGain <- maxV[which(maxV$key >= maxV$key[1]), ]
            #find separating value for first and second max
            minCand <-
                minV[which(minV$key >= maxVGain$key[1] &
                                minV$key <= maxVGain$key[2]), ]
            assumedSepGain <- NA
            if (dim(minCand)[[1]] > 0) {
                assumedSepGain <-
                    minCand[which(minCand$val == min(minCand$val)), "key"]
                if (length(assumedSepGain) > 1) {
                    assumedSepGain <- mean(assumedSepGain)
                }
                if (zeichne) {
                    abline(v = maxVGain$key[2], col = 6)
                    abline(v = assumedSepGain,
                            lwd = 2,
                            col = 4)
                }
            }
            ######################################
            i <- iBak
            
            ##assemble data
            vec <- data.frame(
                rn = rownames(subD)[i],
                i = i,
                cutoffLoss = assumedSepLoss,
                cutoffGain = assumedSepGain,
                cutoffLossWP = assumedSepLossWP,
                cutoffGainWP = assumedSepGainWP,
                baseline = assumedBaseline
            )
            
            groupsDATA <- rbind.fill(groupsDATA, vec)
        }
        cat("\n")
        
        return(groupsDATA)
    }



#' @title
#' Identifying cutoffs
#' 
#' @description
#' Analyzes row-wise histograms by selecting the global maximum (0),
#' as well as one further local maximum right (+) and left (-)
#'
#' @param data Inputdata
#' @param ignoreNAs F / T: if NAs occure, proceed, otherwise fail
#' @param eps difference between x points (internal calculation)
#'
#' @import plyr
#' @import stats
#' @import foreach
#' @import doParallel
#' @import parallel
#'
#' @return identified cutoffs
findCutoffsRunFast <-
    function(data, ignoreNAs, eps) {
        warning("Might be unstable!")
        no_cores <- parallel::detectCores() - 1
        no_cores <- ifelse(no_cores == 0, 1, no_cores)
        doParallel::registerDoParallel(no_cores)
        
        print("### Find Cutoffs (+/-) ... ")
        ##Enable parallelization
        proximity=c(0,0)
        
        subD <- data.matrix(data)
        
        #find differing groups: row wise
        i<-0
        out <- foreach::foreach(i=1:length(subD[,1])) %dopar% {
            epsilon <- eps
            cont <- TRUE
            if (any(is.na(subD[i, ]))) {
                tmpNA <- subD[i, ]
                warning(paste(
                    rownames(subD)[i],
                    ", debug(i:",
                    i,
                    "); ",
                    length(tmpNA[is.na(tmpNA)]),
                    " NAs detected!"
                ))
                #Enough points to proceed?
                if (length(subD[i, ]) - length(tmpNA[is.na(tmpNA)]) <= 2) {
                    warning(
                        paste(
                            "Only 2 datapoints for",
                            rownames(subD)[i],
                            " (debug: i=",
                            i,
                            ") -> skipping!"
                        )
                    )
                    cont <- FALSE
                }
            }
            
            if (cont) {
                minV <- min(subD[i, ], na.rm = TRUE) - epsilon
                maxV <- max(subD[i, ], na.rm = TRUE) + epsilon
                ##fit function
                df <- approxfun(density(subD[i, ], na.rm = TRUE))
                
                #Sollte mit einer etwas stabileren 
                #Loesung ersetzt werden.
                xVal <- seq.int(from=minV, to=maxV, by=eps)
                
                ##find
                d1 <- diff(df(xVal))
                ## Minima
                posD1Min <- .Internal(which(d1[-length(d1)]<0 & d1[-1]>0))+1
                minV <- data.frame(
                    key = xVal[posD1Min], val = df(xVal[posD1Min]))
                
                ## Maxima
                posD1Max <- .Internal(which(d1[-length(d1)]>0 & d1[-1]<0))+1
                maxV <- data.frame(
                    key = xVal[posD1Max], val = df(xVal[posD1Max]))
                
                maxV <- maxV[rev(order(maxV$val)), ]
                assumedBaseline <- maxV$key[1]
                
                ## wendepunkte
                d2 <- diff(df(xVal))
                wendepunkte <- c()
                l=length(d2)
                if (l >= 3) {
                    l1=l-1
                    wendepunkte <- c(.Internal(which(d2[-c(1:2)] > d2[-c(1,l)] & d2[-c(1:2)] < d2[-c(l1,l)]))+1,
                            .Internal(which(d2[-c(1:2)] < d2[-c(1,l)] & d2[-c(1:2)] > d2[-c(l1,l)]))+1)
                }
                
                
                ###########################
                ##potentielle cutoffs: drei wendepunkte zwischen zwei 
                #extremalstellen od extremalstelle u rand
                wp <-
                    data.frame(key = xVal[wendepunkte], 
                                val = df(xVal[wendepunkte]))
                #rechts
                assumedSepGainWP <- NA
                extrR <- rbind(minV[which(minV$key > maxV$key[1]), ])
                if (dim(extrR)[[1]] == 0) {
                    ##rechts kein weiteres extremum -> randwert
                    wpTmp <- wp[which(wp$key > maxV$key[1]), ]
                    wpTmp <- wpTmp[order(wpTmp$key), ]
                    if (length(wpTmp[, 1]) >= 3) {
                        assumedSepGainWP <- mean(wpTmp[2:3, "key"])
                        if (length(assumedSepGainWP) > 1) {
                            assumedSepGainWP <- mean(assumedSepGainWP)
                        }
                    }
                } else {
                    ## rechts weiteres extremum; sind vorher 3 wp vorhanden?
                    wpTmp <-
                        wp[which(wp$key > maxV$key[1] 
                            & wp$key < extrR$key[1]), ]
                    if (length(wpTmp[, 1]) >= 3) {
                        assumedSepGainWP <-  mean(wpTmp[2:3, "key"])
                        if (length(assumedSepGainWP) > 1) {
                            assumedSepGainWP <- mean(assumedSepGainWP)
                        }
                    }
                }
                #links
                assumedSepLossWP <- NA
                extrL <- rbind(minV[which(minV$key < maxV$key[1]), ])
                if (dim(extrL)[[1]] == 0) {
                    ##links kein weiteres extremum -> randwert
                    wpTmp <- wp[which(wp$key < maxV$key[1]), ]
                    wpTmp <- wpTmp[order(wpTmp$key), ]
                    if (length(wpTmp[, 1]) >= 3) {
                        assumedSepLossWP <-  mean(wpTmp[2:3, "key"])
                        if (length(assumedSepLossWP) > 1) {
                            assumedSepLossWP <- mean(assumedSepLossWP)
                        }
                    }
                } else {
                    ## links weiteres extremum; sind vorher 3 wp vorhanden?
                    wpTmp <-
                        wp[which(wp$key < maxV$key[1] 
                            & wp$key > extrL$key[1]), ]
                    if (length(wpTmp[, 1]) >= 3) {
                        assumedSepLossWP <- mean(wpTmp[2:3, "key"])
                        if (length(assumedSepLossWP) > 1) {
                            assumedSepLossWP <- mean(assumedSepLossWP)
                        }
                    }
                }
                
                
                ###########################
                #### find losses
                maxVLoss <- maxV[which(maxV$key <= maxV$key[1]), ]
                #find separating value for first and second max
                minCand <-
                    minV[which(minV$key <= maxVLoss$key[1] &
                                    minV$key >= maxVLoss$key[2]), ]
                assumedSepLoss <- NA
                if (dim(minCand)[[1]] > 0) {
                    assumedSepLoss <-
                        minCand[which(minCand$val == min(minCand$val)), "key"]
                    if (length(assumedSepLoss) > 1) {
                        assumedSepLoss <- mean(assumedSepLoss)
                    }
                }
                
                
                #### find gains
                maxVGain <- maxV[which(maxV$key >= maxV$key[1]), ]
                #find separating value for first and second max
                minCand <-
                    minV[which(minV$key >= maxVGain$key[1] &
                                    minV$key <= maxVGain$key[2]), ]
                assumedSepGain <- NA
                if (dim(minCand)[[1]] > 0) {
                    assumedSepGain <-
                        minCand[which(minCand$val == min(minCand$val)), "key"]
                    if (length(assumedSepGain) > 1) {
                        assumedSepGain <- mean(assumedSepGain)
                    }
                }
                
                ##assemble data
                vec <- data.frame(
                    rn = rownames(subD)[i],
                    i = i,
                    cutoffLoss = assumedSepLoss,
                    cutoffGain = assumedSepGain,
                    cutoffLossWP = assumedSepLossWP,
                    cutoffGainWP = assumedSepGainWP,
                    baseline = assumedBaseline
                )
                
                vec
            }
        }
        
        ret <- do.call(rbind, out)
        ret <- data.frame(ret)
        ret <- ret[order(ret$i),]
        
        doParallel::stopImplicitCluster()
        
        return(ret)
    }

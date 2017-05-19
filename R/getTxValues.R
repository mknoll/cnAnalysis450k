#' @title
#' Get TX CNV values
#'
#' @description
#' Calculates CN values for given list of transcripts.
#'
#' @param op median or mean
#' @param data input data
#' @param ctrl control data
#' @param ctrlAll CN data of all control samples
#' @param statistic statistics test to compare groups (controls
#' vs smp); "t.test" or "wilcoxon" (Mann-Whitney-U Test)
#' @param tx vector of list of transcripts#'
#' @param output "ratio" (Sample/Ctrl) or "diff" (Sample-Ctrl)
#' @param arrayType "auto","450k", "EPIC"; auto -> tries to automatically 
#' determine the array type (450k, EPIC)
#' 
#' @return list with tx-names, p-values and used statistical test
#' 
#' @export
#' 
#' @examples
#' data <- data.frame(
#' smp1=c(-8.12, -5.225, -3.24, -3.62),
#' smp2=c(-5.0, -3.98, -4.06, -4.5),
#' smp3=c(NA, -2.48, -2.27, -2.1)
#' )
#' ctrlAll <- data.frame(
#' ctl1=c(1.0, -3.6, 0.7, -0.73),
#' ctl2=c(-0.4, -4.1, -4.2, -3.9),
#' ctl2=c(0.74, -1.12, -2.8, -1.67)
#' )
#' rownames(data) <- c(
#' "cg05132306", "cg15527168",
#' "cg17434257", "cg17592667"
#' )
#' rownames(ctrlAll) <- rownames(data)
#' ctrl <- c(0.74, -3.6, -2.8, -1.67)
#' #ctrl <- apply(ctrlAll, 1, "median")
#' names(ctrl) <- rownames(data)
#' getTxValues(data, ctrl, ctrlAll, "uc001aih.1", arrayType="450k")
getTxValues <-
    function(data,
            ctrl,
            ctrlAll,
            tx,
            output = "diff",
            statistic = "wilcoxon",
            op = "median",
            arrayType="auto") {
        print("Get CN values for transcripts ...")
        
        ##get annotation
        if (arrayType=="auto") {
            anno <- getAnnoData(determineArrayType(data))
        } else {
            anno <- getAnnoData(arrayType)
        }
        anno <- anno[which(rownames(anno) %in% rownames(ctrlAll)),]
        
        # Get Tx position
        txsel <-
        AnnotationDbi::select(
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
        tx,
        c("CDSCHROM", "TXSTART", "TXEND"),
        "TXNAME"
        )
        txsel <- txsel[which(!is.na(txsel$CDSCHROM)),]
        txsel <- txsel[which(!duplicated(txsel$TXNAME)),]
        
        txVal <- NULL
        p.val <- c()
        for (i in 1:length(txsel[, 1])) {
            cat(".")
            cgs <-
                rownames(anno[which(
                    anno$chr == txsel$CDSCHROM[i] &
                        anno$pos >= txsel$TXSTART[i] &
                        anno$pos <= txsel$TXEND[i]
                ),])
            if (length(cgs) == 0) {
                next
            }
            ct <- ctrl[cgs]
            ctAll <- ctrlAll[cgs,]
            da <- data[cgs,]
            
            if (length(ct) > 1) {
                if (output == "ratio") {
                    vec <- data.frame(apply(da / ct, 2, op))
                } else if (output == "diff") {
                    vec <- data.frame(apply(da - ct, 2, op))
                }
            } else {
                if (output == "ratio") {
                    vec <- data.frame(da / ct)
                } else if (output == "diff") {
                    vec <- data.frame(da - ct)
                }
            }
            
            ##statistics
            if (statistic == "t.test") {
                #standard t-test
                p.val <- c(p.val, t.test(unlist(da), unlist(ctAll))$p.value)
            } else if (statistic == "wilcoxon") {
                #mann-whitney-wilcoxon test
                p.val <- c(p.val, 
                    wilcox.test(unlist(da), unlist(ctAll))$p.value)
            } else {
                stop("Unknown statistic!")
            }
            
            if (length(txVal) == 0) {
                txVal <- vec
            } else {
                txVal <- cbind(txVal, vec)
            }
            colnames(txVal)[length(txVal[1,])] <- txsel$TXNAME[i]
        }
        
        names(p.val) <- colnames(txVal)
        ret <- list(data = t(txVal),
                    p.val = p.val,
                    statistic = statistic)
        return(ret)
    }


#' @title
#' Get TX CNV values
#'
#' @description
#' Calculates CN values for given list of transcripts, uses standard
#' values to increase speed: median; wilcoxon; difference
#'
#' @param data input data
#' @param ctrl control data
#' @param ctrlAll CN data of all control samples
#' @param tx vector of list of transcripts
#' @param arrayType "auto","450k", "EPIC"; auto -> tries to automatically 
#' determine the array type (450k, EPIC)
#'
#' @return list with tx-names, p-values and used statistical test
#'
#' @import foreach
#' @import doParallel
#' @import parallel
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
#' rownames(data) <- c(
#' "cg05132306", "cg15527168",
#' "cg17434257", "cg17592667"
#' )
#' rownames(ctrlAll) <- rownames(data)
#' ctrl <- apply(ctrlAll, 1, "median")
#' names(ctrl) <- rownames(data)
#' getTxValues(data, ctrl, ctrlAll, "uc001aih.1", arrayType="450k")
getTxValuesFast <- function(data, ctrl, ctrlAll, tx, arrayType="auto") {
    no_cores <- parallel::detectCores() - 1
    no_cores <- ifelse(no_cores == 0, 1, no_cores)
    doParallel::registerDoParallel(no_cores)
    
    ##get annotation
    if (arrayType=="auto") {
        anno <- cnAnalysis450k::getAnnoData(
            cnAnalysis450k::determineArrayType(data))
    } else {
        anno <- cnAnalysis450k::getAnnoData(arrayType)
    }
    
    # Get Tx position
    txsel <-
    AnnotationDbi::select(
    TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
    tx,
    c("CDSCHROM", "TXSTART", "TXEND"),
    "TXNAME"
    )
    txsel <- txsel[which(!is.na(txsel$CDSCHROM)),]
    txsel <- txsel[which(!duplicated(txsel$TXNAME)),]
    
    if (length(anno[, 1]) != length(ctrl) ||
        length(ctrl) != length(ctrlAll[, 1]) ||
        length(data[, 1] != length(ctrl))) {
        warning("Incomplete data! Annotation dimensions
                do not match the data dimensions -> missing values!")
    }
    anno <- anno[which(rownames(anno) %in% rownames(ctrlAll)),]
    
    i <- NULL
    outData <- foreach::foreach(i = 1:length(txsel[, 1])) %dopar% {
        cgs <-
            which(
                anno$chr == txsel$CDSCHROM[i] &
                    anno$pos >= txsel$TXSTART[i] &
                    anno$pos <= txsel$TXEND[i]
            )
        
        if (length(cgs) > 0) {
            ct <- ctrl[cgs]
            ctAll <- ctrlAll[cgs,]
            da <- data[cgs,]
            
            if (length(ct) > 1) {
                vec <- data.frame(apply(da - ct, 2, "median"))
            } else {
                vec <- data.frame(da - ct)
            }
            
            ##statistics
            p.val <- wilcox.test(unlist(da), unlist(ctAll))$p.value
            
            list(data = vec,
                pVal = p.val,
                name = txsel$TXNAME[i])
        }
    }
    doParallel::stopImplicitCluster()
    
    txVal <- do.call(cbind, do.call(cbind, outData)["data",])
    colnames(txVal) <- unlist(do.call(cbind, outData)["name",])
    pValAll <- unlist(do.call(cbind, outData)["pVal",])
    
    ret <- list(data = t(txVal),
                p.val = pValAll,
                statistic = "wilcoxon")
    return(ret)
    }

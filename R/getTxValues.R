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
#'
#' @return list with tx-names, p-values and used statistical test
#'
#' @export
#'
#' @examples
#' print("Please refer to the 'completeWorkflow' vignette!")
getTxValues <-
    function(data,
            ctrl,
            ctrlAll,
            tx,
            output = "diff",
            statistic = "wilcoxon",
            op = "median") {
        print("Get CN values for transcripts ...")
        
        # Get annotation
        anno <-
            minfi::getAnnotation(CopyNumber450kData::RGcontrolSetEx)
        
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
                p.val <- c(p.val, t.test(da, ctAll)$p.value)
            } else if (statistic == "wilcoxon") {
                #mann-whitney-wilcoxon test
                p.val <- c(p.val, wilcox.test(da, ctAll)$p.value)
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
#' @param tx vector of list of transcripts#'
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
#' print("Please refer to the 'completeWorkflow' vignette!")
getTxValuesFast <- function(data, ctrl, ctrlAll, tx) {
    warning("Might be unstable!")
    no_cores <- parallel::detectCores() - 1
    doParallel::registerDoParallel(no_cores)
    
    # Get annotation
    anno <- minfi::getAnnotation(CopyNumber450kData::RGcontrolSetEx)
    
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
    
    if (length(anno[, 1]) != length(ctrl) &&
        length(ctrl) != length(ctrlAll[, 1]) &&
        length(data[, 1] != length(ctrl))) {
        stopImplicitCluster()
        stop(
            "Something went terribly wrong :(
            Annotation & controls do not have the same dimensions!"
        )
    }
    
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
            p.val <- wilcox.test(da, ctAll)$p.value
            
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

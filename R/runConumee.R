#' @title
#' Analyze data with conumee
#'
#' @description
#' Analyzes the given 450k data with conumee, extracting 
#' either segments, bins or transcripts.
#'
#' @param data MethylSet: sample data
#' @param ctrl MethylSet: reference
#' @param what type of return data: segments, bins, transcripts
#' @param tx vector of transcripts of interests
#'
#' @return list with $data and $what, containing the CNV.write conumee outputs
#'
#' @import conumee
#' @import plyr
#' @import IlluminaHumanMethylation450kanno.ilmn12.hg19
#'
#' @export
#' 
#' @examples 
#' require(IlluminaHumanMethylation450kanno.ilmn12.hg19)
#' data <- minfi::preprocessRaw(minfiData::RGsetEx[,1:2])
#' ctrl <- minfi::preprocessRaw(minfiData::RGsetEx[,3:4])
#' runConumee(data,ctrl,"segments")
runConumee <- function(data,
                        ctrl,
                        what = "segments",
                        tx = NULL) {
    if (what == "segments") {
        print("## Extracting segments ...")
    } else if (what == "bins") {
        print("## Extracting bins ...")
    } else if (what == "transcripts") {
        print("## Extracting transcripts ...")
    } else {
        stop("Wrong call!")
    }
    
    ## check if all datasets belong to the same array type
    if (cnAnalysis450k::determineArrayType(data) != 
        cnAnalysis450k::determineArrayType(ctrl)) {
        warning(paste("Provided data probably not from the same array: ", 
                    cnAnalysis450k::determineArrayType(data),
                    cnAnalysis450k::determineArrayType(ctrl)))
    }
    
    if (what == "segments" | what == "bins") {
        anno.cnv <- CNV.create_anno(
            array_type=cnAnalysis450k::determineArrayType(data))
    } else if (what == "transcripts") {
        #build GRanges object with transcripts
        anno.cnv <- CNV.create_anno(
                        array_type=cnAnalysis450k::determineArrayType(data),
                        detail_regions = createGRangesObj(tx))
    }
    
    data.cnv <- CNV.load(data)
    ctrl.cnv <- CNV.load(ctrl)
    
    ##collect conumee data
    conumeeData <- NULL
    for (i in 1:length(names(data.cnv))) {
        print(paste("Run conumee analysis for ", names(data.cnv)[i], "..."))
        x <- CNV.fit(query=data.cnv[i], 
                    ref=ctrl.cnv[1:length(colnames(ctrl.cnv@intensity))], 
                    anno.cnv)
        x <- CNV.bin(x)
        
        if (what == "segments") {
            x <- CNV.segment(x)
            out <- CNV.write(x, what = "segments")
            conumeeData <- rbind.fill(conumeeData, out)
        } else if (what == "bins") {
            out <- CNV.write(x, what = "bins")
            if (length(conumeeData) == 0) {
                conumeeData <- out
            } else {
                conumeeData <- cbind(conumeeData, out[, 5])
                colnames(conumeeData)[length(conumeeData[1, ])] <-
                    names(data.cnv)[i]
            }
        } else if (what == "transcripts") {
            x <- CNV.detail(x)
            out <- CNV.write(x, what = "detail")
            if (length(conumeeData) == 0) {
                conumeeData <- out
            } else {
                conumeeData <- cbind(conumeeData, out[, 8])
            }
            colnames(conumeeData)[length(conumeeData[1, ])] <-
                names(data.cnv)[i]
        }
    }
    
    ret <- list(data = conumeeData, what = what)
    return(ret)
}


createGRangesObj <- function(tx) {
    txid <-
        AnnotationDbi::select(
        TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene,
        keys = tx,
        columns = c("TXNAME", "TXSTART", "TXEND", "TXCHROM"),
        keytype = "TXNAME"
        )
    txid <- txid[which(txid$TXCHROM %in% paste("chr", 1:22, sep = "")), ]
    txid <-
        txid[which(!duplicated(paste(
            txid$TXSTART, txid$TXEND, txid$TXCHROM
        ))), ]
    ##build granges object
    gr <- GenomicRanges::GRanges(
        S4Vectors::Rle(c(txid$TXCHROM)),
        IRanges::IRanges(txid$TXSTART, txid$TXEND),
        names = txid$TXNAME,
        genes = txid$TXNAME
    )
    gr$IRanges <- IRanges::ranges(gr)
    
    return(gr)
}

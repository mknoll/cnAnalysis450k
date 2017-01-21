#' @title Determines Array Type
#' 
#' @description
#' Tries to determine the array type (EPIC, 450k)
#' 
#' @param data object, for which the array type should 
#' be determined
#' 
#' @return "EPIC", "450k" or NULL if array type could not 
#' be determined
#' 
#' @export
#' 
#' @examples
#' determineArrayType(minfiData::MsetEx)
determineArrayType <- function(data) {
    array <- NULL

    ## 450k
    # MethylSet: 485512 rows
    # RGSet: 622399 rows

    # Determine by dimension
    if (length(unique(rownames(data))) > 622399) {
        #safe to assume it's EPIC
        array <- "EPIC"
    } else if (length(unique(rownames(data))) == 485512 
                || length(unique(rownames(data)) ==  622399)) {
        #might be a truncated EPIC assay or 450k
        if ("cg08795713" %in% rownames(data)) {
            # cg identifier on Chr9 unique for EPIC
            array <- "EPIC"
        } else {
            # can only be 450k
            array <-"450k"
        }
    } else {
        warning("Array Type could not be determined!")
    }
    return(array)
}

#' @title
#' Get Annotation Object 
#' 
#' @description
#' Returns the getAnnotation() object for 450k or EPIC data
#' 
#' @param type can be "450k" or "EPIC"
#' 
#' @return annotation object, retrieved via getAnnotation()
#' for the respective array design
#' 
#' @import minfiData
#' @import minfiDataEPIC
#' 
#' @export
#'
#' @examples
#' getAnnoData("EPIC")[1:3,]
getAnnoData <- function(type) {
    anno <- NULL
    if (type == "450k") {
        anno <- minfi::getAnnotation(minfiData::RGsetEx)
    } else if (type == "EPIC") {
        anno <- minfi::getAnnotation(minfiDataEPIC::RGsetEPIC)
    } else {
        stop("Unknown array type: can be 'EPIC' or '450k'!")
    }
    return(anno)
}


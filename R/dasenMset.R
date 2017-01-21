#' @title Adapted dasen function
#' 
#' @description
#' dasen function from the wateRmelon package, 
#' returns, however, a MethylSet. Code from bgeqot.R.
#' 
#' @param mns getMethyl() 
#' @param uns getUnmethyl()
#' @param onetwo probe Types ("I" or "II")
#' 
#' @return dasen normalized MethylSet
#' 
#' @import limma
#' @import minfi
#' @import wateRmelon
#' 
#' @export
#' 
#' @examples
#' m <- minfi::getMeth(minfiData::MsetEx[,1])
#' u <- minfi::getUnmeth(minfiData::MsetEx[,1])
#' ot <- minfi::getAnnotation(minfiData::MsetEx)$Type
#' dasenMset(m,u,ot)
dasenMset <- function(mns, uns, onetwo) {
    mnsc <- wateRmelon::dfsfit(mns,  onetwo)  
    unsc <- wateRmelon::dfsfit(uns,  onetwo, roco=NULL)
    mnsc[onetwo=='I' ,] <- limma::normalizeQuantiles(mnsc[onetwo=='I', ])
    mnsc[onetwo=='II',] <- limma::normalizeQuantiles(mnsc[onetwo=='II',])

    unsc[onetwo=='I' ,] <- limma::normalizeQuantiles(unsc[onetwo=='I', ])
    unsc[onetwo=='II',] <- limma::normalizeQuantiles(unsc[onetwo=='II',])
    
    mset <- minfi::MethylSet(mnsc,unsc)
    return(mset)
}

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
#' @import wateRmelon
#' @import minfi
#' 
#' @export
#' 
#' @examples
#' m <- data.frame(
#' smp1=c(22041, 679, 1620, 449),
#' smp2=c(588, 569, 421, 614),
#' smp3=c(20505, 439, 707, 343)
#' )
#' rownames(m) <- c(
#' "cg00050873", "cg00212031",
#' "cg00213748", "cg00214611"
#' )
#' u <- data.frame(
#' smp1=c(1945, 6567,384, 4869),
#' smp2=c(433,300,461,183),
#' smp3=c(1012,2689,295,1655)
#' )
#' rownames(u) <- c(
#' "cg00050873", "cg00212031",
#' "cg00213748", "cg00214611"
#' )
#' ot <- rep("I", 4)
#' dasenMset(m,u,ot)
dasenMset <- function(mns, uns, onetwo) {
  mnsc <- dfsfit(mns,  onetwo)  
  unsc <- dfsfit(uns,  onetwo, roco=NULL)
  mnsc[onetwo=='I' ,] <- normalizeQuantiles(mnsc[onetwo=='I', ])
  mnsc[onetwo=='II',] <- normalizeQuantiles(mnsc[onetwo=='II',])
  
  unsc[onetwo=='I' ,] <- normalizeQuantiles(unsc[onetwo=='I', ])
  unsc[onetwo=='II',] <- normalizeQuantiles(unsc[onetwo=='II',])
  
  mset <- minfi::MethylSet(mnsc,unsc)
  return(mset)
}

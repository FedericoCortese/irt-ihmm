#' Identify best clustering of states in PSBP-iHMM model based on dslo method in salso package 
#'
#' @param ihmm object of type "ihmm"
#'
#' @importFrom salso dlso
#' @return "best" clustering of hidden states z
#' @export
#'
#' 


## this code assumes each time series has equal length ## 

dlso_wrapper <- function(ihmm){
  Z_keep <- ihmm$z.save
  K <- max(unlist(ihmm$K.save))
  Z_newlist <- list()
  for(s in 1:length(Z_keep)){
    Z_newlist[[s]] <- unlist(Z_keep[[s]])
  }
  n <- length(Z_newlist[[1]]) # how many time points to cluster 
  niter <- length(Z_newlist)
  Z_mat <- matrix(unlist(Z_newlist), ncol = n, nrow = niter, byrow = TRUE)
  Zbest = dlso(Z_mat, loss="VI")
  numSub <- length(Z_keep[[1]]); numSub
  t.max <- n/numSub; t.max
  splits <- seq(1,numSub*t.max, t.max)
  z.list <- lapply(1:length(splits), FUN = function(i) Zbest[splits[i]:(splits[i]+t.max-1)])
  # z.list is a list of best clusters for each person 
  return(z.list)
}

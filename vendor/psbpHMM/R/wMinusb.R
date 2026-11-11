#' Calculate W minus Xbeta
#'
#' @param i individual
#' @param t time point
#' @param k state
#' @param w.z w
#' @param beta.k overall beta 
#' @param beta.sk subject specific beta
#' @param X covariates 
#'
#' @return w minus covariate effects
#'

wMinusb <- function(i, t, k, w.z, beta.k, beta.sk=NULL, X){
  if(!is.null(beta.sk)){
    return(w.z[[i]][[t]][k] - crossprod(X[t,], beta.k[[k]]) - crossprod(X[t,], beta.sk[[k]]))
  }else{
    return(w.z[[i]][[t]][k] - crossprod(X[t,], beta.k[[k]]))
  }
}

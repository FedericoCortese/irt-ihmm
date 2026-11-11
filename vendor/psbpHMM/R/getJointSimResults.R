#' Summarize simulation results from jointly fit iHMMs 
#'
#' @param fit1 object of type iHMM
#'
#' @return a list with components 
#' \itemize{
#'     \item hamming: hamming distance
#'     \item mu.mse: MSE for state-specific means
#'     \item mar.mse: MSE for MAR imputations 
#'     \item lod.mse: MSE for imputations below LOD
#' }  
#' @export
#'

getJointSimResults <- function(fit1){
 
  ham1 = mean(fit1$hamming)
  mu.mse1 = mean(fit1$mu.mse)
  mar.mse1 = mean(fit1$mar.mse)
  lod.mse1 = mean(fit1$lod.mse)
  
  list1 = list(hamming = ham1, mu.mse = mu.mse1, mar.mse = mar.mse1, lod.mse = lod.mse1)
  
  return(list1)
  
}
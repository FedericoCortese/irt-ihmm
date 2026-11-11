#' Results from joint cyclical model demonstration for package vignette
#'
#' an object of type "ihmm"
#'
#' @docType data
#'
#' @usage data(fit_cyclical)
#'
#' @format A list with components 
#' \describe{
#'        \item {z.save} {list of estimated hidden states for each time series at each iteration}
#'        \item {K.save} {list of estimated number of hidden states for each time series at each iteration}
#'        \item {ymar} {matrix with len.imp rows of imputed values for MAR data}
#'        \item {ylod} {matrix with len.imp rows of imputed values for data below LOD}
#'        \item {beta.save} {list of posterior estimates of beta_k, state-specific regression coefficients in PSBP} 
#'        \item {gamma.save} {NULL}
#'        \item {mu.save} {list of posterior estimates of mu_k, state-specific means} 
#'        \item {hamming} {posterior hamming distance between true and estimated states, if z.true is given}
#'        \item {mu.mse} {mean squared error for estimated state-specific means, if mu.true is given}
#'        \item {mu.sse} {sum of squared errors for estimated state-specific means, if mu.true is given}
#'        \item {mar.mse} {mean squared error of MAR imputations, if ycomplete is given} 
#'        \item {lod.mse} {mean squared error of imputations below LOD, if ycomplete is given}
#'        \item {mar.sse} {sum of squared errors of MAR imputations, if ycomplete is given} 
#'        \item {lod.sse} {sum of squared errors of imputations below LOD, if ycomplete is given}
#'        \item {mismat} {list, each element is a matrix indicating types of missing data for each time series, 0 = complete, 1 = MAR, 2 = below LOD}
#'        \item {ycomplete} {list of complete data}
#'        \item {MH.arate} {MH acceptance rate for lower triangular elements}
#'        \item {MH.lamrate} {MH acceptance rate for diagonal elements} 
#' }
#'
#' @examples 
#' data(fit_cyclical)
#'
#'
#'
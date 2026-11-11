#' Calculate ALL modes in a vector 
#'
#' @param x vector of factors/integers
#'
#' @return all modes of x 
#'

Mode <- function(x) {
  ux <- unique(x)
  return(ux[which(tabulate(match(x, ux))==max(tabulate(match(x, ux))))])
}

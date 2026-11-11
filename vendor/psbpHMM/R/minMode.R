#' Min Mode: Find the least common value 
#'
#' @param x vector of states 
#'
#' @return least common value 

minMode <- function(x){
  ux <- unique(x)
  return(ux[which(tabulate(match(x, ux))==min(tabulate(match(x, ux))))])
}
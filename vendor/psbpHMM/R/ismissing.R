#' function to indicate missing values 
#'
#' @param x data
#'
#' @return matrix same size as x, recoded to denote missing type. 0 = observed, 1 = NA (missing at random), 2 = -Inf (below LOD)
#'

ismissing <- function(x) {
  if(is.na(x)) return(1) # MAR
  else if(x==-Inf) return(2) # lod
  else return(0) 
  } # observed

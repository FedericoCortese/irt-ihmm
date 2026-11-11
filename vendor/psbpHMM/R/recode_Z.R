#' Relabel z to 1:length(z) while retaining order of z
#'
#' @param z list of states
#'
#' @return list of categorical indicators re-labeled from 1 to the number of unique clusters

recode_Z <- function(z){
  levels <- sort(unique(z))
  for (i in 1:length(levels)){
    z[which(z == levels[i])] <- i
  }
  return(z)
}
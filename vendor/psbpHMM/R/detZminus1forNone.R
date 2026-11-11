#' Determine possible states when no covariates
#'
#' @param i sampling day
#' @param sli state list 
#' @param pizi pi list
#' @param ui u list
#' @param t.max length of time series
#'
#' @return list of current and previous possible states
#'

detZminus1forNone <- function(i, sli, piz, ui, t.max){
  
  return(mclapply(2:t.max, FUN = function(t){
    pass1 <- sli[[t-1]] # possible states for t-1
    pos.zt <- sli[[t]] # possible states for t
    # need loop so we get NULL values where we need them 
    pos.ztminus1 <- list()
    for(k in pos.zt){
      pass2 <- which(piz[[2]][, k] >= ui[t]) # second/backward condition
      pos.ztminus1[[k]] <- intersect(pass1, pass2) # these are the possible values of z_(t-1) given that z_t = k
    }
    return(list(pos.zt = pos.zt, pos.ztminus1 = pos.ztminus1))
  }))
  # pos.zt are possibles states at time t
  # pos.ztminus1 are the specific possible states at t-1 for each k in pos.zt
  # pos.ztminus1 tells me which states to sum over when I calculate prob for z[t] = k
} 
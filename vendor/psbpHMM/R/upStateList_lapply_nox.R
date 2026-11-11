#' update state list when no covariates
#'
#' @param i sampling day
#' @param u list of u
#' @param pi.z list of pi
#' @param K number of states
#' @param t.max time series length
#'
#' @return udpated state list
#'
upStateList_lapply_nox = function(i,u,pi.z,K,t.max){
  for.list <- list()
  for.list[[1]] <- which(pi.z[[1]] >= u[[i]][1])
  for(t in 2:t.max){
    
    for.list[[t]] <- sort(unique(unlist(lapply(intersect(for.list[[t-1]], 1:K),
                                               FUN = function(k) which(pi.z[[2]][k,] >= u[[i]][t])))))
  }
  # backward condition
  back.list <- list()
  back.list[[t.max]] <- for.list[[t.max]]
  for(t in (t.max-1):1){
    back.list[[t]] <- sort(unique(unlist(lapply(intersect(back.list[[t+1]], 1:K),
                                                FUN = function(k) which(pi.z[[2]][,k] >= u[[i]][t+1])))))
  }
  return(lapply(1:t.max, FUN = function(t) {
    return(intersect(for.list[[t]], back.list[[t]]))
  }))
}
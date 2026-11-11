#' Summarize simulation results from indpendently fit iHMMs 
#'
#' @param fit1 list of objects of type iHMM 
#' @param y observed data 

#'
#' @return a list with components 
#' \itemize{
#'     \item hamming: hamming distance
#'     \item mu.mse: MSE for state-specific means
#'     \item mar.mse: MSE for MAR imputations 
#'     \item lod.mse: MSE for imputations below LOD
#' }  
#' 
#' @export
#'
getIndepSimRsults <- function(fit1, y){
  
  t.max = nrow(y[[1]])
  n = length(fit1) # number of independently fit time series 
  
  lodmiss <- sum(unlist(lapply(1:n, FUN = function(i){
    lmiss <- length(which(y[[i]]==-Inf))
    return(lmiss)
  })))
  
  marmiss <- sum(unlist(lapply(1:n, FUN = function(i){
    mmiss <- length(which(is.na(y[[i]])))
    return(mmiss)
  })))
  
  # hamming distance
  ham1 <- unlist(mclapply(1:n, FUN = function(i){
    if(!anyNA(fit1[[i]]$hamming)) return(fit1[[i]]$hamming)
    else return(NULL)
  }))
  
  
  # mu MSE 
  mu.sse1 <- mclapply(1:n, FUN = function(i){
    if(!anyNA(fit1[[i]]$mu.sse)) return(fit1[[i]]$mu.sse)
    else return(NULL)
  })
  muMat <- numeric() # matrix, number of iterations by number of time series 
  for(i in 1:n){
    muMat <- cbind(muMat, mu.sse1[[i]])
  }
  mu.mse1 <- apply(muMat,1, sum)/(n*t.max) # sum over the number of time series before dividing 

  # MAR MSE 
  if(marmiss > 0){
    mar.sse1 <- mclapply(1:n, FUN = function(i){
      if(!anyNA(fit1[[i]]$mar.sse)) return(fit1[[i]]$mar.sse)
      else return(NULL)
    })
    marMat <- numeric()
    for(i in 1:n){
      marMat <- cbind(marMat, mar.sse1[[i]])
    }
    mar.mse1 <- apply(marMat,1,sum)/marmiss # sum the rows to get total sse, then divide by total mar data
  }else{
    mar.mse1 = NULL
  }
  
  # LOD MSE 
  if(lodmiss > 0 ){
    lod.sse1 <- mclapply(1:n, FUN = function(i){
      if(!anyNA(fit1[[i]]$lod.sse)) return(fit1[[i]]$lod.sse)
      else return(NULL)
    })
    lodMat <- numeric()
    for(i in 1:n){
      lodMat <- cbind(lodMat, lod.sse1[[i]])
    }
    lod.mse1 <- apply(lodMat,1,sum)/lodmiss # sum the rows to get total sse, then divide by total mar data
  }else{
    lod.mse1 = NULL
  }
  
  list1 = list(hamming = mean(ham1), mu.mse = mean(mu.mse1), mar.mse = mean(mar.mse1), lod.mse = mean(lod.mse1))
  
  return(list1)
  
}

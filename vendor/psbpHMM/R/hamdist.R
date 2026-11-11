#' Greedily assign estimated states to true states then calculate Hamming distance
#'
#' @param z1 vector of true states
#' @param z2 vector of estimated states 
#'
#' @return hamming distance for most overlap
#'


hamdist <- function(z1, z2){
  
  ## we don't want to replace a state we have already replaced ##
  z1 <- recode_Z(z1) # recode the true states to be 1:K (they should already be like this)

  K <- length(unique(z1))
  x <- list()
  n <- matrix(NA,K,2)
  # loop through true states and find the best est. states to change to k 
  for(k in 1:K){
    x[[k]] <- z2[which(z1==k)]
    n[k,] <- cbind(length(which(x[[k]]==min(Mode(x[[k]])))),k) # how big is the mode, which state
  }
  
  N <- matrix(n[order(n[,1], decreasing = TRUE),],ncol=2)
  # we want to go in the order of N[,2]
  z.temp <- numeric(length(z2))
  idx <- 1
  for(k in N[,2]){
    # values we can't use 
    cantuse <- unique(z2[which(z.temp!=0)]) 
    if(any(x[[k]] %in% cantuse)){
      posz <- x[[k]][-which(x[[k]] %in% cantuse)] 
    }else{
      posz <- x[[k]]
    }
    findz <- Mode(posz)
    # 3 cases:
    # no findz: leave findz at 0, already used everything good
    # only 1 findz: set findz to k
    # more than one findz: choose a findz that minimizes overlap with other states
    if(length(findz)>1){
      if(!any(findz %in% unlist(x[N[-(1:idx),2]]))){
        findz <- findz[1] # doesn't matter which one you use, equally good 
      }else{
        findz <- minMode(unlist(x[N[-(1:idx),2]])) # least overlap with other states not yet used 
      }
    }
    z.temp[which(z2==findz[1])] <- k  # state to replace k with, if findz still > 1 then just use the first one 
    if(! 0 %in% z.temp) break 
    idx <- idx +1
  }
  #return(list(z.temp = z.temp, ham = sum(z1!=z.temp)))
  return(sum(z1!=z.temp))
  
}






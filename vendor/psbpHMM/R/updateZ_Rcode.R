#' Update hidden states
#'
#' @param i sampling day
#' @param yi y
#' @param mu mu
#' @param cholSigma chol Sigma
#' @param detR.stari detRstar
#' @param nu.df nu df
#' @param K number of states
#' @param log.stuff log stuff constant
#' @param t.max length of time series
#' @param detz possible states for t-1
#' @param sli state lists
#'
#' @return updated states
#'

updateZ_Rcode <- function(i, yi, mu, cholSigma, detR.stari, nu.df, K, log.stuff, t.max, detz, sli){
  z.probs <- list()
  # calculate probability for the first time point 
  k.prime <- sli[[1]] # possible states for first time point via slice sampling 
  log.tk.probs <- sapply(intersect(k.prime, 1:K), FUN = function(k) {
    dmvn(yi[1,], mu[[k]], cholSigma[[k]], log = TRUE, isChol = TRUE) # current states
  })
  if((K+1)%in%k.prime){ # new state
    log.tk.probs <- c(log.tk.probs,log.stuff-((nu.df+1)/2)*log(detR.stari[1]))
  }
  tk.probs <- exp(log.tk.probs-logsum(log.tk.probs))
  z.probs[[1]] <-  data.frame(k.prime, tk.probs)
  
  # this has to be a loop because it depends on previous time
  for(t in 2:t.max){ 
    k.prime <- sli[[t]]; k.prime
    j.prime <- z.probs[[t-1]]$k.prime
    j.probs <- z.probs[[t-1]]$tk.probs
    # for each current state k that is possible, sum over the previous possible states and multiply by density of k 
    log.tk.probs <- sapply(k.prime, FUN = function(k){
      log.sumprobs <- log( sum( (j.prime %in% detz[[t-1]]$pos.ztminus1[[k]])*j.probs)) # sum over probs for previous time point 
      if(k <= K){
        loglik_k <- dmvn(yi[t,], mu[[k]], cholSigma[[k]], log = TRUE, isChol = TRUE) # current state
      }else{
        loglik_k <- log.stuff-((nu.df+1)/2)*log(detR.stari[t])
      }
      return(loglik_k + log.sumprobs)
    })
    tk.probs <- exp(log.tk.probs - logsum(log.tk.probs))
    z.probs[[t]] <- data.frame(k.prime, tk.probs)
  }
  
  upz <- mclapply(1:t.max, FUN = function(t){
    k.prime <- z.probs[[t]]$k.prime
    if(length(k.prime) == 1){
      statez <- k.prime
    }else{
      statez <- sample(k.prime, 1, prob = z.probs[[t]]$tk.probs)
    }
    return(statez) 
  })
  return(unlist(upz))
}

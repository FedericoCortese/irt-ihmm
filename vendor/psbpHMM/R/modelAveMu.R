#' Calculate model averaged estimates of state-specific means from a PSBP iHMM fit using the psbpHMM package
#'
#' @param fit object of type "ihmm"
#' @param zbest1 a list of optimal clustering of hidden states
#' @param ymatrix matrix of exposure data for all time series, with n*T rows and p columns 
#'
#' @return a list with components 
#' \itemize{
#'    \item mu_ma: posterior mean of model-averaged state-specific means
#'    \item mu_lwr: .025 quantile of posterior distribution of model-averaged state-specific means
#'    \item mu_upr: .975 quantile of posterior distribution of model-averaged state-specific means
#'    \item exp_lwr: minimum empirical data for each exposure within each state 
#'    \item exp_upr: maximum empirical data for each exposure within each state 
#' } 
#' @export
#'
#' 

modelAveMu <- function(fit, zbest1, ymatrix){
  # model averaged estimates of the harmonic trend (beta only)
  K <- max(unlist(zbest1)) # number of best clusters 
  p <- dim(fit$mu.save[[1]][[1]])[2]
  niter <- length(fit$mu.save);niter
  mu_ma <- list() # model averaged estimate of mu
  mu_upr <- list()
  mu_lwr <- list()
  exp_lwr <- list()
  exp_upr <- list()

  for(k in 1:K){
    whoBest <- which(unlist(zbest1)==k) # the time points assigned to k in the best clustering 
    
    exp_lwr[[k]] = apply(ymatrix[whoBest,], 2, FUN = function(y){
      min(y[which(y>-Inf)], na.rm = TRUE)
    })
    
    exp_upr[[k]] = apply(ymatrix[whoBest,], 2, FUN = function(y){
      max(y[which(y>-Inf)], na.rm = TRUE)
    })
    
    mu.ave <- matrix(NA, p, niter)
    for(s in 1:niter){
      these.states <- unlist(fit$z.save[[s]])[whoBest] # average beta over these states
      mu.matrix <- matrix(0, p, 1) # place holder 
      for(j in 1:length(these.states)){
        mu.matrix <- cbind(mu.matrix, unlist(fit$mu.save[[s]][these.states[j]]))
      }
      mu.ave[,s] <- apply(matrix(mu.matrix[,-1], nrow = p), 1, mean)
    }
    
    # distribution of X%*%beta 
    mu_ma[[k]] <- apply(mu.ave, 1, mean)
    mu_lwr[[k]] <- apply(mu.ave, 1, FUN = function(x) quantile(x, .05))
    mu_upr[[k]] <- apply(mu.ave, 1, FUN = function(x) quantile(x, .95))
    
  }
  return(list(mu_ma = mu_ma, mu_lwr = mu_lwr, mu_upr = mu_upr, exp_lwr = exp_lwr, exp_upr = exp_upr))
}

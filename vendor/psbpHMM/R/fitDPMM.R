#' Fit truncated joint DPMM to multiple time series 
#'
#' @param niter number of total iterations
#' @param nburn number of burn-in iterations
#' @param y list of time series data for each time series 
#' @param ycomplete complete data, if available, for evaluating imputations
#' @param priors list of priors
#' @param K.start maximum allowable number of states
#' @param z.true list of true hidden states, if known
#' @param lod list of lower limits of detection for p exposures for each time series
#' @param mu.true matrix of true exposure means for each true state, if known 
#' @param missing logical; if TRUE then the data set contains missing data, default is FALSE
#' @param tau2 variance tuning parameter for normal proposal in MH update of lower triangular elements in decomposition of Sigma
#' @param a.tune shape tuning parameter for inverse gamma proposal in MH update of diagonal elements in decomposition of Sigma
#' @param b.tune rate tuning parameter for inverse gamma proposal in MH update of diagonal elements in decomposition of Sigma
#' @param resK logical; if TRUE a resolvent kernel is used in MH update for lower triangular elements in decomposition of Sigma
#' @param eta.star resolvent kernel parameter, must be a real value greater than 1. In the resolvent kernel we take a random draw from the geometric distribution with mean (1-p)/p, eta.star = 1/p.
#' @param len.imp number of imputations to save. Imputations will be taken at equally spaced iterations between nburn and niter. 
#' @param holdout list of indicators of missing type in holdout data set, 0 = observed, 1 = MAR, 2 = below LOD, for imputation validation purposes
#'
#' @importFrom parallel mclapply
#' @importFrom stats rnorm runif rgamma rWishart cov cov2cor dnorm rgeom rbeta 
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom tmvmixnorm rtmvn
#' @importFrom matrixcalc matrix.trace
#' @importFrom mvnfast rmvn dmvn
#' @importFrom invgamma dinvgamma
#' @importFrom gdata lowerTriangle<-
#'
#' @return an object of type "dpmm"
#'
#' @return a list with components
#' \itemize{
#'        \item z.save: list of estimated hidden states for each time series at each iteration
#'        \item K.save: list of estimated number of hidden states for each time series at each iteration
#'        \item mu.save: list of posterior estimates of mu_k, state-specific means 
#'        \item ymar: matrix of imputed values for MAR data, number of rows equal to len.imp
#'        \item ylod: matrix of imputed values for data below LOD, number of rows equal to len.imp
#'        \item hamming: posterior hamming distance between true and estimated states, if z.true is given
#'        \item mu.mse: mean squared error for estimated state-specific means, if mu.true is given
#'        \item mar.mse: mean squared error of MAR imputations, if ycomplete is given 
#'        \item lod.mse: mean squared error of imputations below LOD, if ycomplete is given 
#'        \item mar.bias: mean bias for MAR imputations 
#'        \item lod.bias: mean bias for below LOD imputations 
#'        \item mismat: list, each element is a matrix indicating types of missing data for each time series, 0 = complete, 1 = MAR, 2 = below LOD
#'        \item ycomplete: complete data
#'        \item MH.arate: average MH acceptance rate for lower triangular elements
#'        \item MH.lamrate: average MH acceptance rate for diagonal elements 
#' }
#' @export

fitDPMM <- function(niter, nburn, y, ycomplete=NULL,
                          priors=NULL, K.start=NULL, z.true=NULL, lod=NULL,
                          mu.true=NULL, missing = FALSE, 
                          tau2 = NULL, a.tune = NULL, b.tune = NULL,
                          resK = FALSE, eta.star = NULL, len.imp = NULL,
                          holdout = NULL){
  

  # catch problems with parameter input 
  if(niter <= nburn) stop("niter must be greater than nburn")
  if(nburn < 0) stop("nburn must be greater than or equal to 0")

  
  if(!is.numeric(unlist(y))) stop("y must be numeric")
  if(!is.numeric(unlist(ycomplete)) & !is.null(ycomplete)) stop("ycomplete must be numeric")
  
  
  
  # the only prior option here is the jointNIW 
  
  #####################
  ### Initial Setup ###
  #####################
  
  
  if(missing){
    algorithm = "MH"
    SigmaPrior = "wishart"
  }else{
    algorithm = "Gibbs"
    #SigmaPrior = "non-informative"
    SigmaPrior = "wishart"
  }

  
  # how many time series and exposures 
  if(class(y)=="list"){
    p <- ncol(y[[1]]) # number of exposures 
    n <- length(y) # number of time series
    t.max <- nrow(y[[1]]) # number of time points 
  }else if(class(y)=="matrix"){
    p <- ncol(y)
    n <- 1
    t.max <- nrow(y)
    y <- list(y) # make y into a list 
    ycomplete <- list(ycomplete)
    z.true <- list(z.true)
  }else if(class(y)=="numeric"){
    p <- 1
    n <- 1
    t.max <- length(y)
    y <- list(matrix(y, ncol = 1))
    ycomplete <- list(matrix(ycomplete, ncol= 1))
    z.true <- list(z.true)
  }

  
  ##############
  ### Priors ###
  ##############
  
  if(missing(priors)) priors <- NULL
  # mu
  if(is.null(priors$mu0)) priors$mu0 <- matrix(0, p, 1) # mu_k|Sigma_k ~ N(mu0, 1/lambda*Sigma_k) prior mean on p exposures
  # if SigmaPrior = "wish", inverse wishart dist with fixed hyperparams on Sigma_k
  if(SigmaPrior == "wishart"){
    if(is.null(priors$R)) priors$R <- diag(p) # Sigma_k ~ Inv.Wish(nu, R) hyperparameter for Sigma_k
    if(is.null(priors$nu)) priors$nu <- p+2 # Sigma_k ~ Inv.Wish(nu, R) hyperparameter for Sigma_k, nu > p+1
    nu.df <- priors$nu # 
    R.mat <- priors$R 
  }else{
    # if SigmaPrior = "ni", the half-t prior on Sigma_k
    if(is.null(priors$bj)) priors$bj <- rep(1, p) # Huang and Wand advise 10e5
    if(is.null(priors$nu)) priors$nu <- 2 # Huang and Wand advise 2, p+4 for so the variance exists
    nu.df <- priors$nu + p - 1 # Huang and Wand, prior df 
    aj.inv <- rgamma(p, shape = 1/2, rate = 1/(priors$bj^2)) # these are 1/aj 
    # starting value for R.mat, the matrix parameter on Sigma.Inv
    R.mat <- 2*priors$nu*diag(aj.inv) 
    # we model aj.inv with gamma(shape = 1/2, rate = 1/(b_j^2))
  }
  if(is.null(priors$lambda)) priors$lambda <- 1 # concentration on Sigma_k for NIW 
  
  #############################
  ### Indicate Missing Type ###
  #############################
  
  # indicate missing data: obs = 0, mar = 1, lod = 2
  mismat <- list()
  for(i in 1:n){
    mismat[[i]] <- matrix(sapply(y[[i]], ismissing), ncol = p)
  }
  
  if(is.null(holdout)) holdout = mismat
  
  mism <- numeric()
  for(i in 1:n){
    mism <- rbind(mism, mismat[[i]])
  }
  
  # for each i, which time points have any missing data? 
  missingTimes <- lapply(1:n, FUN = function(i){
    which(apply(mismat[[i]],1,sum)>0)
  })
  
  # for each i, which times points are observed? 
  observedTimes <- lapply(1:n, FUN = function(i){
    which(apply(mismat[[i]], 1, sum)==0)
  })
  
  ###############################################
  ### Impute Starting Values for Missing Data ###
  ###############################################
  
  for(i in 1:n){
    if(any(mismat[[i]]==2)){ # lod 
      expLod <- exp(lod[[i]])
      numlod <- apply(mismat[[i]], 1, FUN = function(x) length(which(x==2))) # how many LOD
      for(t in which(numlod>0)){ # loop thru LOD data
        whichlod <- which(mismat[[i]][t,]==2) # which exposures are below LOD 
        y[[i]][t,whichlod] <- log(expLod[whichlod]/sqrt(2)) # impute the LOD with the log(LOD/sqrt(2))
      }
    }
    if(any(mismat[[i]]==1)){ # mar
      nummis <- apply(mismat[[i]], 1, FUN = function(x) length(which(x==1))) # how many missing at each time point 
      for(t in which(nummis>0)){ # only loop thru time points with MAR
        whichmis <- which(mismat[[i]][t,]==1) # which exposures are missing at random 
        if(t == 1){
          for(ws in whichmis){
            lastT = max(which(!is.na(y[[i]][,ws])))
            y[[i]][t,ws] <- y[[i]][lastT, ws] # fill in with the last observed value from the end of the time series
          }
        }else{
          y[[i]][t,whichmis] <- y[[i]][t-1, whichmis] # fill in the missing by LVCF
        }
      }
    }
  } 
  
  #####################
  ## Starting Values ##  
  #####################
  
  if(is.null(K.start)) K = 50 else K = K.start
    
  ### test
  # K.star =20
  # 
  # z <- lapply(1:n, FUN = function(i) {
  #   rcat(matrix(runif(t.max*K.star),t.max,K.star)) 
  # })# initial category assignment  
  
  
  z <- lapply(1:n, FUN = function(i) {
    rcat(matrix(runif(t.max*K),t.max,K))
  })# initial category assignment
  
  
  mu <- list()
  Sigma <- list()
  D <- list()
  L <- list()
  lams <- list()
  al <- list()
  
  ymatrix <- NULL
  for(i in 1:n){
    ymatrix <- rbind(ymatrix, y[[i]])
  }
  
  for(k in 1:K){
    if(missing){
      # we reparameterize Sigma and model L and D instead 
      vj0 <- sapply(1:p, FUN = function(j) priors$nu + j - p); vj0 # fixed for each k 
      deltaj0 <- rep(1,p); deltaj0 # fixed for each k 
      lams[[k]] <- 1/rgamma(3, vj0, rate = deltaj0)
      D[[k]] <- diag(lams[[k]])
      al.list <- list()
      for(j in 2:p){
        al.list[[j-1]] <- rnorm(j-1, 0, lams[[k]][j])
      }
      al[[k]] <- unlist(al.list) # for j = 2 to p
      which.lams <- unlist(sapply(2:p, FUN = function(j) rep(j,j-1))) # which lams to use for each al 
      L[[k]] <- diag(p)
      lowerTriangle(L[[k]]) <- al[[k]]
      Sigma[[k]] <- solve(L[[k]])%*%D[[k]]%*%t(solve(L[[k]])) # Sigma[[k]]
      mu[[k]] <- rmvn(1, priors$mu0, (1/priors$lambda)*Sigma[[k]])
    }else{
      Sigma[[k]] <- chol2inv(chol(matrix(rWishart(1, df = nu.df, Sigma = solve(R.mat)),p,p)))
      mu[[k]] <- rmvn(1, priors$mu0, (1/priors$lambda)*Sigma[[k]])
    }
  }

  ##################################
  ### for stick-breaking process ###
  ##################################
  
  alpha <- 1 # DP parameter
  psi <- rep(1/K, K) # cluster weights
  V <- rbeta(K, 1, alpha) # conditional weights
  V[K] <- 1 # apply DP truncation to C classes
  
  loglik <- list() # log-likelihood of being in each cluster
  for(i in 1:n){
    loglik[[i]] <- matrix(NA, t.max, K)
  }

  nc <- unlist(lapply(1:K, FUN = function(k)  length(which(unlist(z) == k))))
  
  
  ################################################
  ### fixed values for new state probabilities ###
  ################################################
  
  gampp <- mgamma(nu = nu.df, p = p)
  # fixed for "wishart", starting value for "non-informative" because R.mat will change with each iteration as aj.inv changes
  # if "ni" then need to update detR.star and log.stuff each iteration
  detR.star <- mclapply(1:n, FUN = function(i){
    sapply(1:t.max, FUN = function(t){
      x <- R.mat + priors$lambda*tcrossprod(priors$mu0) + tcrossprod(y[[i]][t,]) - 
        (1/(1+priors$lambda))*tcrossprod(priors$lambda*priors$mu0+y[[i]][t,])
      return(det(x))})
  })
  log.stuff <- (p/2)*log(priors$lambda/(pi*(priors$lambda+1)))+log(gampp)+(nu.df/2)*log(det(R.mat)); log.stuff
  

  ########################
  ### For MH algorithm ###
  ########################
  
  # calculate before loop for MH 
  ymatrix <- numeric()
  for(i in 1:n){
    ymatrix <- rbind(ymatrix, y[[i]])
  }
  ycomp <- ymatrix[which(rowSums(mism)==0),]
  
  ##############################
  ### MCMC Storage and Setup ###
  ##############################
  
  z.save <- list()
  mu.save <- list()
  beta.save <- list()
  K.save <- numeric()
  ham <- numeric()
  mu.mse <- numeric()
  mu.sse <- numeric()
  MH.a <- 0 # for MH update a
  MH.lam <- 0 # for MH update lams
  s.save = 1
  
  # missing data sets
  if(!is.null(len.imp)){
    imputes <- ceiling(seq.int(nburn, niter, length.out = len.imp))
    y.mar.save <- matrix(NA, len.imp, length(which(unlist(mismat)==1)))
    y.lod.save <- matrix(NA, len.imp, length(which(unlist(mismat)==2)))
    mar.mse <- numeric()
    lod.mse <- numeric()
    mar.bias = numeric()
    lod.bias = numeric()
    miss.mse <- numeric()
    s.imp <- 1
  }else{
    imputes = 0
    y.mar.save <- NULL
    y.lod.save <- NULL
    mar.mse <- NULL
    lod.mse <- NULL
    mar.bias = NULL
    lod.bias = NULL
    miss.mse <- NULL
    s.imp <- NULL
  }
  
  
  ###############
  ### Sampler ###
  ###############

  for(s in 1:niter){
    
    ##################################
    ### update cluster assignments ###
    ##################################
    
    for(i in 1:n){
      for(k in 1:K){
        loglik[[i]][,k] <- log(psi[k]) +
          dmvnorm(y[[i]], mu[[k]], Sigma[[k]], log = T)
      }
      # calculate likelihood
      lik_k <- exp(loglik[[i]] - apply(loglik[[i]], 1, logsum))
      z[[i]] <- rcat(lik_k) 
    }
    
    nc <- unlist(lapply(1:K, FUN = function(k)  length(which(unlist(z) == k))))
    
    ##############################
    ### update cluster weights ###
    ##############################
    
    # update conditional weights for stick-breaking process
    for (k in 1:(K-1)) V[k] <- rbeta(1, nc[k] + 1, alpha + sum(nc[(k+1):K]))
    
    # handle overflow
    V[which(V == 1)] <- 1-10e-8
    V[K] <- 1
    
    # update mixing weights 
    psi[1] <- V[1]
    cumV <- cumprod(1-V)
    for (k in 2:K) psi[k] <- V[k]*cumV[k-1]
  
    ######################
    ### update theta_k ### 
    ######################
    
    # first update mu and Sigma 
    cholSigma <- lapply(1:K, FUN = function(k) chol(Sigma[[k]]))
    
    for(k in 1:K){
      itimes <- lapply(1:n, FUN = function(i)  which(z[[i]] == k))
      if(sum(unlist(itimes)) == 0){
        # update from prior 
        if(missing){
          # MH sample from prior 
          vj0 <- sapply(1:p, FUN = function(j) priors$nu + j - p); vj0 # fixed for each k 
          deltaj0 <- rep(1,p); deltaj0 # fixed for each k 
          lams[[k]] <- 1/rgamma(3, vj0, rate = deltaj0)
          D[[k]] <- diag(lams[[k]])
          al.list <- list()
          for(j in 2:p){
            al.list[[j-1]] <- rnorm(j-1, 0, lams[[k]][j])
          }
          al[[k]] <- unlist(al.list) # for j = 2 to p
          which.lams <- unlist(sapply(2:p, FUN = function(j) rep(j,j-1))) # which lams to use for each al 
          L[[k]] <- diag(p)
          lowerTriangle(L[[k]]) <- al[[k]]
          Sigma[[k]] <- solve(L[[k]])%*%D[[k]]%*%t(solve(L[[k]])) # Sigma[[k]]
          mu[[k]] <- rmvnorm(1, priors$mu0, (1/priors$lambda)*Sigma[[k]])
        }else{
          # Gibbs sample from prior 
          Sigma[[k]] <- chol2inv(chol(matrix(rWishart(1, df = nu.df, Sigma = solve(R.mat)),p,p)))
          mu[[k]] <- rmvn(1, priors$mu0, (1/priors$lambda)*Sigma[[k]])
        }
      }else{
        nkk.tilde  <- length(unlist(itimes)) # number in state k 
        y.list <- lapply(1:n, FUN = function(i) matrix(y[[i]][itimes[[i]],], ncol = p))
        yk <- numeric()
        for(i in 1:n){
          yk <- rbind(yk, y.list[[i]])
        }
        ybark <- matrix(apply(yk, 2, mean),p,1)
        nu_nk <- nu.df + nkk.tilde 
        if(algorithm == "Gibbs"){
          
          mu_nk <- (priors$lambda*priors$mu0 + nkk.tilde*ybark)/(priors$lambda+nkk.tilde) 
          lambda_nk <- priors$lambda + nkk.tilde 
          if(nkk.tilde == 1){
            M <- R.mat
          }else{
            M <- R.mat + (nkk.tilde-1)*cov(yk) 
          }
          Sigma_nk <- M + (priors$lambda*nkk.tilde)/(nkk.tilde + priors$lambda)*tcrossprod(ybark - priors$mu0) 
          Sigma[[k]] <- chol2inv(chol(matrix(rWishart(1,df=nu_nk, Sigma=invMat(Sigma_nk)),p,p)))
          mu[[k]] <- rmvn(n=1, mu=mu_nk, sigma=chol((1/lambda_nk)*as.matrix(Sigma[[k]], p, p)), isChol = TRUE)  
          
        }else if(algorithm == "MH"){
          # update a
          for(j in 1:length(al[[k]])){
            
            if(resK){
              eta <- rgeom(1, (1/eta.star)) + 1
            }else eta <- 1
            
            if(eta>0){
              for(m in 1:eta){
                al.star <- al[[k]]
                L.star <- L[[k]]
                
                #a.star <- rnorm(1, al[[k]][j], sqrt(tau2)); a.star # proposed value 
                a.star <- rnorm(1, 0, sqrt(tau2)); a.star # proposed value 
                
                al.star[j] <- a.star; al.star
                lowerTriangle(L.star) <- al.star; L.star
                
                SigmaStar <- mhDecomp(L.star, D[[k]]) # cppFunction
                
                #SigmaStar <- solve(L.star)%*%D[[k]]%*%t(solve(L.star)); SigmaStar # function of a.star
                
                # likelihoods
                da.curr <- sum(dmvn(yk, mu = mu[[k]], sigma = cholSigma[[k]], log = TRUE, isChol = TRUE)) +
                  dmvn(mu[[k]], mu = priors$mu0, sigma = chol((1/priors$lambda)*Sigma[[k]]), log = TRUE, isChol = TRUE)
                da.star <- sum(dmvn(yk, mu = mu[[k]], sigma = chol(SigmaStar), log = TRUE, isChol = TRUE)) +
                  dmvn(mu[[k]], mu = priors$mu0, sigma = chol((1/priors$lambda)*SigmaStar), log = TRUE, isChol = TRUE)
                
                # priors 
                pa.curr <- dnorm(al[[k]][j], 0, sqrt(lams[[k]][which.lams[j]]), log = TRUE)
                pa.star <- dnorm(a.star, 0, sqrt(lams[[k]][which.lams[j]]), log = TRUE)
                
                # proposals
                qa.curr <- dnorm(al[[k]][j], 0, sqrt(tau2), log = TRUE)
                qa.star <- dnorm(a.star, 0, sqrt(tau2), log = TRUE)
                
                mh1 <- pa.star + da.star + qa.curr; mh1
                mh2 <- pa.curr + da.curr + qa.star
                # catch error on da.curr
                
                ar <- mh1-mh2
                
                if(runif(1) < exp(ar)){
                  al[[k]][j] <- a.star
                  L[[k]] <- L.star # update this too, fxn of a.star
                  Sigma[[k]] <- SigmaStar # update this too, fxn of a.star
                  MH.a <- MH.a + 1 
                }
              }
            }
          }
          # update lams
          for(j in 1:p){
            D.star <- D[[k]]
            lam.star <- 1/rgamma(1, a.tune, rate = b.tune); lam.star # proposed value 
            D.star[j,j] <- lam.star
            
            SigmaStar <- mhDecomp(L[[k]], D.star) # cppFunction
            
            #SigmaStar <- solve(L[[k]])%*%D.star%*%t(solve(L[[k]]))
            
            # likelihoods
            dlam.curr <- sum(dmvn(yk, mu = mu[[k]], sigma = cholSigma[[k]], log = TRUE, isChol = TRUE)) + 
              dmvn(mu[[k]], mu = priors$mu0, sigma = chol((1/priors$lambda)*Sigma[[k]]), log = TRUE, isChol = TRUE)
            dlam.star <- sum(dmvn(yk, mu = mu[[k]], sigma = chol(SigmaStar), log = TRUE, isChol = TRUE)) +
              dmvn(mu[[k]], mu = priors$mu0, sigma = chol((1/priors$lambda)*SigmaStar), log = TRUE, isChol = TRUE)
            
            # priors
            plam.curr <- dinvgamma(lams[[k]][j], vj0[j]/2, rate = deltaj0[j]/2, log = TRUE)
            plam.star <- dinvgamma(lam.star, vj0[j]/2, rate = deltaj0[j]/2, log = TRUE)
            
            # proposal 
            qlam.curr <- dinvgamma(lams[[k]][j], a.tune, b.tune, log = TRUE)
            qlam.star <- dinvgamma(lam.star, a.tune, b.tune, log = TRUE)
            
            mh1 <- plam.star + dlam.star + qlam.curr; mh1
            mh2 <- plam.curr + dlam.curr + qlam.star; mh2
            ar <- mh1-mh2
            
            if(runif(1) < exp(ar)){
              lams[[k]][j] <- lam.star
              D[[k]] <- D.star # fxn of lam.star
              Sigma[[k]] <- SigmaStar # update this too, fxn of lam.star
              MH.lam <- MH.lam + 1 
            }
          }
          # update mu by Gibbs
          mu_nk <- (priors$lambda*priors$mu0 + nkk.tilde*ybark)/(priors$lambda+nkk.tilde) 
          lambda_nk <- priors$lambda + nkk.tilde 
          mu[[k]] <- rmvn(n=1, mu=mu_nk, sigma=chol((1/lambda_nk)*as.matrix(Sigma[[k]], p, p)), isChol = TRUE)  
        } # end if MH 
      }
    } # end sample theta  
    
    
    ###################################################################################
    ### if SigmaPrior == "non-informative": Update aj.inv, detR.star and log.stuff  ###
    ###################################################################################
    
    if(SigmaPrior == "non-informative"){
      # then update aj.inv for j = 1 to p 
      shape.aj <- (K*(priors$nu+p-1)+1)/2
      sigmajj <- lapply(1:K, FUN = function(k){
        diag(invMat(Sigma[[k]])) # cppFunction
      }) # diagonal elements of each Sigma.Inv_k
      diags <- numeric()
      for(k in 1:K){
        diags <- rbind(diags, sigmajj[[k]])
      }
      sumdiags <- apply(diags, 2, sum) # sum of diagonals of Sigma.Inv for k = 1 to K 
      rate.aj <- 1/(priors$bj^2) + priors$nu*sumdiags
      aj.inv <- rgamma(p, shape = shape.aj, rate = rate.aj) # these are 1/aj 
      R.mat <- 2*priors$nu*diag(aj.inv)
      
      detR.star <- mclapply(1:n, FUN = function(i){
        sapply(1:t.max, FUN = function(t){
          x <- R.mat + priors$lambda*tcrossprod(priors$mu0) + tcrossprod(y[[i]][t,]) - 
            (1/(1+priors$lambda))*tcrossprod(priors$lambda*priors$mu0+y[[i]][t,])
          return(det(x))
        })})
      
      log.stuff <- (p/2)*log(priors$lambda/(pi*(priors$lambda+1)))+log(gampp)+(nu.df/2)*log(det(R.mat))
    }
    
    #################################
    ### Sample New Missing Values ###
    #################################
    
    # Sample new MAR values conditional on observed data and imputed LOD data ###
    for(i in 1:n){
      if(any(mismat[[i]]==1)){ # MAR = 1 
        nummis <- apply(mismat[[i]], 1, FUN = function(x) length(which(x==1))) # how many missing at each time point 
        for(t in which(nummis>0)){ # only loop through time points with missing data
          whichmis <- which(mismat[[i]][t,]==1) # which ones are missing 
          if(length(whichmis)==p){
            y[[i]][t,] <- rmvn(1, mu[[z[[i]][t]]], chol(Sigma[[z[[i]][t]]]), isChol = TRUE)
          }else{
            y.obs <- y[[i]][t,-whichmis]
            mu.obs <- mu[[z[[i]][t]]][,-whichmis]
            mu.miss <- mu[[z[[i]][t]]][,whichmis]
            Sigma.obs <- matrix(Sigma[[z[[i]][t]]][-whichmis, -whichmis], p-length(whichmis), p-length(whichmis))
            Sigma.miss <- matrix(Sigma[[z[[i]][t]]][whichmis, whichmis], length(whichmis), length(whichmis))
            Sigma.obs.miss <-  matrix(Sigma[[z[[i]][t]]][-whichmis, whichmis], p-length(whichmis), length(whichmis))
            Sigma.miss.obs <- t(Sigma.obs.miss)
            Sigma.mis.obs.inv <- Sigma.miss.obs%*%solve(Sigma.obs)
            mu.mgo <- as.numeric(mu.miss + Sigma.mis.obs.inv%*%(y.obs - mu.obs))
            Sigma.mgo <- Sigma.miss + Sigma.mis.obs.inv%*%Sigma.obs.miss
            y[[i]][t,whichmis] <- rmvn(1, mu.mgo, chol(Sigma.mgo), isChol = TRUE)
          }
        }
      }
    }
    
    # Sample new LOD values conditional on observed data and imputed MAR data ###
    for(i in 1:n){
      if(any(mismat[[i]]==2)){ # LOD = 2
        numlod <- apply(mismat[[i]], 1, FUN = function(x) length(which(x==2))) # how many lod at each time point 
        for(t in which(numlod>0)){ # only loop through time points with missing data
          whichlod <- which(mismat[[i]][t,]==2) # which ones are below lod  
          if(length(whichlod)==p){ 
            # before, used int = y[[i]][t, whichlod], trying new thing June 30 2020
            y[[i]][t,] <- rtmvn(1, Mean = as.vector(mu[[z[[i]][t]]]), Sigma = Sigma[[z[[i]][t]]], lower = rep(-Inf, p),
                                upper = lod[[i]], int = y[[i]][t,], burn = 10, thin = 1)
          }else{
            y.obs <- y[[i]][t,-whichlod]
            mu.obs <- mu[[z[[i]][t]]][,-whichlod]
            mu.miss <- mu[[z[[i]][t]]][,whichlod]
            Sigma.obs <- matrix(Sigma[[z[[i]][t]]][-whichlod, -whichlod], p-length(whichlod), p-length(whichlod))
            Sigma.miss <- matrix(Sigma[[z[[i]][t]]][whichlod, whichlod], length(whichlod), length(whichlod))
            Sigma.obs.miss <-  matrix(Sigma[[z[[i]][t]]][-whichlod, whichlod], p-length(whichlod), length(whichlod))
            Sigma.miss.obs <- t(Sigma.obs.miss)
            Sigma.mis.obs.inv <- Sigma.miss.obs%*%chol2inv(chol(Sigma.obs))
            mu.mgo <- as.numeric(mu.miss + Sigma.mis.obs.inv%*%(y.obs - mu.obs))
            Sigma.mgo <- Sigma.miss + Sigma.mis.obs.inv%*%Sigma.obs.miss
            
            # before, used int = y[[i]][t, whichlod], trying new thing June 30 2020
            y[[i]][t,whichlod] <- rtmvn(1, Mean = mu.mgo, Sigma = Sigma.mgo, lower = rep(-Inf, length(whichlod)),
                                        upper = lod[[i]][whichlod], int = y[[i]][t, whichlod], burn = 10, thin = 1)
            
            # tryCatch({
            #   y[[i]][t,whichlod] <- rtmvn(1, Mean = mu.mgo, Sigma = Sigma.mgo, lower = rep(-Inf, length(whichlod)),
            #                               upper = lod[whichlod], int = lod[whichlod]-1, burn = 10, thin = 1) }, 
            #   error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
            
            
          }
        }
      }
    }

    
    #####################
    ### Store Results ###
    #####################
    
    if(s>nburn){
      
      ## Hamming distance ##
      if(!is.null(unlist(z.true))){
        ham.error <- hamdist(unlist(z.true), unlist(z)) 
        ham[s.save] <- ham.error/(n*t.max) # proportion of misplaced states
      }else{
        ham <- NULL
      }    
      
      mu.save[[s.save]] <- mu
      
      ## MSE for mu ##
      if(!is.null(mu.true)){
        sse <- list()
        for(i in 1:n){
          sse[[i]] <- sapply(1:t.max, FUN = function(t){
            as.numeric(crossprod(unlist(mu[z[[i]][t]]) - mu.true[z.true[[i]][t],]))
          })
        }
        mu.mse[s.save] <- mean(unlist(sse))/p # vector mse for mu, divide by # exposures 
      }else{
        mu.mse <- NULL
      }
      z.save[[s.save]] <- z
      K.save[s.save] <- length(unique(unlist(z)))
      #K.save[s.save] <- K 
      
      if(s%in%imputes){
        # imputed values for complete data sets 
        y.mar.save[s.imp,] <- unlist(y)[which(unlist(mismat)==1)] # mar imputations
        y.lod.save[s.imp,] <- unlist(y)[which(unlist(mismat)==2)] # lod imputations
        
        if(!is.null(ycomplete)){
          # separate by the types of missing MSE
          mar.mse[s.imp] <- mean((unlist(ycomplete)[which(unlist(holdout)==1)] - unlist(y)[which(unlist(holdout)==1)])^2)
          lod.mse[s.imp] <- mean((unlist(ycomplete)[which(unlist(holdout)==2)] - unlist(y)[which(unlist(holdout)==2)])^2)

          # mean bias
          mar.bias[s.imp] <- mean((unlist(y)[which(unlist(holdout)==1)] - unlist(ycomplete)[which(unlist(holdout)==1)]))
          lod.bias[s.imp] <- mean((unlist(y)[which(unlist(holdout)==2)] - unlist(ycomplete)[which(unlist(holdout)==2)]))
        }
        s.imp <- s.imp+1
      }
      s.save=s.save+1
    }
  }
  
  
  list1 <- list(z.save = z.save, K.save = K.save, mu.save = mu.save,
                ymar = y.mar.save, ylod = y.lod.save,
                hamming = ham, mu.mse = mu.mse, mar.mse = mar.mse, 
                lod.mse = lod.mse, mar.bias = mar.bias, lod.bias = lod.bias,
                mismat = mismat, ycomplete = ycomplete,
                MH.arate = MH.a/(length(al)*sum(K.save)),
                MH.lamrate = MH.lam/(p*sum(K.save)))
  
  class(list1) <- "dpmm"
  return(list1)
  
  
}
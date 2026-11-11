#' Simulate data to fit methods in the R package `psbpHMM' 
#'
#' Simulate multiple time series data from a fixed-state MVN hidden Markov model 
#'
#' @param n number of time series
#' @param t.max length of each time series
#' @param K number of clusters to be shared across all time series 
#' @param p number of exposures 
#' @param trend "shared" (default) for shared temporal trends or "distinct" for distinct temporal trends among n time series
#' @param missingLevel missing data level, data will be removed and split between missing at random and below LOD
#'
#' @return a list with components
#' \itemize{
#'     \item y: scaled simulated data with missing observations set to NA for MAR or -Inf for below LOD 
#'     \item y.complete: scaled complete data 
#'     \item z.true: true hidden state trajectories for each time series
#'     \item K: true number of hiddens states
#'     \item mu.true: scaled true state-specific means
#'     \item lod: limit of detection
#' }     
#' @export
#'

simData <- function(n=20, t.max=288, K = 6, p = 3, 
                   trend = "shared", missingLevel = 0){
  
  if(n < 1) stop("n must be an integer greater than 0")
  if(t.max < 3) stop("t.max must be an integer greater than 2")
  if(K < 1) stop("K must be an integer greater than 0")  
  if(p < 1) stop("p must be an integer greater than 0")
  if(trend != "shared" & trend != "distinct") stop("must specify trend = shared or trend = distinct")
  if(missingLevel < 0 | missingLevel >= 1) stop("missingLevel must be greater than or equal to 0 and less than 1")
  
  
  if(trend == "shared"){
    z.true <- list()
    for(i in 1:n){
      probs <- rdirichlet(1, rep(10, K))
      nums <- rmultinom(1, t.max, prob = probs)
      num1 <- nums[1]
      numFirsthalf <- floor(nums[1]/2)
      numSecondhalf <- num1 - numFirsthalf
      nums[1] <- numFirsthalf
      nums <- c(nums, numSecondhalf)
      states <- c(1:K, 1)
      z <- unlist(lapply(1:(K+1), FUN = function(j) rep(states[j], nums[j])))
      z.true[[i]] <- z 
    }
  }else{
    z.true <- list()
    for(i in 1:n){
      probs <- rdirichlet(1, rep(10, K))
      ord <- sample(1:K, size = K, replace = FALSE) 
      samps <- sample(1:K, t.max, replace = TRUE, prob = probs) 
      reps <- sapply(1:K, FUN = function(i) length(which(samps==i))) 
      states <- c(unlist(lapply(1:K, FUN = function(k) rep(ord[k], reps[k])))); states
      half1 <- ceiling(0.5*length(which(states==states[1])))
      z <- c(states[-(1:half1)], rep(states[1], half1))
      z.true[[i]] <- z 
    }
  }

  # set covariances for each cluster 
  Sigma.true <- list()
  for(k in 1:K){
    L <- diag(p)
    lowerTriangle(L) <- rnorm(p, 0, .5) 
    Sigma.true[[k]] <- 0.1*solve(L)%*%t(solve(L))
  }
  
  # set means for each cluster 
  corMean <- matrix(c(1, 0.7, 0.4, 0.7, 1, -0.2, 0.4, -0.2, 1), ncol = 3)
  mu.true <- rmvnorm(K, mean = rep(0, p), sigma = corMean)
  
  y <- list()
  for(i in 1:n){
    y[[i]] <- matrix(0, t.max, p) # p exposures at t time points 
    for(k in 1:K){
      if(any(z.true[[i]]==k)){
        which.zk <- which(z.true[[i]] == k)
        y[[i]][which.zk,] <- rmvnorm(length(which.zk), mu.true[k,], Sigma.true[[k]])
      }
    }
  }
  
  # standardize data 
  yMat <- numeric()
  for(i in 1:n){
    yMat <- rbind(yMat, y[[i]])
  }
  
  # scale mu for MSE 
  meanY <- apply(yMat, 2, mean)
  sdY <- apply(yMat, 2, sd)
  mu.true.sc <- cbind((mu.true[,1] - meanY[1])/sdY[1],
                      (mu.true[,2] - meanY[2])/sdY[2],
                      (mu.true[,3] - meanY[3])/sdY[3])
  
  ys <- apply(yMat, 2, scale)
  
  # relist 
  splits <- seq(1,n*t.max, t.max)
  y.list <- lapply(1:length(splits), FUN = function(i) ys[splits[i]:(splits[i]+t.max-1),])
  y <- y.list
  y.complete <- y # list of complete data 
  
  
  ymat <- NULL
  for(i in 1:n){
    ymat <- rbind(ymat, y[[i]])
  }

  
  if(missingLevel > 0){
    # split between MAR and below LOD
    marmis = lodmis = missingLevel/2
    # remove MAR data 
    misnum <- ceiling(marmis*t.max*p)
    numChunks <- marmis*200; numChunks
    chunks <- rmultinom(1, size = misnum, prob = rep(1/numChunks, numChunks)) # generate size of chunks 
    for(i in 1:n){
      for(k in 1:numChunks){
        j <- sample(1:p, 1, prob = rep(1/p,p)) # randomly choose an exposure
        repeat{
          timePoint <- sample(2:(t.max-chunks[k,]),1, prob = rep(1/(t.max-chunks[k,]), (t.max-chunks[k,]-1)) ) # choose a time point to start at
          if(!anyNA(y[[i]][,j][timePoint:(timePoint+chunks[k,]-1)])){
            break
          }
        }
        y[[i]][,j][timePoint:(timePoint+chunks[k,]-1)] <- NA
      }
    }
    # remove LOD data 
    lod = list()
    lod.all = apply(ymat, 2, FUN = function(x) quantile(x, lodmis, na.rm = TRUE))
    for(i in 1:n){
      lod[[i]] = lod.all
      for(j in 1:p){
        y[[i]][which(y[[i]][,j] <= lod[[i]][j]),j] <- -Inf
      }
    }
  }else{
    lod = NULL
  }
  
  list1 <- list(y = y, y.complete = y.complete,
                z.true = z.true, K = K, mu.true = mu.true.sc,
                lod = lod)
  return(list1)
}
#' multivariate gamma function
#'
#' @param nu nu (Wishart df)
#' @param p p (dimension of Sigma; number of exposures)
#'
#' @return Gamma_p((nu+1)/2)/Gamma_p(nu/2)


# nu has to be less than 172 for this to work 
mgamma <- function(nu, p){
  num <- prod(sapply(1:p, FUN = function(j) gamma( (nu-j)/2 + 1)))
  den <- prod(sapply(1:p, FUN = function(j) gamma( (nu-j)/2 + 1/2)))
  return(num/den)
} 

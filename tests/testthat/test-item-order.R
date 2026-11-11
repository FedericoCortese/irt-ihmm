test_that("item maps select emissions by actual item and preserve sequence", {
  x <- matrix(c(0,1,1,0,1,1),2,3); y <- matrix(seq(-.3,.7,length.out=6),2,3)
  map <- rbind(c(2,3,1),c(3,1,2))
  data <- irtihmm_data(x,y,item_index=map)
  e <- list(b=c(-1,0,1),a=c(-.2,.3),d=c(.2,-.1),nu=c(.5,1,2),
    mu=c(-.1,.2),eta=c(.1,-.1),log_sigma=log(c(.4,.6)),rho=c(.2,-.3))
  C <- matrix(c(1L,2L,2L,1L,1L,2L),2,3)
  direct <- matrix(0,2,3)
  for (p in 1:2) for (t in 1:3) {
    k<-map[p,t];s<-C[p,t]
    direct[p,t]<-irtihmm_emission(x[p,t],y[p,t],e$b[k]+e$a[s]+e$d[p],
      e$nu[k]+e$mu[s]+e$eta[p],exp(e$log_sigma[s]),e$rho[s])
  }
  expect_equal(.emission_path_loglik(data,e,C),sum(direct))
  array<-.emission_array(data,e)
  for(p in 1:2) for(t in 1:3) expect_equal(array[p,t,C[p,t]],direct[p,t])
  design<-.original_time_design(C,3,2,2,map)
  expect_equal(as.numeric(design%*%unlist(e[c('nu','mu','eta')])),
    as.numeric(.item_predictor(data,e$eta,e$nu)+matrix(e$mu[C],2)))
  means <- vapply(1:3,function(k) mean(y[map==k]),numeric(1))
  expect_equal(.item_average(data,y),means)
  expect_error(irtihmm_data(x,y,item_index=matrix(1,2,3)), 'permutations')
  fit<-fit_irtihmm(data,iterations=6,burnin=3,seed=4)
  expect_equal(fit$data$item_index, data$item_index)
  expect_true(all(is.finite(fit$diagnostics$conditional_loglik)))
})

test_that("mapped Gaussian fits and item reconstructions use common item identity", {
  map <- rbind(c(2,3,1),c(3,1,2))
  data <- irtihmm_data(matrix(c(0,1,1,0,1,1),2,3),
    matrix(seq(-.3,.7,length.out=6),2,3), item_index=map)
  C <- matrix(c(1L,2L,2L,1L,1L,2L),2,3)
  phi <- seq_len(5)/10
  e <- .effect_from_phi(phi,3,2,2)
  expect_equal(as.numeric(.time_design(C,3,2,2,map)%*%phi),
    as.numeric(matrix(e$item[map],2)+matrix(e$state[C],2)+e$person))
  fit <- fit_irtihmm(data, mode='finite',states=C, copula=FALSE,
    iterations=7,burnin=3,seed=8)
  expect_true(all(is.finite(fit$diagnostics$conditional_loglik)))
  one <- fit; one$draws <- fit$draws[1]
  e <- one$draws[[1]]$emission; C <- one$draws[[1]]$states
  predicted <- plogis(matrix(e$b[map],2)+matrix(e$a[C],2)+e$d)
  expected <- vapply(1:3,function(k)mean(predicted[map==k]),numeric(1))
  expect_equal(.fit_reconstruction(one)$accuracy[,'mean'],expected)
})

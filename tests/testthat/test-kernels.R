test_that("compiled likelihoods and per-margin caches match independent formulas", {
  set.seed(12)
  data <- irtihmm_data(matrix(rbinom(24,1,.6),4,6),matrix(rnorm(24),4,6),
    item_index=t(replicate(4,sample.int(6))))
  e <- .new_emission(data,3,.initial_hyper(irtihmm_priors()),TRUE)
  e$rho <- c(0,-.7,.8); C <- matrix(sample.int(3,24,TRUE),4,6)
  reference <- function(e) {
    a <- .item_predictor(data,e$d,e$b)+matrix(e$a[C],4)
    m <- .item_predictor(data,e$eta,e$nu)+matrix(e$mu[C],4)
    sigma <- exp(e$log_sigma[C]);r <- e$rho[C];z <- (data$y-m)/sigma
    cut <- qnorm(plogis(-a)); b <- (cut-r*z)/sqrt(1-r^2)
    sum(dnorm(z,log=TRUE)-log(sigma)+pnorm(ifelse(data$x==1,-b,b),log.p=TRUE))
  }
  expect_equal(.emission_path_loglik(data,e,C),reference(e),tolerance=1e-12)
  for (time in c(FALSE,TRUE)) {
    ctx <- emission_context_cpp(data,e,C,time)
    candidate <- e
    nm <- if(time)c('nu','mu','eta') else c('b','a','d')
    candidate[nm] <- lapply(candidate[nm],function(v)v+rnorm(length(v),sd=.2))
    actual <- emission_margin_loglik_cpp(ctx,candidate[[nm[1]]],candidate[[nm[2]]],candidate[[nm[3]]])
    expect_equal(actual,reference(candidate),tolerance=1e-12)
  }
  array <- .emission_array(data,e)
  expect_equal(sum(array[cbind(rep(1:4,6),rep(1:6,each=4),as.numeric(C))]),reference(e),tolerance=1e-12)
  a <- as.numeric(.item_predictor(data,e$d,e$b)+matrix(e$a[C],4))
  m <- as.numeric(.item_predictor(data,e$eta,e$nu)+matrix(e$mu[C],4))
  ctx <- emission_scale_context_cpp(as.numeric(data$x),as.numeric(data$y),a,m)
  for(r in c(0,-.9,.8)) for(sigma in c(.2,1.2)) {
    z <- (as.numeric(data$y)-m)/sigma
    b <- (qnorm(plogis(-a))-r*z)/sqrt(1-r^2)
    ref <- sum(dnorm(z,log=TRUE)-log(sigma)+pnorm(ifelse(data$x==1,-b,b),log.p=TRUE))
    expect_equal(emission_scale_loglik_cpp(ctx,log(sigma),r),ref,tolerance=1e-12)
  }
  expect_identical(emission_scale_loglik_cpp(ctx,1000,.2),-Inf)
})

test_that("compiled Gaussian sufficient statistics match dense designs", {
  set.seed(5);P<-3;I<-4;K<-2
  C<-matrix(sample.int(K,P*I,TRUE),P,I);map<-t(replicate(P,sample.int(I)))
  y<-matrix(rnorm(P*I),P,I);w<-runif(P*I,.2,2)
  for(mapped in list(NULL,map)) for(reduced in c(FALSE,TRUE)) {
    X<-if(reduced).time_design(C,I,K,P,mapped) else .original_time_design(C,I,K,P,mapped)
    out<-gaussian_time_stats_cpp(C,mapped,y,w,K,reduced)
    expect_equal(out$precision,crossprod(X,X*w),tolerance=1e-13)
    expect_equal(as.numeric(out$information),as.numeric(crossprod(X,as.numeric(y)*w)),tolerance=1e-13)
  }
  for(dims in list(c(1,1,1),c(4,2,3),c(1,5,1))) {
    I<-dims[1];K<-dims[2];P<-dims[3];v<-c(.2,.7,1.3)
    T<-.location_transform(I,K,P);variance<-rep(v,c(I,K,P))
    expect_equal(location_prior(I,K,P,variances=v)$covariance,
      (T*rep(variance,each=nrow(T)))%*%t(T),tolerance=1e-13)
  }
})

test_that("compiled beam slices and batched FFBS preserve events and RNG order", {
  set.seed(4);tr<-.new_transition(3,1,'infinite');z<-matrix(c(-1,1),2,1)
  pr<-.transition_probs(tr,z);C<-matrix(c(1,2,3,1,2,3),2,3)
  set.seed(8);expected<-matrix(0,2,3)
  for(p in 1:2)for(t in 1:3){mass<-if(t==1)pr$initial[p,C[p,t]] else pr$transition[[p]][C[p,t-1],C[p,t]]
    u<-0;while(u==0)u<-runif(1,0,mass);expected[p,t]<-u}
  set.seed(8);expect_identical(beam_slices_cpp(C,pr$initial,pr$transition),expected)
  loge<-array(rnorm(18),c(2,3,3))
  for(slices in list(matrix(numeric(),0,0),expected)) {
    set.seed(13);ref<-matrix(0L,2,3);ll<-0
    for(p in 1:2){o<-forward_backward_cpp(matrix(loge[p,,],3,3),pr$initial[p,],pr$transition[[p]],
      if(length(slices))slices[p,] else numeric(),TRUE);ref[p,]<-o$path;ll<-ll+o$loglik}
    set.seed(13);out<-hmm_paths_cpp(loge,pr$initial,pr$transition,slices)
    expect_identical(out$C,ref);expect_equal(out$loglik,ll)
  }
  bad<-pr$initial;bad[1,1]<-0
  expect_error(beam_slices_cpp(C,bad,pr$transition),'zero transition mass')
})

test_that("one-sided Gaussian sampler matches truncated-normal CDFs including tails", {
  for(mu in c(-20,-3,0,3)) {
    set.seed(31);x<-positive_normal_cpp(mu,12000)
    # Probability integral transform under N(mu,1) conditional on positivity.
    u<-1-exp(pnorm(x,mu,1,lower.tail=FALSE,log.p=TRUE)-pnorm(0,mu,1,lower.tail=FALSE,log.p=TRUE))
    expect_true(all(x>0 & is.finite(x)))
    expect_lt(abs(mean(u)-.5),5/sqrt(12*length(u)))
    expect_lt(abs(mean(u<.1)-.1),5*sqrt(.1*.9/length(u)))
  }
})

test_that("compiled transition update matches its analytical intercept mean", {
  tr<-.new_transition(2,0,'infinite');tr$alpha[,]<-0;tr$alpha0[]<-0
  tr$mean_diagonal<-.4;tr$var_diagonal<-.8;tr$sigma_alpha2<-1.2
  C<-matrix(c(1L,2L,2L,1L,1L,2L),2,3);z<-matrix(numeric(),2,0)
  before<-tr
  set.seed(42)
  samples<-replicate(2500,as.numeric(.update_transition(tr,C,z)$alpha))
  # Every auxiliary w has E(w)=+/-sqrt(2/pi) at zero linear predictor.
  means<-matrix(0,2,2)
  for(r in 1:2)for(k in 1:2){positive<-negative<-0
    for(p in 1:2)for(t in 2:3)if(C[p,t-1]==r && C[p,t]>=k){
      if(C[p,t]==k)positive<-positive+1 else negative<-negative+1}
    v<-if(r==k).8 else 1.2;m<-if(r==k).4 else 0
    means[r,k]<-(m/v+(positive-negative)*sqrt(2/pi))/(1/v+positive+negative)
  }
  se<-apply(samples,1,sd)/sqrt(ncol(samples))
  expect_lt(max(abs(rowMeans(samples)-as.numeric(means))/se),5)
  expect_identical(tr,before)
})

test_that("compiled covariate and tuned hyperprior updates have the right moments", {
  set.seed(8);P<-3;I<-4;K<-3;L<-3
  z<-cbind(c(-1,0,1),c(.4,-.8,.7));C<-matrix(c(1,2,3,2,3,1,3,1,2,1,2,3),P,I)
  tr<-.new_transition(K,2,'infinite');tr$var_diagonal<-.7;tr$sigma_alpha2<-1.1
  prior<-list(mean=1.2,sd=.25)
  xb<-z%*%tr$beta;wmean<-array(0,c(P,I,L));counts<-sums<-matrix(0,K,L);n0<-s0<-numeric(L)
  for(p in 1:P)for(t in 1:I)for(k in seq_len(C[p,t])) {
    a<-if(t==1)tr$alpha0[k] else tr$alpha[C[p,t-1],k]
    m<-a+xb[p,k];positive<-C[p,t]==k
    ew<-m+if(positive)exp(dnorm(m,log=TRUE)-pnorm(m,log.p=TRUE)) else
      -exp(dnorm(m,log=TRUE)-pnorm(-m,log.p=TRUE))
    wmean[p,t,k]<-ew
    if(t==1){n0[k]<-n0[k]+1;s0[k]<-s0[k]+ew-xb[p,k]} else {
      r<-C[p,t-1];counts[r,k]<-counts[r,k]+1;sums[r,k]<-sums[r,k]+ew-xb[p,k]}
  }
  amean<-matrix(0,K,L)
  for(r in 1:K)for(k in 1:L){v<-if(r==k)tr$var_diagonal else tr$sigma_alpha2
    m<-if(r==k)tr$mean_diagonal else 0
    amean[r,k]<-(m/v+sums[r,k])/(1/v+counts[r,k])}
  a0mean<-s0/(1/tr$sigma_alpha2+n0)
  target<-matrix(0,2,L)
  for(k in 1:L){Q<-diag(2);h<-numeric(2)
    for(p in 1:P)for(t in 1:I)if(C[p,t]>=k){
      a<-if(t==1)a0mean[k] else amean[C[p,t-1],k]
      Q<-Q+tcrossprod(z[p,]);h<-h+z[p,]*(wmean[p,t,k]-a)}
    target[,k]<-solve(Q,h)}
  set.seed(55)
  samples<-replicate(3000,{
    out<-.update_transition(tr,C,z,prior)
    v<-1/(1/prior$sd^2+L/tr$var_diagonal)
    m<-v*(prior$mean/prior$sd^2+sum(diag(out$alpha))/tr$var_diagonal)
    residual<-(out$mean_diagonal-m)/sqrt(v)
    rate<-1+sum((diag(out$alpha)-out$mean_diagonal)^2)/2
    uniform<-pgamma(1/out$var_diagonal,shape=1+L/2,rate=rate)
    c(as.numeric(out$beta),residual,uniform)
  })
  se<-apply(samples[1:6,],1,sd)/sqrt(ncol(samples))
  expect_lt(max(abs(rowMeans(samples[1:6,])-as.numeric(target))/se),5)
  expect_lt(abs(mean(samples[7,])),5/sqrt(ncol(samples)))
  expect_lt(abs(mean(samples[8,])-.5),5/sqrt(12*ncol(samples)))
})

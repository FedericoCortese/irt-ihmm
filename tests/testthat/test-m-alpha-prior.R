test_that("m_alpha prior defaults and legacy objects preserve N(0,1)", {
  prior <- irtihmm_priors()
  expect_equal(prior$m_alpha, list(mean = 0, sd = 1))
  old <- prior; old$m_alpha <- NULL
  expect_identical(.validate_priors(old), prior)
  set.seed(31)
  expected <- .new_transition(3, 2, "infinite")
  set.seed(31)
  expect_identical(.new_transition(3, 2, "infinite", prior$m_alpha), expected)
  d <- c(-.4, 1.5, .7); v <- .8
  variance <- 1/(1+length(d)/v)
  set.seed(20); expected <- rnorm(1, variance*sum(d)/v, sqrt(variance))
  set.seed(20); expect_identical(.draw_m_alpha(d,v,prior$m_alpha),expected)
})

test_that("custom m_alpha prior enters the exact Gaussian full conditional", {
  prior <- irtihmm_priors(m_alpha_mean=1.5,m_alpha_sd=.2)$m_alpha
  d <- c(-.3,.5,1.2); v <- .8
  variance <- 1/(1/.2^2+length(d)/v)
  mean <- variance*(1.5/.2^2+sum(d)/v)
  set.seed(19); expected <- rnorm(1,mean,sqrt(variance))
  set.seed(19); expect_identical(.draw_m_alpha(d,v,prior),expected)
  # Empty diagonal (single-state finite model) must draw from the tuned prior.
  tr <- .new_transition(1,0,'finite',prior)
  set.seed(41); .rinvgamma(1,1); expected <- rnorm(1,1.5,.2)
  set.seed(41)
  out <- .update_transition(tr,matrix(1L,2,3),matrix(numeric(),2,0),prior)
  expect_identical(out$mean_diagonal,expected)
  expect_equal(tr$mean_diagonal,1.5)
  # Initial/off-diagonal priors do not acquire the diagonal mean.
  set.seed(5); a <- .new_transition(3,1,'infinite')
  set.seed(5); b <- .new_transition(3,1,'infinite',prior)
  expect_equal(diag(b$alpha)-diag(a$alpha),rep(1.5,3))
  expect_equal(b$alpha[row(b$alpha)!=col(b$alpha)],a$alpha[row(a$alpha)!=col(a$alpha)])
  expect_identical(b$alpha0,a$alpha0)
})

test_that("invalid m_alpha priors are rejected", {
  for (sd in list(0,-1,NA_real_,Inf,c(1,2),'small'))
    expect_error(irtihmm_priors(m_alpha_sd=sd),'m_alpha')
  for (mean in list(NA_real_,Inf,c(0,1),'zero'))
    expect_error(irtihmm_priors(m_alpha_mean=mean),'m_alpha')
  p <- irtihmm_priors(); p$m_alpha <- list(mean=0,variance=1)
  expect_error(.validate_priors(p),'m_alpha')
})

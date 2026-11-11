test_that("fixed blocks remain exact with and without Gustafson and copulas", {
  sim <- simulate_irtihmm(3, 4, 2, seed = 9)
  for (gustafson in c(TRUE, FALSE)) for (copula in c(TRUE, FALSE)) {
    fit <- fit_irtihmm(sim$data, iterations = 10, burnin = 4, seed = 7,
      gustafson = gustafson, copula = copula,
      fixed_params = list(eta = c(-.1, 0, .2), b = .3, rho = 0, log_sigma = log(.4)))
    for (draw in fit$draws) {
      expect_identical(draw$emission$eta, c(-.1, 0, .2))
      expect_equal(draw$emission$b, rep(.3, 4))
      expect_true(all(draw$emission$rho == 0))
      expect_true(all(draw$emission$log_sigma == log(.4)))
      expect_null(draw$locations$accuracy)
      expect_null(draw$locations$time)
      expect_equal(draw$hyper$time$variances[3], fit$draws[[1]]$hyper$time$variances[3])
    }
    expect_gt(sd(vapply(fit$draws, function(d) d$emission$nu[1], numeric(1))), 0)
  }
  fit <- fit_irtihmm(sim$data, iterations = 8, burnin = 4, seed = 7,
    fixed_params = list(eta = 0))
  expect_true(all(vapply(fit$draws, function(d)
    all(d$emission$eta == 0) && !is.null(d$locations$accuracy) && is.null(d$locations$time), logical(1))))
})

test_that("fixed state effects are enforced during beam expansion", {
  sim <- simulate_irtihmm(2, 3, 2, seed = 9)
  priors <- irtihmm_priors(); hyper <- .initial_hyper(priors)
  hyper$copula <- TRUE; hyper$fixed_params <- list(a = .2, mu = -.1, rho = .3)
  emission <- .apply_fixed_params(.new_emission(sim$data, 2, hyper, TRUE), hyper$fixed_params)
  expanded <- .expand_emission(emission, 7, hyper)
  expect_equal(expanded$a, rep(.2, 7))
  expect_equal(expanded$mu, rep(-.1, 7))
  expect_equal(expanded$rho, rep(.3, 7))
  fit <- fit_irtihmm(sim$data, mode = "finite", K = 2, iterations = 6, burnin = 3,
    fixed_params = list(a = c(-.2, .2), mu = c(-.4, .4), eta = 0), seed = 4)
  expect_true(all(vapply(fit$draws, function(d) identical(d$emission$a, c(-.2, .2)), logical(1))))
})

test_that("invalid fixed blocks fail before sampling", {
  sim <- simulate_irtihmm(2, 3, 2, seed = 9)
  for (fixed in list(list(unknown = 0), list(eta = NA), list(eta = numeric()),
                     list(eta = 1:3), list(rho = 1), list(log_sigma = 1000),
                     list(a = c(-1, 1)), list(eta = 0, eta = 1), list(0)))
    expect_error(fit_irtihmm(sim$data, fixed_params = fixed))
  expect_error(fit_irtihmm(sim$data, copula = FALSE, fixed_params = list(rho = .2)), "rho")
})

test_that("Gaussian constrained updates include fixed nonzero offsets", {
  data <- irtihmm_data(matrix(0, 2, 2), matrix(c(-.2, .4, .7, 1.1), 2))
  C <- matrix(1L, 2, 2)
  e <- list(b = c(0,0), a = 0, d = c(0,0), nu = c(.2,.3), mu = 0,
            eta = c(-.1,.4), log_sigma = log(.5), rho = 0)
  h <- list(time = list(mean = .4, variances = c(.7,.9,.4)))
  design <- .original_time_design(C, 2, 1, 2)
  offset <- as.numeric(design[, -3] %*% c(e$nu, e$eta))
  variance <- 1 / (1/.9 + 4/.5^2)
  mean <- variance * sum((as.numeric(data$y) - offset)/.5^2)
  set.seed(4)
  samples <- replicate(2000, .update_original_margin(data, e, h, C, "time", FALSE,
    list(nu = e$nu, eta = e$eta))$effects$state)
  expect_lt(abs(mean(samples) - mean), 4 * sqrt(variance/2000))
  expect_lt(abs(var(samples) - variance), 4 * variance * sqrt(2/1999))
})

test_that("entire fixed margins skip sampling without shifting their values", {
  sim <- simulate_irtihmm(2, 3, 2, seed = 9)
  fixed <- list(b = .1, a = -.2, d = .3, nu = 2, mu = -.4, eta = 0,
                log_sigma = log(.5), rho = 0)
  fit <- fit_irtihmm(sim$data, mode = "finite", K = 2, iterations = 5,
    burnin = 2, seed = 5, fixed_params = fixed)
  for (draw in fit$draws) for (name in names(fixed))
    expect_equal(draw$emission[[name]], rep(fixed[[name]], length(draw$emission[[name]])))
  expect_true(all(fit$diagnostics$accuracy_evaluations == 0))
  expect_true(all(fit$diagnostics$time_evaluations == 0))
})

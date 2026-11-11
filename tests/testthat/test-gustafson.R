test_that("the default Gustafson sampler matches the explicit TRUE setting", {
  sim <- simulate_irtihmm(3, 4, 2, seed = 9)
  default <- fit_irtihmm(sim$data, iterations = 6, burnin = 3, seed = 13)
  explicit <- fit_irtihmm(sim$data, iterations = 6, burnin = 3, seed = 13, gustafson = TRUE)
  expect_identical(default$draws, explicit$draws)
  expect_identical(default$diagnostics, explicit$diagnostics)
  expect_true(explicit$gustafson)
})

test_that("the direct Gaussian time update has the original-prior conditional moments", {
  data <- irtihmm_data(matrix(0, 2, 2), matrix(c(-.2, .4, .7, 1.1), 2, 2))
  C <- matrix(c(1L, 2L, 2L, 1L), 2, 2)
  emission <- list(b = c(0, 0), a = c(0, 0), d = c(0, 0), nu = c(0, 0),
    mu = c(0, 0), eta = c(0, 0), log_sigma = log(c(.3, .5)), rho = c(0, 0))
  hyper <- list(time = list(mean = .4, variances = c(.7, .9, .4)))
  design <- .original_time_design(C, 2, 2, 2)
  variance <- rep(hyper$time$variances, each = 2)
  mean <- c(.4, .4, 0, 0, 0, 0)
  weight <- as.numeric(exp(-2 * emission$log_sigma[C]))
  precision <- diag(1 / variance) + crossprod(design, design * weight)
  covariance <- solve(precision)
  target <- as.numeric(covariance %*% (mean / variance + crossprod(design, as.numeric(data$y) * weight)))
  set.seed(67)
  values <- replicate(2500, {
    out <- .update_original_margin(data, emission, hyper, C, "time", FALSE)
    unlist(out$effects, use.names = FALSE)
  })
  # Compare Monte Carlo errors with their normal/Wishart standard errors.
  mean_se <- sqrt(diag(covariance) / ncol(values))
  covariance_se <- sqrt((outer(diag(covariance), diag(covariance)) + covariance^2) /
                          (ncol(values) - 1))
  expect_lt(max(abs(rowMeans(values) - target) / mean_se), 4)
  expect_lt(max(abs(cov(t(values)) - covariance) / covariance_se), 4)
})

test_that("direct-effect sampling supports copulas, Gaussian updates and adaptive states", {
  sim <- simulate_irtihmm(3, 5, 2, z = matrix(c(-1, 0, 1), 3, 1), seed = 9)
  for (copula in c(FALSE, TRUE)) for (mode in c("finite", "infinite")) {
    fit <- fit_irtihmm(sim$data, mode = mode, copula = copula, gustafson = FALSE,
      iterations = 8, burnin = 4, seed = 12)
    expect_false(fit$gustafson)
    expect_equal(fit$priors, irtihmm_priors())
    expect_true(all(is.finite(fit$diagnostics$conditional_loglik)))
    for (draw in fit$draws) {
      expect_length(draw$locations, 0)
      expect_equal(draw$reduced$accuracy,
        reduce_effects(draw$emission$b, draw$emission$a, draw$emission$d))
      expect_equal(draw$reduced$time,
        reduce_effects(draw$emission$nu, draw$emission$mu, draw$emission$eta))
      expect_true(all(abs(draw$emission$rho) < 1))
      if (!copula) expect_equal(draw$emission$rho, rep(0, length(draw$emission$rho)))
    }
  }
  expect_error(fit_irtihmm(sim$data, gustafson = NA), "gustafson")
  expect_error(fit_irtihmm(sim$data, gustafson = c(TRUE, FALSE)), "gustafson")
})

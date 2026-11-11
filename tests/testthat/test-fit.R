test_that("known-state, finite and infinite joint samplers produce finite draws", {
  sim <- simulate_irtihmm(3, 5, 2, z = matrix(c(-1, 0, 1), 3, 1), seed = 9)
  configurations <- list(list(mode = "finite", copula = FALSE, states = sim$states),
                         list(mode = "finite", copula = TRUE),
                         list(mode = "infinite", copula = TRUE))
  for (configuration in configurations) {
    fit <- do.call(fit_irtihmm, c(list(data = sim$data, K = 2, iterations = 8,
      burnin = 4, seed = 12), configuration))
    expect_s3_class(fit, "irtihmm_fit")
    expect_length(fit$draws, 4)
    expect_true(all(is.finite(fit$diagnostics$conditional_loglik)))
    for (draw in fit$draws) {
      expect_equal(reduce_effects(draw$emission$b, draw$emission$a, draw$emission$d), draw$reduced$accuracy)
      expect_equal(reduce_effects(draw$emission$nu, draw$emission$mu, draw$emission$eta), draw$reduced$time)
      expect_true(all(abs(draw$emission$rho) < 1))
      if (!fit$copula) expect_true(all(draw$emission$rho == 0))
    }
    if (fit$mode == "infinite") {
      expect_true(all(fit$diagnostics$max_residual < fit$diagnostics$min_slice))
      for (draw in fit$draws)
        expect_equal(draw$transition$K, max(draw$states))
    }
  }
})

test_that("the single-person single-item finite model is supported", {
  data <- irtihmm_data(matrix(1, 1, 1), matrix(0, 1, 1))
  fit <- fit_irtihmm(data, mode = "finite", K = 1, copula = FALSE,
                     iterations = 3, burnin = 1, seed = 7)
  expect_length(fit$draws, 2)
  expect_equal(fit$diagnostics$occupied, rep(1, 3))
})

test_that("unoccupied gaps retain their stick indices", {
  sim <- simulate_irtihmm(2, 3, 3, seed = 15)
  C <- matrix(3L, 2, 3)
  fit <- fit_irtihmm(sim$data, mode = "infinite", K = 3, states = C,
    iterations = 3, burnin = 1, seed = 10)
  for (draw in fit$draws) {
    expect_equal(draw$states, C)
    expect_equal(length(draw$emission$a), 3)
    expect_equal(draw$transition$K, 3)
  }
})

test_that("complete common-order data and prior contracts are enforced", {
  expect_error(irtihmm_data(matrix(NA_real_, 2, 2), matrix(0, 2, 2)), "binary")
  expect_error(irtihmm_data(matrix(0, 2, 2), matrix(0, 2, 2), matrix(1, 2, 1)), "constant")
  sim <- simulate_irtihmm(2, 2, 1, seed = 3)
  expect_error(fit_irtihmm(sim$data, iterations = 2, burnin = 2), "retained")
  expect_error(fit_irtihmm(sim$data, priors = list()), "components")
})

test_that("explicit simulation covariates change the generating transition sticks", {
  z <- matrix(c(0, 1), ncol = 1)
  sim <- simulate_irtihmm(P = 2, I = 4, K = 2, z = z,
                          beta = matrix(-1.2, 1, 1), seed = 19)
  expect_equal(as.numeric(sim$transition$beta), -1.2)
  expect_equal(sim$probabilities$initial[, 1], pnorm(-.5 - 1.2 * z[, 1]))
  expect_gt(sim$probabilities$initial[1, 1], sim$probabilities$initial[2, 1])
  expect_error(simulate_irtihmm(2, 4, 2, z = z, beta = matrix(1, 2, 1)), "beta must")
})

test_that("a reference transition matrix is recovered exactly at zero covariates", {
  z <- cbind(group = c(0, 1), score = c(-1, 1))
  beta <- cbind(c(-1.2, .8), c(.8, -.6))
  target <- matrix(.15, 3, 3); diag(target) <- .7
  sim <- simulate_irtihmm(2, 5, 3, z = z, beta = beta, seed = 19,
                          transition_matrix = target)
  reference <- .transition_probs(sim$transition, matrix(0, 1, 2))$transition[[1]]
  expect_equal(reference, target, tolerance = 1e-12)
  expect_equal(sim$transition$alpha[, 1], qnorm(target[, 1]))
  expect_equal(sim$transition$alpha[, 2], qnorm(target[, 2] / (1 - target[, 1])))
  expect_false(isTRUE(all.equal(sim$probabilities$transition[[1]], target)))
  default <- simulate_irtihmm(2, 5, 3, z = z, beta = beta, seed = 19)
  expect_equal(sim$emission, default$emission)
  expect_equal(sim$probabilities$initial, default$probabilities$initial)
  single <- simulate_irtihmm(2, 3, 1, transition_matrix = matrix(1, 1, 1), seed = 9)
  expect_equal(single$probabilities$transition[[1]], matrix(1, 1, 1))
  expect_error(simulate_irtihmm(2, 3, 3, transition_matrix = diag(3)), "positive K by K")
  expect_error(simulate_irtihmm(2, 3, 3, transition_matrix = matrix(.2, 3, 3)), "unit row sums")
  expect_error(simulate_irtihmm(2, 3, 3, transition_matrix = matrix(1, 2, 2)), "K by K")
})

test_that("reconstruction plots preserve label invariance", {
  sim <- simulate_irtihmm(3, 4, 2, seed = 9)
  fit <- fit_irtihmm(sim$data, mode = "finite", K = 2,
                      iterations = 6, burnin = 3, seed = 4)
  first <- .fit_reconstruction(fit)
  permuted <- fit
  for (j in seq_along(permuted$draws)) {
    d <- permuted$draws[[j]]
    d$states <- matrix(3L - d$states, nrow(d$states))
    for (name in c("a", "mu", "log_sigma", "rho"))
      d$emission[[name]] <- rev(d$emission[[name]])
    permuted$draws[[j]] <- d
  }
  expect_equal(.fit_reconstruction(permuted), first)
  file <- tempfile(fileext = ".png")
  on.exit(unlink(file), add = TRUE)
  plot_simulation(sim, file)
  expect_gt(file.info(file)$size, 0)
  actual <- plot_fit(fit, sim, file)
  expect_equal(actual, first)
  expect_gt(file.info(file)$size, 0)
  plot_fit(fit, file = file)
  expect_error(plot_fit(fit, simulate_irtihmm(3, 4, 2, seed = 10), file), "same data")
})

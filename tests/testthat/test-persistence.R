test_that("omega preserves initial probabilities and scales the entire infinite tail", {
  z <- matrix(c(-1, 1), 2, 1)
  for (mode in c("finite", "infinite")) {
    set.seed(4); tr <- .new_transition(3, 1, mode)
    base <- .transition_probs(tr, z)
    tr$omega <- .65
    # Armadillo stores updated intercept vectors with a column dimension.
    tr$alpha0 <- matrix(tr$alpha0, ncol = 1L)
    actual <- .transition_probs(tr, z)
    expect_identical(actual$initial, base$initial)
    expect_equal(actual$residual[, 1], base$residual[, 1])
    expect_equal(actual$residual[, -1], .35 * base$residual[, -1])
    for (p in 1:2) {
      expect_equal(actual$transition[[p]], .35 * base$transition[[p]] + .65 * diag(3))
      expect_equal(rowSums(actual$transition[[p]]) + actual$residual[p, -1], rep(1, 3))
      expect_equal(.covariate_state_probs(tr, z[p, , drop = FALSE], c(3, 1)),
        rbind(actual$initial[p, c(3, 1)], actual$transition[[p]][c(3, 1), c(3, 1)]))
    }
    tr$omega <- 0
    expect_identical(.transition_probs(tr, z), base)
  }
})

test_that("forced self-transitions do not enter any probit sufficient statistics", {
  set.seed(8); tr <- .new_transition(2, 2, "finite")
  z <- matrix(c(-1, 1, .2, -.2), 2, 2)
  C <- rbind(rep(1, 5), rep(2, 5))
  mask <- matrix(FALSE, 2, 5); mask[, 1] <- TRUE
  set.seed(42)
  actual <- transition_update_cpp(tr, C, z, 0, 1, mask)
  set.seed(42)
  reference <- transition_update_cpp(tr, C[, 1, drop = FALSE], z, 0, 1)
  expect_identical(actual, reference)
  set.seed(14); actual <- transition_update_cpp(tr, C, z, 0, 1)
  set.seed(14); reference <- transition_update_cpp(tr, C, z, 0, 1, matrix(TRUE, 2, 5))
  expect_identical(actual, reference)
  C[1, 2] <- 2
  expect_error(transition_update_cpp(tr, C, z, 0, 1, mask), "Only self-transitions")
})

test_that("auxiliary omega Gibbs updates match a directly integrated posterior", {
  tr <- .new_transition(2, 0, "finite")
  tr$alpha[,] <- qnorm(c(.4, .3)); tr$omega <- .5
  C <- matrix(c(1, 1, 2, 2, 2, 1), 1)
  z <- matrix(numeric(), 1, 0); prior <- c(2, 3)
  # Three self-transitions (.4, .7, .7), two changes (.6, .3).
  density <- function(w) dbeta(w, 2, 3) * (w + (1-w)*.4) *
    (w + (1-w)*.7)^2 * (1-w)^2 * .6 * .3
  norm <- integrate(density, 0, 1)$value
  expected <- integrate(function(w) w*density(w), 0, 1)$value / norm
  set.seed(41)
  samples <- numeric(6000)
  for (j in seq_along(samples)) {
    aux <- .draw_persistence(tr, C, z, prior)
    stopifnot(aux$active[1, 1], aux$active[1, 3], aux$active[1, 6])
    tr$omega <- aux$omega
    samples[j] <- aux$omega
  }
  expect_equal(mean(samples[-(1:1000)]), expected, tolerance = .015)
  # With no transitions the update is exactly the prior draw.
  set.seed(9); actual <- .draw_persistence(tr, C[, 1, drop = FALSE], z, prior)$omega
  set.seed(9); expect_identical(actual, rbeta(1, 2, 3))
})

test_that("disabled persistence preserves seeded fits and enabled fits cover both modes", {
  sim <- simulate_irtihmm(3, 5, 2, seed = 5)
  old <- fit_irtihmm(sim$data, iterations = 8, burnin = 4, seed = 8)
  off <- fit_irtihmm(sim$data, iterations = 8, burnin = 4, seed = 8, persistence = FALSE)
  expect_identical(old$draws, off$draws)
  expect_identical(old$diagnostics, off$diagnostics)
  for (mode in c("finite", "infinite")) {
    fit <- fit_irtihmm(sim$data, mode = mode, iterations = 15, burnin = 5,
      seed = 8, persistence = TRUE, omega_prior = c(2, 1))
    w <- vapply(fit$draws, function(d) d$transition$omega, numeric(1))
    expect_true(all(w > 0 & w < 1))
    expect_true(any(vapply(.parameter_trace_series(fit)$series,
      function(s) s$parameter == "omega", logical(1))))
    if (mode == "infinite") expect_true(all(fit$diagnostics$max_residual < fit$diagnostics$min_slice))
    else expect_true(all(is.finite(fit$diagnostics$finite_marginal_loglik)))
  }
  expect_error(fit_irtihmm(sim$data, persistence = NA), "persistence")
  expect_error(fit_irtihmm(sim$data, omega_prior = c(0, 1)), "omega_prior")
})

test_that("omega beam trajectories agree with the enumerated path posterior", {
  sim <- simulate_irtihmm(1, 3, 2, seed = 17)
  tr <- .new_transition(2, 0, "infinite")
  tr$alpha[, 1] <- qnorm(c(.4, .3)); tr$alpha[, 2] <- 9
  tr$alpha0 <- c(qnorm(.6), 9); tr$omega <- .6
  probs <- .transition_probs(tr, sim$data$z)
  le <- matrix(.emission_array(sim$data, sim$emission), 3, 2)
  paths <- as.matrix(expand.grid(rep(list(1:2), 3)))
  mass <- apply(paths, 1, function(path) probs$initial[1, path[1]] *
    probs$transition[[1]][path[1], path[2]] *
    probs$transition[[1]][path[2], path[3]] * prod(exp(le[cbind(1:3, path)])))
  expected <- sum(mass[paths[, 1] == paths[, 3]]) / sum(mass)
  C <- matrix(1, 1, 3); samples <- logical(2500)
  hyper <- .initial_hyper(irtihmm_priors())
  set.seed(39)
  for (j in seq_along(samples)) {
    draw <- .beam_state_draw(sim$data, sim$emission, tr, C, hyper, 20)
    C <- draw$C
    samples[j] <- C[1, 1] == C[1, 3]
  }
  expect_equal(mean(samples[-(1:500)]), expected, tolerance = .04)
})

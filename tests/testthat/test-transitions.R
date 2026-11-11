test_that("Hoskovec probabilities match the PSBP formula with correct orientation", {
  set.seed(16)
  z <- matrix(c(-1, 1, .3, -.2), 2, 2)
  for (mode in c("finite", "infinite")) {
    tr <- .new_transition(3, 2, mode)
    probs <- .transition_probs(tr, z)
    for (p in 1:2) for (r in 0:3) {
      v <- pnorm((if (r == 0) tr$alpha0 else tr$alpha[r, ]) + z[p, ] %*% tr$beta)
      weights <- v * c(1, head(cumprod(1 - v), -1))
      tail <- prod(1 - v)
      if (mode == "finite") weights <- c(weights, tail)
      got <- if (r == 0) probs$initial[p, ] else probs$transition[[p]][r, ]
      expect_equal(as.numeric(got), as.numeric(weights), tolerance = 1e-12)
      expect_equal(sum(got) + if (mode == "infinite") probs$residual[p, r + 1] else 0, 1)
    }
  }
})

test_that("zero covariates and a single item work through upstream routines", {
  set.seed(3)
  z <- matrix(numeric(), 2, 0)
  tr <- .new_transition(2, 0, "infinite")
  p <- .transition_probs(tr, z)
  expect_equal(dim(p$initial), c(2L, 2L))
  updated <- .update_transition(tr, matrix(c(1L, 2L), 2, 1), z)
  expect_true(all(is.finite(unlist(updated[setdiff(names(updated), "mode")]))))
})

test_that("beam bound aborts rather than accepting a truncated tail", {
  sim <- simulate_irtihmm(2, 3, 1, copula = FALSE, seed = 2)
  expect_error(fit_irtihmm(sim$data, K = 1, max_states = 1,
                         iterations = 2, burnin = 0, seed = 4), "max_states")
})

test_that("the reused Gaussian intercept kernel has the correct conditional moments", {
  set.seed(62)
  K <- 2L; L <- 1L; P <- 2L; I <- 3L
  paths <- list(c(1, 2, 2), c(2, 1, 2))
  X <- list(matrix(c(-1, .4, 1), 3, 1), matrix(c(.7, -.3, 1.2), 3, 1))
  beta <- list(matrix(.4, 1, 1))
  w <- list(lapply(c(.3, -.7, .4), matrix, ncol = 1),
            lapply(c(-.2, .6, -.1), matrix, ncol = 1))
  expected_mean <- expected_var <- numeric(K)
  for (r in 1:K) {
    values <- numeric()
    for (p in 1:P) for (i in 2:I) if (paths[[p]][i - 1] == r)
      values <- c(values, w[[p]][[i]][1] - .4 * X[[p]][i, 1])
    precision <- if (r == 1) 1 / .7 else 1 / 1.2
    prior_mean <- if (r == 1) .3 else 0
    expected_var[r] <- 1 / (precision + length(values))
    expected_mean[r] <- expected_var[r] * (precision * prior_mean + sum(values))
  }
  sampled <- replicate(4000, unlist(hoskovec_alpha(K, L, P, I,
    paths, 1 / .7, 1 / 1.2, w, X, beta, .3, 0)))
  expect_equal(rowMeans(sampled), expected_mean, tolerance = .04)
  expect_equal(apply(sampled, 1, var), expected_var, tolerance = .03)
})

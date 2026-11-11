test_that("forward likelihood and joint backward draws agree with enumeration", {
  le <- log(matrix(c(.9, .3, .7, .1, .7, .3), 3, 2))
  init <- c(.65, .35); A <- matrix(c(.8, .2, .25, .75), 2, byrow = TRUE)
  paths <- as.matrix(expand.grid(rep(list(1:2), 3)))
  mass <- apply(paths, 1, function(path) init[path[1]] * A[path[1], path[2]] *
    A[path[2], path[3]] * prod(exp(le[cbind(1:3, path)])))
  expect_equal(hmm_loglik(le, init, A), log(sum(mass)), tolerance = 1e-12)
  set.seed(5)
  sampled <- replicate(4000, forward_backward_cpp(le, init, A, numeric(), TRUE)$path)
  first_last <- mean(sampled[1, ] == sampled[3, ])
  expect_equal(first_last, sum(mass[paths[, 1] == paths[, 3]]) / sum(mass), tolerance = .03)
})

test_that("slice-restricted trajectories have the exact augmented path law", {
  le <- log(matrix(c(.8, .3, .6, .2, .7, .4), 3, 2))
  init <- c(.6, .4); A <- matrix(c(.8, .2, .3, .7), 2, byrow = TRUE)
  u <- c(.25, .25, .15)
  paths <- as.matrix(expand.grid(rep(list(1:2), 3)))
  mass <- apply(paths, 1, function(path) as.numeric(init[path[1]] > u[1]) *
    as.numeric(A[path[1], path[2]] > u[2]) * as.numeric(A[path[2], path[3]] > u[3]) *
    prod(exp(le[cbind(1:3, path)])))
  set.seed(8)
  sampled <- replicate(4000, forward_backward_cpp(le, init, A, u, TRUE)$path)
  p <- mean(sampled[1, ] == 1 & sampled[3, ] == 2)
  expect_equal(p, sum(mass[paths[, 1] == 1 & paths[, 3] == 2]) / sum(mass), tolerance = .03)
})

test_that("single-state and impossible HMM cases are handled", {
  expect_equal(hmm_loglik(matrix(c(-1, -2), 2, 1), 1, matrix(1, 1, 1)), -3)
  expect_equal(hmm_loglik(matrix(-Inf, 2, 1), 1, matrix(1, 1, 1)), -Inf)
  expect_error(hmm_loglik(matrix(0, 2, 2), c(.2, .3), diag(2)), "normalized")
})

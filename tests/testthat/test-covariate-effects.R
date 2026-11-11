covariate_fixture <- function() {
  data <- irtihmm_data(matrix(0, 2, 2), matrix(0, 2, 2),
                       cbind(group = c(0, 1), score = c(-1, 1)))
  draws <- lapply(1:3, function(j) {
    K <- if (j == 2) 4L else 3L
    tr <- list(K = K, mode = "infinite", alpha = matrix(-.3, K, K),
      alpha0 = rep(-.5, K), beta = rbind(seq(-.5, .5, length.out = K), rep(.7, K)))
    C <- if (j == 2) c(1L, 1L, 4L, 4L) else c(3L, 3L, 1L, 1L)
    if (j == 3) C[1] <- 1L
    list(states = matrix(C, 2, 2, byrow = TRUE), transition = tr,
      emission = list(a = rep(0, K), mu = rep(0, K), log_sigma = rep(0, K), rho = rep(0, K)))
  })
  fit <- structure(list(data = data, draws = draws), class = "irtihmm_fit")
  post <- structure(list(chips = list(subpartition = matrix(c(1L, 1L, 2L, 2L), 2, 2, byrow = TRUE),
                                      probability = 2 / 3)), class = "irtihmm_postprocess")
  list(fit = fit, post = post)
}

test_that("selected probabilities preserve preceding sticks and the original tail", {
  set.seed(7)
  for (mode in c("finite", "infinite")) for (K in c(1L, 4L)) {
    tr <- .new_transition(K, 2, mode)
    z <- c(.3, -.2)
    states <- unique(c(K, 1L))
    selected <- .covariate_state_probs(tr, z, states)
    full <- .transition_probs(tr, matrix(z, 1, 2))
    expected <- rbind(full$initial[1, states], full$transition[[1]][states, states, drop = FALSE])
    expect_equal(selected, expected, tolerance = 1e-12)
    if (mode == "infinite") expect_true(all(rowSums(selected) < 1))
    if (mode == "finite" && K == 1) expect_equal(selected, matrix(1, 2, 1))
  }
})

test_that("coefficients retain stick order and probability changes follow CHIPS matching", {
  f <- covariate_fixture()
  original <- f$fit
  x <- summarize_covariate_effects(f$fit, f$post)
  expect_identical(f$fit, original)
  expect_identical(x$draw_indices, c(1L, 2L))
  expect_equal(x$common_sticks, 3)
  expect_equal(unname(x$baseline), c(0, 0))
  expect_equal(unname(x$change), c(1, 1))
  expect_equal(x$beta_draws[1, "group", ], c(`1` = -.5, `2` = 0, `3` = .5))
  expect_equal(x$beta_draws[2, "group", ],
    setNames(f$fit$draws[[2]]$transition$beta[1, 1:3], 1:3))
  tr <- f$fit$draws[[1]]$transition
  # Core states 1,2 correspond to ordered states 3,1 in the first draw.
  expected <- .covariate_state_probs(tr, c(1, 0), c(3L, 1L)) -
    .covariate_state_probs(tr, c(0, 0), c(3L, 1L))
  expect_equal(unname(x$probability_draws[1, , , "group"]), expected)
  expect_equal(dim(x$probability_draws), c(2, 3, 2, 2))
  expect_equal(x$core_probability, 2 / 3)
  expect_true(all(abs(x$probability_draws) <= 1))
  changed <- summarize_covariate_effects(f$fit, f$post, baseline = c(0, .5), change = c(1, 2))
  expect_false(isTRUE(all.equal(changed$probability_draws, x$probability_draws)))
})

test_that("covariate plots and invalid inputs are handled", {
  f <- covariate_fixture()
  x <- summarize_covariate_effects(f$fit, f$post)
  file <- tempfile(fileext = ".png")
  on.exit(unlink(file), add = TRUE)
  expect_identical(plot_covariate_effects(x, "coefficients", file), x)
  expect_gt(file.info(file)$size, 0)
  plot_covariate_effects(x, "probabilities", file)
  expect_gt(file.info(file)$size, 0)
  expect_error(summarize_covariate_effects(f$fit, f$post, change = c(0, 1)), "nonzero")
  expect_error(summarize_covariate_effects(f$fit, f$post, baseline = 0), "one entry")
  f$fit$data$z <- matrix(numeric(), 2, 0)
  expect_error(summarize_covariate_effects(f$fit, f$post), "respondent covariates")
})

test_that("state and covariate summaries reuse cache without changing results", {
  f <- covariate_fixture()
  cache <- prepare_partition_cache(f$fit)
  state <- summarize_state_parameters(f$fit,f$post)
  expect_identical(summarize_state_parameters(f$fit,f$post,partition_cache=cache),state)
  before <- cache$alignment
  covariates <- summarize_covariate_effects(f$fit,f$post)
  expect_identical(summarize_covariate_effects(f$fit,f$post,partition_cache=cache),covariates)
  expect_identical(cache$alignment,before)
})

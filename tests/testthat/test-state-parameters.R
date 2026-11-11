state_parameter_fixture <- function() {
  data <- irtihmm_data(matrix(0, 2, 3), matrix(0, 2, 3))
  paths <- list(c(2L, 2L, 5L, 5L, 2L, 5L),
                c(4L, 4L, 1L, 1L, 4L, 1L),
                c(2L, 5L, 5L, 5L, 2L, 5L))
  draws <- lapply(seq_along(paths), function(j) {
    e <- list(a = rep(1000, 5), mu = rep(1000, 5),
              log_sigma = rep(log(.8), 5), rho = rep(.1, 5))
    states <- if (j == 2) c(4, 1) else c(2, 5)
    e$a[states] <- c(-1, 2) + .5 * (j - 1)
    e$mu[states] <- c(-.5, .5) + .1 * (j - 1)
    e$log_sigma[states] <- log(c(.2, .5))
    e$rho[states] <- c(-.4, .6)
    list(states = matrix(paths[[j]], 2, 3, byrow = TRUE), emission = e)
  })
  fit <- structure(list(data = data, draws = draws), class = "irtihmm_fit")
  post <- structure(list(chips = list(subpartition =
    matrix(c(1L, 1L, 2L, 2L, NA, NA), 2, 3, byrow = TRUE))),
    class = "irtihmm_postprocess")
  sim <- list(data = data, states = matrix(c(1L, 1L, 2L, 2L, 1L, 2L), 2, 3, byrow = TRUE),
    emission = list(a = c(-1, 2), mu = c(-.5, .5),
                    log_sigma = log(c(.2, .5)), rho = c(-.4, .6)))
  list(fit = fit, post = post, sim = sim)
}

test_that("CHIPS parameter summaries condition jointly and follow draw-specific labels", {
  f <- state_parameter_fixture()
  x <- summarize_state_parameters(f$fit, f$post, level = .9, simulation = f$sim)
  expect_identical(x$compatible_draws, c(1L, 2L))
  expect_equal(unname(x$matched_states), rbind(c(2, 5), c(4, 1)))
  expect_equal(x$core_probability, 2 / 3)
  expect_equal(x$joint_probability, .6)
  expect_equal(x$draws$a[, 1], c(-1, -.5))
  expect_equal(x$draws$a[, 2], c(2, 2.5))
  expect_equal(unname(x$draws$a_contrast), rbind(c(0, 3), c(0, 3)))
  expect_equal(unname(x$draws$mu_contrast), rbind(c(0, 1), c(0, 1)))
  expect_equal(x$summary$truth[x$summary$parameter == "rho"], c(-.4, .6))
  expect_equal(x$mapping$true_state, c(1L, 2L))
  file <- tempfile(fileext = ".png")
  on.exit(unlink(file), add = TRUE)
  expect_identical(plot_state_parameters(x, file), x)
  expect_gt(file.info(file)$size, 0)
  plot_state_parameters(x, file, effects = "absolute")
})

test_that("contrasts survive joint state relabeling and arbitrary location shifts", {
  f <- state_parameter_fixture()
  before <- summarize_state_parameters(f$fit, f$post)
  permutation <- c(5L, 3L, 1L, 2L, 4L)
  for (j in seq_along(f$fit$draws)) {
    d <- f$fit$draws[[j]]
    d$states[] <- permutation[d$states]
    for (name in names(d$emission)) {
      values <- d$emission[[name]]
      d$emission[[name]][permutation] <- values
    }
    d$emission$a <- d$emission$a + 10 * j
    d$emission$mu <- d$emission$mu - 3 * j
    f$fit$draws[[j]] <- d
  }
  after <- summarize_state_parameters(f$fit, f$post)
  expect_equal(after$draws$a_contrast, before$draws$a_contrast)
  expect_equal(after$draws$mu_contrast, before$draws$mu_contrast)
  expect_equal(after$draws$sigma, before$draws$sigma)
  expect_equal(after$draws$rho, before$draws$rho)
})

test_that("mixed true states are not silently matched and invalid inputs fail", {
  f <- state_parameter_fixture()
  f$sim$states[1, 2] <- 2L
  x <- summarize_state_parameters(f$fit, f$post, simulation = f$sim)
  expect_true(is.na(x$mapping$true_state[1]))
  expect_true(is.na(x$summary$truth[x$summary$parameter == "rho"][1]))
  expect_true(all(is.na(x$summary$truth[x$summary$parameter == "a_contrast"])))
  expect_error(summarize_state_parameters(f$fit, f$post, level = 1), "level")
  expect_error(summarize_state_parameters(f$fit, f$post, reference = 9), "reference")
  f$post$chips$probability <- .95
  expect_error(summarize_state_parameters(f$fit, f$post), "do not match")
  f$post$chips$subpartition[,] <- NA
  expect_error(summarize_state_parameters(f$fit, f$post), "positive integer")
})

test_that("a single core state has zero reference contrasts", {
  fit <- structure(list(data = irtihmm_data(matrix(1, 1, 1), matrix(0, 1, 1)),
    draws = list(list(states = matrix(4L, 1, 1), emission =
      list(a = 1:4, mu = 1:4, log_sigma = rep(log(.2), 4), rho = rep(0, 4))))),
    class = "irtihmm_fit")
  post <- structure(list(chips = list(subpartition = matrix(1L, 1, 1))),
                    class = "irtihmm_postprocess")
  x <- summarize_state_parameters(fit, post)
  expect_equal(x$draws$a_contrast, matrix(0, 1, 1, dimnames = list(NULL, "1")))
  expect_equal(x$summary$lower[x$summary$parameter == "sigma"], .2)
  file <- tempfile(fileext = ".png")
  on.exit(unlink(file), add = TRUE)
  plot_state_parameters(x, file)
  expect_gt(file.info(file)$size, 0)
})

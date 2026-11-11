trace_fixture <- function() {
  data <- irtihmm_data(matrix(0, 2, 2), matrix(0, 2, 2), cbind(group = c(0, 1), score = c(-1, 1)))
  draws <- lapply(seq_along(c(1L, 3L, 2L)), function(j) {
    K <- c(1L, 3L, 2L)[j]
    list(iteration = as.integer(10 + j), states = matrix(K, 2, 2),
      emission = list(b = c(-1, 1), nu = c(2, 3), d = c(-.2, .2), eta = c(.1, -.1),
        a = seq_len(K), mu = seq_len(K) / 2, log_sigma = rep(log(.25), K), rho = rep(.2, K)),
      transition = list(alpha = matrix(0, K, K), alpha0 = rep(0, K), beta = matrix(0, 2, K),
        sigma_alpha2 = 1, mean_diagonal = 0, var_diagonal = 1),
      hyper = .initial_hyper(irtihmm_priors()),
      locations = list(accuracy = list(lambda = c(1, 2)), time = list(lambda = c(3, 4))))
  })
  structure(list(data = data, draws = draws, gustafson = TRUE), class = "irtihmm_fit")
}

test_that("traces cover all stored model coordinates and retain variable-prefix gaps", {
  fit <- trace_fixture()
  traces <- .parameter_trace_series(fit)
  names <- vapply(traces$series, `[[`, character(1), "parameter")
  kinds <- vapply(traces$series, `[[`, character(1), "kind")
  expect_equal(sum(kinds == "parameter"), 55)
  expect_equal(sum(kinds == "derived"), 14)
  expect_false(anyDuplicated(names) > 0)
  expect_true(all(c("b[1]", "nu[2]", "d[1]", "eta[2]", "a[state=3]", "mu[state=3]",
    "log_sigma[state=3]", "sigma[state=3]", "rho[state=3]", "beta[group,stick=3]",
    "beta[score,stick=3]", "alpha0[stick=3]", "alpha[origin=3,stick=3]",
    "time.variance_person", "transition.sigma_alpha2", "transition.mean_diagonal",
    "transition.var_diagonal", "lambda_time[2]") %in% names))
  expect_identical(traces$series[[match("a[state=3]", names)]]$values, c(NA_real_, 3, NA_real_))
  expect_equal(traces$series[[match("beta[group,stick=3]", names)]]$values, c(NA, 0, NA))
  expect_equal(traces$series[[match("count[state=3]", names)]]$values, c(0, 4, 0))
})

test_that("trace PDFs have a complete index and CHIPS conditional gaps", {
  fit <- trace_fixture()
  sp <- structure(list(compatible_draws = c(1L, 3L),
    draws = list(a = matrix(c(1, 2), 2, 1, dimnames = list(NULL, "1")))),
    class = "irtihmm_state_parameters")
  traces <- .parameter_trace_series(fit, sp)
  conditional <- traces$series[[length(traces$series)]]
  expect_identical(conditional$values, c(1, NA_real_, 2))
  file <- tempfile(fileext = ".pdf")
  csv <- sub(".pdf", "_index.csv", file, fixed = TRUE)
  on.exit(unlink(c(file, csv)), add = TRUE)
  before <- grDevices::dev.cur()
  index <- plot_parameter_traces(fit, file, sp)
  expect_identical(grDevices::dev.cur(), before)
  expect_equal(nrow(index), 70)
  expect_equal(nrow(utils::read.csv(csv)), nrow(index))
  expect_equal(readChar(file, 5, useBytes = TRUE), "%PDF-")
  expect_equal(index$present_draws[index$parameter == "a[state=3]"], 1)
  expect_true(all(index$page >= 1 & index$panel <= 9))
  single <- fit; single$draws <- single$draws[1]
  single$draws[[1]]$locations <- list()
  expect_silent(plot_parameter_traces(single, file))
  expect_error(plot_parameter_traces(fit, "bad.png"), "PDF")
})

effect_fixture <- function() {
  x <- matrix(0, 3, 4, dimnames = list(paste0("person", 1:3), paste0("item", 1:4)))
  data <- irtihmm_data(x, matrix(0, 3, 4), cbind(group = c(0, 1, 1), score = c(-1, 0, 1)))
  draws <- lapply(1:4, function(j) {
    b <- if (j == 4) c(20, 0, 0, 0) else 1:4 + 10 * j
    C <- matrix(rep(c(1L, 1L, 2L, 2L), 3), 3, 4, byrow = TRUE)
    if (j == 4) C[1, 1] <- 2L
    list(states = C, emission = list(b = b, nu = c(1, 3, 4, 8) - j,
      d = c(-1, 0, 2) + j, eta = c(.1, -.2, .5) + j))
  })
  fit <- structure(list(data = data, draws = draws), class = "irtihmm_fit")
  core <- matrix(NA_integer_, 3, 4)
  core[1, ] <- c(1L, 1L, NA, 2L)
  post <- structure(list(chips = list(subpartition = core, probability = .75)),
                    class = "irtihmm_postprocess")
  sim <- list(data = data, emission = draws[[1]]$emission)
  list(fit = fit, post = post, sim = sim)
}

test_that("effects are centered per draw and retain fixed IDs and centered truth", {
  f <- effect_fixture()
  before <- f$fit
  x <- summarize_effects(f$fit, simulation = f$sim)
  expect_identical(f$fit, before)
  expect_identical(x$draw_indices, 1:4)
  expect_null(x$core_probability)
  for (values in x$draws) expect_equal(rowMeans(values), rep(0, 4), tolerance = 1e-12)
  expect_equal(x$draws$b[1, ], c(-1.5, -.5, .5, 1.5))
  expect_equal(x$items$accuracy_truth, c(-1.5, -.5, .5, 1.5))
  expect_equal(x$items$label, paste0("item", 1:4))
  expect_equal(x$subjects$label, paste0("person", 1:3))
  expect_equal(x$items$accuracy_lower[1], unname(quantile(x$draws$b[, 1], .025)))
  for (j in seq_along(f$fit$draws)) for (name in c("b", "nu", "d", "eta"))
    f$fit$draws[[j]]$emission[[name]] <- f$fit$draws[[j]]$emission[[name]] + 4 * j
  shifted <- summarize_effects(f$fit, simulation = f$sim)
  expect_equal(shifted$draws, x$draws, tolerance = 1e-12)
  expect_equal(shifted$items, x$items, tolerance = 1e-12)
  expect_equal(shifted$subjects, x$subjects, tolerance = 1e-12)
})

test_that("optional CHIPS conditioning uses the joint core event", {
  f <- effect_fixture()
  x <- summarize_effects(f$fit, f$post)
  expect_identical(x$draw_indices, 1:3)
  expect_equal(x$core_probability, .75)
  expect_equal(x$items$accuracy_lower, c(-1.5, -.5, .5, 1.5))
  expect_true(all(is.na(x$items$accuracy_truth)))
  f$post$chips$probability <- .9
  expect_error(summarize_effects(f$fit, f$post), "do not match")
  f$post$chips$subpartition[,] <- NA
  expect_error(summarize_effects(f$fit, f$post), "positive integer")
})

test_that("effect plots support truth, covariates, conditioning and singleton effects", {
  f <- effect_fixture()
  x <- summarize_effects(f$fit, simulation = f$sim)
  file <- tempfile(fileext = ".png")
  on.exit(unlink(file), add = TRUE)
  expect_identical(plot_effects(x, "items", file), x)
  expect_gt(file.info(file)$size, 0)
  plot_effects(x, "subjects", file, covariate = "group")
  plot_effects(x, "subjects", file, covariate = NULL)
  expect_error(plot_effects(x, "subjects", file, covariate = "absent"), "covariate")
  plot_effects(summarize_effects(f$fit, f$post), "subjects", file)
  # Continuous display grouping has no effect on the posterior summaries.
  x$covariates <- matrix(seq_len(8), 8, 1)
  x$subjects <- x$subjects[rep(1:3, length.out = 8), ]
  plot_effects(x, "subjects", file)
  single <- structure(list(data = irtihmm_data(matrix(1, 1, 1), matrix(0, 1, 1)),
    draws = list(list(emission = list(b = 3, nu = 2, d = 1, eta = 0)))),
    class = "irtihmm_fit")
  x <- summarize_effects(single)
  expect_equal(x$items$accuracy_lower, 0)
  expect_equal(x$subjects$time_upper, 0)
  plot_effects(x, "items", file)
  plot_effects(x, "subjects", file)
})

test_that("invalid effect dimensions, levels and simulation data are rejected", {
  f <- effect_fixture()
  expect_error(summarize_effects(f$fit, level = 1), "level")
  f$sim$data$x[1, 1] <- 1
  expect_error(summarize_effects(f$fit, simulation = f$sim), "same data")
  f$fit$draws[[1]]$emission$b <- NA_real_
  expect_error(summarize_effects(f$fit), "finite b, nu, d and eta")
})

test_that("conditional effects reuse partition compatibility without changing summaries", {
  data <- irtihmm_data(matrix(0,2,2),matrix(0,2,2))
  fit <- structure(list(data=data,draws=lapply(1:4,function(j)list(
    states=matrix(c(1L,2L,1L,2L),2,2),emission=list(
      b=c(-.2,.2)+j/10,nu=c(.5,1)+j/10,d=c(-.1,.1)+j/10,eta=c(-.2,.2)+j/10)))),
    class='irtihmm_fit')
  post <- structure(list(chips=list(subpartition=matrix(c(1L,2L,1L,2L),2,2),probability=1)),
    class='irtihmm_postprocess')
  cache <- prepare_partition_cache(fit)
  expect_identical(summarize_effects(fit,post,partition_cache=cache),summarize_effects(fit,post))
})

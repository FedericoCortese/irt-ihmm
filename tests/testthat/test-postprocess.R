partition_fixture <- function(variable = FALSE) {
  data <- irtihmm_data(matrix(0, 2, 3), matrix(0, 2, 3))
  rows <- rep(list(c(1L, 1L, 2L, 2L, 1L, 2L)), 12)
  if (variable) for (j in 1:4) rows[[j]][6] <- 3L
  structure(list(data = data, draws = lapply(rows, function(x)
    list(states = matrix(x, 2, 3, byrow = TRUE)))), class = "irtihmm_fit")
}

test_that("respondent-item mapping and temporal probabilities survive relabeling", {
  skip_if_not_installed("salso", "0.3.80")
  fit <- partition_fixture()
  original <- fit
  set.seed(9)
  rng <- .Random.seed
  # A deterministic two-cluster posterior reaches salso's default search bound.
  expect_warning(post <- postprocess_irtihmm(fit, n_runs = 4, seed = 3),
                 "default maximum possible number of clusters")
  expect_identical(fit, original)
  expect_identical(.Random.seed, rng)
  expect_equal(.posterior_partitions(fit)[1, ], c(1, 1, 2, 2, 1, 2))
  expect_equal(post$change_probability, matrix(c(0, 1, 1, 1), 2, byrow = TRUE))
  expect_equal(post$expected_loss, 0)
  expect_equal(post$chips$probability, 1)
  expect_equal(post$chips$core_fraction, 1)
  for (j in seq_along(fit$draws)) {
    C <- fit$draws[[j]]$states
    fit$draws[[j]]$states <- matrix(c(8L + j, 3L + j)[C], 2, 3)
  }
  expect_warning(relabeled <- postprocess_irtihmm(fit, n_runs = 4, seed = 3),
                 "default maximum possible number of clusters")
  same <- function(C) outer(as.integer(C), as.integer(C), "==")
  expect_equal(same(post$states), same(relabeled$states))
  expect_equal(same(post$chips$states), same(relabeled$chips$states))
  expect_equal(post$change_probability, relabeled$change_probability)
})

test_that("CHIPS reports a simultaneous core with variable posterior counts", {
  skip_if_not_installed("salso", "0.3.80")
  fit <- partition_fixture(TRUE)
  post <- postprocess_irtihmm(fit, credible = .95, n_runs = 8, seed = 6)
  expect_equal(sort(unique(post$occupied_states)), c(2, 3))
  expect_gte(post$chips$probability, .95)
  expect_lt(post$chips$core_fraction, 1)
  core <- as.integer(t(post$chips$subpartition))
  core[is.na(core)] <- -1L
  event <- .core_draws(.posterior_partitions(fit), core)
  expect_equal(post$chips$compatible_draws, which(event))
  expect_equal(post$chips$probability, mean(event))
  complete <- as.integer(t(post$chips$states))
  keep <- !is.na(as.integer(t(post$chips$subpartition)))
  expect_equal(outer(core[keep], core[keep], "=="),
               outer(complete[keep], complete[keep], "=="))
  # Separate clusters in the core must stay separate, not just internally equal.
  expect_identical(.core_draws(rbind(c(3L, 3L, 3L), c(8L, 9L, 8L)),
                              c(1L, 2L, -1L)), c(FALSE, TRUE))
  file <- tempfile(fileext = ".png")
  on.exit(unlink(file), add = TRUE)
  expect_identical(plot_postprocess(post, file), post)
  expect_gt(file.info(file)$size, 0)
  no_chips <- postprocess_irtihmm(fit, chips = FALSE, loss = "binder", n_runs = 2, seed = 1)
  expect_null(no_chips$chips)
  plot_postprocess(no_chips, file)
})

test_that("post-processing handles singleton data and rejects malformed states", {
  skip_if_not_installed("salso", "0.3.80")
  fit <- structure(list(data = irtihmm_data(matrix(1, 1, 1), matrix(0, 1, 1)),
    draws = list(list(states = matrix(4L, 1, 1)))), class = "irtihmm_fit")
  expect_warning(post <- postprocess_irtihmm(fit, n_runs = 1, seed = 4),
                 "default maximum possible number of clusters")
  expect_equal(dim(post$states), c(1, 1))
  expect_equal(dim(post$change_probability), c(1, 0))
  expect_equal(post$chips$probability, 1)
  expect_error(postprocess_irtihmm(fit, credible = 1.1), "credible")
  expect_error(postprocess_irtihmm(fit, n_runs = 0), "positive integers")
  fit$draws[[1]]$states[1, 1] <- NA
  expect_error(postprocess_irtihmm(fit), "positive integer labels")
})

test_that("direct CHIPS preserves all draws and caches whole-core compatibility", {
  skip_if_not_installed('salso','0.3.80')
  fit <- partition_fixture(TRUE)
  cache <- prepare_partition_cache(fit)
  post <- postprocess_irtihmm(fit, seed=6, partition_cache=cache)
  expect_equal(post$settings$n_runs,4)
  expect_identical(post$settings$chips_search,'direct')
  expect_null(post$chips$curve)
  expect_equal(nrow(cache$partitions),length(fit$draws))
  a <- .core_alignment(fit,post,cache)
  expect_identical(a$indices,post$chips$compatible_draws)
  expect_equal(a$matched,cache$partitions[a$indices,a$anchors,drop=FALSE],ignore_attr=TRUE)
  # Reusing the same core preserves the cached alignment rather than re-searching.
  cache$alignment$reuse_marker <- TRUE
  expect_true(.core_alignment(fit,post,cache)$reuse_marker)
  changed <- post
  changed$chips$subpartition[,] <- NA_integer_
  changed$chips$subpartition[1,1] <- 1L
  changed$chips$probability <- NULL
  b <- .core_alignment(fit,changed,cache)
  expect_equal(b$probability,1)
  expect_length(b$indices,length(fit$draws))
  expect_null(b$reuse_marker)
  stale <- fit;stale$draws[[1]]$states[1,1] <- 9L
  expect_error(.core_alignment(stale,post,cache),'cache does not match')
  bad <- prepare_partition_cache(fit);bad$partitions[1,1] <- 9L
  expect_error(.core_alignment(fit,post,bad),'cache does not match')
  expect_error(prepare_partition_cache(fit,parent=cache,rows=rep(1,12)),'Invalid parent')
})

test_that("chain caches share pooled storage and retain local draw indices", {
  fit <- partition_fixture(TRUE)
  pooled <- fit;pooled$draws <- rep(fit$draws,12)  # exercise the 64-draw packing blocks
  cache <- prepare_partition_cache(pooled)
  ref <- do.call(rbind,lapply(pooled$draws,function(d)as.integer(t(d$states))))
  expect_identical(cache$partitions,ref)
  rows <- 13:24
  child <- prepare_partition_cache(fit,parent=cache,rows=rows)
  expect_identical(child$partitions,cache$partitions)
  expect_identical(child$rows,as.integer(rows))
  core <- c(1L,1L,2L,2L,1L,2L)
  expect_identical(core_draws_cpp(child$partitions,core,child$rows),.core_draws(ref[rows,,drop=FALSE],core))
  other <- fit;other$draws[[1]]$states[1,1] <- 7L
  expect_error(prepare_partition_cache(other,parent=cache,rows=rows),'does not match')
  # Numeric state matrices remain supported and are checked before integer packing.
  other <- fit
  for(j in seq_along(other$draws))storage.mode(other$draws[[j]]$states)<-'double'
  expect_identical(.posterior_partitions(other),.posterior_partitions(fit))
})

test_that("deterministic partition kernels do not initialize the caller RNG", {
  fit <- partition_fixture(TRUE)
  previous <- get0('.Random.seed',envir=.GlobalEnv,inherits=FALSE)
  on.exit(if(is.null(previous)){
    if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)
  } else assign('.Random.seed',previous,envir=.GlobalEnv),add=TRUE)
  if(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))rm('.Random.seed',envir=.GlobalEnv)
  cache <- prepare_partition_cache(fit)
  expect_false(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))
  suppressWarnings(postprocess_irtihmm(fit,n_runs=1,seed=5,partition_cache=cache))
  expect_false(exists('.Random.seed',envir=.GlobalEnv,inherits=FALSE))
})

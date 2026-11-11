test_that("ROAR exports omega and diagnostics without any parameter PDF", {
  skip_if_not_installed("posterior")
  skip_if_not_installed("salso")
  skip_if_not("chips" %in% getNamespaceExports("salso"))
  # This helper is an application script, rather than an exported package function.
  root <- if (file.exists("tools/roar_fit_outputs.R")) "." else "../.."
  skip_if_not(file.exists(file.path(root, "R/plots.R")),
    "Application output integration requires the source repository.")
  original_wd <- getwd()
  setwd(root)
  on.exit(setwd(original_wd), add = TRUE)
  env <- new.env(parent = environment())
  sys.source("tools/roar_fit_outputs.R", envir = env)
  sim <- simulate_irtihmm(4, 5, 2,
    z = cbind(z1 = c(0, 1, 0, 1), z2 = c(-1, -.3, .3, 1)), seed = 12)
  chains <- lapply(1:3, function(seed) fit_irtihmm(sim$data, mode = "finite", K = 2,
    iterations = 24, burnin = 4, states = sim$states, seed = seed, persistence = TRUE))
  settings <- list(seeds = 1:3, priors = irtihmm_priors(), persistence = TRUE,
    omega_prior = c(1, 1), excluded_subjects = 2L, time_policy = "Test fixture.",
    transitions = "Known paths for export validation.")
  path <- tempfile("roar_outputs_")
  dir.create(path)
  on.exit(unlink(path, recursive = TRUE), add = TRUE)
  # Fail if a trace exporter is accidentally reached with both flags disabled.
  env$plot_parameter_traces <- function(...) stop("Unexpected parameter PDF export.")
  suppressMessages(env$write_roar_fit_outputs(chains, settings, path,
    traceplots = FALSE, multichain_traceplots = FALSE, post_n_runs = 1L, post_cores = 1L))
  expect_length(list.files(path, pattern = "\\.pdf$", recursive = TRUE), 0L)
  expect_true(file.exists(file.path(path, "omega.png")))
  omega <- read.csv(file.path(path, "omega_summary.csv"))
  expect_equal(nrow(omega), 4L)
  diagnostics <- read.csv(file.path(path, "multichain_diagnostics.csv"))
  expect_true("omega" %in% diagnostics$variable)
  saved <- readRDS(file.path(path, "posterior_summaries.rds"))
  expect_equal(saved$omega_summary$mean[1],
    mean(vapply(chains[[1]]$draws, function(d) d$transition$omega, numeric(1))))
  expect_true(any(grepl("multichain parameter PDF: OFF", readLines(file.path(path, "fit_report.txt")))))

  # OFF must bypass native searches and partition packing, even with existing results.
  cache_file <- file.path(path, "postprocess_direct_runs1_prob0.95.rds")
  cached <- readBin(cache_file, "raw", n = file.info(cache_file)$size)
  previous <- readBin(file.path(path, "posterior_summaries.rds"), "raw",
    n = file.info(file.path(path, "posterior_summaries.rds"))$size)
  for (name in c("prepare_partition_cache", "postprocess_irtihmm",
                 "summarize_state_parameters", "summarize_covariate_effects"))
    assign(name, function(...) stop("Unexpected partition processing in OFF mode."), env)
  off <- suppressMessages(env$write_roar_fit_outputs(chains, settings, path,
    traceplots = FALSE, multichain_traceplots = FALSE, post_processing = FALSE))
  expect_null(off$post)
  saved <- readRDS(file.path(path, "fit_summaries.rds"))
  expect_false(saved$post_processing)
  expect_null(saved$state_parameters)
  expect_true(file.exists(file.path(path, "item_effects.png")))
  expect_true(file.exists(file.path(path, "chain_comparison.png")))
  expect_identical(readBin(cache_file, "raw", n = file.info(cache_file)$size), cached)
  expect_identical(readBin(file.path(path, "posterior_summaries.rds"), "raw",
    n = file.info(file.path(path, "posterior_summaries.rds"))$size), previous)
  expect_true(any(grepl("post-processing: OFF", readLines(file.path(path, "fit_report.txt")), fixed = TRUE)))
})

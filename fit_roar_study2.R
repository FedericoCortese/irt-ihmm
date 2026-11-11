# Full ROAR XYZ fit: 3 independent iHMM chains, then diagnostics and all summaries.
# Run source("fit_roar_study2.R") from the project directory. See overwrite_fits below.
source("load.R")
iterations <- 5000L
burnin <- 1000L
thin <- 1L
seeds <- c(101L, 202L, 303L)
gustafson <- TRUE
estimate_persistence <- FALSE
omega_prior <- c(3, 1)  # Weak preference for persistence: mean 0.75, SD about 0.194.
# Manual prior tuning: m_alpha ~ N(m_alpha_mean, m_alpha_sd^2).
# Defaults reproduce N(0,1); this controls diagonal sticks/persistence.
m_alpha_mean <- 0
m_alpha_sd <- 1
priors <- irtihmm_priors(m_alpha_mean = m_alpha_mean, m_alpha_sd = m_alpha_sd)
fixed_params <- list()
traceplots <- FALSE  # Skip the large per-chain parameter PDFs for this run.
multichain_traceplots <- FALSE  # Independent flag for the large overlay PDF.
# Post-processing settings do not alter or invalidate the sampled chains.
post_processing <- TRUE  # Enable SALSO/CHIPS and their aligned state/covariate summaries.
overwrite_fits <- TRUE  # Set TRUE to replace this working directory; FALSE reuses matching saved chains.
post_n_runs <- 4L
post_credible <- .90
post_cores <- 3L
output_dir <- file.path("output", "roar_study2", "fit")
if (!requireNamespace("posterior", quietly = TRUE)) stop("The posterior package is required for diagnostics.")
if (post_processing && !requireNamespace("salso", quietly = TRUE))
  stop("The salso package is required for partition summaries.")
prepared <- readRDS(file.path("output", "roar_study2", "prepared_data.rds"))
if (basename(prepared$source$csv) != "roar_study2_xyz.csv") stop("Prepare the requested XYZ file first.")
# Complete-time respondents only: keep the actual chronology, never sort by item.
keep <- prepared$availability$complete_log_times
if (!all(keep)) message("Excluding incomplete-time subjects: ", paste(prepared$subject_ids[!keep], collapse = ", "))
map <- prepared$chronological$item_id[keep, , drop = FALSE]
map[] <- match(map, prepared$item_ids)
# Map column names label the common item effects, not the chronological positions.
colnames(map) <- as.character(prepared$item_ids)
data <- irtihmm_data(prepared$chronological$x[keep, , drop = FALSE],
  prepared$chronological$y[keep, , drop = FALSE],
  prepared$chronological$z[keep, , drop = FALSE], item_index = map)
settings <- list(iterations = iterations, burnin = burnin, thin = thin, seeds = seeds,
  gustafson = gustafson, fixed_params = fixed_params, copula = TRUE, initial_K = 2L,
  source = prepared$source, included_subjects = prepared$subject_ids[keep],
  excluded_subjects = prepared$subject_ids[!keep], scaling = prepared$scaling,
  standardized = prepared$standardize_covariates,
  time_policy = "Exclude the respondent with missing log-times; no trimming or imputation.",
  transitions = "One continuous chronological HMM sequence per respondent; no block resets.")
settings$priors <- priors
settings$persistence <- estimate_persistence
settings$omega_prior <- if (estimate_persistence) omega_prior else NULL
# Input validation finishes before replacing any saved output. The source data live elsewhere.
if (overwrite_fits && dir.exists(output_dir)) {
  message("Replacing previous ROAR outputs in ", output_dir)
  unlink(output_dir, recursive = TRUE)
  if (dir.exists(output_dir)) stop("Could not clear previous ROAR outputs.")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
saveRDS(list(data = data, settings = settings), file.path(output_dir, "analysis_input.rds"))
write.csv(prepared$availability, file.path(output_dir, "sample_inclusion.csv"), row.names = FALSE)

run_chain <- function(chain) {
  path <- file.path(output_dir, paste0("chain_", chain))
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  result_file <- file.path(path, "fit.rds")
  if (file.exists(result_file)) {
    old <- readRDS(result_file)
    # Legacy chains used the implicit default prior; never reuse them for tuning.
    expected <- settings
    if (is.null(old$settings$persistence) && !estimate_persistence)
      expected$persistence <- NULL
    if (is.null(old$settings$priors) && identical(priors, irtihmm_priors()))
      expected$priors <- NULL
    if (!identical(old$settings, expected) || !identical(old$fit$data, data) ||
        !identical(.validate_priors(old$fit$priors), priors) ||
        isTRUE(old$fit$persistence) != estimate_persistence ||
        (estimate_persistence && !identical(old$fit$omega_prior, omega_prior)))
      stop("Existing chain does not match current settings/data: ", result_file)
    return(result_file)
  }
  log <- file(file.path(path, "sampling.log"), open = "wt")
  sink(log); sink(log, type = "message")
  on.exit({sink(type = "message"); sink(); close(log)}, add = TRUE)
  fit <- fit_irtihmm(data, mode = "infinite", K = 2,
    iterations = iterations, burnin = burnin, thin = thin, seed = seeds[chain],
    copula = TRUE, progress = TRUE, gustafson = gustafson, fixed_params = fixed_params,
    priors = priors, persistence = estimate_persistence, omega_prior = omega_prior)
  saveRDS(list(fit = fit, settings = settings), paste0(result_file, ".tmp"))
  if (!file.rename(paste0(result_file, ".tmp"), result_file)) stop("Could not finalize chain file.")
  print(fit)
  result_file
}
# Sources and the Rcpp cache are loaded before forks to avoid concurrent compilation.
chain_paths <- parallel::mclapply(seq_along(seeds), run_chain,
  mc.cores = 3L, mc.preschedule = FALSE, mc.set.seed = FALSE)
if (any(vapply(chain_paths, inherits, logical(1), "try-error"))) stop("A chain failed; inspect sampling.log.")
chains <- lapply(chain_paths, function(path) readRDS(path)$fit)
source("tools/roar_fit_outputs.R")
write_roar_fit_outputs(chains, settings, output_dir, traceplots,
  post_n_runs = post_n_runs, post_credible = post_credible, post_cores = post_cores,
  multichain_traceplots = multichain_traceplots, post_processing = post_processing)

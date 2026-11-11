# 3 generating regimes, low reference persistence and explicit Gustafson flag.
# Run from the project directory, optionally one section at a time in RStudio.
source("load.R")
if (!requireNamespace("salso", quietly = TRUE) ||
    !"chips" %in% getNamespaceExports("salso"))
  stop('Install the current salso once: install.packages("salso").')

# SETTINGS: K_true controls the generating process; the fitted iHMM remains infinite.
K_true <- 3L
gustafson <- TRUE
estimate_persistence <- FALSE  # Estimate omega; FALSE recovers the original probit transitions.
omega_prior <- c(1, 1)  # Beta shapes: uniform prior on forced persistence.
compare_persistence <- FALSE  # Reuse or fit an OFF baseline on exactly the same toy.
# Manual prior tuning: m_alpha ~ N(m_alpha_mean, m_alpha_sd^2).
# Defaults reproduce N(0,1); this controls diagonal sticks/persistence.
m_alpha_mean <- 0
m_alpha_sd <- 1
priors <- irtihmm_priors(m_alpha_mean = m_alpha_mean, m_alpha_sd = m_alpha_sd)
prior_tag <- if (m_alpha_mean == 0 && m_alpha_sd == 1) NULL else
  sprintf("m_alpha_mean_%s_sd_%s", format(m_alpha_mean, digits = 15, trim = TRUE),
          format(m_alpha_sd, digits = 15, trim = TRUE))
fixed_params <- list()  # Optional fixed emission blocks, e.g. list(eta = 0).
traceplots <- FALSE  # Set TRUE to export all parameter traceplots and their index.
iterations <- 5000L
burnin <- 1000L
demo_output <- file.path("output", paste0("demo_k", K_true),
                         if (gustafson) "gustafson_on" else "gustafson_off", "eta_fixed_zero")
if (!is.null(prior_tag)) demo_output <- file.path(demo_output, prior_tag)
dir.create(demo_output, recursive = TRUE, showWarnings = FALSE)

# 1. SIMULATE: two respondent covariates and reduced reference persistence.
# Prior tuning changes inference only; generating transition parameters stay fixed.
P <- 30L
I <- 20L
set.seed(18)
z <- cbind(group = rbinom(P, size = 1, prob = .5), score = rnorm(P))
reference_transition <- matrix(.30 / (K_true - 1L), K_true, K_true)
diag(reference_transition) <- .70
beta <- cbind(c(-1.2, .8), c(.8, -.6))
sim <- simulate_irtihmm(
  P = P, I = I, K = K_true, z = z, seed = 19, beta = beta,
  transition_matrix = reference_transition,
  parameters = list(a = seq(-1.2, 1.2, length.out = K_true),
    mu = seq(-.7, .7, length.out = K_true), log_sigma = rep(log(.25), K_true),
    rho = seq(-.4, .4, length.out = K_true))
)
stopifnot(isTRUE(all.equal(
  .transition_probs(sim$transition, matrix(0, 1, ncol(z)))$transition[[1]],
  reference_transition
)))
# At nonzero covariates persistence can differ from 0.70: report it explicitly.
self <- t(vapply(sim$probabilities$transition, diag, numeric(K_true)))
persistence <- data.frame(state = seq_len(K_true),
  responses = as.integer(table(factor(sim$states, levels = seq_len(K_true)))),
  minimum_self = apply(self, 2, min), median_self = apply(self, 2, median),
  maximum_self = apply(self, 2, max))
subject_state_changes <- rowSums(sim$states[, -1, drop = FALSE] !=
                                 sim$states[, -I, drop = FALSE])
utils::write.csv(persistence, file.path(demo_output, "persistence.csv"), row.names = FALSE)
utils::write.csv(data.frame(subject = seq_len(P), changes = subject_state_changes),
  file.path(demo_output, "subject_state_changes.csv"), row.names = FALSE)
print(persistence)
message("Subjects with no true state change: ", sum(subject_state_changes == 0), "/", P)
plot_simulation(sim, file.path(demo_output, "ex_ante.png"))
if (interactive()) plot_simulation(sim)

# Keep the pre-extension result without creating another output directory.
previous_file <- file.path(demo_output, "results.rds")
baseline_file <- file.path(demo_output, "results_persistence_off.rds")
matching_baseline <- function(x) {
  !isTRUE(x$fit$persistence) && identical(x$sim$data, sim$data) &&
    identical(.validate_priors(x$fit$priors), priors) &&
    identical(x$fit$fixed_params, fixed_params) &&
    identical(x$fit$gustafson, gustafson) &&
    nrow(x$fit$diagnostics) == iterations &&
    identical(vapply(x$fit$draws, `[[`, integer(1), "iteration"),
              seq.int(burnin + 1L, iterations))
}
if (estimate_persistence && file.exists(previous_file)) {
  previous <- readRDS(previous_file)
  if (!isTRUE(previous$fit$persistence)) {
    backup_file <- file.path(demo_output, "results_before_persistence.rds")
    if (!file.exists(backup_file)) file.copy(previous_file, backup_file)
    if (matching_baseline(previous)) file.copy(previous_file, baseline_file, overwrite = TRUE)
  }
  rm(previous)
}
if (estimate_persistence && compare_persistence) {
  baseline <- if (file.exists(baseline_file)) readRDS(baseline_file) else NULL
  if (is.null(baseline) || !matching_baseline(baseline)) {
    message("Fitting the persistence-OFF baseline on the same simulated data.")
    baseline_fit <- fit_irtihmm(sim$data, mode = "infinite", K = 2,
      iterations = iterations, burnin = burnin, seed = 23, progress = TRUE,
      gustafson = gustafson, fixed_params = fixed_params, priors = priors,
      persistence = FALSE)
    saveRDS(list(sim = sim, fit = baseline_fit), baseline_file)
    rm(baseline_fit)
  }
  rm(baseline)
}

# 2. FIT: K=2 is an initial prefix, not a fixed number of fitted states.
fit <- fit_irtihmm(sim$data, mode = "infinite", K = 2,
  iterations = iterations, burnin = burnin, seed = 23,
  progress = TRUE, gustafson = gustafson, fixed_params = fixed_params, priors = priors,
  persistence = estimate_persistence, omega_prior = omega_prior)
saveRDS(list(sim = sim, fit = fit), file.path(demo_output, "fit_checkpoint.rds"))
plot_fit(fit, simulation = sim, file = file.path(demo_output, "ex_post.png"))
if (interactive()) plot_fit(fit, simulation = sim)

# 3. PARTITIONS AND STATE PARAMETERS: SALSO and CHIPS-conditional summaries.
partition_cache <- prepare_partition_cache(fit)
post <- postprocess_irtihmm(fit, credible = .95, n_runs = 4, seed = 24,
  partition_cache = partition_cache)
plot_postprocess(post, file.path(demo_output, "post_processing.png"))
state_parameters <- summarize_state_parameters(fit, post, simulation = sim, partition_cache = partition_cache)
plot_state_parameters(state_parameters, file.path(demo_output, "state_parameters.png"))
plot_state_parameters(state_parameters, file.path(demo_output, "state_parameters_absolute.png"),
                      effects = "absolute")
if (interactive()) {
  plot_postprocess(post)
  plot_state_parameters(state_parameters)
}

# 4. ITEM, SUBJECT AND COVARIATE EFFECTS.
effects <- summarize_effects(fit, simulation = sim)
plot_effects(effects, "items", file.path(demo_output, "item_effects.png"))
plot_effects(effects, "subjects", file.path(demo_output, "subject_effects.png"))
covariate_effects <- summarize_covariate_effects(fit, post, partition_cache = partition_cache)
plot_covariate_effects(covariate_effects, "coefficients",
                       file.path(demo_output, "covariate_coefficients.png"))
plot_covariate_effects(covariate_effects, "probabilities",
                       file.path(demo_output, "covariate_probabilities.png"))
if (interactive()) {
  plot_effects(effects, "items")
  plot_effects(effects, "subjects")
  plot_covariate_effects(covariate_effects, "coefficients")
  plot_covariate_effects(covariate_effects, "probabilities")
}

# Omega has its own compact diagnostic, independent of the all-parameter PDF flag.
omega_summary <- NULL
if (estimate_persistence) {
  omega <- vapply(fit$draws, function(d) d$transition$omega, numeric(1))
  omega_summary <- data.frame(mean = mean(omega), median = median(omega),
    lower = unname(quantile(omega, .025)), upper = unname(quantile(omega, .975)))
  utils::write.csv(omega_summary, file.path(demo_output, "omega_summary.csv"), row.names = FALSE)
  png(file.path(demo_output, "omega.png"), width = 1250, height = 550, res = 120)
  par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
  plot(vapply(fit$draws, `[[`, integer(1), "iteration"), omega, type = "l",
    col = "#244b78", xlab = "Iteration (retained draws)", ylab = "Omega",
    ylim = c(0, 1), main = "Forced persistence: trace")
  hist(omega, breaks = 35, probability = TRUE, xlim = c(0, 1),
    col = "#c9dced", border = "white", xlab = "Omega",
    main = "Forced persistence: posterior")
  curve(dbeta(x, omega_prior[1], omega_prior[2]), add = TRUE, col = "#bb4422", lwd = 2)
  legend("topright", "Prior density", col = "#bb4422", lty = 1, bty = "n")
  dev.off()
  print(omega_summary)
}
# Compare occupied counts only when the baseline uses identical data and settings.
if (estimate_persistence && file.exists(baseline_file)) {
  baseline <- readRDS(baseline_file)
  if (matching_baseline(baseline)) {
    count_summary <- function(f, label) {
      k <- summary(f)$occupied_states
      data.frame(model = label, draws = length(k), mean = mean(k), median = median(k),
        lower = unname(quantile(k, .025)), upper = unname(quantile(k, .975)))
    }
    comparison <- rbind(count_summary(baseline$fit, "Persistence OFF"),
                        count_summary(fit, "Persistence ON"))
    utils::write.csv(comparison, file.path(demo_output, "persistence_comparison.csv"), row.names = FALSE)
    print(comparison)
  }
  rm(baseline)
}

# 5. ALL PARAMETERS: paginated traces and an index, including rare stored sticks.
# Burn-in draws are not stored. CHIPS traces are conditional and retain gaps.
trace_index <- NULL
if (traceplots) {
  trace_index <- plot_parameter_traces(fit,
    file.path(demo_output, "traces", "parameters.pdf"), state_parameters = state_parameters)
}
rm(partition_cache)
saveRDS(list(sim = sim, fit = fit, post = post, state_parameters = state_parameters,
  effects = effects, covariate_effects = covariate_effects,
  reference_transition = reference_transition, persistence = persistence,
  subject_state_changes = subject_state_changes, omega_summary = omega_summary,
  trace_index = trace_index),
  file.path(demo_output, "results.rds"))
unlink(file.path(demo_output, "fit_checkpoint.rds"))
print(fit)
print(c(SALSO_states = length(unique(as.integer(post$states))),
  CHIPS_probability = post$chips$probability, CHIPS_core_fraction = post$chips$core_fraction,
  trace_panels = if (is.null(trace_index)) 0L else nrow(trace_index),
  trace_pages = if (is.null(trace_index)) 0L else max(trace_index$page)))
message("Plots and results saved in ", demo_output, "/. Parameter trace export: ",
        if (traceplots) "ON" else "OFF", ".")

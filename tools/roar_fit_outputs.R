# Outputs for the ROAR application. Sampling is performed in fit_roar_study2.R.
# Reload display helpers at export time, also when resuming saved chain fits.
source("R/plots.R")
source("R/effects.R")
write_roar_fit_outputs <- function(chains, settings, output_dir, traceplots = TRUE,
  post_n_runs = 4L, post_credible = .95, post_cores = 3L,
  multichain_traceplots = traceplots, post_processing = TRUE) {
  if (!is.logical(post_processing) || length(post_processing) != 1L || is.na(post_processing))
    stop("post_processing must be TRUE or FALSE.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  message("Sampling finished for all chains. Building multichain diagnostics.")
  M <- length(chains[[1]]$draws); n_chains <- length(chains)
  stopifnot(all(vapply(chains, function(f) length(f$draws) == M, logical(1))))
  # Original ordered coordinates and invariant contrasts are both reported.
  series <- lapply(chains, .parameter_trace_series)
  for (k in seq_along(chains)) {
    kept <- match(series[[k]]$iterations, chains[[k]]$diagnostics$iteration)
    for (name in c("occupied", "conditional_loglik")) {
      series[[k]]$series[[length(series[[k]]$series) + 1L]] <- list(
        family = "Global diagnostics", parameter = name,
        values = chains[[k]]$diagnostics[[name]][kept],
        kind = "derived", scope = "state-label invariant")
    }
  }
  ids <- lapply(series, function(s) vapply(s$series, `[[`, character(1), "parameter"))
  common <- Reduce(intersect, ids)
  finite <- vapply(common, function(name) all(vapply(seq_along(series), function(k)
    all(is.finite(series[[k]]$series[[match(name, ids[[k]])]]$values)), logical(1))), logical(1))
  common <- common[finite]
  values <- array(NA_real_, c(M, n_chains, length(common)),
    dimnames = list(iteration = NULL, chain = paste0("chain_", seq_len(n_chains)), variable = common))
  for (k in seq_len(n_chains)) for (j in seq_along(common))
    values[, k, j] <- series[[k]]$series[[match(common[j], ids[[k]])]]$values
  diagnostic_table <- as.data.frame(posterior::summarise_draws(posterior::as_draws_array(values)))
  first <- match(diagnostic_table$variable, ids[[1]])
  diagnostic_table$scope <- vapply(series[[1]]$series[first], `[[`, character(1), "scope")
  diagnostic_table$kind <- vapply(series[[1]]$series[first], `[[`, character(1), "kind")
  write.csv(diagnostic_table, file.path(output_dir, "multichain_diagnostics.csv"), row.names = FALSE)
  all_names <- Reduce(union, ids)
  unavailable <- setdiff(all_names, common)
  write.csv(data.frame(parameter = unavailable,
    reason = rep("Not instantiated in every retained draw of every chain; no padded R-hat or ESS.", length(unavailable))),
    file.path(output_dir, "diagnostics_unavailable.csv"), row.names = FALSE)
  chain_summary <- do.call(rbind, lapply(seq_along(chains), function(k) {
    f <- chains[[k]]
    data.frame(chain = k, seed = settings$seeds[k], draws = M, elapsed_seconds = f$elapsed,
      mean_occupied_states = mean(summary(f)$occupied_states),
      min_occupied_states = min(summary(f)$occupied_states),
      max_occupied_states = max(summary(f)$occupied_states))
  }))
  write.csv(chain_summary, file.path(output_dir, "chain_summary.csv"), row.names = FALSE)
  # Compare to the neighboring original-model run only when all inputs/settings match.
  if (isTRUE(settings$persistence)) {
    comparison <- list()
    expected <- settings; expected$persistence <- NULL; expected$omega_prior <- NULL
    for (k in seq_along(chains)) {
      baseline_file <- file.path(dirname(output_dir), paste0("chain_", k), "fit.rds")
      if (!file.exists(baseline_file)) next
      baseline <- readRDS(baseline_file)
      previous <- baseline$settings
      previous$persistence <- NULL; previous$omega_prior <- NULL
      if (!isTRUE(baseline$fit$persistence) && identical(previous, expected) &&
          identical(baseline$fit$data, chains[[k]]$data)) {
        counts <- summary(baseline$fit)$occupied_states
        comparison[[length(comparison) + 1L]] <- data.frame(chain = k,
          mean_occupied_OFF = mean(counts),
          mean_occupied_ON = chain_summary$mean_occupied_states[k],
          baseline = baseline_file)
      }
      rm(baseline)
    }
    if (length(comparison) == n_chains)
      write.csv(do.call(rbind, comparison), file.path(output_dir, "persistence_comparison.csv"),
        row.names = FALSE)
  }
  png(file.path(output_dir, "chain_comparison.png"), width = 1800, height = 1200, res = 150)
  tryCatch({
    par(mfrow = c(2,2), mar = c(4.5,4.5,3,1))
    colors <- c("#244b78", "#c04b39", "#009988")
    for (name in c("occupied", "conditional_loglik")) {
      value <- do.call(cbind, lapply(chains, function(f) f$diagnostics[[name]]))
      matplot(chains[[1]]$diagnostics$iteration, value, type = "l", col = colors, lty = 1,
        xlab = "Iteration", ylab = if (name == "occupied") "Occupied states" else "Conditional log-likelihood",
        main = if (name == "occupied") "Occupied states: all chains" else "Emission likelihood: all chains")
      abline(v = min(series[[1]]$iterations)-1L, lty = 3, col = "grey50")
      legend("bottomright", paste("Chain",seq_len(n_chains)), col = colors, lty = 1, bty = "n", cex = .75)
    }
    boxplot(lapply(chains, function(f) summary(f)$occupied_states),
      names = paste("Chain",seq_len(n_chains)), col = colors,
      main = "Retained occupied-state counts", ylab = "Occupied states")
    finite_rhat <- diagnostic_table$rhat[is.finite(diagnostic_table$rhat)]
    if (length(finite_rhat)) {
      hist(finite_rhat, breaks = 30, col = "#648eae", border = "white",
        main = "Rank-normalized split R-hat", xlab = "R-hat (shared complete coordinates)")
      abline(v = 1.01, col = "#c04b39", lty = 2)
    } else { plot.new(); title("No finite R-hat values") }
  }, finally = dev.off())

  # One separate PDF overlays all 3 chains on every shared, complete coordinate.
  if (multichain_traceplots) {
    grDevices::cairo_pdf(file.path(output_dir, "multichain_traces.pdf"), width = 14, height = 10)
    tryCatch({
      par(mfrow = c(3,3), mar = c(3,4,3,.5), oma = c(2.5,0,2,0), mgp = c(1.8,.5,0))
      capacity <- 9L; pages <- split(seq_along(common), ceiling(seq_along(common)/capacity))
      index <- list(); colors <- c("#244b78", "#c04b39", "#009988")
      for (page in seq_along(pages)) {
        for (slot in seq_len(capacity)) {
          if (slot > length(pages[[page]])) {plot.new(); next}
          j <- pages[[page]][slot]
          matplot(series[[1]]$iterations, values[,,j], type = "l", lty = 1, lwd = .5,
            col = colors, main = common[j], cex.main = .7, xlab = "Iteration (retained draws)", ylab = "Value")
          index[[length(index)+1L]] <- data.frame(parameter = common[j], page = page, panel = slot)
        }
        mtext("ROAR XYZ: three independent chains", outer = TRUE, side = 3, cex = 1, line = .5)
        mtext("Chain 1: blue | Chain 2: red | Chain 3: green. State/stick coordinates keep their order.",
          outer = TRUE, side = 1, cex = .7, line = .5)
      }
      write.csv(do.call(rbind,index), file.path(output_dir,"multichain_traces_index.csv"),row.names=FALSE)
    }, finally = dev.off())
  }
  rm(values, series); invisible(gc())

  # Keep chain identities separately. Pooling is for exploratory posterior summaries;
  # it does not assert that independent chains have mixed or reached convergence.
  pooled <- chains[[1]]
  pooled$draws <- unlist(lapply(chains, `[[`, "draws"), recursive = FALSE)
  pooled$chain_id <- rep(seq_len(n_chains), each = M)
  omega_summary <- NULL
  if (isTRUE(pooled$persistence)) {
    omega <- lapply(chains, function(f)
      vapply(f$draws, function(d) d$transition$omega, numeric(1)))
    summarize_omega <- function(x, label) data.frame(chain = label,
      mean = mean(x), median = median(x), lower = unname(quantile(x, .025)),
      upper = unname(quantile(x, .975)))
    omega_summary <- do.call(rbind, c(lapply(seq_along(omega), function(k)
      summarize_omega(omega[[k]], paste0("chain_", k))),
      list(summarize_omega(unlist(omega), "pooled_exploratory"))))
    write.csv(omega_summary, file.path(output_dir, "omega_summary.csv"), row.names = FALSE)
    png(file.path(output_dir, "omega.png"), width = 1500, height = 650, res = 130)
    tryCatch({
      par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
      boxplot(omega, names = paste("Chain", seq_len(n_chains)), ylim = c(0, 1),
        col = c("#244b78", "#c04b39", "#009988"),
        main = "Forced persistence by chain", ylab = "Omega")
      hist(unlist(omega), breaks = 45, probability = TRUE, xlim = c(0, 1),
        col = "#c9dced", border = "white", xlab = "Omega",
        main = "Pooled posterior (exploratory)")
      curve(dbeta(x, settings$omega_prior[1], settings$omega_prior[2]),
        add = TRUE, col = "#bb4422", lwd = 2)
      legend("topright", "Prior density", col = "#bb4422", lty = 1, bty = "n")
    }, finally = dev.off())
  }
  post <- state_parameters <- covariate_effects <- NULL
  post_settings <- partition_cache <- NULL
  post_warnings <- character()
  if (post_processing) {
    message("Computing pooled SALSO and CHIPS on all retained allocations.")
    partition_cache <- prepare_partition_cache(pooled)
    post_settings <- list(n_runs = post_n_runs, credible = post_credible,
      cores = post_cores, method = "direct", seed = 404L)
    post_file <- file.path(output_dir, sprintf("postprocess_direct_runs%d_prob%s.rds",
      post_n_runs, format(post_credible, digits = 15, trim = TRUE)))
    post_warnings <- character()
    if (file.exists(post_file)) {
      cached <- readRDS(post_file)
      expected <- settings
      if (is.null(cached$settings$persistence) && !isTRUE(settings$persistence))
        expected$persistence <- NULL
      if (is.null(cached$settings$priors) && identical(settings$priors, irtihmm_priors()))
        expected$priors <- NULL
      if (!identical(cached$settings, expected) || !identical(cached$post_settings, post_settings))
        stop("Post-processing settings mismatch.")
      post <- cached$post
      if (!is.null(cached$post_warnings)) post_warnings <- cached$post_warnings
    } else {
      post <- withCallingHandlers(postprocess_irtihmm(pooled, credible = post_credible,
        n_runs = post_n_runs, cores = post_cores, seed = 404, partition_cache = partition_cache),
        warning = function(w) post_warnings <<- unique(c(post_warnings, conditionMessage(w))))
      saveRDS(list(settings = settings, post_settings = post_settings, post = post,
        post_warnings = post_warnings), post_file)
    }
    message("SALSO/CHIPS finished. Exporting partitions and state parameters.")
    plot_postprocess(post, file.path(output_dir, "post_processing.png"))
    state_parameters <- summarize_state_parameters(pooled, post, partition_cache = partition_cache)
    plot_state_parameters(state_parameters, file.path(output_dir, "state_parameters.png"))
    plot_state_parameters(state_parameters, file.path(output_dir, "state_parameters_absolute.png"), effects = "absolute")
    write.csv(state_parameters$summary, file.path(output_dir, "state_parameters.csv"), row.names = FALSE)
  } else message("SALSO/CHIPS OFF: skipping partition packing, searches and aligned summaries.")
  message("Summarizing item and respondent effects.")
  effects <- summarize_effects(pooled)
  plot_effects(effects, "items", file.path(output_dir, "item_effects.png"))
  plot_effects(effects, "subjects", file.path(output_dir, "subject_effects.png"))
  if (post_processing) {
    message("Summarizing covariate effects.")
    covariate_effects <- summarize_covariate_effects(pooled, post, partition_cache = partition_cache)
    plot_covariate_effects(covariate_effects, "coefficients", file.path(output_dir, "covariate_coefficients.png"))
    plot_covariate_effects(covariate_effects, "probabilities", file.path(output_dir, "covariate_probabilities.png"))
  }
  write.csv(effects$items, file.path(output_dir, "item_effects.csv"), row.names = FALSE)
  write.csv(effects$subjects, file.path(output_dir, "subject_effects.csv"), row.names = FALSE)
  saveRDS(list(settings = settings, post_processing = post_processing,
    post_settings = post_settings, post = post, state_parameters = state_parameters,
    effects = effects, covariate_effects = covariate_effects, omega_summary = omega_summary,
    chain_id = pooled$chain_id),
    file.path(output_dir, if (post_processing) "posterior_summaries.rds" else "fit_summaries.rds"))
  for (k in seq_along(chains)) {
    message("Exporting plots for chain ", k, if (traceplots) " (with parameter traces)." else " (parameter traces OFF).")
    path <- file.path(output_dir, paste0("chain_", k))
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    plot_fit(chains[[k]], file = file.path(path, "ex_post.png"))
    # Use the same pooled core in every chain, with chain-specific compatibility.
    if (traceplots) {
      aligned <- NULL
      if (post_processing) {
        chain_post <- post; chain_post$chips$probability <- NULL
        chain_cache <- prepare_partition_cache(chains[[k]], parent = partition_cache,
          rows = ((k - 1L) * M + 1L):(k * M))
        aligned <- summarize_state_parameters(chains[[k]], chain_post, partition_cache = chain_cache)
        rm(chain_cache)
      }
      plot_parameter_traces(chains[[k]], file.path(path, "traces", "parameters.pdf"), aligned)
    }
  }
  rm(partition_cache)
  png(file.path(output_dir, "observed_data.png"), width = 1600, height = 1000, res = 140)
  tryCatch({
    par(mfrow = c(2,2), mar = c(4,4,3,1))
    data <- pooled$data; order <- order(data$z[,1])
    .item_heatmap(data$x, order, c("#f2f4f8", "#244b78"), "Observed accuracy by actual trial position", c(0,1))
    .item_heatmap(data$y, order, hcl.colors(64,"YlOrRd",rev=TRUE), "Observed log-times by actual trial position")
    for (j in 1:2) hist(data$z[,j], main = paste("Covariate",colnames(data$z)[j]),
      xlab = "Standardized value", col = "#648eae", border = "white")
  }, finally = dev.off())
  bad_rhat <- sum(diagnostic_table$rhat > 1.01, na.rm = TRUE)
  m_alpha_prior <- if (is.null(settings$priors)) irtihmm_priors()$m_alpha else settings$priors$m_alpha
  report <- c("ROAR XYZ: three-chain iHMM analysis", "",
    sprintf("Respondents: %d; items: %d. Excluded subject IDs: %s.",
      pooled$data$P, pooled$data$I, paste(settings$excluded_subjects, collapse=", ")),
    "Each chain: 5000 iterations, 1000 burn-in, thin=1; 4000 retained draws.",
    "Three chains with seeds 101, 202, 303; full model, copula and Gustafson ON.",
    sprintf("Diagonal mean prior: m_alpha ~ N(%g, %g^2); second value is SD.",
      m_alpha_prior$mean, m_alpha_prior$sd),
    if (isTRUE(pooled$persistence)) sprintf("Forced persistence ON; omega ~ Beta(%g, %g).",
      settings$omega_prior[1], settings$omega_prior[2]) else "Forced persistence OFF.",
    "Respondent-specific item order is retained; item effects use the supplied item map.",
    settings$time_policy, settings$transitions,
    "Covariates use the preparation scaling based on all 85 source respondents.",
    "No simulated truth is available; truth-only panels are replaced or omitted.",
    "multichain_diagnostics.csv: posterior rank-normalized split R-hat and bulk/tail ESS.",
    sprintf("Coordinates with finite availability in all chains: %d; R-hat > 1.01: %d.",
      nrow(diagnostic_table), bad_rhat),
    "Missing state/stick coordinates are never padded for R-hat/ESS.",
    "Raw ordered state coordinates are not CHIPS-relabeled; interpret with occupancy and invariant contrasts.",
    if (post_processing) c(
    sprintf("Post-processing: direct CHIPS threshold %.3f; %d search restarts; all retained allocations used.",
      post_credible, post_n_runs),
    "No complete probability/core-size path is computed; the plot shows the selected core.",
    "Pooled SALSO/CHIPS summaries are exploratory when chain diagnostics indicate poor mixing.",
    sprintf("Pooled SALSO clusters: %d; CHIPS core probability %.4f; core fraction %.4f.",
      length(unique(as.integer(post$states))), post$chips$probability, post$chips$core_fraction),
    "Conditional core intervals do not establish convergence or resolve all state/person confounding.") else
      c("SALSO/CHIPS post-processing: OFF. No partition cache or aligned state/covariate summaries computed.",
        "Any pre-existing post-processing files are preserved and were not refreshed by this run."),
    sprintf("Per-chain parameter PDFs: %s; multichain parameter PDF: %s.",
      if (traceplots) "ON" else "OFF", if (multichain_traceplots) "ON" else "OFF"),
    if (length(post_warnings)) paste("Post-processing warning:", post_warnings))
  writeLines(report, file.path(output_dir, "fit_report.txt"))
  message(paste(report, collapse = "\n"))
  invisible(list(post = post, diagnostics = diagnostic_table))
}

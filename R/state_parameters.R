#' Summarize state emission parameters conditional on a CHIPS core
#' @param fit An irtihmm_fit object with retained emission draws.
#' @param postprocess Output from postprocess_irtihmm with CHIPS enabled.
#' @param level Probability of each equal-tailed interval conditional on the core.
#' @param reference Core state label used for accuracy and log-time contrasts.
#'   NULL selects the largest core cluster (first label breaks ties).
#' @param partition_cache Optional reusable cache from prepare_partition_cache.
#' @param simulation Optional generating simulation for truth comparisons. A core
#'   cluster is mapped to truth only when all its responses share one true state.
#' @return An irtihmm_state_parameters object with matched parameter draws,
#'   a summary data frame, core-cluster mapping, compatible draw indices, core
#'   probability, and the joint core-and-interval probability for each interval.
#' @details In each compatible draw, every core cluster has a unique containing
#'   state. Its emission parameters are taken from the existing joint posterior
#'   draw, including information from responses outside the core. No model is
#'   refitted and no numeric correspondence of labels across draws is assumed.
#'   Absolute a and mu retain the original location-prior dependence; contrasts
#'   relative to the reference core cluster remove additive location ambiguity.
#'   Intervals are marginal per parameter and state, not simultaneous across all
#'   parameters. Their joint probability with the core is level * core_probability
#'   in the target distribution, subject to Monte Carlo quantile error.
#' @export
summarize_state_parameters <- function(fit, postprocess, level = .95,
                                       reference = NULL, simulation = NULL, partition_cache = NULL) {
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) stop("level must lie in (0, 1).")
  alignment <- .core_alignment(fit, postprocess, partition_cache)
  core <- postprocess$chips$subpartition
  labels <- alignment$labels; included <- alignment$included
  indices <- alignment$indices; probability <- alignment$probability
  ids <- alignment$ids; anchors <- alignment$anchors; sizes <- alignment$sizes
  if (is.null(reference)) reference <- ids[which.max(sizes)]
  if (length(reference) != 1L || is.na(reference) || !reference %in% ids)
    stop("reference must identify a state in the CHIPS core.")
  matched <- alignment$matched
  parameters <- setNames(lapply(c("a", "mu", "sigma", "rho"), function(name)
    matrix(NA_real_, length(indices), length(ids), dimnames = list(NULL, as.character(ids)))),
    c("a", "mu", "sigma", "rho"))
  for (j in seq_along(indices)) {
    e <- fit$draws[[indices[j]]]$emission
    states <- matched[j, ]
    if (!all(c("a", "mu", "log_sigma", "rho") %in% names(e)) ||
        any(lengths(e[c("a", "mu", "log_sigma", "rho")]) < max(states)))
      stop("Compatible draws must contain emission parameters for every matched state.")
    parameters$a[j, ] <- e$a[states]
    parameters$mu[j, ] <- e$mu[states]
    parameters$sigma[j, ] <- exp(e$log_sigma[states])
    parameters$rho[j, ] <- e$rho[states]
  }
  if (any(!is.finite(unlist(parameters))) || any(parameters$sigma <= 0) ||
      any(abs(parameters$rho) >= 1)) stop("Invalid matched emission parameters.")
  ref <- match(reference, ids)
  parameters$a_contrast <- sweep(parameters$a, 1L, parameters$a[, ref], "-")
  parameters$mu_contrast <- sweep(parameters$mu, 1L, parameters$mu[, ref], "-")
  mapping <- data.frame(state = ids, core_size = sizes,
    respondent = (anchors - 1L) %/% fit$data$I + 1L,
    item = (anchors - 1L) %% fit$data$I + 1L, true_state = NA_integer_)
  truth <- setNames(lapply(parameters, function(x) rep(NA_real_, length(ids))), names(parameters))
  if (!is.null(simulation)) {
    if (!identical(dim(simulation$states), dim(core)) ||
        !isTRUE(all.equal(simulation$data$x, fit$data$x)) ||
        !isTRUE(all.equal(simulation$data$y, fit$data$y)) ||
        !isTRUE(all.equal(simulation$data$z, fit$data$z)))
      stop("simulation must describe the same data and covariates used by fit.")
    true_labels <- as.integer(t(simulation$states))
    mapping$true_state <- vapply(ids, function(s) {
      values <- unique(true_labels[labels == s])
      if (length(values) == 1L) values else NA_integer_
    }, integer(1))
    e <- simulation$emission
    truth$a <- e$a[mapping$true_state]
    truth$mu <- e$mu[mapping$true_state]
    truth$sigma <- exp(e$log_sigma[mapping$true_state])
    truth$rho <- e$rho[mapping$true_state]
    truth$a_contrast <- truth$a - truth$a[ref]
    truth$mu_contrast <- truth$mu - truth$mu[ref]
  }
  tail <- (1 - level) / 2
  summary <- do.call(rbind, lapply(names(parameters), function(name) {
    values <- parameters[[name]]
    bounds <- vapply(seq_along(ids), function(k)
      quantile(values[, k], c(tail, .5, 1 - tail), names = FALSE), numeric(3))
    data.frame(parameter = name, state = ids, mean = colMeans(values),
      median = bounds[2, ], lower = bounds[1, ], upper = bounds[3, ], truth = truth[[name]],
      row.names = NULL)
  }))
  structure(list(summary = summary, draws = parameters, mapping = mapping,
    matched_states = matched, compatible_draws = indices, reference = reference,
    level = level, core_probability = probability,
    joint_probability = level * probability), class = "irtihmm_state_parameters")
}

#' Plot CHIPS-conditional state emission parameter intervals
#' @param parameters Output from summarize_state_parameters.
#' @param file Optional PNG filename. NULL plots on the current graphics device.
#' @param effects Plot identifiable "contrasts" (default) or "absolute" a and mu.
#' @return Invisibly, parameters. Points show posterior medians and bars show
#'   marginal equal-tailed intervals conditional on the CHIPS core. Red crosses
#'   show generating values only for unambiguously mapped true states. Absolute
#'   effects depend on the location priors. State IDs refer to CHIPS core clusters.
#' @export
plot_state_parameters <- function(parameters, file = NULL,
                                  effects = c("contrasts", "absolute")) {
  if (!inherits(parameters, "irtihmm_state_parameters"))
    stop("Expected output from summarize_state_parameters.")
  effects <- match.arg(effects)
  x <- parameters
  names <- if (effects == "contrasts") c("a_contrast", "mu_contrast", "sigma", "rho") else
    c("a", "mu", "sigma", "rho")
  titles <- if (effects == "contrasts") c("State accuracy contrasts", "State log-time contrasts") else
    c("Absolute state accuracy effects", "Absolute state log-time effects")
  titles <- c(titles, "State log-time standard deviations", "State copula correlations")
  ylabels <- if (effects == "contrasts") c("a[s] - a[reference]", "mu[s] - mu[reference]") else
    c("a[s] (location-prior dependent)", "mu[s] (location-prior dependent)")
  ylabels <- c(ylabels, "sigma[s]", "rho[s]")
  device <- .plot_start(file, rows = 2, width = 1400, height = 1050)
  on.exit(.plot_end(device), add = TRUE)
  graphics::par(oma = c(2, 0, 3, 0))
  K <- nrow(x$mapping)
  palette <- grDevices::hcl.colors(max(K, 2), "Dark 3")[seq_len(K)]
  for (j in seq_along(names)) {
    values <- x$summary[x$summary$parameter == names[j], ]
    limits <- range(values$lower, values$upper, values$truth, finite = TRUE)
    if (diff(limits) == 0) limits <- limits + c(-.1, .1)
    limits <- limits + c(-.08, .12) * diff(limits)
    if (names[j] == "rho") limits <- c(-1, 1)
    if (names[j] == "sigma") limits[1] <- max(0, limits[1])
    graphics::plot(seq_len(K), values$median, type = "n", xlim = c(.5, K + .5),
      ylim = limits, xaxt = "n", xlab = "CHIPS core state", ylab = ylabels[j], main = titles[j])
    graphics::axis(1, at = seq_len(K), labels = values$state)
    if (names[j] != "sigma") graphics::abline(h = 0, lty = 3, col = "#bbbbbb")
    graphics::segments(seq_len(K), values$lower, seq_len(K), values$upper,
      col = palette, lwd = 3)
    graphics::points(seq_len(K), values$median, pch = 16, col = palette)
    valid <- which(is.finite(values$truth))
    if (length(valid)) graphics::points(valid, values$truth[valid], pch = 4, col = "#bb4422", lwd = 2)
    graphics::legend("topright", legend = c("Posterior median / interval", if (length(valid)) "Generating value"),
      col = c("#444444", if (length(valid)) "#bb4422"),
      pch = c(16, if (length(valid)) 4), bty = "n", cex = .7)
  }
  graphics::mtext(sprintf("CHIPS: core probability %.3f; %.1f%% conditional intervals; %d compatible draws",
    x$core_probability, 100 * x$level, length(x$compatible_draws)), side = 3, outer = TRUE, line = 1)
  graphics::mtext(sprintf("Reference core state: %s. Nominal joint core + interval probability: %.3f. Truth omitted for mixed cores.",
    x$reference, x$joint_probability), side = 1, outer = TRUE, line = .5, cex = .8)
  invisible(x)
}

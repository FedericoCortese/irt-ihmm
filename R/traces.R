# Build trace series without truncating to the common stored state prefix.
# Uninstantiated parameter coordinates stay NA; their lines must have gaps.
.parameter_trace_series <- function(fit, state_parameters = NULL) {
  if (!inherits(fit, "irtihmm_fit") || !length(fit$draws))
    stop("Expected a fitted model with retained draws.")
  draws <- fit$draws; M <- length(draws)
  iterations <- vapply(draws, `[[`, integer(1), "iteration")
  K <- max(vapply(draws, function(d) length(d$emission$a), integer(1)))
  L <- max(vapply(draws, function(d) ncol(d$transition$beta), integer(1)))
  R <- max(vapply(draws, function(d) nrow(d$transition$alpha), integer(1)))
  A <- max(vapply(draws, function(d) ncol(d$transition$alpha), integer(1)))
  series <- list()
  add <- function(family, parameter, getter = NULL, values = NULL,
                  kind = "parameter", scope = "original coordinates") {
    if (is.null(values)) values <- vapply(draws, function(d) {
      value <- getter(d)
      if (length(value) == 0L) NA_real_ else as.numeric(value)
    }, numeric(1))
    if (length(values) != M || any(!is.finite(values[!is.na(values)])))
      stop("Invalid trace values for ", parameter, ".")
    series[[length(series) + 1L]] <<- list(family = family, parameter = parameter,
      values = values, kind = kind, scope = scope)
  }
  for (name in c("b", "nu", "d", "eta")) {
    item <- name %in% c("b", "nu")
    n <- if (item) fit$data$I else fit$data$P
    family <- if (item) "Item effects" else "Respondent effects"
    for (j in seq_len(n)) {
      add(paste(family, "- absolute"), sprintf("%s[%d]", name, j),
          function(d) d$emission[[name]][j])
      add(paste(family, "- centered"), sprintf("%s[%d] - mean(%s)", name, j, name),
          function(d) d$emission[[name]][j] - mean(d$emission[[name]]),
          kind = "derived", scope = "centered within every draw")
    }
  }
  for (name in c("a", "mu", "log_sigma", "sigma", "rho")) for (k in seq_len(K)) {
    add("State emissions - ordered indices", sprintf("%s[state=%d]", name, k),
      function(d) if (name == "sigma") exp(d$emission$log_sigma[k]) else d$emission[[name]][k],
      kind = if (name == "sigma") "derived" else "parameter",
      scope = "ordered state index; NA when not instantiated")
  }
  for (k in seq_len(K)) add("State occupancy - ordered indices", sprintf("count[state=%d]", k),
    function(d) sum(d$states == k), kind = "derived", scope = "ordered state index")
  q <- ncol(fit$data$z)
  covariates <- colnames(fit$data$z)
  if (is.null(covariates)) covariates <- paste0("z", seq_len(q))
  for (k in seq_len(L)) for (j in seq_len(q)) {
    add("Covariate coefficients - ordered sticks", sprintf("beta[%s,stick=%d]", covariates[j], k),
      function(d) if (k <= ncol(d$transition$beta)) d$transition$beta[j, k] else NA_real_,
      scope = "ordered stick index; NA when not instantiated")
  }
  for (k in seq_len(A)) {
    add("Initial distribution intercepts", sprintf("alpha0[stick=%d]", k),
        function(d) d$transition$alpha0[k], scope = "ordered stick index; NA when not instantiated")
    for (r in seq_len(R)) add("Transition intercepts - ordered indices",
      sprintf("alpha[origin=%d,stick=%d]", r, k),
      function(d) if (r <= nrow(d$transition$alpha) && k <= ncol(d$transition$alpha))
        d$transition$alpha[r, k] else NA_real_,
      scope = "ordered origin/stick indices; NA when not instantiated")
  }
  for (margin in c("accuracy", "time")) {
    add("Hyperparameters", paste0(margin, ".item_mean"), function(d) d$hyper[[margin]]$mean)
    for (j in 1:3) add("Hyperparameters",
      paste0(margin, ".variance_", c("item", "state", "person")[j]),
      function(d) d$hyper[[margin]]$variances[j])
    if (any(vapply(draws, function(d) !is.null(d$locations[[margin]]), logical(1)))) {
      for (j in 1:2) add("Recovered locations", sprintf("lambda_%s[%d]", margin, j),
        function(d) d$locations[[margin]]$lambda[j], scope = "Gaussian conditional recovery")
    }
  }
  for (name in c("mean", "variance")) add("Hyperparameters", paste0("log_sigma.", name),
    function(d) d$hyper$scale[[name]])
  for (name in c("sigma_alpha2", "mean_diagonal", "var_diagonal"))
    add("Hyperparameters", paste0("transition.", name), function(d) d$transition[[name]])
  if (isTRUE(fit$persistence))
    add("Persistence", "omega", function(d) d$transition$omega,
        scope = "global probability of forced persistence")
  if (!is.null(state_parameters)) {
    if (!inherits(state_parameters, "irtihmm_state_parameters"))
      stop("state_parameters must come from summarize_state_parameters.")
    indices <- state_parameters$compatible_draws
    if (!length(indices) || anyNA(indices) || any(indices < 1L | indices > M) ||
        anyDuplicated(indices)) stop("Invalid CHIPS-compatible draw indices.")
    for (name in names(state_parameters$draws)) {
      values <- state_parameters$draws[[name]]
      if (nrow(values) != length(indices)) stop("CHIPS draw dimensions do not match.")
      for (k in seq_len(ncol(values))) {
        full <- rep(NA_real_, M); full[indices] <- values[, k]
        add("State emissions - CHIPS conditional",
          sprintf("%s[core_state=%s]", name, colnames(values)[k]), values = full,
          kind = if (name %in% c("sigma", "a_contrast", "mu_contrast")) "derived" else "parameter",
          scope = "CHIPS conditional; gaps at incompatible draws, not a full-chain diagnostic")
      }
    }
  }
  list(iterations = iterations, series = series)
}

#' Export paginated traceplots of all saved model parameters
#' @param fit An irtihmm_fit object with retained draws.
#' @param file Output PDF path. A companion _index.csv lists all plotted coordinates.
#' @param state_parameters Optional output from summarize_state_parameters for
#'   this fit, adding conditional CHIPS-aligned emission traces to the raw traces.
#' @param panels Number of rows and columns per page; defaults to 3 by 3.
#' @return Invisibly, the trace index with family, parameter, scope, page, panel,
#'   number of present draws and first/last saved iterations where it is present.
#' @details Every saved ordered state/stick coordinate is included up to the
#'   largest instantiated prefix, including rare coordinates. Missing coordinates
#'   produce gaps, never zero or prior padding. Burn-in was discarded by the
#'   sampler and cannot be reconstructed. Item/person centered traces and sigma
#'   are included alongside original parameters. Fixed prior constants and latent
#'   allocations are not model-parameter traces; ordered state counts are supplied.
#'   CHIPS traces condition on a core and cannot establish full-chain convergence.
#' @export
plot_parameter_traces <- function(fit, file, state_parameters = NULL, panels = c(3L, 3L)) {
  if (!is.character(file) || length(file) != 1L || is.na(file) ||
      !grepl("\\.pdf$", file, ignore.case = TRUE)) stop("file must be a PDF path.")
  if (!is.numeric(panels) || length(panels) != 2L || any(!is.finite(panels)) ||
      any(panels < 1 | panels != floor(panels))) stop("panels must contain two positive integers.")
  traces <- .parameter_trace_series(fit, state_parameters)
  n <- length(traces$iterations); index <- list(); page <- 0L
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  # Cairo embeds fonts and preserves axis/text clipping on dense multipage traces.
  if (isTRUE(capabilities("cairo"))) {
    grDevices::cairo_pdf(file, width = 14, height = 10, onefile = TRUE)
  } else {
    grDevices::pdf(file, width = 14, height = 10, onefile = TRUE, compress = TRUE)
  }
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::par(mfrow = panels, mar = c(3, 4, 2.8, .8), oma = c(3, 0, 3, 0),
                mgp = c(1.8, .5, 0), cex = .85, las = 1)
  families <- vapply(traces$series, `[[`, character(1), "family")
  capacity <- prod(panels)
  for (family in unique(families)) {
    coordinates <- which(families == family)
    chunks <- split(coordinates, ceiling(seq_along(coordinates) / capacity))
    for (chunk in chunks) {
      page <- page + 1L
      for (slot in seq_len(capacity)) {
        if (slot > length(chunk)) { graphics::plot.new(); next }
        trace <- traces$series[[chunk[slot]]]; values <- trace$values
        present <- !is.na(values); count <- sum(present)
        bounds <- if (count) range(values[present]) else c(-1, 1)
        if (diff(bounds) == 0) bounds <- bounds + c(-1, 1) * .05 * max(1, abs(bounds))
        bounds <- bounds + c(-.04, .04) * diff(bounds)
        graphics::plot(traces$iterations, values, type = "l", col = "#244b78", lwd = .6,
          xlim = range(traces$iterations) + if (n == 1L) c(-.5, .5) else c(0, 0),
          ylim = bounds, main = trace$parameter, cex.main = .9,
          xlab = "Iteration (saved draws)", ylab = "Value")
        isolated <- present & !c(FALSE, head(present, -1)) & !c(tail(present, -1), FALSE)
        if (any(isolated)) graphics::points(traces$iterations[isolated], values[isolated],
                                            pch = 16, cex = .4, col = "#244b78")
        graphics::mtext(sprintf("Present: %d / %d", count, n), side = 3, adj = 1, line = .1, cex = .55)
        index[[length(index) + 1L]] <- data.frame(family = family, parameter = trace$parameter,
          kind = trace$kind, scope = trace$scope, page = page, panel = slot,
          present_draws = count, saved_draws = n,
          first_present_iteration = if (count) min(traces$iterations[present]) else NA_integer_,
          last_present_iteration = if (count) max(traces$iterations[present]) else NA_integer_)
      }
      graphics::mtext(family, side = 3, outer = TRUE, line = 1, font = 2, cex = 1.1)
      graphics::mtext(sprintf("Retained draws only | Gustafson %s | Page %d",
        if (isFALSE(fit$gustafson)) "OFF" else "ON", page), side = 1, outer = TRUE, line = .5, cex = .8)
      note <- if (family == "State emissions - CHIPS conditional")
        "CHIPS-conditional traces: incompatible draws are gaps. Not a full-chain convergence diagnostic." else
        "Raw state/stick indices retain their order; absent parameters are gaps. Burn-in is not stored."
      graphics::mtext(note, side = 1, outer = TRUE, line = 1.6, cex = .75)
    }
  }
  index <- do.call(rbind, index)
  utils::write.csv(index, sub("\\.pdf$", "_index.csv", file, ignore.case = TRUE), row.names = FALSE)
  invisible(index)
}

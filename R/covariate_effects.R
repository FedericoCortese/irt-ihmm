# Evaluate probabilities of selected destinations, preserving all preceding
# sticks, the finite terminal remainder, and the infinite unallocated tail.
.covariate_state_probs <- function(transition, z, states) {
  tr <- transition
  eta <- rbind(as.numeric(tr$alpha0), tr$alpha[states, , drop = FALSE])
  eta <- sweep(eta, 2L, as.numeric(z %*% tr$beta), "+")
  survive <- pnorm(-eta, log.p = TRUE)
  out <- matrix(0, length(states) + 1L, length(states))
  for (j in seq_along(states)) {
    k <- states[j]
    before <- if (k == 1L) rep(0, nrow(eta)) else
      rowSums(survive[, seq_len(k - 1L), drop = FALSE])
    out[, j] <- if (tr$mode == "finite" && k == tr$K) exp(before) else
      exp(before + pnorm(eta[, k], log.p = TRUE))
  }
  if (!is.null(tr$omega) && tr$omega != 0) {
    out[-1L, ] <- (1 - tr$omega) * out[-1L, , drop = FALSE]
    out[cbind(seq_along(states) + 1L, seq_along(states))] <-
      diag(out[-1L, , drop = FALSE]) + tr$omega
  }
  out
}

#' Summarize covariate coefficients and CHIPS-aligned probability contrasts
#' @param fit An irtihmm_fit object with respondent covariates.
#' @param postprocess CHIPS post-processing for the same retained draws.
#' @param partition_cache Optional reusable cache from prepare_partition_cache.
#' @param baseline Numeric reference profile, ordered like columns of data$z.
#'   NULL uses zero for binary 0/1 covariates and the sample mean otherwise.
#' @param change Numeric changes, one per covariate, applied one at a time.
#'   NULL uses +1 for every covariate. Other covariates stay at baseline.
#' @param level Probability of each marginal equal-tailed conditional interval.
#' @return An irtihmm_covariate_effects object containing coefficient and
#'   probability-contrast summaries and draw arrays, profile settings, CHIPS state
#'   labels, selected draw indices and the core probability.
#' @details Coefficients are summarized at fixed ordered stick indices over the
#'   stored prefix common to every compatible draw. They are not relabeled by
#'   CHIPS. A positive beta increases its probit stick but need not increase the
#'   corresponding full destination probability, which also uses earlier sticks.
#'   Probability contrasts compare the changed and baseline profiles in each
#'   draw, then summarize changes after matching origins and destinations to core
#'   states. All preceding sticks are included and probabilities are not
#'   renormalized over the core. The first origin row is the initial distribution.
#'   Intervals condition on the whole core and are marginal per coefficient or
#'   probability contrast. No linear coefficient on an aligned core state is
#'   inferred. Covariates enter transitions, not the emission predictors directly.
#' @export
summarize_covariate_effects <- function(fit, postprocess, baseline = NULL,
                                        change = NULL, level = .95, partition_cache = NULL) {
  if (!inherits(fit, "irtihmm_fit") || !ncol(fit$data$z))
    stop("Expected a fitted model with respondent covariates.")
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) || level <= 0 || level >= 1)
    stop("level must lie in (0, 1).")
  alignment <- .core_alignment(fit, postprocess, partition_cache)
  indices <- alignment$indices
  states <- alignment$ids
  z <- fit$data$z; q <- ncol(z); K <- length(states); M <- length(indices)
  covariates <- colnames(z)
  if (is.null(covariates)) covariates <- paste0("z", seq_len(q))
  if (is.null(baseline)) baseline <- vapply(seq_len(q), function(j)
    if (all(z[, j] %in% c(0, 1))) 0 else mean(z[, j]), numeric(1))
  if (is.null(change)) change <- rep(1, q)
  for (value in list(baseline, change))
    if (!is.numeric(value) || length(value) != q || any(!is.finite(value)))
      stop("baseline and change must be finite numeric vectors with one entry per covariate.")
  if (any(change == 0)) stop("Each covariate change must be nonzero.")
  names(baseline) <- names(change) <- covariates
  transitions <- lapply(indices, function(m) fit$draws[[m]]$transition)
  for (tr in transitions) {
    if (!is.matrix(tr$beta) || nrow(tr$beta) != q || any(!is.finite(tr$beta)))
      stop("Compatible draws must contain finite transition coefficients for each covariate.")
  }
  L <- min(vapply(transitions, function(tr) ncol(tr$beta), integer(1)))
  beta <- array(NA_real_, c(M, q, L),
    dimnames = list(NULL, covariates, as.character(seq_len(L))))
  origins <- c("Initial", paste("State", states))
  destinations <- paste("State", states)
  delta <- array(NA_real_, c(M, K + 1L, K, q),
    dimnames = list(NULL, origins, destinations, covariates))
  for (m in seq_len(M)) {
    tr <- transitions[[m]]
    if (L) beta[m, , ] <- tr$beta[, seq_len(L), drop = FALSE]
    matched <- alignment$matched[m, ]
    base <- .covariate_state_probs(tr, baseline, matched)
    for (j in seq_len(q)) {
      profile <- baseline; profile[j] <- profile[j] + change[j]
      delta[m, , , j] <- .covariate_state_probs(tr, profile, matched) - base
    }
  }
  summarize <- function(x) {
    bounds <- quantile(x, c((1 - level) / 2, .5, (1 + level) / 2), names = FALSE)
    c(mean = mean(x), median = bounds[2], lower = bounds[1], upper = bounds[3],
      probability_positive = mean(x > 0))
  }
  coefficients <- data.frame()
  if (L) {
    coefficients <- expand.grid(covariate = covariates, stick = seq_len(L),
                                 stringsAsFactors = FALSE)
    coefficients <- cbind(coefficients, as.data.frame(t(vapply(seq_len(nrow(coefficients)),
      function(i) summarize(beta[, match(coefficients$covariate[i], covariates),
                                  coefficients$stick[i]]), numeric(5)))))
  }
  probabilities <- expand.grid(origin = origins, destination = destinations,
    covariate = covariates, stringsAsFactors = FALSE)
  probabilities <- cbind(probabilities, as.data.frame(t(vapply(seq_len(nrow(probabilities)),
    function(i) summarize(delta[, match(probabilities$origin[i], origins),
      match(probabilities$destination[i], destinations),
      match(probabilities$covariate[i], covariates)]), numeric(5)))))
  structure(list(coefficients = coefficients, probabilities = probabilities,
    beta_draws = beta, probability_draws = delta, common_sticks = L,
    covariates = covariates, origins = origins, destinations = destinations,
    baseline = baseline, change = change, level = level, draw_indices = indices,
    core_probability = alignment$probability), class = "irtihmm_covariate_effects")
}

#' Plot ordered-stick coefficients or CHIPS-aligned covariate probability effects
#' @param effects Output from summarize_covariate_effects.
#' @param type Either "coefficients" or "probabilities".
#' @param file Optional PNG filename. NULL plots on the current graphics device.
#' @return Invisibly, effects. Coefficient plots retain ordered stick indices.
#'   Probability heatmaps show changes (changed profile minus baseline) and
#'   marginal equal-tailed intervals, with origins/destinations aligned by CHIPS.
#' @export
plot_covariate_effects <- function(effects, type = c("coefficients", "probabilities"),
                                   file = NULL) {
  if (!inherits(effects, "irtihmm_covariate_effects"))
    stop("Expected output from summarize_covariate_effects.")
  type <- match.arg(type)
  q <- length(effects$covariates)
  if (type == "coefficients" && effects$common_sticks == 0)
    stop("There are no stochastic sticks shared by all selected draws.")
  device <- .plot_start(file, rows = ceiling(q / 2),
                        width = if (type == "probabilities") 2000 else 1500,
                        height = 650 * ceiling(q / 2) + 160)
  on.exit(.plot_end(device), add = TRUE)
  graphics::par(oma = c(if (type == "probabilities") 5 else 2, 0, 3, 0),
                mar = c(4.5, if (type == "probabilities") 7 else 5, 3, 1),
                mfrow = c(ceiling(q / 2), min(q, 2)))
  color_limit <- max(.01, abs(effects$probabilities$median))
  for (name in effects$covariates) {
    if (type == "coefficients") {
      x <- effects$coefficients[effects$coefficients$covariate == name, ]
      limits <- range(x$lower, x$upper, 0)
      padding <- max(.1, .1 * diff(limits))
      graphics::plot(x$stick, x$median, type = "n", xaxt = "n",
        xlim = range(x$stick) + c(-.5, .5), ylim = limits + c(-padding, padding),
        xlab = "Ordered stick index (not a CHIPS state)", ylab = "Probit stick coefficient beta[j,k]",
        main = paste("Covariate:", name))
      graphics::axis(1, at = x$stick)
      graphics::abline(h = 0, lty = 3, col = "#999999")
      graphics::segments(x$stick, x$lower, x$stick, x$upper, col = "#244b78", lwd = 3)
      graphics::points(x$stick, x$median, col = "#244b78", pch = 16)
    } else {
      x <- effects$probabilities[effects$probabilities$covariate == name, ]
      nr <- length(effects$origins); nc <- length(effects$destinations)
      median <- matrix(x$median, nr, nc)[nr:1, , drop = FALSE]
      text <- matrix(sprintf("%+.2g\n[%+.2g, %+.2g]", x$median, x$lower, x$upper), nr, nc)
      text <- text[nr:1, , drop = FALSE]
      limit <- color_limit
      palette <- grDevices::colorRampPalette(c("#2166ac", "#ffffff", "#b2182b"))(101)
      graphics::image(seq(.5, nc + .5, length.out = nc + 1),
        seq(.5, nr + .5, length.out = nr + 1), t(median), zlim = c(-limit, limit),
        col = palette, axes = FALSE, xlab = "Destination core state", ylab = "",
        main = sprintf("%s: %.2f -> %.2f", name, effects$baseline[name],
          effects$baseline[name] + effects$change[name]))
      graphics::axis(1, at = seq_len(nc), labels = effects$destinations)
      graphics::axis(2, at = seq_len(nr), labels = rev(effects$origins), las = 1, cex.axis = .8)
      graphics::mtext("Origin", side = 2, line = 5.2, las = 0)
      graphics::text(rep(seq_len(nc), each = nr), rep(seq_len(nr), nc),
        labels = as.vector(text), cex = min(.85, 4 / max(nr, nc)),
        col = ifelse(abs(as.vector(median)) > .65 * limit, "white", "black"))
      graphics::box()
    }
  }
  graphics::mtext(sprintf("CHIPS core probability %.3f | %.1f%% conditional intervals | %d draws",
    effects$core_probability, effects$level * 100, length(effects$draw_indices)),
    side = 3, outer = TRUE, line = 1)
  if (type == "coefficients") {
    note <- sprintf("Fixed ordered sticks shared by every compatible draw: 1-%d. Coefficients are not relabeled by CHIPS.",
                      effects$common_sticks)
    graphics::mtext(note, side = 1, outer = TRUE, line = .5, cex = .8)
  } else {
    profile <- paste(sprintf("%s = %.2f", effects$covariates, effects$baseline), collapse = "; ")
    graphics::mtext(paste("Baseline:", profile), side = 1, outer = TRUE, line = .5, cex = .85)
    graphics::mtext("One covariate changes at a time. Cell: median change [interval]. Red = increase; blue = decrease.",
      side = 1, outer = TRUE, line = 1.8, cex = .85)
    graphics::mtext("Original probabilities, including earlier sticks; no renormalization over the core. Intervals are marginal.",
      side = 1, outer = TRUE, line = 3.1, cex = .85)
  }
  invisible(effects)
}

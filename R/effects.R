#' Summarize centered item and respondent effects
#' @param fit An irtihmm_fit object with retained draws.
#' @param postprocess Optional CHIPS post-processing. NULL uses all retained draws;
#'   otherwise summaries condition on agreement with the entire CHIPS core.
#' @param partition_cache Optional reusable cache from prepare_partition_cache.
#' @param simulation Optional generating simulation for truth comparisons.
#' @param level Probability of each marginal equal-tailed interval.
#' @return An irtihmm_effects object containing item and subject summary tables,
#'   centered draw matrices, selected draw indices, covariates, interval level,
#'   and the core probability when conditioning on CHIPS.
#' @details b, nu, d and eta are each centered across their respective items or
#'   respondents separately within every draw. This removes additive location
#'   shifts without changing the model or refitting it. Item and respondent IDs
#'   are fixed across draws, so state-label alignment is not required. Intervals
#'   are marginal, not simultaneous across effects. CHIPS-conditioned intervals
#'   have nominal joint probability level times core_probability with the core.
#' @export
summarize_effects <- function(fit, postprocess = NULL, simulation = NULL, level = .95,
                              partition_cache = NULL) {
  if (!inherits(fit, "irtihmm_fit") || !length(fit$draws))
    stop("Expected a fitted model with retained draws.")
  if (!is.numeric(level) || length(level) != 1L || !is.finite(level) ||
      level <= 0 || level >= 1) stop("level must lie in (0, 1).")
  indices <- seq_along(fit$draws)
  probability <- NULL
  if (!is.null(postprocess)) {
    alignment <- .core_alignment(fit, postprocess, partition_cache)
    indices <- alignment$indices; probability <- alignment$probability
  }
  if (!is.null(simulation) &&
      (!isTRUE(all.equal(simulation$data$x, fit$data$x)) ||
       !isTRUE(all.equal(simulation$data$y, fit$data$y)) ||
       !isTRUE(all.equal(simulation$data$z, fit$data$z))))
    stop("simulation must describe the same data and covariates used by fit.")
  sizes <- c(b = fit$data$I, nu = fit$data$I, d = fit$data$P, eta = fit$data$P)
  truth <- list()
  draws <- lapply(names(sizes), function(name) {
    values <- do.call(rbind, lapply(indices, function(j) {
      value <- fit$draws[[j]]$emission[[name]]
      if (!is.numeric(value) || length(value) != sizes[[name]] || any(!is.finite(value)))
        stop("Retained emission draws must contain finite b, nu, d and eta vectors matching the data.")
      value
    }))
    sweep(values, 1L, rowMeans(values), "-")
  })
  names(draws) <- names(sizes)
  for (name in names(sizes)) {
    truth[[name]] <- rep(NA_real_, sizes[[name]])
    if (!is.null(simulation)) {
      value <- simulation$emission[[name]]
      if (!is.numeric(value) || length(value) != sizes[[name]] || any(!is.finite(value)))
        stop("simulation must contain finite item and respondent effects matching the data.")
      truth[[name]] <- value - mean(value)
    }
  }
  summarize <- function(name, prefix) {
    values <- draws[[name]]
    tail <- (1 - level) / 2
    bounds <- vapply(seq_len(ncol(values)), function(k)
      quantile(values[, k], c(tail, .5, 1 - tail), names = FALSE), numeric(3))
    out <- data.frame(mean = colMeans(values), median = bounds[2, ],
      lower = bounds[1, ], upper = bounds[3, ], truth = truth[[name]], row.names = NULL)
    names(out) <- paste0(prefix, "_", names(out))
    out
  }
  labels <- function(x, n, prefix) if (is.null(x)) paste0(prefix, seq_len(n)) else x
  items <- data.frame(id = seq_len(fit$data$I),
    label = labels(if (is.null(fit$data$item_index)) colnames(fit$data$x) else
      colnames(fit$data$item_index), fit$data$I, "i"),
    summarize("b", "accuracy"), summarize("nu", "time"))
  subjects <- data.frame(id = seq_len(fit$data$P),
    label = labels(rownames(fit$data$x), fit$data$P, "p"),
    summarize("d", "accuracy"), summarize("eta", "time"))
  structure(list(items = items, subjects = subjects, draws = draws,
    draw_indices = indices, level = level, core_probability = probability,
    covariates = fit$data$z), class = "irtihmm_effects")
}

#' Plot item or respondent effects with marginal intervals and a joint scatterplot
#' @param effects Output from summarize_effects.
#' @param type Either "items" or "subjects".
#' @param file Optional PNG filename. NULL plots on the current graphics device.
#' @param covariate Column index or name used to color subjects, defaulting to the
#'   first covariate. NULL disables grouping. More than six distinct values are
#'   split at the median for display only. Ignored for items or absent covariates.
#' @return Invisibly, effects. The first two panels share an ordering by median
#'   accuracy effect. Scatterplot crosses show marginal intervals, not a joint
#'   credible region. Red crosses show centered generating values when available.
#' @export
plot_effects <- function(effects, type = c("items", "subjects"), file = NULL,
                         covariate = 1L) {
  if (!inherits(effects, "irtihmm_effects")) stop("Expected output from summarize_effects.")
  type <- match.arg(type)
  x <- effects[[type]]; n <- nrow(x)
  colors <- rep("#244b78", n)
  group <- NULL; group_name <- NULL
  if (type == "subjects" && !is.null(covariate) && ncol(effects$covariates)) {
    z <- effects$covariates
    index <- if (is.character(covariate)) match(covariate, colnames(z)) else covariate
    if (!is.numeric(index) || length(index) != 1L || is.na(index) ||
        index != floor(index) || index < 1 || index > ncol(z))
      stop("covariate must identify one column of the respondent covariates.")
    value <- z[, index]
    group <- if (length(unique(value)) > 6L)
      factor(ifelse(value <= median(value), "<= median", "> median")) else factor(value)
    group_name <- colnames(z)[index]
    if (is.null(group_name) || is.na(group_name) || !nzchar(group_name))
      group_name <- paste0("z[,", index, "]")
    palette <- grDevices::hcl.colors(max(nlevels(group), 2), "Dark 3")[seq_len(nlevels(group))]
    colors <- palette[as.integer(group)]
  }
  device <- .plot_start(file, rows = 1, width = 1900, height = max(950, 500 + 15 * n))
  on.exit(.plot_end(device), add = TRUE)
  graphics::par(mfrow = c(1, 3), oma = c(3, 0, 3, 0),
                mar = c(5.5, min(9, max(4, .35 * max(nchar(x$label)) + 1.5)), 3, 1))
  order <- order(x$accuracy_median, x$id)
  limits <- function(...) {
    value <- range(..., finite = TRUE)
    pad <- if (diff(value) == 0) .1 else .12 * diff(value)
    value + c(-pad, pad)
  }
  has_truth <- any(is.finite(x$accuracy_truth))
  first_label <- if (type == "items") "Centered easiness: b[i] - mean(b)" else
    "Centered ability: d[p] - mean(d)"
  second_label <- if (type == "items") "Centered time intensity: nu[i] - mean(nu)" else
    "Centered log-time: eta[p] - mean(eta)"
  titles <- if (type == "items") c("Item easiness", "Item time intensity", "Easiness and time intensity") else
    c("Respondent ability", "Respondent log-time", "Ability and log-time")
  legend <- function() graphics::legend("topleft",
    legend = c("Posterior median / interval", if (has_truth) "Generating value"),
    col = c("#444444", if (has_truth) "#bb4422"),
    pch = c(16, if (has_truth) 4), bty = "n", cex = .72)
  for (j in 1:2) {
    prefix <- c("accuracy", "time")[j]
    lower <- x[[paste0(prefix, "_lower")]]
    upper <- x[[paste0(prefix, "_upper")]]
    middle <- x[[paste0(prefix, "_median")]]
    truth <- x[[paste0(prefix, "_truth")]]
    graphics::plot(middle[order], seq_len(n), type = "n", yaxt = "n",
      xlim = limits(lower, upper, truth), ylim = c(.5, n + 2.5),
      xlab = c(first_label, second_label)[j], ylab = "", main = titles[j])
    graphics::abline(v = 0, lty = 3, col = "#aaaaaa")
    graphics::abline(h = seq_len(n), col = "#eeeeee")
    # The image height already grows with n; keep row labels readable.
    graphics::axis(2, at = seq_len(n), labels = x$label[order], cex.axis = .7)
    graphics::segments(lower[order], seq_len(n), upper[order], seq_len(n),
      col = colors[order], lwd = 2)
    graphics::points(middle[order], seq_len(n), pch = 16, col = colors[order], cex = .8)
    if (has_truth) graphics::points(truth[order], seq_len(n), pch = 4, col = "#bb4422", lwd = 1.5)
    legend()
  }
  graphics::plot(x$accuracy_median, x$time_median, type = "n",
    xlim = limits(x$accuracy_lower, x$accuracy_upper, x$accuracy_truth),
    ylim = limits(x$time_lower, x$time_upper, x$time_truth),
    xlab = first_label, ylab = second_label, main = titles[3])
  graphics::abline(h = 0, v = 0, lty = 3, col = "#aaaaaa")
  faint <- grDevices::adjustcolor(colors, alpha.f = .2)
  graphics::segments(x$accuracy_lower, x$time_median, x$accuracy_upper,
                      x$time_median, col = faint)
  graphics::segments(x$accuracy_median, x$time_lower, x$accuracy_median,
                      x$time_upper, col = faint)
  if (has_truth) graphics::points(x$accuracy_truth, x$time_truth, pch = 4, col = "#bb4422", cex = .8)
  graphics::points(x$accuracy_median, x$time_median, pch = 16, col = colors)
  labeled <- if (n <= 60L) seq_len(n) else unique(c(
    head(order(x$accuracy_median), 6), tail(order(x$accuracy_median), 6),
    head(order(x$time_median), 6), tail(order(x$time_median), 6)))
  graphics::text(x$accuracy_median[labeled], x$time_median[labeled], labels = x$label[labeled],
                 pos = 3, cex = .65, col = colors[labeled], offset = .4)
  if (!is.null(group)) graphics::legend("topleft", title = group_name,
    legend = levels(group), col = palette, pch = 16, bty = "n", cex = .75) else legend()
  scope <- if (is.null(effects$core_probability)) "All retained posterior draws" else
    sprintf("Conditional on CHIPS core (probability %.3f)", effects$core_probability)
  graphics::mtext(sprintf("%s: centered effects | %s | %d draws",
    if (type == "items") "Items" else "Respondents", scope, length(effects$draw_indices)),
    side = 3, outer = TRUE, line = 1, cex = 1.05)
  note <- sprintf("%.1f%% marginal intervals; scatterplot crosses are not joint credible regions. Higher time effects mean slower responses.",
                    100 * effects$level)
  graphics::mtext(note, side = 1, outer = TRUE, line = 1, cex = .85)
  invisible(effects)
}

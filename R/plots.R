# Plot helpers are deliberately base R: no plotting dependencies or state relabeling.
.plot_start <- function(file, rows, width = 1250, height = 950) {
  if (!is.null(file)) {
    dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
    grDevices::png(file, width = width, height = height, res = 120, bg = "white")
  }
  old <- graphics::par(no.readonly = TRUE)
  graphics::par(mfrow = c(rows, 2), mar = c(3.5, 4, 3, 1),
                mgp = c(2.2, .6, 0), las = 1, bg = "white")
  list(old = old, close = !is.null(file))
}

.plot_end <- function(device) {
  graphics::par(device$old)
  if (device$close) grDevices::dev.off()
  invisible(NULL)
}

.item_heatmap <- function(value, order, colors, main, zlim = range(value)) {
  I <- ncol(value); P <- nrow(value)
  graphics::image(x = seq(.5, I + .5, length.out = I + 1),
    y = seq(.5, P + .5, length.out = P + 1),
    z = t(value[order, , drop = FALSE]), col = colors, zlim = zlim,
    xlab = "Item / position", ylab = "Respondent", main = main,
    yaxt = "n", useRaster = TRUE)
  ticks <- unique(round(seq(1, P, length.out = min(P, 6))))
  labels <- if (is.null(rownames(value))) order[ticks] else rownames(value)[order[ticks]]
  graphics::axis(2, at = ticks, labels = labels)
}

#' Plot simulated data before fitting
#' @param simulation Output from simulate_irtihmm.
#' @param file Optional PNG filename. NULL plots on the current graphics device.
#' @return Invisibly, the simulation, after plotting responses, log-times, true
#'   regimes and regime occupancy by the first covariate (when supplied).
#' @export
plot_simulation <- function(simulation, file = NULL) {
  data <- simulation$data
  if (!inherits(data, "irtihmm_data") || !identical(dim(simulation$states), dim(data$x)))
    stop("Expected a simulation from simulate_irtihmm.")
  device <- .plot_start(file, rows = 2)
  on.exit(.plot_end(device), add = TRUE)
  order <- if (ncol(data$z)) order(data$z[, 1]) else seq_len(data$P)
  .item_heatmap(data$x, order, c("#f2f4f8", "#244b78"),
                "Responses: white = 0, blue = 1", c(0, 1))
  .item_heatmap(data$y, order, grDevices::hcl.colors(64, "YlOrRd", rev = TRUE), "Observed log-response times")
  graphics::legend("bottomright", legend = format(range(data$y), digits = 2),
    fill = grDevices::hcl.colors(64, "YlOrRd", rev = TRUE)[c(1, 64)], bg = "white", cex = .75)
  K <- length(simulation$emission$a)
  palette <- grDevices::hcl.colors(max(K, 2), "Dark 3")[seq_len(K)]
  .item_heatmap(simulation$states, order, palette, "True states", c(.5, K + .5))
  graphics::legend("bottomright", legend = paste("State", seq_len(K)),
                    fill = palette, bg = "white", cex = .75)
  if (ncol(data$z)) {
    covariate <- data$z[, 1]
    grouped <- length(unique(covariate)) > 6L
    group <- if (grouped) factor(ifelse(covariate <= median(covariate),
      "<= median", "> median")) else factor(covariate)
    name <- colnames(data$z)[1]
    if (is.null(name) || is.na(name) || !nzchar(name)) name <- "z[,1]"
    xlab <- paste(name, if (grouped) "(grouped)" else "")
  } else {
    group <- factor(rep("All", data$P))
    xlab <- "Respondents"
  }
  counts <- table(factor(as.numeric(simulation$states), levels = seq_len(K)),
                  rep(group, data$I))
  graphics::barplot(prop.table(counts, 2), col = palette, ylim = c(0, 1),
    xlab = xlab, ylab = "Proportion of responses", main = "State occupancy by group")
  graphics::legend("topright", legend = paste("State", seq_len(K)),
                    fill = palette, bty = "n", cex = .75)
  invisible(simulation)
}

.fit_reconstruction <- function(fit) {
  data <- fit$data; M <- length(fit$draws)
  accuracy <- time <- matrix(0, data$I, M)
  for (j in seq_len(M)) {
    draw <- fit$draws[[j]]; e <- draw$emission; C <- draw$states
    accuracy[, j] <- .item_average(data, plogis(.item_predictor(data, e$d, e$b) + matrix(e$a[C], data$P)))
    time[, j] <- .item_average(data, .item_predictor(data, e$eta, e$nu) + matrix(e$mu[C], data$P))
  }
  summarize <- function(x) {
    intervals <- t(apply(x, 1, quantile, probs = c(.025, .975), names = FALSE))
    cbind(mean = rowMeans(x), lower = intervals[, 1], upper = intervals[, 2])
  }
  list(accuracy = summarize(accuracy), time = summarize(time))
}

.reconstruction_panel <- function(summary, observed, truth, main, ylab,
                                   bounded = FALSE) {
  x <- seq_along(observed)
  ylim <- if (bounded) c(0, 1) else range(summary, observed, truth, finite = TRUE)
  graphics::plot(x, observed, type = "n", ylim = ylim,
                  xlab = "Item", ylab = ylab, main = main)
  graphics::polygon(c(x, rev(x)), c(summary[, "lower"], rev(summary[, "upper"])),
    col = grDevices::adjustcolor("#4477aa", alpha.f = .2), border = NA)
  graphics::lines(x, summary[, "mean"], col = "#244b78", lwd = 2)
  graphics::points(x, observed, pch = 16, cex = .6)
  if (!is.null(truth)) graphics::lines(x, truth, col = "#bb4422", lwd = 2, lty = 2)
  labels <- c("Observed", "Posterior mean", if (!is.null(truth)) "True mean")
  graphics::legend("topleft", legend = labels,
    col = c("black", "#244b78", if (!is.null(truth)) "#bb4422"),
    pch = c(16, NA, if (!is.null(truth)) NA),
    lty = c(NA, 1, if (!is.null(truth)) 2), bty = "n", cex = .7)
}

.partition_panel <- function(value, people, I, main, respondent_labels = NULL) {
  n <- nrow(value)
  graphics::image(seq(.5, n + .5, length.out = n + 1),
    seq(.5, n + .5, length.out = n + 1), value,
    col = grDevices::hcl.colors(64, "Blues 3", rev = TRUE), zlim = c(0, 1),
    xlab = "Responses grouped by respondent", ylab = "Responses", axes = FALSE,
    main = main, useRaster = TRUE)
  centers <- (seq_along(people) - .5) * I + .5
  labels <- if (is.null(respondent_labels)) paste0("p", people) else respondent_labels
  graphics::axis(1, at = centers, labels = labels)
  graphics::axis(2, at = centers, labels = labels)
  if (length(people) > 1) {
    boundaries <- seq_len(length(people) - 1) * I + .5
    graphics::abline(v = boundaries, h = boundaries, col = "#bbbbbb", lwd = .5)
  }
  graphics::box()
}

#' Plot posterior diagnostics and reconstruction without matching state labels
#' @param fit An irtihmm_fit object.
#' @param simulation Optional generating simulation, to compare against truth.
#' @param file Optional PNG filename. NULL plots on the current graphics device.
#' @return Invisibly, posterior item-level reconstruction summaries. Bands are
#'   95 percent posterior intervals for conditional mean accuracy and log-time,
#'   not predictive intervals for new observations. The state panels compare
#'   co-membership of responses for up to three selected respondents, a quantity
#'   invariant to state-label changes. A single-chain trace does not prove convergence.
#' @export
plot_fit <- function(fit, simulation = NULL, file = NULL) {
  if (!inherits(fit, "irtihmm_fit") || !length(fit$draws)) stop("Expected a fitted model with retained draws.")
  data <- fit$data
  if (!is.null(simulation) &&
      (!isTRUE(all.equal(simulation$data$x, data$x)) ||
       !isTRUE(all.equal(simulation$data$y, data$y)) ||
       !isTRUE(all.equal(simulation$data$z, data$z))))
    stop("simulation must describe the same data and covariates used by fit.")
  reconstruction <- .fit_reconstruction(fit)
  device <- .plot_start(file, rows = 3, height = 1300)
  on.exit(.plot_end(device), add = TRUE)
  diagnostics <- fit$diagnostics
  saved <- vapply(fit$draws, `[[`, integer(1), "iteration")
  burnin <- min(saved) - 1L  # Discarded prefix; includes any initial thinning gap.
  true_K <- if (!is.null(simulation)) length(unique(as.numeric(simulation$states))) else NULL
  graphics::plot(diagnostics$iteration, diagnostics$occupied, type = "l",
    xlab = "Iteration", ylab = "Occupied states", main = "Number of occupied states",
    ylim = range(diagnostics$occupied, true_K))
  graphics::abline(v = burnin, lty = 3, col = "#777777")
  if (!is.null(simulation)) graphics::abline(h = true_K,
    lty = 2, col = "#bb4422")
  graphics::legend("topright", legend = c("End of burn-in", if (!is.null(simulation)) "True count"),
    lty = c(3, if (!is.null(simulation)) 2), col = c("#777777", if (!is.null(simulation)) "#bb4422"),
    bty = "n", cex = .7)
  graphics::plot(diagnostics$iteration, diagnostics$conditional_loglik, type = "l",
    xlab = "Iteration", ylab = "Emission log-likelihood",
    main = "Likelihood conditional on the state path")
  graphics::abline(v = burnin, lty = 3, col = "#777777")
  truth_accuracy <- truth_time <- NULL
  if (!is.null(simulation)) {
    e <- simulation$emission; C <- simulation$states
    truth_accuracy <- .item_average(data, plogis(.item_predictor(data, e$d, e$b) + matrix(e$a[C], data$P)))
    truth_time <- .item_average(data, .item_predictor(data, e$eta, e$nu) + matrix(e$mu[C], data$P))
  }
  .reconstruction_panel(reconstruction$accuracy, .item_average(data, data$x), truth_accuracy,
    "Reconstruction: mean item accuracy", "Probability of a correct response", TRUE)
  .reconstruction_panel(reconstruction$time, .item_average(data, data$y), truth_time,
    "Reconstruction: mean item log-time", "Mean log-response time")
  people <- unique(round(seq(1, data$P, length.out = min(3, data$P))))
  n <- length(people) * data$I
  membership <- matrix(0, n, n)
  for (draw in fit$draws) {
    labels <- as.vector(t(draw$states[people, , drop = FALSE]))
    membership <- membership + outer(labels, labels, "==")
  }
  membership <- membership / length(fit$draws)
  if (!is.null(simulation)) {
    labels <- as.vector(t(simulation$states[people, , drop = FALSE]))
    .partition_panel(outer(labels, labels, "==") * 1, people, data$I,
                      "True co-membership: white = no, blue = yes", rownames(data$x)[people])
  } else {
    graphics::hist(summary(fit)$occupied_states, main = "Posterior occupied-state distribution",
      xlab = "Number of states", col = "#b7cfdf", ylab = "Frequency")
  }
  .partition_panel(membership, people, data$I,
                    "Posterior co-membership probability (0-1)", rownames(data$x)[people])
  invisible(reconstruction)
}

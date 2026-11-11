# Partitions use respondent-major order, as in Hoskovec's dlso_wrapper:
# (p=1,i=1), ..., (p=1,i=I), (p=2,i=1), ... .
.posterior_partitions <- function(fit) {
  if (!inherits(fit, "irtihmm_fit") || !length(fit$draws))
    stop("Expected a fitted model with retained draws.")
  posterior_partitions_cpp(fit$draws, fit$data$P, fit$data$I)
}

# Compare the whole subpartition, preserving separation between core clusters.
.core_draws <- function(partitions, core) {
  if (!is.integer(partitions)) {
    if (any(!is.finite(partitions) | partitions < 1 | partitions != floor(partitions) |
            partitions > .Machine$integer.max)) stop("Invalid posterior state labels.")
    storage.mode(partitions) <- "integer"
  }
  core_draws_cpp(partitions, as.integer(core), seq_len(nrow(partitions)))
}

#' Prepare a reusable in-memory posterior partition cache
#' @param fit A fitted irtihmm model.
#' @param parent Optional cache for a pooled fit, to share its matrix with a chain.
#' @param rows Rows of the parent partition matrix matching this fit's draws.
#' @return A process-local cache environment. Pass it as partition_cache to
#'   post-processing and summaries. Do not save it with results: it holds the
#'   large allocation matrix and references to states. No draw thinning is used.
#' @export
prepare_partition_cache <- function(fit, parent = NULL, rows = NULL) {
  if (!inherits(fit, "irtihmm_fit") || !length(fit$draws)) stop("Expected a fitted model with retained draws.")
  cache <- new.env(parent = emptyenv())
  if (is.null(parent)) {
    if (!is.null(rows)) stop("rows requires a parent cache.")
    cache$partitions <- .posterior_partitions(fit)
    cache$rows <- seq_along(fit$draws)
  } else {
    if (!inherits(parent, "irtihmm_partition_cache")) stop("Invalid parent cache.")
    if (!is.numeric(rows) || length(rows) != length(fit$draws) || anyNA(rows) ||
        any(!is.finite(rows) | rows < 1 | rows > nrow(parent$partitions) | rows != floor(rows)) ||
        anyDuplicated(rows)) stop("Invalid parent partition rows.")
    cache$partitions <- parent$partitions; cache$rows <- as.integer(rows)
    if (!partitions_match_cpp(cache$partitions, fit$draws, fit$data$P, fit$data$I, cache$rows))
      stop("The partition cache does not match the fitted draws.")
  }
  cache$dimensions <- c(fit$data$P, fit$data$I)
  cache$verified_partitions <- cache$partitions
  cache$verified_rows <- cache$rows
  cache$verified_states <- lapply(fit$draws, `[[`, "states")
  cache$alignment <- NULL
  class(cache) <- "irtihmm_partition_cache"
  cache
}

.checked_partition_cache <- function(fit, cache = NULL) {
  if (is.null(cache)) return(prepare_partition_cache(fit))
  if (!inherits(cache, "irtihmm_partition_cache") ||
      !identical(cache$dimensions, c(fit$data$P, fit$data$I)))
    stop("The partition cache does not match the fitted data dimensions.")
  states <- lapply(fit$draws, `[[`, "states")
  # Unmodified matrices share references, so identical() avoids repeated scans.
  # R copy-on-modify keeps these snapshots valid when a caller edits a draw.
  if (!identical(cache$verified_states, states) ||
      !identical(cache$verified_partitions, cache$partitions) ||
      !identical(cache$verified_rows, cache$rows)) {
    if (!partitions_match_cpp(cache$partitions, fit$draws, fit$data$P, fit$data$I, cache$rows))
      stop("The partition cache does not match the fitted draws.")
    cache$alignment <- NULL
    cache$verified_states <- states; cache$verified_partitions <- cache$partitions
    cache$verified_rows <- cache$rows
  }
  cache
}

.core_alignment <- function(fit, postprocess, partition_cache = NULL) {
  if (!inherits(postprocess, "irtihmm_postprocess") || is.null(postprocess$chips))
    stop("Expected post-processing with CHIPS enabled.")
  core <- postprocess$chips$subpartition
  if (!identical(dim(core), dim(fit$data$x)) || all(is.na(core)) ||
      any(!is.finite(core[!is.na(core)]) | core[!is.na(core)] < 1 |
          core[!is.na(core)] != floor(core[!is.na(core)])))
    stop("The CHIPS core must contain positive integer state labels and match the fitted data.")
  labels <- as.integer(t(core)); included <- !is.na(labels); labels[!included] <- -1L
  cache <- .checked_partition_cache(fit, partition_cache)
  if (is.null(cache$alignment) || !identical(cache$alignment$labels, labels)) {
    indices <- which(core_draws_cpp(cache$partitions, labels, cache$rows))
    if (!length(indices)) stop("No fitted draws are compatible with the CHIPS core.")
    ids <- sort(unique(labels[included]))
    anchors <- vapply(ids, function(s) which(labels == s)[1], integer(1))
    cache$alignment <- list(labels = labels, included = included, indices = indices,
      probability = length(indices) / length(fit$draws), ids = ids, anchors = anchors,
      sizes = vapply(ids, function(s) sum(labels == s), integer(1)),
      matched = cache$partitions[cache$rows[indices], anchors, drop = FALSE])
    colnames(cache$alignment$matched) <- as.character(ids)
  }
  out <- cache$alignment
  if (!is.null(postprocess$chips$probability) && abs(out$probability - postprocess$chips$probability) > 1e-12)
    stop("The fitted draws do not match the saved CHIPS probability; recompute post-processing for this fit.")
  out
}

#' Summarize posterior state partitions with SALSO and optionally CHIPS
#' @param fit An irtihmm_fit object with retained draws.
#' @param loss Partition loss, either "VI" or "binder", with equal costs.
#' @param chips Whether to also compute the CHIPS core and CHIPS-SALSO completion.
#' @param credible Minimum posterior probability of the simultaneous CHIPS core.
#' @param n_runs Number of randomized search restarts for each algorithm.
#' @param cores Number of parallel search cores. Defaults to one.
#' @param seed Optional random seed. The caller's RNG state is restored afterwards.
#' @param partition_cache Optional output from prepare_partition_cache, reused by summaries.
#' @return A list of class irtihmm_postprocess containing SALSO states (P by I),
#'   posterior occupied-state counts, change probabilities (P by I-1), and optional
#'   CHIPS results. CHIPS includes a subpartition (NA outside its core), a completed
#'   partition, the empirical joint core probability and compatible draw indices.
#'   CHIPS runs directly at credible; no full probability/core-size curve is built.
#'   Original draws and labels are unchanged.
#' @details All respondent-item responses are partitioned jointly. Numeric labels
#'   of the summaries do not index the original emission or transition arrays.
#'   VI does not put a special cost on misplaced temporal change points. CHIPS
#'   credible probability concerns the whole core jointly, not individual labels.
#'   No fixed number of clusters is imposed beyond the package's default search
#'   bound (the largest occupied count in the retained draws). These summaries
#'   neither change the posterior nor resolve location identifiability.
#' @references Dahl, Johnson and Mueller (2022), doi:10.1080/10618600.2022.2069779.
#'   Page, Barrientos, Dahl and Dunson, Uncertainty Quantification in Bayesian
#'   Clustering, arXiv:2511.16040.
#' @export
postprocess_irtihmm <- function(fit, loss = c("VI", "binder"), chips = TRUE,
                               credible = .95, n_runs = 4L, cores = 1L,
                               seed = NULL, partition_cache = NULL) {
  loss <- match.arg(loss)
  cache <- .checked_partition_cache(fit, partition_cache)
  Z <- if (identical(cache$rows, seq_len(nrow(cache$partitions)))) cache$partitions else
    cache$partitions[cache$rows, , drop = FALSE]
  if (!is.logical(chips) || length(chips) != 1L || is.na(chips))
    stop("chips must be TRUE or FALSE.")
  if (!is.numeric(credible) || length(credible) != 1L || !is.finite(credible) ||
      credible <= 0 || credible > 1) stop("credible must lie in (0, 1].")
  for (value in list(n_runs, cores))
    if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
        value < 1 || value != floor(value) || value > .Machine$integer.max)
      stop("n_runs and cores must be positive integers.")
  if (!requireNamespace("salso", quietly = TRUE))
    stop('Post-processing requires salso. Install it once with install.packages("salso").')
  if (chips && !"chips" %in% getNamespaceExports("salso"))
    stop("Update salso: this version does not export chips().")
  if (!is.null(seed)) {
    previous <- get0(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    on.exit(if (is.null(previous)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    } else assign(".Random.seed", previous, envir = .GlobalEnv), add = TRUE)
    set.seed(seed)
  }
  P <- fit$data$P; I <- fit$data$I
  reshape <- function(x) matrix(as.integer(x), P, I, byrow = TRUE,
                                dimnames = dimnames(fit$data$x))
  estimate <- salso::salso(Z, loss = loss, nRuns = as.integer(n_runs),
                            nCores = as.integer(cores))
  # No quadratic PSM is required for VI or for these temporal summaries.
  change <- matrix(0, P, max(0L, I - 1L))
  for (draw in fit$draws)
    change <- change + (draw$states[, -1L, drop = FALSE] !=
                          draw$states[, -I, drop = FALSE])
  result <- list(states = reshape(estimate), loss = loss,
    expected_loss = as.numeric(salso::partition.loss(Z, as.integer(estimate), loss)),
    occupied_states = apply(Z, 1L, function(x) length(unique(x))),
    change_probability = change / nrow(Z),
    chips = NULL, covariates = fit$data$z,
    settings = list(credible = credible, n_runs = n_runs, cores = cores,
                    seed = seed, chips_search = "direct", salso_version = as.character(utils::packageVersion("salso"))))
  if (chips) {
    page <- salso::chips(Z, threshold = credible, intermediateResults = FALSE,
      loss = loss, nRuns = as.integer(n_runs), nCores = as.integer(cores), andSALSO = TRUE)
    core <- as.integer(page$chips_partition)
    included <- core >= 0L
    core[!included] <- NA_integer_
    result$chips <- list(subpartition = reshape(core),
      states = reshape(page$chips_and_salso_partition),
      probability = NULL, core_fraction = mean(included), curve = NULL)
    aligned <- .core_alignment(fit, structure(result, class = "irtihmm_postprocess"), cache)
    if (aligned$probability + 1e-12 < credible)
      stop("CHIPS returned a core below the requested empirical credible probability.")
    result$chips$probability <- aligned$probability
    result$chips$compatible_draws <- aligned$indices
  }
  structure(result, class = "irtihmm_postprocess")
}

#' Plot representative state partitions and posterior uncertainty
#' @param postprocess Output from postprocess_irtihmm.
#' @param file Optional PNG filename. NULL plots on the current graphics device.
#' @return Invisibly, postprocess. State colors are local to each partition panel;
#'   grey cells are outside the CHIPS core. Respondents are ordered by their first
#'   covariate when available. The last panels show posterior changes and counts.
#' @export
plot_postprocess <- function(postprocess, file = NULL) {
  if (!inherits(postprocess, "irtihmm_postprocess"))
    stop("Expected output from postprocess_irtihmm.")
  x <- postprocess; P <- nrow(x$states); I <- ncol(x$states)
  order <- if (ncol(x$covariates)) order(x$covariates[, 1]) else seq_len(P)
  device <- .plot_start(file, rows = if (is.null(x$chips)) 2 else 3,
                        height = if (is.null(x$chips)) 950 else 1300)
  on.exit(.plot_end(device), add = TRUE)
  states_panel <- function(C, title, core = FALSE) {
    labels <- sort(unique(as.integer(C[!is.na(C)])))
    encoded <- matrix(match(C, labels), P, I)
    palette <- grDevices::hcl.colors(max(length(labels), 2), "Dark 3")[seq_along(labels)]
    if (core) {
      encoded <- encoded + 1L
      encoded[is.na(encoded)] <- 1L
      palette <- c("#dddddd", palette)
    }
    .item_heatmap(encoded, order, palette, title, c(.5, length(palette) + .5))
    graphics::legend("bottomright", legend = c(if (core) "Outside core", paste("State", labels)),
      fill = palette, bg = "white", cex = .65)
  }
  states_panel(x$states, paste("SALSO -", x$loss))
  if (!is.null(x$chips)) {
    states_panel(x$chips$subpartition,
      sprintf("CHIPS core: prob. %.3f, responses %.1f%%",
              x$chips$probability, 100 * x$chips$core_fraction), TRUE)
    states_panel(x$chips$states, "CHIPS + SALSO (own labels)")
    curve <- x$chips$curve
    if (is.null(curve)) {
      graphics::plot(x$chips$core_fraction, x$chips$probability, pch = 19,
        col = "#244b78", ylim = c(0, 1), xlim = c(0, 1),
        xlab = "Proportion of responses in the core", ylab = "Joint probability",
        main = "Selected CHIPS core (direct search)")
    } else {
      indices <- order(curve$n_items)
      graphics::plot(curve$n_items[indices] / (P * I), curve$probability[indices],
        type = "l", ylim = c(0, 1), xlim = c(0, 1),
        xlab = "Proportion of responses in the core", ylab = "Joint probability",
        main = "CHIPS curve: partition uncertainty")
    }
    graphics::abline(h = x$settings$credible, lty = 2, col = "#bb4422")
  }
  if (I > 1L) {
    .item_heatmap(x$change_probability, order,
      grDevices::hcl.colors(64, "Blues 3", rev = TRUE),
      "Change probability: position j means j -> j+1", c(0, 1))
  } else {
    graphics::plot.new()
    graphics::title("One item: no observable state changes")
  }
  counts <- table(x$occupied_states)
  graphics::barplot(counts / sum(counts), xlab = "Occupied states per draw",
    ylab = "Posterior probability", col = "#b7cfdf",
    main = paste("Original posterior; SALSO:", length(unique(as.integer(x$states)))))
  invisible(x)
}

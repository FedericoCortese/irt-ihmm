.new_emission <- function(data, K, hyper, copula) {
  empirical_sd <- sd(as.numeric(data$y))
  if (!is.finite(empirical_sd)) empirical_sd <- .5
  list(b = rep(0, data$I), a = rnorm(K, sd = .3), d = rep(0, data$P),
       nu = .item_average(data, data$y), mu = rnorm(K, sd = .1), eta = rep(0, data$P),
       log_sigma = rep(log(max(empirical_sd, .1)), K),
       rho = if (copula) 2 * rbeta(K, hyper$rho_shape, hyper$rho_shape) - 1 else rep(0, K))
}

.expand_emission <- function(emission, K, hyper) {
  old <- length(emission$a)
  if (K <= old) return(emission)
  n <- K - old
  emission$a <- c(emission$a, rnorm(n, sd = sqrt(hyper$accuracy$variances[2])))
  emission$mu <- c(emission$mu, rnorm(n, sd = sqrt(hyper$time$variances[2])))
  emission$log_sigma <- c(emission$log_sigma, rnorm(n, hyper$scale$mean, sqrt(hyper$scale$variance)))
  emission$rho <- c(emission$rho, if (hyper$copula)
    2 * rbeta(n, hyper$rho_shape, hyper$rho_shape) - 1 else rep(0, n))
  .apply_fixed_params(emission, hyper$fixed_params)
}

.set_margin <- function(emission, margin, effects) {
  names <- if (margin == "accuracy") c("b", "a", "d") else c("nu", "mu", "eta")
  emission[names] <- effects[c("item", "state", "person")]
  emission
}

.original_time_design <- function(C, I, K, P, item_index = NULL) {
  # Dense mathematical reference for tests; the sampler accumulates sparse statistics.
  # Rows follow as.numeric(data$y), i.e. respondent varies fastest.
  original <- matrix(0, P * I, I + K + P)
  p <- rep(seq_len(P), I); i <- rep(seq_len(I), each = P)
  if (!is.null(item_index)) i <- as.numeric(item_index)
  original[cbind(seq_len(P * I), i)] <- 1
  original[cbind(seq_len(P * I), I + as.numeric(C))] <- 1
  original[cbind(seq_len(P * I), I + K + p)] <- 1
  original
}

.time_design <- function(C, I, K, P, item_index = NULL) {
  original <- .original_time_design(C, I, K, P, item_index) 
  # Canonical inverse with lambda=(0,0): g=G, h=theta+S, u=R.
  inverse <- matrix(0, I + K + P, I + K + P - 2)
  inverse[I + seq_len(K), 1] <- 1
  if (I > 1) inverse[cbind(2:I, 2:I)] <- 1
  if (K > 1) inverse[cbind(I + 2:K, I + seq_len(K - 1))] <- 1
  if (P > 1) inverse[cbind(I + K + 2:P, I + K + seq_len(P - 1) - 1)] <- 1
  original %*% inverse
}

.update_original_margin <- function(data, emission, hyper, C, margin, copula, fixed_params = list()) {
  I <- data$I; P <- data$P; K <- length(emission$a)
  nm <- if (margin == "accuracy") c("b", "a", "d") else c("nu", "mu", "eta")
  unpack <- function(value) list(item = value[seq_len(I)],
    state = value[I + seq_len(K)], person = value[I + K + seq_len(P)])
  current <- unlist(emission[nm], use.names = FALSE)
  prior_mean <- c(rep(hyper[[margin]]$mean, I), rep(0, K + P))
  prior_variance <- rep(hyper[[margin]]$variances, c(I, K, P))
  free <- !rep(nm %in% names(fixed_params), c(I, K, P))
  if (!any(free)) return(list(effects = unpack(current), evaluations = 0L))
  template <- current
  expand <- function(value) { out <- template; out[free] <- value; out }
  current <- current[free]
  prior_mean <- prior_mean[free]; prior_variance <- prior_variance[free]
  evaluations <- 0L
  if (margin == "time" && !copula) {
    weight <- as.numeric(exp(-2 * emission$log_sigma[C]))
    stats <- gaussian_time_stats_cpp(C, data$item_index, data$y, weight, K, FALSE)
    precision <- diag(1 / prior_variance, nrow = length(prior_variance)) +
      stats$precision[free, free, drop = FALSE]
    information <- prior_mean / prior_variance + stats$information[free] -
      stats$precision[free, !free, drop = FALSE] %*% template[!free]
    value <- .rmvn_precision(precision, information)
  } else {
    context <- emission_context_cpp(data, emission, C, margin == "time")
    loglik <- function(value) {
      effects <- unpack(expand(value))
      emission_margin_loglik_cpp(context, effects$item, effects$state, effects$person)
    }
    result <- .elliptical_slice(current, prior_mean, NULL, loglik,
      direction = rnorm(length(current), sd = sqrt(prior_variance)))
    value <- result$value; evaluations <- result$evaluations
  }
  list(effects = unpack(expand(value)), evaluations = evaluations)
}

.update_emission <- function(data, emission, hyper, C, priors, copula, gustafson = TRUE, fixed_params = list()) {
  I <- data$I; P <- data$P; K <- length(emission$a)
  evaluations <- c(accuracy = 0L, time = 0L)
  reduced <- list()
  for (margin in c("accuracy", "time")) {
    nm <- if (margin == "accuracy") c("b", "a", "d") else c("nu", "mu", "eta")
    if (!gustafson || any(nm %in% names(fixed_params))) {
      result <- .update_original_margin(data, emission, hyper, C, margin, copula, fixed_params)
      emission <- .set_margin(emission, margin, result$effects)
      evaluations[margin] <- result$evaluations
      # Saved reduced coordinates are diagnostics only in the direct sampler.
      reduced[[margin]] <- reduce_effects(result$effects$item,
        result$effects$state, result$effects$person)
      next
    }
    phi <- reduce_effects(emission[[nm[1]]], emission[[nm[2]]], emission[[nm[3]]])
    if (margin == "time" && !copula) {
      prior <- location_prior(I, K, P, hyper[[margin]]$mean, hyper[[margin]]$variances)
      weight <- as.numeric(exp(-2 * emission$log_sigma[C]))
      stats <- gaussian_time_stats_cpp(C, data$item_index, data$y, weight, K, TRUE)
      precision_prior <- chol2inv(chol(prior$covariance))
      precision <- precision_prior + stats$precision
      information <- precision_prior %*% prior$mean + stats$information
      phi <- .rmvn_precision(precision, information)
    } else {
      context <- emission_context_cpp(data, emission, C, margin == "time")
      location_loglik <- function(value) {
        effects <- .effect_from_phi(value, I, K, P)
        emission_margin_loglik_cpp(context, effects$item, effects$state, effects$person)
      }
      # T times independent original-prior normals has the exact reduced prior
      # covariance, avoiding a dense covariance matrix and Cholesky for ESS.
      v <- hyper[[margin]]$variances
      direction <- reduce_effects(rnorm(I, sd = sqrt(v[1])),
        rnorm(K, sd = sqrt(v[2])), rnorm(P, sd = sqrt(v[3])))
      mean <- c(hyper[[margin]]$mean, rep(0, length(phi) - 1L))
      result <- .elliptical_slice(phi, mean, NULL, location_loglik, direction)
      phi <- result$value; evaluations[margin] <- result$evaluations
    }
    reduced[[margin]] <- phi
    emission <- .set_margin(emission, margin, .effect_from_phi(phi, I, K, P))
  }
  # Recover positions only for the collapsed sampler. Both variants update the
  # same original hyperpriors using their current absolute effect vectors.
  locations <- list()
  for (margin in c("accuracy", "time")) {
    nm <- if (margin == "accuracy") c("b", "a", "d") else c("nu", "mu", "eta")
    if (gustafson && !any(nm %in% names(fixed_params))) {
      recovered <- recover_effects(reduced[[margin]], I, K, P, hyper[[margin]]$mean,
                                hyper[[margin]]$variances)
      emission <- .set_margin(emission, margin, recovered)
      locations[[margin]] <- recovered[c("lambda", "location_mean", "location_covariance")]
    } else {
      nm <- if (margin == "accuracy") c("b", "a", "d") else c("nu", "mu", "eta")
      recovered <- setNames(emission[nm], c("item", "state", "person"))
    }
    hyper[[margin]] <- .update_location_hyper(recovered$item, recovered$state,
      recovered$person, hyper[[margin]], priors[[margin]],
      fixed = nm %in% names(fixed_params))
  }
  accuracy <- .item_predictor(data, emission$d, emission$b) + matrix(emission$a[C], P)
  time_mean <- .item_predictor(data, emission$eta, emission$nu) + matrix(emission$mu[C], P)
  for (s in seq_len(K)) {
    idx <- which(C == s)
    if (!length(idx)) {
      if (!"log_sigma" %in% names(fixed_params)) emission$log_sigma[s] <- rnorm(1, hyper$scale$mean, sqrt(hyper$scale$variance))
      if (!"rho" %in% names(fixed_params)) emission$rho[s] <- if (copula) 2 * rbeta(1, priors$rho_shape, priors$rho_shape) - 1 else 0
      next
    }
    context <- emission_scale_context_cpp(data$x[idx], data$y[idx], accuracy[idx], time_mean[idx])
    scale_loglik <- function(log_sigma, rho)
      emission_scale_loglik_cpp(context, log_sigma, rho)
    if (!"log_sigma" %in% names(fixed_params)) emission$log_sigma[s] <- .slice_scalar(emission$log_sigma[s], function(value)
      scale_loglik(value, emission$rho[s]) +
        dnorm(value, hyper$scale$mean, sqrt(hyper$scale$variance), log = TRUE))
    if (copula && !"rho" %in% names(fixed_params)) {
      density <- function(value) {
        rho <- tanh(value)
        if (!is.finite(rho) || abs(rho) >= 1) return(-Inf)
        scale_loglik(emission$log_sigma[s], rho) + priors$rho_shape * log1p(-rho^2)
      }
      emission$rho[s] <- tanh(.slice_scalar(atanh(emission$rho[s]), density))
    }
  }
  if (!"log_sigma" %in% names(fixed_params)) {
    scale_prior <- priors$log_scale
    hyper$scale$variance <- .rinvgamma(scale_prior$shape + K / 2,
      scale_prior$rate + sum((emission$log_sigma - hyper$scale$mean)^2) / 2)
    variance <- 1 / (1 / scale_prior$mean_sd^2 + K / hyper$scale$variance)
    mean <- variance * (scale_prior$mean / scale_prior$mean_sd^2 +
                         sum(emission$log_sigma) / hyper$scale$variance)
    hyper$scale$mean <- rnorm(1, mean, sqrt(variance))
  }
  list(emission = emission, hyper = hyper, locations = locations,
       reduced = reduced, evaluations = evaluations)
}

#' Fit the joint IRT and response-time model
#' @param data A validated irtihmm_data object.
#' @param mode Infinite beam sampling or a finite-state reference model.
#' @param K Initial prefix length (infinite) or fixed number of states (finite).
#' @param iterations Total number of MCMC iterations.
#' @param burnin Number of initial iterations discarded.
#' @param thin Keep every thin-th iteration after burnin.
#' @param copula Whether to estimate state-specific correlations; FALSE fixes rho=0.
#' @param states Optional known state-path matrix; freezes the state block for validation.
#' @param initial_states Optional initial path matrix when states are unknown.
#' @param priors Proper prior settings from irtihmm_priors.
#' @param seed Optional integer seed passed to set.seed.
#' @param max_states Safety bound for beam expansion; reaching it aborts, never truncates.
#' @param progress Whether to print progress every 100 iterations.
#' @param gustafson TRUE uses location-reduced updates and exact Gaussian recovery.
#'   FALSE samples original additive effects directly under the same proper priors.
#'   Both target the same model posterior; the flag does not turn prior correction
#'   weights on/off or impose location constraints. Defaults to TRUE.
#' @param fixed_params Named list fixing complete emission blocks: b, a, d, nu,
#'   mu, eta, log_sigma, rho. Scalars broadcast; full vectors are also accepted
#'   for item/person effects and finite-state effects. Infinite-state effects
#'   require scalars. Hyperparameters belonging only to fixed blocks are inactive
#'   and kept at initialization; they are not posterior estimates.
#'   Margins containing fixed blocks use constrained original-coordinate updates;
#'   Gustafson reduction remains active for unconstrained margins.
#' @param persistence Whether to estimate a global forced-persistence probability
#'   omega. FALSE recovers the original probit model without extra RNG draws.
#' @param omega_prior Positive length-two vector of Beta shape parameters for omega.
#'   Used only when persistence is TRUE; defaults to c(1,1).
#' @return An irtihmm_fit object containing complete retained draws, iteration
#'   diagnostics, data, priors and upstream commit provenance. State labels are
#'   ordered stick indices and are not compacted. The first release supports
#'   complete responses, respondent-specific item permutations and respondent covariates.
#' @export
fit_irtihmm <- function(data, mode = c("infinite", "finite"), K = 2L,
                        iterations = 1000L, burnin = 500L, thin = 1L,
                        copula = TRUE, states = NULL, initial_states = NULL,
                        priors = irtihmm_priors(), seed = NULL,
                        max_states = 200L, progress = FALSE, gustafson = TRUE,
                        fixed_params = list(), persistence = FALSE,
                        omega_prior = c(1, 1)) {
  call <- match.call(); mode <- match.arg(mode)
  if (!inherits(data, "irtihmm_data")) stop("Use irtihmm_data to construct data.")
  # Revalidate fields so edited data objects cannot bypass the public contract.
  data <- irtihmm_data(data$x, data$y, data$z, data$item_index)
  K <- .integer_scalar(K, "K"); iterations <- .integer_scalar(iterations, "iterations")
  burnin <- .integer_scalar(burnin, "burnin", 0L); thin <- .integer_scalar(thin, "thin")
  max_states <- .integer_scalar(max_states, "max_states")
  if (burnin >= iterations || thin > iterations - burnin) stop("No retained iterations.")
  if (!is.logical(copula) || length(copula) != 1 || is.na(copula)) stop("copula must be TRUE or FALSE.")
  if (!is.logical(gustafson) || length(gustafson) != 1 || is.na(gustafson))
    stop("gustafson must be TRUE or FALSE.")
  fixed_params <- .validate_fixed_params(fixed_params, data, K, mode, copula)
  priors <- .validate_priors(priors)
  if (!is.logical(persistence) || length(persistence) != 1L || is.na(persistence))
    stop("persistence must be TRUE or FALSE.")
  if (!is.numeric(omega_prior) || length(omega_prior) != 2L ||
      any(!is.finite(omega_prior)) || any(omega_prior <= 0))
    stop("omega_prior must contain two positive finite Beta shapes.")
  if (!is.null(seed)) set.seed(seed)
  if (!is.null(states) && !is.null(initial_states)) stop("Supply states or initial_states, not both.")
  C <- if (!is.null(states)) states else initial_states
  if (is.null(C)) C <- matrix(sample.int(K, data$P * data$I, replace = TRUE), data$P, data$I)
  if (!is.matrix(C) || !is.numeric(C) || !identical(dim(C), dim(data$x)) ||
      any(!is.finite(C)) || any(C < 1) || any(C != floor(C))) stop("Invalid state matrix.")
  if (max(C) > K) stop("State labels exceed K.")
  if (mode == "infinite" && K > max_states) stop("K exceeds max_states.")
  storage.mode(C) <- "integer"
  hyper <- .initial_hyper(priors); hyper$copula <- copula; hyper$fixed_params <- fixed_params
  emission <- .new_emission(data, K, hyper, copula)
  emission <- .apply_fixed_params(emission, fixed_params)
  tr <- .new_transition(K, ncol(data$z), mode, priors$m_alpha)
  if (persistence) tr$omega <- omega_prior[1] / sum(omega_prior)
  save_iterations <- seq.int(burnin + thin, iterations, by = thin)
  draws <- vector("list", length(save_iterations)); saved <- 0L
  diagnostics <- data.frame(iteration = seq_len(iterations), occupied = NA_integer_,
    retained_prefix = NA_integer_, beam_instantiated = NA_integer_,
    max_residual = NA_real_, min_slice = NA_real_, conditional_loglik = NA_real_,
    finite_marginal_loglik = NA_real_, accuracy_evaluations = NA_integer_, time_evaluations = NA_integer_)
  started <- proc.time()[[3]]
  for (iteration in seq_len(iterations)) {
    if (is.null(states)) {
      state_result <- if (mode == "finite") .finite_state_draw(data, emission, tr) else
        .beam_state_draw(data, emission, tr, C, hyper, max_states)
      C <- state_result$C; emission <- state_result$emission; tr <- state_result$transition
    } else state_result <- list(instantiated = K, max_residual = NA_real_, min_slice = NA_real_)
    # State slices have been discarded before the unconstrained probit update.
    tr <- .update_transition(tr, C, data$z, priors$m_alpha, omega_prior)
    result <- .update_emission(data, emission, hyper, C, priors, copula, gustafson, fixed_params)
    emission <- result$emission; hyper <- result$hyper
    diagnostics[iteration, c("occupied", "retained_prefix", "beam_instantiated")] <-
      c(length(unique(as.numeric(C))), tr$K, state_result$instantiated)
    diagnostics$max_residual[iteration] <- state_result$max_residual
    diagnostics$min_slice[iteration] <- state_result$min_slice
    diagnostics$conditional_loglik[iteration] <- .emission_path_loglik(data, emission, C)
    if (mode == "finite") {
      probabilities <- .transition_probs(tr, data$z)
      loge <- .emission_array(data, emission)
      diagnostics$finite_marginal_loglik[iteration] <- sum(vapply(seq_len(data$P), function(p)
        hmm_loglik(matrix(loge[p, , ], data$I, tr$K), probabilities$initial[p, ],
                   probabilities$transition[[p]]), numeric(1)))
    }
    diagnostics$accuracy_evaluations[iteration] <- result$evaluations[1]
    diagnostics$time_evaluations[iteration] <- result$evaluations[2]
    if (iteration %in% save_iterations) {
      saved <- saved + 1L
      draws[[saved]] <- list(iteration = iteration, states = C, emission = emission,
        transition = tr, hyper = hyper, reduced = result$reduced, locations = result$locations)
    }
    if (progress && (iteration %% 100L == 0L || iteration == iterations))
      message("Iteration ", iteration, ": ", diagnostics$occupied[iteration],
              " occupied states; prefix ", tr$K)
  }
  structure(list(call = call, data = data, mode = mode, copula = copula,
    known_states = !is.null(states), gustafson = gustafson, fixed_params = fixed_params,
    persistence = persistence, omega_prior = if (persistence) omega_prior else NULL,
    draws = draws, diagnostics = diagnostics,
    priors = priors, elapsed = proc.time()[[3]] - started,
    upstream = list(repository = "https://github.com/lvhoskovec/psbpHMM",
                    commit = "1240b4da72a504dd5c97112f0b06f4b4a40f602a")), class = "irtihmm_fit")
}

#' @export
print.irtihmm_fit <- function(x, ...) {
  cat("IRT-iHMM research fit (", x$mode, ")\n", sep = "")
  cat(nrow(x$data$x), "respondents;", ncol(x$data$x), "items;", length(x$draws), "retained draws\n")
  cat("Copula:", x$copula, "; known states:", x$known_states, "\n")
  if (isTRUE(x$persistence)) cat("Forced persistence: estimated omega (Beta prior)\n")
  cat("Location update:", if (isFALSE(x$gustafson)) "direct original effects" else
    "Gustafson-type reduction and Gaussian recovery", "\n")
  if (length(x$fixed_params)) cat("Fixed emission blocks:", paste(names(x$fixed_params), collapse = ", "),
    "(constrained margins use original-coordinate updates)\n")
  cat("Elapsed:", round(x$elapsed, 2), "seconds\n")
  invisible(x)
}

#' Summarize location-invariant fit quantities
#' @param object A fitted irtihmm model.
#' @param ... Unused.
#' @return Means and 95 percent equal-tail intervals for item and person contrasts,
#'   and retained occupied-state counts. State-specific effects require alignment.
#' @export
summary.irtihmm_fit <- function(object, ...) {
  summarize <- function(values) t(apply(do.call(rbind, values), 2L, function(x)
    c(mean = mean(x), lower = quantile(x, .025, names = FALSE), upper = quantile(x, .975, names = FALSE))))
  list(item_accuracy = summarize(lapply(object$draws, function(d) d$emission$b - d$emission$b[1])),
       item_logtime = summarize(lapply(object$draws, function(d) d$emission$nu - d$emission$nu[1])),
       person_accuracy = summarize(lapply(object$draws, function(d) d$emission$d - d$emission$d[1])),
       person_logtime = summarize(lapply(object$draws, function(d) d$emission$eta - d$emission$eta[1])),
       occupied_states = vapply(object$draws, function(d) length(unique(as.numeric(d$states))), integer(1)))
}

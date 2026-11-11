#' Evaluate a finite HMM log likelihood
#' @param log_emission Matrix with time points in rows and states in columns.
#' @param initial Initial state probability vector summing to one.
#' @param transition Row-stochastic state transition matrix.
#' @return The marginal log likelihood computed by the forward algorithm.
#' @export
hmm_loglik <- function(log_emission, initial, transition) {
  .validate_hmm(log_emission, initial, transition)
  forward_backward_cpp(log_emission, initial, transition, numeric(), FALSE)$loglik
}

.validate_hmm <- function(log_emission, initial, transition) {
  if (!is.matrix(log_emission) || !is.numeric(log_emission) ||
      !length(log_emission) || anyNA(log_emission) || any(log_emission == Inf))
    stop("Invalid log-emission matrix.")
  K <- ncol(log_emission)
  if (length(initial) != K || !is.matrix(transition) ||
      !identical(dim(transition), c(K, K)) ||
      any(!is.finite(c(initial, transition))) || any(c(initial, transition) < 0) ||
      abs(sum(initial) - 1) > 1e-10 || any(abs(rowSums(transition) - 1) > 1e-10))
    stop("Require a normalized initial vector and row-stochastic transition matrix.")
}

.finite_state_draw <- function(data, emission, tr) {
  probs <- .transition_probs(tr, data$z)
  loge <- .emission_array(data, emission)
  result <- hmm_paths_cpp(loge, probs$initial, probs$transition, matrix(numeric(), 0, 0))
  C <- result$C; loglik <- result$loglik
  list(C = C, emission = emission, transition = tr, instantiated = tr$K,
       max_residual = 0, min_slice = NA_real_, loglik = loglik)
}

.beam_state_draw <- function(data, emission, tr, C, hyper, max_states) {
  probs <- .transition_probs(tr, data$z)
  slices <- beam_slices_cpp(C, probs$initial, probs$transition)
  bound <- min(slices)
  while (max(probs$residual) >= bound) {
    if (tr$K >= max_states)
      stop("Beam expansion reached max_states before exhausting the slice tail; ",
           "increase max_states. No truncated draw was accepted.")
    tr <- .expand_transition(tr, tr$K + 1L)
    emission <- .expand_emission(emission, tr$K, hyper)
    probs <- .transition_probs(tr, data$z)
  }
  instantiated <- tr$K; residual <- max(probs$residual)
  loge <- .emission_array(data, emission)
  C <- hmm_paths_cpp(loge, probs$initial, probs$transition, slices)$C
  # The slices are discarded here. Keep prefix labels; do not compact gaps.
  retained <- max(C)
  tr <- .trim_transition(tr, retained)
  for (name in c("a", "mu", "log_sigma", "rho")) emission[[name]] <- emission[[name]][seq_len(retained)]
  list(C = C, emission = emission, transition = tr, instantiated = instantiated,
       max_residual = residual, min_slice = bound, loglik = NA_real_)
}

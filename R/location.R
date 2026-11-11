#' Reduce additive effects to location-invariant coordinates
#' @param item Item effects.
#' @param state State effects; the first entry is the reference state.
#' @param person Person effects.
#' @return The vector (overall predictor, item contrasts, state contrasts, person contrasts).
#' @export
reduce_effects <- function(item, state, person) {
  if (any(lengths(list(item, state, person)) < 1L) ||
      any(!is.finite(c(item, state, person)))) stop("Effects must be finite and nonempty.")
  c(item[1] + state[1] + person[1], item[-1] - item[1],
    state[-1] - state[1], person[-1] - person[1])
}

.location_parts <- function(phi, I, K, P) {
  if (length(phi) != I + K + P - 2L) stop("Incorrect reduced coordinate length.")
  G <- c(0, phi[seq_len(I - 1L) + 1L])
  S <- c(0, phi[seq_len(K - 1L) + I])
  R <- c(0, phi[seq_len(P - 1L) + I + K - 1L])
  list(theta = phi[1], G = G, S = S, R = R)
}

.location_transform <- function(I, K, P) {
  T <- matrix(0, I + K + P - 2L, I + K + P)
  T[1, c(1L, I + 1L, I + K + 1L)] <- 1
  row <- 1L
  for (block in list(c(0L, I), c(I, K), c(I + K, P))) {
    if (block[2] > 1L) for (j in 2:block[2]) {
      row <- row + 1L
      T[row, block[1] + c(1L, j)] <- c(-1, 1)
    }
  }
  T
}

#' Induced Gaussian prior on reduced additive effects
#' @param I Number of items.
#' @param K Number of retained states.
#' @param P Number of respondents.
#' @param mean_item Population mean of item effects.
#' @param variances Positive variances for item, state, person effects, in that order.
#' @return A list with mean and covariance of the exact induced prior.
#' @export
location_prior <- function(I, K, P, mean_item = 0, variances = c(1, 1, 1)) {
  I <- .integer_scalar(I, "I"); K <- .integer_scalar(K, "K")
  P <- .integer_scalar(P, "P")
  if (length(variances) != 3L || any(!is.finite(variances)) || any(variances <= 0) ||
      length(mean_item) != 1L || !is.finite(mean_item)) stop("Invalid location hyperparameters.")
  list(mean = c(mean_item, rep(0, I + K + P - 3L)),
       covariance = location_covariance_cpp(I, K, P, variances))
}

#' Exact conditional recovery of absolute additive effects
#' @param phi Reduced coordinates from reduce_effects.
#' @param I Number of items.
#' @param K Number of retained states.
#' @param P Number of respondents.
#' @param mean_item Population mean of item effects.
#' @param variances Item, state, person prior variances.
#' @param draw If TRUE draw the absolute effects; otherwise return their conditional means.
#' @return A list with item, state, person effects and the bivariate location moments.
#' @export
recover_effects <- function(phi, I, K, P, mean_item = 0, variances = c(1, 1, 1),
                            draw = TRUE) {
  if (!is.numeric(phi) || any(!is.finite(phi)) || !is.logical(draw) ||
      length(draw) != 1 || is.na(draw)) stop("Invalid reduced coordinates or draw flag.")
  # Validate dimensions and hyperparameters using the same public prior contract.
  .integer_scalar(I, "I"); .integer_scalar(K, "K"); .integer_scalar(P, "P")
  if (length(variances) != 3L || any(!is.finite(variances)) || any(variances <= 0) ||
      length(mean_item) != 1L || !is.finite(mean_item)) stop("Invalid location hyperparameters.")
  parts <- .location_parts(phi, I, K, P)
  vg <- variances[1]; vh <- variances[2]; vu <- variances[3]
  Q <- matrix(c(I / vg + K / vh, K / vh, K / vh, P / vu + K / vh), 2L)
  h <- c((I * mean_item - sum(parts$G)) / vg + sum(parts$theta + parts$S) / vh,
         -sum(parts$R) / vu + sum(parts$theta + parts$S) / vh)
  mean <- as.numeric(solve(Q, h))
  lambda <- if (draw) .rmvn_precision(Q, h) else mean
  list(item = lambda[1] + parts$G, person = lambda[2] + parts$R,
       state = parts$theta - sum(lambda) + parts$S,
       lambda = lambda, location_mean = mean, location_covariance = solve(Q))
}

.effect_from_phi <- function(phi, I, K, P) {
  part <- .location_parts(phi, I, K, P)
  list(item = part$G, state = part$theta + part$S, person = part$R)
}

.update_location_hyper <- function(item, state, person, current, prior, fixed = rep(FALSE, 3L)) {
  values <- list(item - current$mean, state, person)
  current$variances <- vapply(seq_len(3), function(j) {
    if (fixed[j]) return(current$variances[j])
    .rinvgamma(prior$shape[j] + length(values[[j]]) / 2,
               prior$rate[j] + sum(values[[j]]^2) / 2)
  }, numeric(1))
  if (!fixed[1]) {
    variance <- 1 / (1 / prior$mean_sd^2 + length(item) / current$variances[1])
    mean <- variance * (prior$mean / prior$mean_sd^2 + sum(item) / current$variances[1])
    current$mean <- rnorm(1, mean, sqrt(variance))
  }
  current
}

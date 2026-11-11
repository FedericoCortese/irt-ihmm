# PSBP probabilities and Gaussian probit updates adapted from psbpHMM.
# Original sources and differences are recorded in UPSTREAM.md.
.new_transition <- function(K, q, mode, m_alpha_prior = list(mean = 0, sd = 1)) {
  L <- if (mode == "infinite") K else K - 1L
  alpha <- matrix(rnorm(K * L), K, L)
  if (L) alpha[cbind(seq_len(L), seq_len(L))] <- rnorm(L, m_alpha_prior$mean)
  list(K = K, mode = mode, alpha = alpha, alpha0 = rnorm(L),
       beta = matrix(rnorm(q * L), q, L), sigma_alpha2 = 1,
       mean_diagonal = m_alpha_prior$mean, var_diagonal = 1)
}

.expand_transition <- function(tr, K) {
  if (tr$mode != "infinite" || K <= tr$K) return(tr)
  old <- tr$K; q <- nrow(tr$beta)
  alpha <- matrix(rnorm(K * K, sd = sqrt(tr$sigma_alpha2)), K, K)
  diag(alpha) <- rnorm(K, tr$mean_diagonal, sqrt(tr$var_diagonal))
  alpha[seq_len(old), seq_len(old)] <- tr$alpha
  beta <- matrix(rnorm(q * K), q, K)
  beta[, seq_len(old)] <- tr$beta
  alpha0 <- rnorm(K, sd = sqrt(tr$sigma_alpha2))
  alpha0[seq_len(old)] <- tr$alpha0
  tr$alpha <- alpha; tr$alpha0 <- alpha0; tr$beta <- beta; tr$K <- K
  tr
}

.trim_transition <- function(tr, K) {
  if (tr$mode == "finite" || K == tr$K) return(tr)
  idx <- seq_len(K)
  tr$alpha <- tr$alpha[idx, idx, drop = FALSE]
  tr$alpha0 <- tr$alpha0[idx]
  tr$beta <- tr$beta[, idx, drop = FALSE]
  tr$K <- K
  tr
}

.transition_probs <- function(tr, z) {
  out <- transition_probabilities_cpp(tr, z)
  if (!is.null(tr$omega) && tr$omega != 0) {
    out$transition <- lapply(out$transition, function(row) {
      row <- (1 - tr$omega) * row
      diag(row) <- diag(row) + tr$omega
      row
    })
    # Column 1 is the initial-distribution tail, which is unchanged.
    out$residual[, -1L] <- (1 - tr$omega) * out$residual[, -1L, drop = FALSE]
  }
  out
}

.draw_persistence <- function(tr, C, z, prior) {
  # Initial observations always belong to the probit component.
  active <- matrix(TRUE, nrow(C), ncol(C))
  forced <- 0L
  if (ncol(C) > 1L) {
    same <- which(C[, -1L, drop = FALSE] == C[, -ncol(C), drop = FALSE], arr.ind = TRUE)
    if (nrow(same)) {
      base <- transition_probabilities_cpp(tr, z)$transition
      probability <- vapply(seq_len(nrow(same)), function(j) {
        p <- same[j, 1]; s <- C[p, same[j, 2]]
        tr$omega / (tr$omega + (1 - tr$omega) * base[[p]][s, s])
      }, numeric(1))
      H <- rbinom(length(probability), 1L, probability)
      active[cbind(same[, 1], same[, 2] + 1L)] <- H == 0L
      forced <- sum(H)
    }
  }
  list(active = active, omega = rbeta(1L, prior[1] + forced,
    prior[2] + nrow(C) * (ncol(C) - 1L) - forced))
}

.draw_m_alpha <- function(diagonal, var_diagonal, prior) {
  # Scalar reference conditional for validation; production uses the compiled block.
  variance <- 1 / (1 / prior$sd^2 + length(diagonal) / var_diagonal)
  # Keep the original arithmetic for the default prior to preserve seeded draws.
  mean <- if (prior$mean == 0 && prior$sd == 1)
    variance * sum(diagonal) / var_diagonal else
    variance * (prior$mean / prior$sd^2 + sum(diagonal) / var_diagonal)
  rnorm(1, mean, sqrt(variance))
}

.update_transition <- function(tr, C, z, m_alpha_prior = list(mean = 0, sd = 1),
                               omega_prior = c(1, 1)) {
  if (is.null(tr$omega))
    return(transition_update_cpp(tr, C, z, m_alpha_prior$mean, m_alpha_prior$sd))
  auxiliary <- .draw_persistence(tr, C, z, omega_prior)
  tr <- transition_update_cpp(tr, C, z, m_alpha_prior$mean, m_alpha_prior$sd,
                              auxiliary$active)
  tr$omega <- auxiliary$omega
  tr
}

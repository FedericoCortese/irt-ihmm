#' Validate item responses in presentation order
#' @param x Numeric matrix of binary responses; respondents in rows and actual
#'   presentation positions in columns (common item order when item_index is NULL).
#' @param y Numeric matrix of log-response times, with the same dimensions as x.
#' @param z Optional numeric respondent-covariate matrix, without an intercept.
#' @param item_index Optional P by I integer matrix mapping presentation positions
#'   to common item indices 1:I. Each row must be a permutation of 1:I.
#' @return An object of class `irtihmm_data`.
#' @export
irtihmm_data <- function(x, y, z = NULL, item_index = NULL) {
  if (!is.matrix(x) || !is.numeric(x) || !length(x) ||
      anyNA(x) || any(!x %in% c(0, 1)))
    stop("x must be a complete numeric binary matrix.")
  if (!is.matrix(y) || !is.numeric(y) || !identical(dim(x), dim(y)) ||
      any(!is.finite(y)))
    stop("y must be a complete finite matrix of log-times matching x.")
  P <- nrow(x); I <- ncol(x)
  if (is.null(z)) z <- matrix(numeric(), P, 0L)
  if (!is.matrix(z) || !is.numeric(z) || nrow(z) != P ||
      any(!is.finite(z))) stop("z must be a finite numeric matrix with P rows.")
  if (ncol(z) && any(apply(z, 2, function(v) length(unique(v)) == 1L)))
    stop("z must not contain an intercept or constant covariate column.")
  out <- list(x = x, y = y, z = z, P = P, I = I)
  if (!is.null(item_index)) {
    if (!is.matrix(item_index) || !is.numeric(item_index) ||
        !identical(dim(item_index), dim(x)) || any(!is.finite(item_index)) ||
        any(item_index != floor(item_index)) ||
        any(!apply(item_index, 1L, function(v) setequal(v, seq_len(I)))))
      stop("item_index rows must be permutations of 1:I matching x dimensions.")
    storage.mode(item_index) <- "integer"
    out$item_index <- item_index
  }
  structure(out, class = "irtihmm_data")
}

.item_indices <- function(data) {
  if (!is.null(data$item_index)) data$item_index else
    matrix(rep(seq_len(data$I), each = data$P), data$P, data$I)
}

.item_predictor <- function(data, person, item) {
  if (is.null(data$item_index)) return(outer(person, item, "+"))
  matrix(item[data$item_index], data$P, data$I) + person
}

.item_average <- function(data, value) {
  if (is.null(data$item_index)) return(colMeans(value))
  as.numeric(rowsum(as.numeric(value), as.numeric(data$item_index), reorder = TRUE)) / data$P
}

#' Proper default priors for the research sampler
#' @param m_alpha_mean Mean of the Gaussian prior for the mean diagonal
#'   transition intercept m_alpha. Defaults to 0.
#' @param m_alpha_sd Standard deviation (not variance) of that prior.
#'   Must be positive; defaults to 1. This controls diagonal sticks, not a
#'   direct prior on the total number of occupied states.
#' @return A named list of location, log-scale and copula hyperpriors. Location
#'   group variances have inverse-gamma(shape, rate) densities proportional to
#'   v^(-shape-1) exp(-rate/v). Transition hyperpriors follow Hoskovec:
#'   off-diagonal and diagonal variances IG(1,1), diagonal mean N(0,1), and
#'   regression coefficients N(0,I).
#' @export
irtihmm_priors <- function(m_alpha_mean = 0, m_alpha_sd = 1) {
  .validate_m_alpha_prior(list(mean = m_alpha_mean, sd = m_alpha_sd))
  location <- list(mean = 0, mean_sd = 2, shape = rep(2, 3), rate = rep(1, 3))
  list(accuracy = location, time = location,
       log_scale = list(mean = 0, mean_sd = 1, shape = 2, rate = 1),
       rho_shape = 2,
       m_alpha = list(mean = m_alpha_mean, sd = m_alpha_sd))
}

.validate_m_alpha_prior <- function(prior) {
  if (!is.list(prior) || !setequal(names(prior), c("mean", "sd")) ||
      !is.numeric(prior$mean) || length(prior$mean) != 1L || !is.finite(prior$mean) ||
      !is.numeric(prior$sd) || length(prior$sd) != 1L || !is.finite(prior$sd) ||
      prior$sd <= 0 || !is.finite(prior$sd^2) || prior$sd^2 == 0)
    stop("m_alpha requires a finite numeric mean and a positive finite sd.")
  prior
}

.validate_priors <- function(priors) {
  # Older saved prior objects used the implicit N(0,1) diagonal-mean prior.
  if (is.list(priors) && !"m_alpha" %in% names(priors))
    priors$m_alpha <- list(mean = 0, sd = 1)
  if (!is.list(priors) || !setequal(names(priors), names(irtihmm_priors())))
    stop("priors must have accuracy, time, log_scale, rho_shape and m_alpha components.")
  .validate_m_alpha_prior(priors$m_alpha)
  for (name in c("accuracy", "time")) {
    p <- priors[[name]]
    if (!is.list(p) || !setequal(names(p), c("mean", "mean_sd", "shape", "rate")) ||
        length(p$mean) != 1 || !is.finite(p$mean) || length(p$mean_sd) != 1 ||
        !is.finite(p$mean_sd) || p$mean_sd <= 0 || length(p$shape) != 3 ||
        length(p$rate) != 3 || any(!is.finite(c(p$shape, p$rate))) ||
        any(c(p$shape, p$rate) <= 0)) stop("Invalid ", name, " location prior.")
  }
  p <- priors$log_scale
  if (!is.list(p) || !setequal(names(p), c("mean", "mean_sd", "shape", "rate")) ||
      any(lengths(p) != 1) || any(!is.finite(unlist(p))) ||
      any(unlist(p[c("mean_sd", "shape", "rate")]) <= 0) ||
      length(priors$rho_shape) != 1 || !is.finite(priors$rho_shape) || priors$rho_shape <= 0)
    stop("Invalid log-scale or copula prior.")
  priors
}

.initial_hyper <- function(priors) {
  list(accuracy = list(mean = priors$accuracy$mean, variances = rep(1, 3)),
       time = list(mean = priors$time$mean, variances = rep(1, 3)),
       scale = list(mean = priors$log_scale$mean, variance = 1),
       rho_shape = priors$rho_shape)
}

#' Gaussian-copula binary and normal joint emission
#' @param x Binary response vector.
#' @param y Log-response-time vector.
#' @param accuracy Accuracy predictor on the logit scale.
#' @param time_mean Normal mean of log-response time.
#' @param sigma Positive standard deviation of log-response time.
#' @param rho Gaussian-copula correlation, strictly between -1 and 1.
#' @param log Return log densities if TRUE.
#' @return Vector of joint densities or log densities.
#' @export
irtihmm_emission <- function(x, y, accuracy, time_mean, sigma, rho = 0,
                             log = TRUE) {
  args <- list(x, y, accuracy, time_mean, sigma, rho)
  n <- max(lengths(args))
  if (!n || any(!lengths(args) %in% c(1L, n)) ||
      any(!vapply(args, is.numeric, logical(1))) ||
      any(!is.finite(unlist(args)))) stop("Emission arguments must be finite and conformable.")
  args <- lapply(args, rep_len, length.out = n)
  x <- args[[1]]; y <- args[[2]]; accuracy <- args[[3]]
  time_mean <- args[[4]]; sigma <- args[[5]]; rho <- args[[6]]
  if (any(!x %in% c(0, 1)) || any(sigma <= 0) || any(abs(rho) >= 1))
    stop("Require binary x, sigma > 0 and abs(rho) < 1.")
  ans <- emission_density_cpp(x, y, accuracy, time_mean, sigma, rho)
  if (log) ans else exp(ans)
}

.emission_array <- function(data, emission, states = seq_along(emission$a)) {
  emission_array_cpp(data, emission, states)
}

.emission_path_loglik <- function(data, emission, C) {
  context <- emission_context_cpp(data, emission, C, TRUE)
  emission_margin_loglik_cpp(context, emission$nu, emission$mu, emission$eta)
}

#' Simulate common-order binary responses and log-response times
#' @param P Number of respondents.
#' @param I Number of items.
#' @param K Number of states in the finite data-generating model.
#' @param z Optional respondent covariates, without an intercept.
#' @param copula Whether to generate nonzero Gaussian-copula correlations.
#' @param seed Optional integer random seed; see set.seed.
#' @param parameters Optional named emission overrides (b,a,d,nu,mu,eta,log_sigma,rho).
#' @param beta Optional covariate effects on stochastic transition sticks: a
#'   numeric matrix with ncol(z) rows and K-1 columns. Shared with the initial
#'   distribution. NULL retains the small random coefficients used by default.
#' @param transition_matrix Optional positive K by K finite generating transition
#'   matrix with unit row sums at the zero covariate profile. Its probabilities
#'   are converted to probit stick intercepts. NULL keeps the default persistent
#'   intercepts. Initial-distribution intercepts and coefficient settings are unchanged.
#' @return A list containing validated data, true state paths, emissions and transitions.
#' @export
simulate_irtihmm <- function(P = 30L, I = 20L, K = 3L, z = NULL,
                             copula = TRUE, seed = NULL, parameters = list(),
                             beta = NULL, transition_matrix = NULL) {
  P <- .integer_scalar(P, "P"); I <- .integer_scalar(I, "I")
  K <- .integer_scalar(K, "K")
  if (!is.null(seed)) set.seed(seed)
  data <- irtihmm_data(matrix(0, P, I), matrix(0, P, I), z)
  state_sequence <- if (K == 1L) 0 else seq(-1, 1, length.out = K)
  emission <- list(b = seq(-1, 1, length.out = I), a = state_sequence,
    d = rnorm(P, sd = .7), nu = log(20) + seq(-.3, .3, length.out = I),
    mu = .5 * state_sequence, eta = rnorm(P, sd = .2),
    log_sigma = rep(log(.3), K), rho = if (copula) .4 * state_sequence else rep(0, K))
  if (!is.list(parameters) || any(!names(parameters) %in% names(emission)))
    stop("Unknown emission parameter override.")
  emission <- utils::modifyList(emission, parameters)
  lengths_expected <- c(b = I, a = K, d = P, nu = I, mu = K, eta = P, log_sigma = K, rho = K)
  if (any(lengths(emission) != lengths_expected) || any(!is.finite(unlist(emission))) ||
      any(abs(emission$rho) >= 1) || any(!is.finite(exp(emission$log_sigma))) ||
      any(exp(emission$log_sigma) <= 0)) stop("Invalid emission parameter dimensions or values.")
  if (!copula && any(emission$rho != 0)) stop("Nonzero rho requires copula=TRUE.")
  tr <- .new_transition(K, ncol(data$z), "finite")
  if (K > 1L) {
    tr$alpha[,] <- -1
    tr$alpha[cbind(seq_len(K - 1L), seq_len(K - 1L))] <- 2
    tr$alpha0[] <- -.5
    tr$beta <- tr$beta * .2
  }
  if (!is.null(beta)) {
    if (!is.matrix(beta) || !is.numeric(beta) || any(!is.finite(beta)) ||
        !identical(dim(beta), c(ncol(data$z), K - 1L)))
      stop("beta must be a finite matrix with ncol(z) rows and K-1 columns.")
    tr$beta <- beta
  }
  if (!is.null(transition_matrix)) {
    if (!is.matrix(transition_matrix) || !is.numeric(transition_matrix) ||
        !identical(dim(transition_matrix), c(K, K)) ||
        any(!is.finite(transition_matrix)) || any(transition_matrix <= 0) ||
        any(abs(rowSums(transition_matrix) - 1) > 1e-10))
      stop("transition_matrix must be a positive K by K matrix with unit row sums.")
    remaining <- rep(1, K)
    if (K > 1L) for (k in seq_len(K - 1L)) {
      fraction <- transition_matrix[, k] / remaining
      if (any(fraction <= 0 | fraction >= 1))
        stop("transition_matrix must yield finite probit stick intercepts.")
      tr$alpha[, k] <- qnorm(fraction)
      remaining <- remaining - transition_matrix[, k]
    }
  }
  probs <- .transition_probs(tr, data$z)
  C <- matrix(0L, P, I)
  for (p in seq_len(P)) for (i in seq_len(I)) {
    weights <- if (i == 1L) probs$initial[p, ] else probs$transition[[p]][C[p, i - 1L], ]
    C[p, i] <- sample.int(K, 1, prob = weights)
    s <- C[p, i]
    predictor <- emission$b[i] + emission$a[s] + emission$d[p]
    log_notq <- -.softplus(predictor); log_q <- -.softplus(-predictor)
    threshold <- if (predictor >= 0) qnorm(log_notq, log.p = TRUE) else
      qnorm(log_q, lower.tail = FALSE, log.p = TRUE)
    time_residual <- rnorm(1)
    accuracy_residual <- emission$rho[s] * time_residual + sqrt(1 - emission$rho[s]^2) * rnorm(1)
    data$x[p, i] <- as.numeric(accuracy_residual > threshold)
    data$y[p, i] <- emission$nu[i] + emission$mu[s] + emission$eta[p] +
      exp(emission$log_sigma[s]) * time_residual
  }
  list(data = data, states = C, emission = emission, transition = tr,
       probabilities = probs)
}

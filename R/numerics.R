.softplus <- function(x) pmax(x, 0) + log1p(exp(-abs(x)))

.logsumexp <- function(x) {
  a <- max(x)
  if (!is.finite(a)) return(a)
  a + log(sum(exp(x - a)))
}

.draw_log <- function(x) {
  norm <- .logsumexp(x)
  if (!is.finite(norm)) stop("No finite probability in state recursion.")
  sample.int(length(x), 1L, prob = exp(x - norm))
}

.rmvn_precision <- function(precision, information) {
  U <- chol(precision)
  as.numeric(backsolve(U, forwardsolve(t(U), information)) +
               backsolve(U, rnorm(nrow(precision))))
}

.rinvgamma <- function(shape, rate) 1 / rgamma(1L, shape, rate = rate)

# Elliptical slice update under an exact Gaussian prior (no tuning parameter).
.elliptical_slice <- function(current, mean, covariance, loglik, direction = NULL) {
  centered <- current - mean
  if (is.null(direction))
    direction <- as.numeric(t(chol(covariance)) %*% rnorm(length(current)))
  log_height <- loglik(current) - rexp(1)
  angle <- runif(1, 0, 2 * pi)
  lower <- angle - 2 * pi
  upper <- angle
  for (j in seq_len(10000L)) {
    proposal <- mean + centered * cos(angle) + direction * sin(angle)
    value <- loglik(proposal)
    if (is.finite(value) && value >= log_height)
      return(list(value = proposal, evaluations = j))
    if (angle < 0) lower <- angle else upper <- angle
    angle <- runif(1, lower, upper)
  }
  stop("Elliptical slice bracket did not close; chain aborted.")
}

# Neal-style stepping-out slice, with a randomized finite expansion budget.
.slice_scalar <- function(current, logdensity, width = 1, max_steps = 50L) {
  height <- logdensity(current) - rexp(1)
  if (!is.finite(height)) stop("Non-finite current scalar log density.")
  left <- current - width * runif(1)
  right <- left + width
  j <- sample.int(max_steps, 1L) - 1L
  k <- max_steps - 1L - j
  while (j > 0L && logdensity(left) > height) {
    left <- left - width; j <- j - 1L
  }
  while (k > 0L && logdensity(right) > height) {
    right <- right + width; k <- k - 1L
  }
  for (attempt in seq_len(10000L)) {
    proposal <- runif(1, left, right)
    if (logdensity(proposal) >= height) return(proposal)
    if (proposal < current) left <- proposal else right <- proposal
  }
  stop("Scalar slice bracket did not close; chain aborted.")
}

.integer_scalar <- function(x, name, minimum = 1L) {
  if (length(x) != 1L || !is.numeric(x) || !is.finite(x) ||
      x != floor(x) || x < minimum)
    stop(name, " must be an integer >= ", minimum, ".")
  as.integer(x)
}

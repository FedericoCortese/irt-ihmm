test_that("Gaussian recovery agrees with conditioning the original joint normal", {
  set.seed(14)
  I <- 3; K <- 2; P <- 4; v <- c(.4, 1.2, .8); m <- .7
  original <- rnorm(I + K + P)
  T <- .location_transform(I, K, P)
  phi <- as.numeric(T %*% original)
  out <- recover_effects(phi, I, K, P, m, v, draw = FALSE)
  A <- rbind(T, diag(I + K + P)[c(1, I + K + 1), ])
  expect_equal(abs(det(A)), 1)
  fullcov <- A %*% diag(rep(v, c(I, K, P))) %*% t(A)
  fullmean <- A %*% c(rep(m, I), rep(0, K + P))
  n <- length(phi)
  cross <- fullcov[n + 1:2, seq_len(n), drop = FALSE]
  conditional_mean <- fullmean[n + 1:2] + cross %*%
    solve(fullcov[seq_len(n), seq_len(n)], phi - fullmean[seq_len(n)])
  conditional_cov <- fullcov[n + 1:2, n + 1:2] - cross %*%
    solve(fullcov[seq_len(n), seq_len(n)], t(cross))
  expect_equal(out$location_mean, as.numeric(conditional_mean))
  expect_equal(out$location_covariance, conditional_cov)
  expect_equal(reduce_effects(out$item, out$state, out$person), phi)
  recovered <- replicate(4000, recover_effects(phi, I, K, P, m, v)$lambda)
  expect_equal(rowMeans(recovered), out$location_mean, tolerance = .04)
  expect_equal(cov(t(recovered)), out$location_covariance, tolerance = .04)
})

test_that("location shifts leave joint copula emissions unchanged", {
  sim <- simulate_irtihmm(4, 5, 2, seed = 4)
  e <- sim$emission; shifted <- e
  shifted$b <- e$b + 3; shifted$a <- e$a - 1; shifted$d <- e$d - 2
  shifted$nu <- e$nu - 4; shifted$mu <- e$mu + 1; shifted$eta <- e$eta + 3
  expect_equal(.emission_array(sim$data, e),
               .emission_array(sim$data, shifted), tolerance = 1e-12)
})

test_that("independent time regression reproduces original additive predictors", {
  for (dims in list(c(3, 2, 4), c(1, 1, 1), c(4, 1, 2))) {
    I <- dims[1]; K <- dims[2]; P <- dims[3]
    set.seed(1); g <- rnorm(I); h <- rnorm(K); u <- rnorm(P)
    C <- matrix(sample.int(K, P * I, TRUE), P, I)
    phi <- reduce_effects(g, h, u)
    design <- .time_design(C, I, K, P)
    expect_equal(as.numeric(design %*% phi),
                 as.numeric(outer(u, g, "+") + matrix(h[C], P)))
  }
})

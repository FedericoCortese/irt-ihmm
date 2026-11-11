test_that("the copula density normalizes to both requested margins", {
  q <- .7; mean <- .4; sigma <- .8
  for (rho in c(-.6, 0, .6)) {
    for (x in 0:1) {
      mass <- integrate(function(y) irtihmm_emission(x, y, qlogis(q), mean, sigma,
        rho, log = FALSE), -Inf, Inf)$value
      expect_equal(mass, if (x == 1) q else 1 - q, tolerance = 1e-7)
    }
    y <- seq(-2, 2, length.out = 20)
    total <- irtihmm_emission(0, y, qlogis(q), mean, sigma, rho, FALSE) +
      irtihmm_emission(1, y, qlogis(q), mean, sigma, rho, FALSE)
    expect_equal(total, dnorm(y, mean, sigma), tolerance = 1e-13)
  }
})

test_that("independence is exactly Bernoulli times normal, including logit tails", {
  predictor <- c(-1000, -40, 0, 40, 1000)
  for (x in 0:1) {
    got <- irtihmm_emission(x, 0, predictor, 0, 1, 0)
    expected <- dnorm(0, log = TRUE) + if (x == 1)
      -.softplus(-predictor) else -.softplus(predictor)
    expect_equal(got, expected)
    expect_true(all(is.finite(irtihmm_emission(x, 0, predictor, 0, 1, .5))))
  }
})

test_that("positive rho associates correctness with longer response times", {
  y <- c(-1, 1)
  conditional <- irtihmm_emission(1, y, 0, 0, 1, .6, FALSE) / dnorm(y)
  expect_gt(conditional[2], conditional[1])
})

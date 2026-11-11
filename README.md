# irt-ihmm

`irt-ihmm` is an R implementation of a Bayesian item response model with an
infinite hidden Markov process. It jointly models response accuracy and response
time, allowing a respondent to move between latent regimes over the course of
an assessment.

For respondent $p$ and presentation position $i$, let $j_{pi}$ be the item
presented, $X_{pi}$ the binary response, and $Y_{pi}$ the log response time.
Conditional on regime $C_{pi}=s$, the marginal models are

$$
\operatorname{logit}\Pr(X_{pi}=1\mid C_{pi}=s)
  = b_{j_{pi}} + d_p + a_s,
$$

$$
Y_{pi}\mid C_{pi}=s
  \sim \mathcal N(\nu_{j_{pi}} + \eta_p + \mu_s,\sigma_s^2).
$$

The item effects $b_j$ and $\nu_j$ describe ease and time intensity. The
respondent effects $d_p$ and $\eta_p$ describe ability and log-time offsets.
Each regime has its own accuracy shift $a_s$, log-time shift $\mu_s$, and
log-time standard deviation $\sigma_s$. An optional Gaussian copula with
regime-specific correlation $\rho_s$ links the two margins.

A covariate-dependent probit stick-breaking process specifies the initial
distribution and transitions among a countable set of shared regimes. Each
respondent has a separate trajectory in presentation order. An optional
persistence parameter $\omega$ adds probability to staying in the current
regime. Inference uses beam sampling for the infinite model, with a finite
model available for comparison. The sampler also supports location-reduced
updates with exact Gaussian recovery under the original hierarchical priors.
See the [identification and inference notes](irt_ihmm_identification_gustafson.tex)
and the [upstream adaptation record](inst/UPSTREAM.md).

## Quick start

Use R 4.1 or later with a C++17 compiler. Install the R dependencies once, then
run the example from the repository root. The repository name is `irt-ihmm`
and the R package name is `irtihmm`.

```r
install.packages(c("Rcpp", "RcppArmadillo"))
source("load.R")

set.seed(18)
z <- cbind(group = rbinom(12, 1, 0.5), score = rnorm(12))
sim <- simulate_irtihmm(P = 12, I = 10, K = 2, z = z, seed = 19)

fit <- fit_irtihmm(sim$data, mode = "infinite", K = 2,
                   iterations = 300, burnin = 150, seed = 23)
plot_fit(fit, simulation = sim)
table(summary(fit)$occupied_states)
```

This short run illustrates the interface. In the simulator, `K` sets the number
of generating regimes. In the infinite fit, it sets the initial state prefix.
The number of occupied regimes can vary across draws.

Use `irtihmm_data(x, y, z, item_index = map)` to prepare your own data.

| Input | Format |
| --- | --- |
| `x` | Respondent-by-position matrix of binary responses. |
| `y` | Matching matrix of log response times. |
| `z` | Optional respondent covariate matrix in the same row order, without an intercept column. |
| `map` | Optional respondent-by-position matrix of item indices. Each row is a permutation of `1:I`. Supply it when respondents see the same items in different orders. |

The fitting interface expects complete responses and finite log response times.
The [two-regime demo](demo.R) and [three-regime demo](demo_three_regimes.R) show
the full simulation workflow. The [ROAR preparation script](prepare_roar_study2.R)
and [fitting script](fit_roar_study2.R) illustrate an application with separately
obtained data. The ROAR fitting script also requires the `posterior` R package.

Transition code is adapted from Lauren Hoskovec's [psbpHMM](https://github.com/lvhoskovec/psbpHMM).
The project is distributed under GPL version 2 or later. [License text](COPYING).

## Post-processing

SALSO uses the sampled state partitions to find a representative grouping of all
responses, using variation of information as the default loss. CHIPS looks for a
subset of responses whose grouping holds jointly with a specified posterior
probability. That subset anchors comparisons of state-specific parameters
across compatible draws. [SALSO paper](https://doi.org/10.1080/10618600.2022.2069779)
· [CHIPS paper](https://arxiv.org/abs/2511.16040).

```r
install.packages("salso")  # Version 0.3.80 or later
post <- postprocess_irtihmm(fit, credible = 0.90, seed = 24)
plot_postprocess(post)
```

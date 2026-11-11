# Upstream source and adaptation record

Original author: Lauren Hoskovec and psbpHMM contributors.
Repository: https://github.com/lvhoskovec/psbpHMM
Commit: 1240b4da72a504dd5c97112f0b06f4b4a40f602a
Original package version in DESCRIPTION: 1.0.
Original license: GPL (>= 2).

Reference: Hoskovec et al. (2023), *Infinite hidden Markov models for multiple
multivariate time series with missing data*, Biometrics 79, 2592--2604,
https://doi.org/10.1111/biom.13715.

## Unmodified source snapshot

The development repository's `vendor/psbpHMM/` contains the upstream R and
C++ sources and root DESCRIPTION, NAMESPACE and README.md. The files are
byte-for-byte copies, checked against `vendor/manifest.json`. Large datasets,
cached vignette outputs and simulation results are not imported. The vendor
snapshot is excluded from the installed package; derived compiled source is
included in `src/` and retains this provenance record.

## Imported compiled reference functions

From upstream `src/tester.cpp`:

* `updatePi` -> `hoskovec_pi`: body unchanged. The original R adapter transposed the
  origin-by-destination intercept matrix to upstream's destination-by-origin
  representation, evaluated two covariate-identical positions, and reused the
  initial vector and transition matrix for the common-order sequence. The
  adapter recomputed tail masses with log survival probabilities to avoid
  subtractive cancellation.
* `up_ajk` -> `hoskovec_alpha`: adds `n_sticks`, using it for the destination
  loop and output-vector size; this supports K-1 stochastic sticks in the
  finite reference model. Replaces Boolean `&` by short-circuit `&&`. The
  Gaussian update formulas, loops and posterior sufficient statistics are
  otherwise retained. No original helper functions are needed by these kernels.

The package uses C++17 and generated Rcpp registration. These two imported
functions remain compiled and tested as references, but the production hot path
now uses `src/sampler_kernels.cpp`. It pools the same reached-stick Gaussian
sufficient statistics in one compiled pass, avoiding the original nested R lists
and repeated scans per origin/stick. Probabilities are evaluated directly with
normal log survival functions, including the original orientation and finite
terminal residual mass.

Probit auxiliaries now use a new exact one-sided normal rejection sampler:
normal rejection when the mean is nonnegative and an exponential tail envelope
otherwise. The original custom upstream rejection helpers are not imported.
Tests check probability-integral transforms at central and tail means and the
analytical transition-intercept conditional means. This replaces the R-level
truncnorm calls without changing the augmentation distribution. Manual m_alpha
mean/SD tuning enters the same Gaussian full conditional.

## R transition blocks adapted from mciHMM.R

`R/transitions.R` follows the same truncated-normal stick augmentation,
Gaussian initial-intercept update, Gaussian beta update, and inverse-gamma
variance updates as upstream `R/mciHMM.R`. Changes:

* beta updates pair covariates, latent variables and previous-state intercepts
  by explicit (respondent,position) indices;
* the initial intercepts contribute to the shared off-diagonal variance
  conditional, as required by their shared N(0,sigma_alpha^2) prior;
* a finite model has no stochastic terminal stick and therefore no probit
  latent variable at that terminal stick;
* a zero-column covariate matrix and single-state models are supported;
* no individual random regression effects are added, matching the proposal.

## State sampling and emissions

The slice construction and indicator-restricted forward recursion follow
upstream `mciHMM.R`, `upStateList_lapply.R` and `upZ` in `tester.cpp`.
The new `src/forward_backward.cpp` accepts generic supplied log emissions
and includes an explicit backward draw of the whole trajectory. It does not
reuse the upstream conjugate multivariate-Gaussian new-state predictive density:
that density is inappropriate for binary/continuous copula emissions.

All states permitted by the slice bounds are explicitly instantiated with
their conditional-prior parameters. Expansion checks every current origin
row and the initial row. After the path draw the unused suffix is integrated
out, and slices are discarded before unconstrained transition updates. State
indices are never automatically compacted. These are intentional changes
from the original orchestration, not claims that the two samplers have
identical paths or numerically identical output.

The Rasch/time/copula likelihood, simulator, reduced Gaussian priors,
elliptical slice emission updates, exact location recovery and complete-data
validation interface are new model-specific code. The emitted prior remains
the original independent normal hierarchy, with no independent convenience
prior on contrasts. The original Gaussian exposure/imputation model and
its dependencies are not incorporated.

This initial validation reproduces the imported probability and intercept
kernels, not the original paper's complete exposure-data analysis. Several
optional imputation dependencies needed by that original analysis are not
installed in the development environment.

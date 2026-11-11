# Initial validation record

Date: 2026-10-19. Package version: 0.0.1.
Environment: R 4.5.2, macOS arm64, Apple clang 21, C++17.

This is the historical validation record for the initial version. Since
2026-10-20 the simplified root `demo.R` uses a different simulation with a
binary covariate and one 2,000-iteration infinite-model fit. The results below
refer to the earlier experiment, not to that current demo.

## Provenance

The 38 preserved psbpHMM source/metadata files match their recorded SHA-256
hashes at commit 1240b4da72a504dd5c97112f0b06f4b4a40f602a. The two compiled
upstream kernels also match the documented source transformations exactly.
Run `python3 tools/verify_upstream.py` to reproduce these checks.

## Numerical and software checks

* Numerical integration of each binary copula cell recovers its Bernoulli
  probability; summing cells recovers the normal marginal density.
* At rho=0, emissions reduce to the independent Bernoulli/normal product.
* Extreme positive/negative logit predictors remain finite on the log scale
  in the tested cases, including nonzero correlation.
* The additive coordinate map preserves predictors, has unit absolute
  Jacobian, and its bivariate Gaussian recovery matches an independent joint
  multivariate-normal conditioning calculation.
* The reduced time design reconstructs the original predictor, including
  singleton dimensions.
* Finite forward likelihoods equal explicit path enumeration. Monte Carlo
  joint trajectory probabilities agree with enumeration for both ordinary
  and slice-restricted paths.
* The imported PSBP routine agrees with a direct formula with nonzero
  covariates. The reused Gaussian intercept kernel agrees with its analytic
  conditional moments.
* Known-state, finite-state and infinite-state model fits complete with
  finite emission log likelihoods. Every recovered draw preserves its reduced
  coordinates; no-copula draws keep rho exactly zero.
* Slice residual bounds are satisfied before accepting infinite-state path
  draws. A deliberately insufficient max_states aborts instead of accepting
  a truncated path. Known trajectories with gaps preserve their stick labels.
* Invalid/missing observations, constant covariate columns and invalid prior
  settings are rejected.

The testthat suite passed. R CMD check --no-manual --no-vignettes returned
Status: OK (no errors, warnings or notes) using the documented local Makevars
override. Network access to package-index servers was unavailable during the
dependency lookup, but all declared dependencies were installed and checked.
The override addresses this host's absent Fortran linker paths and an
unsupported warning name in R system headers, not a statistical code change.

## Reproducible demonstration

`scripts/run_demo.R` generates 12 respondents by 12 items from a finite
two-state model with Gaussian-copula dependence, using seed 19. It runs 600
iterations and discards 300 in each of three configurations, with seeds
21, 22 and 23. The initial observed results were:

| Fit | Saved draws | Occupied states, saved draws | Largest saved prefix | Largest temporary prefix |
| --- | ---: | ---: | ---: | ---: |
| Known trajectories | 300 | 2 | 2 | 2 |
| Finite K=2, unknown trajectories | 300 | 2 | 2 | 2 |
| Infinite, unknown trajectories | 300 | 3--9 | 14 | 22 |

These are execution and tail-expansion checks. In particular, the infinite
demonstration does not recover two occupied regimes reliably in this short
run. Its larger state counts must not be presented as successful recovery
of the generating state count. Neither a conditional emission likelihood
increase nor the absence of software errors proves posterior convergence.

## Work still required for statistical validation

Longer independent chains, posterior calibration, state alignment and
parameter/trajectory recovery across multiple simulated datasets remain
necessary. In particular, evaluate the state-count behavior against the
transition and emission priors, posterior state occupancy, run lengths and
sensitivity to rare states. No validated ESS/R-hat, evidence estimate or
global identification theorem is supplied by the initial demonstration.

## Compiled-kernel validation and timing

The optimized backend is covered by independent checks of joint copula densities,
per-margin and per-state likelihood caches, sparse Gaussian sufficient statistics
against dense reference designs (including individual item orders), and structured
reduced-prior covariances against T diag(v) T'. Batched FFBS and beam slice draws
are compared to the original respondent loops with identical RNG seeds. The
one-sided normal sampler is checked through probability-integral transforms,
including mean -20. Transition intercept and covariate means and tuned m_alpha/
diagonal-variance conditionals are checked against analytic moments. The complete
suite, source provenance check and compiled-package registration/load checks pass.
Fixed emission blocks, both Gustafson modes, finite/infinite models and manual
m_alpha controls remain covered. CHIPS/SALSO code is unchanged.

Warm-call timings on this development host, excluding compilation:

| Operation | Earlier R orchestration | Compiled backend | Approximate speedup |
| --- | ---: | ---: | ---: |
| Synthetic transition update (30 x 80, K=6) | 22.2 ms | 0.25 ms | 89x |
| Synthetic emission update | 19.2 ms | 4.3 ms | 4.5x |
| Synthetic beam block | 74.5 ms | 25.5 ms | 2.9x |
| Synthetic short fit, seconds per iteration | 0.0730 | 0.0180 | 4.0x |
| ROAR transition update (84 x 252, saved prefix K=20) | 376 ms | 3.5 ms | 108x |
| ROAR emission update | 141 ms | 34.5 ms | 4.1x |
| ROAR beam block | 126 ms | 56.3 ms | 2.2x |

The ROAR measurements reuse a saved allocation/parameter configuration, not a
new long fit. Repeated stochastic update calls target the same conditionals but
consume RNG streams differently. The short-fit comparison uses 30 synthetic
iterations and is not an end-to-end ROAR runtime or an ESS-per-second comparison.
The timing resolution and host load limit the precision of very fast calls.
`tools/benchmark_sampler.R` measures the current kernels and optionally compares
a saved prechange timing/input object. CSV timings are under output/performance/.


## Direct CHIPS and shared partition-cache validation

Post-processing now retains standalone SALSO plus the CHIPS-SALSO completion,
but calls CHIPS directly at the requested threshold (0.95 by default) with four
search restarts. It no longer computes the full threshold-zero intermediate path
and then invokes a second CHIPS search. All retained draws and responses remain
included; no subsampling or post-processing thinning is used.

Tests cover simultaneous core probability, separation between distinct clusters,
label invariance, singleton and no-CHIPS cases, RNG restoration, selected-core
plots, blocked matrix packing, numeric/integer state matrices, cached versus
uncached summaries and stale/mismatched cache rejection. Per-chain caches share
pooled storage while retaining chain-local compatibility indices. The large
matrix/cache is not serialized in result objects.

A bounded synthetic comparison (150 draws, 600 responses, one search core)
measured 0.641 seconds for the earlier standalone SALSO + full CHIPS path +
threshold completion pattern with 32 restarts, versus 0.042 seconds for the new
four-restart direct workflow (about 15x in that example). The heuristically
selected core/partition need not be identical when reducing restarts. This
measurement is not a prediction of end-to-end ROAR runtime. No long application
chains were rerun for this post-processing change.

## Optional global persistence (2026-11-09)

The full test suite and a temporary-library package installation pass with the
optional omega mixture. Checks include normalized finite/infinite rows, unchanged
initial probabilities and tail, covariate probabilities under reordered selected
states, masked probit sufficient statistics (including covariate coefficients),
seeded default/OFF equivalence, auxiliary Gibbs moments against numerical
integration, and beam trajectory frequencies against an enumerated posterior.
The installed compiled wrapper was exercised with persistence enabled.

The existing two-regime toy was run with the same generated data, seeds, priors,
free emission effects, 5,000 iterations and 1,000 burn-in for both OFF and ON.
The older saved result fixed eta=0 and was therefore backed up rather than used
as the comparison baseline. All outputs remain in the existing demo directory.
With a Beta(1,1) omega prior, its posterior mean was 0.282, with a 95% equal-tail
interval [0.026, 0.555]. Mean occupied-state counts were 5.198 (OFF) and 5.698
(ON). Thus this single-chain comparison does not demonstrate fewer occupied
states. SALSO on the ON draws selected 3 states; CHIPS probability was 0.9505
and the core fraction was 0.235. These summaries do not establish convergence;
omega and the base diagonal sticks can both explain self-transitions.

## ROAR application with global persistence (2026-11-09)

Three independent chains completed 5,000 iterations each, discarding 1,000 and
retaining every subsequent draw. The analysis retained 84 complete-time subjects
and 252 items in their actual respondent-specific order. Omega had a Beta(1,1)
prior; the existing m_alpha ~ N(1,0.25^2), free emission effects, copula, covariate
scaling and Gustafson settings were preserved. The original neighboring run's
settings and data matched after excluding the new persistence fields.

Mean occupied counts OFF/ON by chain were 16.433/22.691, 17.444/21.015 and
16.744/25.916. The pooled sample mean of omega was 0.4893. Rank-normalized
split R-hat was 1.3626 for omega, 2.4215 for occupied-state count and 2.4719 for
conditional emission log-likelihood. These chains have not established convergence;
the pooled summaries are exploratory and do not demonstrate state-count reduction.
All chain fits are stored separately from the original model under persistence_on.

A source-repository integration test verifies omega summaries, R-hat/ESS export,
and the absence of any parameter PDFs when both trace flags are disabled. The
full test suite passed (the small SALSO fixture issues its expected maximum-cluster
search warning). The real run disables both per-chain and multichain parameter PDFs.

The full ROAR post-processing and graphical export completed successfully on all
12,000 retained allocations, with no parameter PDF files produced. SALSO selected
15 clusters; CHIPS had empirical joint probability 0.951083 and covered about 1.16%
of responses in 6 core clusters. SALSO warned that the maximum number of zealous
attempts had been reached. The warning is recorded in the post-processing cache
and fit report, and the application helper now retains search warnings on future
runs. The selected partitions are heuristic summaries, with the convergence and
small-core limitations above. Final graph readability was checked; the covariate
heatmap's width and left-axis spacing were adjusted and regenerated from saved
summaries, without rerunning any MCMC.

## ROAR fit-only rerun with Beta(3,1) (2026-11-09)

At the user's request, the previous output/roar_study2/fit working directory was
replaced rather than creating another experiment directory. Input data were
validated before replacement. The earlier ROAR artifact locations above are now
historical. Three fresh chains ran for 5,000 iterations each, with 1,000 burn-in,
thin=1, seeds 101/202/303, persistence enabled and omega ~ Beta(3,1). The existing
m_alpha ~ N(1,0.25^2), free emission effects, copula and Gustafson were unchanged.

Mean occupied-state counts by chain were 22.556, 24.740 and 21.880; mean omega
values were 0.4039, 0.4318 and 0.5716. R-hat was 2.1385 for omega, 1.3916 for
occupied-state count and 1.8119 for conditional emission log-likelihood. This fit
has not established convergence and cannot support definitive prior comparisons.

SALSO/CHIPS, partition packing and aligned summaries were skipped entirely. No
parameter PDFs were generated. Fit-only summaries, item/respondent plots, omega
summaries, per-chain reconstruction plots and R-hat/ESS were exported successfully.
After the completed replacement run, overwrite_fits was restored to FALSE so
future post-processing can reuse these saved chains. Set it to TRUE explicitly
for another fresh replacement experiment in the same working directory.

## SALSO/CHIPS on the saved Beta(3,1) ROAR fit at 90% (2026-11-09)

Post-processing was run directly on the three existing 4,000-draw chain files,
with no sampler invocation. All chain settings and data matched analysis_input.
File sizes, modification times and SHA-256 fingerprints were unchanged afterward.
The run used all 12,000 allocations, 4 search restarts, 3 cores and direct CHIPS
with threshold 0.90. SALSO selected a representative partition with 9 clusters;
CHIPS retained 285 responses (1.3464%) in 6 core clusters, with empirical joint
probability 10814/12000 = 0.9011667. The core-constrained SALSO completion had
11 clusters. These point summaries do not change the sampled occupied-state counts.
Parameter summaries retain 95% conditional marginal intervals.

The export completed successfully, including aligned state parameters, covariate
coefficients/probability contrasts and the existing item/respondent diagnostics.
Graphs were checked visually and no parameter PDFs were generated. SALSO's
maximum-zealous-attempts warning is retained in its cache and fit_report.txt.
The pre-existing MCMC convergence limitations remain; pooled summaries are exploratory.
The launcher now enables post-processing at 0.90 with overwrite_fits still FALSE.

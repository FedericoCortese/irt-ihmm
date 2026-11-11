# irtihmm 0.0.1

* Import pinned Hoskovec source and reuse PSBP probability/intercept kernels.
* Add common-order joint Rasch/log-time Gaussian-copula emissions.
* Preserve Gaussian location priors through exact reduction and recovery.
* Add known-state, finite-state FFBS and adaptive infinite beam inference.
* Add simulation, analytical tests, reproducible examples and provenance.
* Support fixed emission blocks through fixed_params, including fixed eta=0,
  with constrained updates and consistent infinite-state expansion.
* Support respondent-specific item permutations in likelihoods, Gaussian designs,
  initialization and item reconstruction, retaining chronological HMM paths.
* Add a reproducible three-chain ROAR XYZ application and multichain diagnostics.
* Expose manual Gaussian m_alpha prior mean/SD in both demos and the ROAR fit;
  defaults retain N(0,1), with separate output paths for nondefault priors.
* Compile transition sufficient statistics/auxiliaries, copula likelihood caches,
  beam slices and Gaussian design accumulation in one new kernel file.
* Batch FFBS paths, precompute finite transition logs and use exact structured
  Gaussian prior directions; retain CHIPS/SALSO and all existing model controls.
* Run CHIPS directly at the target probability with four search restarts by
  default; retain standalone SALSO and CHIPS-SALSO completion, without a full path.
* Pack partitions and test core compatibility in C++, sharing an in-memory cache
  and chain views across summaries while checking for stale draw/core pairs.
# Optional forced persistence

* Use one ROAR working directory with an explicit `overwrite_fits` switch for
  repeated fit experiments. Overwrite mode validates inputs first, then clears
  old results and partition caches before sampling fresh chains.

* Add a ROAR `post_processing` flag, off by default, to defer SALSO/CHIPS and
  aligned summaries while keeping fit diagnostics and unaligned effect plots.
  Re-enabling it reuses matching chain files; fit-only summaries preserve older
  partition caches and posterior summaries.

* Add `persistence` and `omega_prior` arguments to `fit_irtihmm`, with a global
  Beta-distributed omega and a self-transition mixture. The default remains off.
* Exclude forced transitions from both intercept and coefficient probit updates;
  use the mixture throughout finite and infinite state sampling and covariate plots.
* Include omega in parameter traces. Run the existing two-state toy with the flag
  on and keep outputs and a compatible previous baseline in the same directory.
* Enable omega for the ROAR application in a separate persistence-on output
  directory, with cache checks for the model flag and Beta prior. Export omega
  summaries, diagnostics and a compatible original-model comparison. Disable
  parameter PDFs through independent per-chain and multichain flags.


# Fixed emission blocks replace their priors with point masses. Scalars broadcast;
# state vectors are allowed only when the number of states is fixed.
.validate_fixed_params <- function(fixed_params, data, K, mode, copula) {
  allowed <- c(b = data$I, a = K, d = data$P, nu = data$I, mu = K,
               eta = data$P, log_sigma = K, rho = K)
  if (!is.list(fixed_params) || (length(fixed_params) &&
      (is.null(names(fixed_params)) || anyNA(names(fixed_params)) ||
       any(!nzchar(names(fixed_params))) || anyDuplicated(names(fixed_params)) ||
       any(!names(fixed_params) %in% names(allowed)))))
    stop("fixed_params must be a uniquely named list of emission parameters: ",
         paste(names(allowed), collapse = ", "), ".")
  for (name in names(fixed_params)) {
    value <- fixed_params[[name]]
    if (!is.numeric(value) || !is.null(dim(value)) || !length(value) ||
        any(!is.finite(value)) || !length(value) %in% c(1L, allowed[[name]]))
      stop("Invalid fixed_params value for ", name, ". Supply a finite scalar or full vector.")
    if (mode == "infinite" && name %in% c("a", "mu", "log_sigma", "rho") && length(value) != 1L)
      stop("Infinite-state fixed parameters require scalar state effects: ", name, ".")
    if (name == "rho" && (any(abs(value) >= 1) || (!copula && any(value != 0))))
      stop("Fixed rho must lie in (-1,1) and equal zero when copula=FALSE.")
    if (name == "log_sigma" && any(!is.finite(exp(value)) | exp(value) <= 0))
      stop("Fixed log_sigma must imply a finite positive standard deviation.")
  }
  fixed_params
}

.apply_fixed_params <- function(emission, fixed_params) {
  for (name in names(fixed_params))
    emission[[name]] <- rep(fixed_params[[name]], length.out = length(emission[[name]]))
  emission
}

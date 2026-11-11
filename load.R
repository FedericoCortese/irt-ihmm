# Development entry point: source("load.R") from RStudio or Rscript.
# No package installation and no changes to .libPaths() are required.
local({
  source_paths <- vapply(sys.frames(), function(frame) {
    path <- frame$ofile
    if (is.character(path) && length(path) == 1L) path else ""
  }, character(1))
  source_paths <- source_paths[nzchar(source_paths)]
  entry <- if (length(source_paths)) tail(source_paths, 1L) else "load.R"
  root <- dirname(normalizePath(entry, mustWork = TRUE))
  target <- parent.env(environment())

  load_sources <- function(root, target) {
    dependencies <- c("Rcpp", "RcppArmadillo")
    missing <- dependencies[!vapply(dependencies, requireNamespace, logical(1), quietly = TRUE)]
    if (length(missing))
      stop("Missing dependencies: ", paste(missing, collapse = ", "),
           ". Install these dependencies once; irtihmm itself need not be installed.")

    # Scope this host-specific linker workaround to compilation, restoring the
    # caller's environment afterwards. A user-specified Makevars takes precedence.
    previous_makevars <- Sys.getenv("R_MAKEVARS_USER", unset = NA_character_)
    local_makevars <- file.path(root, "tools", "Makevars.macos-no-fortran")
    if (Sys.info()[["sysname"]] == "Darwin" && is.na(previous_makevars) &&
        !dir.exists("/opt/gfortran/lib") && file.exists(local_makevars)) {
      Sys.setenv(R_MAKEVARS_USER = local_makevars)
      on.exit(Sys.unsetenv("R_MAKEVARS_USER"), add = TRUE)
    }

    cache <- file.path(root, "output", "sourceCpp")
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    for (file in c("hoskovec.cpp", "forward_backward.cpp", "sampler_kernels.cpp")) {
      Rcpp::sourceCpp(file.path(root, "src", file), env = target,
        cacheDir = cache, rebuild = FALSE, showOutput = FALSE, verbose = FALSE)
    }
    files <- sort(list.files(file.path(root, "R"), pattern = "\\.R$", full.names = TRUE))
    # Package registration wrappers would overwrite the sourceCpp wrappers with
    # references to an installed DLL. Package roxygen metadata is not executable.
    files <- files[!basename(files) %in% c("RcppExports.R", "package.R")]
    for (file in files) sys.source(file, envir = target)
    message("Loaded irtihmm sources (", length(files),
            " R files and 3 cached Rcpp modules); package installation is not required.")
    invisible(NULL)
  }
  load_sources(root, target)
})

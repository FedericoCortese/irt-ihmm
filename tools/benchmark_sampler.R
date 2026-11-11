# Short kernel benchmark, not an application fit. Run from the project root:
# Rscript tools/benchmark_sampler.R [optional_prechange_baseline.rds]
source("load.R")
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  baseline <- readRDS(args[1])
  sim <- baseline$sim; tr <- baseline$tr; hyper <- baseline$hyper
} else {
  set.seed(44)
  sim <- simulate_irtihmm(30,80,6,z=cbind(z1=rnorm(30),z2=rnorm(30)),seed=17)
  tr <- .new_transition(6,2,"infinite")
  hyper <- .initial_hyper(irtihmm_priors()); hyper$copula <- TRUE
}
data <- sim$data; C <- sim$states; emission <- sim$emission
elapsed <- function(fun,n=20L) {
  invisible(fun())
  unname(system.time(for(j in seq_len(n))fun())["elapsed"])/n
}
set.seed(29)
timing <- c(
  transition=elapsed(function().update_transition(tr,C,data$z)),
  probabilities=elapsed(function().transition_probs(tr,data$z),100),
  emission_array=elapsed(function().emission_array(data,emission),50),
  path_likelihood=elapsed(function().emission_path_loglik(data,emission,C),100),
  beam=elapsed(function().beam_state_draw(data,emission,tr,C,hyper,200)),
  emission_update=elapsed(function().update_emission(data,emission,hyper,C,irtihmm_priors(),TRUE)),
  short_fit=unname(system.time(fit_irtihmm(data,iterations=30,burnin=15,seed=11))["elapsed"])/30)
report <- data.frame(block=names(timing),seconds_per_call=as.numeric(timing))
if(length(args)) {
  report$previous_seconds <- unname(baseline$timing[report$block])
  report$speedup <- report$previous_seconds/report$seconds_per_call
}
dir.create("output/performance",recursive=TRUE,showWarnings=FALSE)
write.csv(report,"output/performance/sampler_benchmark.csv",row.names=FALSE)
print(report,row.names=FALSE)
cat("Synthetic 30 x 80, 6 generating states, 2 covariates, copula ON. Compilation excluded.\n",
    "Short-fit timings are per iteration, not a claim about ESS per second.\n")

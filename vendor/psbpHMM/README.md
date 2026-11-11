

This folder contains code for the R package psbpHMM. 

The main purpose of this package is to fit a covariate-dependent infinite hidden Markov model to multiple time series with or without missing data. The hidden state transition distribution is modeled via the probit stick-breaking process (psbp) which allows covariates to inform the transitions. 


This package relies on the following R packages. Install these packages in the R environment by using the install.packages("") command.  


Rcpp

RcppArmadillo

parallel

gdata

invgamma

gtools

mvtnorm

matrixcalc

tmvmixnorm

rlist

truncnorm

mvnfast


Then you can install psbpHMM by running the following lines in the R console: 

library(devtools)

install_github( "lvhoskovec/psbpHMM", build_vignettes = TRUE)

library(psbpHMM)

vignette("psbpHMM_vignette")


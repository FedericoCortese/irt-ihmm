// Forward recursion follows the slice indicator construction in psbpHMM/upZ.
// Emission-specific code is replaced by supplied log densities; an explicit
// backward draw produces a joint trajectory rather than successive filters.
// [[Rcpp::plugins(cpp17)]]
#include <Rcpp.h>
#include <cmath>
#include <limits>
using namespace Rcpp;

static double lse(const NumericVector& x) {
  double a = R_NegInf;
  for (double v : x) if (v > a) a = v;
  if (!R_finite(a)) return a;
  double s = 0;
  for (double v : x) s += std::exp(v - a);
  return a + std::log(s);
}

static int categorical_log(const NumericVector& weights) {
  double norm = lse(weights);
  if (!R_finite(norm)) stop("No feasible path in forward-backward recursion.");
  double u = R::runif(0, 1), total = 0;
  int last = 0;
  for (int k = 0; k < weights.size(); ++k) {
    double p = std::exp(weights[k] - norm);
    if (p > 0) last = k;
    total += p;
    if (u < total) return k;
  }
  return last;
}

// [[Rcpp::export]]
List forward_backward_cpp(NumericMatrix log_emission, NumericVector initial,
                          NumericMatrix transition, NumericVector slices,
                          bool draw = true) {
  int T = log_emission.nrow(), K = log_emission.ncol();
  if (T < 1 || K < 1 || initial.size() != K ||
      transition.nrow() != K || transition.ncol() != K ||
      (slices.size() != 0 && slices.size() != T)) stop("Invalid HMM dimensions.");
  bool beam = slices.size() > 0;
  NumericMatrix log_transition(K, K);
  if (!beam) for (int r=0;r<K;++r) for(int k=0;k<K;++k)
    log_transition(r,k)=std::log(transition(r,k));
  NumericMatrix forward(T, K);
  NumericVector row(K), work(K);
  for (int k = 0; k < K; ++k) {
    double value = beam ? (initial[k] > slices[0] ? 0 : R_NegInf) : std::log(initial[k]);
    row[k] = log_emission(0, k) + value;
  }
  double norm = lse(row), total_norm = norm;
  if (!R_finite(norm)) {
    if (draw) stop("Initial state has zero likelihood.");
    return List::create(_["path"] = IntegerVector(), _["loglik"] = norm);
  }
  for (int k = 0; k < K; ++k) forward(0, k) = row[k] - norm;
  for (int t = 1; t < T; ++t) {
    if (t % 100 == 0) checkUserInterrupt();
    for (int k = 0; k < K; ++k) {
      for (int r = 0; r < K; ++r) {
        double value = beam ? (transition(r, k) > slices[t] ? 0 : R_NegInf) : log_transition(r, k);
        work[r] = forward(t - 1, r) + value;
      }
      row[k] = log_emission(t, k) + lse(work);
    }
    norm = lse(row);
    if (!R_finite(norm)) {
      if (draw) stop("No feasible state trajectory.");
      return List::create(_["path"] = IntegerVector(), _["loglik"] = norm);
    }
    total_norm += norm;
    for (int k = 0; k < K; ++k) forward(t, k) = row[k] - norm;
  }
  IntegerVector path;
  if (draw) {
    path = IntegerVector(T);
    for (int k = 0; k < K; ++k) row[k] = forward(T - 1, k);
    path[T - 1] = categorical_log(row) + 1;
    for (int t = T - 2; t >= 0; --t) {
      int next = path[t + 1] - 1;
      for (int r = 0; r < K; ++r) {
        double value = beam ? (transition(r, next) > slices[t + 1] ? 0 : R_NegInf) : log_transition(r, next);
        row[r] = forward(t, r) + value;
      }
      path[t] = categorical_log(row) + 1;
    }
  }
  return List::create(_["path"] = path, _["loglik"] = total_norm);
}

// Batch trajectories in compiled code; RNG order matches the respondent loop.
// [[Rcpp::export]]
List hmm_paths_cpp(NumericVector emission, NumericMatrix initial, List transition,
                   NumericMatrix slices) {
  IntegerVector dims=emission.attr("dim");
  if(dims.size()!=3)stop("Expected a P by I by K emission array.");
  int P=dims[0],I=dims[1],K=dims[2];bool beam=slices.size()>0;
  if(initial.nrow()!=P||initial.ncol()!=K||transition.size()!=P ||
    (beam&&(slices.nrow()!=P||slices.ncol()!=I)))stop("Invalid batch HMM dimensions.");
  IntegerMatrix C(P,I);double loglik=0;
  for(int p=0;p<P;++p) {
    checkUserInterrupt();NumericMatrix e(I,K);NumericVector a(K),u(beam?I:0);
    for(int k=0;k<K;++k){a[k]=initial(p,k);for(int t=0;t<I;++t)e(t,k)=emission[p+P*t+P*I*k];}
    if(beam)for(int t=0;t<I;++t)u[t]=slices(p,t);
    List out=forward_backward_cpp(e,a,as<NumericMatrix>(transition[p]),u,true);
    IntegerVector path=out["path"];for(int t=0;t<I;++t)C(p,t)=path[t];
    loglik+=as<double>(out["loglik"]);
  }
  return List::create(_["C"]=C,_["loglik"]=loglik);
}

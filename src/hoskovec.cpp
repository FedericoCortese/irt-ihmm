// Derived from Lauren Hoskovec, psbpHMM, GPL (>= 2).
// Upstream commit: 1240b4da72a504dd5c97112f0b06f4b4a40f602a
// See UPSTREAM.md for exact provenance and modifications.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List hoskovec_pi(List beta, List X, arma::vec a0, arma::mat ajk, int tmax){
  
  // X is a List for each i 
  int n = X.size(); // number of people
  int K = beta.size();
  arma::vec first(K);
  arma::vec second(K);
  second[0] = 1; 
  arma::mat xi; 
  arma::mat x; 
  arma::mat xb; 
  arma::mat y; 
  NumericVector fourth; 
  double s; 
  arma::vec third;
  arma::vec pieces; 
  
  List probMat(n); // each element is a T-length list of matrices 
  
  for(int ind = 0; ind < n; ++ind){
    List pimat(tmax);
    xi = as<arma::mat>(X[ind]);
    
    for(int t = 0; t < tmax; ++t){
      x = xi.row(t);
      if(t == 0){
        for(int k = 0; k < K; ++k){
          xb = as<arma::rowvec>(beta[k])*x.t();
          y = a0[k] + xb; 
          first[k] = R::pnorm(as_scalar(y), 0, 1, 1, 0);
          if (k > 0) second[k] = 1-first[k-1];
        }
        second[0] = 1;
        third = cumprod(second);
        pieces = first % third;
        s = sum(pieces);
        fourth = wrap(pieces);
        fourth.push_back(1-s);
        pimat[t] = fourth;
      }else {
        arma::mat Z;
        for(int j = 0; j < K; ++j){
          for(int k = 0; k < K; ++k){
            xb = as<arma::rowvec>(beta[k])*x.t();
            double a = ajk(k,j);
            y = a + xb; 
            first[k] = R::pnorm(as_scalar(y), 0, 1, 1, 0);
            if (k > 0) second[k] = 1-first[k-1];
          }
          second[0] = 1;
          third = cumprod(second);
          pieces = first % third;
          s = sum(pieces);
          fourth = wrap(pieces);
          fourth.push_back(1-s);
          arma::vec four1 = as<arma::vec>(fourth);
          Z = join_rows(Z, four1);
        }
        pimat[t] = Z.t();
      }
    }
    probMat[ind] = pimat;
  }
  return(probMat);
}

// [[Rcpp::export]]
List hoskovec_alpha(int K, int n_sticks, int n, int tmax, List z, double vinv_alpha, double sig2inv_alpha, List w, List X, List beta_k, 
               double m_alpha, double mu_alpha){
  
  // loop thru j and k
  int i;
  List wi; 
  arma::mat xi; 
  arma::mat xit;
  arma::mat wk; 
  arma::mat wpiece; 
  double wdouble; 
  double wsum = 0.0; 
  double sigjk;
  NumericVector ivec;
  double sig2jk;
  double mujk;
  arma::vec isave;
  isave.set_size(0);
  double ans;
  int njk;
  List aj(K);
  for(int j = 0; j < K; ++j){
    NumericVector ajk(n_sticks);
    for(int k = 0; k < n_sticks; ++k){
      njk = 0;
      wsum = 0.0;
      for(i = 0; i < n; ++i){
        ivec = wrap(isave);
        arma::vec zi = as<arma::vec>(z[i]);
        wi = w[i];
        xi = as<arma::mat>(X[i]);
        
        for(int t = 1; t < tmax; ++t){
          if(zi[t]>=(k+1) && zi[t-1]==(j+1)){
            
            // add the w part 
            xit = xi.row(t);
            arma::mat xbl = as<arma::rowvec>(beta_k[k])*xit.t();
            arma::mat wpart = as<arma::mat>(wi[t]);
            wk = wpart[k]; // indexing issue  
            wpiece = wk - xbl; 
            wdouble = as_scalar(wpiece);
            wsum = wsum + wdouble;
            
            ivec.push_back(t);
          }
        }
        njk = njk + ivec.length();
      }
      
      if(j!=k){
        sig2jk = 1/(sig2inv_alpha + njk);
        mujk = sig2jk*(mu_alpha*sig2inv_alpha + wsum);
      }else{
        sig2jk = 1/(vinv_alpha + njk); 
        mujk = sig2jk*(m_alpha*vinv_alpha + wsum);
      }
      sigjk = pow(sig2jk, 0.5);
      ans = R::rnorm(mujk, sigjk);
      ajk[k] = ans;
    }
    aj[j] = ajk;
  }
  return aj;
}

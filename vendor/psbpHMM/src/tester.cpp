#include <RcppArmadillo.h>
#include <Rcpp.h>
// [[Rcpp::depends(RcppArmadillo]]
#include <RcppArmadilloExtensions/sample.h>
using namespace Rcpp;


namespace help{

  // which equal 
  NumericVector getWhich_equal(NumericVector a, int k){
    
    sugar::SeqLen seqK = seq_len(a.length()); // 1 to n
    NumericVector seq_vec;
    seq_vec = seqK; // remove 0 so its 1 to n
    arma::vec seq1 = as<arma::vec>(seq_vec); 
    
    NumericVector ans;
    ans = (a == k); // 1s and 0s
    arma::vec x = as<arma::vec>(ans);
    
    arma::vec y = x%seq1; // 0s and indices
    NumericVector z = wrap(y);
    z = sort_unique(z);
    if(z[0]==0) z.erase(0);
    return z; 
  }



  // log density 
  double logmvndensity(arma::vec y, arma::vec mu, arma::mat Sigma, double d){
    double detSig = log_det_sympd(Sigma); 
    double norm_part = -0.5*detSig + (-d/2)*log(2*M_PI);
    double exp_part = as_scalar(-0.5 * (y - mu).t() * inv(Sigma) * (y - mu));
    double ans = norm_part + exp_part;
    return(ans);
  }

  // mean of are these two vectors the same? 
  double vectorSame(NumericVector a, NumericVector b){
    int n = a.length();
    LogicalVector c(n); 
    for(int i=0; i<n; ++i){
      c[i] = (a[i]==b[i]);
    }
    double d = mean(c); 
    return(d); 
  }


  // log sum trick 
  arma::vec logTrick(arma::vec probs){
    
    // log trick 
    double maxlog = max(probs);
    double thresh = log(.00000000000000000001) - log(probs.size());
    double sum_denom = 0; 
    for(int s = 0; s < probs.size(); ++s){
      if(probs[s] - maxlog > thresh){
        sum_denom += exp(probs[s]-maxlog);
        probs[s] = exp(probs[s]-maxlog);
      }else{
        probs[s] = 0; 
      }
    }
    // normalize probability vector
    arma::vec probVec = probs/sum_denom; // real normalized probabilities 
    if(sum(probVec)<0.9999999999999){
      throw std::range_error("sum probs neq 1");
    }
    
    return probVec;
    
  }

  //multivariate normal density 
  double mvndensity(arma::vec y, arma::vec mu, arma::mat Sigma, double d){
    arma::vec x1; 
    double x2 = det(Sigma); 
    x1 = x2; 
    arma::vec x3 = pow( x1, -0.5 );
    arma::vec x4 = pow( (2*M_PI), (-d/2))*x3*exp(-0.5 * (y - mu).t() * inv(Sigma) * (y - mu));
    double x5 = as_scalar(x4);
    return(x5);
  }


  // get index of vector elements of a greater than or equal to num
  NumericVector getWhich_greater(NumericVector a, double num, int K){
    sugar::SeqLen seqK = seq_len(K+1);
    NumericVector seqxK;
    seqxK = seqK-1;
    
    NumericVector ans;
    ans = a >= num;
    arma::vec x = as<arma::vec>(ans);
    int n = a.length();
    sugar::SeqLen seq = seq_len(n);
    NumericVector seqx;
    seqx = seq;
    arma::vec seqa = as<arma::vec>(seqx);
    arma::vec mult = x%seqa;
    NumericVector ret = wrap(mult);
    
    ret = sort_unique(ret - 1);
    ret = intersect(ret, seqxK);
    ret = sort_unique(ret);
    
    return ret;
  }


  // call sample function from R 
  double sample_cpp( IntegerVector x, NumericVector p ){
    Function f("sample");
    IntegerVector sampled = f(x, 1, 0, p);
    return sampled[0];
  }

  // call the which function from R
  // call anyNA from R 
  bool anyNA_cpp(NumericVector x){
    Function f("anyNA");
    bool y = f(x);
    return y;
  }

  double rtnorm1(double a) {
    
    double x = 0.0;
    double y;
    
    if(a<0){
      // normal rejection sampling
      while(x==0){
        y = R::rnorm(0,1);
        if(y>a){
          x=y;
        }
      }
    }else if(a<0.2570){
      // half-normal rejection sampling
      while(x==0){
        y = fabs(R::rnorm(0,1));
        if(y>a){
          x=y;
        }
      }
    }else{
      // one-sided translated-exponential rejection sampling
      while(x==0){
        double lambdastar = (a + sqrt(a*a+4.0))/2.0;
        y = R::rexp(1)/lambdastar + a;
        if(R::runif(0,1) < exp(-0.5*y*y + lambdastar*y - 0.5*lambdastar + log(lambdastar))){
          x=y;
        }
      }
    }
    
    return x;
  }

}




// [[Rcpp::export]]
List up_ajk_rm(int K, int n, int tmax, List z, double vinv_alpha, double sig2inv_alpha, List w, List X, List beta_k, 
           List beta_sk, double m_alpha, double mu_alpha){
  
  // loop thru j and k
  int i;
  List wi; 
  List beta_ski;
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
    NumericVector ajk(K);
    for(int k = 0; k < K; ++k){
      njk = 0;
      wsum = 0.0;
      for(i = 0; i < n; ++i){
        ivec = wrap(isave);
        arma::vec zi = as<arma::vec>(z[i]);
        wi = w[i];
        beta_ski = beta_sk[i];
        xi = as<arma::mat>(X[i]);

        for(int t = 1; t < tmax; ++t){
          if(zi[t]>=(k+1) & zi[t-1]==(j+1)){
            
            // add the w part 
            xit = xi.row(t);
            arma::mat xbl = as<arma::rowvec>(beta_k[k])*xit.t();
            arma::mat xbli_rm = as<arma::rowvec>(beta_ski[k])*xit.t();
            arma::mat wpart = as<arma::mat>(wi[t]);
            wk = wpart[k]; // indexing issue  
            wpiece = wk - xbl - xbli_rm; 
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

// [[Rcpp::export]]
List up_ajk(int K, int n, int tmax, List z, double vinv_alpha, double sig2inv_alpha, List w, List X, List beta_k, 
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
    NumericVector ajk(K);
    for(int k = 0; k < K; ++k){
      njk = 0;
      wsum = 0.0;
      for(i = 0; i < n; ++i){
        ivec = wrap(isave);
        arma::vec zi = as<arma::vec>(z[i]);
        wi = w[i];
        xi = as<arma::mat>(X[i]);
        
        for(int t = 1; t < tmax; ++t){
          if(zi[t]>=(k+1) & zi[t-1]==(j+1)){
            
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

// [[Rcpp::export]]
List up_ajk_nox(int K, int n, int tmax, List z, double vinv_alpha, double sig2inv_alpha, List w, double m_alpha, double mu_alpha){
  
  // loop thru j and k
  int i;
  List wi; 

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
    NumericVector ajk(K);
    for(int k = 0; k < K; ++k){
      njk = 0;
      wsum = 0.0;
      for(i = 0; i < n; ++i){
        ivec = wrap(isave);
        arma::vec zi = as<arma::vec>(z[i]);
        wi = w[i];
        for(int t = 1; t < tmax; ++t){
          if(zi[t]>=(k+1) & zi[t-1]==(j+1)){
            arma::mat wpart = as<arma::mat>(wi[t]);
            wk = wpart[k]; // indexing issue  
            wpiece = wk; 
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



// [[Rcpp::export]]
List upW_rm(arma::vec alpha0, List X, List beta, List beta_rm, arma::mat ajk, List z, int tmax, int n){
  
  int l; 
  double m1; 
  double out;
  arma::vec w; 
  int zitminus1 = 0; 
  double apart; 
  List wall(n); 
  arma::vec zi; 
  List beta_i; 
  arma::mat xi; 
  int zit; 
  arma::mat xit; 
  arma::mat m; 
  
  // loop through i 
  for(int i = 0; i < n; ++i){
    zi = as<arma::vec>(z[i]);
    beta_i = beta_rm[i]; 
    xi = as<arma::mat>(X[i]);
    
    // loop through t
    List wi(tmax); 
    for(int t = 0; t < tmax; ++t){
      zit = zi[t]; // cpp indexing 
      xit = xi.row(t); 
      w.set_size(zit);
      
      if(t > 0) {
        zitminus1 = zi[t-1]-1; 
      }
      
      // loop through l 
      for(l = 0; l < zit; ++l){
        //begin loop
        arma::mat xbl = as<arma::rowvec>(beta[l])*xit.t();
        arma::mat xbli_rm = as<arma::rowvec>(beta_i[l])*xit.t();
        if(t == 0) {
          apart = alpha0[l];
        }else{
          apart = ajk(zitminus1, l);
        }
        m = apart + xbl + xbli_rm; 
        m1 = as_scalar(m); 
        if(l<(zit-1)){
          // negative 
          out = -1*(help::rtnorm1(m1)) + m1;
          //out = -1;
        }else{
          // positive 
          out = help::rtnorm1(-m1) + m1;  // check this add m1 
          //out = 1; 
        } 
        w[l] = out; 
        // end loop 
      }
      wi[t] = w;  
    }
    wall[i] = wi; 
  }
  return(wall); 
}


// [[Rcpp::export]]
List upW(arma::vec alpha0, List X, List beta, arma::mat ajk, List z, int tmax, int n){
  
  int l; 
  double m1; 
  double out;
  arma::vec w; 
  int zitminus1 = 0; 
  double apart; 
  List wall(n); 
  arma::vec zi; 
  arma::mat xi; 
  int zit; 
  arma::mat xit; 
  arma::mat m; 
  
  // loop through i 
  for(int i = 0; i < n; ++i){
    zi = as<arma::vec>(z[i]);
    xi = as<arma::mat>(X[i]);
    
    // loop through t
    List wi(tmax); 
    for(int t = 0; t < tmax; ++t){
      zit = zi[t]; // cpp indexing 
      xit = xi.row(t); 
      w.set_size(zit);
      
      if(t > 0) {
        zitminus1 = zi[t-1]-1; 
      }
      
      // loop through l 
      for(l = 0; l < zit; ++l){
        //begin loop
        arma::mat xbl = as<arma::rowvec>(beta[l])*xit.t();
        if(t == 0) {
          apart = alpha0[l];
        }else{
          apart = ajk(zitminus1, l);
        }
        m = apart + xbl;
        m1 = as_scalar(m); 
        if(l<(zit-1)){
          // negative 
          out = -1*(help::rtnorm1(m1)) + m1;
          //out = -1;
        }else{
          // positive 
          out = help::rtnorm1(-m1) + m1;  // check this add m1 
          //out = 1; 
        } 
        w[l] = out; 
        // end loop 
      }
      wi[t] = w;  
    }
    wall[i] = wi; 
  }
  return(wall); 
}

// [[Rcpp::export]]
List upW_nox(arma::vec alpha0, arma::mat ajk, List z, int tmax, int n){
  
  int l; 
  double out;
  arma::vec w; 
  int zitminus1 = 0; 
  double apart; 
  List wall(n); 
  arma::vec zi; 
  int zit; 
  double m; 
  
  // loop through i 
  for(int i = 0; i < n; ++i){
    zi = as<arma::vec>(z[i]);
    
    // loop through t
    List wi(tmax); 
    for(int t = 0; t < tmax; ++t){
      zit = zi[t]; // cpp indexing 
      w.set_size(zit);
      
      if(t > 0) {
        zitminus1 = zi[t-1]-1; 
      }
      
      // loop through l 
      for(l = 0; l < zit; ++l){
        //begin loop
        if(t == 0) {
          apart = alpha0[l];
        }else{
          apart = ajk(zitminus1, l);
        }
        m = apart;
        if(l<(zit-1)){
          // negative 
          out = -1*(help::rtnorm1(m)) + m;
          //out = -1;
        }else{
          // positive 
          out = help::rtnorm1(-m) + m;  
          //out = 1; 
        } 
        w[l] = out; 
        // end loop 
      }
      wi[t] = w;  
    }
    wall[i] = wi; 
  }
  return(wall); 
}


// [[Rcpp::export]]
arma::mat invMat(arma::mat x) {
  return( inv(x));
}

// [[Rcpp::export]]
arma::mat mhDecomp(arma::mat L, arma::mat D) {
  return( inv(L) * D * inv(L).t() );
}

// [[Rcpp::export]]
List updatePi(List beta, List X, arma::vec a0, arma::mat ajk, int tmax){
  
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
List updatePi_rm(List beta, List beta_sk, List X, arma::vec a0, arma::mat ajk, int tmax){
  
  // X is a List for each i 
  int n = X.size(); // number of people
  int K = beta.size();
  arma::vec first(K);
  arma::vec second(K);
  second[0] = 1; 
  arma::mat xi; 
  arma::mat x; 
  arma::mat xb; 
  arma::mat xbs; 
  List bs; 
  arma::mat y; 
  NumericVector fourth; 
  double s; 
  arma::vec third;
  arma::vec pieces; 
  
  List probMat(n); // each element is a T-length list of matrices 
  
  for(int ind = 0; ind < n; ++ind){
    List pimat(tmax);
    xi = as<arma::mat>(X[ind]);
    bs = beta_sk[ind]; 
    
    for(int t = 0; t < tmax; ++t){
      x = xi.row(t);
      if(t == 0){
        for(int k = 0; k < K; ++k){
          xb = as<arma::rowvec>(beta[k])*x.t();
          xbs = as<arma::rowvec>(bs[k])*x.t();
          y = a0[k] + xb + xbs; 
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
            xbs = as<arma::rowvec>(bs[k])*x.t();
            double a = ajk(k,j);
            y = a + xb + xbs; 
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
double mvndensity(arma::vec y, arma::vec mu, arma::mat Sigma, double d){
  arma::vec x1; 
  double x2 = det(Sigma); 
  x1 = x2; 
  arma::vec x3 = pow( x1, -0.5 );
  arma::vec x4 = pow( (2*M_PI), (-d/2))*x3*exp(-0.5 * (y - mu).t() * inv(Sigma) * (y - mu));
  double x5 = as_scalar(x4);
  return(x5);
}


// [[Rcpp::export]]
double logmvndensity(arma::vec y, arma::vec mu, arma::mat Sigma, double d){
  double detSig = log_det_sympd(Sigma); 
  double norm_part = -0.5*detSig + (-d/2)*log(2*M_PI);
  double exp_part = as_scalar(-0.5 * (y - mu).t() * inv(Sigma) * (y - mu));
  double ans = norm_part + exp_part;
  return(ans);
}

// [[Rcpp::export]]
List upZ(List stateList, List y, List mu, List Sigma, double logStuff, double nudf, List detRstar, List piz,
         List u, int tmax, int K, int n, double d){
  
  // initialization 
  List z(n); // set to n 
  List sli; // state list at i 
  arma::mat yi; // y at i 
  arma::vec yit; // y_i at t 
  IntegerVector kprime; 
  int nprime;  
  arma::vec logTKprobs; 
  arma::vec tkProbs;
  //double maxlog; 
  //double logsum; 
  IntegerVector jprime; 
  arma::vec jprobs; 
  NumericVector probVec; 
  IntegerVector ret; 
  int state; 
  arma::vec zt;
  int idx; 
  int k; 
  arma::mat sigK;
  arma::vec muK;

  
  arma::vec detRi; 
  double logR; 
  double newState; 
  double loglik; 
  int s1; 
  int jindex; 
  int s2; 
  double probj;
  List pizi; 
  arma::mat pizit; 
  arma::vec pizitk; 
  arma::vec ui; 
  double uit; 
  IntegerVector slit;
  IntegerVector prevK; 
  IntegerVector whichK(0);
  arma::vec whichKtemp;

  double logsumall; 
  int kstar; 
  int t; 
  
  double sumAllProbs; 
  
  for(int i = 0; i < n; ++i){
    sli = stateList[i];
    yi = as<arma::mat>(y[i]);
    zt.set_size(tmax); 
    pizi = piz[i]; 
    ui = as<arma::vec>(u[i]);
    
    // first time point 
    t = 0; 
    yit = (yi.row(t)).t(); // redefined for each t 
    kprime = sli[t];
    nprime = kprime.length(); 
    
    logTKprobs.set_size(nprime); 
    for(idx = 0; idx < nprime; ++idx){
      k = kprime[idx]; 
      if(k <= K){
        sigK = as<arma::mat>(Sigma[k-1]);
        muK = as<arma::mat>(mu[k-1]).t();
        loglik = help::logmvndensity(yit, muK, sigK, d);
        logTKprobs[idx] = loglik;
      }else{
        detRi = as<arma::vec>(detRstar[i]);
        logR = log(detRi[t]); 
        newState = logStuff - ((nudf + 1)/2)*logR;
        logTKprobs[idx] = newState;
      }
    }
    tkProbs = help::logTrick(logTKprobs); // arma::vec tkProbs are true normalized probabilities 
    if(nprime == 1){
      zt[t] = kprime[0]; // only one element in kprime
    }else{
      probVec = wrap(tkProbs);
      ret = RcppArmadillo::sample(kprime, 1, 0, probVec);
      state = ret[0]; 
      zt[t] = state;
    }
    // the rest of the time points 
    for(t = 1; t < tmax; ++t){
      pizit = as<arma::mat>(pizi[t]); 
      uit = ui[t]; 
      yit = (yi.row(t)).t(); // y_it 
      jprime = kprime; // previous possible states
      jprobs = tkProbs; // previous state probabilities 
      kprime = sli[t]; // update kprime, current possible states
      nprime = kprime.length(); 
      logTKprobs.set_size(nprime); // reset 
      slit = sli[t-1];
      
      // log likelihood for all current possible states kprime 
      for(idx = 0; idx < nprime; ++idx){
        // reset whichK
        whichKtemp = as<arma::vec>(whichK); 
        whichKtemp.set_size(0);
        whichK = wrap(whichKtemp); 
        // for each k we sum over previous possible states for logsumprobs 
        k = kprime[idx]; 
        // because the indexing starts at 0 for mu and Sigma 
        if(k <= K){
          sigK = as<arma::mat>(Sigma[k-1]);
          muK = as<arma::mat>(mu[k-1]).t();
          loglik = help::logmvndensity(yit, muK, sigK, d); // make this the log likelihood
        }else{
          detRi = as<arma::vec>(detRstar[i]);
          logR = log(detRi[t]); 
          loglik = logStuff - ((nudf + 1)/2)*logR; // loglik
        }
        pizitk = pizit.col(k-1); 
        for(kstar = 0; kstar < K; ++kstar){
          if(pizitk[kstar] >= uit) {
            whichK.push_back(kstar+1); // indexing 
          }
        }
        prevK = intersect(whichK, slit); // previous possible states for the current possible state k 
        sumAllProbs = 0.0; 
        
        for(s1 = 0; s1 < jprime.length(); s1++){
          jindex = jprime[s1];
          for(s2 = 0; s2 < prevK.length(); s2++){
            if(prevK[s2]==jindex){
              probj = jprobs[s1];
              sumAllProbs += probj;
              break;
            }
          }
        }
        
        logsumall = log(sumAllProbs); 
        logTKprobs[idx] = loglik + logsumall;
      }

      tkProbs = help::logTrick(logTKprobs); // arma::vec tkProbs are true normalized probabilities 

      if(nprime == 1){
        zt[t] = kprime[0]; // only one element in kprime
      }else{
        probVec = wrap(tkProbs);
        ret = RcppArmadillo::sample(kprime, 1, 0, probVec);
        state = ret[0]; 
        zt[t] = state;
      }
    }
    z[i] = zt; 
  }
  //return(logsumall);
  
  return(z); 
}

// [[Rcpp::export]]
List upZnox(List stateList, List y, List mu, List Sigma, double logStuff, double nudf, List detRstar, List piz, 
            List u, int tmax, int K, int n, double d){
  
  // no X so piz is the same for all time points 
  
  // initialization 
  List z(n); // set to n 
  List sli; // state list at i 
  arma::mat yi; // y at i 
  arma::vec yit; // y_i at t 
  IntegerVector kprime; 
  int nprime;  
  arma::vec logTKprobs; 
  arma::vec tkProbs;

  IntegerVector jprime; 
  arma::vec jprobs; 
  NumericVector probVec; 
  IntegerVector ret; 
  int state; 
  arma::vec zt;
  int idx; 
  int k; 
  arma::mat sigK;
  arma::vec muK;

  arma::vec detRi; 
  double logR; 
  double newState; 
  double loglik; 
  int s1; 
  int jindex; 
  int s2; 
  double probj;
  //List pizi; 
  arma::mat pizit; 
  arma::vec pizitk; 
  arma::vec ui; 
  double uit; 
  IntegerVector slit;
  IntegerVector prevK; 
  IntegerVector whichK(0);
  arma::vec whichKtemp;
  double logsumall; 
  int idx2; 
  int kstar; 
  
  
  double sumAllProbs; 
  
  for(int i = 0; i < n; ++i){
    sli = stateList[i];
    yi = as<arma::mat>(y[i]);
    zt.set_size(tmax); 
    ui = as<arma::vec>(u[i]);
    
    // first time point 
    int t = 0; 
    yit = (yi.row(t)).t(); // redefined for each t 
    kprime = sli[t];
    nprime = kprime.length(); 
    logTKprobs.set_size(nprime); 
    for(idx = 0; idx < nprime; ++idx){
      k = kprime[idx]; // because the indexing starts at 0 for mu and Sigma 
      if(k <= K){
        sigK = as<arma::mat>(Sigma[k-1]);
        muK = as<arma::mat>(mu[k-1]).t();
        loglik = help::logmvndensity(yit, muK, sigK, d);
        logTKprobs[idx] = loglik;
      }else{
        detRi = as<arma::vec>(detRstar[i]);
        logR = log(detRi[t]); 
        newState = logStuff - ((nudf + 1)/2)*logR;
        logTKprobs[idx] = newState;
      }
    }
    tkProbs = help::logTrick(logTKprobs); // arma::vec tkProbs are true normalized probabilities 
    if(nprime == 1){
      zt[t] = kprime[0]; // only one element in kprime
    }else{
      probVec = wrap(tkProbs);
      ret = RcppArmadillo::sample(kprime, 1, 0, probVec);
      state = ret[0]; 
      zt[t] = state;
    }
    
    // the rest of the time points 
    for(t = 1; t < tmax; ++t){
      pizit = as<arma::mat>(piz[1]); // second index 
      uit = ui[t]; 
      yit = (yi.row(t)).t(); // y_it 
      jprime = kprime; // previous possible states
      jprobs = tkProbs; // previous state probabilities 
      kprime = sli[t]; // update kprime, current possible states
      nprime = kprime.length(); 
      logTKprobs.set_size(nprime); // reset 
      slit = sli[t-1];
      
      // log likelihood for all current possible states kprime 
      for(idx = 0; idx < nprime; ++idx){
        // reset whichK
        whichKtemp = as<arma::vec>(whichK); 
        whichKtemp.set_size(0);
        whichK = wrap(whichKtemp); 
        // for each k we sum over previous possible states for logsumprobs 
        k = kprime[idx]; 
        // because the indexing starts at 0 for mu and Sigma 
        if(k <= K){
          sigK = as<arma::mat>(Sigma[k-1]);
          muK = as<arma::mat>(mu[k-1]).t();
          loglik = help::logmvndensity(yit, muK, sigK, d);
          logTKprobs[idx] = loglik;
        }else{
          detRi = as<arma::vec>(detRstar[i]);
          logR = log(detRi[t]); 
          loglik = logStuff - ((nudf + 1)/2)*logR; // loglik
        }
        pizitk = pizit.col(k-1); 
        idx2 = 0; 
        for(kstar = 0; kstar < K; ++kstar){
          if(pizitk[kstar] >= uit) {
            whichK.push_back(kstar+1); // indexing 
            idx2++; 
          }
        }
        prevK = intersect(whichK, slit);
        sumAllProbs = 0.0; 
        
        for(s1 = 0; s1 < jprime.length(); s1++){
          jindex = jprime[s1];
          for(s2 = 0; s2 < prevK.length(); s2++){
            if(prevK[s2]==jindex){
              probj = jprobs[s1];
              sumAllProbs += probj;
              break;
            }
          }
        }
        
        logsumall = log(sumAllProbs); 
        logTKprobs[idx] = loglik + logsumall;
      }
      tkProbs = help::logTrick(logTKprobs); // arma::vec tkProbs are true normalized probabilities 
      if(nprime == 1){
        zt[t] = kprime[0]; // only one element in kprime
      }else{
        probVec = wrap(tkProbs);
        ret = RcppArmadillo::sample(kprime, 1, 0, probVec);
        state = ret[0]; 
        zt[t] = state;
      }
    }
    z[i] = zt; 
  }
  return(z); 
}

// [[Rcpp::export]]
NumericVector getWhich_greater(NumericVector a, double num, int K){
  sugar::SeqLen seqK = seq_len(K+1); // 1 to K+1
  NumericVector seqxK;
  seqxK = seqK-1; // 0 to K 
  
  NumericVector ans;
  ans = a >= num; // 1's and 0's
  arma::vec x = as<arma::vec>(ans);
  int n = a.length();
  sugar::SeqLen seq = seq_len(n); // 1 to length of a 
  NumericVector seqx;
  seqx = seq;
  arma::vec seqa = as<arma::vec>(seqx);
  arma::vec mult = x%seqa; // get the indices
  NumericVector ret = wrap(mult);
  ret = sort_unique(ret - 1); // the indices are off by one 
  ret = intersect(ret, seqxK); // remove the -1 if it is there 
  ret = sort_unique(ret);
  
  return ret;
}


// [[Rcpp::export]]
List upStateList(List piz, List u, int K, int tmax, int n){
  
  List stateListAll(n);
  
  List pizi;
  arma::vec ui;
  arma::mat pizit; 
  double uit;
  int k;
  int t; 
  int idx; 
  int tempK; 
  NumericVector pizitk;
  NumericVector pizitk0; 
  
  double uitplus1; 
  arma::mat pizitplus1; 
  
  NumericVector whichK;
  NumericVector setzero(0);
  // loop through individuals 
  for(int i = 0; i < n; ++i){
    
    List stateList(tmax); 
    List forList(tmax); 
    List backList(tmax);
    
    ui = as<arma::vec>(u[i]);
    pizi = piz[i]; // list
    
    // first time point 
    t = 0; 
    pizitk0 = pizi[t]; 
    uit = ui[t]; 
    whichK = setzero;
    for(k = 0; k < K+1; ++k){
      if(pizitk0[k] >= uit){
        whichK.push_back(k); 
      }
    }
    
    forList[t] = whichK; 
    // rest of time points 
    for(t = 1; t < tmax; ++t){
      NumericVector whichKprev = forList[t-1];
      uit = ui[t]; 
      whichK = setzero; 
      pizit = as<arma::mat>(pizi[t]); // prob matrix for current time point 
      for(idx = 0; idx < whichKprev.length(); ++idx){
        tempK = whichKprev[idx]; // one previous state at a time 
        if(tempK < K){ // can't go from a new state 
          pizitk = (pizit.row(tempK)).t(); // get the transition probs for the one previous state 
          for(k = 0; k < K+1; ++k){
            if(pizitk[k] >= uit){
              whichK.push_back(k); 
            }
          }
        }
      }
      forList[t] = sort_unique(whichK); 
    }
    
    // last time point 
    backList[tmax-1] = forList[tmax-1]; 
    for(t = tmax-2; t>=0; --t){
      NumericVector whichKnext = backList[t+1];
      uitplus1 = ui[t+1]; 
      whichK = setzero; 
      pizitplus1 = as<arma::mat>(pizi[t+1]); // prob matrix for current time point 
      for(idx = 0; idx < whichKnext.length(); ++idx){
        tempK = whichKnext[idx]; // one previous state at a time 
        pizitk = (pizitplus1.col(tempK)); // get the transition probs for the one previous state 
        // whichK = help::getWhich_greater(pizitk, uitplus1, K);
        for(k = 0; k < K; ++k){
          if(pizitk[k] >= uitplus1){
            whichK.push_back(k); 
          }
        }
      }
      backList[t] = sort_unique(whichK); 
    }
    
    // intersect
    for(t = 0; t < tmax; ++t){
      IntegerVector f3 = forList[t];
      IntegerVector f4 = backList[t];
      stateList[t] = sort_unique(intersect(f3, f4))+1;
    }
    
    stateListAll[i] = stateList;
    
  }
  
  return(stateListAll);
}


// [[Rcpp::export]]
List upStateListnox(List piz, List u, int K, int tmax, int n){
  
  List stateListAll(n);
  
  //List pizi;
  arma::vec ui;
  arma::mat pizit; 
  double uit;
  int k;
  int t; 
  arma::vec whichKtemp; 
  IntegerVector whichK(0);
  int idx; 
  arma::vec whichKminus1; 
  int tempK; 
  arma::vec pizitk;
  IntegerVector whichKprev; 
  arma::vec pizitk0; 
  double uitplus1; 
  arma::vec whichKplus1; 
  IntegerVector whichKnext; 
  arma::mat pizitplus1; 
  arma::vec f1; 
  arma::vec f2; 
  IntegerVector f3; 
  IntegerVector f4; 
  
  
  // loop through individuals 
  for(int i = 0; i < n; ++i){
    
    List stateList(tmax); 
    List forList(tmax); 
    List backList(tmax);
    ui = as<arma::vec>(u[i]);
    // first time point 
    t = 0; 
    pizitk0 = as<arma::vec>(piz[0]); 
    uit = ui[t]; 
    whichKtemp = as<arma::vec>(whichK); 
    whichKtemp.set_size(0);
    whichK = wrap(whichKtemp); 
    for(k = 0; k < K+1; ++k){
      if(pizitk0[k] >= uit){
        whichK.push_back(k); 
      }
    }
    forList[t] = whichK; 
    // rest of time points 
    
    pizit = as<arma::mat>(piz[1]); // prob matrix for current time point 
    for(t = 1; t < tmax; ++t){
      whichKminus1 = as<arma::vec>(forList[t-1]); // vector of previous states 
      whichKprev = wrap(whichKminus1); // IntegerVector
      uit = ui[t]; 
      whichKtemp = as<arma::vec>(whichK); 
      whichKtemp.set_size(0);
      whichK = wrap(whichKtemp); 
      
      for(idx = 0; idx < whichKprev.length(); ++idx){
        tempK = whichKprev[idx]; // one previous state at a time 
        if(tempK < K){ // can't go from a new state 
          pizitk = (pizit.row(tempK)).t(); // get the transition probs for the one previous state 
          for(k = 0; k < K+1; ++k){
            if(pizitk[k] >= uit){
              whichK.push_back(k); 
            }
          }
        }
      }
      forList[t] = sort_unique(whichK); 
    }
    
    // last time point 
    backList[tmax-1] = forList[tmax-1]; 
    for(t = tmax-2; t>=0; --t){
      whichKplus1 = as<arma::vec>(backList[t+1]); // vector of previous states 
      whichKnext = wrap(whichKplus1); // IntegerVector
      uitplus1 = ui[t+1]; 
      
      whichKtemp = as<arma::vec>(whichK); 
      whichKtemp.set_size(0);
      whichK = wrap(whichKtemp); 
      
      pizitplus1 = as<arma::mat>(piz[1]); // prob matrix for current time point 
      for(idx = 0; idx < whichKnext.length(); ++idx){
        tempK = whichKnext[idx]; // one previous state at a time 
        pizitk = (pizitplus1.col(tempK)); // get the transition probs for the one previous state 
        for(k = 0; k < K; ++k){
          if(pizitk[k] >= uitplus1){
            whichK.push_back(k); 
          }
        }
      }
      backList[t] = sort_unique(whichK); 
    }
    
    // intersect
    for(t = 0; t < tmax; ++t){
      f1 = as<arma::vec>(forList[t]); 
      f2 = as<arma::vec>(backList[t]); 
      f3 = wrap(f1);
      f4 = wrap(f2); 
      // needs to be a vector 
      stateList[t] = sort_unique(intersect(f3, f4))+1;
    }
    
    stateListAll[i] = stateList;
    
  }
  
  return(stateListAll);
}


// [[Rcpp::export]]
NumericVector getWhich_equal(NumericVector a, int k){
  
  sugar::SeqLen seqK = seq_len(a.length()); // 1 to n
  NumericVector seq_vec;
  seq_vec = seqK; // remove 0 so its 1 to n
  arma::vec seq1 = as<arma::vec>(seq_vec); 
  
  NumericVector ans;
  ans = (a == k); // 1s and 0s
  arma::vec x = as<arma::vec>(ans);
  
  arma::vec y = x%seq1; // 0s and indices
  NumericVector z = wrap(y);
  z = sort_unique(z);
  if(z[0]==0) z.erase(0);

  
  return z; 
}


// [[Rcpp::export]]
double vectorSame(NumericVector a, NumericVector b){
  int n = a.length();
  LogicalVector c(n); 
  for(int i=0; i<n; ++i){
    c[i] = (a[i]==b[i]);
  }
  double d = mean(c); 
  return(d); 
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/

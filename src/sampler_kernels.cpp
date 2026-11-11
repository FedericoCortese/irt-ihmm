// Compiled numerical kernels for the existing PSBP Gibbs/beam and emission updates.
// The model, Gaussian conditionals and slice targets are unchanged.
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(cpp17)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <unordered_map>
#include <climits>
using namespace Rcpp;

namespace {
double softplus(double x) { return std::max(x, 0.0) + std::log1p(std::exp(-std::abs(x))); }
double threshold(double a) {
  return a >= 0 ? R::qnorm(-softplus(a), 0, 1, true, true) :
    R::qnorm(-softplus(-a), 0, 1, false, true);
}
double log_density(double x, double y, double a, double m, double sigma, double rho) {
  if (!R_finite(a) || !R_finite(m) || !R_finite(sigma) || sigma <= 0 ||
      !R_finite(rho) || std::abs(rho) >= 1) return R_NegInf;
  double z = (y-m)/sigma;
  double ans = R::dnorm(z, 0, 1, true) - std::log(sigma);
  if (rho == 0) return ans - softplus(x == 1 ? -a : a);
  double cut = (threshold(a)-rho*z)/std::sqrt(1-rho*rho);
  return ans + R::pnorm(x == 1 ? -cut : cut, 0, 1, true, true);
}
int index(double value, int size) {
  if (!R_finite(value) || value < 1 || value > size || value != std::floor(value))
    stop("Invalid parameter index in numerical kernel.");
  return static_cast<int>(value)-1;
}
arma::vec precision_draw(const arma::mat& Q, const arma::vec& h) {
  arma::mat U;
  if (!arma::chol(U,Q)) stop("Gaussian precision is not positive definite.");
  arma::vec noise(Q.n_rows);
  for (unsigned j=0; j<noise.n_elem; ++j) noise[j]=R::rnorm(0,1);
  return arma::solve(arma::trimatu(U), arma::solve(arma::trimatl(U.t()),h)) +
    arma::solve(arma::trimatu(U),noise);
}
// Exact one-sided normal rejection sampling; exponential tails use Robert's envelope.
double positive_normal(double mean) {
  if (mean >= 0) {
    double x; do {x=R::rnorm(mean,1);} while (x <= 0);
    return x;
  }
  double a=-mean, rate=.5*a+.5*std::hypot(a,2.0);
  for (;;) {
    double excess=R::rexp(1/rate);
    if (excess > 0 && std::log(R::runif(0,1)) <= -.5*std::pow(a-rate+excess,2)) return excess;
  }
}
struct PathContext {
  int P,I,K; bool time;
  std::vector<double> x,y,base,ls,rho,cut;
  std::vector<int> item,state;
  PathContext(List data,List e,NumericMatrix C,bool update_time):time(update_time) {
    NumericMatrix X=data["x"],Y=data["y"];
    P=X.nrow();I=X.ncol(); NumericVector a=e["a"],mu=e["mu"],b=e["b"],nu=e["nu"],d=e["d"],eta=e["eta"];
    NumericVector log_sigma=e["log_sigma"],r=e["rho"];K=a.size();
    if (C.nrow()!=P || C.ncol()!=I || Y.nrow()!=P || Y.ncol()!=I ||
        b.size()!=I || nu.size()!=I || d.size()!=P || eta.size()!=P ||
        mu.size()!=K || log_sigma.size()!=K || r.size()!=K) stop("Invalid emission dimensions.");
    bool mapped=data.containsElementNamed("item_index");
    NumericMatrix map; if(mapped) { map=as<NumericMatrix>(data["item_index"]); if(map.nrow()!=P||map.ncol()!=I)stop("Invalid item map."); }
    int N=P*I;x.resize(N);y.resize(N);base.resize(N);ls.resize(N);rho.resize(N);cut.resize(N);item.resize(N);state.resize(N);
    for(int n=0;n<N;++n) {
      int p=n%P,i=mapped?index(map[n],I):n/P,s=index(C[n],K);
      x[n]=X[n];y[n]=Y[n];item[n]=i;state[n]=s;ls[n]=log_sigma[s];rho[n]=r[s];
      if(time) { double ac=b[i]+a[s]+d[p];base[n]=-softplus(x[n]==1?-ac:ac);cut[n]=r[s]==0?0:threshold(ac); }
      else { double z=(Y[n]-nu[i]-mu[s]-eta[p])/std::exp(ls[n]);
        base[n]=R::dnorm(z,0,1,true)-ls[n];cut[n]=r[s]*z; }
    }
  }
};
struct ScaleContext {
  std::vector<double> x,delta,cut,lq;
  ScaleContext(NumericVector X,NumericVector Y,NumericVector A,NumericVector M) {
    int N=X.size(); if(Y.size()!=N||A.size()!=N||M.size()!=N) stop("Invalid scale dimensions.");
    x.resize(N);delta.resize(N);cut.resize(N);lq.resize(N);
    for(int n=0;n<N;++n) { x[n]=X[n];delta[n]=Y[n]-M[n];cut[n]=threshold(A[n]);lq[n]=-softplus(X[n]==1?-A[n]:A[n]); }
  }
};
}

// [[Rcpp::export]]
NumericVector emission_density_cpp(NumericVector x,NumericVector y,NumericVector a,
 NumericVector m,NumericVector sigma,NumericVector rho) {
  int N=x.size(); if(y.size()!=N||a.size()!=N||m.size()!=N||sigma.size()!=N||rho.size()!=N)stop("Nonconformable emissions.");
  NumericVector out(N);
  for(int n=0;n<N;++n) out[n]=log_density(x[n],y[n],a[n],m[n],sigma[n],rho[n]);
  return out;
}

// [[Rcpp::export]]
SEXP emission_context_cpp(List data,List emission,NumericMatrix C,bool time) {
  XPtr<PathContext> ptr(new PathContext(data,emission,C,time),true);return ptr;
}

// [[Rcpp::export]]
double emission_margin_loglik_cpp(SEXP pointer,NumericVector item,NumericVector state,NumericVector person) {
  XPtr<PathContext> c(pointer);
  if(item.size()!=c->I||state.size()!=c->K||person.size()!=c->P)stop("Invalid margin effects.");
  double total=0;
  for(int n=0;n<c->P*c->I;++n) {
    double v=item[c->item[n]]+state[c->state[n]]+person[n%c->P], ans;
    if(!R_finite(v)) return R_NegInf;
    double r=c->rho[n];
    if(c->time) {
      double z=(c->y[n]-v)/std::exp(c->ls[n]);
      ans=R::dnorm(z,0,1,true)-c->ls[n];
      if(r==0) ans+=c->base[n];
      else { double cut=(c->cut[n]-r*z)/std::sqrt(1-r*r);
        ans+=R::pnorm(c->x[n]==1?-cut:cut,0,1,true,true); }
    } else {
      ans=c->base[n];
      if(r==0) ans-=softplus(c->x[n]==1?-v:v);
      else { double cut=(threshold(v)-c->cut[n])/std::sqrt(1-r*r);
        ans+=R::pnorm(c->x[n]==1?-cut:cut,0,1,true,true); }
    }
    total+=ans;
  }
  return total;
}

// [[Rcpp::export]]
NumericVector emission_array_cpp(List data,List e,IntegerVector states) {
  NumericMatrix X=data["x"],Y=data["y"];int P=X.nrow(),I=X.ncol(),N=P*I;
  NumericVector b=e["b"],a=e["a"],d=e["d"],nu=e["nu"],mu=e["mu"],eta=e["eta"],ls=e["log_sigma"],rho=e["rho"];
  bool mapped=data.containsElementNamed("item_index");NumericMatrix map;
  if(mapped)map=as<NumericMatrix>(data["item_index"]);
  if(Y.nrow()!=P||Y.ncol()!=I||b.size()!=I||nu.size()!=I||d.size()!=P||eta.size()!=P||
      mu.size()!=a.size()||ls.size()!=a.size()||rho.size()!=a.size()||
      (mapped&&(map.nrow()!=P||map.ncol()!=I)))stop("Invalid emission dimensions.");
  NumericVector out(N*states.size());out.attr("dim")=IntegerVector::create(P,I,states.size());
  std::vector<double> ac(N),tm(N);
  for(int n=0;n<N;++n) { int p=n%P,i=mapped?index(map[n],I):n/P;ac[n]=b[i]+d[p];tm[n]=nu[i]+eta[p]; }
  for(int k=0;k<states.size();++k) {
    checkUserInterrupt();int s=index(states[k],a.size());double sigma=std::exp(ls[s]);
    for(int n=0;n<N;++n)out[n+N*k]=log_density(X[n],Y[n],ac[n]+a[s],tm[n]+mu[s],sigma,rho[s]);
  }
  return out;
}

// [[Rcpp::export]]
SEXP emission_scale_context_cpp(NumericVector x,NumericVector y,NumericVector a,NumericVector m) {
  XPtr<ScaleContext> ptr(new ScaleContext(x,y,a,m),true);return ptr;
}

// [[Rcpp::export]]
double emission_scale_loglik_cpp(SEXP pointer,double ls,double rho) {
  XPtr<ScaleContext> c(pointer);double sigma=std::exp(ls);
  if(!R_finite(sigma)||sigma<=0||!R_finite(rho)||std::abs(rho)>=1)return R_NegInf;
  double total=0,denom=std::sqrt(1-rho*rho);
  for(unsigned n=0;n<c->x.size();++n) {
    double z=c->delta[n]/sigma,ans=R::dnorm(z,0,1,true)-ls;
    if(rho==0)ans+=c->lq[n];
    else {double cut=(c->cut[n]-rho*z)/denom;ans+=R::pnorm(c->x[n]==1?-cut:cut,0,1,true,true);}
    total+=ans;
  }
  return total;
}

// [[Rcpp::export]]
List transition_probabilities_cpp(List tr,NumericMatrix z) {
  NumericMatrix alpha=tr["alpha"],beta=tr["beta"];NumericVector a0=tr["alpha0"];
  int K=as<int>(tr["K"]),P=z.nrow(),L=alpha.ncol(),q=z.ncol();bool finite=as<std::string>(tr["mode"])=="finite";
  if(alpha.nrow()!=K||beta.nrow()!=q||beta.ncol()!=L||a0.size()!=L||L!=(finite?K-1:K))stop("Invalid transition dimensions.");
  NumericMatrix initial(P,K),residual(P,K+1);List transitions(P);
  for(int p=0;p<P;++p) {
    checkUserInterrupt();NumericMatrix row(K,K);std::vector<double> xb(L,0);
    for(int k=0;k<L;++k)for(int j=0;j<q;++j)xb[k]+=z(p,j)*beta(j,k);
    for(int r=-1;r<K;++r) {
      double log_remaining=0;
      for(int k=0;k<L;++k) {
        double a=(r<0?a0[k]:alpha(r,k))+xb[k];
        double mass=std::exp(log_remaining+R::pnorm(a,0,1,true,true));
        if(r<0)initial(p,k)=mass;else row(r,k)=mass;
        log_remaining+=R::pnorm(a,0,1,false,true);
      }
      double tail=std::exp(log_remaining);
      if(finite) {if(r<0)initial(p,K-1)=tail;else row(r,K-1)=tail;}
      else residual(p,r+1)=tail;
    }
    transitions[p]=row;
  }
  return List::create(_["initial"]=initial,_["transition"]=transitions,_["residual"]=residual);
}

// [[Rcpp::export]]
List transition_update_cpp(List input,NumericMatrix C,arma::mat z,double prior_mean,double prior_sd,
                          Nullable<LogicalMatrix> active = R_NilValue) {
  List tr=clone(input);int P=C.nrow(),I=C.ncol(),K=as<int>(tr["K"]);
  arma::mat alpha=as<arma::mat>(tr["alpha"]),beta=as<arma::mat>(tr["beta"]);
  arma::vec a0=as<arma::vec>(tr["alpha0"]);int L=alpha.n_cols,q=z.n_cols;
  double va=as<double>(tr["sigma_alpha2"]),vd=as<double>(tr["var_diagonal"]),md=as<double>(tr["mean_diagonal"]);
  if(z.n_rows!=(unsigned)P||alpha.n_rows!=(unsigned)K||beta.n_rows!=(unsigned)q||beta.n_cols!=(unsigned)L||a0.n_elem!=(unsigned)L)
    stop("Invalid transition-update dimensions.");
  if(!R_finite(va)||va<=0||!R_finite(vd)||vd<=0||!R_finite(md)||
      !R_finite(prior_mean)||!R_finite(prior_sd)||prior_sd<=0||
      !alpha.is_finite()||!a0.is_finite()||!beta.is_finite()||!z.is_finite())
    stop("Transition parameters must be finite with positive variances.");
  if(!L) {
    tr["sigma_alpha2"]=1/R::rgamma(1,1);tr["mean_diagonal"]=R::rnorm(prior_mean,prior_sd);
    tr["var_diagonal"]=1/R::rgamma(1,1);return tr;
  }
  arma::mat xb=z*beta,counts(K,L,arma::fill::zeros),sums(K,L,arma::fill::zeros);
  bool masked=active.isNotNull();LogicalMatrix mask;
  if(masked) {
    mask=LogicalMatrix(active.get());
    if(mask.nrow()!=P||mask.ncol()!=I)stop("Invalid transition mask dimensions.");
    for(int p=0;p<P;++p)for(int t=0;t<I;++t) {
      if(mask(p,t)==NA_LOGICAL)stop("Transition mask cannot contain NA.");
      if(t==0&&!mask(p,t))stop("Initial observations must remain active.");
      if(t>0&&!mask(p,t)&&C(p,t)!=C(p,t-1))stop("Only self-transitions can be forced.");
    }
  }
  arma::vec n0(L,arma::fill::zeros),s0(L,arma::fill::zeros);
  NumericMatrix w(P*I,L);
  // Same reached-stick event as Hoskovec: C[p,t] >= k, with sign C[p,t] == k.
  for(int p=0;p<P;++p) {
    checkUserInterrupt();
    for(int t=0;t<I;++t) {
      if(masked&&!mask(p,t))continue;
      int dest=index(C(p,t),K),origin=t?index(C(p,t-1),K):-1;
      for(int k=0;k<std::min(dest+1,L);++k) {
        double mean=(origin<0?a0[k]:alpha(origin,k))+xb(p,k);
        double value=k==dest?positive_normal(mean):-positive_normal(-mean);
        w(p+P*t,k)=value;
        if(origin<0){n0[k]++;s0[k]+=value-xb(p,k);}
        else {counts(origin,k)++;sums(origin,k)+=value-xb(p,k);}
      }
    }
  }
  for(int r=0;r<K;++r)for(int k=0;k<L;++k) {
    double pv=r==k?vd:va,pm=r==k?md:0,v=1/(1/pv+counts(r,k));
    alpha(r,k)=R::rnorm(v*(pm/pv+sums(r,k)),std::sqrt(v));
  }
  for(int k=0;k<L;++k){double v=1/(1/va+n0[k]);a0[k]=R::rnorm(v*s0[k],std::sqrt(v));}
  for(int k=0;k<L && q;++k) {
    arma::mat Q(q,q,arma::fill::eye);arma::vec h(q,arma::fill::zeros);
    for(int p=0;p<P;++p) {
      double n=0,sum=0;
      for(int t=0;t<I;++t)if((!masked||mask(p,t))&&C(p,t)>=k+1) {
        n++;sum+=w(p+P*t,k)-(t?alpha(index(C(p,t-1),K),k):a0[k]);
      }
      Q+=n*z.row(p).t()*z.row(p);h+=sum*z.row(p).t();
    }
    beta.col(k)=precision_draw(Q,h);
  }
  double sum_diag=0;for(int k=0;k<L;++k)sum_diag+=alpha(k,k);
  double v=1/(1/(prior_sd*prior_sd)+L/vd);
  md=R::rnorm(v*(prior_mean/(prior_sd*prior_sd)+sum_diag/vd),std::sqrt(v));
  double ss=0;for(int k=0;k<L;++k)ss+=std::pow(alpha(k,k)-md,2);
  vd=1/R::rgamma(1+L/2.0,1/(1+ss/2));
  ss=arma::dot(a0,a0);int n=L;
  for(int r=0;r<K;++r)for(int k=0;k<L;++k)if(r!=k){ss+=alpha(r,k)*alpha(r,k);n++;}
  va=1/R::rgamma(1+n/2.0,1/(1+ss/2));
  tr["alpha"]=alpha;tr["alpha0"]=a0;tr["beta"]=beta;
  tr["mean_diagonal"]=md;tr["var_diagonal"]=vd;tr["sigma_alpha2"]=va;return tr;
}

// [[Rcpp::export]]
NumericVector positive_normal_cpp(double mean,int n) {
  if(!R_finite(mean)||n<0)stop("Invalid truncated-normal request.");
  NumericVector out(n);for(int j=0;j<n;++j){if(j%4096==0)checkUserInterrupt();out[j]=positive_normal(mean);}return out;
}

// [[Rcpp::export]]
NumericMatrix beam_slices_cpp(NumericMatrix C,NumericMatrix initial,List transitions) {
  int P=C.nrow(),I=C.ncol(),K=initial.ncol();
  if(initial.nrow()!=P||transitions.size()!=P)stop("Invalid slice dimensions.");
  NumericMatrix slices(P,I);
  for(int p=0;p<P;++p){NumericMatrix tr=transitions[p];checkUserInterrupt();
    if(tr.nrow()!=K||tr.ncol()!=K)stop("Invalid transition matrix.");
    for(int t=0;t<I;++t) {
      int s=index(C(p,t),K);double mass=t?tr(index(C(p,t-1),K),s):initial(p,s);
      if(!R_finite(mass)||mass<=0)stop("Current path has zero transition mass.");
      double u=0;while(u==0)u=R::runif(0,mass);slices(p,t)=u;
    }
  }
  return slices;
}

// [[Rcpp::export]]
List gaussian_time_stats_cpp(NumericMatrix C,Nullable<NumericMatrix> item_index,
 NumericMatrix y,NumericVector weights,int K,bool reduced) {
  int P=C.nrow(),I=C.ncol(),D=I+K+P-(reduced?2:0);
  if(y.nrow()!=P||y.ncol()!=I||weights.size()!=P*I)stop("Invalid Gaussian dimensions.");
  bool mapped=item_index.isNotNull();NumericMatrix map;if(mapped)map=as<NumericMatrix>(item_index);
  if(mapped&&(map.nrow()!=P||map.ncol()!=I))stop("Invalid item map.");
  arma::mat Q(D,D,arma::fill::zeros);arma::vec h(D,arma::fill::zeros);
  for(int n=0;n<P*I;++n) {
    if(n%4096==0)checkUserInterrupt();int p=n%P,i=mapped?index(map[n],I):n/P,s=index(C[n],K);
    int cols[4],count=0;
    if(reduced){cols[count++]=0;if(i)cols[count++]=i;if(s)cols[count++]=I+s-1;if(p)cols[count++]=I+K+p-2;}
    else {cols[count++]=i;cols[count++]=I+s;cols[count++]=I+K+p;}
    for(int a=0;a<count;++a) {h[cols[a]]+=weights[n]*y[n];for(int b=0;b<count;++b)Q(cols[a],cols[b])+=weights[n];}
  }
  return List::create(_["precision"]=Q,_["information"]=h);
}

// [[Rcpp::export]]
NumericMatrix location_covariance_cpp(int I,int K,int P,NumericVector variance) {
  int D=I+K+P-2;if(I<1||K<1||P<1||variance.size()!=3)stop("Invalid location dimensions.");
  NumericMatrix V(D,D);V(0,0)=sum(variance);int offset=1;
  int sizes[3]={I,K,P};
  for(int block=0;block<3;++block){int n=sizes[block]-1;
    for(int j=0;j<n;++j){V(0,offset+j)=V(offset+j,0)=-variance[block];
      for(int k=0;k<n;++k)V(offset+j,offset+k)=variance[block]*(1+(j==k));}
    offset+=n;
  }
  return V;
}

namespace {
struct StateView {
  SEXP value; int type;
  StateView(SEXP x,int P,int I):value(x),type(TYPEOF(x)) {
    if(!Rf_isMatrix(x)||(type!=INTSXP&&type!=REALSXP))stop("Invalid states matrix.");
    IntegerVector dims=Rf_getAttrib(x,R_DimSymbol);
    if(dims[0]!=P||dims[1]!=I)stop("States dimensions do not match data.");
  }
  int at(int n) const {
    if(type==INTSXP) {
      int x=INTEGER(value)[n];
      if(x<1)stop("Each retained states matrix must contain positive integer labels.");
      return x;
    }
    double x=REAL(value)[n];
    if(!R_finite(x)||x<1||x!=std::floor(x)||x>INT_MAX)
      stop("Each retained states matrix must contain positive integer labels.");
    return static_cast<int>(x);
  }
};
}

// [[Rcpp::export(rng = false)]]
IntegerMatrix posterior_partitions_cpp(List draws,int P,int I) {
  int M=draws.size();if(P<1||I<1||M<1)stop("Expected nonempty posterior draws.");
  IntegerMatrix Z(M,P*I);
  // Block draws so the destination writes are contiguous without a second big matrix.
  for(int start=0;start<M;start+=64) {
    checkUserInterrupt();int end=std::min(start+64,M);std::vector<StateView> state;
    for(int m=start;m<end;++m){List draw=draws[m];state.emplace_back(draw["states"],P,I);}
    for(int p=0;p<P;++p)for(int t=0;t<I;++t)for(int m=start;m<end;++m)
      Z(m,p*I+t)=state[m-start].at(p+P*t);
  }
  return Z;
}

// [[Rcpp::export(rng = false)]]
bool partitions_match_cpp(IntegerMatrix Z,List draws,int P,int I,IntegerVector rows) {
  if(Z.ncol()!=P*I||rows.size()!=draws.size())return false;
  for(int m=0;m<draws.size();++m) {
    if(m%128==0)checkUserInterrupt();int r=index(rows[m],Z.nrow());
    List draw=draws[m];StateView C(draw["states"],P,I);
    for(int p=0;p<P;++p)for(int t=0;t<I;++t)if(Z(r,p*I+t)!=C.at(p+P*t))return false;
  }
  return true;
}

// [[Rcpp::export(rng = false)]]
LogicalVector core_draws_cpp(IntegerMatrix Z,IntegerVector core,IntegerVector rows) {
  if(core.size()!=Z.ncol())stop("Core dimensions do not match partitions.");
  std::vector<int> keep,targets;
  for(int n=0;n<core.size();++n)if(core[n]>=0){keep.push_back(n);targets.push_back(core[n]);}
  LogicalVector out(rows.size(),true);
  for(int m=0;m<rows.size();++m) {
    if(m%128==0)checkUserInterrupt();int r=index(rows[m],Z.nrow());
    std::unordered_map<int,int> forward,reverse;
    for(unsigned j=0;j<keep.size();++j) {
      int a=targets[j],b=Z(r,keep[j]);
      if(b==NA_INTEGER||b<1)stop("Invalid posterior state label.");
      auto f=forward.find(a),v=reverse.find(b);
      if((f!=forward.end()&&f->second!=b)||(v!=reverse.end()&&v->second!=a)){out[m]=false;break;}
      forward[a]=b;reverse[b]=a;
    }
  }
  return out;
}

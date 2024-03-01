#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
double
hazdist_cpp(double lambda0,
            double sigma,
            double d) {
  double
  haz = lambda0* std::exp(-((d * d)/(2 * (sigma * sigma))));
  return haz;
}

// [[Rcpp::export]]
double 
  distkxt_cpp(int k,
              int x,
              Rcpp::DatetimeVector t,
              Rcpp::List dist_dat,
              double timesnap = 600
  ){
    //referencing a distance matrix object
    arma::cube
    distmat = dist_dat["distmat"];
    Rcpp::DatetimeVector
    times = dist_dat["times"];
    Rcpp::NumericVector
    timediffs = Rcpp::na_omit(abs(times - Rcpp::rep(t, times.length())));
    int
    closesttime = Rcpp::which_min(timediffs);
    int 
    tindex = NA_INTEGER;
    if(timediffs(closesttime) <= timesnap){
      tindex = closesttime;
    } else {
      Rcpp::stop("Error: Closest time in more than timesnap seconds away");
    }
    double
      distout = distmat(k - 1, x - 1, tindex); //0 index to 1 index
    return distout;
  }

// [[Rcpp::export]]
double
  haz_cpp(Rcpp::DatetimeVector t,
      int x, 
      int k,
      double lambda0,
      double sigma,
      Rcpp::List dist_dat,
      double timesnap = 600
  ){
    double
    distkxt_cpp_ = distkxt_cpp(k, x, t, dist_dat, timesnap);
    double
      hazout = hazdist_cpp(lambda0, sigma, distkxt_cpp_);
    return(hazout);
  }

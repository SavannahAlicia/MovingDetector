#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
Rcpp::NumericVector
hazdist_cpp(double lambda0,
            double sigma,
            double d) {
  Rcpp::NumericVector
  haz(1);
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
    timediffs = abs(times - t);
    int
    closesttime = which_min(na_omit(timediffs));
    int 
    tindex = NA_INTEGER;
    if(timediffs(closesttime) <= timesnap){
      tindex = closesttime;
    } else {
      Rcpp::stop("Error: Closest time in more than timesnap seconds away");
    }
    double
      distout = distmat(k, x, tindex);
    return distout;
  }


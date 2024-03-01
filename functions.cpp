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
              Rcpp::Datetime t,
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
  haz_cpp(Rcpp::Datetime t,
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

// [[Rcpp::export]]
double
  surv_cpp(Rcpp::Datetime timestart,
           Rcpp::Datetime timeend,
           int timeincr,
           int x,
           int k,
           double lambda0,
           double sigma, 
           Rcpp::List dist_dat,
           double timesnap = 600
  ){
    double survout;
    if(timeend < timestart) {
      survout = NA_REAL;
    } else {
      //figure out how many steps of length timeincr between timestart and timeend
      Rcpp::NumericVector
        steps = {0, ((timeend - timestart)/timeincr)};
      //create sequence of 1:number of steps * timeincr
      Rcpp::IntegerVector 
        timesteps = (Rcpp::seq(Rcpp::min(steps), Rcpp::max(steps)) * timeincr);
      //add that to timestart. These are timesopen
      Rcpp::DatetimeVector 
        timesopen(timesteps.length());
      for(int i = 0; i < timesteps.length(); i++){
        timesopen(i) = timestart + timesteps(i);
      }
      Rcpp::NumericVector 
        hazs(timesopen.length());
      for(int tt = 0; tt < timesopen.length(); tt++){
        Rcpp::Datetime
        timet = timesopen(tt);
        hazs(tt) = haz_cpp(timet, x, k, lambda0, sigma, dist_dat);
      }
      double 
        integ = Rcpp::sum(hazs);
      survout = std::exp(-1 * integ);
    }
    return(survout);
  }
  

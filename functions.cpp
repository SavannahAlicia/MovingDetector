#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::plugins(cpp11)]]

//----------------- Supporting functions ---------------------------------------
// [[Rcpp::export]]
Rcpp::NumericVector which_is_not_na(Rcpp::NumericVector x) {
  Rcpp::NumericVector indices;
  for(int i = 0; i < x.size(); ++i) {
    if(Rcpp::NumericVector::is_na(x[i])) {}else{
      indices.push_back(i); 
    }
  }
  return indices;
}


// [[Rcpp::export]]
Rcpp::NumericVector seqC(int &first, int &last) {
  Rcpp::NumericVector y(abs(last - first) + 1);
  if (first < last) 
    std::iota(y.begin(), y.end(), first);
  else {
    std::iota(y.begin(), y.end(), last);
    std::reverse(y.begin(), y.end());
  }
  return y;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix cubeRowToNumericMatrix(arma::cube cube,
                                           int row) {
  // Get dimensions of the cube
  int n_cols = cube.n_cols;
  int n_slices = cube.n_slices;
  
  // Create a NumericMatrix to store the row
  Rcpp::NumericMatrix result(n_cols, n_slices);
  
  // Copy the row from the cube to the NumericMatrix
  for (int col = 0; col < n_cols; ++col) {
    for (int slice = 0; slice < n_slices; ++slice){
      result(col, slice) = cube(row, col, slice);
    }
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector Sugar_colSums(const Rcpp::NumericMatrix x) {
  return colSums(x);
}

// [[Rcpp::export]]
double product(Rcpp::NumericVector x) {
  double prod = 1.0;
  for(int i = 0; i < x.size(); ++i) {
    prod *= x[i];
  }
  return prod;
}

// [[Rcpp::export]]
double sumC(Rcpp::NumericVector x) {
  int n = x.size();
  double total = 0;
  for(int i = 0; i < n; ++i) {
    total += x[i];
  }
  return total;
}

// [[Rcpp::export]]
Rcpp::NumericVector logvec(Rcpp::NumericVector x) {
  int n = x.size();
  Rcpp::NumericVector newvec(n);
  for(int i = 0; i < n; i++){
    newvec(i) = log(x(i));
  }
  return newvec;
}

// [[Rcpp::export]]
double vec_min(Rcpp::NumericVector x) {
  // Rcpp supports STL-style iterators
  Rcpp::NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference 
  return *it;
}

// [[Rcpp::export]]
double vec_max(Rcpp::NumericVector x) {
  // Rcpp supports STL-style iterators
  Rcpp::NumericVector::iterator it = std::max_element(x.begin(), x.end());
  // we want the value so dereference 
  return *it;
}

//----------------- Likelihood related functions -------------------------------
// [[Rcpp::export]]
double hazdist_cpp(double lambda0,
                   double sigma,
                   double d) {
  double haz = lambda0* std::exp(-((d * d)/(2 * (sigma * sigma))));
  return haz;
}

// [[Rcpp::export]]
double distkxt_cpp(int k,
                   int x,
                   Rcpp::Datetime t,
                   Rcpp::List dist_dat,
                   double timesnap = 600) {
  //referencing a distance matrix object
  arma::cube distmat = dist_dat["distmat"];
  Rcpp::DatetimeVector times = dist_dat["times"];
  Rcpp::NumericVector timediffs = Rcpp::na_omit(Rcpp::abs(times - Rcpp::rep(t, times.length())));
  int closesttime = Rcpp::which_min(timediffs);
  int tindex = NA_INTEGER;
  if(timediffs(closesttime) <= timesnap){
    tindex = closesttime;
  } else {
    Rcpp::stop("Error: Closest time in more than timesnap seconds away");
  }
  double distout = distmat(k, x, tindex); //remember 0 index
  return distout;
}

// [[Rcpp::export]]
double haz_cpp(Rcpp::Datetime t,
               int x, 
               int k,
               double lambda0,
               double sigma,
               Rcpp::List dist_dat,
               double timesnap = 600) {
  double distkxt_cpp_ = distkxt_cpp(k, x, t, dist_dat, timesnap);
  double hazout = hazdist_cpp(lambda0, sigma, distkxt_cpp_);
  return(hazout);
}

// [[Rcpp::export]]
double surv_cpp(Rcpp::Datetime timestart,
                Rcpp::Datetime timeend,
                int timeincr,
                int x,
                int k,
                double lambda0,
                double sigma, 
                Rcpp::List dist_dat,
                double timesnap = 600) {
  double survout;
  if(timeend < timestart) {
    survout = NA_REAL;
  } else {
    //figure out how many steps of length timeincr between timestart and timeend
    Rcpp::NumericVector steps = {0, ((timeend - timestart)/timeincr)};
    //create sequence of 1:number of steps * timeincr
    Rcpp::IntegerVector timesteps = (Rcpp::seq(Rcpp::min(steps), Rcpp::max(steps)) * timeincr);
    //add that to timestart. These are timesopen
    Rcpp::DatetimeVector timesopen(timesteps.length());
    for(int i = 0; i < timesteps.length(); i++){
      timesopen(i) = timestart + timesteps(i);
    }
    Rcpp::NumericVector hazs(timesopen.length());
    for(int tt = 0; tt < timesopen.length(); tt++){
      Rcpp::Datetime timet = timesopen(tt);
      hazs(tt) = haz_cpp(timet, x, k, lambda0, sigma, dist_dat);
    }
    //set intervals of time over which to approx integrate
    Rcpp::NumericVector opentimediffs(timesteps.length());
    for(int d = 1; d < opentimediffs.length(); d++){
      //double timeincr_db = static_cast<double>(Rcpp::clone(timeincr));
      opentimediffs(d) = timeincr;
    }
    //the gap between first time and itself is 0
    opentimediffs(0) = 0;
    
     // hazard times the time interval from the previous t to current t (0 for first time point)
    double integ = Rcpp::sum(hazs * opentimediffs);
    survout = std::exp(-1 * integ);
  }
  return(survout);
}



//-----------------Likelihood --------------------------------------------------
// [[Rcpp::export]]
double
negloglikelihood_cpp( //add log link 
  double lambda0,
  double sigma, 
  Rcpp::NumericVector D_mesh,
  int timeincr,
  Rcpp::NumericMatrix capthist,
  Rcpp::List dist_dat,
  double timesnap = 600) {//specify objects
  Rcpp::List trapspre = dist_dat["traps"];
  Rcpp::List traps = trapspre[0];
  arma::cube distmat = dist_dat["distmat"]; //traps x mesh x times
  Rcpp::List meshpre = dist_dat["mesh"];
  Rcpp::DatetimeVector times = dist_dat["times"];
  times.attr("tzn") = "UTC";
  //calculate mesh area (note this is for a rectangular mesh grid, will need to
  //be recalculated for irregular shaped mesh)
  Rcpp::NumericVector meshx = meshpre["x"];
  Rcpp::NumericVector meshy = meshpre["y"];
  Rcpp::NumericVector meshxsorted = Rcpp::sort_unique(meshx);
  Rcpp::NumericVector meshysorted = Rcpp::sort_unique(meshy);
  double mesharea = ((meshxsorted(2) - meshxsorted(1)) * (meshysorted(2) - meshysorted(1)))/10000; //hectares
  //begin for loops for lambdan calculation
  Rcpp::NumericVector Dx_pdotxs(meshx.length());
  for(int m = 0; m < meshx.length(); m++){
    double Dx = D_mesh(m);
    Rcpp::NumericVector surv_eachtrap((traps.length()));
    for(int trapk = 0; trapk < traps.length(); trapk++){ 
      Rcpp::NumericMatrix distmatslicek = cubeRowToNumericMatrix(distmat, trapk); //mesh x times
      Rcpp::NumericVector opentimeindx = which_is_not_na(Sugar_colSums(distmatslicek));
      double openidx = vec_min(opentimeindx);
      double closeidx = vec_max(opentimeindx);
      Rcpp::Datetime topentime = times[openidx];
      Rcpp::Datetime tclosetime = times[closeidx];
      surv_eachtrap(trapk) = surv_cpp(topentime, tclosetime, timeincr, m, trapk, lambda0, sigma, dist_dat); 
    }
    double survalltraps = product(surv_eachtrap);
    double pdot = 1 - survalltraps;
    double Dx_pdotx = Dx * pdot;
    Dx_pdotxs(m) = Dx_pdotx;
  }
  double lambdan = sum(Dx_pdotxs) * mesharea;
  //rest of likelihood
  double n = capthist.nrow();
  Rcpp::NumericVector integral_eachi(n);
  for(int i = 0; i < n; i++){
    Rcpp::NumericVector DKprod_eachx(meshx.length());
    for(int x = 0; x < meshx.length(); x++){
      Rcpp::NumericVector Sxhx_eachtrap(traps.length());
      Rcpp::NumericVector captik(1);
      for(int trapk = 0; trapk < traps.length(); trapk++){
        Rcpp::NumericMatrix distmatslicek = cubeRowToNumericMatrix(distmat, trapk);
        Rcpp::NumericVector opentimeindx = which_is_not_na(Sugar_colSums(distmatslicek));
        double openidx = vec_min(opentimeindx);
        double closeidx = vec_max(opentimeindx);
        Rcpp::Datetime starttime = times[openidx];
        Rcpp::Datetime endtime = times[closeidx];
        captik = capthist(i,trapk);
        bool ikcaught = all(Rcpp::is_na(captik)).is_false();//isna returns false, so a capture time exists
        if(ikcaught){
          Rcpp::Datetime dettime = capthist(i,trapk);
          double Sx = surv_cpp(starttime, dettime, timeincr, x, trapk, lambda0, sigma, dist_dat);
          double hx = haz_cpp(dettime, x, trapk, lambda0, sigma, dist_dat);
          Sxhx_eachtrap(trapk) = Sx * hx;
        }else{
          double Sx = surv_cpp(starttime, endtime, timeincr, x, trapk, lambda0, sigma, dist_dat);
          Sxhx_eachtrap(trapk) = Sx;
        }
      }
      double Sxhx_alltraps = product(Sxhx_eachtrap);
      DKprod_eachx(x) = D_mesh(x) * Sxhx_alltraps;
    }
    double DKprod_sum = sumC(DKprod_eachx);
    integral_eachi(i) = DKprod_sum * mesharea;
  }
  int n_int = std::round(n);
  int one = 1;
  Rcpp::NumericVector ns(n);
  ns = seqC(one, n_int);
  Rcpp::NumericVector logns(n);
  logns = log(ns);
  double lognfact = sum(logns);
  Rcpp::NumericVector logint = logvec(integral_eachi);
  double sumlogint = sumC(logint);
  double out = -1 * (-lambdan - lognfact + sumlogint); 
  return(out);
}

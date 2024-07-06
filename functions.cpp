#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppClock)]]
//[[Rcpp::plugins(cpp11)]]
#include <RcppClock.h>

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

// [[Rcpp::export]]
Rcpp::Datetime addSecondsToDatetime(Rcpp::Datetime dt, int secondsToAdd) {
  // Extract the double value from the Rcpp::Datetime object
  double dt_double = dt;
  
  // Convert secondsToAdd to double (in case there's any implicit conversion needed)
  double seconds = static_cast<double>(secondsToAdd);
  
  // Add the integer seconds to the double value
  double new_dt_double = dt_double + seconds;
  
  // Create a new Rcpp::Datetime object from the modified double value
  Rcpp::Datetime new_dt(new_dt_double);
  
  return new_dt;
}

// [[Rcpp::export]]
double extract_index(Rcpp::DatetimeVector x, Rcpp::Datetime v) {
  // Find the indices where (x - v) == 0
  int index = -1;
  // Find the indices where (x - v) == 0
  for (int i = 0; i < x.size(); ++i) {
    double vtimestamp = v.getFractionalTimestamp();
    Rcpp::Datetime xi = x[i];
    double xitimestamp = xi.getFractionalTimestamp();
    if (xitimestamp  == vtimestamp) {
      index = i;
      break;
    }
  }
  if (index != -1) {
    return index;
  } else {
    return NA_REAL; // Use NA_REAL for NA_DATETIME
  }
}

//----------------- Likelihood related functions -------------------------------
// [[Rcpp::export]]
double hazdist_cpp(double lambda0,
                   double sigma,
                   double d,
                   int timeincr) {
  double g0 = lambda0* std::exp(-((d * d)/(2 * (sigma * sigma))));
  //experimenting, not sure what to do for 1/T, seems somewhat arbitrary
  double haz  = log(1 - g0) * -1 *  1/timeincr ; //
  return haz;
}

// [[Rcpp::export]]
double disttxk_cpp(
                   Rcpp::Datetime t,
                   int x,
                   int k,
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
               int timeincr,
               double timesnap = 600) {
  double distkxt_cpp_ = disttxk_cpp(t, x, k, dist_dat, timesnap);
  double hazout = hazdist_cpp(lambda0, sigma, distkxt_cpp_, timeincr);
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
      hazs(tt) = haz_cpp(timet, x, k, lambda0, sigma, dist_dat, timeincr);
    }
    //set intervals of time over which to approx integrate
    Rcpp::NumericVector opentimediffs(timesteps.length());
    //the gap between first time and itself is 0
    opentimediffs(0) = 0;
    for(int d = 1; d < opentimediffs.length(); d++){
      opentimediffs(d) = timeincr;
    }
     // hazard times the time interval from the previous t to current t (0 for first time point)
    double integ = Rcpp::sum(hazs * opentimediffs);//will this preserve small numbers?
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
  Rcpp::Nullable<Rcpp::NumericVector> haz_kmt_ = R_NilValue,//traps x mesh x times
  Rcpp::Nullable<Rcpp::NumericVector> surv_kmt_ = R_NilValue,//traps x mesh x times
  double timesnap = 600) {//specify objects
  Rcpp::Clock clock;
  clock.tick("whole_enchilada");
  clock.tick("setup");
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

  Rcpp::NumericVector haz_kmt(traps.length() * meshx.length() * distmat.n_slices);
  Rcpp::NumericVector surv_kmt(traps.length() * meshx.length() * distmat.n_slices);
  Rcpp::NumericVector sumhaz_kmt(traps.length() * meshx.length() * distmat.n_slices);
  Rcpp::NumericVector incsurv_kmt(traps.length() * meshx.length() * distmat.n_slices);
  clock.tock("setup");
  //start with matrix of hazards
  clock.tick("matrices_creation");
  if(haz_kmt_.isNotNull()){
    haz_kmt = haz_kmt_;
  } else{
      for(int trapk = 0; trapk < traps.length(); trapk++){
        Rcpp::NumericMatrix distmatslicek = cubeRowToNumericMatrix(distmat, trapk); //mesh x times
        Rcpp::NumericVector opentimeindx = which_is_not_na(Sugar_colSums(distmatslicek));
        double openidx = vec_min(opentimeindx);
        double closeidx = vec_max(opentimeindx);
        for (int m = 0; m < meshx.length(); m++){
          for(int t = 0; t < distmat.n_slices; t++){
            if(t >= openidx & t <= closeidx){ 
              haz_kmt[trapk + traps.length() * (m + meshx.length() * t)] = haz_cpp(times[t], m, trapk, lambda0, sigma, dist_dat, timeincr);
              if(t == openidx){
                //this is the sum of the hazard to approx the integral for survival.
                sumhaz_kmt[trapk + traps.length() * (m + meshx.length() * t)] = 0;
              } else {
                sumhaz_kmt[trapk + traps.length() * (m + meshx.length() * t)] = haz_kmt[trapk + traps.length() * (m + meshx.length() * (t))] + sumhaz_kmt[trapk + traps.length() * (m + meshx.length() * (t - 1))];
              }
            } else{//hazard when trap isn't open is 0
              haz_kmt[trapk + traps.length() * (m + meshx.length() * t)] = 0;
              sumhaz_kmt[trapk + traps.length() * (m + meshx.length() * t)] = 0;
            }
          }
      }
    }
      haz_kmt.attr("dim") = Rcpp::Dimension(traps.length(), meshx.length(), distmat.n_slices);
    sumhaz_kmt.attr("dim") = Rcpp::Dimension(traps.length(), meshx.length(), distmat.n_slices);
  }
  if (surv_kmt_.isNotNull()){
    surv_kmt = surv_kmt_;       // casting to underlying type NumericMatrix?
  } else {
    for(int trapk = 0; trapk < traps.length(); trapk++){
      Rcpp::NumericMatrix distmatslicek = cubeRowToNumericMatrix(distmat, trapk); //mesh x times
      Rcpp::NumericVector opentimeindx = which_is_not_na(Sugar_colSums(distmatslicek));
      double openidx = vec_min(opentimeindx);
      double closeidx = vec_max(opentimeindx);
      for (int m = 0; m < meshx.length(); m++){
        for(int t = 0; t < distmat.n_slices; t++){
          if(t >= openidx & t <= closeidx){ //survival when trap isn't open is 1
            if(t == openidx) {
              surv_kmt[trapk + traps.length() * (m + meshx.length() * t)] = 1;
              incsurv_kmt[trapk + traps.length() * (m + meshx.length() * t)] = 1;
            } else {
              surv_kmt[trapk + traps.length() * (m + meshx.length() * t)] = exp(-1 * sumhaz_kmt[trapk + traps.length() * (m + meshx.length() * t)] * timeincr);
              incsurv_kmt[trapk + traps.length() * (m + meshx.length() * t)] = surv_cpp(times[(t-1)], times[t], timeincr, m, trapk, lambda0, sigma, dist_dat);
            }
          } else {
            //survival when trap is closed is 1
            surv_kmt[trapk + traps.length() * (m + meshx.length() * t)] = 1;
            incsurv_kmt[trapk + traps.length() * (m + meshx.length() * t)] = 1;
          }
        }
      }
    }
    surv_kmt.attr("dim") = Rcpp::Dimension(traps.length(), meshx.length(), distmat.n_slices);
    //then need to add references to this matrix below instead of calling haz and surv
  }
  clock.tock("matrices_creation");
  //begin for loops for lambdan calculation
  clock.tick("lambdan");
  Rcpp::NumericVector Dx_pdotxs(meshx.length());
  for(int m = 0; m < meshx.length(); m++){
    double Dx = D_mesh(m);
    Rcpp::NumericVector surv_eachtrap((traps.length()));
    for(int trapk = 0; trapk < traps.length(); trapk++){
      Rcpp::NumericMatrix distmatslicek = cubeRowToNumericMatrix(distmat, trapk); //mesh x times
      Rcpp::NumericVector opentimeindx = which_is_not_na(Sugar_colSums(distmatslicek));
      double closeidx = vec_max(opentimeindx);
      surv_eachtrap(trapk) = surv_kmt[trapk + traps.length() * (m + meshx.length() * closeidx)]; //trap x mesh x time
    }
    double survalltraps = product(surv_eachtrap);
    double pdot = 1 - survalltraps;
    double Dx_pdotx = Dx * pdot;
    Dx_pdotxs(m) = Dx_pdotx;
  }
  double lambdan = sum(Dx_pdotxs) * mesharea;
  clock.tock("lambdan");

  //rest of likelihood
  clock.tick("loopllk");
  double n = capthist.nrow();
  Rcpp::NumericVector integral_eachi(n);
  for(int i = 0; i < n; i++){
    clock.tick("inteachi" + std::to_string(i));
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
        captik = capthist(i,trapk);
        bool ikcaught = all(Rcpp::is_na(captik)).is_false();//isna returns false, so a capture time exists
        if(ikcaught){
          Rcpp::Datetime dettime = capthist(i,trapk);
          // suvival up to time interval in which capture occurs, then hazard times time interval
          //time of trap before dettime
          //Rcpp::Datetime beforedettime = addSecondsToDatetime(dettime, -timeincr);
          if (dettime == starttime){
            //if no time has passed, multiply by time increment of 0
            Sxhx_eachtrap(trapk) = 0 + 1e-16;
          } else {
            //#probability the first detection happens at t given at least one detection happens
            //double Sx = surv_cpp(starttime, dettime, timeincr, x, trapk, lambda0, sigma, dist_dat) + 1e-16;
            //double hx = haz_cpp(dettime, x, trapk, lambda0, sigma, dist_dat, timeincr);
            //Sxhx_eachtrap(trapk) = Sx * hx * timeincr;

            //#probability that at least one detection happens in increment ending at t given at least one detection happens
            //need to find index of beforedettime (or of dettime and just subtract 1)
            double dettimeindx = extract_index(times, dettime);
            double Sx = surv_kmt(trapk + traps.length() * (x + meshx.length() * (dettimeindx - 1)));
           // double Sx = surv_cpp(starttime, beforedettime, timeincr, x, trapk, lambda0, sigma, dist_dat) + 1e-16;
            double atleastoneseen = 1 - incsurv_kmt(trapk + traps.length() * (x + meshx.length() * dettimeindx));
            Sxhx_eachtrap(trapk) = Sx * atleastoneseen;
          }
        }else{
          double Sx = surv_kmt(trapk + traps.length() * (x + meshx.length() * closeidx));
          Sxhx_eachtrap(trapk) = Sx;
        }
      }
      double Sxhx_alltraps = product(Sxhx_eachtrap);
      DKprod_eachx(x) = D_mesh(x) * Sxhx_alltraps;
    }
    double DKprod_sum = sumC(DKprod_eachx);
    integral_eachi(i) = DKprod_sum * mesharea;
    clock.tock("inteachi" + std::to_string(i));
  }
  clock.tock("loopllk");
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
  clock.tock("whole_enchilada");
  clock.stop("llktimes");
  return(out);
}


//#---------------------Stationary Detector multi-catch LLK
//-----------------Likelihood --------------------------------------------------

// [[Rcpp::export]]
double halfnormdet_cpp(double lambda0,
                   double sigma,
                   double d) {
  double g0 = lambda0* std::exp(-((d * d)/(2 * (sigma * sigma))));
  return g0;
}

// [[Rcpp::export]]
double
  negloglikelihood_stationary_cpp( //add log link 
    double lambda0,
    double sigma, 
    Rcpp::NumericVector D_mesh,
    arma::cube capthist,
    Rcpp::NumericMatrix usage, //traps by occ
    Rcpp::NumericMatrix distmat, //traps x mesh 
    Rcpp::List dist_dat
    ) {//specify objects
    int one = 1;
    int captrap = capthist.n_slices;
    Rcpp::NumericVector traps = seqC(one, captrap); //needs to still just be a list of trap ID number
    int capocc = capthist.n_cols;
    Rcpp::NumericVector occs = seqC(one, capocc); //seq(1:n_occasions) but 0 base
    Rcpp::List meshpre = dist_dat["mesh"]; //still just an xy mesh
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
      Rcpp::NumericVector notseen_eachocc((occs.size()));
      for(int occk = 0; occk < occs.size(); occk++){
        Rcpp::NumericVector notseen_eachtrap((traps.size()));
        for(int trapj = 0; trapj < traps.size(); trapj++){
          if(usage(trapj, occk) == 0){
            notseen_eachtrap(trapj) = 1; //if the trap isn't used, can't be seen at it
          } else {
            double thisdist = distmat(trapj, m);
            notseen_eachtrap(trapj) = 1 - halfnormdet_cpp(lambda0, sigma, thisdist); //1 minus prob seen
          }
        }
        double notseenalltraps = product(notseen_eachtrap);
        notseen_eachocc(occk) = notseenalltraps;
      }
      double notseen_alloccs = product(notseen_eachocc);
      double pdot = 1 - notseen_alloccs;
      double Dx_pdotx = Dx * pdot;
      Dx_pdotxs(m) = Dx_pdotx;
     }
    double lambdan = sum(Dx_pdotxs) * mesharea;
    //rest of likelihood
    double n = capthist.n_rows;
    Rcpp::NumericVector integral_eachi(n);
    for(int i = 0; i < n; i++){
      Rcpp::NumericVector DKprod_eachx(meshx.length());
      for(int x = 0; x < meshx.length(); x++){
        Rcpp::NumericVector probcapthist_eachocc(occs.length());
        for(int occk = 0; occk < occs.length(); occk++){
          bool ikcaught = FALSE;
          int trapijk;
          Rcpp::NumericVector pj_xis(traps.length()); //need to calculate prod(1-pji)
          Rcpp::NumericVector oneminuspj_xis(traps.length());  //and sum(pjk)
          //need to limit this to traps used in the occasion
          for(int trapj = 0; trapj < traps.length(); trapj++){
            double captik = capthist(i, occk, trapj);
            if(usage(trapj, occk) == 0){
              pj_xis(trapj) = 0;
              oneminuspj_xis(trapj) = 1;
            } else{
              double thisdist2 = distmat(trapj, x);
              pj_xis(trapj) = halfnormdet_cpp(lambda0, sigma, thisdist2);
              oneminuspj_xis(trapj) = 1 - halfnormdet_cpp(lambda0, sigma, thisdist2);
            }
            //assign if i is caught and which trap caught it
            if(captik == 1){ // this could be within the above else (only captures if used)?
              ikcaught = TRUE;
              trapijk = trapj;
            }
          }
          double prod_oneminuspj_xis = product(oneminuspj_xis);
          if(ikcaught){
            double sum_pj_xis = sum(pj_xis);
            double probj = halfnormdet_cpp(lambda0, sigma, distmat(trapijk, x));
            probcapthist_eachocc(occk) = (probj/sum_pj_xis) * (1 - prod_oneminuspj_xis);
          }else{
            probcapthist_eachocc(occk) = prod_oneminuspj_xis;
          }
        }
        double probcapthist_alloccs = product(probcapthist_eachocc);
        DKprod_eachx(x) = D_mesh(x) * probcapthist_alloccs;
      }
      double DKprod_sum = sumC(DKprod_eachx);
      integral_eachi(i) = DKprod_sum * mesharea;
    }
    int n_int = std::round(n);
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


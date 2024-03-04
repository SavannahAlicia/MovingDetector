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
      double integ = Rcpp::sum(hazs);
      survout = std::exp(-1 * integ);
    }
    return(survout);
  }

//-----------------Likelihood --------------------------------------------------
// [[Rcpp::export]]
double negloglikelihood_cpp(
    double lambda0,
    double sigma, 
    Rcpp::NumericVector D_mesh,
    int timeincr,
    Rcpp::NumericMatrix capthist,
    Rcpp::List dist_dat,
    double timesnap = 600) {//specify objects
    Rcpp::List trapspre = dist_dat["traps"];
    Rcpp::List traps = trapspre[0];
    arma::cube distmat = dist_dat["distmat"];
    Rcpp::List meshpre = dist_dat["mesh"];
    Rcpp::DatetimeVector times = dist_dat["times"];
    times.attr("tzn") = "UTC";
    //calculate mesh area
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
        Rcpp::NumericMatrix distmatslicek = cubeRowToNumericMatrix(distmat, trapk);
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
    return(lambdan);
  }

// double out;

// negloglikelihood_RTMB <- function(pars){
//function(lambda0, sigma, D_mesh, timeincr, capthist, dist_dat){
//   getAll(dist_dat)
//   lambda0 <- exp(pars$loglambda0)
//   sigma <- exp(pars$logsigma)
//   D_mesh <- rep(exp(pars$logD), 144)
//   out <- as.numeric(0)
//   
//   meshx_array <- as.array(1:nrow(mesh))
//     Dx_pdotxs <- rep(0, length(meshx_array))
//     for (meshx in meshx_array){
//       Dx <- as.numeric(D_mesh[meshx])
//       surv_eachtrap <- rep(0, nrow(traps))
//       for(trapk in 1:nrow(traps)){
//         opentimeindx <- which(!is.na(colSums(distmat[trapk,,], na.rm = F)))
//         topentime <- times[min(opentimeindx)]
//         tclosetime <- times[max(opentimeindx)]
//         surv_eachtrap[trapk] <- as.numeric(surv_rtmb(timestart = topentime, 
//                                                      timeend = tclosetime, 
//                                                      timeincr = timeincr, 
//                                                      x = meshx, k = trapk, 
//                                                      lambda0 = lambda0, 
//                                                      sigma = sigma,
//                                                      distmat = distmat,
//                                                      times = times))
//       }
//       surv_alltraps <- prod(surv_eachtrap)
//         pdot <- 1 - surv_alltraps
//       Dx_pdotx <- Dx * pdot
//       Dx_pdotxs[meshx] <- (Dx_pdotx)
//     }
//     lambdan_ <- as.numeric(sum(Dx_pdotxs) * attr(mesh, "area"))
//       n <- nrow(capthist)
//       integral_eachi <- rep(0, n)
//       for (i in 1:n){
//         DKprod_eachx <- rep(0, nrow(mesh))
//         for (x in 1:nrow(mesh)){
//           Sxhx_eachtrap <- rep(0, nrow(traps))
//           for(trapk in 1:nrow(traps)){
//             opentimeindx <- which(!is.na(colSums(distmat[trapk,,], na.rm = F)))
//             starttime <- times[min(opentimeindx)]
//             if(!is.na(capthist[i,trapk])){
//               dettime <- capthist[i, trapk]
//               Sx <- as.numeric(surv_rtmb(starttime, dettime, timeincr, x, trapk, lambda0, sigma, distmat, times))
//               hx <- as.numeric(haz_rtmb(dettime, x, trapk, lambda0, sigma, distmat, times))
//               Sxhxout <- as.numeric(Sx * hx)
//             } else { #individual never detected by this trap
//               endtime <- times[max(opentimeindx)]
//               Sx <- as.numeric(surv_rtmb(starttime, endtime, timeincr, x, trapk, lambda0, sigma, distmat, times))
//               Sxhxout <- Sx
//             }
//             Sxhx_eachtrap[trapk] <- Sxhxout
//           }
//           Sxhx_alltraps <- prod(Sxhx_eachtrap)
//             DKprod_out <- as.numeric(D_mesh[x]) * Sxhx_alltraps 
//           DKprod_eachx[x] <- DKprod_out
//         }
//         integral_eachi[i] <- sum(DKprod_eachx) * attr(mesh, "area") #sum * mesh area
//       }
//       ns <- 1:n
//       logns <- log(ns)
//         lognfact <- sum(logns)
//         out <- (-lambdan_ - lognfact + n * sum(log(integral_eachi)))
//         ADREPORT(sigma) ## If I want estimates on the real scale.
//       ADREPORT(lambda0)
//         ADREPORT(D_mesh)
//         return(-out)
// }
// 


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
//----------------- Likelihood related functions -------------------------------

// [[Rcpp::export]]
double halfnormdet_cpp(double lambda0,
                       double sigma,
                       double d) {
  double g0 = lambda0* std::exp(-((d * d)/(2 * (sigma * sigma))));
  return g0;
}


// [[Rcpp::export]]
  double hazdist_cpp(double lambda0,
    double sigma,
    double d,
    int timeincr) {
double g0 = halfnormdet_cpp(lambda0, sigma, d);
//timeincr is the unit for hazard rate (could be time or distance (eg meters))
double haz  = log(1 - g0) * -1 *  1/timeincr ; //
return haz;
}
//#---------------------Moving Detector multi-catch LLK (approximated with traps)
//-----------------Likelihood --------------------------------------------------

// [[Rcpp::export]]
double
negloglikelihood_moving_cpp( //add log link 
                  double lambda0,
                  double sigma, 
                  int timeincr,
                  Rcpp::NumericVector D_mesh,
                  arma::cube capthist,
                  Rcpp::NumericMatrix usage, //traps by occ, could be calculated from indusage
                  arma::cube indusage, // ind by traps by occ
                  Rcpp::NumericMatrix distmat, //traps x mesh 
                  Rcpp::NumericMatrix mesh //first column x, second y
) {//specify objects
Rcpp::Clock clock;
clock.tick("wholeenchilada");
clock.tick("setup");
int one = 1;
int captrap = capthist.n_slices;
Rcpp::NumericVector traps = seqC(one, captrap); //needs to still just be a list of trap ID number
int capocc = capthist.n_cols;
Rcpp::NumericVector occs = seqC(one, capocc); //seq(1:n_occasions) but 0 base
//calculate mesh area (note this is for a rectangular mesh grid, will need to
          //be recalculated for irregular shaped mesh)
Rcpp::NumericVector meshx = mesh.column(0);
Rcpp::NumericVector meshy = mesh.column(1);
Rcpp::NumericVector meshxsorted = Rcpp::sort_unique(meshx);
Rcpp::NumericVector meshysorted = Rcpp::sort_unique(meshy);
double mesharea = ((meshxsorted(2) - meshxsorted(1)) * (meshysorted(2) - meshysorted(1)))/10000; //hectares
clock.tock("setup");
//begin for loops for lambdan calculation
clock.tick("lambdan");
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
double thisdist = distmat(trapj, m);  //note this can't be recycled below except for individuals never detected. detected inds will have different induse
notseen_eachtrap(trapj) = exp(-hazdist_cpp(lambda0, sigma, thisdist, timeincr) * (usage(trapj, occk)/timeincr)) +  1e-16; //survival
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
clock.tock("lambdan");
//rest of likelihood
clock.tick("loopllk");
double n = capthist.n_rows;
Rcpp::NumericVector integral_eachi(n);
for(int i = 0; i < n; i++){
Rcpp::NumericVector DKprod_eachx(meshx.length());
for(int x = 0; x < meshx.length(); x++){
Rcpp::NumericVector probcapthist_eachocc(occs.length());
for(int occk = 0; occk < occs.length(); occk++){
bool ikcaught = FALSE;
int trapijk;
Rcpp::NumericVector hu_js(traps.length()); // hazard times use for m, j, k 
for(int trapj = 0; trapj < traps.length(); trapj++){//could limit this to traps used in the occasion
double captik = capthist(i, occk, trapj);
  hu_js(trapj) = hazdist_cpp(lambda0, sigma, distmat(trapj, x), timeincr) * (indusage(i, trapj, occk)/timeincr);//hazard times individual usage for each trap. 
if(captik == 1){ // this could be within the above else (only captures if used)?
ikcaught = TRUE; //assign if i is caught and which trap caught it
trapijk = trapj;
}
}
if(ikcaught){
probcapthist_eachocc(occk) = exp(-sum(hu_js)) * hu_js(trapijk) ; //survived to time of detection (usage) and detected
}else{
probcapthist_eachocc(occk) = exp(-sum(hu_js)) ; //survived all traps
}
}
double probcapthist_alloccs = product(probcapthist_eachocc);
DKprod_eachx(x) = D_mesh(x) * probcapthist_alloccs;
}
double DKprod_sum = sumC(DKprod_eachx) + 1e-16;
integral_eachi(i) = DKprod_sum * mesharea;
}
clock.tock("loopllk");
int n_int = std::round(n);
Rcpp::NumericVector ns(n);
ns = seqC(one, n_int);
Rcpp::NumericVector logns(n);
logns = log(ns);
double lognfact = sum(logns);
Rcpp::NumericVector logint = logvec(integral_eachi);
double sumlogint = sumC(logint);
double out = -1 * (-lambdan - lognfact + sumlogint);
clock.tock("wholeenchilada");
clock.stop("approxllktimes");
return(out);
}

//---------------Stationary

// [[Rcpp::export]]
double
  negloglikelihood_stationary_cpp( //add log link 
    double lambda0,
    double sigma, 
    int timeincr,
    Rcpp::NumericVector D_mesh,
    arma::cube capthist,
    Rcpp::NumericMatrix usage, //traps by occ, could be calculated from indusage
    Rcpp::NumericMatrix distmat, //traps x mesh 
    Rcpp::NumericMatrix mesh //first column x, second y
  ) {//specify objects
    Rcpp::Clock clock;
    clock.tick("wholeenchilada");
    clock.tick("setup");
    int one = 1;
    int captrap = capthist.n_slices;
    Rcpp::NumericVector traps = seqC(one, captrap); //needs to still just be a list of trap ID number
    int capocc = capthist.n_cols;
    Rcpp::NumericVector occs = seqC(one, capocc); //seq(1:n_occasions) but 0 base
    //calculate mesh area (note this is for a rectangular mesh grid, will need to
    //be recalculated for irregular shaped mesh)
    Rcpp::NumericVector meshx = mesh.column(0);
    Rcpp::NumericVector meshy = mesh.column(1);
    Rcpp::NumericVector meshxsorted = Rcpp::sort_unique(meshx);
    Rcpp::NumericVector meshysorted = Rcpp::sort_unique(meshy);
    double mesharea = ((meshxsorted(2) - meshxsorted(1)) * (meshysorted(2) - meshysorted(1)))/10000; //hectares
    clock.tock("setup");
    //begin for loops for lambdan calculation
    clock.tick("lambdan");
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
            double thisdist = distmat(trapj, m);  //note this can't be recycled below except for individuals never detected. detected inds will have different induse
            notseen_eachtrap(trapj) = exp(-hazdist_cpp(lambda0, sigma, thisdist, timeincr) * (usage(trapj, occk)/timeincr)); //survival
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
    clock.tock("lambdan");
    //rest of likelihood
    clock.tick("loopllk");
    double n = capthist.n_rows;
    Rcpp::NumericVector integral_eachi(n);
    for(int i = 0; i < n; i++){
      Rcpp::NumericVector DKprod_eachx(meshx.length());
      for(int x = 0; x < meshx.length(); x++){
        Rcpp::NumericVector probcapthist_eachocc(occs.length());
        for(int occk = 0; occk < occs.length(); occk++){
          bool ikcaught = FALSE;
          int trapijk;
          Rcpp::NumericVector hu_js(traps.length()); // hazard times use for m, j, k 
          for(int trapj = 0; trapj < traps.length(); trapj++){//could limit this to traps used in the occasion
            double captik = capthist(i, occk, trapj);
            hu_js(trapj) = hazdist_cpp(lambda0, sigma, distmat(trapj, x), timeincr) * (usage(trapj, occk)/timeincr);//hazard times usage for each trap. 
            if(captik == 1){ // this could be within the above else (only captures if used)?
              ikcaught = TRUE; //assign if i is caught and which trap caught it
              trapijk = trapj;
            }
          }
          double sum_hujs = sum(hu_js);
          if(ikcaught){
            //prob that trap j made the detection given i was detected in k
            probcapthist_eachocc(occk) = hu_js(trapijk)/(sum_hujs) * (1 - exp(-sum_hujs)); 
          }else{
            //prob i wasn't detected in k
            probcapthist_eachocc(occk) = exp(-sum_hujs) ; //survived all traps
          }
        }
        double probcapthist_alloccs = product(probcapthist_eachocc);
        DKprod_eachx(x) = D_mesh(x) * probcapthist_alloccs;
      }
      double DKprod_sum = sumC(DKprod_eachx);
      integral_eachi(i) = DKprod_sum * mesharea;
    }
    clock.tock("loopllk");
    int n_int = std::round(n);
    Rcpp::NumericVector ns(n);
    ns = seqC(one, n_int);
    Rcpp::NumericVector logns(n);
    logns = log(ns);
    double lognfact = sum(logns);
    Rcpp::NumericVector logint = logvec(integral_eachi);
    double sumlogint = sumC(logint);
    double out = -1 * (-lambdan - lognfact + sumlogint);
    clock.tock("wholeenchilada");
    clock.stop("approxstatllktimes");
    return(out);
  }
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppClock)]]
//[[Rcpp::plugins(cpp11)]]
#include <RcppClock.h>
#include <algorithm>
#include <string>
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
Rcpp::NumericMatrix calc_dist_matC(Rcpp::NumericMatrix t,//traps (j) xy
                                   Rcpp::NumericMatrix h //hrcs (i) xy
                                     ) {
  // Distance matrix has i rows and j columns
  Rcpp::NumericMatrix dist_mat(t.nrow(), h.nrow());
  for(int i = 0; i < h.nrow(); ++i) {
    for(int j = 0; j < t.nrow(); ++j) {
      double dx = t(j,0) - h(i,0);
      double dy = t(j, 1) - h(i, 1);
      dist_mat(j,i) = std::sqrt(dx * dx + dy * dy);
    }
  }
  return dist_mat;
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

// ----- Functions for calculating tracklength effort ------

// [[Rcpp::export]]

Rcpp::List create_line_list_C(Rcpp::DataFrame tracksdf, 
                             std::string scenario = "everything") {
  
  // Check that required columns exist
  Rcpp::CharacterVector colnames = tracksdf.names();
  
  // Check required columns
  if (std::find(colnames.begin(), colnames.end(), "x") == colnames.end() ||
      std::find(colnames.begin(), colnames.end(), "y") == colnames.end() ||
      std::find(colnames.begin(), colnames.end(), "time") == colnames.end() ||
      std::find(colnames.begin(), colnames.end(), "occ") == colnames.end()) {
    Rcpp::stop("tracksdf must contain columns: x, y, occ, time");
  }
  
  //Check if "effort" column exists
  bool has_effort = (std::find(colnames.begin(), colnames.end(), "effort") != colnames.end());
  
  // If "effort" exists, extract it
  Rcpp::CharacterVector effort;
  if (has_effort) {
    effort = tracksdf["effort"];
  }
  
  Rcpp::NumericVector x = tracksdf["x"];
  Rcpp::NumericVector y = tracksdf["y"];
  Rcpp::IntegerVector occ = tracksdf["occ"];
  Rcpp::NumericVector time = tracksdf["time"];
  
  // Find unique track IDs
  Rcpp::IntegerVector uniq_occ = sort_unique(occ);
  int n_occ = uniq_occ.size();
  
  // Initialize output
  Rcpp::List out_list(n_occ);
  Rcpp::CharacterVector out_names(n_occ);
  
  // Column names for numeric matrices
  Rcpp::CharacterVector mat_colnames(5);
  mat_colnames[0] = "x1";
  mat_colnames[1] = "y1";
  mat_colnames[2] = "x2";
  mat_colnames[3] = "y2";
  mat_colnames[4] = "time";
  
  // Iterate over each unique track ID
  for (int k = 0; k < n_occ; ++k) {
    int track_id = uniq_occ[k];
    
    // Indices for this occasion from tracksdf (of points, not segments)
    Rcpp::LogicalVector mask = occ == track_id;
    Rcpp::IntegerVector idx = Rcpp::seq(0, occ.size() - 1); // indexes rows of entire tracksdf
    idx = idx[mask]; // subsetted to rows that are in this occasion
    
    // If there are less than two points in a track
    int n_pt = idx.size();
    if (n_pt < 2) {
      out_list[k] = R_NilValue;
      out_names[k] = std::to_string(track_id);
      continue;
    }
    
    // --- Determine number of segments to allocate ---
    int seg_count = 0;
    
    if (scenario == "onison") {
      if (!has_effort)
        Rcpp::stop("Onison scenario invalid without effort column");
      
      for (int i = 0; i < n_pt - 1; ++i) { // -1 makes sure there is a point to complete the segment
        if (effort[idx[i]] == "OnEffort") seg_count++; 
      }
    } else {
      seg_count = n_pt - 1; // one segment between each consecutive point
    }
    
    // --- Preallocate numeric vectors ---
    Rcpp::NumericMatrix seg_mat(seg_count, 4);
    
    // --- Fill them ---
    int seg_index = 0;
    for (int i = 0; i < n_pt - 1; ++i) { // n_pt in occasion
      bool keep_seg = (scenario == "everything") ||
        (scenario == "onison" && effort[idx[i]] == "OnEffort");
      
      if (keep_seg) {
        seg_mat(seg_index, 0) = x[idx[i]];       // x1
        seg_mat(seg_index, 1) = y[idx[i]];       // y1
        seg_mat(seg_index, 2) = x[idx[i + 1]];   // x2
        seg_mat(seg_index, 3) = y[idx[i + 1]];   // y2
        seg_mat(seg_index, 4) = time[idk[i]];
        seg_index++;
      }
    }
    
    Rcpp::List dimnames(2);
    dimnames[0] = Rcpp::CharacterVector(seg_count); // empty row names
    dimnames[1] = mat_colnames;                     // column names
    Rf_setAttrib(seg_mat, R_DimNamesSymbol, dimnames);
    
    out_list[k] = seg_mat;
    out_names[k] = std::to_string(track_id);
  }
  
  out_list.attr("names") = out_names;
  return out_list;
}
 
 // [[Rcpp::export]]
 Rcpp::NumericMatrix create_grid_bboxes_C(Rcpp::DataFrame grid,
                                             double spacing) {
   
   // Extract coordinates
   Rcpp::NumericVector x = grid["x"];
   Rcpp::NumericVector y = grid["y"];
   
   int n = x.size();
   if (spacing <= 0) Rcpp::stop("spacing must be positive");
   
   double h = spacing / 2.0;
   
   // Prepare output matrix: columns = left, right, bottom, top
   Rcpp::NumericMatrix bboxes(n, 4);
   
   for (int i = 0; i < n; ++i) {
     bboxes(i,0) = x[i] - h; // left
     bboxes(i,1) = x[i] + h; // right
     bboxes(i,2) = y[i] - h; // bottom
     bboxes(i,3) = y[i] + h; // top
   }
   
   Rcpp::CharacterVector rownames(bboxes.nrow()); 
   Rcpp::CharacterVector colnames(4);
   colnames[0] = "left";
   colnames[1] = "right";
   colnames[2] = "bottom";
   colnames[3] = "top";
   
   Rcpp::List dimnames(2);
   dimnames[0] = rownames; // row names
   dimnames[1] = colnames; // column names
   
   Rf_setAttrib(bboxes, R_DimNamesSymbol, dimnames);
   
   return bboxes;
 }
 
 // [[Rcpp::export]]
 double get_length_C(Rcpp::NumericVector line, 
                              Rcpp::NumericVector box) {
   
   // Constants for region codes
   const int INSIDE = 0; // 0000
   const int LEFT   = 1; // 0001
   const int RIGHT  = 2; // 0010
   const int BOTTOM = 4; // 0100
   const int TOP    = 8; // 1000
   
   // Extract coordinates
   double x1 = line[0];
   double y1 = line[1];
   double x2 = line[2];
   double y2 = line[3];
   
   double xmin = box[0];
   double xmax = box[1];
   double ymin = box[2];
   double ymax = box[3];
   
   // Function to compute region code
   auto computeCode = [&](double x, double y) {
     int code = INSIDE;
     if (x < xmin) {
       code = code | LEFT;
     }
     if (x > xmax) {
       code = code | RIGHT;
     }
     if (y < ymin) {
       code = code | BOTTOM;
     }
     if (y > ymax) {
       code = code | TOP;
     }
     return code;
   };
   
   int code1 = computeCode(x1, y1);
   int code2 = computeCode(x2, y2);
   bool accept = false;
   
   while (true) {
     if (code1 == 0 && code2 == 0) {
       // Both points inside
       accept = true;
       break;
     }
     else if ((code1 & code2) != 0) { //bitwise comparison &
       // Both points share an outside region
       break;
     }
     else {
       // Line needs clipping
       int codeOut;
       double x = 0.0;
       double y = 0.0;
       
       if (code1 != 0) {
         codeOut = code1;
       } else {
         codeOut = code2;
       }
       
       // Find intersection point
       if ((codeOut & TOP) != 0) {
         x = x1 + (x2 - x1) * (ymax - y1) / (y2 - y1);
         y = ymax;
       } 
       else if ((codeOut & BOTTOM) != 0) {
         x = x1 + (x2 - x1) * (ymin - y1) / (y2 - y1);
         y = ymin;
       } 
       else if ((codeOut & RIGHT) != 0) {
         y = y1 + (y2 - y1) * (xmax - x1) / (x2 - x1);
         x = xmax;
       } 
       else if ((codeOut & LEFT) != 0) {
         y = y1 + (y2 - y1) * (xmin - x1) / (x2 - x1);
         x = xmin;
       }
       
       // Replace the point outside the box
       if (codeOut == code1) {
         x1 = x;
         y1 = y;
         code1 = computeCode(x1, y1);
       } 
       else {
         x2 = x;
         y2 = y;
         code2 = computeCode(x2, y2);
       }
     }
   }
   
   if (accept) {
     // Compute the length of the clipped line segment
     double dx = x2 - x1;
     double dy = y2 - y1;
     double length = std::sqrt(dx * dx + dy * dy);
     return length;
   } 
   else {
     return 0.0; // No part inside
   }
 }
//#' Create individual use matrix
//#' 
//#' @param ch capture history matrix i x k x j of times of detection
//#' @param trapcells polygons of grid cells of each trap
//#' @param tracksdf dataframe of track locations, occ, and times
//#' 
//#' @return list of use matrices for each individual
Rcpp::NumericVector create_ind_use(Rcpp::NumericVector ch,
                          Rcpp::DataFrame traps,
                          double spacing,
                          Rcpp::DataFrame tracksdf,
                          Rcpp::NumericMatrix useall,
                          std::string scenario = "everything"
                                     ){
  Rcpp::NumericVector use(ch.size());
  Rcpp::IntegerVector newdims(3);
  Rcpp::IntegerVector olddims = ch.attr("dim");
  newdims[0] = olddims[0]; //i
  newdims[1] = olddims[2]; //j
  newdims[2] = olddims[1]; //k
  use.attr("dim") = newdims;
  //create tracklines
  Rcpp::List lines = create_line_list_C(tracksdf, scenario);
  
  //create trap grid bboxes
  Rcpp::NumericMatrix trap_cells = create_grid_bboxes_C(traps, spacing);

  auto indexify3D = [&](int i, 
                        int j,
                        int k,
                        int I,
                        int J) {
    int index = i + I * j + I * J * k;

    return index;
  };
  
  Rcpp::IntegerVector occ  = tracksdf["occ"];
  Rcpp::NumericVector time = tracksdf["time"];
  
  auto min_time = [&](int k) {
    double minval = NA_REAL;
    bool found = false;
    for (int i = 0; i < occ.size(); i++) {
      if (occ[i] == k) {
        if (!found || time[i] < minval) {
          minval = time[i];
          found = true;
        }
      }
    }
    return minval;
  };
  
  for (int i = 0; i < newdims[0]; ++i) {
    for (int k = 0; k < newdims[2]; ++k) {
      //create output 
      Rcpp::NumericVector lengths_in_traps(newdims[1]);
      
      //check if all ch[i,k,] are NA
      bool i_det_k = FALSE;
      double det_time_ik = NA_REAL;
      for (int j = 0; j < newdims[1]; ++j) {
        int indx3D = indexify3D(i, j, k, newdims[0], newdims[1]); //get vector index
        if(ch[indx3D] != NA_REAL) {
          i_det_k = TRUE;
          det_time_ik = min_time(k) + ch[indx3D];
          Rcpp::NumericMatrix linek = lines[k];
          Rcpp::NumericVector linetimes = linek.column(4);
          
          for (int step = 0; step < linek.nrow(); ++step) {
            if (linetimes[step] <= det_time_ik) {
              // now add this line length to total length for each trap
              lengths_in_traps[j] =  lengths_in_traps[j] + get_length_C(linek[step], trap_cells.row(j));
            }
          }
          break;
        }
      }
      if(i_det_k == FALSE){
        //if i not detected just useall[,k]
        lengths_in_traps = useall.column(k);
      }
    }//now just need to organize output
}
  
//  use <- mclapply(as.list(1:dim(ch)[1]), FUN = function(i){
//    apply(as.array(1:dim(ch)[2]), 1, function(k){
//      trackoccdf <- tracksdf[tracksdf$occ == k,]
//      if(!all(is.na(ch[i,k,]))){
//#if i detected k, discard survey pts after detection in occasion k
//        dettime <- min(trackoccdf$time) + ch[i,k,which(!is.na(ch[i,k,]))]
//        trackoccdf <- trackoccdf[trackoccdf$time <= dettime,]
//#add up length of trackline in each grid cell
//        useik <- lengths_in_grid(create_line_spatlines(trackoccdf), k, trap_cells)
//      } else {#if i wasn't detected in k, all traps used full amount, so skip subsetting
//        useik <- useall[,k]
//      }
//#(if there are covariates, perhaps could create matrix for them here, but for now ignore)
//    })
//  }, mc.cores = 3)
//  return(use)
return use;
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
  double haz  = g0; //log(1 - g0) * -1 *  1/timeincr ; //
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
    Rcpp::NumericMatrix mesh, //first column x, second y
    bool linear = false
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
    double mesharea;
    if(linear){
      mesharea = (meshxsorted(2) - meshxsorted(1))/1000; //km     
    } else {
      mesharea = ((meshxsorted(2) - meshxsorted(1)) * (meshysorted(2) - meshysorted(1)))/1000000; //km^2
    }
    clock.tock("setup");
    //begin for loops for lambdan calculation
    clock.tick("lambdan");
    Rcpp::NumericMatrix notseen_mk(meshx.length(), occs.size());
    Rcpp::NumericVector Dx_pdotxs(meshx.length());
    for(int m = 0; m < meshx.length(); m++){
      double Dx = D_mesh(m);
      //Rcpp::NumericVector notseen_eachocc((occs.size()));
      for(int occk = 0; occk < occs.size(); occk++){
        Rcpp::NumericVector hu_eachtrap((traps.size()));
        for(int trapj = 0; trapj < traps.size(); trapj++){
          if(usage(trapj, occk) == 0){
            hu_eachtrap(trapj) = 0; //if the trap isn't used, can't be seen at it
          } else {
            double thisdist = distmat(trapj, m);  //note this can't be recycled below except for individuals never detected. detected inds will have different induse
            hu_eachtrap(trapj) = hazdist_cpp(lambda0, sigma, thisdist, timeincr) * (usage(trapj, occk)/timeincr); //survival
          }
        }
        notseen_mk(m,occk) = exp(-sum(hu_eachtrap)); 
      }
      double notseen_alloccs = product(notseen_mk.row(m));
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
          double sumtoj_ind_ijk;
          double hu_ind_ijk;
          double sumtoj_ind = 0;
          for(int trapj = 0; trapj < traps.length(); trapj++){//could limit this to traps used in the occasion
            double captik = capthist(i, occk, trapj);
            //hazard for ind at trap
            double hu_ind_j = hazdist_cpp(lambda0, sigma, distmat(trapj, x), timeincr) * (indusage(i, trapj, occk)/timeincr);//hazard times individual usage for each trap. 
            
            if(captik == 1){ // this could be within the above else (only captures if used)?
              ikcaught = TRUE; //assign if i is caught
              //if(indusage(i,trapj,occk) <= timeincr){
                //do what i've been doing
                hu_ind_ijk = hu_ind_j; //hazard times use at trap/time of capture
                sumtoj_ind_ijk = sumtoj_ind; //sum of hazards up to trap before capture
              //} else { //else survive up to current trap - hazdenom and don't survive an interval of hazdenom
              //  hu_ind_ijk = hazdist_cpp(lambda0, sigma, distmat(trapj, x), timeincr) * timeincr; //hu for the last timeincr in this trap
              //  sumtoj_ind_ijk = sumtoj_ind + (hazdist_cpp(lambda0, sigma, distmat(trapj, x), timeincr) * ((indusage(i, trapj, occk) - timeincr)/timeincr)); //survive up to last time increment in this trap
            //  }
              break; //don't need to keep summing hazard after detection
            }
            //add cumulative hazard for ind at trap if not detected
            sumtoj_ind += hu_ind_j;
          }
          double prob_notseenk = notseen_mk(x,occk);
          if(ikcaught){
            probcapthist_eachocc(occk) = exp(-sumtoj_ind_ijk) * (1 - exp(-hu_ind_ijk));
          }else{
            probcapthist_eachocc(occk) = prob_notseenk; //survived all traps
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
    Rcpp::NumericMatrix mesh, //first column x, second y
    bool linear = false
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
    double mesharea;
    if(linear){
      mesharea = (meshxsorted(2) - meshxsorted(1))/1000; //km
    } else {
      mesharea = ((meshxsorted(2) - meshxsorted(1)) * (meshysorted(2) - meshysorted(1)))/1000000; //km^2 instead of ha
    }
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
          Rcpp::NumericVector hu_js(traps.length());
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
           if(sum_hujs < 1e-16){//would be odd if individual detected if sumhujs was almost 0...
             probcapthist_eachocc(occk) = (1/traps.length())   * (1 - exp(-sum_hujs)); 
           } else {
             //prob that trap j made the detection given i was detected in k
             probcapthist_eachocc(occk) = exp(log(hu_js(trapijk))-log(sum_hujs))   * (1 - exp(-sum_hujs)); 
           }
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
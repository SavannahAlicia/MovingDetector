#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppClock)]]
//[[Rcpp::plugins(cpp11)]]
#include <RcppClock.h>
#include <algorithm>
#include <string>
using namespace Rcpp;
//----------------- Supporting functions ---------------------------------------
// [[Rcpp::export]]
Rcpp::NumericVector which_is_not_na(Rcpp::NumericVector x) {
  NumericVector indices;
  for(int i = 0; i < x.size(); ++i) {
    if(Rcpp::NumericVector::is_na(x[i])) {}else{
      indices.push_back(i); 
    }
  }
  return indices;
}

// [[Rcpp::export]]
Rcpp::NumericMatrix calc_dist_matC(Rcpp::NumericMatrix t,//first pts (i) xy
                                   NumericMatrix h //second pts (j) xy
                                     ) {
  // Distance matrix has i rows and j columns
  NumericMatrix dist_mat(t.nrow(), h.nrow());
  for(int j = 0; j < h.nrow(); ++j) {
    for(int i = 0; i < t.nrow(); ++i) {
      double dx = t(i,0) - h(j,0);
      double dy = t(i, 1) - h(j, 1);
      dist_mat(i,j) = std::sqrt(dx * dx + dy * dy);
    }
  }
  return dist_mat; // [i, j] = distance(pt i, pt j)
}

// [[Rcpp::export]]
Rcpp::NumericVector seqC(int &first, int &last) {
  NumericVector y(abs(last - first) + 1);
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
  NumericMatrix result(n_cols, n_slices);
  
  // Copy the row from the cube to the NumericMatrix
  for (int col = 0; col < n_cols; ++col) {
    for (int slice = 0; slice < n_slices; ++slice){
      result(col, slice) = cube(row, col, slice);
    }
  }
  return result;
}

// [[Rcpp::export]]
Rcpp::NumericVector Sugar_colSums(const NumericMatrix x) {
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
  NumericVector newvec(n);
  for(int i = 0; i < n; i++){
    newvec(i) = log(x(i));
  }
  return newvec;
}

// [[Rcpp::export]]
double vec_min(Rcpp::NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::min_element(x.begin(), x.end());
  // we want the value so dereference 
  return *it;
}

// [[Rcpp::export]]
double vec_max(Rcpp::NumericVector x) {
  // Rcpp supports STL-style iterators
  NumericVector::iterator it = std::max_element(x.begin(), x.end());
  // we want the value so dereference 
  return *it;
}

// ----- Functions for calculating tracklength effort ------

// [[Rcpp::export]]

Rcpp::List create_line_list_C(Rcpp::DataFrame tracksdf, 
                             std::string scenario) {
  
  // Check that required columns exist
  CharacterVector colnames = tracksdf.names();
  
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
  CharacterVector effort;
  if (has_effort) {
    effort = tracksdf["effort"];
  }
  
  NumericVector x = tracksdf["x"];
  NumericVector y = tracksdf["y"];
  IntegerVector occ = tracksdf["occ"];
  NumericVector time = tracksdf["time"];
  
  // Find unique track IDs
  IntegerVector uniq_occ = sort_unique(occ);
  int n_occ = uniq_occ.size();
  
  // Initialize output
  List out_list(n_occ);
  CharacterVector out_names(n_occ);
  
  // Column names for numeric matrices
  CharacterVector mat_colnames(6);
  mat_colnames[0] = "x1";
  mat_colnames[1] = "y1";
  mat_colnames[2] = "x2";
  mat_colnames[3] = "y2";
  mat_colnames[4] = "time";
  mat_colnames[5] = "timeend";
  
  // Iterate over each unique track ID
  for (int k = 0; k < n_occ; ++k) {
    int track_id = uniq_occ[k];
    
    // Indices for this occasion from tracksdf (of points, not segments)
    LogicalVector mask = occ == track_id;
    IntegerVector idx = Rcpp::seq(0, occ.size() - 1); // indexes rows of entire tracksdf
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
    NumericMatrix seg_mat(seg_count, 6);
    
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
        seg_mat(seg_index, 4) = time[idx[i]]; //time of first pt
        seg_mat(seg_index, 5) = time[idx[i+1]];
        seg_index++;
      }
    }
    
    List dimnames(2);
    dimnames[0] = CharacterVector(seg_count); // empty row names
    dimnames[1] = mat_colnames;                     // column names
    Rf_setAttrib(seg_mat, R_DimNamesSymbol, dimnames);
    
    out_list[k] = seg_mat;
    out_names[k] = std::to_string(track_id);
  }
  
  out_list.attr("names") = out_names;
  return out_list;
}
 
 // [[Rcpp::export]]
 NumericMatrix create_grid_bboxes_C(Rcpp::DataFrame grid,
                                             double spacing) {
   
   // Extract coordinates
   NumericVector x = grid["x"];
   NumericVector y = grid["y"];
   
   int n = x.size();
   if (spacing <= 0) Rcpp::stop("spacing must be positive");
   
   double h = spacing / 2.0;
   
   // Prepare output matrix: columns = left, right, bottom, top
   NumericMatrix bboxes(n, 4);
   
   for (int i = 0; i < n; ++i) {
     bboxes(i,0) = x[i] - h; // left
     bboxes(i,1) = x[i] + h; // right
     bboxes(i,2) = y[i] - h; // bottom
     bboxes(i,3) = y[i] + h; // top
   }
   
   CharacterVector rownames(bboxes.nrow()); 
   CharacterVector colnames(4);
   colnames[0] = "left";
   colnames[1] = "right";
   colnames[2] = "bottom";
   colnames[3] = "top";
   
   List dimnames(2);
   dimnames[0] = rownames; // row names
   dimnames[1] = colnames; // column names
   
   Rf_setAttrib(bboxes, R_DimNamesSymbol, dimnames);
   
   return bboxes;
 }
 
 // [[Rcpp::export]]
 // Need to check error handling for vertical/horizontal lines
 double get_length_C(Rcpp::NumericVector line, 
                              NumericVector box) {
   
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

// [[Rcpp::export]]
arma::cube create_ind_use_C(arma::cube ch,
                                   DataFrame traps,
                                   double spacing,
                                   DataFrame tracksdf,
                                   std::string scenario
){
  int N = ch.n_rows; //first dim
  int J = ch.n_slices; //third dim of ch is traps
  int K = ch.n_cols; //second dim of ch is occ
  
  arma::cube use(N,J,K, arma::fill::zeros);

  //create tracklines
  List lines = create_line_list_C(tracksdf, scenario);
  
  //create trap grid bboxes
  NumericMatrix trap_cells = create_grid_bboxes_C(traps, spacing);
  
  IntegerVector occ  = tracksdf["occ"];
  IntegerVector uniq_occs = sort_unique(occ);
  
  for (int k = 0; k < K; ++k) {
    int track_id = uniq_occs[k];
    std::string track_id_str = std::to_string(track_id);

    // setup effort for if i not detected in k
    bool useallk_exists = false;
    NumericVector useallk(J, 0.0);
    
    // line segments for occasion k
    if (lines[track_id_str] == R_NilValue){
      Rcpp::warning("Line %d has no length, all effort 0 for this occasion.", k);
    } else {
     
     NumericMatrix linek = lines[track_id_str];
      if(!Rf_isMatrix(linek)){
        stop("Line " + track_id_str + " is not a matrix");
      }
      NumericVector linetimes = linek(_,4); //time of first pt in line
      
     // double min_timek = vec_min(linetimes);
      
      for (int i = 0; i < N; ++i) {
        
        //check if all ch[i,k,] are NA
        //bool i_det_k = false;
        double det_time_ik = NA_REAL;
        
        //this loop just searches ch for detection time
        for (int j = 0; j < J; ++j) {
          
          
          if (!Rcpp::NumericVector::is_na(ch(i,k,j))){ //trap where detection occurs
            //i_det_k = true;
            det_time_ik =  ch(i,k,j);
            
            break;
          }
        }
        
        // if i not detected in k, check if useallk calculated
        if(Rcpp::NumericVector::is_na(det_time_ik)){
          // i_det_k == false
          if(useallk_exists){
            
            for(int j = 0; j < J; ++j){
              
              use(i,j,k) = useallk[j];
            }
            
          } else { //useallk doesn't exist so calculated it
            
            //this loop adds up length of each segment within trap grid
            for (int j = 0; j < J; ++ j ){

              for (int step = 0; step < linek.nrow(); ++step) {
                // now add this line length to total length for each trap
                useallk[j] =  useallk[j] + get_length_C(linek(step,_), trap_cells.row(j));
              }
              
              use(i,j,k) = useallk[j];
            }
            
            // it does exist now
            useallk_exists = true;
            
          } 
        } else { //det_time_ik does exist
          
          //this loop adds up length of each segment within trap grid
          for (int j = 0; j < J; ++ j ){
            
            for (int step = 0; step < linek.nrow(); ++step) {
              if (linetimes[step] <= det_time_ik) {
                // now add this line length to total length for each trap
                use(i,j,k) =  use(i,j,k) + get_length_C(linek(step,_), trap_cells.row(j));
              }
            }
          }
          
        }
        
      }
      
    }
    
  }
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
                   int haz_denom) {
  double g0 = halfnormdet_cpp(lambda0, sigma, d);
  //haz_denom is the unit for hazard rate (could be time or distance (eg meters))
  double haz  = g0; //log(1 - g0) * -1 *  1/haz_denom ; //
  return haz;
}

//----------------- Simulating capture histories -------------------------------
// (still need these functions for linear meshes)

// [[Rcpp::export]]
NumericMatrix sim_pop_C( NumericVector D_mesh,
                               NumericMatrix mesh,
                               double meshspacing) {
  // N ~ Poisson(sum(D_mesh))
  double Dlambda = 0.0;
  for (int i = 0; i < D_mesh.size(); ++i) {
    Dlambda += D_mesh[i];
  }

  int N = R::rpois(Dlambda * meshspacing * meshspacing);  //square ms
  
  // sample mesh according to D_mesh
  NumericVector meshprobs = D_mesh/Rcpp::sum(D_mesh);
  IntegerVector Dmeshsample = seq(0, D_mesh.size()-1);
  IntegerVector popmeshAC = Rcpp::sample(Dmeshsample, N, true, meshprobs);
  
   // choose locations within each mesh cell uniformly
   NumericMatrix pop(N,2);
   for(int i = 0; i < N; ++i) {
     double x = mesh(popmeshAC[i],0);
     double y = mesh(popmeshAC[i],1);
     pop(i,0) = R::runif(x - (meshspacing/2), x + (meshspacing/2));
     pop(i,1) = R::runif(y - (meshspacing/2), y + (meshspacing/2));
   }
  return pop;
}

// [[Rcpp::export]]
NumericVector
sim_capthist_C(NumericMatrix traps,
                             DataFrame tracksdf,
                             double lambda0,
                             double sigma,
                             NumericVector D_mesh,
                             NumericMatrix mesh,
                             double meshspacing,
                             int hazdenom,
                             NumericMatrix pop,
                             Nullable<NumericMatrix> dist_dat_pop = R_NilValue,
                             bool report_probseenxk = false
){
  //dist_dat_pop needs to be distances to midpoints of increments if precalculated
  bool recalc_distdatpop = false;
  NumericMatrix dist_dat_pop_r;
  if(dist_dat_pop.isNotNull()){
    dist_dat_pop_r = NumericMatrix(dist_dat_pop);
  } else {
    dist_dat_pop_r = NumericMatrix(pop.nrow(), tracksdf.nrows());
    recalc_distdatpop = true;
  }

  if(D_mesh.size() == 0) Rcpp::stop("D_mesh cannot be empty or NULL");
  if(mesh.nrow() == 0) Rcpp::stop("mesh cannot be empty");
  if(pop.nrow() != dist_dat_pop_r.nrow()) Rcpp::stop("dist_dat_pop dimension does not match pop");
  if(tracksdf.nrows() != dist_dat_pop_r.ncol()) Rcpp::stop("dist_dat_pop ncol does not match tracksdf nrow");
  
  // Check that required columns exist
  CharacterVector colnames = tracksdf.names();
  
  // Check required columns
  if (std::find(colnames.begin(), colnames.end(), "x") == colnames.end() ||
      std::find(colnames.begin(), colnames.end(), "y") == colnames.end() ||
      std::find(colnames.begin(), colnames.end(), "time") == colnames.end() ||
      std::find(colnames.begin(), colnames.end(), "trapno") == colnames.end() ||
      std::find(colnames.begin(), colnames.end(), "occ") == colnames.end()) {
    Rcpp::stop("tracksdf must contain columns: x, y, occ, time, trapno");
  }
  
  NumericVector x = tracksdf["x"];
  NumericVector y = tracksdf["y"];
  NumericVector tracktimes = tracksdf["time"];
  IntegerVector occ = tracksdf["occ"];
  IntegerVector trapnos = tracksdf["trapno"];
  
  IntegerVector uniq_occs = sort_unique(Rcpp::as<Rcpp::IntegerVector>(tracksdf["occ"]));
  int nocc = uniq_occs.size();
  
  int N = pop.nrow();

  IntegerVector dims = {N, nocc, traps.nrow()};
  NumericVector ch(N * nocc * traps.nrow(), NA_REAL); //needs to filled with NA
  ch.attr("dim") = dims;

  NumericMatrix tracksxy(tracksdf.nrows(), 2);
  tracksxy = Rcpp::cbind(x,y);
 
  if(recalc_distdatpop){
    //need to calculate midpoints of tracksdf pts for each occasion for 
    // centered Riemann sum
    NumericMatrix tracksdfmidptxy(tracksdf.nrows(),2);
    for (int k = 0; k < nocc; ++k){
      
      int track_id = uniq_occs[k];
      
      // Indices for this occasion from tracksdf (of points, not segments)
      LogicalVector mask = occ == track_id;
      IntegerVector idx = Rcpp::seq(0, tracksdf.nrows() - 1); // indexes rows of entire tracksdf
      idx = idx[mask]; // subsetted to rows that are in this occasion
      int nsteps = idx.size();
      if(nsteps < 2){
        Rcpp::warning("Occ %d has less than 2 increments", k);
      }
      if(nsteps < 1) continue;
      
      // first midpoint is just the first point
      tracksdfmidptxy(idx[0], 0) = tracksxy(idx[0], 0);
      tracksdfmidptxy(idx[0], 1) = tracksxy(idx[0], 1);
      
      // Remaining midpoints are the averages of consecutive points
      for (int i = 1; i < nsteps; ++i) {
        tracksdfmidptxy(idx[i], 0) = 0.5 * (tracksxy(idx[i-1], 0) + tracksxy(idx[i], 0));
        tracksdfmidptxy(idx[i], 1) = 0.5 * (tracksxy(idx[i-1], 1) + tracksxy(idx[i], 1));
      }
    }
    
  //recalculate distances between hrcs and midpts for hazard later
    dist_dat_pop_r = calc_dist_matC(pop, tracksdfmidptxy);
  }
  
  //temporary output for testing
  NumericVector firsttprob(N * tracksdf.nrows());
  IntegerVector fpdims = {N, tracksdf.nrows()};
  firsttprob.attr("dim") = fpdims;

  //capthist loops
  NumericVector probseen(N * nocc);
  IntegerVector probseendims = {N, nocc};
  probseen.attr("dim") = probseendims;
   for (int k = 0; k < nocc; ++k){

      int track_id = uniq_occs[k];

      // Indices for this occasion from tracksdf (of points, not segments)
      LogicalVector mask = occ == track_id;
      IntegerVector idx = Rcpp::seq(0, tracksdf.nrows() - 1); // indexes rows of entire tracksdf
      idx = idx[mask]; // subsetted to rows that are in this occasion
      int nsteps = idx.size();
      if(nsteps < 2){
        Rcpp::warning("Occ %d has less than 2 increments", k);
      }
      // calculate 'increments' (step sizes along track, can be time or dist)
      //if hazard is per unit distance
      // calculate distance of line segments associated with tracksdf
      // (currently no way to subset line segments by effort label)
      // time denominated hazard not yet implemented
     // double begintime = vec_min(tracktimes[idx]);
      NumericVector increments(nsteps);

      for(int d = 1; d <  nsteps; d ++){ //first one is 0, so start at 1
        int x = idx[d];
        int xbefore = idx[d-1]; 
        increments[d] = std::sqrt((tracksxy(xbefore, 0) - tracksxy(x,0)) * (tracksxy(xbefore, 0) - tracksxy(x,0)) +
          (tracksxy(xbefore,1) - tracksxy(x,1)) * (tracksxy(xbefore,1) - tracksxy(x,1)));

      }

      for (int i = 0; i < N; ++i){
        NumericVector hus(nsteps);
        double cumhus = 0.0;
        NumericVector seenfirstatt(nsteps + 1); //last entry is not seen any t
        for(int e = 0; e < nsteps; e ++){
          int x = idx[e]; //index of row in tracksdf
          double d = dist_dat_pop_r(i,x);
          double haz = hazdist_cpp(lambda0, sigma, d, hazdenom);
          hus[e] = haz * increments[e];
          if(hus[e] < 0){
            Rcpp::stop("Error: negative hazard * effort");
          }
          double survive_until_tminus1 = exp(-cumhus);
          cumhus += hus[e];
         // double survive_until_t = exp(-cumhus);
          double survive_t_inc = exp(-hus[e]);
          seenfirstatt[e] = survive_until_tminus1 * (1.0 - survive_t_inc); //zero for e = 0
          firsttprob[i + idx[e] * N] = seenfirstatt[e];
        }
        double survk = exp(-cumhus);
        probseen[i + k * N] = 1 - survk;
        if(!report_probseenxk){
          seenfirstatt[nsteps] = survk;
          IntegerVector nstepsample = seq(0, nsteps); //length nsteps+1
          int stepdet = Rcpp::sample(nstepsample, 1, false, seenfirstatt)[0];
          if(stepdet < nsteps){ //if detection happened
            int trapdet = trapnos[idx[stepdet]] - 1; //trapno column is index based 1 in R
            if(trapdet < 0 || trapdet >= traps.nrow()) {
              Rcpp::stop("trapdet index out of range: %d (max allowed %d)", trapdet, traps.nrow()-1);
            }
            
            double timedet = tracktimes[idx[stepdet]];// - begintime;
            
            ch[i + N * k + N * nocc * trapdet] = timedet; //ch is ind x occ x trap, ixkxj 
          } 
          
        }
    }
  }
  if(report_probseenxk){
    return probseen;
  } else {
    return ch;
  }
}

//#---------------------Moving Detector multi-catch LLK (approximated with traps)
//-----------------Likelihood --------------------------------------------------

// [[Rcpp::export]]
double
  negloglikelihood_moving_cpp( //add log link 
    double lambda0,
    double sigma, 
    int haz_denom,
    NumericVector D_mesh,
    arma::cube capthist,
    arma::mat usage, //traps by occ, could be calculated from indusage
    arma::cube indusage, // ind by traps by occ
    NumericMatrix distmat, //traps x mesh 
    NumericMatrix mesh, //first column x, second y
    double mesharea
  ) {//specify objects
    Rcpp::Clock clock;
    clock.tick("wholeenchilada");
    clock.tick("setup");
    int captrap = capthist.n_slices;
    int capocc = capthist.n_cols;
    //calculate mesh area (note this is for a rectangular mesh grid, will need to
    //be recalculated for irregular shaped mesh)
    NumericVector meshx = mesh.column(0);
    NumericVector meshy = mesh.column(1);
    clock.tock("setup");
    //begin for loops for lambdan calculation
    clock.tick("lambdan");
    NumericMatrix notseen_mk_log(meshx.length(), capocc);
    NumericVector Dx_pdotxs(meshx.length());
    for(int m = 0; m < meshx.length(); m++){
      double Dx = D_mesh(m);
      //Rcpp::NumericVector notseen_eachocc((occs.size()));
      for(int occk = 0; occk < capocc; occk++){
        NumericVector hu_eachtrap(captrap);
        for(int trapj = 0; trapj < captrap; trapj++){
          if(usage(trapj, occk) == 0){
            hu_eachtrap(trapj) = 0; //if the trap isn't used, can't be seen at it
          } else {
            double thisdist = distmat(trapj, m);  //note this can't be recycled
            //below except for individuals never detected. detected inds will 
            //have different induse
            hu_eachtrap(trapj) = 
              hazdist_cpp(lambda0, sigma, thisdist, haz_denom) * 
              (usage(trapj, occk)/haz_denom); 
          }
        }
        notseen_mk_log(m,occk) = -sum(hu_eachtrap); // survival is exp(-sum(x))
      }
      double notseen_alloccs_log = sum(notseen_mk_log.row(m));
      double pdot = 1 - exp(notseen_alloccs_log);
      double Dx_pdotx = Dx * pdot;
      Dx_pdotxs(m) = Dx_pdotx;
    }
    double lambdan = sum(Dx_pdotxs) * mesharea;
    clock.tock("lambdan");
    //rest of likelihood
    clock.tick("loopllk");
    double n = capthist.n_rows;
    NumericVector integral_eachi_log(n);
    for(int i = 0; i < n; i++){
      NumericVector DKprod_eachx_log(meshx.length());
      for(int x = 0; x < meshx.length(); x++){
        NumericVector probcapthist_eachocc_log(capocc);
        for(int occk = 0; occk < capocc; occk++){
          bool ikcaught = false;
          double sumtoj_ind_ijk;
          double hu_ind_ijk;
          double sumtoj_ind = 0;
          arma::subview_row<double> induseik = indusage.slice(occk).row(i);
          for (arma::uword trapj = 0; trapj < captrap; trapj++) {
            //limit to used traps
            if (induseik[trapj] > 0) {
              //hazard times individual usage for each trap. 
              double hu_ind_j = hazdist_cpp(lambda0, sigma, distmat(trapj, x), haz_denom) * (induseik[trapj]/haz_denom);
              if(capthist(i, occk, trapj) == 1){
                ikcaught = true; //assign if i is caught
                hu_ind_ijk = hu_ind_j; //hazard times use at trap/time of capture
                //traps are not ordered by time, so I need to keep summing 
                //cumulative hazard after the detecting trap
              } else {
                //add cumulative hazard for ind at trap if not detected
                sumtoj_ind += hu_ind_j;
              }
            }
          }
          sumtoj_ind_ijk = sumtoj_ind; //sum of all hazards for ind i EXCEPT 
          // trap detected
          double prob_notseenk_log = notseen_mk_log(x,occk);
          if(ikcaught){
            probcapthist_eachocc_log(occk) = -sumtoj_ind_ijk + log(1 - exp(-hu_ind_ijk));
          } else {
            probcapthist_eachocc_log(occk) = prob_notseenk_log; //survived all traps
          }
        }
        double probcapthist_alloccs_log = sum(probcapthist_eachocc_log);
        DKprod_eachx_log(x) = log(D_mesh(x)) + probcapthist_alloccs_log;
      }
      double maxv = max(DKprod_eachx_log); //attempt to fix underflow
      double DKprod_sum_log = maxv  + log(sum(exp(DKprod_eachx_log - maxv)));

      integral_eachi_log(i) = DKprod_sum_log + log(mesharea);
    }
    clock.tock("loopllk");
    int n_int = std::round(n);
    NumericVector ns(n);
    int one = 1;
    ns = seqC(one, n_int);
    NumericVector logns(n);
    logns = log(ns);
    double lognfact = sum(logns);
    double sumlogint = sum(integral_eachi_log);
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
    int haz_denom,
    NumericVector D_mesh,
    arma::cube capthist,
    arma::mat usage, //traps by occ, could be calculated from indusage
    NumericMatrix distmat, //traps x mesh 
    NumericMatrix mesh, //first column x, second y
    double mesharea //area of single mesh cell
  ) {//specify objects
    Rcpp::Clock clock;
    clock.tick("wholeenchilada");
    clock.tick("setup");
    int one = 1;
    int captrap = capthist.n_slices;
    int capocc = capthist.n_cols;
    //calculate mesh area (note this is for a rectangular mesh grid, will need to
    //be recalculated for irregular shaped mesh)
    NumericVector meshx = mesh.column(0);
    NumericVector meshy = mesh.column(1);
    clock.tock("setup");
    //begin for loops for lambdan calculation
    clock.tick("lambdan");
    NumericMatrix notseen_mk_log(meshx.length(), capocc);
    NumericVector Dx_pdotxs(meshx.length());
    for(int m = 0; m < meshx.length(); m++){
      double Dx = D_mesh(m);
      //NumericVector notseen_eachocc_log((capocc));
      for(int occk = 0; occk < capocc; occk++){
        NumericVector hu_eachtrap((captrap));
        for(int trapj = 0; trapj < captrap; trapj++){
          if(usage(trapj, occk) == 0){
            hu_eachtrap(trapj) = 0; //if the trap isn't used, can't be seen at it
          } else {
            double thisdist = distmat(trapj, m);  
            hu_eachtrap(trapj) =
              hazdist_cpp(lambda0, sigma, thisdist, haz_denom) * 
              usage(trapj, occk)/haz_denom; 
          }
        }
        notseen_mk_log(m,occk) = -sum(hu_eachtrap); // survival is exp(-sum(x))
      }
      double notseen_alloccs_log = sum(notseen_mk_log.row(m));
      double pdot = 1 - exp(notseen_alloccs_log);
      double Dx_pdotx = Dx * pdot;
      Dx_pdotxs(m) = Dx_pdotx;
    }
    double lambdan = sum(Dx_pdotxs) * mesharea;
    clock.tock("lambdan");
    //rest of likelihood
    clock.tick("loopllk");
    double n = capthist.n_rows;
    NumericVector integral_eachi(n);
    for(int i = 0; i < n; i++){
      NumericVector DKprod_eachx_log(meshx.length());
      for(int x = 0; x < meshx.length(); x++){
        NumericVector probcapthist_eachocc(capocc);
        for(int occk = 0; occk < capocc; occk++){
          bool ikcaught = false;
          int trapijk;
          arma::vec hu_js(captrap, arma::fill::zeros);
          arma::vec usek = usage.col(occk);
          for(int trapj = 0; trapj < captrap; trapj++){
            //limit to traps used in occasion
            if(usek(trapj) >0 ){
            double captik = capthist(i, occk, trapj);
              hu_js(trapj) = hazdist_cpp(lambda0, sigma, distmat(trapj, x), haz_denom) * (usage(trapj, occk)/haz_denom);//hazard times usage for each trap. 
            if(captik == 1){ // this could be within the above else (only captures if used)?
              ikcaught = true; //assign if i is caught and which trap caught it
              trapijk = trapj;
            }
            }
          }
          double sum_hujs = sum(hu_js);
          if(ikcaught){
           //if(sum_hujs < 1e-16){//would be odd if individual detected if sumhujs was almost 0...
          //   probcapthist_eachocc(occk) = (1/captrap)   * (1 - exp(-sum_hujs)); 
           //} else {
             //prob that trap j made the detection given i was detected in k
             probcapthist_eachocc(occk) = exp(log(hu_js(trapijk))-log(sum_hujs))   * (1 - exp(-sum_hujs)); 
          // }
           
          } else {
            //prob i wasn't detected in k
            probcapthist_eachocc(occk) = exp(-sum_hujs) ; //survived all traps
          }
          //prevent underflow
          probcapthist_eachocc(occk) = std::max(probcapthist_eachocc(occk), 1e-16);
        }
        double probcapthist_alloccs_log = sum(log(probcapthist_eachocc));
        DKprod_eachx_log(x) = log(D_mesh(x)) + probcapthist_alloccs_log;
      }
      double maxv = max(DKprod_eachx_log); //prevent underflow
      double DKprod_sum = exp(maxv) * sum(exp(DKprod_eachx_log - maxv));
      
      integral_eachi(i) = DKprod_sum * mesharea;
      integral_eachi(i) = std::max(integral_eachi(i),  1e-16);
    }
    clock.tock("loopllk");
    int n_int = std::round(n);
    NumericVector ns(n);
    ns = seqC(one, n_int);
    NumericVector logns(n);
    logns = log(ns);
    double lognfact = sum(logns);
    NumericVector logint = logvec(integral_eachi);
    double sumlogint = sumC(logint);
    double out = -1 * (-lambdan - lognfact + sumlogint);
    clock.tock("wholeenchilada");
    clock.stop("approxstatllktimes");
    
    return(out);
  }





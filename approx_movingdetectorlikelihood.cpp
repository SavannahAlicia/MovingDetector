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
Rcpp::NumericMatrix calc_dist_matC(Rcpp::NumericMatrix t,//traps (j) xy
                                   NumericMatrix h //hrcs (i) xy
                                     ) {
  // Distance matrix has i rows and j columns
  NumericMatrix dist_mat(t.nrow(), h.nrow());
  for(int i = 0; i < h.nrow(); ++i) {
    for(int j = 0; j < t.nrow(); ++j) {
      double dx = t(j,0) - h(i,0);
      double dy = t(j, 1) - h(i, 1);
      dist_mat(j,i) = std::sqrt(dx * dx + dy * dy);
    }
  }
  return dist_mat; // t by h
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
                             std::string scenario = "everything") {
  
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
  CharacterVector mat_colnames(5);
  mat_colnames[0] = "x1";
  mat_colnames[1] = "y1";
  mat_colnames[2] = "x2";
  mat_colnames[3] = "y2";
  mat_colnames[4] = "time";
  
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
    NumericMatrix seg_mat(seg_count, 5);
    
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
        seg_mat(seg_index, 4) = time[idx[i]];
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
Rcpp::NumericVector create_ind_use_C(Rcpp::NumericVector ch,
                                   DataFrame traps,
                                   double spacing,
                                   DataFrame tracksdf,
                                   std::string scenario = "everything"
){
  NumericVector use(ch.size(), 0.0);
  IntegerVector newdims(3);
  IntegerVector olddims = ch.attr("dim");
  newdims[0] = olddims[0]; //i
  newdims[1] = olddims[2]; //j
  newdims[2] = olddims[1]; //k
  use.attr("dim") = newdims;
  //create tracklines
  List lines = create_line_list_C(tracksdf, scenario);
  
  //create trap grid bboxes
  NumericMatrix trap_cells = create_grid_bboxes_C(traps, spacing);
  
  auto indexify3D = [&](int i, 
                        int j,
                        int k,
                        int I,
                        int J) {
    int index = i + I * j + I * J * k;
    
    return index;
  };
  
  IntegerVector occ  = tracksdf["occ"];
  NumericVector time = tracksdf["time"];
  
  auto min_time = [&](int k) {
    double minval = NA_REAL;
    bool found = false;
    for (int i = 0; i < occ.size(); i++) {
      if (occ[i] == occ[k]) {
        if (!found || time[i] < minval) {
          minval = time[i];
          found = true;
        }
      }
    }
    return minval;
  };
  
  for (int k = 0; k < newdims[2]; ++k) {
    
    // setup effort for if i not detected in k
    bool useallk_exists = false;
    NumericVector useallk(newdims[1], 0.0);
    
    
    // line segments for occasion k
    if (lines[k] == R_NilValue){
      Rcpp::warning("Line %d has no length, all effort 0 for this occasion.", k);
    } else {
     
     NumericMatrix linek = lines[k];
      NumericVector linetimes = linek(_,4); 
      
      double min_timek = min_time(k);
      
      for (int i = 0; i < newdims[0]; ++i) {
        //create output 
        // NumericVector lengths_in_traps(newdims[1]);
        
        //check if all ch[i,k,] are NA
        //bool i_det_k = false;
        double det_time_ik = NA_REAL;
        
        //this loop just searches ch for detection time
        for (int j = 0; j < newdims[1]; ++j) {
          
          //get vector index for capture history
          int ix_ch = indexify3D(i, k, j, olddims[0], olddims[1]); 
          
          if (!Rcpp::NumericVector::is_na(ch[ix_ch])){ //trap where detection occurs
            //i_det_k = TRUE;
            det_time_ik = min_timek + ch[ix_ch];
            
            break;
          }
        }
        
        // if i not detected in k, check if useallk calculated
        if(Rcpp::NumericVector::is_na(det_time_ik)){
          // i_det_k == false
          if(useallk_exists){
            
            for(int j = 0; j < newdims[1]; ++j){
              
              int ix_use = indexify3D(i,j,k,newdims[0], newdims[1]);
              use[ix_use] = useallk[j];
            }
            
          } else { //useallk doesn't exist so calculated it
            
            //this loop adds up length of each segment within trap grid
            for (int j = 0; j < newdims[1]; ++ j ){
              
              int ix_use = indexify3D(i,j,k,newdims[0], newdims[1]);
              
              for (int step = 0; step < linek.nrow(); ++step) {
                // now add this line length to total length for each trap
                useallk[j] =  useallk[j] + get_length_C(linek(step,_), trap_cells.row(j));
              }
              
              use[ix_use] = useallk[j];
            }
            
            // it does exist now
            useallk_exists = TRUE;
            
          } 
        } else { //det_time_ik does exist
          
          //this loop adds up length of each segment within trap grid
          for (int j = 0; j < newdims[1]; ++ j ){
            
            int ix_use = indexify3D(i,j,k,newdims[0], newdims[1]);
            
            for (int step = 0; step < linek.nrow(); ++step) {
              if (linetimes[step] <= det_time_ik) {
                // now add this line length to total length for each trap
                use[ix_use] =  use[ix_use] + get_length_C(linek(step,_), trap_cells.row(j));
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
  double mesharea = meshspacing * meshspacing/ 1000000; //square kms
  int N = R::rpois(Dlambda * mesharea);
  
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
                             NumericMatrix dist_dat_pop,
                             bool report_probseenxk = false
){
  if(D_mesh.size() == 0) Rcpp::stop("D_mesh cannot be empty or NULL");
  if(mesh.nrow() == 0) Rcpp::stop("mesh cannot be empty");
  if(pop.nrow() != dist_dat_pop.nrow()) Rcpp::stop("dist_dat_pop dimension does not match pop");
  if(tracksdf.nrows() != dist_dat_pop.ncol()) Rcpp::stop("dist_dat_pop ncol does not match tracksdf nrow");
  
  // Check that required columns exist
  CharacterVector colnames = tracksdf.names();
  
  // Check required columns
  if (std::find(colnames.begin(), colnames.end(), "x") == colnames.end() ||
      std::find(colnames.begin(), colnames.end(), "y") == colnames.end() ||
      std::find(colnames.begin(), colnames.end(), "time") == colnames.end() ||
      std::find(colnames.begin(), colnames.end(), "trapno") == colnames.end() ||
      std::find(colnames.begin(), colnames.end(), "occ") == colnames.end()) {
    Rcpp::stop("tracksdf must contain columns: x, y, occ, time");
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
  tracksxy(_, 0) = x;
  tracksxy(_, 1) = y;

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
      double begintime = vec_min(tracktimes[idx]);
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
          double d = dist_dat_pop(i,x);
          double haz = hazdist_cpp(lambda0, sigma, d, hazdenom);
          hus[e] = haz * increments[e];
          if(hus[e] < 0){
            Rcpp::stop("Negative hazard * effort");
          }
          double survive_until_tminus1 = exp(-cumhus);
          cumhus += hus[e];
         // double survive_until_t = exp(-cumhus);
          double survive_t_inc = exp(-hus[e]);
          if(e == 0){
            seenfirstatt[e] = 0;
          } else{
            seenfirstatt[e] = survive_until_tminus1 * (1.0 - survive_t_inc);
          }
        }
        double integ = Rcpp::sum(hus);
        double survk = exp(-integ);
        probseen[i + k * N] = 1 - survk;
        if(!report_probseenxk){
          seenfirstatt[nsteps] = survk;
          IntegerVector nstepsample = seq(0, nsteps); //length nsteps+1
          int stepdet = Rcpp::sample(nstepsample, 1, false, seenfirstatt)[0];
          if(stepdet < nsteps){
            int trapdet = trapnos[idx[stepdet]] - 1; //trapno column is index based 1 in R
            if(trapdet < 0 || trapdet >= traps.nrow()) {
              Rcpp::stop("trapdet index out of range: %d (max allowed %d)", trapdet, traps.nrow()-1);
            }
            
            double timedet = tracktimes[idx[stepdet]] - begintime;
            
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
    NumericMatrix usage, //traps by occ, could be calculated from indusage
    arma::cube indusage, // ind by traps by occ
    NumericMatrix distmat, //traps x mesh 
    NumericMatrix mesh, //first column x, second y
    bool linear = false
  ) {//specify objects
    Rcpp::Clock clock;
    clock.tick("wholeenchilada");
    clock.tick("setup");
    int one = 1;
    int captrap = capthist.n_slices;
    NumericVector traps = seqC(one, captrap); //needs to still just be a list of trap ID number
    int capocc = capthist.n_cols;
    NumericVector occs = seqC(one, capocc); //seq(1:n_occasions) but 0 base
    //calculate mesh area (note this is for a rectangular mesh grid, will need to
    //be recalculated for irregular shaped mesh)
    NumericVector meshx = mesh.column(0);
    NumericVector meshy = mesh.column(1);
    NumericVector meshxsorted = Rcpp::sort_unique(meshx);
    NumericVector meshysorted = Rcpp::sort_unique(meshy);
    double mesharea;
    if(linear){
      mesharea = (meshxsorted(2) - meshxsorted(1))/1000; //km     
    } else {
      mesharea = ((meshxsorted(2) - meshxsorted(1)) * (meshysorted(2) - meshysorted(1)))/1000000; //km^2
    }
    clock.tock("setup");
    //begin for loops for lambdan calculation
    clock.tick("lambdan");
    NumericMatrix notseen_mk(meshx.length(), occs.size());
    NumericVector Dx_pdotxs(meshx.length());
    for(int m = 0; m < meshx.length(); m++){
      double Dx = D_mesh(m);
      //Rcpp::NumericVector notseen_eachocc((occs.size()));
      for(int occk = 0; occk < occs.size(); occk++){
        NumericVector hu_eachtrap((traps.size()));
        for(int trapj = 0; trapj < traps.size(); trapj++){
          if(usage(trapj, occk) == 0){
            hu_eachtrap(trapj) = 0; //if the trap isn't used, can't be seen at it
          } else {
            double thisdist = distmat(trapj, m);  //note this can't be recycled below except for individuals never detected. detected inds will have different induse
            hu_eachtrap(trapj) = hazdist_cpp(lambda0, sigma, thisdist, haz_denom) * (usage(trapj, occk)/haz_denom); //survival
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
    NumericVector integral_eachi(n);
    for(int i = 0; i < n; i++){
      NumericVector DKprod_eachx(meshx.length());
      for(int x = 0; x < meshx.length(); x++){
        NumericVector probcapthist_eachocc(occs.length());
        for(int occk = 0; occk < occs.length(); occk++){
          bool ikcaught = false;
          double sumtoj_ind_ijk;
          double hu_ind_ijk;
          double sumtoj_ind = 0;
          for(int trapj = 0; trapj < traps.length(); trapj++){//could limit this to traps used in the occasion
            double captik = capthist(i, occk, trapj);
            //hazard for ind at trap
            double hu_ind_j = hazdist_cpp(lambda0, sigma, distmat(trapj, x), haz_denom) * (indusage(i, trapj, occk)/haz_denom);//hazard times individual usage for each trap. 
            
            if(captik == 1){ // this could be within the above else (only captures if used)?
              ikcaught = TRUE; //assign if i is caught
              //if(indusage(i,trapj,occk) <= haz_denom){
                //do what i've been doing
                hu_ind_ijk = hu_ind_j; //hazard times use at trap/time of capture
                sumtoj_ind_ijk = sumtoj_ind; //sum of hazards up to trap before capture
              //} else { //else survive up to current trap - hazdenom and don't survive an interval of hazdenom
              //  hu_ind_ijk = hazdist_cpp(lambda0, sigma, distmat(trapj, x), haz_denom) * haz_denom; //hu for the last haz_denom in this trap
              //  sumtoj_ind_ijk = sumtoj_ind + (hazdist_cpp(lambda0, sigma, distmat(trapj, x), haz_denom) * ((indusage(i, trapj, occk) - haz_denom)/haz_denom)); //survive up to last time increment in this trap
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
    NumericVector ns(n);
    ns = seqC(one, n_int);
    NumericVector logns(n);
    logns = log(ns);
    double lognfact = sum(logns);
    NumericVector logint = logvec(integral_eachi);
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
    int haz_denom,
    NumericVector D_mesh,
    arma::cube capthist,
    NumericMatrix usage, //traps by occ, could be calculated from indusage
    NumericMatrix distmat, //traps x mesh 
    NumericMatrix mesh, //first column x, second y
    bool linear = false
  ) {//specify objects
    Rcpp::Clock clock;
    clock.tick("wholeenchilada");
    clock.tick("setup");
    int one = 1;
    int captrap = capthist.n_slices;
    NumericVector traps = seqC(one, captrap); //needs to still just be a list of trap ID number
    int capocc = capthist.n_cols;
    NumericVector occs = seqC(one, capocc); //seq(1:n_occasions) but 0 base
    //calculate mesh area (note this is for a rectangular mesh grid, will need to
    //be recalculated for irregular shaped mesh)
    NumericVector meshx = mesh.column(0);
    NumericVector meshy = mesh.column(1);
    NumericVector meshxsorted = Rcpp::sort_unique(meshx);
    NumericVector meshysorted = Rcpp::sort_unique(meshy);
    double mesharea;
    if(linear){
      mesharea = (meshxsorted(2) - meshxsorted(1))/1000; //km
    } else {
      mesharea = ((meshxsorted(2) - meshxsorted(1)) * (meshysorted(2) - meshysorted(1)))/1000000; //km^2 instead of ha
    }
    clock.tock("setup");
    //begin for loops for lambdan calculation
    clock.tick("lambdan");
    NumericVector Dx_pdotxs(meshx.length());
    for(int m = 0; m < meshx.length(); m++){
      double Dx = D_mesh(m);
      NumericVector notseen_eachocc((occs.size()));
      for(int occk = 0; occk < occs.size(); occk++){
        NumericVector notseen_eachtrap((traps.size()));
        for(int trapj = 0; trapj < traps.size(); trapj++){
          if(usage(trapj, occk) == 0){
            notseen_eachtrap(trapj) = 1; //if the trap isn't used, can't be seen at it
          } else {
            double thisdist = distmat(trapj, m);  //note this can't be recycled below except for individuals never detected. detected inds will have different induse
            notseen_eachtrap(trapj) = exp(-hazdist_cpp(lambda0, sigma, thisdist, haz_denom) * (usage(trapj, occk)/haz_denom)); //survival
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
    NumericVector integral_eachi(n);
    for(int i = 0; i < n; i++){
      NumericVector DKprod_eachx(meshx.length());
      for(int x = 0; x < meshx.length(); x++){
        NumericVector probcapthist_eachocc(occs.length());
        for(int occk = 0; occk < occs.length(); occk++){
          bool ikcaught = false;
          int trapijk;
          NumericVector hu_js(traps.length());
          for(int trapj = 0; trapj < traps.length(); trapj++){//could limit this to traps used in the occasion
            double captik = capthist(i, occk, trapj);
              hu_js(trapj) = hazdist_cpp(lambda0, sigma, distmat(trapj, x), haz_denom) * (usage(trapj, occk)/haz_denom);//hazard times usage for each trap. 
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





#simulate data from fitted CHS

#-------------- Read in objects from CHS ----------

Xmat <- readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/Xmats.Rds")[2][[1]]
fit <- readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/m_moves(x, y, k = 6).Rds")
#use fitted values from moving detector model 
lambda0 <- exp(fit$movdet_est$value[7])
sigma <- exp(fit$movdet_est$value[8])
D_mesh <- exp(Xmat %*% fit$movdet_est$value[1:6])

traps <- readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHSinput/traps.Rds")
tracksdf <- readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHSinput/tracksdf.Rds")
mesh <- readRDS("data/onison/all_occasions/meshscr_NSbuff_2000.Rds")
meshspacing <- 2000
hazdenom <- 1
dist_trapmesh <- readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHSinput/distmat_trapmesh.Rds")

#-------------- Read in functions --------------------

Rcpp::sourceCpp("~/Documents/UniStAndrews/MovingDetector/approx_movingdetectorlikelihood.cpp")

simulate_popandcapthist <- function(traps,
                                    tracksdf, 
                                    lambda0,
                                    sigma,
                                    D_mesh,
                                    mesh,
                                    meshspacing,
                                    hazdenom){
  start.time.sim <- Sys.time()
  #capthist dim(inds, traps)
  pop <- sim_pop_C(D_mesh, 
                   as.matrix(mesh), 
                   meshspacing)
  if(is.null(dim(pop))){
    stop("No population simulated, try larger density.")
  } else if(nrow(pop) < 2){
    warning("Only one individual simulated.")
  } else if(nrow(pop) <20) {
    warning("Less than 20 individuals simualted.")
  }
  dist_dat_pop <- calc_dist_matC(pop, 
                                 as.matrix(tracksdf[,c("x","y")]))
  
  capthist_full <- sim_capthist_C(as.matrix(traps),
                                  tracksdf, 
                                  lambda0,
                                  sigma,
                                  D_mesh,
                                  as.matrix(mesh),
                                  meshspacing,
                                  hazdenom,
                                  pop,
                                  dist_dat_pop,
                                  report_probseenxk = F)
  
  #standard scr likelihood
  nocc <- length(unique(tracksdf$occ))
  trapspacing <- sort(unique(traps$x))[2]- sort(unique(traps$x))[1]
  
  if(length(which(apply((!is.na(capthist_full)), 1, sum)>0)) == 0){
    warning("Empty capture history.")
    capthist <- array(NA, dim = c(1, nocc, nrow(traps)))
  } else {
    capthist <- capthist_full[which(apply((!is.na(capthist_full)), 1, sum)>0),,]
  }
  
  #use
  induse <- create_ind_use_C(capthist, as.matrix(traps), trapspacing, tracksdf,
                             scenario = "everything")
  
  #convert capthist to 1s and 0s
  capthist[is.na(capthist)] <- 0
  capthist[capthist!=0] <- 1
  
  fit.time.sim <- difftime(Sys.time(), start.time.sim, units = "secs")
  out_ls <- list(capthist = capthist,
                 induse = induse,
                 fit.time.sim = fit.time.sim)
}

fit_capthist <- function(dist_trapmesh,
                         useall,
                         lambda0, 
                         sigma, 
                         D_mesh, 
                         beta1, 
                         beta2,
                         beta3,
                         hazdenom, 
                         mesh, 
                         capthistout
){

  
  capthist = capthistout$capthist
  induse = capthistout$induse
  fit.time.sim = capthistout$fit.time.sim
  
  #in case mesh is df
  mesh_mat <- as.matrix(mesh)
  
    start0 <- c(
      log(lambda0),
      log(sigma), fit$movdet_est$value[1:6])
    scaling_factors <- 10^round(log10(abs(start0)))
    start <- start0/scaling_factors
    
    #quadratic density function
    stat_nll <- function(v_scaled){
      v <- v_scaled * scaling_factors 
      lambda0_ <- exp(v[1])
      sigma_ <- exp(v[2])
      D_mesh_ <- exp(Xmat %*% v[3:8])
      out <- negloglikelihood_stationary_cpp(lambda0_, sigma_,
                                             hazdenom, D_mesh_, 
                                             capthist, useall,
                                             dist_trapmesh, mesh_mat)
      return(out)
    }
    #moving detector likelihood
    nll <- function(v_scaled){
      v <- v_scaled * scaling_factors 
      lambda0_ <- exp(v[1])#invlogit(v[1])
      sigma_ <- exp(v[2])
      D_mesh_ <- exp(Xmat %*% v[3:8])
      out <- negloglikelihood_moving_cpp(lambda0_, sigma_,  
                                         hazdenom, D_mesh_,
                                         capthist, useall,
                                         induse, dist_trapmesh, mesh_mat)
      return(out)
    }
  
  
  start.time.sd <- Sys.time()
  fit_sd <- optim(par = start,
                  fn = stat_nll,
                  hessian = F, method = "Nelder-Mead") #NM is best at this likelihood, even though slower
  fit_sd$hessian <- numDeriv::hessian(stat_nll, x = fit_sd$par,method = "Richardson",
                                      method.args = list(eps = 1e-5, d = 1e-3, r = 3))
  fit.time.sd <- difftime(Sys.time(), start.time.sd, units = "secs") #includes hessian
  
  start.time.md <- Sys.time()
  fit_md <- optim(par = start,
                  fn = nll,
                  hessian = F, method = "Nelder-Mead")
  fit_md$hessian <- numDeriv::hessian(nll, x = fit_md$par,method = "Richardson",
                                      method.args = list(eps = 1e-6, d = 1e-4, r = 4))
  fit.time.md <- difftime(Sys.time(), start.time.md, units = "secs")
  
    outnames <- c("lambda0", "sigma", 
                  "D1", "D2", "D3", "D4", "D5", "D6")
  
  assemble_CIs <- function(fit){
    fisher_info <- MASS::ginv(fit$hessian)
    prop_sigma <- sqrt(diag(fisher_info))
    prop_sigma <- diag(prop_sigma)
    upper <- fit$par+1.96*prop_sigma
    lower <- fit$par-1.96*prop_sigma
    interval <- data.frame(name = outnames,
                           value = fit$par * scaling_factors, 
                           upper = diag(upper) * scaling_factors, 
                           lower = diag(lower) * scaling_factors
    )
    
    return(interval)
  } 
  
  out <- list(statdet_est = assemble_CIs(fit_sd), 
              movdet_est = assemble_CIs(fit_md),
              statdet_time = fit.time.sd,
              movdet_time = fit.time.md,
              sim_time = fit.time.sim,
              n = dim(capthist)[1])
  return(out)
}

sim_fit <- function(traps,
                    tracksdf, 
                    mesh, 
                    meshspacing,
                    dist_trapmesh,
                    useall,
                    lambda0, 
                    sigma, 
                    D_mesh, 
                    beta1, 
                    beta2,
                    beta3,
                    hazdenom){
  capthistout <- simulate_popandcapthist(traps,
                                         tracksdf, 
                                         lambda0,
                                         sigma,
                                         D_mesh,
                                         mesh,
                                         meshspacing,
                                         hazdenom)
  out <- fit_capthist(dist_trapmesh,
                      useall,
                      lambda0, 
                      sigma, 
                      D_mesh, 
                      beta1, 
                      beta2,
                      beta3,
                      hazdenom, 
                      mesh, 
                      capthistout)
  return(out)
}


#--------------- Do some simulations ----------------
nsims = 50 
start.time.all_q <- Sys.time()
all_sim_fits_q <- mclapply(X = as.list(1:nsims),
                           FUN = function(sim){
                             return(sim_fit(traps,
                                            tracksdf, 
                                            mesh, 
                                            meshspacing,
                                            dist_trapmesh,
                                            useall,
                                            lambda0, sigma, D_mesh, 
                                            beta1, beta2, beta3,
                                            hazdenom, 
                                            Dmod = "~x^2"))
                           },
                           mc.cores = 6
)
tot.time.all_q <- difftime(Sys.time(), start.time.all_q, units = "secs")
print(tot.time.all_q)
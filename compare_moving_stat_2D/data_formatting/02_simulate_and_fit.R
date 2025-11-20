## 2 dimensional space
###-----------simulate data for moving and stationary detectors-----------------
###-------and fit model with stationary and moving detector---------------------

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
  capthist <- capthist_full[which(apply((!is.na(capthist_full)), 1, sum)>0),,]
  
  #standard scr likelihood
  nocc <- length(unique(tracksdf$occ))
  trapspacing <- sort(unique(traps$x))[2]- sort(unique(traps$x))[1]
  
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
                    Dmod = "~1",
                    capthistout){

  capthist = capthistout$capthist
  induse = capthistout$induse
  
  #in case mesh is df
  mesh_mat <- as.matrix(mesh)
  
  #grid with stationary detector likelihood
  if (Dmod == "~1"){
    start0 <- c( 
      log(lambda0),
      log(sigma), log(D_mesh[1]))
    scaling_factors <- 10^round(log10(abs(start0)))
    start <- start0/scaling_factors
    
    stat_nll <- function(v_scaled ){
      v <- v_scaled * scaling_factors 
      lambda0_ <- exp(v[1])
      sigma_ <- exp(v[2])
      D_mesh_ <- rep(exp(v[3]), nrow(mesh_mat))
      out <- negloglikelihood_stationary_cpp(lambda0_, sigma_,
                                             hazdenom, D_mesh_, 
                                             capthist, useall,
                                             dist_trapmesh, mesh_mat)
      return(out)
    }
    #moving detector likelihood
    nll <- function(v_scaled ){
      v <- v_scaled * scaling_factors 
      lambda0_ <- exp(v[1])#invlogit(v[1]) #logit link?
      sigma_ <- exp(v[2])
      D_mesh_ <- rep(exp(v[3]), nrow(mesh))
      out <- negloglikelihood_moving_cpp(lambda0_, sigma_,  
                                         hazdenom, D_mesh_,
                                         capthist, useall,
                                         induse, dist_trapmesh, mesh_mat)
      return(out)
    }

    
  }else if(Dmod == "~x^2"){  
    start0 <- c(
      log(lambda0),
      log(sigma), beta1, beta2, beta3)
    scaling_factors <- 10^round(log10(abs(start0)))
    start <- start0/scaling_factors
    
    #quadratic density function
    stat_nll <- function(v_scaled ){
      v <- v_scaled * scaling_factors 
      lambda0_ <- exp(v[1])#invlogit(v[1])
      sigma_ <- exp(v[2])
      D_mesh_ <- exp(v[3]*(mesh_mat[,1] + v[4])^2 + v[5])
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
      D_mesh_ <- exp(v[3]*(mesh_mat[,1] + v[4])^2 + v[5])#exp(beta1*(mesh$x + beta2)^2 + beta3)
      out <- negloglikelihood_moving_cpp(lambda0_, sigma_,  
                                         hazdenom, D_mesh_,
                                         capthist, useall,
                                         induse, dist_trapmesh, mesh_mat)
      return(out)
    }
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
  
  if (Dmod == "~1"){
    outnames <- c("lambda0", "sigma", "D")
  } else if(Dmod == "~x^2"){
    outnames <- c("lambda0", "sigma", "beta1", "beta2")
  }
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

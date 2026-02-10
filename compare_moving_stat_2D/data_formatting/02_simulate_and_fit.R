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
  if(is.null(dim(pop))){
    stop("No population simulated, try larger density.")
  } else if(nrow(pop) < 2){
    warning("Only one individual simulated.")
  } else if(nrow(pop) <20) {
    warning("Less than 20 individuals simualted.")
  }
  
  
  capthist_full <- sim_capthist_C(as.matrix(traps),
                                  tracksdf, 
                                  lambda0,
                                  sigma,
                                  D_mesh,
                                  as.matrix(mesh),
                                  meshspacing,
                                  hazdenom,
                                  pop,
                                  dist_dat_pop = NULL,
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
                    N,
                    hazdenom, 
                    mesh, 
                    capthistout,
                    Dmod 
                    ){
  if(!Dmod %in% c("~1", "~x^2")){
    stop("Dmod must be specified ~1 or ~x^2")
  }

  capthist = capthistout$capthist
  induse = capthistout$induse
  fit.time.sim = capthistout$fit.time.sim
  
  #in case mesh is df
  mesh_mat <- as.matrix(mesh)
  
  #grid with stationary detector likelihood
  if (Dmod == "~1"){
    start0 <- c( 
      log(lambda0),
      log(sigma), log(D_mesh[1]))
    scaling_factors <-  rep(1, length(start0)) #10^round(log10(pmax(abs(start0), 1e-8)))
    start <- start0/scaling_factors
    
    stat_nll <- function(v_scaled ){
      v <- v_scaled * scaling_factors 
      lambda0_ <- exp(v[1])
      if (!is.finite(lambda0_)) return(1e12)
      
      sigma_ <- exp(v[2])
      if (!is.finite(sigma_)) return(1e12)
      
      D_mesh_ <- rep(exp(v[3]), nrow(mesh_mat))
      if (any(!is.finite(D_mesh_))) return(1e12)
      
      out <- negloglikelihood_stationary_cpp(lambda0_, sigma_,
                                             hazdenom, D_mesh_, 
                                             capthist, useall,
                                             dist_trapmesh, mesh_mat,
                                             mesharea = meshspacing^2)
      
      return(out)
    }
    #moving detector likelihood
    nll <- function(v_scaled ){
      v <- v_scaled * scaling_factors 
      lambda0_ <- exp(v[1])
      if (!is.finite(lambda0_)) return(1e12)
      
      sigma_ <- exp(v[2])
      if (!is.finite(sigma_)) return(1e12)
      
      D_mesh_ <- rep(exp(v[3]), nrow(mesh))
      if (any(!is.finite(D_mesh_))) return(1e12)
      
      out <- negloglikelihood_moving_cpp(lambda0_, sigma_,  
                                         hazdenom, D_mesh_,
                                         capthist, useall,
                                         induse, dist_trapmesh, mesh_mat,
                                         mesharea = meshspacing^2)
      return(out)
    }

    
  }else if(Dmod == "~x^2"){  
    start0 <- c(
      log(lambda0),
      log(sigma),
      #beta1 * x_sd^2, 
      beta2, 
      log(N))
    scaling_factors <- rep(1, length(start0)) #10^round(log10(abs(start0)))
    start <- start0/scaling_factors
    
    #quadratic density function
    stat_nll <- function(v_scaled){
      v <- v_scaled * scaling_factors 
      lambda0_ <- exp(v[1])
      if (!is.finite(lambda0_)) return(1e12)
      
      sigma_ <- exp(v[2])
      if (!is.finite(sigma_)) return(1e12)
      
      eta <- #v[3]*
        fixed_beta1 * (mesh_mat[,1] + v[3])^2
      if (any(!is.finite(eta))) return(1e12)
      exp_eta <- exp(eta)
      if (any(!is.finite(exp_eta))) return(1e12)
      
      Z   <- sum(exp_eta) * meshspacing^2
      if (!is.finite(Z) || Z <= 0) return(1e12)
      
      D_mesh_  <- exp(v[4]) * exp_eta / Z
      if (any(!is.finite(D_mesh_))) return(1e12)
      
      out <- negloglikelihood_stationary_cpp(lambda0_, sigma_,
                                             hazdenom, D_mesh_, 
                                             capthist, useall,
                                             dist_trapmesh, mesh_mat,
                                             mesharea = meshspacing^2)
      return(out)
    }
    #moving detector likelihood
    nll <- function(v_scaled){
      v <- v_scaled * scaling_factors 
      lambda0_ <- exp(v[1])
      if (!is.finite(lambda0_)) return(1e12)
      
      sigma_ <- exp(v[2])
      if (!is.finite(sigma_)) return(1e12)
      
      eta <- #v[3]*
        fixed_beta1 * (mesh_mat[,1] + v[3])^2
      if (any(!is.finite(eta))) return(1e12)
      exp_eta <- exp(eta)
      if (any(!is.finite(exp_eta))) return(1e12)
      
      Z   <- sum(exp_eta) * meshspacing^2
      if (!is.finite(Z) || Z <= 0) return(1e12)
      
      D_mesh_ <- exp(v[4]) * exp_eta / Z
      if (any(!is.finite(D_mesh_))) return(1e12)
      
      out <- negloglikelihood_moving_cpp(lambda0_, sigma_,  
                                         hazdenom, D_mesh_,
                                         capthist, useall,
                                         induse, dist_trapmesh, mesh_mat,
                                         mesharea = meshspacing^2)
      return(out)
    }
  }  
  
  start.time.sd <- Sys.time()
  fit_sd <- optim(par = start,
                  fn = stat_nll,
                  hessian = F, method = "Nelder-Mead") #NM is best at this likelihood, even though slower
  fit_sd_con = fit_sd$convergence
  if(fit_sd_con == 0){
    fit_sd$hessian <- numDeriv::hessian(stat_nll, x = fit_sd$par,method = "Richardson",
                                        method.args = list(eps = 1e-5, d = 1e-3, r = 3))
    fit.time.sd <- difftime(Sys.time(), start.time.sd, units = "secs") #includes hessian
    
  } else { #if model did not converge, don't bother hessian 
    fit_sd$hessian = NA
    fit.time.sd = NA
  }
  
  start.time.md <- Sys.time()
  fit_md <- optim(par = start,
                  fn = nll,
                  hessian = F, method = "Nelder-Mead")
  fit_md_con = fit_md$convergence
  if(fit_md_con == 0){
    fit_md$hessian <- numDeriv::hessian(nll, x = fit_md$par,method = "Richardson",
                                        method.args = list(eps = 1e-6, d = 1e-4, r = 4))
    
    fit.time.md <- difftime(Sys.time(), start.time.md, units = "secs")
    
  } else {
    fit_md$hessian <- NA
    fit.time.md <- NA
  }
  
  if (Dmod == "~1"){
    outnames <- c("lambda0", "sigma", "D")
  } else if(Dmod == "~x^2"){
    outnames <- c("lambda0", "sigma", #"beta1", 
                  "beta2", "N")
  }
  assemble_CIs <- function(fit){
    if(fit$convergence == 0){
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
    } else { #if model did not converge, don't return values
      interval <- data.frame(name = outnames,
                             value = NA, 
                             upper = NA, 
                             lower = NA
      )
    }
    
    return(interval)
  } 
  
  out <- list(statdet_est = assemble_CIs(fit_sd), 
              movdet_est = assemble_CIs(fit_md),
              statdet_time = fit.time.sd,
              movdet_time = fit.time.md,
              stat_conv = fit_sd_con, #return convergence codes
              mov_conv = fit_md_con,
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
                    N,
                    hazdenom, 
                    Dmod){
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
                                  N,
                                  hazdenom, 
                                  mesh, 
                                  capthistout,
                                  Dmod)
  return(out)
}

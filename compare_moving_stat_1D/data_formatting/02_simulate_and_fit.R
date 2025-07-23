## 1 dimensional space

###-----------simulate data for moving and stationary detectors-----------------
###-------and fit model with stationary and moving detector---------------------

sim_fit <- function(tracksdf, 
                    traps,
                    trapcells,
                    dist_trapmesh,
                    tracksmeshdistmat,
                    useall,
                    lambda0, sigma, D_mesh, beta1, beta2,
                    hazdenom, 
                    mesh, 
                    meshunit,
                    Dmod = "~1"){
  start.time.sim <- Sys.time()
  #capthist dim(inds, traps)
  capthist_full <- sim_capthist(pop = NULL, 
                                traps, 
                                tracksdf,
                                lambda0, 
                                sigma, 
                                D_mesh,
                                hazdenom, #for hazard rate
                                mesh,
                                meshunit,
                                tracksmeshdistmat) 
  capthist <- capthist_full[which(apply((!is.na(capthist_full)), 1, sum)>0),,]
  
  #in case mesh is df
  mesh_mat <- as.matrix(mesh)
  
  #standard scr likelihood
  #occasions will be each survey (so for stationary, single trap per occasion)
  nocc <- length(unique(tracksdf$occ))
  
  #use
  induse_ls <- create_ind_use(capthist, trapcells, tracksdf)
  induse <- aperm(
    array(unlist(induse_ls), 
          dim =c(nrow(traps), nocc, dim(capthist)[1])), 
    c(3, 1, 2))
  
  
  #convert capthist to 1s and 0s
  capthist[is.na(capthist)] <- 0
  capthist[capthist!=0] <- 1
  
  fit.time.sim <- difftime(Sys.time(), start.time.sim, units = "secs")
  
  #grid with stationary detector likelihood
  if (Dmod == "~1"){
    start0 <- c( #logit(lambda0), 
      log(lambda0),
      log(sigma), log(D_mesh[1]))
    scaling_factors <- 10^round(log10(abs(start0)))
    start <- start0/scaling_factors
    
    stat_nll <- function(v_scaled){
      v <- v_scaled * scaling_factors 
      lambda0_ <- exp(v[1])#invlogit(v[1])
      sigma_ <- exp(v[2])
      D_mesh_ <- rep(exp(v[3]), nrow(mesh_mat))
      out <- negloglikelihood_stationary_cpp(lambda0_, sigma_,
                                             hazdenom, D_mesh_, 
                                             capthist, useall,
                                             dist_trapmesh, mesh_mat,
                                             linear = T)
      return(out)
    }
    #moving detector likelihood
    nll <- function(v_scaled){
      v <- v_scaled * scaling_factors 
      lambda0_ <- exp(v[1])#invlogit(v[1]) #logit link?
      sigma_ <- exp(v[2])
      D_mesh_ <- rep(exp(v[3]), nrow(mesh))
      out <- negloglikelihood_moving_cpp(lambda0_, sigma_,  
                                         hazdenom, D_mesh_,
                                         capthist, useall,
                                         induse, dist_trapmesh, mesh_mat,
                                         linear = T)
      return(out)
    }

    
  }else if(Dmod == "~x^2"){
    start0 <- c(#logit(lambda0),
      log(lambda0),
      log(sigma), beta1, beta2)
    scaling_factors <- 10^round(log10(abs(start0)))
    start <- start0/scaling_factors
    #quadratic density function
    stat_nll <- function(v_scaled){
      v <- v_scaled * scaling_factors 
      lambda0_ <- exp(v[1])#invlogit(v[1])
      sigma_ <- exp(v[2])
      D_mesh_ <- exp(v[3]*(mesh_mat[,1] + v[4])^2)
      out <- negloglikelihood_stationary_cpp(lambda0_, sigma_,
                                             hazdenom, D_mesh_, 
                                             capthist, useall,
                                             dist_trapmesh, mesh_mat,
                                             linear = T)
      return(out)
    }
    #moving detector likelihood
    nll <- function(v_scaled){
      v <- v_scaled * scaling_factors 
      lambda0_ <- exp(v[1])#invlogit(v[1])
      sigma_ <- exp(v[2])
      D_mesh_ <- exp(v[3]*(mesh_mat[,1] + v[4])^2)#exp(beta1*(mesh$x + beta2)^2)
      out <- negloglikelihood_moving_cpp(lambda0_, sigma_,  
                                         hazdenom, D_mesh_,
                                         capthist, useall,
                                         induse, dist_trapmesh, mesh_mat,
                                         linear = T)
      return(out)
    }
  }  
  
  start.time.sd <- Sys.time()
  fit_sd <- optim(par = start,
                  fn = stat_nll,
                  hessian = F, method = "Nelder-Mead") #NM is best at this likelihood, even though slower
  fit_sd$hessian <- numDeriv::hessian(stat_nll, x = fit_sd$par,
                                      method = "Richardson",
                                      method.args = list(eps = 1e-5, d = 1e-3, r = 3))
  fit.time.sd <- difftime(Sys.time(), start.time.sd, units = "secs") #includes hessian
  
  start.time.md <- Sys.time()
  fit_md <- optim(par = start,
                  fn = nll,
                  hessian = F, method = "Nelder-Mead")
  fit_md$hessian <- numDeriv::hessian(nll, x = fit_md$par,
                                      method = "Richardson",
                                      method.args = list(eps = 1e-5, d = 1e-3, r = 3))
  fit.time.md <- difftime(Sys.time(), start.time.md, units = "secs")
  
  if (Dmod == "~1"){
    outnames <- c("lambda0", "sigma", "D")
  } else if(Dmod == "~x^2"){
    outnames <- c("lambda0", "sigma", "beta1", "beta2")
  }
  assemble_CIs <- function(fit){
    fisher_info <- solve(fit$hessian)
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
              sim_time = fit.time.sim)
  return(out)
}


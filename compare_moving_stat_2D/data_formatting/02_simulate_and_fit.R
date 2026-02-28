## 2 dimensional space
###-----------simulate data for moving and stationary detectors-----------------
###-------and fit model with stationary and moving detector---------------------

simulate_popandcapthist <- function(tracksdf,
                                                             D_mesh,
                                                             lambda0,
                                                             sigma,
                                                             mesh,
                                                             traps,
                                                             trapspacing){
  
  start.time.sim <- Sys.time()
  
  # use 
  usetracksdf <- as.data.frame(tidyr::pivot_wider(
    tracksdf[,c("occ", "inc", "midx", "midy", "trapno")], 
    names_from = occ, 
    values_from = inc))
  usetracksdf[(is.na(usetracksdf))] <- 0
  colnames(usetracksdf)[c(1,2)] <- c("x","y")
  usetracksdf <- usetracksdf[
    rowSums(
      usetracksdf[,-c(which(colnames(usetracksdf) %in% 
                              c("x","y","trapno")))]) > 0,
  ]
  
  trackstrapscr <- read.traps(data = usetracksdf,
                              detector = "proximity",
                              binary.usage = FALSE)
  usage(trackstrapscr) <- usetracksdf[,-c(1,2,3)]
  
  popscr <- sim.popn(D = D_mesh * (10000), 
                     core = mesh, 
                     model2D = "IHP")
  
  # simulate with secr and proximity detectors 
  capthist_full <- sim.capthist(trackstrapscr,
                                popn = popscr,
                                detectfn = "HHN",
                                detectpar = list("lambda0" = lambda0,
                                                 "sigma" = sigma),
                                noccasions = ncol(usage(trackstrapscr)),
                                renumber = F)

  # delete detections after first 
  capthist <- apply(as.array(1:(dim(capthist_full)[1])), 1, 
                    function(i){
                      apply(as.array(1:(dim(capthist_full)[2])), 1, 
                            function(k){
                              capik <- capthist_full[i,k,]
                              #if detection happened that occasion
                              if(sum(capik) > 0){
                                if(sum(capik) == 1){
                                  thej <- which(capik > 0)
                                  dettime <- as.numeric(tracksdf[which(tracksdf$occ == k &
                                                              tracksdf$midx == trackstrapscr[thej,1] &
                                                              tracksdf$midy == trackstrapscr[thej,2]),
                                                      "time"])
                                } else {
                                  #which traps detected
                                  js <- which(capik > 0)
                                  #what were the coordinates
                                  jxys <- trackstrapscr[js,]
                                  jtimes <- array(NA, dim = nrow(jxys))
                                  #check rows of tracksdf to see times of those traps
                                  for(j in 1:nrow(jxys)){
                                    jtimes[j] <- tracksdf[which(tracksdf$occ == k &
                                                                  tracksdf$midx == jxys[j,1] &
                                                                  tracksdf$midy == jxys[j,2]),
                                                          "time"]
                                  }
                                  #keep only first one
                                  thej <- js[which.min(jtimes)] #index
                                  dettime <- min(jtimes)
                                }
                                
                                capik <- array(NA, dim = length(capik))
                                capik[thej] <- dettime
                                
                              } else {
                                #no detections
                                capik <- array(NA, dim = length(capik))
                              }
                              return(capik)
                            })
                    })
  dim(capthist) <- (dim(capthist_full)[c(3,2,1)])
  capthist <- aperm(capthist, c(3,2,1))
  
  capthist_times <- capthist
  
  # traps
  trapno_step <- usetracksdf$trapno
  
  # now scale back down to desired traps
  capthist_attrap_ls <- apply(as.array(1:(dim(capthist)[2])), 1, function(k){
    apply(as.array(1:(dim(capthist)[1])), 1, function(i){
      apply(as.array(1:length(unique(trapno_step))), 1, function(j){
        # index steps based on trap
        c_ikjs <- capthist_times[i,k,][trapno_step == j]
        # if all NA no detection happened
        if(all(is.na(c_ikjs))){
          ijk <- NA
          # otherwise keep detection time from steps for that trap
        } else {
          ijk <- c_ikjs[which(!is.na(c_ikjs))]
        }
      }, simplify = F)
    })
  })
  
  # rearrange from apply structure
  capthist_attrap <- aperm(array(
    unlist(capthist_attrap_ls),
    dim = c(length(unique(trapno_step)),
            (dim(capthist)[1]), 
            (dim(capthist)[2]))
  ), c(2,3,1))
  
  # calculate ind use for traps
  induse <- create_ind_use_C(capthist_attrap,
                             as.matrix(traps),
                             trapspacing,
                             tracksdf,
                             scenario = "everything")
  
  
  # convert capthist to 1s and 0s
  capthist_attrap[is.na(capthist_attrap)] <- 0
  capthist_attrap[capthist_attrap!=0] <- 1
  
  fit.time.sim <- difftime(Sys.time(), start.time.sim, units = "secs")
  out_ls <- list(capthist = capthist_attrap,
                 induse = induse,
                 fit.time.sim = fit.time.sim)
  return(out_ls)
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
                    Dmod,
                    meshspacing,
                    meanstepsize,
                    fitstat = FALSE
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
                                         induse, 
                                         dist_trapmesh, 
                                         mesh_mat,
                                         mesharea = meshspacing^2,
                                         meanstepsize)
      return(out)
    }

    
  }else if(Dmod == "~x^2"){  
    start0 <- c(
      log(lambda0),
      log(sigma),
      beta1, 
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
      
      D_mesh_  <- calcDv(mesh[,1] ,
                                     mesh[,2],
                                     v[3],
                                     v[4],
                                    exp(v[5]),
                                     meshspacing)
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
      
      D_mesh_  <- calcDv(mesh[,1] ,
                         mesh[,2],
                         v[3],
                         v[4],
                         exp(v[5]),
                         meshspacing)
      
      out <- negloglikelihood_moving_cpp(lambda0_, sigma_,  
                                         hazdenom, D_mesh_,
                                         capthist, useall,
                                         induse, 
                                         dist_trapmesh,
                                         mesh_mat,
                                         mesharea = meshspacing^2,
                                         meanstepsize)
      return(out)
    }
  }  
  
  if(fitstat){
    
    start.time.sd <- Sys.time()
    fit_sd <- nlm(p = start,
                    f = stat_nll,
                  hessian = FALSE) 
    fit_sd_con = fit_sd$code
    if(fit_sd_con == 1){
      fit_sd$hessian <- numDeriv::hessian(stat_nll, x = fit_sd$estimate,method = "Richardson",
                                          method.args = list(eps = 1e-5, d = 1e-3, r = 3))
      fit.time.sd <- difftime(Sys.time(), start.time.sd, units = "secs") #includes hessian
      
    } else { #if model did not converge, don't bother hessian 
      fit_sd$hessian = NA
      fit.time.sd = NA
    }
  } else {
    fit_sd <- list(par = NA, convergence = 4)
    fit.time.sd <- 99999999
    fit_sd_con <- 4
    fit_sd$hessian <- NA
  }
  
  start.time.md <- Sys.time()
  fit_md <- nlm(p = start,
                  f = nll,
                hessian = FALSE)
  fit_md_con = fit_md$code
  if(fit_md_con == 1){
    fit_md$hessian <- numDeriv::hessian(nll, x = fit_md$estimate,method = "Richardson",
                                        method.args = list(eps = 1e-6, d = 1e-4, r = 4))
    
    fit.time.md <- difftime(Sys.time(), start.time.md, units = "secs")
    
  } else {
    fit_md$hessian <- NA
    fit.time.md <- NA
  }
  
  if (Dmod == "~1"){
    outnames <- c("lambda0", "sigma", "D")
  } else if(Dmod == "~x^2"){
    outnames <- c("lambda0", "sigma", "beta1", 
                  "beta2", "N")
  }
  assemble_CIs <- function(fit){
    if(fit$code == 1){
      fisher_info <- MASS::ginv(fit$hessian)
      prop_sigma <- sqrt(diag(fisher_info))
      prop_sigma <- diag(prop_sigma)
      upper <- fit$estimate+1.96*prop_sigma
      lower <- fit$estimate-1.96*prop_sigma
      interval <- data.frame(name = outnames,
                             value = fit$estimate * scaling_factors, 
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

sim_fit <- function(traps_,
                    tracksdf_, 
                    mesh_, 
                    meshspacing_,
                    dist_trapmesh_,
                    useall_,
                    lambda0_, 
                    sigma_, 
                    D_mesh_, 
                    beta1_, 
                    beta2_,
                    N_,
                    hazdenom_, 
                    Dmod_,
                    trapspacing_,
                    meanstepsize_,
                    fitstat_ = FALSE){
  
  
  capthistout <- simulate_popandcapthist(tracksdf_,
                                         D_mesh_,
                                         lambda0_,
                                         sigma_,
                                         mesh_,
                                         traps_,
                                         trapspacing_)
  out <- fit_capthist(dist_trapmesh_,
                      useall_,
                      lambda0_,
                      sigma_,
                      D_mesh_,
                      beta1_,
                      beta2_,
                      N_,
                      hazdenom_,
                      mesh_,
                      capthistout,
                      Dmod_,
                      meshspacing_,
                      meanstepsize_,
                      fitstat_)
  return(out)
}

calc_nll <- function(dist_trapmesh,
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
                         Dmod,
                         meshspacing,
                         meanstepsize
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
                                         induse, 
                                         dist_trapmesh, 
                                         mesh_mat,
                                         mesharea = meshspacing^2,
                                         meanstepsize)
      return(out)
    }
    
    
  }else if(Dmod == "~x^2"){  
    start0 <- c(
      log(lambda0),
      log(sigma),
      beta1, 
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
      
      D_mesh_  <- calcDv(mesh[,1] ,
                         mesh[,2],
                         v[3],
                         v[4],
                         exp(v[5]),
                         meshspacing)
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
      
      D_mesh_  <- calcDv(mesh[,1] ,
                         mesh[,2],
                         v[3],
                         v[4],
                         exp(v[5]),
                         meshspacing)
      
      out <- negloglikelihood_moving_cpp(lambda0_, sigma_,  
                                         hazdenom, D_mesh_,
                                         capthist, useall,
                                         induse, 
                                         dist_trapmesh,
                                         mesh_mat,
                                         mesharea = meshspacing^2,
                                         meanstepsize)
      return(out)
    }
  }  
  
  out_sd  <- stat_nll(start)
  out_md <- nll(start)
  return(list(sd = out_sd,
              md = out_md))
    
} 


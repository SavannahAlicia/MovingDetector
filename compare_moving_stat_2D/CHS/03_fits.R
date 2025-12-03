
#------------------------------Moving detector --------------------------------


move_fit <- function(capthist,
                     tracksdf, 
                     trapcells,
                     dist_trapmesh,
                     useall,
                     induse,
                     DdesignX, 
                     hazdenom, 
                     mesh, 
                     startpar0){
  
  mesh_mat <- as.matrix(mesh)
  par0 <- unlist(startpar0)
  #rescale for easy hessian
  scaling_factors <- 10^round(log10(abs(par0)))
  par <- par0/scaling_factors
  lambda0parindex <- which(names(par) == "lambda0")
  sigmaparindex <- which(names(par) == "sigma")
  Dparindex <- grep("D", names(par))
  
   #stationary detector likelihood
   stat_nll <- function(v_scaled){
     v <- v_scaled * scaling_factors 
    lambda0_ <- exp(v[lambda0parindex])
    sigma_ <- exp(v[sigmaparindex])
    D_mesh_ <- exp(DdesignX %*% v[Dparindex]) 
     out <- negloglikelihood_stationary_cpp(lambda0_, sigma_,
                                            hazdenom, D_mesh_, 
                                            capthist, useall,
                                            dist_trapmesh, mesh_mat)
     return(out)
   }
  #moving detector likelihood
  nll <- function(v_scaled){
    v <- v_scaled * scaling_factors 
    lambda0_ <- exp(v[lambda0parindex])
    sigma_ <- exp(v[sigmaparindex])
    D_mesh_ <- exp(DdesignX %*% v[Dparindex]) 
    out <- negloglikelihood_moving_cpp(lambda0_, sigma_,  
                                       hazdenom, D_mesh_,
                                       capthist, useall,
                                       induse, dist_trapmesh, mesh_mat)
    return(out)
  }
  
   start.time.sd <- Sys.time()
   fit_sd <- optim(par = par,
                   fn = stat_nll,
                   hessian = F, method = "Nelder-Mead")
   fit_sd$hessian <- numDeriv::hessian(stat_nll, x = fit_sd$par,
                                       method = "Richardson",
                                       method.args = list(eps = 1e-5, d = 1e-3, r = 3))
   fit.time.sd <- difftime(Sys.time(), start.time.sd, units = "secs")
  
  start.time.md <- Sys.time()
  fit_md <- optim(par = par,
                  fn = nll,
                  hessian = F, method = "Nelder-Mead")
  fit_md$hessian <- numDeriv::hessian(nll, x = fit_md$par,
                                      method = "Richardson",
                                      method.args = list(eps = 1e-5, d = 1e-3, r = 3))
  fit.time.md <- difftime(Sys.time(), start.time.md, units = "secs")
  

  
  assemble_CIs <- function(fit){
    fisher_info <- MASS::ginv(fit$hessian)
    prop_sigma <- sqrt(diag(fisher_info))
    prop_sigma <- diag(prop_sigma)
    upper <- fit$par+1.96*prop_sigma
    lower <- fit$par-1.96*prop_sigma
    interval <- data.frame(name = names(par),
                           value = fit$par * scaling_factors, 
                           sd = diag(prop_sigma),
                           upper = diag(upper) * scaling_factors, 
                           lower = diag(lower) * scaling_factors
    )
    
    return(interval)
  } 
  
  out <- list(statdet_est = assemble_CIs(fit_sd), 
              movdet_est = assemble_CIs(fit_md),
              statdet_time = fit.time.sd,
              movdet_time = fit.time.md)
  return(out)
}




formulas <- list(D~s(x,y,k=5),
                 D~s(x,y,k=6),
                 D~s(x,y,k=7),
                 D~s(x,y,k=8),
                 D~s(x,y,k=9),
                 D~s(x, y, bs = "so", xt = list(bnd = bound)),
                 D~s(x, y, bs = "so", xt = list(bnd = bound)))

get_X_mat <- function(f, knots = NULL){
  formula <- formulas[[f]]
  split <- interpret.gam(formula)
  sml =  mgcv::smoothCon(split$smooth.spec[[1]], data = scrmesh, 
                         knots = knots, absorb.cons = T, 
                         scale.penalty = T, 
                         null.space.penalty = F,
                         sparse.cons = 0, 
                         diagonal.penalty = F, 
                         apply.by = T, 
                         modCon = 0)
  DdesignX = cbind( rep(1, nrow(scrmesh)), sml[[1]]$X)
  colnames(DdesignX) = c("(Intercept)", paste0("s(x,y).", 1:(ncol(DdesignX)-1)))
  return(DdesignX)
  }

fit_smooth <- function(f, startother = NULL, addtl_name = ""){
  formula = formulas[[f]]
  if(is.null(startother)){
    Dpar =  exp(fits[[f+1]]$fit$par[fits[[f+1]]$parindx$D]) 
    lambda0 = exp(fits[[f+1]]$fit$par[fits[[f+1]]$parindx$lambda0])
    sigma = exp(fits[[f+1]]$fit$par[fits[[f+1]]$parindx$sigma])
  } else {
    Dpar =  exp(fits[[startother]]$fit$par[fits[[startother]]$parindx$D]) 
    lambda0 = exp(fits[[startother]]$fit$par[fits[[startother]]$parindx$lambda0])
    sigma = exp(fits[[startother]]$fit$par[fits[[startother]]$parindx$sigma])
  }
  
  startparf = list(D = log(Dpar), 
                   lambda0 = log(lambda0), 
                   sigma = log(sigma))
  DdesignX <- Xmats[[f]]
  
  m_move <- move_fit(capthist = ch_10,
                     tracksdf = tracksdf, 
                     trapcells = trapcells,
                     dist_trapmesh = distmatscr,
                     useall = usage(trapscr),
                     induse = induse,
                     DdesignX = DdesignX, 
                     hazdenom = 1, 
                     mesh = scrmesh, 
                     startpar0 = startparf)
  saveRDS(m_move, file = paste("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/m_move", paste(formula)[3], addtl_name, ".Rds", sep = ""))
  return(m_move)
}

Xmats <- lapply(as.list(1:5), get_X_mat)
Xmats[[6]] <- get_X_mat(f = 6, knots = knots_soap)
Xmats[[7]] <- get_X_mat(f = 6, knots = knots_soap2)

myfits <- lapply(as.list(1:5), fit_smooth)
myfits[[6]] <- fit_smooth(6, startother = 4)
myfits[[7]] <- fit_smooth(7, startother = 8, addtl_name = "moreknots")
myfits <- list(readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/m_moves(x, y, k = 5).Rds"),
               readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/m_moves(x, y, k = 6).Rds"),
               readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/m_moves(x, y, k = 7).Rds"),
               readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/m_moves(x, y, k = 8).Rds"),
               readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/m_moves(x, y, k = 9).Rds"),
               readRDS('~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/m_moves(x, y, bs = "so", xt = list(bnd = bound)).Rds'),
               readRDS('~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/m_moves(x, y, bs = "so", xt = list(bnd = bound))moreknots.Rds'))

AICs <- apply(as.array(1:length(myfits)), 1, function(i){
  AIC_fn <- function(n,L){2*n + 2*L}
  parloc <- list(lambda0 = which(myfits[[i]]$movdet_est$name == "lambda0"),
                 sigma = which(myfits[[i]]$movdet_est$name == "sigma"),
                 D = which(grepl("D", myfits[[i]]$movdet_est$name)))

  movL <- negloglikelihood_moving_cpp(lambda0 = exp(myfits[[i]]$movdet_est$value[parloc$lambda0]), 
                                      sigma = exp(myfits[[i]]$movdet_est$value[parloc$sigma]),  
                              timeincr = 1, 
                              D_mesh = exp(Xmats[[i]] %*% myfits[[i]]$movdet_est$value[parloc$D]),
                              capthist= ch_10, 
                              usage = usage(trapscr),
                              indusage = induse, 
                              distmat =  distmatscr,
                              mesh = as.matrix(scrmesh))
  statL <- negloglikelihood_stationary_cpp(lambda0 = exp(myfits[[i]]$movdet_est$value[parloc$lambda0]), 
                                           sigma = exp(myfits[[i]]$movdet_est$value[parloc$sigma]),  
                                           timeincr = 1, 
                                           D_mesh = exp(Xmats[[i]] %*% myfits[[i]]$movdet_est$value[parloc$D]),
                                           capthist= ch_10, 
                                           usage = usage(trapscr),
                                           distmat =  distmatscr,
                                           mesh = as.matrix(scrmesh))
  statAIC <- AIC_fn(n = length(unlist(parloc)),
                    L = statL)
  movAIC <- AIC_fn(n = length(unlist(parloc)),
                   L = movL)
  return(
    data.frame(stat = statAIC, 
                    mov = movAIC, 
                    mod = paste(formulas[[i]])[c(3)])
         )
})
myAICs <- do.call(rbind, AICs)


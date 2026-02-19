#compare fits between my script and Abinand's
source("compare_moving_stat_2D/data_formatting/00_functions_and_parameters.R")
source("compare_moving_stat_2D/data_formatting/01_data_setup.R")
source("compare_moving_stat_2D/data_formatting/02_simulate_and_fit.R")
source("inst/BananSim.R")
trap_n_vert = round(ntrapsish/trap_n_horiz)
tracksteplength = round(trapspacing/nsteps_pertrap)

#create input objects for Abinand's
survObj <- simulateScrTrapsMask(
  nxTraps = trap_n_horiz,
  nyTraps = trap_n_vert,
  trapSpacing = trapspacing,
  nSteps = round(trapspacing/tracksteplength),
  maskSpacing = meshspacing,
  sigma = sigma, 
  N = N,
  b1 = beta1,
  b2 = beta2
)

trapSteps <- survObj$trapSteps
mask <- survObj$mask

surv_obj <- setup_data(sigma,
                         N,
                         beta1,
                         beta2,
                         ntrapsish,
                         trackxmin,
                         trapspacing,
                         meshspacing,
                         trap_n_horiz,
                         nsteps_pertrap,
                         occreps
  )
tracksdf <- surv_obj$tracksdf
traps <- surv_obj$traps
mesh <- surv_obj$mesh
D_mesh <- surv_obj$D_mesh_v

ch_out <- simulate_popandcapthist(traps,
                                    tracksdf, 
                                    lambda0,
                                    sigma,
                                    D_mesh,
                                    mesh,
                                    meshspacing = surv_obj$meshspacing,
                                    hazdenom = 1)

ch <- ch_out$capthist
induse <- ch_out$induse

p <- sim.popn(D = D_mesh * 100^2*2, 
              core = mask,
              model2D = 'IHP')
CH <- simCapthist(pop = p,
            trapSteps = trapSteps, 
            mask = mask, 
            lambda0 = lambda0*100, 
            sigma = sigma,
            nOccasionsTransect = occreps)
CH <- subset(CH, subset = 1:(dim(ch)[1]))

#replace CH with data from ch (or collapse ch by transect)
occkey <- unique(tracksdf[,c("occ", "transect", "rep")])

ch_sm <- sapply(as.list(1:occreps),  
      function(t){
        apply(ch[,occkey$rep == t,], c(1,3), sum)
      }, simplify = "array")
ch_sm <- aperm(ch_sm, c(1, 3, 2))

meanstepsize = mean(tracksdf$inc[tracksdf$inc !=0])

for(i in 1:(dim(CH)[1])){
  for(k in 1:(dim(CH)[2])){
    for(j in 1:(dim(CH)[3])){
      CH[i,k,j] <- ch_sm[i,k,j]
    }
  }
}

fit_abinand <-  scrFitMov(CH,
                                                mask, 
                                                trapSteps, 
                                                model = NULL, 
                                                startparams = NULL,
                                                hessian = F)

start0 <- c(
  log(lambda0),
  log(sigma),
  beta1, 
  beta2, 
  log(N))
scaling_factors <- rep(1, length(start0)) #10^round(log10(abs(start0)))
start <- start0/scaling_factors

#moving detector likelihood
nll <- function(v_scaled){
  v <- v_scaled * 1
  lambda0_ <- exp(v[1])
  sigma_ <- exp(v[2])
  eta <- v[3] * (mesh[,1] + v[4])^2
  exp_eta <- exp(eta)
  Z   <- sum(exp_eta) * meshspacing^2
  D_mesh_ <- exp(v[5]) * exp_eta / Z
  
  out <- negloglikelihood_moving_cpp(lambda0 = lambda0_, 
                                     sigma = sigma_,  
                                     haz_denom = 1,
                                     D_mesh = D_mesh_,
                                     capthist = ch, 
                                     usage = surv_obj$useall,
                                     indusage = induse, 
                                     distmat = surv_obj$dist_trapmesh,
                                     mesh = as.matrix(mesh),
                                     mesharea = meshspacing^2,
                                     meanstepsize = meanstepsize)
  return(out)
}
fit_md <- optim(par = start,
                fn = nll,
                hessian = F, method = "Nelder-Mead")

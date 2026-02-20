
Rcpp::sourceCpp("approx_movingdetectorlikelihood.cpp")
source("compare_moving_stat_2D/data_formatting/00_functions_and_parameters.R")
source("compare_moving_stat_2D/data_formatting/01_data_setup.R")
source("compare_moving_stat_2D/data_formatting/02_simulate_and_fit.R")
#calculate for one ind one occasion all fixed

nsims = 30
trap_n_vert = round(ntrapsish/trap_n_horiz)
tracksteplength = round(trapspacing/nsteps_pertrap)
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
meanstepsize <- mean(tracksdf$inc)
useall <- surv_obj$useall

pop <- matrix(c(-500, 500), nrow = 1, ncol = 2)

capthist_full <- sim_capthist_C(as.matrix(traps),
                                tracksdf, 
                                lambda0,
                                sigma,
                                D_mesh,
                                as.matrix(mesh),
                                meshspacing,
                                hazdenom=1,
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

capthist <- array(capthist, dim = c(1, (dim(capthist)[1]),(dim(capthist)[2])))
capthist_times <- capthist

#use
induse <- create_ind_use_C(capthist_times,
                           as.matrix(traps),
                           trapspacing, 
                           tracksdf,
                           scenario = "everything")

#convert capthist to 1s and 0s
capthist[is.na(capthist)] <- 0
capthist[capthist!=0] <- 1

out <- negloglikelihood_moving_cpp(lambda0 = lambda0, 
                                       sigma = sigma,  
                                       haz_denom = 1,
                                       D_mesh = D_mesh,
                                       capthist = capthist, 
                                       usage = surv_obj$useall,
                                       indusage = induse, 
                                       distmat = surv_obj$dist_trapmesh,
                                       mesh = as.matrix(mesh),
                                       mesharea = meshspacing^2,
                                       meanstepsize = meanstepsize)
    


#calculate by hand the probability of the capture history
Dprob_log <- array(NA, dim = nrow(mesh))
notseen_log <- sumallhuexcept <- didntsurvivej_log <- array(NA, dim = c(3800,4))
for(x in 1:nrow(mesh)){
  prob_ch_log <- rep(NA, nocc)
  for(k in 1:nocc){
    ikcap = F
    sum_hu <- 0
    hu_trapdet <- NA
    for(j in 1:nrow(traps)){
      d <- sqrt((traps[j,1] - mesh[x,1])^2 + (traps[j,2] - mesh[x,2])^2)
      hazard <- hazdist_cpp(lambda0, sigma, d, 1)
      hu <- hazard * induse[1, j, k]
      if(capthist[1,k,j] == 1){
        ikcap = T
        if(induse[1,j,k] > meanstepsize){
          hu_trapdet <- hazard * meanstepsize
          sum_hu = sum_hu + (hu - hu_trapdet)
        } else {
          hu_trapdet <- hu
        }
      } else {
        sum_hu = sum_hu + hu
      }
    }
    sumallhuexcept[x,k] <- sum_hu
    if(ikcap){
      #probability of capture history for this occasion given x
      # exp(-sum_hu) * (1 - exp(-hu_trapdet))
      prob_ch_log[k] <- -sum_hu + (log(1 - exp(-hu_trapdet)))
      notseen_log[x,k] <- -sum_hu - hu_trapdet
      didntsurvivej_log[x,k] <- (log(1 - exp(-hu_trapdet)))
    } else {
      prob_ch_log[k] <- -sum_hu 
      notseen_log[x,k] <- -sum_hu
    }

  }
  sum_prob_ch_log <- sum(prob_ch_log)
  Dprob_log[x] <- log(D_mesh[x]) + sum_prob_ch_log 
}

DKdat <- data.frame(x = mesh$x,
                    y = mesh$y,
                    Dprobch = exp(out$DKprod_eachx_log[,1]),
                    myprobch = exp(Dprob_log),
                    D = D_mesh)
ggplot() +
  geom_raster(DKdat, 
              mapping = aes(x = x, y = y, fill = Dprobch/D_mesh)) +
  geom_point(data = data.frame(x = traps$x,
                               y = traps$y,
                               dets = apply(capthist, 3, sum)),
             mapping = aes(x = x,
                           y =y,
                           color = dets)) +
  coord_equal() +
  scale_color_viridis_c()+
  geom_point(data.frame(x = pop[1,1], y = pop[1,2]),
             mapping = aes(x = x, y = y), color= "red")
ggplot() +
  geom_raster(DKdat, 
              mapping = aes(x = x, y = y, fill = myprobch/D_mesh)) +
  geom_point(data = data.frame(x = traps$x,
                               y = traps$y,
                               dets = apply(capthist, 3, sum)),
             mapping = aes(x = x,
                           y =y,
                           color = dets)) +
  coord_equal() +
  scale_color_viridis_c()+
  geom_point(data.frame(x = pop[1,1], y = pop[1,2]),
             mapping = aes(x = x, y = y), color= "red")

outDpdot <- negloglikelihood_moving_cpp(lambda0 = lambda0, 
                                           sigma = sigma,  
                                           haz_denom = 1,
                                           D_mesh = D_mesh,
                                           capthist = capthist, 
                                           usage = surv_obj$useall,
                                           indusage = induse, 
                                           distmat = surv_obj$dist_trapmesh,
                                           mesh = as.matrix(mesh),
                                           mesharea = meshspacing^2,
                                           meanstepsize = meanstepsize)


Dpdot <- array(NA, nrow(mesh))
for(x in 1:nrow(mesh)){
  prob_notcap_log <- rep(NA, nocc)
  for(k in 1:nocc){
    sum_hu <- 0
    for(j in 1:nrow(traps)){
      d <- sqrt((traps[j,1] - mesh[x,1])^2 + (traps[j,2] - mesh[x,2])^2)
      hazard <- hazdist_cpp(lambda0, sigma, d, 1)
      hu <- hazard * useall[j, k]
      sum_hu = sum_hu + hu
    }
      prob_notcap_log[k] <- -sum_hu 
  }
  sum_prob_notcap_log <- sum(prob_notcap_log)
  pdot <- 1 - exp(sum_prob_notcap_log)
  Dpdot[x] <- (D_mesh[x]) * pdot
}

DKdat$myDpdot <- Dpdot
DKdat$outDpdot <- outDpdot

ggplot() +
  geom_raster(DKdat, 
              mapping = aes(x = x, y = y, fill = outDpdot/D_mesh)) +
  geom_point(data = data.frame(x = traps$x,
                               y = traps$y,
                               dets = apply(capthist, 3, sum)),
             mapping = aes(x = x,
                           y =y,
                           color = dets)) +
  scale_color_viridis_c()+
  geom_point(data.frame(x = pop[1,1], y = pop[1,2]),
             mapping = aes(x = x, y = y), color= "red")
ggplot() +
  geom_raster(DKdat, 
              mapping = aes(x = x, y = y, fill = myDpdot/D_mesh)) +
  geom_point(data = data.frame(x = traps$x,
                               y = traps$y,
                               dets = apply(capthist, 3, sum)),
             mapping = aes(x = x,
                           y =y,
                           color = dets)) +
  scale_color_viridis_c()+
  geom_point(data.frame(x = pop[1,1], y = pop[1,2]),
             mapping = aes(x = x, y = y), color= "red")

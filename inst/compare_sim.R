######### trying different simulator
#compare fits between my script and Abinand's
source("compare_moving_stat_2D/data_formatting/00_functions_and_parameters.R")
source("compare_moving_stat_2D/data_formatting/01_data_setup.R")
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
meanstepsize = mean(tracksdf$inc[tracksdf$inc !=0])


sim_capthist_bysubset <- function(tracksdf,
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
  
  # ggplot() +
  #   geom_raster(data = data.frame(x = mesh$x,
  #                                 y = mesh$y,
  #                                 D = D_mesh),
  #               mapping = aes(x = x, y = y,
  #                             fill = D)) +
  #   scale_fill_viridis_c() +
  #   geom_point(data = tracksdf,
  #              mapping = aes(x = x, 
  #                            y = y),
  #              color = "white", shape = "+", size = 5) +
  #   geom_point(data = data.frame(x = trackstrapscr$x, 
  #                                y = trackstrapscr$y),
  #              mapping = aes(x = x, y = y),
  #              color = "red",
  #              size = 1) +
  #   scale_x_continuous(
  #     #limits = c(-1100,-900)
  #   )
  
  # simulate with secr and proximity detectors 
  capthist_full <- sim.capthist(trackstrapscr,
                                popn = popscr,
                                detectfn = "HHN",
                                detectpar = list("lambda0" = lambda0,
                                                 "sigma" = sigma),
                                noccasions = ncol(usage(trackstrapscr)),
                                renumber = F)
  
  # add zero capthists back on for visualization
  zero_ch <- array(0, dim = c(nrow(popscr) - (dim(capthist_full)[1]), 
                              dim(capthist_full)[2],
                              dim(capthist_full)[3]))
  rownames(zero_ch) <- (1:nrow(popscr))[which(!1:nrow(popscr) %in% rownames(capthist_full))]
  capthist_full <- abind::abind(capthist_full, zero_ch, along = 1)
  capthist_full <- capthist_full[order(as.numeric(rownames(capthist_full))),,]
  
  # 
  # ggplot() +
  #   scale_fill_viridis_c() +
  #   geom_raster(data = data.frame(x = mesh$x,
  #                                 y = mesh$y,
  #                                 D = D_mesh),
  #               mapping = aes(x = x, y = y,
  #                             fill = D),
  #               alpha = 0.3) +
  #   new_scale_fill() +
  #   geom_point(data = data.frame(x = popscr$x,
  #                                y = popscr$y,
  #                                dets = apply(capthist_full, 1, sum)),
  #              mapping = aes(x = x, 
  #                            y = y, 
  #                            fill = dets),
  #              shape = 21) +
  #   scale_fill_viridis_c(option = "magma") +
  #   geom_point(data = data.frame(x = trackstrapscr$x, 
  #                                y = trackstrapscr$y),
  #              mapping = aes(x = x, y = y),
  #              color = "red", 
  #              shape = "+",
  #              size = 1) +
  #   scale_x_continuous(
  #     #limits = c(-1100,-900)
  #   )
  
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
                                  dettime <- tracksdf[which(tracksdf$occ == k &
                                                              tracksdf$midx == trackstrapscr[thej,1] &
                                                              tracksdf$midy == trackstrapscr[thej,2]),
                                                      "time"]
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
  
  # dim(capthist)
  # summary(apply(capthist_full, 1, function(x){sum(x>0)}))     # detections per individual
  # 
  # ggplot() +
  #   scale_fill_viridis_c() +
  #   geom_raster(data = data.frame(x = mesh$x,
  #                                 y = mesh$y,
  #                                 D = D_mesh),
  #               mapping = aes(x = x, y = y,
  #                             fill = D),
  #               alpha = 0.3) +
  #   new_scale_fill() +
  #   geom_point(data = data.frame(x = trackstrapscr$x, 
  #                                y = trackstrapscr$y),
  #              mapping = aes(x = x, y = y),
  #              color = "white", 
  #              shape = 19,
  #              size = 1) +
  #   geom_point(data = data.frame(x = popscr$x,
  #                                y = popscr$y,
  #                                dets = apply(capthist, 1, sum)),
  #              mapping = aes(x = x, 
  #                            y = y, 
  #                            fill = dets),
  #              shape = 21) +
  #   scale_fill_viridis_c(option = "magma") +
  #   geom_point(data = traps,
  #              mapping = aes(x = x, y = y),
  #              color = "red") +
  #   scale_x_continuous(
  #     #limits = c(-1100,-900)
  #   )
  # 
  
  # delete zero capture histories
  if(length(which(apply((!is.na(capthist)), 1, sum)>0)) == 0){
    warning("Empty capture history.")
    capthist <- array(NA, dim = c(1, nocc, nrow(traps)))
  } else {
    capthist <- capthist[which(apply((!is.na(capthist)), 1, sum)>0),,]
  }
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

ch_bysub <- sim_capthist_bysubset(tracksdf,
                                  D_mesh,
                                  lambda0,
                                  sigma,
                                  mesh,
                                  traps,
                                  trapspacing)

saveRDS(ch_bysub$capthist, "inst/ch.Rds")
saveRDS(ch_bysub$induse, "inst/induse.Rds")

# ch_bysub_ls <- lapply(1:10, function(x){
#   sim_capthist_bysubset(tracksdf,
#                         D_mesh,
#                         lambda0,
#                         sigma,
#                         mesh,
#                         traps,
#                         trapspacing)$capthist})
# 
# ch_bystep_ls <- lapply(1:10, function(x){
#   simulate_popandcapthist(traps,
#                         tracksdf, 
#                         lambda0,
#                         sigma,
#                         D_mesh,
#                         mesh,
#                         meshspacing,
#                         hazdenom)$capthist})
# 
# summary(unlist(lapply(ch_bysub_ls, function(ch){
#   (nrow(ch))
#   })))
# summary(unlist(lapply(ch_bystep_ls, function(ch){
#   (nrow(ch))
# })))
# bigtrackstrapscr <-  read.traps(data = traps,
#                           detector = "multi",
#                           binary.usage = FALSE)
# usage(bigtrackstrapscr) <- useall
# 
# pdot <- pdot(as.matrix(mesh),
#      traps = bigtrackstrapscr,
#      detectfn = "HHN",
#      detectpar = list(lambda0 = lambda0, 
#                       sigma = sigma),
#      noccasions = length(unique(tracksdf$occ))
#      )
# sum(pdot * D_mesh)*meshspacing^2
# 
# 
# dim(ch_bysub$capthist)
# summary(apply(ch_bysub$capthist, 1, function(x){sum((x))}))    
# dim(ch_bystep$capthist)
# summary(apply(ch_bystep$capthist, 1, function(x){sum((x))}))    
# 


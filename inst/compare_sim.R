######### trying different simulator


sim_capthist_bysubset <- function(tracksdf,
                               D_mesh,
                               lambda0,
                               sigma,
                               mesh,
                               traps,
                               trapspacing){
  
  start.time.sim <- Sys.time()
  
  #tracksdf[,c("midx", "midy")] <- calc_trackmidpts(tracksdf)
  usetracksdf <- as.data.frame(tidyr::pivot_wider(tracksdf[,c("occ", "inc", "midx", "midy")], 
                                                  names_from = occ, values_from = inc))
  usetracksdf[(is.na(usetracksdf))] <- 0
  colnames(usetracksdf)[c(1,2)] <- c("x","y")
  usetracksdf <- usetracksdf[rowSums(usetracksdf[,-c(1,2)]) > 0,]
  
  trackstrapscr <- read.traps(data = usetracksdf,
                        detector = "proximity",
                        binary.usage = FALSE)
  usage(trackstrapscr) <- usetracksdf[,-c(1,2)]
  
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
  
  capthist_full <- sim.capthist(trackstrapscr,
                                popn = popscr,
                                detectfn = "HHN",
                                detectpar = list("lambda0" = lambda0,
                                                 "sigma" = sigma),
                                noccasions = ncol(usage(trackstrapscr)),
                                renumber = F)
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
  
  #delete detections after first 
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
  
  if(length(which(apply((!is.na(capthist)), 1, sum)>0)) == 0){
    warning("Empty capture history.")
    capthist <- array(NA, dim = c(1, nocc, nrow(traps)))
  } else {
    capthist <- capthist[which(apply((!is.na(capthist)), 1, sum)>0),,]
  }
  capthist_times <- capthist
  
  #convert capthist to 1s and 0s
  capthist[is.na(capthist)] <- 0
  capthist[capthist!=0] <- 1

  
  # allocate track records to grid
  trapno_step <- rep(0, nrow(trackstrapscr))
  for (tr in 1:length(trapno_step)) {
    d <- sqrt((trackstrapscr[tr, ]$x - traps[, 1])^2 + (trackstrapscr[tr, ]$y - traps[, 2])^2)
    dmin <- min(d)
    #this needs to assign to the first trap
    candidatetrap <- which(d==dmin) #but this returns first trap
    if(length(candidatetrap)==1){
      trapno_step[tr] <- candidatetrap
    } else { #if there are two candidate traps
      #we need to choose the first one
      #check if there is an earlier tracksdfpt
      if(tr>1){ #if its not the first point in tracksdf
        #calculate the midpoint between two trackpoints
        midx <- mean(trackstrapscr[tr,"x"],trackstrapscr[tr-1, "x"])
        #return the trap of two candidates closest to that
        trapno_step[tr] <- candidatetrap[which.min(c(gr[candidatetrap[1],1], gr[candidatetrap[2],1])-midx)]
      } else {
        #if this if the first point in tracksdf, it doesn't matter
        trapno_step[tr] <- candidatetrap[2]
      }
      #
    }
    
  }
  
  #now scale back down to desired trap
  apply(as.array(1:(dim(capthist)[2])), 1, function(k){
    apply(as.array(1:(dim(capthist)[1])), 1, function(i){
      
      
    })
  })
  
  #condense capture history to traps
  X <- matrix(capthist, nrow = (dim(capthist)[1]) * (dim(capthist)[2]),
              ncol = (dim(capthist)[3]))
  
  X_new <- rowsum(t(X), trapno_step)
  X_newt <- t(X_new)
  capthist_attrap <- array(X_newt, dim = c((dim(capthist)[1]),
                                           (dim(capthist)[2]),
                                           (nrow(traps))))
  #I need to condense before calculating induse!!!
  #use
  # induse <- create_ind_use_C(capthist,
  #                            as.matrix(traps),
  #                            trapspacing, 
  #                            tracksdf,
  #                            scenario = "everything")

  fit.time.sim <- difftime(Sys.time(), start.time.sim, units = "secs")
  out_ls <- list(capthist = capthist_attrap,
                 induse = induse,
                 fit.time.sim = fit.time.sim)
  return(out_ls)
}

ch_bysub <- sim_capthist_bysubset(tracksdf, D_mesh_v,
                                  lambda0, sigma, mesh, traps,
                                  trapspacing)

ch_bysub_ls <- lapply(1:300, function(x){
  sim_capthist_bysubset(tracksdf, D_mesh_v,
                      lambda0, sigma, mesh, traps)$capthist})

ch_bystep_ls <- lapply(1:300, function(x){
  simulate_popandcapthist(traps,
                        tracksdf, 
                        lambda0,
                        sigma,
                        D_mesh_v,
                        mesh,
                        meshspacing,
                        hazdenom)$capthist})

summary(unlist(lapply(ch_bysub_ls, function(ch){
  (nrow(ch))
  })))
summary(unlist(lapply(ch_bystep_ls, function(ch){
  (nrow(ch))
})))
bigtrackstrapscr <-  read.traps(data = traps,
                          detector = "multi",
                          binary.usage = FALSE)
usage(bigtrackstrapscr) <- useall

pdot <- pdot(as.matrix(mesh),
     traps = bigtrackstrapscr,
     detectfn = "HHN",
     detectpar = list(lambda0 = lambda0, 
                      sigma = sigma),
     noccasions = length(unique(tracksdf$occ))
     )
sum(pdot * D_mesh_v)*meshspacing^2


dim(ch_bysub$capthist)
summary(apply(ch_bysub$capthist, 1, function(x){sum((x))}))    
dim(ch_bystep$capthist)
summary(apply(ch_bystep$capthist, 1, function(x){sum((x))}))    



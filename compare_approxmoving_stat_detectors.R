#compare moving and stationary detector fits
library(secr)
library(lubridate)
library(parallel)
library(ggplot2)
library(sp)
library(sf)
setwd("~/Documents/UniStAndrews/MovingDetector")
#source("movingdetectorlikelihood.R")
Rcpp::sourceCpp("approx_movingdetectorlikelihood.cpp")


set.seed(12345)
nsims = 100

#-------------------------------functions---------------------------------------
#' Create polygons for each box centered at grid pt
#' 
#' @param grid dataframe with x and y columns specifying grid pts
#' @param spacing the spacing between grid pts, if not specified will calculate
#' @return SpatialPolygons object in list form 
#' @export
create_grid_polygons <- function(grid, spacing = NULL){
  #so we only have to create gridbboxes and lines list once

  #specify corners of box centered on grid pt
  grid_spacing <- spacing
  grid_bboxes <- data.frame(
    left_bound = grid$x - (grid_spacing/2),
    right_bound = grid$x + (grid_spacing/2), 
    lower_bound = grid$y - (grid_spacing/2),
    upper_bound = grid$y + (grid_spacing/2)
  )
  #create a polygon for box centered on grid pt
  create_polygon_for_grid_row <- function(grid_bboxes_row){
    upperleft <- as.numeric(c(grid_bboxes_row[1], grid_bboxes_row[4]))
    lowerleft <- as.numeric(c(grid_bboxes_row[1], grid_bboxes_row[3]))
    upperright <- as.numeric(c(grid_bboxes_row[2], grid_bboxes_row[4]))
    lowerright <- as.numeric(c(grid_bboxes_row[2], grid_bboxes_row[3]))
    boxpts <- rbind(upperleft, upperright, lowerright, lowerleft)
    colnames(boxpts) <- c("x", "y")
    rownames(boxpts) <- NULL
    p = sp::Polygon(boxpts)
    ps = sp::Polygons(list(p),1)
    sps = sp::SpatialPolygons(list(ps))
    proj4string(sps) <- sp::CRS(sf::st_crs(26916)$input)
    return(sps)
  }
  #create polygon box for all grid pts
  grid_bboxes$polygons  <- apply(X = grid_bboxes, MARGIN = 1, 
                                 FUN = create_polygon_for_grid_row)
  return(grid_bboxes)
}

#' Create named list of SpatialLines where name is track ID (occ)
#' @param tracksdf dataframe with xy, occ, trapno, and time
#' @return named list of SpatialLines objects
#' @export
create_line_spatlines <- function(tracksdf, scenario = "everything",
                                  projto = sp::CRS(sf::st_crs(26916)$input)){
  #create a line for between each pt of track, put into named list by ID
  #check that it's not the last point in a track
  #(we don't want a line between tracks)
  check_pt <- function(row, tracksdf. = tracksdf){
    if (tracksdf.$occ[row] == tracksdf.$occ[row+1]){
      newline <- Line(tracksdf.[c(row, row+1),c("x", "y")])
    } else {
      newline <- NULL
    }
    return(newline)
  }
  linelist <- lapply(X = as.list(1:(nrow(tracksdf)-1)), FUN = check_pt)
  names(linelist) <- tracksdf$occ[1:(length(tracksdf$occ)-1)]
  linelist <- linelist[which(!lapply(linelist, is.null) == TRUE)]
  #Testing: check there are lines for each track segment, but not between tracksdf
  #dim(tracksdf)[1] - length(unique(tracksdf$occ)) == length(linelist)
  #for each track ID, create a Lines object by scenario
  #append that Lines object to a list of them
  big_Lines_list <- list()
  nocc <- length(unique(tracksdf$occ))
  for (trackindex in 1:nocc){
    #print(trackindex)
    trackIDi <- unique(tracksdf$occ)[trackindex]
    thetrack <- tracksdf[which(tracksdf$occ == trackIDi),]
    thelines <- linelist[which(names(linelist) == paste(trackIDi))]
    
    #check scenario
    if (scenario == "onison"){ #on effort is on effort
      newlines <- list(Lines(thelines[which(
        thetrack[1:(nrow(thetrack)-1),]$effort == "OnEffort")], 
        ID = paste(trackIDi)))
    } else if (scenario == "everything"){ #everything is on effort
      newlines <- list(Lines(thelines, ID = paste(trackIDi)))
    }
    
    #just creating desired class of object
    newlines_spat <- SpatialLines(newlines)
    #specify projection
    proj4string(newlines_spat) <- projto
    #adding new objects to list
    big_Lines_list <- append(big_Lines_list, newlines_spat)
    #naming list for track ID
    names(big_Lines_list)[trackindex] <- paste(trackIDi)
  }
  return(big_Lines_list)
}

#' Get length line that cross polygon
#' @param polygon a spatial polygon object
#' @param sp_line a Spatial Lines object
#' @return returns the length of the line within polygon
#' @export
get_length <- function(polygon, sp_line){
  boxlength = 0
  #attempts to chop line to just the bit within the polygon
  theintersection <- st_intersection(st_as_sf(polygon), st_as_sf(sp_line))
  #only calculate length of the line if there is a line that exists
  if(length(theintersection$geometry) > 0){
    boxlength <- as.numeric(st_length(theintersection))
  } #otherwise the length remains 0
  return(boxlength)
}

#'Calculate track lengths in grid boxes
#' function to find intersection between lines and grid box 
#' @param tracks_sp_lines list of lines labelled by track ID in name field
#' @param grid_polygons dataframe with a polygon in each row corresponding to 
#'                      that grid pt
#' @return column of summed intersections, each row is grid pt (trap)
#'         and sum is for all tracks specified in function
lengths_in_grid <- function(tracks_sp_lines, tracksinoc, grid_polygons){
  trackspl_to_use <- tracks_sp_lines[which(names(tracks_sp_lines) %in%
                                             as.character(tracksinoc))]
  
  if (length(trackspl_to_use) == 1){#just one track per occ now 
    inters <-  sapply(X = grid_polygons$polygons, 
                      FUN = function(x){get_length(x, trackspl_to_use[[1]])})
  } else {
    warning(paste("track length", paste0(tracksinoc, collapse = " "), "not equal to 1"))
    inters = rep(NA, length(grid_polygons$polygons))
  }
  return(inters)
}

#' Simulate moving detector capture history
#' 
#' @param pop defaults to null, population locations from sim.popn
#' @param traps dataframe of x y trap coordinates for trap grid
#' @param tracksdf dataframe with track locations (xy), occ, time, and trap 
#'                  number for corresponding trap grid row
#' @param lambda0 
#' @param sigma 
#' @param D_mesh density at each mesh pt (for simulating population if not specified)
#' 
#' @return capture history array of dimension individual by trap by occ that 
#' contains NA if not captured or time of first capture (seconds since the first
#' time logged for that occasion) if captured
#' @export
sim_capthist <- function(pop = NULL, 
                         traps, 
                         tracksdf,
                         lambda0, 
                         sigma, 
                         D_mesh,
                         hazdenom, #for hazard rate
                         report_probseenxk = FALSE){
  if(is.null(pop)){
    pop <- sim.popn(D = D_mesh, core = mesh, model2D = "IHP", 
                    Ndist = "poisson", buffertype = "rect")
    rownames(pop) <- NULL
  }
  #this is now distances between pop hrcs and track locations
  dist_dat_pop <- apply(as.matrix(tracksdf[,c("x","y")]), 1, function(trapj){
    apply(as.matrix(pop), 1, function(popi){
    dist(rbind(trapj,
               popi), method = "euclidean")
  }) }) #ind by track location
  
  capthist <- lapply(as.list(1:nrow(pop)), #for each individual
                     FUN = function(i){
                       lapply(as.list(unique(tracksdf$occ)),
                              FUN = function(occk){
                                #prob individual is seen this occasion (survey)
                                trackoccdf <- tracksdf[tracksdf$occ == occk,]
                                begintimek <-min(trackoccdf$time)
                                #NOTE if you want a distance denominated hazard, use length diffs for hazdenom
                                #if you want rate per unit time, use time increment in seconds
                                increments <- #c(0, 
                                              # difftime(trackoccdf$time[-1], 
                                              #          trackoccdf$time[-nrow(trackoccdf)], 
                                              #          unit = "secs"))
                                            c(0, apply(as.array(1:(nrow(trackoccdf)-1)), 1, function(x){
                                              dist(trackoccdf[x:(x+1), c("x", "y")])}))
                                hazs <- apply(as.array(dist_dat_pop[i,tracksdf$occ == occk]),
                                              1, FUN = function(d){
                                                hazdist_cpp(lambda0, sigma, d, hazdenom)
                                                })
                                integ <- sum(hazs * (increments/hazdenom)) #hazard is per time incr, need how many increments (or if dist, per dist incr)
                                survk <- exp(-1 * integ)
                                probseenxk <- 1 - survk 
                                if (report_probseenxk) { 
                                  return(probseenxk)
                                } else {
                                  seenxk_bool <- rbinom(1,1, probseenxk)
                                  if (seenxk_bool){
                                    survive_until_t <- exp(-1 * (cumsum(hazs * (increments/hazdenom))))
                                    survive_t_inc <- exp(-1 * hazs * (increments/hazdenom))
                                    seenfirstat_t <- survive_until_t[-length(survive_until_t)] * (1 - survive_t_inc[-1])
                                    seenfirstatt_givenseen <- c(0, seenfirstat_t/probseenxk)
                                    #need to record time of detection at traps as defined in grid
                                    det <- sample(x = c(1:nrow(trackoccdf)), size = 1, replace = T,  prob = seenfirstatt_givenseen)
                                    trapdet <- trackoccdf[det,"trapno"]
                                    timedet <- difftime(trackoccdf[det, "time"], begintimek, units = "secs")
                                    capik <- rep((NA), nrow(traps))
                                    capik[trapdet] <- timedet
                                    
                                  } else {
                                    capik <- rep((NA), nrow(traps))
                                  }
                                  
                                  return(capik) 
                                }
                        
                              })
                     }) 
  if(report_probseenxk){
    capthist_array <- t(
      array(unlist(capthist), 
            dim = c(nocc,nrow(pop))) 
    ) #ind by occasion
  } else {
    capthist_array <- aperm(
        array(unlist(capthist), 
              dim =c(nrow(traps), nocc,  nrow(pop))), 
        c(3, 2, 1)) #ind x occ x traps

  }
  return(capthist_array)
}

#----------------------------------data setup-----------------------------------
lambda0 = .4
sigma = 300
#multiple tracklines, keep seperate occasions since we only take first detection
#per occasion

#each trackline is a series of points with x, y, and time
tracksteps = 55 #intervals 
trackint = 360 #seconds
tracksdf <- rbind(
  data.frame(occ = 1,
             x = seq(from = 1500, to = 3000, length.out = tracksteps+1),
             y = 1250, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
  data.frame(occ = 2,
             x = seq(from = 1500, to = 3000, length.out = tracksteps+1),
             y = 1450, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
  data.frame(occ = 3,
             x = seq(from = 1500, to = 3000, length.out = tracksteps+1),
             y = 1650, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
  data.frame(occ = 4,
             x = seq(from = 1500, to = 3000, length.out = tracksteps+1),
             y = 1850, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                        by = trackint))
  
)
nocc <- length(unique(tracksdf$occ))

#mesh grid
meshspacing = 250
mesh <- make.mask(tracksdf[,c("x","y")], buffer = 3*sigma, spacing = meshspacing)

D_mesh <- rep(.4, nrow(mesh))
beta1 <- -(1/40000)
beta2 <- -2250
D_mesh_q <- exp(beta1*(mesh$x + beta2)^2)
hazdenom <- 1 #hazard is per time or distance, currently specified as distance

#trap grid
trapspacing = meshspacing
xgr <- seq(min(tracksdf$x), max(tracksdf$x), by = trapspacing)
ygr <- seq(min(tracksdf$y), max(tracksdf$y), by = trapspacing)
gr <- expand.grid(xgr, ygr)
# allocate track records to grid
trapno <- rep(0, nrow(tracksdf))
for (tr in 1:length(trapno)) {
  cat(tr, " / ", length(trapno), "\r")
  d <- sqrt((tracksdf[tr, ]$x - gr[, 1])^2 + (tracksdf[tr, ]$y - gr[, 2])^2)
  dmin <- min(d)
  trapno[tr] <- which.min(d)
}

# pick out possible traps as used cells
un <- sort(unique(trapno))
tracksdf$trapno <- apply(as.array(1:length(trapno)), 1, FUN = function(x){which(un == trapno[x])})
traps <- data.frame(x = gr[un, 1],
                    y = gr[un, 2])

#example population
expop <- sim.popn(D = D_mesh_q, core = mesh, model2D = "IHP", 
                Ndist = "poisson", buffertype = "rect")
excapthist <- sim_capthist(pop = expop, traps, tracksdf, lambda0, sigma, D_mesh, hazdenom, report_probseenxk = F)
exch <- excapthist[which(apply((!is.na(excapthist)), 1, sum)>0),,]

ggplot() +
  geom_point(data.frame(x = mesh$x, y = mesh$y, D = D_mesh_q), mapping = aes(x = x, y = y, alpha = D), shape = 21) + 
  geom_point(data.frame(x = expop$x, y = expop$y, dets = apply((!is.na(excapthist)), 1, sum)),
             mapping = aes(x = x, y = y, color = as.factor(dets)), size = 2) +
  scale_color_viridis_d() +
  geom_line(tracksdf, mapping = aes(x = x, y = y, group = occ)) +
  geom_point(data.frame(x = traps$x, y = traps$y, dets = as.factor(apply((!is.na(excapthist)), 3, sum))),
             mapping = aes(x = x, y = y), shape = 19)


#create trapgrid (will break tracklines into trap grid, keeping times)
trap_cells <- create_grid_polygons(traps, spacing = trapspacing)
#all lines
tracklines <- create_line_spatlines(tracksdf)

#use for undetected inds
getuse <- function(oc){
  usecol <- lengths_in_grid(tracklines, oc, trap_cells)
  return(usecol)
}
useall <- matrix(0, nr = nrow(trap_cells), nc = nocc)
colnames(useall) <- 1:nocc
useall[,c(1:ncol(useall))] <- do.call(cbind,
                                   mclapply(X= as.list(1:ncol(useall)), 
                                            FUN = getuse, mc.cores = 3))

#for each individual, create a use matrix
create_ind_use <- function(ch, trapcells, tracksdf){
  use <- mclapply(as.list(1:dim(ch)[1]), FUN = function(i){
   apply(as.array(1:dim(ch)[2]), 1, function(k){
     trackoccdf <- tracksdf[tracksdf$occ == k,]
     if(!all(is.na(ch[i,k,]))){
       #if i detected k, discard survey pts after detection in occasion k
       dettime <- min(trackoccdf$time) + ch[i,k,which(!is.na(ch[i,k,]))]
       trackoccdf <- trackoccdf[trackoccdf$time <= dettime,]
       #add up length of trackline in each grid cell
       useik <- lengths_in_grid(create_line_spatlines(trackoccdf), k, trap_cells)
     } else {#if i wasn't detected in k, all traps used full amount, so skip subsetting
      useik <- useall[,k]
      }
       #(if there are covariates, perhaps could create matrix for them here, but for now ignore)
   })
}, mc.cores = 3)
  return(use)
}

#assemble into a i,j,k use matrix
induse_ls <- create_ind_use(exch, trapcells, tracksdf) #takes about 30 seconds
induse <- aperm(
  array(unlist(induse_ls), 
        dim =c(nrow(traps), nocc, nrow(expop))), 
  c(3, 1, 2))

#calculate distance matrix for all trap cells and mesh cells
dist_trapmesh <- apply(as.matrix(mesh), 1, function(meshm){
  apply(as.matrix(traps), 1, function(trapj){
    dist(rbind(trapj,
               meshm), method = "euclidean")
  }) })




###-----------simulate data for moving and stationary detectors-----------------
###-------and fit model with stationary and moving detector---------------------

sim_fit <- function(tracksdf, 
                    traps,
                    trapcells,
                    dist_trapmesh,
                    useall,
                    lambda0, sigma, D_mesh,
                    hazdenom, 
                    mesh, 
                    Dmod = "~1"){
  start.time.sim <- Sys.time()
  #capthist dim(inds, traps)
  capthist_full<- sim_capthist(pop = NULL, traps,
                               tracksdf, lambda0, sigma, D_mesh,
                               hazdenom) #function in movingdetectorlikelihood.R
  capthist <- capthist_full[which(apply((!is.na(capthist_full)), 1, sum)>0),,]
  
  #in case mesh is df
  mesh <- as.matrix(mesh)
  
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
    stat_nll <- function(v){
      lambda0_ <- exp(v[1])
      sigma_ <- exp(v[2])
      D_mesh_ <- rep(exp(v[3]), nrow(mesh))
      out <- negloglikelihood_stationary_cpp(lambda0_, sigma_,
                                             hazdenom, D_mesh_, 
                                             capthist, useall,
                                             dist_trapmesh, mesh)
      return(out)
    }
    #moving detector likelihood
    nll <- function(v){
      lambda0_ <- exp(v[1])
      sigma_ <- exp(v[2])
      D_mesh_ <- rep(exp(v[3]), nrow(mesh))
      out <- negloglikelihood_moving_cpp(lambda0_, sigma_,  
                                         hazdenom, D_mesh_,
                                         capthist, useall,
                                         induse, dist_trapmesh, mesh)
      return(out)
    }
    start <- c( log(.4), log(300), log(.4))
    
  }else if(Dmod == "~x^2"){
    #quadratic density function
    stat_nll <- function(v){
      lambda0_ <- exp(v[1])
      sigma_ <- exp(v[2])
      D_mesh_ <- exp(v[3]*(mesh[,1] + v[4])^2)
      out <- negloglikelihood_stationary_cpp(lambda0_, sigma_,
                                             hazdenom, D_mesh_, 
                                             capthist, useall,
                                             dist_trapmesh, mesh)
      return(out)
    }
    #moving detector likelihood
    nll <- function(v){
      lambda0_ <- exp(v[1])
      sigma_ <- exp(v[2])
      D_mesh_ <- exp(v[3]*(mesh[,1] + v[4])^2)#exp(beta1*(mesh$x + beta2)^2)
      out <- negloglikelihood_moving_cpp(lambda0_, sigma_,  
                                         hazdenom, D_mesh_,
                                         capthist, useall,
                                         induse, dist_trapmesh, mesh)
      return(out)
    }
    start <- c(log(.4), log(300), -(1/40000), -1250)
  }  

  start.time.sd <- Sys.time()
  fit_sd <- optim(par = start,
                  fn = stat_nll,
                  hessian = F, method = "Nelder-Mead")
  fit.time.sd <- difftime(Sys.time(), start.time.sd, units = "secs")
  
  start.time.md <- Sys.time()
  fit_md <- optim(par = start,
               fn = nll,
               hessian = F, method = "Nelder-Mead")
  fit.time.md <- difftime(Sys.time(), start.time.md, units = "secs")
  
  if (Dmod == "~1"){
    outnames <- c("lambda0", "sigma", "D")
  } else if(Dmod == "~x^2"){
    outnames <- c("lambda0", "sigma", "beta1", "beta2")
  }
  assemble_CIs <- function(fit){
    #fisher_info <- solve(fit$hessian)
    #prop_sigma <- sqrt(diag(fisher_info))
    #prop_sigma <- diag(prop_sigma)
    #upper <- fit$par+1.96*prop_sigma
    #lower <- fit$par-1.96*prop_sigma
    interval <- data.frame(name = outnames,
                           value = fit$par#, 
                           #upper = diag(upper), 
                           #lower = diag(lower)
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

#test one of each
fit1 <- sim_fit(tracksdf, 
                traps,
                trapcells,
                dist_trapmesh,
                useall,
                lambda0, sigma, D_mesh,
                hazdenom, 
                mesh, 
                Dmod = "~1")
fit2 <- sim_fit(tracksdf, 
                traps,
                trapcells,
                dist_trapmesh,
                useall,
                lambda0, sigma, D_mesh_q,
                hazdenom, 
                mesh, Dmod = "~x^2")

start.time.all_q <- Sys.time()
all_sim_fits_q <- mclapply(X = as.list(1:nsims),
                           FUN = function(sim){
                             return(sim_fit(tracksdf, 
                                            traps,
                                            trapcells,
                                            dist_trapmesh,
                                            useall,
                                            lambda0, sigma, D_mesh_q,
                                            hazdenom, 
                                            mesh, 
                                            Dmod = "~x^2"))
                           },
                           mc.cores = 6
)
tot.time.all_q <- difftime(Sys.time(), start.time.all_q, units = "secs")



start.time.all <- Sys.time()
all_sim_fits <- mclapply(X = as.list(1:nsims),
                   FUN = function(sim){
                     return(sim_fit(tracksdf, 
                                    traps,
                                    trapcells,
                                    dist_trapmesh,
                                    useall,
                                    lambda0, sigma, D_mesh,
                                    hazdenom, 
                                    mesh, 
                                    Dmod = "~1"))
                   },
                   mc.cores = 6
)
tot.time.all <- difftime(Sys.time(), start.time.all, units = "secs")


###------------------------compare computation time-----------------------------



###--------------------------compare precision ---------------------------------
make_plot_dat<- function(all_sim_fits){
  
stat_outs <- do.call(rbind,lapply(as.list(1:length(all_sim_fits)), FUN = function(x){
  df <- all_sim_fits[[x]]$statdet_est
  df$sim = rep(x,nrow(all_sim_fits[[1]]$statdet_est))
  return(df)
  }))
move_outs <-  do.call(rbind,lapply(as.list(1:length(all_sim_fits)), FUN = function(x){
  df <- all_sim_fits[[x]]$movdet_est
  df$sim = rep(x,nrow(all_sim_fits[[1]]$statdet_est))
  return(df)
}))

all_outs <- rbind(cbind(stat_outs, data.frame(model = rep("stationary", 
                                                          nrow(stat_outs)))),
      cbind(move_outs, data.frame(model = rep("moving",nrow(move_outs)))))

library(dplyr)
all_outs2 <- all_outs %>%
  group_by(name, model) %>%
  summarize(mean = mean(value), meanupper = quantile(value, probs = .975), meanlower = quantile(value, probs = .025))


out <- list(all_outs= all_outs, all_outs2 =all_outs2)

}

plotdat <- make_plot_dat(all_sim_fits_q)
plotdat <- make_plot_dat(all_sim_fits)
all_outs <- plotdat$all_outs
all_outs2 <- plotdat$all_outs2
plotcols <- c("#178A28", "#D81B60")

ggplot() +
  geom_density(all_outs[all_outs$name == "lambda0",], mapping = aes(x = exp(value), col = model), size = 1.3) +
  geom_vline(data = all_outs2[all_outs2$name == "lambda0",], aes(xintercept = exp(c(mean)), col = model), size = 1.3) +
  geom_vline(data = all_outs2[all_outs2$name == "lambda0",], aes(xintercept = exp(c(meanlower)), col = model), linetype = "dashed", size = 1.3) +
  geom_vline(data = all_outs2[all_outs2$name == "lambda0",], aes(xintercept = exp(c(meanupper)), col = model), linetype = "dashed", size = 1.3) +
  geom_vline(xintercept = lambda0, size = 1.3, col = "black") +
  xlab("lambda0") +
  scale_color_manual(values = plotcols) +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))

  
ggplot() +
  geom_density(all_outs[all_outs$name == "sigma",], mapping = aes(x = exp(value), col = model), size = 1.3) +
  geom_vline(data = all_outs2[all_outs2$name == "sigma",], aes(xintercept = exp(mean), col = model), size = 1.3) +
  geom_vline(data = all_outs2[all_outs2$name == "sigma",], aes(xintercept = exp(meanlower), col = model), linetype = "dashed", size = 1.3) +
  geom_vline(data = all_outs2[all_outs2$name == "sigma",], aes(xintercept = exp(meanupper), col = model), linetype = "dashed", size = 1.3) +
  geom_vline(xintercept = sigma, size = 1.3, col = "black")+
  scale_color_manual(values = plotcols) +
  xlab("sigma") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
ggplot() +
  geom_density(all_outs[all_outs$name == "beta1",], mapping = aes(x = value, col = model), size = 1.3) +
  geom_vline(data = all_outs2[all_outs2$name == "beta1",], aes(xintercept = c(mean), col = model), size = 1.3) +
  geom_vline(data = all_outs2[all_outs2$name == "beta1",], aes(xintercept = c(meanlower), col = model), linetype = "dashed", size = 1.3) +
  geom_vline(data = all_outs2[all_outs2$name == "beta1",], aes(xintercept = c(meanupper), col = model), linetype = "dashed", size = 1.3) +
  geom_vline(xintercept = beta1) +
  scale_color_manual(values = plotcols) +
  xlab("beta1") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))
ggplot() +
  geom_density(all_outs[all_outs$name == "beta2",], mapping = aes(x = value, col = model)) +
  geom_vline(data = all_outs2[all_outs2$name == "beta2",], aes(xintercept = c(mean), col = model)) +
  geom_vline(data = all_outs2[all_outs2$name == "beta2",], aes(xintercept = c(meanlower), col = model), linetype = "dashed") +
  geom_vline(data = all_outs2[all_outs2$name == "beta2",], aes(xintercept = c(meanupper), col = model), linetype = "dashed") +
  geom_vline(xintercept = beta2) +
  scale_color_manual(values = plotcols) +
  xlab("beta2") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))

D_plotdat <- data.frame(x = rep(seq(min(mesh$x), max(mesh$x), 50),3),
                        y = rep(rep(mesh$y[1], 61),3),
                        trueD = rep(exp(beta1*(seq(min(mesh$x), max(mesh$x), 50) + beta2)^2), 3),
                        stationarydets = c(exp(as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                         all_outs2$name == "beta1", "mean"]) * 
                                               (seq(min(mesh$x), max(mesh$x), 50) + 
                                                  as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                     all_outs2$name == "beta2", "mean"]))^2),
                                           exp(as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                      all_outs2$name == "beta1", "meanupper"]) * 
                                                 (seq(min(mesh$x), max(mesh$x), 50) + 
                                                    as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                           all_outs2$name == "beta2", "meanupper"]))^2),
                                           exp(as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                      all_outs2$name == "beta1", "meanlower"]) * 
                                                 (seq(min(mesh$x), max(mesh$x), 50) + 
                                                    as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                           all_outs2$name == "beta2", "meanlower"]))^2)),
                        movingdets = c(exp(as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                     all_outs2$name == "beta1", "mean"]) * 
                                           (seq(min(mesh$x), max(mesh$x), 50) + 
                                              as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                 all_outs2$name == "beta2", "mean"]))^2),
                                       exp(as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                  all_outs2$name == "beta1", "meanupper"]) * 
                                             (seq(min(mesh$x), max(mesh$x), 50) + 
                                                as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                       all_outs2$name == "beta2", "meanupper"]))^2),
                                       exp(as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                  all_outs2$name == "beta1", "meanlower"]) * 
                                             (seq(min(mesh$x), max(mesh$x), 50) + 
                                                as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                       all_outs2$name == "beta2", "meanlower"]))^2)),
                        quantile = c(rep("mean", 61), rep("2.5%",61), rep("97.5%",61))
                        )

D_plotdat <- data.frame(x = rep(mesh$x, 3),
                        y = rep(mesh$y, 3),
                        trueD = rep(D_mesh, 3),
                        stationarydets = c(rep(exp(as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                           all_outs2$name == "D", "mean"])), nrow(mesh)),
                                           rep(exp(as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                      all_outs2$name == "D", "meanupper"])), nrow(mesh)),
                                           rep(exp(as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                      all_outs2$name == "D", "meanlower"])), nrow(mesh))),
                        movingdets = c(rep(exp(as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                       all_outs2$name == "D", "mean"])), nrow(mesh)),
                                       rep(exp(as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                  all_outs2$name == "D", "meanupper"])), nrow(mesh)),
                                       rep(exp(as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                  all_outs2$name == "D", "meanlower"])), nrow(mesh)) 
                                       ),
                        quantile = c(rep("mean", nrow(mesh)), rep("2.5%", nrow(mesh)), rep("97.5%",nrow(mesh)))
)

ggplot() + 
  geom_tile(data = D_plotdat[,c(1,2,3)], aes(x = x, y= y, fill = trueD)) +
  scale_color_manual(values = "black",  name = "D") +
  theme_classic() +
  labs(title = "True density") +
  theme(axis.text.y = element_blank()
        ,
        legend.position = "none") 
ggplot() + 
  geom_tile(data = D_plotdat[,c(1,2,4)], aes(x = x, y= y, fill = stationarydets)) +
  scale_fill_viridis_c(name = "D") +
  theme_classic() +
  labs(title = "Stationary detectors") +
  theme(axis.text.y = element_blank(),
        legend.position = "none") 
ggplot() + 
  geom_tile(data = D_plotdat[,c(1,2,5)], aes(x = x, y= y, fill = movingdets)) +
  scale_fill_viridis_c(name = "D") +
  theme_classic() +
  labs(title = "Moving detectors") +
  theme(axis.text.y = element_blank(),
        legend.position = "none") 


D_plotdatlong <- tidyr::pivot_longer(D_plotdat, cols = c("trueD", "stationarydets", "movingdets"))

ggplot() + 
  geom_line(data = D_plotdatlong[D_plotdatlong$quantile == "mean",], mapping = aes(x = x, y = value, col = name), size = 1.3) +
  geom_line(data = D_plotdatlong[D_plotdatlong$quantile == "2.5%",], mapping = aes(x = x, y = value, col = name), linetype = "dashed", size = 1.3) +
  geom_line(data = D_plotdatlong[D_plotdatlong$quantile == "97.5%",], mapping = aes(x = x, y = value, col = name), linetype = "dashed", size = 1.3) +
  scale_color_manual(values = c(plotcols, "black"), labels = c("moving", "stationary", "true D")) +
  guides(col=guide_legend(title="model")) +
  #ylim(0,.5) +
  xlim(500,1750)+
  ylab("Density") +
  theme_classic() +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20))


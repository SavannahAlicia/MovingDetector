#compare moving and stationary detector fits
library(secr)
library(lubridate)
library(ggplot2)
setwd("~/Documents/UniStAndrews/MovingDetector")
Rcpp::sourceCpp("functions.cpp")
source("movingdetectorlikelihood.R")

set.seed(12345)
nsims = 100

#----------------------------------data setup-----------------------------------
lambda0 = .4
sigma = 300
meshgrid <- expand.grid(x = seq(0, 2600, 300), y = seq(0,2600, 300))
mesh <- make.mask(meshgrid, buffer = 800, spacing = 250)
D_mesh <- rep(.4, nrow(mesh))

timeincr = 6*60 #specifies time increments for integration AND for distance matrix 

tracksteps = 10

#-----------------------------detectors moving left to right--------------------

trapsdf <- rbind(
  data.frame(trapID = 1, 
             occ = 1,
             x = seq(from = 750, to = 1750, length.out = tracksteps),
             y = 750, 
             time = seq(ymd_hms("2024-01-01 8:00:00"), (ymd_hms("2024-01-01 8:00:00") + (tracksteps-1)*360), 
                        by = 360)),
  data.frame(trapID = 2, 
             occ = 1,
             x = seq(from = 750, to = 1750, length.out = tracksteps),
             y = 1050, 
             time = seq(ymd_hms("2024-01-01 8:00:00"), (ymd_hms("2024-01-01 8:00:00") + (tracksteps-1)*360), 
                        by = 360)),
  data.frame(trapID = 3, 
             occ = 1,
             x = seq(from = 750, to = 1750, length.out = tracksteps),
             y = 1375, 
             time = seq(ymd_hms("2024-01-01 8:00:00"), (ymd_hms("2024-01-01 8:00:00") + (tracksteps-1)*360), 
                        by = 360)),
  data.frame(trapID = 4, 
             occ = 1,
             x = seq(from = 750, to = 1750, length.out = tracksteps),
             y = 1750, 
             time = seq(ymd_hms("2024-01-01 8:00:00"), (ymd_hms("2024-01-01 8:00:00") + (tracksteps-1)*360), 
                        by = 360)),
  data.frame(trapID = 5, 
             occ = 1,
             x = seq(from = 750, to = 1750, length.out = tracksteps),
             y = 1750, 
             time = seq(ymd_hms("2024-01-02 8:00:00"), (ymd_hms("2024-01-02 8:00:00") + (tracksteps-1)*360), 
                        by = 360))
  
)

dist_dat <- create_distdat(trapsdf, mesh) #calls function in movingdetectorlikelihood.R to create data object
ggplot() +
  geom_point(data = mesh, aes(x = x, y = y), col = "grey") +
  geom_point(data = trapsdf, aes(x = x, y = y, col = as.factor(trapID)), shape = "+", size = 5) +
  theme_classic()

###-----------simulate data for moving and stationary detectors-----------------
###-------and fit model with stationary and moving detector---------------------
sim_fit <- function(trapsdf, dist_dat, lambda0, sigma, D_mesh, timeincr, mesh){
  #capthist dim(inds, traps)
  capthist_full<- sim_capthist(pop = NULL, trapsdf, timeincr, lambda0, sigma, D_mesh) #function in movingdetectorlikelihood.R
  capthist <- capthist_full[which(rowSums(capthist_full, na.rm = T) > 0),]
  
  #standard scr likelihood
  #occasions will be each detector (so for stationary, single trap per occasion)
  occn <- nrow(dist_dat$traps)
  
  #traps are every location
  trapxy <- trapsdf[!duplicated(trapsdf[,c("x","y")]), c("x","y")]
  rownames(trapxy) <-  paste("trap", rownames(trapxy))
  secrtraps <- read.traps(data = trapxy, detector = "multi")
  usage(secrtraps) <- t(apply(as.array(1:nrow(trapxy)), 1,
                     FUN = function(traprow){
                       #usage for traps is 1 for occasions that have a time
                       occT <- unique(trapsdf[which(trapsdf$x == secrtraps[traprow,"x"] & 
                                              trapsdf$y == secrtraps[traprow,"y"]) ,
                                      "trapID"]) 
                       userow <- rep(0, occn)
                       userow[occT] <- 1
                       return(userow)
                       }))
  
  
  #secr capthist must be session, ID, occasion, trap
  #each time in capthist corresponds to an individual and occasion
  #and the time/occasion can be used to find the xy in trapsdf

  secrcapdata <- data.frame(session = numeric(),
                        ID = numeric(),
                        occasion = numeric(),
                        trap = numeric())
  for(ind in 1:nrow(capthist)){
    for(trapo in 1:ncol(capthist)){
      if(!is.na(capthist[ind,trapo])){
        trapxyid <- trapsdf[which(trapsdf$trapID == trapo &
                                  trapsdf$time == capthist[ind,trapo]), c("x","y","occ")]
        secrtrap <- which(secrtraps$x == trapxyid$x & 
                            secrtraps$y == trapxyid$y)
        out <- data.frame(session = 1,
                          ID = ind,
                          occasion = trapo,
                          trap = secrtrap)
        secrcapdata <- rbind(secrcapdata, out)
      }
      
    }
  }
  secrcapdata$trap <- paste("trap", secrcapdata$trap )
  secrcapthist <- make.capthist(secrcapdata, secrtraps, fmt = c("trapID"))
  distmat <-   #distance data object (created from traps and mesh)
    traptomesh <- apply(as.array(1:nrow(mesh)), 1, FUN = 
                          function(meshrow){
                            apply(as.array(1:nrow(trapxy)), 
                                  1, FUN = 
                                    function(traprow){
                                      dist(rbind(trapxy[traprow, c("x", "y")],
                                                 mesh[meshrow,c("x","y")]), method = "euclidean")
                                    })})
  

  #grid with stationary detector likelihood
  stat_nll <- function(v){
    lambda0_ <- v[1]
    sigma_ <- v[2]
    D_mesh_ <- rep(v[3], nrow(mesh))
    out <- negloglikelihood_stationary_cpp(lambda0_, sigma_, D_mesh_, secrcapthist, usage(secrtraps), distmat, dist_dat)
    return(out)
  }
  start <- c( .3, 200, .3)
  start.time.sd <- Sys.time()
  fit_sd <- optim(par = start,
               fn = stat_nll,
               hessian = T)
  fit.time.sd <- difftime(Sys.time(), start.time.sd, units = "secs")
  
  #moving detector likelihood
  nll <- function(v){
    lambda0_ <- v[1]
    sigma_ <- v[2]
    D_mesh_ <- rep(v[3], nrow(mesh))
    out <- negloglikelihood_cpp(lambda0_, sigma_, D_mesh_, timeincr, capthist, dist_dat)
    return(out)
  }
  
  start <- c( .3, 200, .3)
  start.time.md <- Sys.time()
  fit_md <- optim(par = start,
               fn = nll,
               hessian = T)
  fit.time.md <- difftime(Sys.time(), start.time.md, units = "secs")
  
  assemble_CIs <- function(fit){
    fisher_info <- solve(fit$hessian)
    prop_sigma <- sqrt(diag(fisher_info))
    prop_sigma <- diag(prop_sigma)
    upper <- fit$par+1.96*prop_sigma
    lower <- fit$par-1.96*prop_sigma
    interval <- data.frame(name = c("lambda0", "sigma", "D"),
                           value = fit$par, 
                           upper = diag(upper), 
                           lower = diag(lower))
    
    return(interval)
  } 
  
  out <- list(statdet_est = assemble_CIs(fit_sd), 
              movdet_est = assemble_CIs(fit_md),
              statdet_time = fit.time.sd,
              movdet_time = fit.time.md)
  return(out)
}

start.time.all <- Sys.time()
all_sim_fits <- lapply(as.list(1:nsims),
                   FUN = function(sim){
                     return(sim_fit(trapsdf, dist_dat, lambda0, sigma, D_mesh, timeincr, mesh))
                   }
)
tot.time.all <- difftime(Sys.time(), start.time.all, units = "secs")


###------------------------compare computation time-----------------------------



###--------------------------compare precision ---------------------------------

stat_outs <- do.call(rbind,lapply(as.list(1:length(all_sim_fits)), FUN = function(x){
  df <- all_sim_fits[[x]]$statdet_est
  df$sim = rep(x,3)
  return(df)
  }))
move_outs <-  do.call(rbind,lapply(as.list(1:length(all_sim_fits)), FUN = function(x){
  df <- all_sim_fits[[x]]$movdet_est
  df$sim = rep(x,3)
  return(df)
}))

ggplot() +
  geom_density(stat_outs[stat_outs$name == "lambda0",], mapping = aes(x = value)) +
  geom_vline(xintercept = mean(stat_outs[stat_outs$name == "lambda0","value"])) +
  geom_vline(xintercept = quantile(stat_outs[stat_outs$name == "lambda0","value"], 
                                   probs = c(0.025, .975)), linetype = "dashed") +
  geom_density(move_outs[move_outs$name == "lambda0",], mapping = aes(x = value),
               col = "blue") +
  geom_vline(xintercept = mean(move_outs[move_outs$name == "lambda0","value"]),
             col = "blue") +
  geom_vline(xintercept = quantile(move_outs[move_outs$name == "lambda0","value"], 
                                   probs = c(0.025, .975)), linetype = "dashed",
             col = "blue") +
  geom_vline(xintercept = lambda0, col = "red")
  
               



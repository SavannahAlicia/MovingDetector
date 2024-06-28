#compare moving and stationary detector fits
library(secr)
library(lubridate)
library(ggplot2)
setwd("~/Documents/UniStAndrews/MovingDetector")
Rcpp::sourceCpp("functions.cpp")
source("movingdetectorlikelihood.R")

set.seed(12345)
numsims = 100

#----------------------------------data setup-----------------------------------
lambda0 = .4
sigma = 300
meshgrid <- expand.grid(x = seq(0, 2600, 300), y = seq(0,2600, 300))
mesh <- make.mask(meshgrid, buffer = 800, spacing = 250)
D_mesh <- rep(.4, nrow(mesh))

timeincr = 6*60 #specifies time increments for integration AND for distance matrix 

tracksteps = 10
#-------------------------stationary detectors-----------------------------
# #dataframe with column for trapid, x, y, and time
# stat_trapsdf <- rbind(
#   data.frame(trapID = 1, 
#              occ = 1,
#              x = rep(750, tracksteps),
#              y = 1000, 
#              time = seq(ymd_hms("2024-01-01 8:00:00"), (ymd_hms("2024-01-01 8:00:00") + (tracksteps-1)*360), 
#                         by = 360)),
#   data.frame(trapID = 2, 
#              occ = 1,
#              x = rep(750, tracksteps),
#              y = 1250, 
#              time = seq(ymd_hms("2024-01-01 8:00:00"), (ymd_hms("2024-01-01 8:00:00") + (tracksteps-1)*360), 
#                         by = 360)),
#   data.frame(trapID = 3, 
#              occ = 1,
#              x = rep(750, tracksteps),
#              y = 1500, 
#              time = seq(ymd_hms("2024-01-01 8:00:00"), (ymd_hms("2024-01-01 8:00:00") + (tracksteps-1)*360), 
#                         by = 360)),
#   data.frame(trapID = 4, 
#              occ = 1,
#              x = rep(1250, tracksteps),
#              y = 1000, 
#              time = seq(ymd_hms("2024-01-01 8:00:00"), (ymd_hms("2024-01-01 8:00:00") + (tracksteps-1)*360), 
#                         by = 360)),
#   data.frame(trapID = 5, 
#              occ = 1,
#              x = rep(1250, tracksteps),
#              y = 1250, 
#              time = seq(ymd_hms("2024-01-01 8:00:00"), (ymd_hms("2024-01-01 8:00:00") + (tracksteps-1)*360), 
#                         by = 360)),
#   data.frame(trapID = 6, 
#              occ = 1,
#              x = rep(1250, tracksteps),
#              y = 1500, 
#              time = seq(ymd_hms("2024-01-01 8:00:00"), (ymd_hms("2024-01-01 8:00:00") + (tracksteps-1)*360), 
#                         by = 360)),
#   data.frame(trapID = 7, 
#              occ = 1,
#              x = rep(1750, tracksteps),
#              y = 1000, 
#              time = seq(ymd_hms("2024-01-01 8:00:00"), (ymd_hms("2024-01-01 8:00:00") + (tracksteps-1)*360), 
#                         by = 360)),
#   data.frame(trapID = 8, 
#              occ = 2,
#              x = rep(1750, tracksteps),
#              y = 1250, 
#              time = seq(ymd_hms("2024-01-02 8:00:00"), (ymd_hms("2024-01-02 8:00:00") + (tracksteps-1)*360), 
#                         by = 360)),
#   data.frame(trapID = 9, 
#              occ = 2,
#              x = rep(1750, tracksteps),
#              y = 1500, 
#              time = seq(ymd_hms("2024-01-02 8:00:00"), (ymd_hms("2024-01-02 8:00:00") + (tracksteps-1)*360), 
#                         by = 360))
# )
# 
# stat_dist_dat <- create_distdat(stat_trapsdf, mesh) #calls function in movingdetectorlikelihood.R to create data object
# ggplot() +
#   geom_point(data = mesh, aes(x = x, y = y), col = "grey") +
#   geom_point(data = stat_trapsdf, aes(x = x, y = y, col = as.factor(trapID)), shape = "+", size = 5) +
#   theme_classic()
  
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
# #add 4 more occasions, each a day later
# add <- data.frame(trapID = rep(5:20, each = 4*tracksteps),
#                   occ = rep(2:5, each = nrow(ltr_trapsdf)),
#                   x = rep(ltr_trapsdf$x, 4),
#                   y = rep(ltr_trapsdf$y, 4),
#                   time = c(ltr_trapsdf$time+ 60*60*24,
#                            ltr_trapsdf$time+ 60*60*24*2,
#                            ltr_trapsdf$time+ 60*60*24*3,
#                            ltr_trapsdf$time+ 60*60*24*4))
# ltr_trapsdf <- rbind(ltr_trapsdf, add)


dist_dat <- create_distdat(ltr_trapsdf, mesh) #calls function in movingdetectorlikelihood.R to create data object
ggplot() +
  geom_point(data = mesh, aes(x = x, y = y), col = "grey") +
  geom_point(data = trapsdf, aes(x = x, y = y, col = as.factor(trapID)), shape = "+", size = 5) +
  theme_classic()

###-----------simulate data for moving and stationary detectors-----------------
###-------and fit model with stationary and moving detector---------------------
sim_fit <- function(trapsdf, dist_dat, lambda0, sigma, D_mesh, timeincr){
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
        trapxy <- trapsdf[which(trapsdf$trapID == trapo &
                                  trapsdf$time == capthist[ind,trapo]), c("x","y","occ")]
        secrtrap <- which(secrtraps$x == trapxy$x & 
                            secrtraps$y == trapxy$y)
        out <- data.frame(session = 1,
                          ID = ind,
                          occasion = trapo,
                          trap = secrtrap)
        secrcapdata <- rbind(secrcapdata, out)
      }
      
    }
  }
  make.capthist(secrcapdata, secrtraps, fmt = c("trapID"))
                
  
  #moving detector likelihood
  nll <- function(v){
    lambda0_ <- v[1]
    sigma_ <- v[2]
    D_mesh_ <- rep(v[3], nrow(mesh))
    out <- negloglikelihood_cpp(lambda0_, sigma_, D_mesh_, timeincr, capthist, dist_dat)
    return(out)
  }
  
  start <- c( .3, 200, .3)
  fit <- optim(par = start,
               fn = nll,
               hessian = T)
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

start.time <- Sys.time()
stat_fit <- sim_fit(stat_trapsdf, stat_dist_dat, lambda0, sigma, D_mesh, timeincr)
tot.time <- Sys.time()- start.time


start.time2 <- Sys.time()
ltr_fit <- sim_fit(ltr_trapsdf, ltr_dist_dat, lambda0, sigma, D_mesh, timeincr)
tot.time2 <- Sys.time()- start.time2


###------------------------compare computation time-----------------------------



###--------------------------compare precision ---------------------------------
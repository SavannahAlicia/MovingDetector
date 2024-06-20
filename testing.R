library(secr)
library(lubridate)
library(ggplot2)
setwd("~/Documents/UniStAndrews/MovingDetector")
Rcpp::sourceCpp("functions.cpp")
source("movingdetectorlikelihood.R")

set.seed(555)

######data setup
lambda0 = .5
sigma = 300
timeincr = 6*60 #specifies time increments for integration AND for distance matrix 
meshgrid <- expand.grid(x = seq(0, 2600, 300), y = seq(0,2600, 300))
mesh <- make.mask(meshgrid, buffer = 800, spacing = 250)
D_mesh <- rep(.4, nrow(mesh))

#traps
#specify K moving traps
#dataframe with column for trapid, x, y, and time
trapsdf <- rbind(
  data.frame(trapID = 1, 
             x = seq(750, 1750, 250),
             #x = rep(750, 5),
             y = 1000, 
             time = seq(ymd_hms("2024-01-01 8:00:00"), (ymd_hms("2024-01-01 8:00:00") + (5-1)*360), 
                        by = 360)),
  data.frame(trapID = 2, 
             x = seq(750, 1750, 250), 
             #x = rep(750, 5),
             y = 1500, 
             time = seq(ymd_hms("2024-01-02 8:00:00"), (ymd_hms("2024-01-02 8:00:00") + (5-1)*360), 
                        by = 360)),
  data.frame(trapID = 3, 
             #x = seq(750, 1750, 250),
             x = rep(1750, 5),
             y = 1000, 
             time = seq(ymd_hms("2024-01-01 8:00:00"), (ymd_hms("2024-01-01 8:00:00") + (5-1)*360), 
                        by = 360)),
  data.frame(trapID = 4, 
             #x = seq(750, 1750, 250), 
             x = rep(1750, 5),
             y = 1500, 
             time = seq(ymd_hms("2024-01-02 8:00:00"), (ymd_hms("2024-01-02 8:00:00") + (5-1)*360), 
                        by = 360))
)

dist_dat <- create_distdat(trapsdf, mesh) #calls function in movingdetectorlikelihood.R to create data object
pop <- sim.popn(D = D_mesh, core = mesh, model2D = "IHP", 
                Ndist = "poisson", buffertype = "rect")
rownames(pop) <- NULL
dist_dat_pop <- create_distdat(trapsdf, pop)

#capthist should be dim(inds, traps)
capthist_full <- sim_capthist(pop, trapsdf, timeincr, lambda0, sigma, D_mesh) #function in movingdetectorlikelihood.R
capthist <- capthist_full[which(rowSums(capthist_full, na.rm = T) > 0),]


########################################################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
library(raster)
library(ggplot2)

#check distance matrix
#trap by mesh by time
trapx <- 2
timex <- 9
plotdatx <- data.frame(x = mesh$x, y = mesh$y, dist = dist_dat$distmat[trapx, ,timex])

ggplot() +
  geom_tile(data = plotdatx, aes(x = x, y = y, fill = dist)) +
  geom_point(data = trapsdf[trapsdf$trapID == dist_dat$traps[trapx, ] &
                            trapsdf$time == dist_dat$times[timex],], aes(x = x,y = y))



#must choose a ind/trap that has a detection
ex_j <- 1
#example individual (from pop and capthist_full) that was detected
ex_i <- min(which(!is.na(capthist_full[,ex_j])))
#calculate hazard of detection at and survival up to observed detection time for all mesh pts
plotdat1 <- data.frame(x = mesh$x,
                       y = mesh$y,
                       hazard = apply(as.array(0:(nrow(mesh)-1)), 1, FUN = function(meshx){haz_cpp(capthist_full[ex_i,ex_j], meshx, (ex_j-1), lambda0, sigma, dist_dat)}),
                       survival = apply(as.array(0:(nrow(mesh)-1)), 1, FUN = function(x){surv_cpp(dist_dat$times[min(which(!is.na(colSums(dist_dat$distmat[ex_j,,], na.rm = F))))],capthist_full[ex_i,ex_j], timeincr, x, (ex_j-1), lambda0, sigma, dist_dat)})
)

ggplot() +
  geom_raster(data = plotdat1,
              mapping = aes(x = x, y = y, fill = hazard)) +
  geom_point(data = trapsdf[trapsdf$trapID == ex_j & trapsdf$time == capthist_full[ex_i,ex_j],], aes(x = x, y = y), shape = "+") +
  geom_point(data = pop[ex_i,], aes(x = x, y = y)) +
  ggtitle(paste("hazard for each mesh at time", capthist_full[ex_i,ex_j], "by trap ", ex_j))

ggplot() +
  geom_raster(data = plotdat1,
              mapping = aes(x = x, y = y, fill = survival)) +
  geom_point(data = trapsdf[trapsdf$trapID == ex_j & trapsdf$time == capthist_full[ex_i,ex_j],], aes(x = x, y = y), shape = "+") +
  geom_point(data = pop[ex_i,], aes(x = x, y = y)) +
  ggtitle(paste("Survival prob for each mesh at time", capthist_full[ex_i,ex_j], "by trap", ex_j))

# #find mesh pt closest to individual i (pop[ex_i,]) to use dist_dat for hazard over time
# dist_imesh <- do.call(rbind, apply(as.array(1:nrow(mesh)), 1, function(x)abs(pop[ex_i,]-mesh[x,])))
# dist_imesh_minx <- dist_imesh[dist_imesh$x == min(dist_imesh$x),]
# dist_imesh_minxy <- dist_imesh_minx[dist_imesh_minx$y == min(dist_imesh_minx$y),]
# closestmesh_m <- which(apply(as.array(1:nrow(dist_imesh)), 1, function(x)all(dist_imesh[x,] == dist_imesh_minxy)))

plotdat2 <- data.frame(times = trapsdf[trapsdf$trapID == ex_j, "time"],
                       haz = sapply(1:length(trapsdf[trapsdf$trapID == ex_j, "time"]),
                                    FUN = function(timej){
                                      t = trapsdf[trapsdf$trapID == ex_j, "time"][timej]
                                      haz_cpp(t, (ex_i - 1), (ex_j - 1), lambda0, sigma, dist_dat_pop) }),
                       surv = sapply(1:length(trapsdf[trapsdf$trapID == ex_j, "time"]),
                                     FUN = function(timej){
                                       t = trapsdf[trapsdf$trapID == ex_j, "time"][timej]
                                       tstart = dist_dat$times[min(which(!is.na(colSums(dist_dat$distmat[ex_j,,], na.rm = F))))]
                                       surv_cpp(tstart, t, timeincr, (ex_i - 1), (ex_j-1), lambda0, sigma, dist_dat_pop)
                                     }),
                       prob1stdett = sapply(1:length(trapsdf[trapsdf$trapID == ex_j, "time"]),
                                            FUN = function(timej){
                                              t = trapsdf[trapsdf$trapID == ex_j, "time"][timej]
                                              tstart = dist_dat$times[min(which(!is.na(colSums(dist_dat$distmat[ex_j,,], na.rm = F))))]
                                              haz_cpp(t, (ex_i - 1), (ex_j - 1), lambda0, sigma, dist_dat_pop) *
                                              surv_cpp(tstart, t, timeincr, (ex_i - 1), (ex_j-1), lambda0, sigma, dist_dat_pop)
                                            }))
ggplot()+
  geom_line(data = plotdat2, aes(x = times, y = haz)) +
  geom_vline(xintercept = capthist_full[ex_i,ex_j], col = "blue")
ggplot()+
  geom_line(data = plotdat2, aes(x = times, y = surv))+
  geom_vline(xintercept = capthist_full[ex_i,ex_j], col = "blue")

ggplot()+
  geom_line(data = plotdat2, aes(x = times, y = log(surv)))+
  geom_vline(xintercept = capthist_full[ex_i,ex_j], col = "blue")
ggplot()+
  geom_line(data = plotdat2, aes(x = times, y = prob1stdett))+
  geom_vline(xintercept = capthist_full[ex_i,ex_j], col = "blue")
ggplot()+
  geom_line(data = plotdat2, aes(x = times, y = log(prob1stdett)))+
  geom_vline(xintercept = capthist_full[ex_i,ex_j], col = "blue")


#trying to minimize
trylambda0s <- parallel::mclapply(as.list(seq(0,.01,.001)),  FUN = function(lambda0_){
  negloglikelihood_cpp(lambda0_, sigma, D_mesh, timeincr, capthist, dist_dat)
}, mc.cores = 4)
plotlambdas <- data.frame(lambda0 = seq(0,.01,.001),
                          llk = unlist(trylambda0s))
ggplot() +
  geom_point(data = plotlambdas, aes(x = lambda0, y = llk)) +
#since lambda is currently on a completely different scale than specified
  geom_vline(xintercept = lambda0)

trysigmas <- parallel::mclapply(as.list(seq(50,500,50)),  FUN = function(sigma_){
  negloglikelihood_cpp(lambda0, sigma_, D_mesh, timeincr, capthist, dist_dat)
}, mc.cores = 4)
plotsigmas <- data.frame(sigma = seq(50,500,50),
                         llk = unlist(trysigmas))
ggplot() +
  geom_point(data = plotsigmas, aes(x = sigma, y = llk))+
  geom_vline(xintercept = sigma)


#check pdot
pdot_m <- rep(NA, nrow(mesh))
for (m in 1:nrow(mesh)){
  surv_eachtrap <- rep(NA, nrow(dist_dat$traps)) 
  for (trapk in 1:nrow(dist_dat$traps)){
    opentimeindx <- which(!is.na(colSums(dist_dat$distmat[trapk,,], na.rm = F)))
    topentime <- dist_dat$times[min(opentimeindx)]
    tclosetime <- dist_dat$times[max(opentimeindx)]
    
    surv_eachtrap[trapk] <-  surv_cpp(topentime, tclosetime, timeincr, m-1, trapk-1, lambda0, sigma, dist_dat)
  }
  survalltraps = product(surv_eachtrap)
  pdot = 1 - survalltraps
  pdot_m[m] <- pdot
}

gettrapwcap <- function(trapk){
  trapwcap <- trapsdf[trapsdf$trapID == trapk & as.numeric(trapsdf$time) %in% as.numeric(capthist[,trapk]),]
  timecapcounts <- as.data.frame(table(as.numeric(capthist[,trapk])))
  if (nrow(timecapcounts) == 0){
    timecapcounts <- data.frame(var1 = NA, Freq = NA)
  }
  colnames(timecapcounts) <- c("timenum", "count")
  trapwcap$numdets <- timecapcounts[timecapcounts$time == as.numeric(trapwcap$time),"count"]
  return(trapwcap)
}
trapwcap <- do.call(rbind,
                    apply(as.array(1:nrow(dist_dat$traps)), 1, 
                          gettrapwcap))
popcap <- cbind(pop, apply(capthist_full, 1, FUN = function(x){length(which(!is.na(x)))}))
colnames(popcap)[3] <- "captures"
ggplot() +
  geom_point(data = cbind(mesh, pdot_m), aes(x = x, y = y, col = pdot_m)) +
  geom_point(data = popcap, aes(x = x, y = y, size = as.factor(captures)), shape = 5) +
  geom_point(data = trapwcap, aes(x = x, y = y#, col = numdets
                                  ), size = 5) +
  scale_color_viridis_c() +
  
  #May 16 so this all works sometimes, but not all the time. Sometimes, lambda0
  #needs to be much smaller to minimize the neg llk. I think this means
  #that in my simulation, I sometimes miss detections?
  #but also sometimes it looks like a detection happens when the hazard is small
  #and survival is high. Which is problematic
  theme_classic()

library(secr)
library(lubridate)

Rcpp::sourceCpp("functions.cpp")

######data setup
lambda0 = .9
sigma = 300
timeincr = 6*60 #specifies time increments for integration AND for distance matrix 
meshgrid <- expand.grid(x = seq(000, 6000, 300), y = seq(000,6000, 300))
mesh <- make.mask(meshgrid, buffer = 0, spacing = 500)
D_mesh <- rep(.1, nrow(mesh))

#traps
#specify K moving traps
#dataframe with column for trapid, x, y, and time
trapsdf <- rbind(
  data.frame(trapID = 1, x = seq(1500, 4500, 50), y = 1500, 
             time = seq(ymd_hms("2024-01-01 8:00:00"), ymd_hms("2024-01-01 14:00:00"), 
                        by = 360)),
  data.frame(trapID = 2, x = seq(1500, 4500, 50), y = 3500, 
             time = seq(ymd_hms("2024-01-02 8:00:00"), ymd_hms("2024-01-02 14:00:00"), 
                        by = 360)),
  data.frame(trapID = 3, x = seq(1500, 4500, 50), y = 2500, 
             time = seq(ymd_hms("2024-01-03 8:00:00"), ymd_hms("2024-01-03 14:00:00"), 
                        by = 360)),
  data.frame(trapID = 4, x = seq(1500, 4500, 50), y = 4500, 
             time = seq(ymd_hms("2024-01-01 9:00:00"), ymd_hms("2024-01-01 15:00:00"), 
                        by = 360))
)

dist_dat <- create_distdat(trapsdf, mesh) #calls function in movingdetectorlikelihood.R to create data object
pop <- sim.popn(D = D_mesh, core = mesh, model2D = "IHP", 
                Ndist = "poisson", buffertype = "rect")
rownames(pop) <- NULL

capthist <- sim_capthist(pop, trapsdf, timeincr, lambda0, sigma, D_mesh) #function in movingdetectorlikelihood.R



########################################################################################################################################
######################################################################################################
######################################################################################################
######################################################################################################
library(raster)
library(ggplot2)
plotdat1 <- cbind(mesh, apply(as.array(0:(nrow(mesh)-1)), 1, FUN = function(x){haz_cpp(capthist[37,1], x, 0, lambda0, sigma, dist_dat)}))
colnames(plotdat1)[3] <- "hazard"
plotdat1$survival <- apply(as.array(0:(nrow(mesh)-1)), 1, FUN = function(x){surv_cpp(dist_dat$times[1],capthist[37,1], timeincr/100,x, 0, lambda0, sigma, dist_dat)})
ggplot() +
  geom_raster(data = plotdat1,
              mapping = aes(x = x, y = y, fill = log(hazard))) +
  geom_point(data = pop[37,], aes(x = x, y = y)) +
  ggtitle(paste("log hazard for each mesh at time", capthist[37,1], "by trap 1"))

ggplot() +
  geom_raster(data = plotdat1,
              mapping = aes(x = x, y = y, fill = survival)) +
  geom_point(data = pop[37,], aes(x = x, y = y)) + 
  ggtitle(paste("Survival prob for each mesh at time", capthist[37,1], "by trap 1"))

plotdat2 <- data.frame(times = trapsdf[trapsdf$trapID == 1, "time"],
                       haz = sapply(1:length(trapsdf[trapsdf$trapID == 1, "time"]),
                                    FUN = function(i){
                                      t = trapsdf[trapsdf$trapID == 1, "time"][i]
                                      haz_cpp(t, 16, 0, lambda0, sigma, dist_dat) }),
                       surv = sapply(1:length(trapsdf[trapsdf$trapID == 1, "time"]),
                                     FUN = function(i){
                                       t = trapsdf[trapsdf$trapID == 1, "time"][i]
                                         surv_cpp(dist_dat$times[1], t, timeincr, 16, 0, lambda0, sigma, dist_dat)
                                     }),
                       prob1stdett = sapply(1:length(trapsdf[trapsdf$trapID == 1, "time"]),
                                            FUN = function(i){
                                              t = trapsdf[trapsdf$trapID == 1, "time"][i]
                                              haz_cpp(t, 16, 0, lambda0, sigma, dist_dat) *
                                                surv_cpp(dist_dat$times[1], t, timeincr, 16, 0, lambda0, sigma, dist_dat)
                                            }))
ggplot()+
  geom_line(data = plotdat2, aes(x = times, y = haz))
ggplot()+
  geom_line(data = plotdat2, aes(x = times, y = log(surv)))
ggplot()+
  geom_line(data = plotdat2, aes(x = times, y = prob1stdett))
ggplot()+
  geom_line(data = plotdat2, aes(x = times, y = log(prob1stdett)))

trylambda0s <- parallel::mclapply(as.list(seq(0,1,.1)),  FUN = function(lambda0_){
  negloglikelihood_cpp(lambda0_, sigma, D_mesh, timeincr, capthist, dist_dat)
}, mc.cores = 4)
trysigmas <- parallel::mclapply(as.list(seq(50,500,50)),  FUN = function(sigma_){
  negloglikelihood_cpp(lambda0, sigma_, D_mesh, timeincr, capthist, dist_dat)
}, mc.cores = 4)
plotlambdas <- data.frame(lambda0 = seq(0,1,.1),
                          llk = unlist(trylambda0s))
plotsigmas <- data.frame(sigma = seq(50,500,50),
                         llk = unlist(trysigmas))
ggplot() +
  geom_point(data = plotlambdas, aes(x = lambda0, y = llk)) +
  geom_vline(xintercept = lambda0)
ggplot() +
  geom_point(data = plotsigmas, aes(x = sigma, y = llk))+
  geom_vline(xintercept = sigma)


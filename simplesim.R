library(secr)
library(ggplot2)
library(tidyr)
source("movingdetectorlikelihood.R")
lambda0 = .9
sigma = 300
timeincrement = 6*60*10
meshgrid <- expand.grid(x = seq(000, 6000, 300), y = seq(000,6000, 300))
mesh <- make.mask(meshgrid, buffer = 0, spacing = 500)
D_mesh <- rep(.1, nrow(mesh))


#specify K moving traps
#traps should be dataframe with column for trapid, x, y, and time
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

#distance data object (created from traps and mesh)
dist_dat <- create_distdat(trapsdf, mesh)

#simulate capture history
pop <- sim.popn(D = D_mesh, core = mesh, model2D = "IHP", 
                Ndist = "poisson", buffertype = "rect")
rownames(pop) <- NULL

capthist <- sim_capthist(pop, trapsdf, timeincrement, lambda0, sigma, D_mesh)
gettrapwcap <- function(trapk){
  trapwcap <- trapsdf[trapsdf$trapID == trapk & as.numeric(trapsdf$time) %in% as.numeric(capthist[,trapk]),]
  timecapcounts <- as.data.frame(table(as.numeric(capthist[,trapk])))
  colnames(timecapcounts) <- c("timenum", "count")
  trapwcap$numdets <- timecapcounts[timecapcounts$time == as.numeric(trapwcap$time),"count"]
  return(trapwcap)
}
trapwcap <- rbind(gettrapwcap(1), 
                  gettrapwcap(2),
                  gettrapwcap(3),
                  gettrapwcap(4))
popcap <- cbind(pop, apply(capthist, 1, FUN = function(x){length(which(!is.na(x)))}))
colnames(popcap)[3] <- "captures"

ggplot() +
  geom_point(data = mesh, aes(x = x, y = y)) +
  geom_point(data = popcap, aes(x = x, y = y, size = as.factor(captures)), shape = 5) +
  geom_point(data = trapwcap, aes(x = x, y = y, col = numdets), size = 5) +
  scale_color_viridis_c() +
  theme_classic()


llk <- negloglikelihood(lambda0, sigma, D_mesh, timeincrement, capthist, dist_dat)

trylambda0s <- parallel::mclapply(as.list(seq(0,1,.1)),  FUN = function(lambda0_){
  negloglikelihood(lambda0_, sigma, D_mesh, timeincrement, capthist, dist_dat)
}, mc.cores = 4)
trysigmas <- parallel::mclapply(as.list(seq(50,500,50)),  FUN = function(sigma_){
  negloglikelihood(lambda0, sigma_, D_mesh, timeincrement, capthist, dist_dat)
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


#inspect
#example individual 1, trap 1
testk = 3
testx = 54
xi <- mesh[testx,]
dxi <- dist_dat$distmat[testk,testx,]
hazxi <- apply(as.array(dist_dat$times), 1, FUN = function(t){ 
  haz(t, testx, testk, lambda0, sigma, dist_dat)})
survxi <- apply(as.array(1:length(dist_dat$times)), 1, FUN = function(et){
  endt <- dist_dat$times[et]
  opentimeindx <- which(!is.na(colSums(dist_dat$distmat[testk,,], na.rm = F)))
  starttime <- dist_dat$times[min(opentimeindx)]
  surv(timestart = starttime, timeend = endt, timeincr = timeincrement, testx, testk, lambda0, sigma, dist_dat)
})
cumhaz <- -log(survxi)
detxi <- survxi * hazxi
plotdat <- data.frame(time = dist_dat$times, 
                      distance = dxi/2000*3,
                      hazard = hazxi,
                      cumulativehazard = cumhaz,
                      survival = survxi,
                      detection = detxi)
plotdat_long <- pivot_longer(plotdat, cols = c("hazard", "survival", "detection", "distance","cumulativehazard"))
ggplot() +
  geom_point(data = mesh, mapping = aes(x = x, y = y)) +
  geom_line(data = traps[traps$trapID == testk,], mapping = aes(x =x, y =y)) +
  geom_point(data = mesh[testx,], mapping = aes(x = x, y = y), col = "red")

ggplot() +
  geom_line(data = plotdat_long, mapping = aes(x = time, y = value, col = name)) +
  #geom_line(data = plotdat_long, mapping = aes(x = time, y = distance/2000*3)) +
  xlim(min(plotdat$time[!is.na(plotdat$distance)]), max(plotdat$time[!is.na(plotdat$distance)])) +
  scale_y_continuous(
    name = "",
    sec.axis = sec_axis(~.*2000/3, name="Distance (m)")
  )

ggplot()+
  geom_line(data = plotdat, mapping = aes(x = time, y = distance*2000/3))+
  xlim(min(plotdat$time[!is.na(plotdat$distance)]), max(plotdat$time[!is.na(plotdat$distance)]))

ggplot()+
  geom_line(data = plotdat, mapping = aes(x = time, y = hazard))+
  xlim(min(plotdat$time[!is.na(plotdat$distance)]), max(plotdat$time[!is.na(plotdat$distance)]))

ggplot()+
  geom_line(data = plotdat, mapping = aes(x = time, y = cumulativehazard))+
  xlim(min(plotdat$time[!is.na(plotdat$distance)]), max(plotdat$time[!is.na(plotdat$distance)]))

ggplot()+
  geom_line(data = plotdat, mapping = aes(x = time, y = survival))+
  xlim(min(plotdat$time[!is.na(plotdat$distance)]), max(plotdat$time[!is.na(plotdat$distance)]))

ggplot()+
  geom_line(data = plotdat, mapping = aes(x = time, y = detection))+
  xlim(min(plotdat$time[!is.na(plotdat$distance)]), max(plotdat$time[!is.na(plotdat$distance)]))


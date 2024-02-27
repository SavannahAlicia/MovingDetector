library(secr)
source("movingdetectorlikelihood.R")
lambda0 = .9
sigma = 300
timeincrement = 6*60*10
meshgrid <- expand.grid(x = seq(000, 6000, 300), y = seq(000,6000, 300))
mesh <- make.mask(meshgrid, buffer = 0, spacing = 500)
D_mesh <- rep(.1, nrow(mesh))

#simulate N home range centers
# simulate population
pop <- sim.popn(D = D_mesh, core = mesh, model2D = "IHP", 
                Ndist = "poisson", buffertype = "rect")
rownames(pop) <- NULL

#specify K moving traps
#traps should be dataframe with column for trapid, x, y, and time
traps <- rbind(
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
traptomesh <- apply(as.array(1:nrow(mesh)), 1, FUN = 
                      function(meshrow){
                        apply(as.array(1:nrow(traps)), 
                              1, FUN = 
                                function(traprow){
                                  dist(rbind(traps[traprow, c("x", "y")],
                                             mesh[meshrow,c("x","y")]), method = "euclidean")
                                  })})
dist_dat <- list(traps = data.frame(traps = unique(traps$trapID)),
                 mesh = as.data.frame(mesh),
                 times = sort(unique(traps$time)))
dist_dat$distmat <- sapply(as.array(dist_dat$times),  FUN = function(timet){ 
  apply(as.array(1:nrow(dist_dat$mesh)), 1, FUN = function(meshcol){
    apply(as.array(1:nrow(dist_dat$traps)), 1, FUN = function(trapid){
      dist <- traptomesh[traps$trapID == trapid & traps$time == timet,meshcol]
      if (length(dist) == 0){
        dist <- NA
      }
      return(dist)
      })
    })
  }, simplify = "array")

#simulate capture history
#i think i can just use a new dist_dat that uses pop as hrcs? instead of mesh?
#for each trap, for each ind
#check if detected while trap was open 
#distance data object (created from traps and mesh)


  traptopop <- apply(as.array(1:nrow(pop)), 1, FUN = 
                       function(meshrow){
                         apply(as.array(1:nrow(traps)), 
                               1, FUN = 
                                 function(traprow){
                                   dist(rbind(traps[traprow, c("x", "y")],
                                              pop[meshrow,c("x","y")]), method = "euclidean")
                                 })})
  
  dist_dat_pop <- list(traps = data.frame(traps = unique(traps$trapID)),
                       mesh = as.data.frame(pop),
                       times = sort(unique(traps$time)))
  dist_dat_pop$distmat <- sapply(as.array(dist_dat_pop$times),  FUN = function(timet){ 
    apply(as.array(1:nrow(dist_dat_pop$mesh)), 1, FUN = function(meshcol){
      apply(as.array(1:nrow(dist_dat_pop$traps)), 1, FUN = function(trapid){
        dist <- traptopop[traps$trapID == trapid & traps$time == timet,meshcol]
        if (length(dist) == 0){
          dist <- NA
        }
        return(dist)
      })
    })
  }, simplify = "array")
  capthist <- lapply(as.list(1:nrow(dist_dat_pop$mesh)), 
                     FUN = function(hrcx){
                       lapply(as.list(1:nrow(dist_dat_pop$traps)),
                              FUN = function(trapk){
                                #can do a dominating process of rate lambda0*T and thin
                                opentimeindx <- which(!is.na(colSums(dist_dat_pop$distmat[trapk,,], na.rm = F)))
                                topentime <- dist_dat_pop$times[min(opentimeindx)]
                                tclosetime <- dist_dat_pop$times[max(opentimeindx)]
                                probseenxk <- 1 - surv(topentime, tclosetime, timeincrement, hrcx, trapk, lambda0, sigma, dist_dat_pop)
                                seenxk_bool <- rbinom(1,1, probseenxk)
                                if (seenxk_bool){
                                  timesopen <- seq(topentime, tclosetime, timeincrement)
                                  detattime <- apply(as.array(1:length(timesopen)), 1, FUN = function(t){
                                    notseenuntil <- surv(topentime, timesopen[t], timeincrement, hrcx, trapk, lambda0, sigma, dist_dat_pop)
                                    seenat <- haz(timesopen[t], hrcx, trapk, lambda0, sigma, dist_dat_pop)
                                    probdettime <- notseenuntil * seenat
                                  })
                                  capik <- timesopen[sample(length(detattime), 1, prob = detattime/sum(detattime))]
                                } else {
                                  capik <- ymd_hms(NA)
                                }
                                return(ymd_hms(capik))
                              })
                     }) 

  capthist_array <- t(structure(array(unlist(capthist), dim = c(4,nrow(dist_dat_pop$mesh))), class = c("POSIXct", "POSIXt")))
gettrapwcap <- function(trapk){
  trapwcap <- traps[traps$trapID == trapk & as.numeric(traps$time) %in% as.numeric(capthist_array[,trapk]),]
  timecapcounts <- as.data.frame(table(as.numeric(capthist_array[,trapk])))
  colnames(timecapcounts) <- c("timenum", "count")
  trapwcap$numdets <- timecapcounts[timecapcounts$time == as.numeric(trapwcap$time),"count"]
  return(trapwcap)
}
trapwcap <- rbind(gettrapwcap(1), 
                  gettrapwcap(2),
                  gettrapwcap(3),
                  gettrapwcap(4))
popcap <- cbind(pop, apply(capthist_array, 1, FUN = function(x){length(which(!is.na(x)))}))
colnames(popcap)[3] <- "captures"

ggplot() +
  geom_point(data = mesh, aes(x = x, y = y)) +
  geom_point(data = popcap, aes(x = x, y = y, size = as.factor(captures)), shape = 5) +
  geom_point(data = trapwcap, aes(x = x, y = y, col = numdets), size = 5) +
  scale_color_viridis_c() +
  theme_classic()


llk <- negloglikelihood(lambda0, sigma, D_mesh, timeincrement, capthist_array, dist_dat)

trylambda0s <- parallel::mclapply(as.list(seq(0,1,.1)),  FUN = function(lambda0_){
  negloglikelihood(lambda0_, sigma, D_mesh, timeincrement, capthist_array, dist_dat)
}, mc.cores = 4)
trysigmas <- parallel::mclapply(as.list(seq(50,500,50)),  FUN = function(sigma_){
  negloglikelihood(lambda0, sigma_, D_mesh, timeincrement, capthist_array, dist_dat)
}, mc.cores = 4)
plotlambdas <- data.frame(lambda0 = seq(0,1,.1),
                          llk = unlist(trylambda0s))
plotsigmas <- data.frame(sigma = seq(50,500,50),
                          llk = unlist(trysigmas))
ggplot() +
  geom_point(data = plotlambdas, aes(x = lambda0, y = llk))
ggplot() +
  geom_point(data = plotsigmas, aes(x = sigma, y = llk))


#inspect
#example individual 1, trap 1
testk = 3
testx = 30
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


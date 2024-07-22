#compare moving and stationary detector fits
library(secr)
library(lubridate)
library(parallel)
library(ggplot2)
setwd("~/Documents/UniStAndrews/MovingDetector")
Rcpp::sourceCpp("functions.cpp")
source("movingdetectorlikelihood.R")

set.seed(12345)
nsims = 100

#----------------------------------data setup-----------------------------------
lambda0 = .7
sigma = 300
meshgrid <- expand.grid(x = seq(0, 2600, 100), y = seq(0,2600, 100))
mesh <- make.mask(meshgrid, buffer = 800, spacing = 250)
D_mesh <- rep(.4, nrow(mesh))
beta1 <- -(1/40000)
beta2 <- -1250
D_mesh_q <- exp(beta1*(mesh$x + beta2)^2)
truepar <- data.frame(name = c("lambda0", "sigma", "beta1", "beta2"),
                      value = c(lambda0, sigma, beta1, beta2))
##REMOVE THIS AFTER DEBUG
#D_mesh <- D_mesh_q
#Dmod = "~x^2"
##
ggplot() + 
  geom_line(data = data.frame(x = mesh$x, y = D_mesh_q), mapping = aes(x = x, y = y)) + 
  geom_line(data = data.frame(x = mesh$x, y = D_mesh_q), mapping = aes(x = x, y = y)) + 
  geom_vline(xintercept = c(750, 1750))

timeincr = 6*60 #specifies time increments for integration AND for distance matrix 

tracksteps = 10

#-----------------------------detectors moving left to right--------------------

trapsdf <- rbind(
  data.frame(trapID = 1, 
             occ = 1,
             x = seq(from = 750, to = 1750, length.out = tracksteps),
             y = 750, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps-1)*360), 
                        by = 360)),
  data.frame(trapID = 2, 
             occ = 1,
             x = seq(from = 750, to = 1750, length.out = tracksteps),
             y = 1050, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps-1)*360), 
                        by = 360)),
  data.frame(trapID = 3, 
             occ = 1,
             x = seq(from = 750, to = 1750, length.out = tracksteps),
             y = 1375, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps-1)*360), 
                        by = 360)),
  data.frame(trapID = 4, 
             occ = 1,
             x = seq(from = 750, to = 1750, length.out = tracksteps),
             y = 1750, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps-1)*360), 
                        by = 360)),
  data.frame(trapID = 5, 
             occ = 1,
             x = seq(from = 750, to = 1750, length.out = tracksteps),
             y = 1750, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps-1)*360), 
                        by = 360))
  
)

dist_dat <- create_distdat(trapsdf, mesh) #calls function in movingdetectorlikelihood.R to create data object
ggplot() +
  geom_tile(data = data.frame(x = mesh$x, y = mesh$y, D = D_mesh_q), aes(x= x, y = y, fill = D)) +
  geom_point(data = mesh, aes(x = x, y = y), col = "grey") +
  geom_point(data = trapsdf, aes(x = x, y = y, col = as.factor(trapID)), shape = "+", size = 5) +
  theme_classic()

###-----------simulate data for moving and stationary detectors-----------------
###-------and fit model with stationary and moving detector---------------------
sim_fit <- function(trapsdf, dist_dat, lambda0, sigma, D_mesh, timeincr, mesh, Dmod = "~1"){
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
  if (Dmod == "~1"){
    stat_nll <- function(v){
      lambda0_ <- exp(v[1])
      sigma_ <- exp(v[2])
      D_mesh_ <- rep(exp(v[3]), nrow(mesh))
      out <- negloglikelihood_stationary_cpp(lambda0_, sigma_, D_mesh_, secrcapthist, usage(secrtraps), distmat, dist_dat)
      return(out)
    }
    #moving detector likelihood
    nll <- function(v){
      lambda0_ <- exp(v[1])
      sigma_ <- exp(v[2])
      D_mesh_ <- rep(exp(v[3]), nrow(mesh))
      out <- negloglikelihood_cpp(lambda0_, sigma_, D_mesh_, timeincr, capthist, dist_dat)
      return(out)
    }
    start <- c( log(.4), log(300), log(.4))
    
  }else if(Dmod == "~x^2"){
    #quadratic density function
    stat_nll <- function(v){
      lambda0_ <- exp(v[1])
      sigma_ <- exp(v[2])
      D_mesh_ <- exp(v[3]*(mesh$x + v[4])^2)
      out <- negloglikelihood_stationary_cpp(lambda0_, sigma_, D_mesh_, secrcapthist, usage(secrtraps), distmat, dist_dat)
      return(out)
    }
    #moving detector likelihood
    nll <- function(v){
      lambda0_ <- exp(v[1])
      sigma_ <- exp(v[2])
      D_mesh_ <- exp(v[3]*(mesh$x + v[4])^2)#exp(beta1*(mesh$x + beta2)^2)
      out <- negloglikelihood_cpp(lambda0_, sigma_, D_mesh_, timeincr, capthist, dist_dat)
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
              movdet_time = fit.time.md)
  return(out)
}

#test one of each
#fit1 <- sim_fit(trapsdf, dist_dat, lambda0, sigma, D_mesh, timeincr, mesh)
#fit2 <- sim_fit(trapsdf, dist_dat, lambda0, sigma, D_mesh_q, timeincr, mesh, Dmod = "~x^2")

start.time.all_q <- Sys.time()
all_sim_fits_q <- mclapply(X = as.list(1:nsims),
                           FUN = function(sim){
                             return(sim_fit(trapsdf, dist_dat, lambda0, sigma, D_mesh_q, timeincr, mesh,Dmod = "~x^2"))
                           },
                           mc.cores = 6
)
tot.time.all_q <- difftime(Sys.time(), start.time.all_q, units = "secs")



start.time.all <- Sys.time()
all_sim_fits <- mclapply(X = as.list(1:nsims),
                   FUN = function(sim){
                     return(sim_fit(trapsdf, dist_dat, lambda0, sigma, D_mesh, timeincr, mesh))
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
                        y = rep(rep(mesh$y[1], 81),3),
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
                        quantile = c(rep("mean", 81), rep("2.5%",81), rep("97.5%",81))
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


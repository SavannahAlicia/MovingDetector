library(secr)
source("movingdetectorlikelihood.R")
lambda0 = .9
sigma = 300
timeincrement = 6*60*10
meshgrid <- expand.grid(x = seq(000, 6000, 300), y = seq(000,6000, 300))
mesh <- make.mask(meshgrid, buffer = 0, spacing = 500)
D_mesh <- rep(.1, nrow(mesh))

#traps
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

dist_dat <- create_distdat(traps, mesh)
#simulate population
pop <- sim.popn(D = D_mesh, core = mesh, model2D = "IHP", 
                Ndist = "poisson", buffertype = "rect")
rownames(pop) <- NULL

#simulate capture histories
ch_ls <- lapply(list(1:1), FUN = function(x){
  sim_capthist(pop, traps, timeincrement, lambda0, sigma, D_mesh)})

system.time(foo <- optim(init, nllk_pars, timeincr = timeincrement, capthist = capthist_array, dist_dat = dist_dat))

#estimate parameters based on capture histories
startsystime <- Sys.time()
estimates <- lapply(ch_ls, FUN = function(capthist_arrays){
  optim(init, nllk_pars, timeincr = timeincrement, capthist = capthist_arrays, dist_dat = dist_dat, control = list(trace = 1))
})
runtime <- Sys.time() - startsystime
                    

init <- log(c(lambda0 = 0.8, sigma = 500, D = .001))
nlm(likelihood, init, print.level = 2)
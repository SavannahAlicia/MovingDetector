#compare fits between my script and Abinand's
source("compare_moving_stat_2D/data_formatting/00_functions_and_parameters.R")
source("compare_moving_stat_2D/data_formatting/01_data_setup.R")
source("compare_moving_stat_2D/data_formatting/01.5_visualize.R")
source("compare_moving_stat_2D/data_formatting/02_simulate_and_fit.R")
source("inst/abifunct.R")

#--------------Try calculating one llk value------------------------------------
simdch <- simulate_popandcapthist(tracksdf, D_mesh_v, lambda0, sigma,
                                  mesh, traps, trapspacing)
ch <- simdch$capthist

#me
myllk <- calc_nll(dist_trapmesh,
                  useall,
                  lambda0, 
                  sigma, 
                  D_mesh = D_mesh_v, 
                  beta1, 
                  beta2,
                  N,
                  hazdenom, 
                  mesh, 
                  capthistout = simdch,
                  Dmod = "~x^2",
                  meshspacing,
                  meanstepsize = mean(tracksdf$inc[tracksdf$inc != 0]))
#Abinand
eta = beta1*((mesh$x/meshspacing + beta2)^2 )
Z = sum(exp(eta)) * meshspacing
b0 = log(N/Z) 
startparams <- c(b0,beta2,beta1,log(lambda0),log(sigma))
#make trapSteps
trapSteps <- data.frame(matrix(nrow = nrow(tracksdf), ncol = 0))
trapSteps$x <- traps$x[tracksdf$trapno]
trapSteps$y <- traps$y[tracksdf$trapno]
trapSteps$xStep <- tracksdf$midx
trapSteps$yStep <- tracksdf$midy
trapSteps$effort <- tracksdf$inc/100
trapSteps$transect <- tracksdf$transect
trapSteps$TrapID <- paste0("Trap",tracksdf$trapno)
trapSteps <- trapSteps[which(tracksdf$rep == 1 & tracksdf$inc != 0),]
trapSteps$trapOrder <- rep((1:(trap_n_horiz * nsteps_pertrap)),2)
trapSteps$stepOrder <- rep(c(1:10), (nrow(trapSteps)/10))
trapSteps$StepID <- paste0("Step", rownames(trapSteps))

#make secr capthist
mask <- mesh
covariates(mask) <- data.frame(cov = D_mesh_v,
                               x = mesh$x,
                               x2 = mesh$x^2/meshspacing^2)
trapscr <- read.traps(data = data.frame(TrapID = paste0("Trap", 1:nrow(traps)),
                                        x = traps$x,
                                        y = traps$y),
                      detector = "proximity")
rownames(trapscr) <- paste0("Trap", 1:nrow(traps))
detinch <- which(ch == 1, arr.ind = T) #ikj
detinch <- data.frame(detinch) %>%
  arrange(as.character(dim1), dim3, dim2)
trap_n_vert <- ntrapsish/trap_n_horiz
CH <- make.capthist(captures = data.frame(
  Session = 1,
  ID = paste0("ind_", detinch[,1]),
  Occasion = ((detinch[,2] - 1) %/% trap_n_vert) + 1,
  TrapID = paste0("Trap",detinch[,3])
),
traps = trapscr,
fmt = "trapID"
)

stepOrder <- apply(as.array(1:nrow(detinch)), 1, function(r){
  induse[detinch[r,1], detinch[r,3], detinch[r,2]]/10
})

datobj <- prep_dat_for_lik(trapSteps, CH)
abillk <- lik_opt(startparams,
                  capthistscr= CH,
                  mask = mesh,
                  distmat = ,
                  nOccasion = nreps,
                  dets,
                  noDets,
                  effort,
                  model = NULL
)

#Try one fit

myfit <- fit_capthist(dist_trapmesh,
                      useall,
                      lambda0, 
                      sigma, 
                      D_mesh, 
                      beta1, 
                      beta2,
                      N,
                      hazdenom, 
                      mesh, 
                      capthistout,
                      Dmod,
                      meshspacing,
                      meanstepsize,
                      fitstat = FALSE
)

abifit <- nlm(lik_opt, startparams, print = 1)



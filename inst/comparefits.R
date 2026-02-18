#compare fits between my script and Abinand's
source("compare_moving_stat_2D/data_formatting/00_functions_and_parameters.R")
source("compare_moving_stat_2D/data_formatting/01_data_setup.R")
source("compare_moving_stat_2D/data_formatting/02_simulate_and_fit.R")
source("inst/BananSim.R")
trap_n_vert = round(ntrapsish/trap_n_horiz)
tracksteplength = round(trapspacing/nsteps_pertrap)

#create input objects for Abinand's
survObj <- simulateScrTrapsMask(
  nxTraps = trap_n_horiz,
  nyTraps = trap_n_vert,
  trapSpacing = trapspacing,
  nSteps = round(trapspacing/tracksteplength),
  maskSpacing = meshspacing,
  sigma = sigma, 
  N = N,
  b1 = beta1,
  b2 = beta2
)

trapSteps <- survObj$trapSteps
mask <- survObj$mask

surv_obj <- setup_data(sigma,
                         N,
                         beta1,
                         beta2,
                         ntrapsish,
                         trackxmin,
                         trapspacing,
                         meshspacing,
                         trap_n_horiz,
                         nsteps_pertrap,
                         occreps
  )
tracksdf <- surv_obj$tracksdf
traps <- surv_obj$traps
mesh <- surv_obj$mesh

ch_out <- simulate_popandcapthist(traps,
                                    tracksdf, 
                                    lambda0,
                                    sigma,
                                    D_mesh = surv_obj$D_mesh_v,
                                    mesh,
                                    meshspacing = surv_obj$meshspacing,
                                    hazdenom = 1)

ch <- ch_out$capthist
induse <- ch_out$induse






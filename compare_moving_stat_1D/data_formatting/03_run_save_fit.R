#1 dimensional
set.seed(4994)
start.time.all_q <- Sys.time()
all_sim_fits_q <- mclapply(X = as.list(1:nsims),
                           FUN = function(sim){
                             set.seed(sim)
                             return(sim_fit(tracksdf = tracksdf, 
                                            traps = traps,
                                            trapcells = trapcells,
                                            dist_trapmesh = dist_trapmesh_lin,
                                            tracksmeshdistmat = tracksmeshdistmat_lin,
                                            useall = useall,
                                            lambda0 = lambda0,
                                            sigma = sigma,
                                            D_mesh = D_meshlin_q,
                                            beta1 = beta1, 
                                            beta2 = beta2, 
                                            hazdenom = hazdenom, 
                                            mesh = meshlin, 
                                            meshunit = meshunit_lin,
                                            Dmod = "~x^2",
                                            linear = T))
                           },
                           mc.cores = 6
)
tot.time.all_q <- difftime(Sys.time(), start.time.all_q, units = "secs")

saveRDS(all_sim_fits_q, file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/simulation_results/1D/all_sim_fitsq.Rds")

start.time.all1 <- Sys.time()
all_sim_fits1 <- mclapply(X = as.list(1:nsims),
                          FUN = function(sim){
                            return(sim_fit(tracksdf = tracksdf, 
                                           traps = traps,
                                           trapcells = trapcells,
                                           dist_trapmesh = dist_trapmesh_lin,
                                           tracksmeshdistmat = tracksmeshdistmat_lin,
                                           useall = useall,
                                           lambda0 = lambda0,
                                           sigma = sigma,
                                           D_mesh = D_meshlin,
                                           beta1 = beta1, 
                                           beta2 = beta2, 
                                           hazdenom = hazdenom, 
                                           mesh = meshlin, 
                                           meshunit = meshunit_lin,
                                           Dmod = "~1",
                                           linear = T))
                          },
                          mc.cores = 6
)
tot.time.all1 <- difftime(Sys.time(), start.time.all1, units = "secs")

saveRDS(all_sim_fits1, file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/simulation_results/1D/all_sim_fits1.Rds")

#-----------------------------------------------------------------------------
#2 dimensional
start.time.all_q2D <- Sys.time()
all_sim_fits_q2D <- mclapply(X = as.list(1:nsims),
                           FUN = function(sim){
                             return(sim_fit(tracksdf = tracksdf, 
                                            traps = traps,
                                            trapcells = trapcells,
                                            dist_trapmesh = dist_trapmesh_2D,
                                            tracksmeshdistmat = tracksmeshdistmat_2D,
                                            useall = useall,
                                            lambda0 = lambda0,
                                            sigma = sigma,
                                            D_mesh = D_mesh2D_q,
                                            beta1 = beta1, 
                                            beta2 = beta2, 
                                            hazdenom = hazdenom, 
                                            mesh = mesh2D, 
                                            meshunit = meshunit_2D,
                                            Dmod = "~x^2",
                                            linear = F))
                           },
                           mc.cores = 6
)
tot.time.all_q2D <- difftime(Sys.time(), start.time.all_q2D, units = "secs")

saveRDS(all_sim_fits_q2D, file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/simulation_results/2D/all_sim_fitsq.Rds")

start.time.all12D <- Sys.time()
all_sim_fits12D <- mclapply(X = as.list(1:nsims),
                          FUN = function(sim){
                            return(sim_fit(tracksdf = tracksdf, 
                                           traps = traps,
                                           trapcells = trapcells,
                                           dist_trapmesh = dist_trapmesh_2D,
                                           tracksmeshdistmat = tracksmeshdistmat_2D,
                                           useall = useall,
                                           lambda0 = lambda0,
                                           sigma = sigma,
                                           D_mesh = D_mesh2D,
                                           beta1 = beta1, 
                                           beta2 = beta2, 
                                           hazdenom = hazdenom, 
                                           mesh = mesh2D, 
                                           meshunit = meshunit_2D,
                                           Dmod = "~1",
                                           linear = F))
                          },
                          mc.cores = 6
)
tot.time.all12D <- difftime(Sys.time(), start.time.all12D, units = "secs")

saveRDS(all_sim_fits12D, file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/simulation_results/2D/all_sim_fits1.Rds")

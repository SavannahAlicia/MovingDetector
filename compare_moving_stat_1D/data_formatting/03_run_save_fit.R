start.time.all_q <- Sys.time()
all_sim_fits_q <- mclapply(X = as.list(1:nsims),
                           FUN = function(sim){
                             return(sim_fit(tracksdf, 
                                            traps,
                                            trapcells,
                                            dist_trapmesh,
                                            useall,
                                            lambda0, sigma, D_mesh_q,
                                            beta1, beta2, beta3,
                                            hazdenom, 
                                            mesh, 
                                            Dmod = "~x^2"))
                           },
                           mc.cores = 6
)
tot.time.all_q <- difftime(Sys.time(), start.time.all_q, units = "secs")

saveRDS(all_sim_fits_q, file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/simulation_results/all_sim_fitsq.Rds")

start.time.all1 <- Sys.time()
all_sim_fits1 <- mclapply(X = as.list(1:nsims),
                          FUN = function(sim){
                            return(sim_fit(tracksdf, 
                                           traps,
                                           trapcells,
                                           dist_trapmesh,
                                           useall,
                                           lambda0, sigma, D_mesh,
                                           beta1, beta2, beta3,
                                           hazdenom, 
                                           mesh, 
                                           Dmod = "~1"))
                          },
                          mc.cores = 6
)
tot.time.all1 <- difftime(Sys.time(), start.time.all1, units = "secs")

saveRDS(all_sim_fits1, file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/simulation_results/all_sim_fits1.Rds")


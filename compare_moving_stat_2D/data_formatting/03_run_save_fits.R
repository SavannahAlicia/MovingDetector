#fit multiple models and save output

start.time.all_q <- Sys.time()
all_sim_fits_q <- mclapply(X = as.list(1:nsims),
                           FUN = function(sim){
                             return(sim_fit(tracksdf, 
                                            traps,
                                            trapcells,
                                            dist_trapmesh,
                                            useall,
                                            lambda0, sigma, D_mesh_q,
                                            beta1, beta2, 
                                            hazdenom, 
                                            mesh, 
                                            Dmod = "~x^2"))
                           },
                           mc.cores = 6
)
tot.time.all_q <- difftime(Sys.time(), start.time.all_q, units = "secs")
print(tot.time.all_q)


start.time.all1 <- Sys.time()
all_sim_fits <- mclapply(X = as.list(1:nsims),
                          FUN = function(sim){
                            return(sim_fit(tracksdf, 
                                           traps,
                                           trapcells,
                                           dist_trapmesh,
                                           useall,
                                           lambda0, sigma, D_mesh,
                                           beta1, beta2,
                                           hazdenom, 
                                           mesh, 
                                           Dmod = "~1"))
                          },
                          mc.cores = 6
)
tot.time.all1 <- difftime(Sys.time(), start.time.all1, units = "secs")
print(tot.time.all1)

saveRDS(all_sim_fits_q, "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/simulation_results/variable_dens.Rds")
saveRDS(all_sim_fits, "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/simulation_results/flat_dens.Rds")



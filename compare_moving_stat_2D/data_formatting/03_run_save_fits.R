#fit multiple models and save output

start.time.all_q <- Sys.time()
all_sim_fits_q <- mclapply(X = as.list(1:nsims),
                           FUN = function(sim){
                             return(sim_fit(traps,
                                            tracksdf, 
                                            mesh, 
                                            meshspacing,
                                            dist_trapmesh,
                                            useall,
                                            lambda0, sigma, D_mesh_v, 
                                            beta1, beta2, beta3,
                                            hazdenom, 
                                            Dmod = "~x^2"))
                           },
                           mc.cores = 6
)
tot.time.all_q <- difftime(Sys.time(), start.time.all_q, units = "secs")
print(tot.time.all_q)


start.time.all1 <- Sys.time()
all_sim_fits <- mclapply(X = as.list(1:nsims),
                          FUN = function(sim){
                            return(sim_fit(traps,
                                           tracksdf, 
                                           mesh, 
                                           meshspacing,
                                           dist_trapmesh,
                                           useall,
                                           lambda0, sigma, D_mesh_f, 
                                           beta1, beta2, beta3,
                                           hazdenom,
                                           Dmod = "~1"))
                          },
                          mc.cores = 6
)
tot.time.all1 <- difftime(Sys.time(), start.time.all1, units = "secs")
print(tot.time.all1)

#directory names
dirstart <- paste("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/simulation_results/",
                  "l", lambda0,
                  "D", round(sum(D_mesh_f)*meshspacing^2),
                  "/",
                  sep = "")

# Check and create the directory
if (!dir.exists(dirstart)) {
  dir.create(dirstart, recursive = TRUE)
}

saveRDS(all_sim_fits_q, paste(dirstart,"variable_dens.Rds", sep = ""))
saveRDS(all_sim_fits, paste(dirstart, "flat_dens.Rds", sep = ""))



#fit multiple models and save output
ncores <- 30
cl <- makeCluster(ncores)

# Load packages on workers
clusterEvalQ(cl, 
             {library(Rcpp)
             sourceCpp("approx_movingdetectorlikelihood.cpp")})

# Export all objects/functions needed by sim_fit
clusterExport(cl, varlist = c(
  "sim_fit", "simulate_popandcapthist", "fit_capthist",
  "traps", "tracksdf", "mesh", "meshspacing",
  "dist_trapmesh", "useall", "lambda0", "sigma",
  "D_mesh_v", "D_mesh_f", "beta1", "beta2", "N", "hazdenom"
))

# Run first batch (inhomogeneous)
start.time.all_q <- Sys.time()
all_sim_fits_q <- parLapply(cl, 1:nsims, function(sim) {
  tryCatch({
    sim_fit(traps, tracksdf, mesh, meshspacing,
            dist_trapmesh, useall,
            lambda0, sigma, D_mesh_v, 
            beta1, beta2, N,
            hazdenom, 
            Dmod = "~x^2")
  }, error = function(e) {
    message("Simulation ", sim, " failed: ", e$message)
    return(NULL)
  })
})
tot.time.all_q <- difftime(Sys.time(), start.time.all_q, units = "secs")
print(tot.time.all_q)

# Run second batch (~1)
start.time.all1 <- Sys.time()
all_sim_fits <- parLapply(cl, 1:nsims, function(sim) {
  tryCatch({
    sim_fit(traps, tracksdf, mesh, meshspacing,
            dist_trapmesh, useall,
            lambda0, sigma, D_mesh_f, 
            beta1, beta2, N,
            hazdenom,
            Dmod = "~1")
  }, error = function(e) {
    message("Simulation ", sim, " failed: ", e$message)
    return(NULL)
  })
})
tot.time.all1 <- difftime(Sys.time(), start.time.all1, units = "secs")
print(tot.time.all1)

stopCluster(cl)


saveRDS(all_sim_fits_q, paste(dirstart,"variable_dens.Rds", sep = ""))
saveRDS(all_sim_fits, paste(dirstart, "flat_dens.Rds", sep = ""))



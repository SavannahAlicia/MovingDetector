#2D
#----------------------------------data setup-----------------------------------
#multiple tracklines, keep separate occasions since we only take first detection
#per occasion

#each trackline is a series of points with x, y, and time
trackxmin = 1500
trackxmax = 3500
tracksteplength = 125/5
tracksteps = (trackxmax - trackxmin)/tracksteplength #intervals 
trackint = 360 #seconds (doesn't really matter for length based hazard as long as its positive)
tracksdf <- rbind(
  data.frame(occ = 1,
             x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
             y = 1250, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
  data.frame(occ = 2,
             x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
             y = 1500, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
  data.frame(occ = 3,
             x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
             y = 1750, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
  data.frame(occ = 4,
             x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
             y = 2000, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                        by = trackint))
)
tracksteplength <- abs(tracksdf[2,"x"] - tracksdf[1,"x"])
nocc <- length(unique(tracksdf$occ))

#mesh grid
meshspacing = tracksteplength * 5
mesh <- make.mask(tracksdf[,c("x","y")], buffer = 5*sigma, spacing = meshspacing)
D_mesh <- rep(flatD, nrow(mesh))
D_mesh_q <- exp(beta1*(mesh$x + beta2)^2)

hazdenom <- 1 #hazard is per time or distance, currently specified as distance

#trap grid
trapspacing = meshspacing
xgr <- seq(min(tracksdf$x)-(.5*trapspacing), max(tracksdf$x)+(.5*trapspacing), by = trapspacing)
ygr <- seq(min(tracksdf$y), max(tracksdf$y), by = trapspacing)
gr <- expand.grid(xgr, ygr)
# allocate track records to grid
trapno <- rep(0, nrow(tracksdf))
for (tr in 1:length(trapno)) {
  cat(tr, " / ", length(trapno), "\r")
  d <- sqrt((tracksdf[tr, ]$x - gr[, 1])^2 + (tracksdf[tr, ]$y - gr[, 2])^2)
  dmin <- min(d)
  trapno[tr] <- which.min(d)
}

# pick out possible traps as used cells
un <- sort(unique(trapno))
tracksdf$trapno <- apply(as.array(1:length(trapno)), 1, FUN = function(x){which(un == trapno[x])})
traps <- data.frame(x = gr[un, 1],
                    y = gr[un, 2])

#create trapgrid (will break tracklines into trap grid, keeping times)
trap_cells <- create_grid_bboxes_C(traps, trapspacing)
#all lines
tracklines <- create_line_list_C(tracksdf, scenario = "everything")

emptych <- array(NA, dim = c(1, nocc, nrow(traps)))
useallC <- create_ind_use_C(emptych, traps, trapspacing, tracksdf, scenario = "everything")

#use for undetected inds
useall <- as.matrix(useallC[1,,])

#calculate distance matrix for all trap cells and mesh cells
dist_trapmesh <- apply(as.matrix(mesh), 1, function(meshm){
  apply(as.matrix(traps), 1, function(trapj){
    dist(rbind(trapj,
               meshm), method = "euclidean")
  }) })
dist_trapmesh <- calc_dist_matC(as.matrix(traps),(as.matrix(mesh)))

ggplot() +
  geom_raster(data.frame(x = mesh$x, y = mesh$y, D = D_mesh_q), 
             mapping = aes(x = x, y = y, fill = D)) +
  geom_point(data = tracksdf, mapping = aes(x = x, y = y, group = occ),
             size = 3,shape = "+", color= "white") +
  scale_fill_viridis_c() +
  #xlim(1000, 4000) +
  geom_point(data.frame(x = traps$x,
                        y = traps$y),
             mapping = aes(x = x, y = y), shape = 21, , color = "pink", fill = "red", alpha = .7, size = 2) +
  theme_bw()


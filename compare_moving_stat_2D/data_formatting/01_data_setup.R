#2D
#----------------------------------data setup-----------------------------------
#multiple tracklines, keep separate occasions since we only take first detection
#per occasion

#each trackline is a series of points with x, y, and time
trackxmin = 1500
trackxmax = 3500
trapspacing = 25
tracksteplength = trapspacing/5

tracksteps = (trackxmax - trackxmin)/tracksteplength #intervals 
trackint = 360 #seconds (doesn't really matter for length based hazard as long as its positive)
tracksdf <- rbind(
  data.frame(occ = 1,
             x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
             y = 1300, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
  data.frame(occ = 2,
             x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
             y = 1400, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
   data.frame(occ = 3,
              x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
              y = 1500, 
              time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                         by = trackint)),
   data.frame(occ = 4,
              x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
              y = 1600, 
              time = seq(ymd_hms("2024-01-02 0:00:00"), (ymd_hms("2024-01-02 0:00:00") + (tracksteps)*trackint), 
                         by = trackint)),
   data.frame(occ = 5,
              x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
              y = 1300, 
              time = seq(ymd_hms("2024-01-02 0:00:00"), (ymd_hms("2024-01-02 0:00:00") + (tracksteps)*trackint), 
                         by = trackint)),
   data.frame(occ = 6,
              x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
              y = 1400, 
              time = seq(ymd_hms("2024-01-02 0:00:00"), (ymd_hms("2024-01-02 0:00:00") + (tracksteps)*trackint), 
                         by = trackint)),
  data.frame(occ = 7,
             x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
             y = 1500, 
             time = seq(ymd_hms("2024-01-03 0:00:00"), (ymd_hms("2024-01-03 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
  data.frame(occ = 8,
             x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
             y = 1600, 
             time = seq(ymd_hms("2024-01-03 0:00:00"), (ymd_hms("2024-01-03 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
  data.frame(occ = 9,
             x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
             y = 1300, 
             time = seq(ymd_hms("2024-01-03 0:00:00"), (ymd_hms("2024-01-03 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
   data.frame(occ = 10,
              x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
              y = 1400, 
              time = seq(ymd_hms("2024-01-04 0:00:00"), (ymd_hms("2024-01-04 0:00:00") + (tracksteps)*trackint), 
                         by = trackint)),
   data.frame(occ = 11,
              x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
              y = 1500, 
              time = seq(ymd_hms("2024-01-04 0:00:00"), (ymd_hms("2024-01-04 0:00:00") + (tracksteps)*trackint), 
                         by = trackint)),
   data.frame(occ = 12,
              x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
              y = 1600, 
              time = seq(ymd_hms("2024-01-04 0:00:00"), (ymd_hms("2024-01-04 0:00:00") + (tracksteps)*trackint), 
                         by = trackint))
)
tracksteplength <- abs(tracksdf[2,"x"] - tracksdf[1,"x"])
nocc <- length(unique(tracksdf$occ))

#mesh grid
meshspacing = trapspacing
mesh <- make.mask(tracksdf[,c("x","y")], buffer = 3*sigma, spacing = meshspacing)
D_mesh_f <- rep(flatD, nrow(mesh))
D_mesh_v <- exp(beta1*(mesh$x + beta2)^2 + beta3)
#rescale Dmeshv 
beta3 <- log(sum(D_mesh_f)/sum(D_mesh_v))
D_mesh_v <- exp(beta1*(mesh$x + beta2)^2 + beta3)


hazdenom <- 1 #hazard is per time or distance, currently specified as distance

#trap grid

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
dist_trapmesh <- calc_dist_matC(as.matrix(traps),(as.matrix(mesh)))

layoutplot <- ggplot() +
  geom_raster(data.frame(x = mesh$x, y = mesh$y, D = D_mesh_v), 
             mapping = aes(x = x, y = y, fill = D)) +
  geom_point(data = tracksdf, mapping = aes(x = x, y = y, group = occ),
             size = 3,shape = "+", color= "white") +
  scale_fill_viridis_c()
layoutplot +
  #xlim(1000, 4000) +
  geom_point(data.frame(x = traps$x,
                        y = traps$y),
             mapping = aes(x = x, y = y), shape = 21, , color = "pink", fill = "red", alpha = .7, size = 2) +
  theme_bw()

# visualize a capture history
testpop <- sim_pop_C(D_mesh_v, 
                     as.matrix(mesh), 
                     meshspacing)
testdist_dat_pop <- calc_dist_matC(testpop, 
                               as.matrix(tracksdf[,c("x","y")]))

testcapthist_full <- sim_capthist_C(as.matrix(traps),
                                tracksdf, 
                                lambda0,
                                sigma,
                                D_mesh_v,
                                as.matrix(mesh),
                                meshspacing,
                                hazdenom,
                                testpop,
                                testdist_dat_pop,
                                report_probseenxk = F)
testcapthist <- testcapthist_full[which(apply((!is.na(testcapthist_full)), 1, sum)>0),,]
testinduse <- create_ind_use_C(testcapthist_full, as.matrix(traps),
                               trapspacing, tracksdf, scenario = "everything")

tocck = sample(nocc,1)
tindi = sample(nrow(testpop),1)
ggplot() + 
  scale_fill_viridis_c()+
  geom_point(data.frame(x = traps$x,
                        y = traps$y,
                        induse1 = testinduse[tindi,,tocck]),
             mapping = aes(x = x ,y = y, fill = induse1),
                      shape = 21, size = 2) +
  new_scale_fill() +
  geom_point(data.frame(x = testpop[,1],
                        y = testpop[,2],
                        detocci = apply(testcapthist_full[,tocck,], 1, sum, na.rm = T)>0,
                        interest = seq(1:nrow(testpop)) == tindi
  ), mapping = aes(x = x, y = y, fill = detocci, color = interest), shape = 21, stroke = 1.5) +
  scale_color_manual(values = c("transparent", "red")) +
  scale_fill_discrete() +
  geom_point(data.frame(x = traps$x[which(!is.na(testcapthist_full[tindi,tocck,]))], 
                        y = traps$y[which(!is.na(testcapthist_full[tindi,tocck,]))]
                        ), mapping = aes(x = x, y = y), color = "red", shape = 1,
             stroke = 1.5) 

popdf <- data.frame(x = testpop[,1],
                    y = testpop[,2],
                    totdets = apply(testcapthist_full, 1, function(i){sum(!is.na(i))})
                    )

layoutplot + geom_point(popdf, mapping = aes(x = x, y = y, color = totdets), size = 3) +
  scale_color_viridis_c(option = "magma")

layoutplot +   
  geom_point(data.frame(x = traps$x,
                                     y = traps$y,
                                     dets = apply(testcapthist_full, 3, function(j){sum(!is.na(j))})),
                          mapping = aes(x = x, y = y, fill = dets, color = dets), 
              size = 3) +
  scale_color_viridis_c(option = "magma")

  
    

#each trackline is a series of points with x, y, and time
trackxmin = 1500
trackxmax = 3500
tracksteplength = 125/5
diaglength = sqrt(tracksteplength^2/2)
tracksteps = (trackxmax - trackxmin)/tracksteplength #intervals 
trackint = 360 #seconds (doesn't really matter for length based hazard as long as its positive)
swlen = 90
#create winding stream
streamdf <- data.frame(x = c(trackxmin,
                             trackxmin+tracksteplength,
                             c(trackxmin+tracksteplength + tracksteplength*1:swlen),
                             rep(trackxmin+tracksteplength + tracksteplength*swlen, 5),
                             trackxmin+tracksteplength + tracksteplength*(swlen-1):1,
                             rep(trackxmin+tracksteplength + tracksteplength*1, 4),
                             c(trackxmin+tracksteplength + tracksteplength*1:swlen),
                             rep(trackxmin+tracksteplength + tracksteplength*swlen, 5),
                             trackxmin+tracksteplength + tracksteplength*(swlen-1):1,
                             rep(trackxmin+tracksteplength + tracksteplength*1, 4),
                             c(trackxmin+tracksteplength + tracksteplength*1:swlen),
                             rep(trackxmin+tracksteplength + tracksteplength*swlen, 5),
                             trackxmin+tracksteplength + tracksteplength*(swlen-1):1,
                             rep(trackxmin+tracksteplength + tracksteplength*1, 4),
                             c(trackxmin+tracksteplength + tracksteplength*1:7)
),
y = c(0,
      0,
      0*1:swlen,
      0*swlen + tracksteplength*1:5,
      0*(swlen-1):1 + tracksteplength*5,
      0*1 + tracksteplength*6:9,
      0*1:swlen + tracksteplength*10,
      0*swlen + tracksteplength*11:15,
      0*(swlen-1):1 + tracksteplength*15,
      0*1 + tracksteplength*16:19,
      0*1:swlen + tracksteplength*20,
      0*swlen + tracksteplength*21:25,
      0*(swlen-1):1 + tracksteplength*25,
      0*1 + tracksteplength*26:29,
      0*1:7 + tracksteplength*30
))
#turn it into tracksdf
tracksteps = nrow(streamdf)-1
tracksdf <- rbind(
  data.frame(occ = 1,
             x = streamdf$x,
             y = streamdf$y, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
  data.frame(occ = 2,
             x = streamdf$x,
             y = streamdf$y, 
             time = seq(ymd_hms("2024-01-02 0:00:00"), (ymd_hms("2024-01-02 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
  data.frame(occ = 3,
             x = streamdf$x,
             y = streamdf$y,  
             time = seq(ymd_hms("2024-01-03 0:00:00"), (ymd_hms("2024-01-03 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
  data.frame(occ = 4,
             x = streamdf$x,
             y = streamdf$y, 
             time = seq(ymd_hms("2024-01-04 0:00:00"), (ymd_hms("2024-01-04 0:00:00") + (tracksteps)*trackint), 
                        by = trackint))
)
tracksteplength <- abs(tracksdf[2,"x"] - tracksdf[1,"x"])
nocc <- length(unique(tracksdf$occ))

#mesh grid
meshspacing = tracksteplength * meshstepmult
mesh <- secrlinear::read.linearmask(data = rbind(data.frame(x = min(streamdf$x)-3*sigma, y = 0),
                                                 streamdf[,c("x","y")],
                                                 data.frame(x = streamdf$x[(tracksteps+1)] + sqrt((3*sigma)^2/2),
                                                            y = streamdf$y[(tracksteps+1)] + sqrt((3*sigma)^2/2))), 
                                    spacing = meshspacing)

D_mesh <- rep(flatD, nrow(mesh))
D_mesh_q <- beta3*exp(beta1*(mesh$x + beta2)^2)
hazdenom <- 1 #hazard is per time or distance, currently specified as distance



#trap grid
trapspacing = meshspacing
xgr <- seq(min(tracksdf$x)+50, max(tracksdf$x)+50, by = trapspacing)
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
trap_cells <- create_grid_polygons(traps, spacing = trapspacing)
#all lines
tracklines <- create_line_spatlines(tracksdf)

#use for undetected inds
getuse <- function(oc){
  usecol <- lengths_in_grid(tracklines, oc, trap_cells)
  return(usecol)
}
useall <- matrix(0, nr = nrow(trap_cells), nc = nocc)
colnames(useall) <- 1:nocc
useall[,c(1:ncol(useall))] <- do.call(cbind,
                                      mclapply(X= as.list(1:ncol(useall)), 
                                               FUN = getuse, mc.cores = 3))
# #assemble into a i,j,k use matrix
# induse_ls <- create_ind_use(exch, trapcells, tracksdf) #takes about 30 seconds
# induse <- aperm(
#   array(unlist(induse_ls), 
#         dim =c(nrow(traps), nocc, length(induse_ls))), 
#   c(3, 1, 2))
# saveRDS(induse, file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/data_objs/induse.Rds")
induse <- readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/data_objs/induse.Rds")

#calculate non Euclidean distance matrix for all trap cells and mesh cells

polypts <- rbind(data.frame(x = min(streamdf$x)-3*sigma, y = 0, time = 0, occ = 1:4 ),
                 cbind(tracksdf[,c("occ","x","y")], data.frame(time= tracksdf$time + 1)),
                 data.frame(x = streamdf$x[(tracksteps+1)] + sqrt((3*sigma)^2/2),
                            y = streamdf$y[(tracksteps+1)] + sqrt((3*sigma)^2/2),
                            occ = 1:4, time = (max(tracksdf$time)+1)))
polypts <- polypts[order(polypts[,"occ"], polypts[,"time"]),]
riverpoly <- st_buffer(st_as_sfc(do.call(rbind, 
                                         create_line_spatlines(polypts)), crs = 26916), 
                       dist = 20)
# trappts <- st_as_sf(x = traps, coords = c("x","y"), crs = 26916)
# connects <- nngeo::st_connect(trappts, riverpoly)
# connects <- connects[ which(as.numeric(st_length(connects)) > 0.001)]
# conbuf <- st_buffer(connects, 20)
# all_poly <- as_Spatial(st_union(x = riverpoly,y = conbuf))
# r <- raster(ncol = 1000, nrow = 1000)
# extent(r) <- extent(all_poly)
# rp <- rasterize(all_poly, r)
# rp2 <- rasterize(as(all_poly, "SpatialLines"), r)
# crs(rp) <- crs(rp2) <- crs(all_poly)
# values(rp)[!is.na(values(rp))] = 1
# values(rp2)[!is.na(values(rp2))] = 1
# rp_df <- as.data.frame(as(rp, "SpatialPixelsDataFrame"))
# rp_df2 <- as.data.frame(as(rp2, "SpatialPixelsDataFrame"))
# colnames(rp_df) <- colnames(rp_df2) <- c("value", "x", "y")
# rp_dfm <- rbind(rp_df, rp_df2)
# rp_m <- rasterFromXYZ(rp_dfm[,c("x", "y", "value")], crs = crs(all_poly))
# trans <- gdistance::transition(rp_m, mean, directions = 16)
# trans.c <- gdistance::geoCorrection(trans, type = "c")
# saveRDS(trans.c, file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/data_objs/transc.Rds")
trans.c <- readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/data_objs/transc.Rds")
userdfn1 <- function (xy1, xy2, trans.c) {
  gdistance::costDistance(trans.c, as.matrix(xy1), as.matrix(xy2))
}

dist_trapmesh <- userdfn1(traps[,1:2], mesh[,1:2], trans.c)

ggplot() +
  geom_sf(riverpoly, mapping = aes(), fill = "lightblue") +
  geom_point(data.frame(x = mesh$x, y = mesh$y, D = D_mesh_q), mapping = aes(x = x, y = y, alpha = D), shape = 21) +
  geom_point(data = tracksdf, mapping = aes(x = x, y = y, group = occ), shape = "+") +
  scale_color_viridis_d() +
  geom_sf(st_as_sfc(do.call(rbind, create_line_spatlines(tracksdf)), crs = 26916),
          mapping = aes()) +
  geom_point(data.frame(x = traps$x,
                        y = traps$y, 
                        dets = as.factor(apply((!is.na(excapthist)), 3, sum))),
             mapping = aes(x = x, y = y), shape = 21, fill = "transparent", color = "red") +
  coord_sf(crs = 26916)


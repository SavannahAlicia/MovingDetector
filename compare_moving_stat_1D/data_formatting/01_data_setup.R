## 1 D
#each trackline is a series of points with x, y, and time
trackxmin = 11500
trackxmax = 3500
tracksteplength = 1250/5
diaglength = sqrt(tracksteplength^2/2)
tracksteps = 80 #intervals 
trackint = 360 #seconds (doesn't really matter for length based hazard as long as its positive)
swlen = 91
ystretch = 3 #must be integer
#create winding stream

streamdf <- data.frame(x = c(trackxmin, #1
                             trackxmin+tracksteplength, #1
                             c(trackxmin+tracksteplength + tracksteplength*1:swlen),#91
                             rep(trackxmin+tracksteplength + tracksteplength*swlen, 5*ystretch), #5
                             trackxmin+tracksteplength + tracksteplength*(swlen-1):1, #90
                             rep(trackxmin+tracksteplength + tracksteplength*1, (5*ystretch-1)), #4
                             c(trackxmin+tracksteplength + tracksteplength*1:swlen), #91
                             rep(trackxmin+tracksteplength + tracksteplength*swlen, 5*ystretch),#5
                             trackxmin+tracksteplength + tracksteplength*(swlen-1):1, #90
                             rep(trackxmin+tracksteplength + tracksteplength*1, (5*ystretch-1)),#4
                             c(trackxmin+tracksteplength + tracksteplength*1:swlen), #91
                             rep(trackxmin+tracksteplength + tracksteplength*swlen, 5*ystretch),#5
                             trackxmin+tracksteplength + tracksteplength*(swlen-1):1, #90
                             rep(trackxmin+tracksteplength + tracksteplength*1, (5*ystretch-1)),
                             c(trackxmin+tracksteplength + tracksteplength*1:7)
),
y = c(0, #1
      0, #1
      0*1:swlen, #91
      0*swlen + tracksteplength*(0 + 1:(5*ystretch)), #5
      0*(swlen-1):1 + tracksteplength*1*(5*ystretch), #90 #same as last of 5
      0*1 + tracksteplength*((5 * ystretch * 1) + 1:(5*ystretch-1)),#4
      0*1:swlen + tracksteplength*(2*5*ystretch), #91 increase from last of 4
      0*swlen + tracksteplength*((5 * ystretch * 2) + 1:(5*ystretch)), #5
      0*(swlen-1):1 + tracksteplength*3*(5*ystretch), #90
      0*1 + tracksteplength*((5 * ystretch * 3) + 1:(5*ystretch-1)), #4
      0*1:swlen + tracksteplength*(4*5*ystretch),#91
      0*swlen + tracksteplength*((5 * ystretch * 4) + 1:(5*ystretch)), #5
      0*(swlen-1):1 + tracksteplength*5*(5*ystretch), #90
      0*1 + tracksteplength*((5 * ystretch * 5) + 1:(5*ystretch-1)), #4
      0*1:7 + tracksteplength*6*(5*ystretch) #7
))
#turn it into tracksdf
tracksteps = nrow(streamdf)/3 -1
tracksdf <- rbind(
  data.frame(occ = 1,
             x = streamdf$x[1:(tracksteps+1)],
             y = streamdf$y[1:(tracksteps+1)], 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
  data.frame(occ = 2,
             x = streamdf$x[(1:(tracksteps+1) + (tracksteps+1))],
             y = streamdf$y[(1:(tracksteps+1) + (tracksteps+1))], 
             time = seq(ymd_hms("2024-01-02 0:00:00"), (ymd_hms("2024-01-02 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
  data.frame(occ = 3,
             x = streamdf$x[(1:(tracksteps+1) + 2*(tracksteps+1))],
             y = streamdf$y[(1:(tracksteps+1) + 2*(tracksteps+1))],  
             time = seq(ymd_hms("2024-01-03 0:00:00"), (ymd_hms("2024-01-03 0:00:00") + (tracksteps)*trackint), 
                        by = trackint)),
  data.frame(occ = 4,
             x = streamdf$x[1:(tracksteps+1)],
             y = streamdf$y[1:(tracksteps+1)], 
             time = seq(ymd_hms("2024-01-04 0:00:00"), (ymd_hms("2024-01-04 0:00:00") + (tracksteps)*trackint), 
                        by = trackint))
)
tracksteplength <- abs(tracksdf[2,"x"] - tracksdf[1,"x"])
nocc <- length(unique(tracksdf$occ))

#mesh grid
meshspacing = tracksteplength * meshstepmult
meshlin <- secr::read.mask(data = unique(rbind(data.frame(x = seq(streamdf$x[3]-floor((3*sigma)/meshspacing)*meshspacing, 
                                                           to = streamdf$x[3],
                                                           by = meshspacing),
                                                   y = 0),
                                        streamdf[which(streamdf$x %in% c(seq(streamdf$x[3], max(streamdf$x)+meshspacing, by = meshspacing)) &
                                                         streamdf$y %in% c(seq(streamdf$y[1], max(streamdf$y)+meshspacing, by = meshspacing))),],
                                        data.frame(x = streamdf$x[(3*(tracksteps+1))] + sqrt(meshspacing^2/2) * 1:ceiling((3*sigma)/meshspacing),
                                                   y = streamdf$y[(3*(tracksteps+1))] + sqrt(meshspacing^2/2) * 1:ceiling((3*sigma)/meshspacing)))), 
                           spacing = meshspacing)
streamwidth = (meshspacing*1.2)*2
D_meshlin <- rep(flatD, nrow(meshlin))*(streamwidth/1000)
D_meshlin_q <- exp(beta1*(meshlin$x + beta2)^2)*(streamwidth/1000)
hazdenom <- 1 #hazard is per time or distance, currently specified as distance

ggplot() +
  geom_point(data.frame(x = meshlin$x, y = meshlin$y, D = D_meshlin_q), 
             mapping = aes(x = x, y = y), shape = 21) +
  scale_color_viridis_d() +
  geom_sf(st_as_sfc(do.call(rbind, create_line_spatlines(tracksdf)), crs = 26916),
                    mapping = aes())  +
  geom_point(data = tracksdf, mapping = aes(x = x, y = y, group = occ, color = as.factor(occ)), 
             size = 3,shape = "+") +
  coord_sf(crs = 26916)

#trap grid
trapspacing = meshspacing
xgr <- seq(min(meshlin$x), max(meshlin$x), by = trapspacing)
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

ggplot() +
  geom_point(data.frame(x = meshlin$x, y = meshlin$y, D = D_meshlin_q), 
             mapping = aes(x = x, y = y), shape = 21) +
  geom_point(data = tracksdf, mapping = aes(x = x, y = y, group = occ), shape = "+") +
  geom_point(data = traps, mapping = aes(x = x, y = y), color = "red") +
  geom_sf(st_as_sfc(do.call(rbind, create_line_spatlines(tracksdf)), crs = 26916),
          mapping = aes())  +
  coord_sf(crs = 26916)

#create trapgrid (will break tracklines into trap grid, keeping times)
trapcells <- create_grid_polygons(traps, spacing = trapspacing)
#all lines
tracklines <- create_line_spatlines(tracksdf)

#use for undetected inds
getuse <- function(oc){
  usecol <- lengths_in_grid(tracklines, oc, trapcells)
  return(usecol)
}
useall <- matrix(0, nr = nrow(trapcells), nc = nocc)
colnames(useall) <- 1:nocc
useall[,c(1:ncol(useall))] <- do.call(cbind,
                                      mclapply(X= as.list(1:ncol(useall)), 
                                               FUN = getuse, mc.cores = 3))

#calculate non Euclidean distance matrix for all trap cells and mesh cells
#both graph distance and 2D
polypts <- data.frame(x = meshlin$x,
                      y = meshlin$y,
                      occ = 1, 
                      t = 1:nrow(meshlin))
riverpoly <- st_buffer(st_as_sfc(do.call(rbind, 
                                         create_line_spatlines(polypts)), crs = 26916), 
                       dist = streamwidth/2)
ggplot() +
  geom_sf(riverpoly, mapping = aes(), fill = "lightblue") +
  geom_sf(st_as_sfc(do.call(rbind, create_line_spatlines(tracksdf)), crs = NULL),
          mapping = aes())  +
  coord_sf(datum = NULL) +
  theme_bw()
trappts <- st_as_sf(x = traps, coords = c("x","y"), crs = 26916)
# connects <- nngeo::st_connect(trappts, riverpoly)
# connects <- connects[ which(as.numeric(st_length(connects)) > 0.001)]
# conbuf <- st_buffer(connects, 20)
all_poly <- as_Spatial(riverpoly)#as_Spatial(st_union(x = riverpoly,y = conbuf))
r <- raster(ncol = 500, nrow = 500)
extent(r) <- extent(all_poly)
rp <- rasterize(all_poly, r)
rp2 <- rasterize(as(all_poly, "SpatialLines"), r)
crs(rp) <- crs(rp2) <- crs(all_poly)
values(rp)[!is.na(values(rp))] = 1
values(rp2)[!is.na(values(rp2))] = 1
rp_df <- as.data.frame(as(rp, "SpatialPixelsDataFrame"))
rp_df2 <- as.data.frame(as(rp2, "SpatialPixelsDataFrame"))
colnames(rp_df) <- colnames(rp_df2) <- c("value", "x", "y")
rp_dfm <- rbind(rp_df, rp_df2)
rp_m <- rasterFromXYZ(rp_dfm[,c("x", "y", "value")], crs = crs(all_poly))
trans <- gdistance::transition(rp_m, mean, directions = 16)
trans.c <- gdistance::geoCorrection(trans, type = "c")
saveRDS(trans.c, file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/data_objs/transc.Rds")
trans.c <- readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/data_objs/transc.Rds")
userdfn1 <- function (xy1, xy2, trans.c) {
  gdistance::costDistance(trans.c, as.matrix(xy1), as.matrix(xy2))
}

#distances along graph
unique(tracksdf$trapno[tracksdf$occ == 1])
meshgraph <- make_graph(edges = c(rep(1:nrow(meshlin), each = 2)[-c(1, 2*nrow(meshlin))]), n = nrow(meshlin), directed = FALSE)
meshgraph$layout <- as.matrix(meshlin)
dist_meshmesh_lin <- distances(meshgraph) * meshspacing

#which mesh indices are traps
meshistraps_lin <- which(do.call(paste, meshlin[,c("x","y")]) %in% do.call(paste, traps[,c("x","y")]))
dist_trapmesh_lin <- dist_meshmesh_lin[meshistraps_lin, ]

#also create 2D mesh for comparison with 2D
xgr <- seq(st_bbox(riverpoly)[1]-1000,st_bbox(riverpoly)[3]+trapspacing, by = trapspacing)
ygr <- seq(st_bbox(riverpoly)[2]-1000, st_bbox(riverpoly)[4]+trapspacing, by = trapspacing)
gr <- expand.grid(xgr, ygr)

mesh2Dxy <- gr[st_intersects(riverpoly, st_as_sf(gr, coords = c("Var1", "Var2"), crs = crs(riverpoly)))[[1]],]
colnames(mesh2Dxy) <- c("x", "y")
mesh2D <-  secr::read.mask(data = mesh2Dxy,
                           spacing = meshspacing)

#take density along the line and divide that by area

D_mesh2D <- rep(flatD, nrow(mesh2D))
D_mesh2D_q <- exp(beta1*(mesh2D$x + beta2)^2)
meshunit_lin <- (meshspacing)/1000
meshunit_2D <- meshspacing^2/1000^2 #
tracksmeshdistmat_2D <- userdfn1(tracksdf[,c("x","y")], mesh2D[,1:2], trans.c)
tracksmeshdistmat_lin <- userdfn1(tracksdf[,c("x","y")], meshlin[,1:2], trans.c)

exch <- sim_capthist(pop = NULL, 
                                        traps, 
                                        tracksdf,
                                        lambda0, 
                                        sigma, 
                                        D_mesh2D,
                                        hazdenom, #for hazard rate
                                        mesh2D,
                                        meshunit_2D,
                                        tracksmeshdistmat_2D) 
excapthist <- exch[apply(exch, 1, function(x){!all(is.na(x))}),,]
excapthist[is.na(excapthist)] <- 0
excapthist[excapthist!=0] <- 1

trapcaps <- apply(excapthist, 3, sum)

ggplot() +
  geom_sf(riverpoly, mapping = aes(), fill = "lightblue") +
  geom_sf(st_as_sfc(do.call(rbind, create_line_spatlines(tracksdf)), crs = NULL),
          mapping = aes())  +
  geom_point(data.frame(x = meshlin$x, y = meshlin$y, D = D_meshlin_q), 
             mapping = aes(x = x, y = y, alpha = D), shape = 21, size = 2) +
  geom_point(data = tracksdf, mapping = aes(x = x, y = y, group = occ), size = .5, 
             alpha = .5, shape = "+") +
  geom_point(data.frame(x = mesh2D$x, y = mesh2D$y, D = D_mesh2D_q), shape = 2,
             mapping = aes(x = x, y = y, alpha = D)) +
  geom_point(data = data.frame(x = traps$x, 
                               y = traps$y,
                               capts = trapcaps),
             mapping = aes(x = x, y = y, color = as.factor(capts)),
             size = 3, shape = 1, stroke = 1) +
  scale_color_viridis_d() +
  coord_sf(datum = NULL) +
  theme_bw()

dist_meshmesh_2D <- userdfn1(mesh2D[,1:2], mesh2D[,1:2], trans.c)
meshistraps_2D <- which(do.call(paste, mesh2D[,c("x","y")]) %in% do.call(paste, traps[,c("x","y")]))
dist_trapmesh_2D <- dist_meshmesh_2D[meshistraps_2D,]


ggsave(file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/plots/setup1Driver.png",
       plot = ggplot() +
         geom_sf(riverpoly, mapping = aes(), fill = "lightblue") +
         geom_sf(st_as_sfc(do.call(rbind, create_line_spatlines(tracksdf)), crs = NULL),
                 mapping = aes())  +
         geom_point(data.frame(x = meshlin$x, y = meshlin$y, D = D_meshlin_q), 
                    mapping = aes(x = x, y = y, alpha = D), shape = 21, size = 2) +
         geom_point(data = tracksdf, mapping = aes(x = x, y = y, group = occ), size = .5, 
                    alpha = .5, shape = "+") +
         geom_point(data = traps, mapping = aes(x = x, y = y), color = "red",
                    size = 1.5, shape = "+") +
         scale_color_viridis_d() +
         geom_point(data.frame(x = mesh2D$x, y = mesh2D$y, D = D_mesh2D_q), shape = 2,
                    mapping = aes(x = x, y = y, alpha = D)) +
         coord_sf(datum = NULL) +
         theme_bw(),
       width = 169,
       height = 169,
       units = c("mm"),
       dpi = 300)

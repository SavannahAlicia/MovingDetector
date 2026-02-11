#2D
#----------------------------------data setup-----------------------------------
#multiple tracklines, keep separate occasions since we only take first detection
#per occasion

#each trackline is a series of points with x, y, and time
ntrapsish = 350 #98/2 #it'll be the first number if there's two types of tracks
trackxmin = -1000
trapspacing = round(sigma/3)
trap_n_horiz = 15 #round(sqrt(ntrapsish))
trap_n_vert = round(ntrapsish/trap_n_horiz)
trackxmax = trackxmin + trapspacing * trap_n_horiz #roughly ntraps x
tracksteplength = round(trapspacing/10)
occreps = 6


tracksteps = (trackxmax - trackxmin)/tracksteplength #intervals 
trackint = 360 #seconds (doesn't really matter for length based hazard as long as its positive)
tracksdf <- rbind(data.frame(occ = 1,
             x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
             y = 0, 
             time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                        by = trackint))#,
            # data.frame(occ = 2,
            #            x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
            #            y = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1)-900, 
            #            time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
            #                       by = trackint))
             )
uniquetracktypes <- length(unique(tracksdf$occ))
nexttrack <- tracksdf
for(i in 2:trap_n_vert){ #total, so including existing df
  nexttrack[,"y"] <- nexttrack[,"y"] + trapspacing
  nexttrack[,"occ"] <- nexttrack[,"occ"] + uniquetracktypes
  tracksdf <- rbind(tracksdf, nexttrack)
}

df2 <- tracksdf
occ_1rep <- length(unique(tracksdf$occ))
for(i in 2:occreps){
  rep <- 1
  df2$occ <- df2$occ + (occ_1rep * rep)
  df2$time <- df2$time + 24*60*60
  tracksdf <- rbind(tracksdf, df2)
  rep <- rep + 1
}


#tracksteplength <- abs(tracksdf[2,"x"] - tracksdf[1,"x"])
nocc <- length(unique(tracksdf$occ))

#mesh grid
meshspacing = trapspacing/2
mesh <- make.mask(tracksdf[,c("x","y")], buffer = 3*sigma, spacing = meshspacing)
#eta <- beta1*((mesh$x + beta2)^2 + (mesh$y + beta2)^2)
#Z   <- sum(exp(eta)) * meshspacing^2
D_mesh_v   <- calcDv(mesh$x,
                     mesh$y,
                     beta1,
                     beta2,
                     N,
                     meshspacing)#N * exp(eta) / Z
flatD <- N/(nrow(mesh) * meshspacing^2)
D_mesh_f <- rep(flatD, nrow(mesh))


hazdenom <- 1 #hazard is per time or distance, currently specified as distance

#trap grid

xgr <- seq(min(tracksdf$x)+(0.5 * trapspacing), max(tracksdf$x), by = trapspacing)
ygr <- seq(min(tracksdf$y),#-(0.5 * trapspacing),
           max(tracksdf$y), by = trapspacing)
gr <- expand.grid(xgr, ygr)
# allocate track records to grid
trapno <- rep(0, nrow(tracksdf))
for (tr in 1:length(trapno)) {
  cat(tr, " / ", length(trapno), "\r")
  d <- sqrt((tracksdf[tr, ]$x - gr[, 1])^2 + (tracksdf[tr, ]$y - gr[, 2])^2)
  dmin <- min(d)
  #this needs to assign to the first trap
  candidatetrap <- which(d==dmin) #but this returns first trap
  if(length(candidatetrap)==1){
    trapno[tr] <- candidatetrap
  } else { #if there are two candidate traps
    #we need to choose the first one
    #check if there is an earlier tracksdfpt
    if(tr>1){ #if its not the first point in tracksdf
      #calculate the midpoint between two trackpoints
      midx <- mean(tracksdf[tr,"x"],tracksdf[tr-1, "x"])
      #return the trap of two candidates closest to that
      trapno[tr] <- candidatetrap[which.min(c(gr[candidatetrap[1],1], gr[candidatetrap[2],1])-midx)]
    } else {
      #if this if the first point in tracksdf, it doesn't matter
      trapno[tr] <- candidatetrap[2]
    }
    #
  }

}

# pick out possible traps as used cells
un <- sort(unique(trapno))
tracksdf$trapno <- apply(as.array(1:length(trapno)), 1, 
                         FUN = function(x){which(un == trapno[x])})
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


#-------visualize
# pick which mesh to check
testmesh <- D_mesh_v

layoutplot <- ggplot() +
  geom_raster(data.frame(x = mesh$x, y = mesh$y, D = testmesh), 
             mapping = aes(x = x, y = y, fill = D)) +
  geom_point(data = tracksdf, mapping = aes(x = x, y = y, group = occ),
             size = 3,shape = "+", color= "white"
             ) +
  scale_fill_viridis_c()
seetrap_plot <- layoutplot +
  #xlim(1000, 4000) +
  geom_point(data.frame(x = traps$x,
                        y = traps$y),
             mapping = aes(x = x, y = y), 
             shape = 21, color = "pink", 
             fill = "red", 
             alpha = .7, size = 2) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = "none",
        text = element_text(size = fontsize))


# visualize a capture history
testpop <- sim_pop_C(testmesh, 
                     as.matrix(mesh), 
                     meshspacing)

testcapthist_full <- sim_capthist_C(as.matrix(traps),
                                tracksdf, 
                                lambda0,
                                sigma,
                                testmesh,
                                as.matrix(mesh),
                                meshspacing,
                                hazdenom,
                                testpop,
                                dist_dat_pop = NULL,
                                report_probseenxk = F)
testcapthist <- testcapthist_full[which(apply((!is.na(testcapthist_full)), 1, sum)>0),,]
testinduse <- create_ind_use_C(testcapthist_full, as.matrix(traps),
                               trapspacing, tracksdf, scenario = "everything")

tocck = sample(nocc,1)
tindi = sample(nrow(testpop),1)
example_ch_plot <- ggplot() + 
  scale_fill_viridis_c(guide = "none")+
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
  scale_color_manual(values = c("transparent", "red"), 
                     labels = NULL,
                     name = "Focus\nindividual\nand trap") +
  scale_fill_discrete(name = "Detected\nthis survey?") +
  geom_point(data.frame(x = traps$x[which(!is.na(testcapthist_full[tindi,tocck,]))], 
                        y = traps$y[which(!is.na(testcapthist_full[tindi,tocck,]))]
                        ), mapping = aes(x = x, y = y), color = "red", shape = 1,
             stroke = 1.5) +
  guides(
    fill = guide_legend(
      order = 1,
      override.aes = list(
        shape = 21,
        colour = "transparent",
        stroke = 1.5
      )
    ),
    color = "none"#guide_legend(
      #order = 2,
      #override.aes = list(
      #  shape = 21,
       # fill = "white",
      #  stroke = 1.5
      #)
   # )
  ) +
  theme_bw() +
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size = fontsize))

detmax <- max(c(apply(testcapthist_full, 1, function(i){sum(!is.na(i))}), 
                apply(testcapthist_full, 3, function(j){sum(!is.na(j))})))

popdf <- data.frame(x = testpop[,1],
                    y = testpop[,2],
                    totdets = apply(testcapthist_full, 1, function(i){sum(!is.na(i))})
                    )

popdet_plot <- layoutplot + 
  geom_point(popdf, mapping = aes(x = x, y = y, 
                                             color = totdets), 
                        size = 3) +
  scale_color_viridis_c(option = "magma", name = "dets",
                        limits = c(0, detmax)) +
  guides(fill = "none") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        legend.position = "none",
        text = element_text(size = fontsize))

trapdet_plot <- layoutplot +   
  geom_point(data.frame(x = traps$x,
                                     y = traps$y,
                                     dets = apply(testcapthist_full, 3, function(j){sum(!is.na(j))})),
                          mapping = aes(x = x, y = y, color = dets), 
              size = 3) +
  geom_vline(xintercept = -beta2 - 3*sigma, color = "white", linetype = "dashed") +
  scale_color_viridis_c(option = "magma", name = "Detections"#,
                       # limits = c(0, detmax)
                        ) +
  guides(fill = "none") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_continuous(expand = c(0,0)) +
  theme_bw() +
  theme(axis.text.y = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank(),
        text = element_text(size = fontsize))

grid.arrange(seetrap_plot,
                  example_ch_plot,
                  popdet_plot,
                  trapdet_plot)
dim(testcapthist)
summary(apply(testcapthist_full, 1, function(x){sum(!is.na(x))}))     # detections per individual

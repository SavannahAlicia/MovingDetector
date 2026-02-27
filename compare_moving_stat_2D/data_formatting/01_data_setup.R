#2D
#----------------------------------data setup-----------------------------------

setup_data <- function(sigma,
                       N,
                       beta1,
                       beta2,
                       ntrapsish,
                       trackxmin,
                       trapspacing,
                       meshspacing,
                       trap_n_horiz,
                       nsteps_pertrap,
                       occreps
){
  
  #multiple tracklines, keep separate occasions since we only take first detection
  #per occasion
  
  #each trackline is a series of points with x, y, and time
  
  trap_n_vert = round(ntrapsish/trap_n_horiz)
  trackxmax = trackxmin + trapspacing * trap_n_horiz #roughly ntraps x
  tracksteplength = round(trapspacing/nsteps_pertrap)
  
  tracksteps = (trackxmax - trackxmin)/tracksteplength #intervals 
  trackint = 360 #seconds (doesn't really matter for length based hazard as long as its positive)
  tracksdf <- rbind(data.frame(occ = 1,
                               x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
                               y = 0, 
                               time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                                          by = trackint),
                               transect = 1,
                               rep = 1)#,
                    # data.frame(occ = 2,
                    #            x = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1),
                    #            y = seq(from = trackxmin, to = trackxmax, length.out = tracksteps+1)-900, 
                    #            time = seq(ymd_hms("2024-01-01 0:00:00"), (ymd_hms("2024-01-01 0:00:00") + (tracksteps)*trackint), 
                    #                       by = trackint))
  )
  uniquetracktypes <- length(unique(tracksdf$occ))
  nexttrack <- tracksdf
  if(trap_n_vert > 1){
    for(i in 2:(trap_n_vert)){ #total, so including existing df
      nexttrack[,"y"] <- nexttrack[,"y"] + trapspacing
      nexttrack[,"occ"] <- nexttrack[,"occ"] + uniquetracktypes
      nexttrack[,"transect"] <- nexttrack[,"transect"] + 1
      tracksdf <- rbind(tracksdf, nexttrack)
    }
  }

  df2 <- tracksdf
  occ_1rep <- length(unique(tracksdf$occ))
  for(i in 2:occreps){
    df2$occ <- df2$occ + (occ_1rep * (i-1))
    df2$time <- df2$time + 24*60*60
    df2$rep <- i
    tracksdf <- rbind(tracksdf, df2)
    
  }
  
  tracksdf[,c("midx", "midy")] <- calc_trackmidpts(tracksdf)
  tracksdf$inc <- c(0,sqrt((tracksdf$x[2:nrow(tracksdf)] - tracksdf$x[1:(nrow(tracksdf)-1)])^2 + 
                             (tracksdf$y[2:nrow(tracksdf)] - tracksdf$y[1:(nrow(tracksdf)-1)])^2))
  #set first inc of each occ to 0
  tracksdf$inc[apply(as.array(unique(tracksdf$occ)), 1, function(x){min(which(tracksdf$occ == x))})] <- 0
  
  nocc <- length(unique(tracksdf$occ))
  
  hazdenom <- 1 #hazard is per time or distance, currently specified as distance
  
  #trap grid
  
  xgr <- seq(min(tracksdf$x)+(0.5 * trapspacing), max(tracksdf$x), by = trapspacing)
  ygr <- seq(min(tracksdf$y),#-(0.5 * trapspacing),
             max(tracksdf$y), by = trapspacing)
  #candidate traps
  gr <- expand.grid(xgr, ygr)
  # allocate track records to grid
  trapno <- rep(0, nrow(tracksdf))
  for (tdf_row in 1:nrow(tracksdf)) {
    cat(tdf_row, " / ", length(trapno), "\r")
    #distance from this tracksdf midpoint (between it and prev) to all grid points
    d <- sqrt((tracksdf[tdf_row, ]$midx - gr[, 1])^2 + (tracksdf[tdf_row, ]$midy - gr[, 2])^2)
    dmin <- min(d)
    #candidate traps are those traps at min distance
    candidatetrap <- which(d==dmin) # of gr
    if(length(candidatetrap)==1){
      trapno[tdf_row] <- candidatetrap
    } else { #if there are two candidate traps
      #we need to choose the first one
      #check if there is an earlier tracksdfpt
      if(tdf_row == 1){ 
        #if this if the first point in tracksdf, it doesn't matter
        trapno[tdf_row] <- candidatetrap[2]
      } else {
        # associate the step with whichever trap is closest to the xy in row tdf_row
        # (which is the endpoint of the step) since this is where detection is 
        # attributed (and effort will be cut off between traps)
        td <- sqrt((gr[candidatetrap,1] - tracksdf[tdf_row,"x"])^2 + 
                     (gr[candidatetrap,2] - tracksdf[tdf_row,"y"])^2)
        candidatetrap <- candidatetrap[which(td == min(td))]
        
        trapno[tdf_row] <- candidatetrap
                
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
  
  #mesh grid
  mesh <- make.mask(traps[,c("x","y")], 
                    buffer = 4*sigma, 
                    spacing = meshspacing)
  
  D_mesh_v  <- calcDv(mesh$x,
                      mesh$y,
                      beta1, #will be 0 for flat D
                      beta2,
                      N,
                      meshspacing)
  D_mesh_f  <- calcDv(mesh$x,
                      mesh$y,
                      0, 
                      0,
                      N,
                      meshspacing)
  
  
  
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
  
  meanstepsize = mean(tracksdf$inc[tracksdf$inc != 0])
  
  out_ls <- list(traps = traps,
                 tracksdf = tracksdf, 
                 mesh = mesh,
                 meshspacing = meshspacing,
                 dist_trapmesh = dist_trapmesh,
                 useall = useall,
                 D_mesh_f = D_mesh_f,
                 D_mesh_v = D_mesh_v,
                 meanstepsize = meanstepsize
                 )
  return(out_ls)
}


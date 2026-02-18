#compare fits between my script and Abinand's

#data setup
#--------------------------------------true parameters -------------------------
lambda0 = .008 #expected number of detections per m of trackline at AC
sigma = 300
N <- 65
beta1 <-  -0.015
beta2 <- 0
calcDv <- function(xs, 
                   ys, 
                   beta1_,
                   beta2_, 
                   N_,
                   meshspacing) {
  eta = beta1_*((xs/meshspacing + beta2_)^2 )#+ (ys + beta2_)^2)
  Z = sum(exp(eta)) * meshspacing^2
  D = N * exp(eta) / Z
  return(D)
}
#2D
#----------------------------------data setup-----------------------------------
#multiple tracklines, keep separate occasions since we only take first detection
#per occasion

#each trackline is a series of points with x, y, and time
ntrapsish = 200 #98/2 #it'll be the first number if there's two types of tracks
trackxmin = -1000
trapspacing = round(sigma/3)
meshspacing = trapspacing/2
trap_n_horiz = 15 #round(sqrt(ntrapsish))
trap_n_vert = round(ntrapsish/trap_n_horiz)
trackxmax = trackxmin + trapspacing * trap_n_horiz #roughly ntraps x
nsteps_pertrap = 10
tracksteplength = round(trapspacing/nsteps_pertrap)
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
mesh <- make.mask(tracksdf[,c("x","y")], buffer = 3*sigma, spacing = meshspacing)

D_mesh_v   <- calcDv(mesh$x,
                     mesh$y,
                     beta1,
                     beta2,
                     N,
                     meshspacing)
D_mesh_f <- calcDv(mesh$x,
                   mesh$y,
                   beta1_ = 0,
                   beta2_ = 0,
                   N,
                   meshspacing)


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



#create input objects for Abinand's
survObj <- simulateScrTrapsMask(
  nxTraps = trap_n_horiz,
  nyTraps = trap_n_horiz,
  trapSpacing = trapspacing,
  nSteps = round(trapspacing/tracksteplength),
  maskSpacing = meshspacing,
  sigma = sigma, 
  N = N,
  b1 = beta1,
  b2 = beta2
)

trapSteps <- survObj$trapSteps

mask <- survObj$mask

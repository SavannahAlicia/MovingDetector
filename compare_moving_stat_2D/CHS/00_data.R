library(ggplot2)
library(sf)
library(parallel)
library(sp)
library(openpopscr)
library(secr)
library(mgcv)
library(dplyr)
library(tidyr)
library(gridExtra)

#-------------------------------functions---------------------------------------

setwd("~/Documents/UniStAndrews/MovingDetector")
Rcpp::sourceCpp("approx_movingdetectorlikelihood.cpp")
setwd("~/Documents/UniStAndrews/Dolphins/Charleston")


spreadD <- function(mesh_dist_mat, lambda0, sigma, D){
  rowSums(apply(as.array(1:nrow(mesh_dist_mat)), 1, function(meshx){
    #this is the probability of detection at x
    probdet = apply(as.array(1:nrow(mesh_dist_mat)), 1, function(x){
      lambda0 * exp(-(mesh_dist_mat[meshx,x]^2)/(2*sigma^2))
    })
    probdet = probdet/sum(probdet)
    D[meshx] * probdet #* area #keep density per square km, not per mesh
  }))
}
#-------------Read files -------------------------------------------------------

lpoly <- readRDS("data/all_scenarios/larger_poly.Rds")
tracks <- readRDS("data/all_scenarios/lbl_trackdat.Rds")
olddat <- readRDS("data/onison/all_occasions/model_objs/chs_noneuclidean_extraprims2000.Rds")
sight <- readRDS("data/onison/all_occasions/sight_for_capthist_noopen_2000.Rds")
sightwn <-  readRDS("data/all_scenarios/sightings_data.Rds")
wpt_summaries <- readRDS("data/all_scenarios/wpt_summaries.Rds")
capthist <- readRDS("data/onison/all_occasions/capthistscr_noopen_2000.Rds")
prim <- readRDS("data/all_scenarios/all_occasions/primary.Rds")
traps <- readRDS("data/onison/all_occasions/trapscr_noopen_2000.Rds")
mesh <- readRDS("data/onison/all_occasions/meshscr_NSbuff_2000.Rds")
meshpoly <- readRDS("data/onison/all_occasions/meshpoly_2000.Rds")
distmat <- readRDS("data/onison/all_occasions/user_ne_dist_mat_2000.Rds")
meshdistmat <- readRDS("data/onison/all_occasions/mesh_dist_mat_2000.Rds")
all_poly_2000 <- readRDS("data/onison/all_occasions/all_poly_2000.Rds")


#-------------Tracks dataframe -------------------------------------------------
#format one primary of data 
olddat$time()[11]+2004
surveys <-  unique(sight[which(sight$occ_key %in% which(prim$primary == 11)),"survey"])
oldoccs <- unique(sight[which(sight$occ_key %in% which(prim$primary == 11)),"occ_key"])
tracks <- tracks[tracks$ID %in% surveys,]
tracks <- tracks[order(tracks$ID, tracks$t),]


#sort by track ID, then t
tracksdf <- data.frame(occ = apply(as.array(1:nrow(tracks)), 1, 
                                   function(x){which(unique(tracks$ID) == tracks$ID[x])}),
                       x = tracks$x,
                       y = tracks$y,
                       effort = tracks$effort,
                       time = tracks$t) #time right now only determines order, not relative time
#note SIghting 1 in survey 710 (occasion 4) happens immediately, and the two trackpts aren't labelled on effort
tracksdf[which(tracksdf$occ == 4 & tracksdf$time < 101),"effort"] <- "OnEffort"
# midpoint is between previous pt and current row pt
tracksdf[,c("midx", "midy")] <- calc_trackmidpts(tracksdf)
# step size from previous to current row pt
tracksdf$inc <- c(0,sqrt((tracksdf$x[2:nrow(tracksdf)] - tracksdf$x[1:(nrow(tracksdf)-1)])^2 + 
                           (tracksdf$y[2:nrow(tracksdf)] - tracksdf$y[1:(nrow(tracksdf)-1)])^2))
#set first inc of each occ to 0
tracksdf$inc[apply(as.array(unique(tracksdf$occ)), 1, function(x){min(which(tracksdf$occ == x))})] <- 0

# set inc to 0 if previous was off effort
for(row in 2:nrow(tracksdf)){
  if(is.na(tracksdf[(row-1), "effort"])){
    tracksdf[row,"inc"] <- 0
  }
  else if(tracksdf[(row-1), "effort"] == "OffEffort"){
    tracksdf[row,"inc"] <- 0
  }
}

traps <- traps[which(rowSums(usage(traps)[,c(oldoccs)])>0),]
trapscr <- traps
saveRDS(trapscr, "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHSinput/traps.Rds")

# allocate track records to grid
trapno <- rep(0, nrow(tracksdf))
for (tdf_row in 1:nrow(tracksdf)) {
  cat(tdf_row, " / ", length(trapno), "\r")
  #distance from this tracksdf midpoint (between it and prev) to all grid points
  d <- sqrt((tracksdf[tdf_row, ]$midx - traps[, 1])^2 + (tracksdf[tdf_row, ]$midy - traps[, 2])^2)
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
tracksdf$trapno <- apply(as.array(1:length(trapno)), 1, FUN = function(x){which(un == trapno[x])})
saveRDS(tracksdf, "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHSinput/tracksdf.Rds")


nocc <- length(unique(tracksdf$occ))
dx = 2000

#-------------Subset capthist -------------------------------------------------
capthist <- subset(capthist, occasions = oldoccs) #for comparison, not used later
sight <- sight[sight$survey %in% surveys,]
sight$occ_key <- apply(as.array(sight$survey), 1, function(x){which(surveys == x)})

emptych <- array(NA, dim = c(1, nocc, nrow(traps)))
useallC <- create_ind_use_C(emptych, traps, dx, tracksdf, scenario = "onison")
useall <- as.matrix(useallC[1,,])
usage(trapscr) <- useall
saveRDS(usage(trapscr), "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHSinput/useall.Rds")

distmatscr <- distmat[which(rowSums(usage(traps)[,c(oldoccs)])>0),]
saveRDS(distmatscr, "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHSinput/distmat_trapmesh.Rds")

#-------------Times of detections ----------------------------------------------
#need to identify time from survey start that each individual is detected in each survey for induse
sightwn <- sightwn[which(sightwn$SurveyNumber %in% surveys &
                           sightwn$CatalogID %in% unique(sight$ID) &
                           sightwn$OnEffort == "Yes"),]
sight <- left_join(sight, sightwn[,c("CatalogID", "SurveyNumber", "Sighting", "lat")], 
                   by=c("ID"="CatalogID","survey" = "SurveyNumber", "lat" = "lat"))

sight$t <- NA

wpt_summaries <- wpt_summaries[which(as.numeric(wpt_summaries$ID) %in% surveys),]
for (i in 1:nrow(sight)){
  survi <- sight[i, "survey"]
  sightnumi <- sight[i, "Sighting"]
  wpt <- wpt_summaries[which(
    wpt_summaries$ID == survi & 
      wpt_summaries$SightingNum == sightnumi), 
  ]
  sight[i, "t"] <- wpt$t
  sight[i, "x"] <- wpt$x
  sight[i, "y"] <- wpt$y #note waypoint logs sighting at different point than 
                    #sightings spreadsheet. Suspect they used "start" instead of
                    # "boat"
}
#discard any second detections
pcomboIDsurv <- do.call(paste, sight[,c("ID", "survey")])
sight$freq <- apply(as.array(pcomboIDsurv), 1, function(x){length(which(pcomboIDsurv == x))})
repeaters <- unique(sight[sight$freq>1, c("ID", "survey")])
sight_single <- sight
for(row in 1:nrow(repeaters)){
  rID = repeaters[row, "ID"]
  rsurvey = repeaters[row, "survey"]
  sightrows = which(sight_single$ID == rID & sight_single$survey == rsurvey)
  print(length(sightrows))
  mults <- sight_single[sightrows,]
  sight_single <- sight_single[-sightrows,]
  sight_single <- rbind(sight_single, mults[order(mults$t),][1,])
}

sight_single$trap <- NA
#add 1 second to sighting times (this is how it is added to tracksdat)
sight_single$t <- sight_single$t + 1

#closest trap
for (tr in 1:nrow(sight_single)) {
  d <- sqrt((sight_single[tr, ]$x - traps[, 1])^2 + (sight_single[tr, ]$y - traps[, 2])^2)
  dmin <- min(d)
  if (dmin > 5000) {
    sight_single[tr, "trap"] <- NA
    next
  }
  if(length(which.min(d))==1){
    sight_single[tr, "trap"] <- which.min(d)
    
  } else {
    print(paste(tr))
  }
}


#format this into a ch of times
cht <- apply(as.array(unique(sight_single$ID)), 1, 
             function(i){
               apply(as.array(surveys), 1, function(k){
                 chik <- rep(NA, nrow(traps))
                 s <- sight_single[which(sight_single$ID == i & 
                                           sight_single$survey ==k),]
                 if(nrow(s) > 0){
                   chik[s$trap] <- s$t
                 }
                 return(chik)
               })
             })

ch_t <- aperm(array(cht, dim = c(nrow(traps), nocc, length(unique(sight_single$ID)))),
              c(3,2,1))
saveRDS(ch_t, "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHSinput/ch_t.Rds")

induse <- create_ind_use_C(ch_t, trapscr, 2000, tracksdf, scenario = "onison")
saveRDS(induse, "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHSinput/induse.Rds")
#convert capthist back to 1s and 0s
ch_10 <- ch_t
ch_10[is.na(ch_10)] <- 0
ch_10[ch_10!=0] <- 1
saveRDS(ch_10, "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHSinput/ch.Rds")

exi <- 200
exk <- which(ch_10[exi,,] > 0, arr.ind = T)[1,1]
exline <- as.data.frame(create_line_list_C(tracksdf,"onison")[[exk]])
testbbx <- create_grid_bboxes_C(traps, 2000)
exj <- which(ch_10[exi,exk,]>0)
exsight <- sight_single[which(sight_single$ID == unique(sight_single$ID)[exi] &
                                sight_single$occ_key == exk),]
exwpts <- wpt_summaries[wpt_summaries$ID == unique(tracks$ID)[exk],]
exwpts <- exwpts[exwpts$SightingNumber == exsight$Sighting,]
ggplot() +
  geom_point(data.frame(x = traps$x, 
                        y = traps$y,
                        use = induse[exi,,exk],
                        cap = ch_10[exi,exk,]),
             mapping = aes(x = x, 
                           y =y, 
                           color = use,
                           size = cap)) +
  geom_vline(xintercept = testbbx[exj,c(1,2)], linetype = "dashed") +
  geom_hline(yintercept = testbbx[exj,c(3,4)], linetype = "dashed") +
  scale_color_viridis_c() +
  geom_segment(exline[which(exline$timeend <= ch_t[exi,exk,exj]),], 
               mapping = aes(x = x1,
                             y = y1,
                             xend = x2, 
                             yend = y2
                                     )) +
  geom_point(exsight,
             mapping = aes(x = x, y = y),
             color = "red") +
  geom_point(exwpts, 
             mapping = aes(
    x = x,
    y =y
  ), color = "green", size = .5) #+
 # scale_x_continuous(limits = c(1165550, 1167500)) +
 # scale_y_continuous(limits = c(3649700, 3650000)) 
 # coord_fixed()
  



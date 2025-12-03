library(ggplot2)
library(sf)
library(parallel)
library(sp)
library(openpopscr)
library(secr)
library(mgcv)
library(dplyr)
library(tidyr)

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
nocc <- length(unique(tracksdf$occ))
dx = 2000


#-------------Subset capthist -------------------------------------------------
capthist <- subset(capthist, occasions = oldoccs)
sight <- sight[sight$survey %in% surveys,]
sight$occ_key <- apply(as.array(sight$survey), 1, function(x){which(surveys == x)})
traps <- traps[which(rowSums(usage(traps)[,c(oldoccs)])>0),]
trapscr <- traps
usage(trapscr) <- usage(traps)[which(rowSums(usage(traps)[,c(oldoccs)])>0),c(oldoccs)]
distmatscr <- distmat[which(rowSums(usage(traps)[,c(oldoccs)])>0),]
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

#closest trap
for (tr in 1:nrow(sight_single)) {
  d <- sqrt((sight_single[tr, ]$x - traps[, 1])^2 + (sight_single[tr, ]$y - traps[, 2])^2)
  dmin <- min(d)
  if (dmin > 5000) {
    sight_single[tr, "trap"] <- NA
    next
  }
  sight_single[tr, "trap"] <- which.min(d)
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

induse <- create_ind_use_C(ch_t, trapscr, 2000, tracksdf, scenario = "onison")
#convert capthist back to 1s and 0s
ch_10 <- ch_t
ch_10[is.na(ch_10)] <- 0
ch_10[ch_10!=0] <- 1

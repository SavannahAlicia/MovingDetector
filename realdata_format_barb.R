setwd("~/Documents/UniStAndrews/Dolphins/barb")
scenario = "onison"
# load secr objects
capthist <- readRDS(paste("data/", scenario, "/capthistscr.Rds", sep = ""))
traps <- readRDS(paste("data/", scenario, "/trapscr.Rds", sep =""))
mesh <- readRDS(paste("data/", scenario, "/meshscr.Rds", sep = ""))
prim <- readRDS("data/all_scenarios/primary.Rds")
alltracks <- readRDS(paste("data/onison/trackdat.Rds"))
survsum <- read.csv("data/all_scenarios/survey_summary.csv")
allsightings <- readRDS("data/onison/sightings.Rds")
sightingwnum <- readRDS("data/all_scenarios/sightings_data.Rds")
wpts <- readRDS("data/onison/wpt_summaries.Rds")

#just subset last primary (11)
relsurvsum <- survsum[survsum$Primary.Period == 11,]
tracks <- alltracks[alltracks$ID %in% relsurvsum$SurveyNum,]
sightings <- allsightings[allsightings$survey %in% relsurvsum$SurveyNum,]
sightings$SightingNumber <- apply(as.array(1:nrow(sightings)),1, FUN = function(row){
  sightingwnum[which(sightingwnum$SurveyNumber == sightings[row, "survey"] & 
                       sightingwnum$lat == sightings[row, "lat"] &
                       sightingwnum$lon == sightings[row, "lon"]),"Sighting"][1]
})
wpts <- wpts[which(!is.na(wpts$SightingNum) & wpts$ID %in% relsurvsum$SurveyNum),]
sightings$datetime <- as_datetime(apply(as.array(1:nrow(sightings)),1, FUN = function(row){
  wpts[wpts$ID == sightings[row, "survey"] &
         wpts$SightingNum == sightings[row, "SightingNumber"],"datetime"]
}))



#find average amount of time spent in a mesh cell and use that as time increment
m_per_supersurvey <- colSums(usage(traps))
mesh_visited_ss <- apply(as.array(1:34), 1, FUN = function(x){
  length(which(usage(traps)[,x] > 0))
})
tau_supersurveys <- apply(as.array(1:34), 1, FUN = function(x){
  start_ss <- min(alltracks[alltracks$occasion == x,"datetime"])
  end_ss <- max(alltracks[alltracks$occasion == x,"datetime"])
  ssdur <- difftime(end_ss, start_ss, units = "secs")
  return(ssdur)
})
#time per mesh visited
timeincr <- median(tau_supersurveys/mesh_visited_ss)

#to make a smaller matrix, set beginning of all traps to the same time
smallest_t <- min(tracks[,"datetime"])
diff_surveys <- data.frame(ID = unique(tracks$ID),
                           diff = apply(as.array(unique(tracks$ID)), 1, FUN = function(x){
  start_ss <- min(tracks[tracks$ID == x,"datetime"])
  diff_from_smallest <- difftime(start_ss, smallest_t, units = "secs")
  return(diff_from_smallest)
}))
for (s in unique(tracks$ID)){
  #subtract diff_from_smallest from all times in each survey
  tracks[tracks$ID == s,"datetime"] <- tracks[tracks$ID == s,"datetime"] - diff_surveys[diff_surveys == s,"diff"]
  sightings[sightings$survey == s, "datetime"] <-  sightings[sightings$survey == s, "datetime"] - diff_surveys[diff_surveys == s,"diff"]
}

incrementsoftime <- seq.POSIXt(smallest_t, max(tracks[,"datetime"]), by = timeincr)


#now I need to snap track locations to times on timeincr
tracks$timestep <-  as_datetime(apply(as.array(1:nrow(tracks)), 1, FUN = function(row){
  realtime <- tracks[row,"datetime"]
  d <- abs(difftime(realtime, incrementsoftime, units = "secs"))
  dmin <- min(d)
  if (dmin > timeincr){
    timestep <- NA
    next
  }
  timestep <- incrementsoftime[which.min(d)]
  return(timestep)
}))

sightings$timestep <- as_datetime(apply(as.array(1:nrow(sightings)), 1, FUN = function(row){
    realtime <- sightings[row,"datetime"]
    d <- abs(difftime(realtime, incrementsoftime, units = "secs"))
    dmin <- min(d)
    if (dmin > timeincr){
      timestep <- NA
      next
    }
    timestep <- incrementsoftime[which.min(d)]
    return(timestep)
  }))

#take average location at that time point
snappedtracks <- do.call(rbind, apply(as.array(unique(tracks$ID)), 1, FUN = function(trackID){
  timesfortrack <- unique(tracks[tracks$ID == trackID,"timestep"])
  do.call(rbind, apply(as.array(1:length(timesfortrack)), 1, FUN = function(t){
    timestep <- timesfortrack[t]
    dat <- tracks[which(tracks$ID == trackID &
                          tracks$timestep == timestep),]
    outrow <- dat[1,]
    outrow$x <- mean(dat$x)
    outrow$y <- mean(dat$y)
    return(outrow)
  }))
})
)
#trapID, x, y, time
trapsdf <- snappedtracks[,c("ID", "x", "y", "timestep")]
dist_dat <- create_distdat(trapsdf, mesh)

#and snap detection times to those times
#merge wpts datetime and sightings 
capthist <- matrix(nrow = length(unique(sightings$ID)),
                   ncol = length(unique(tracks$ID)),
                   dimnames = list(unique(sightings$ID), 1:length(unique(tracks$ID))))


for(c in 1:length(unique(sightings$ID))){
                    catalogID <- unique(sightings$ID)[c]
                    print(paste(catalogID))
  out <- as_datetime(apply(as.array(1:length(unique(tracks$ID))), 1, 
                    FUN = function(occ){
                      surveyID = unique(tracks$ID)[occ]
    relsight <- sightings[sightings$survey == surveyID & 
                            sightings$ID == catalogID,]
    if(nrow(relsight) == 0){
      det = NA
    } else {
      det = relsight$timestep[1]
    }
    return(det)
  }))
  capthist[c,] <- as_datetime(out)
}

negloglikelihood_cpp(-12, 8, rep(2.5, nrow(mesh)), timeincr, capthist, dist_dat)

nll <- function(v){
  lambda0_ <- exp(v[1])
  sigma_ <- exp(v[2])
  D_mesh_ <- exp(v[3]*(mesh$x + v[4])^2)#exp(beta1*(mesh$x + beta2)^2)
  out <- negloglikelihood_cpp(lambda0_, sigma_, D_mesh_, timeincr, capthist, dist_dat)
  return(out)
}

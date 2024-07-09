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
wpts <- readRDS("data/onison/wpt_summaries.Rds")

#just subset last primary (11)
relsurvsum <- survsum[survsum$Primary.Period == 11,]
tracks <- alltracks[alltracks$ID %in% relsurvsum$SurveyNum,]
sightings <- allsightings[allsightings$survey %in% relsurvsum$SurveyNum,]
wpts <- wpts[which(!is.na(wpts$SightingNum) & wpts$ID %in% relsurvsum$SurveyNum),]


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

#just do this to the surveys (not supersurveys) in final primary
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
}

incrementsoftime <- seq.POSIXt(smallest_t, max(tracks[,"datetime"]), by = timeincr)


#now I need to snap track locations to times on timeincr
tracks$timestep <- NA
for (t in 1:length(unique(tracks$ID))){
  trackID <- unique(tracks$ID)[t]
  for (tr in 1:nrow(tracks[tracks$ID == trackID,])) {
    realtime <- tracks[tracks$ID == trackID,"datetime"][tr]
    d <- abs(difftime(realtime, incrementsoftime, unit = "secs"))
    dmin <- min(d)
    if (dmin > 2000) {
      trapno[tr] <- NA
      next
    }
    tracks[tracks$ID == trackID, "timestep"][tr] <- incrementsoftime[which.min(d)]
  }
}
tracks$timestep <- as_datetime(tracks$timestep)

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
trapsdf <- snappedtracks[,c("ID", "x", "y", "datetime")]

#and snap detection times to those times

#merge wpts datetime and sightings 

dist_dat <- create_distdat(trapsdf, mesh)

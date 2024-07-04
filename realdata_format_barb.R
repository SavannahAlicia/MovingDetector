setwd("~/Documents/UniStAndrews/Dolphins/barb")
# load secr objects
capthist <- readRDS(paste("data/", scenario, "/capthistscr.Rds", sep = ""))
traps <- readRDS(paste("data/", scenario, "/trapscr.Rds", sep =""))
mesh <- readRDS(paste("data/", scenario, "/meshscr.Rds", sep = ""))
prim <- readRDS("data/all_scenarios/primary.Rds")
tracks <- readRDS(paste("data/onison/trackdat.Rds"))
survsum <- read.csv("data/all_scenarios/survey_summary.csv")
sightings <- readRDS("data/onison/sightings.Rds")

#just subset last primary (11)
relsurvsum <- survsum[survsum$Primary.Period == 11,]
tracks <- tracks[tracks$ID %in% relsurvsum$SurveyNum,]
sightings <- sightings[sightings$survey %in% relsurvsum$SurveyNum,]


#trapID, x, y, time
trapsdf <- tracks[,c("ID", "x", "y", "datetime")]
dist_dat <- create_distdat(trapsdf, mesh)

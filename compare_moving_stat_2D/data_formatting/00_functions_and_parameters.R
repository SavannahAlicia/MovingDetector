#compare moving and stationary detector fits
library(secr)
library(lubridate)
library(parallel)
library(ggplot2)
library(sp)
library(sf)
library(gridExtra)
library(ggnewscale)
library(RcppClock)
library(dplyr)
library(tidyr)
#setwd("~/Documents/UniStAndrews/MovingDetector")
Rcpp::sourceCpp("approx_movingdetectorlikelihood.cpp")


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
  D = N_ * exp(eta) / Z
    return(D)
}


#----------------------------setup parameters ----------------------------------
ntrapsish = 150 #98/2 #it'll be the first number if there's two types of tracks
trackxmin = -1000
trapspacing = round(sigma/3)
meshspacing = trapspacing/2
trap_n_horiz = 15 #n vertical determined by total traps
nsteps_pertrap = 10
occreps = 5

nsims = 30
n_cores = 30
theseed = 1994
set.seed(theseed)
fontsize = 22
printplots = F
hazdenom = 1

#directory names
dirstart <- paste("compare_moving_stat_2D/simulation_results/",
                  "l", lambda0,
                  "D", N,
                  "/",
                  sep = "")

# Check and create the directory
if (!dir.exists(dirstart)) {
  dir.create(dirstart, recursive = TRUE)
}


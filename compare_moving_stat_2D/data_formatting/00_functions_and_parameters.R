#compare moving and stationary detector fits
library(secr)
library(lubridate)
library(parallel)
library(ggplot2)
library(sp)
library(sf)
library(gridExtra)
library(ggnewscale)
library(dplyr)
#setwd("~/Documents/UniStAndrews/MovingDetector")
Rcpp::sourceCpp("approx_movingdetectorlikelihood.cpp")


#--------------------------------------true parameters -------------------------
lambda0 = .0012 #expected number of detections per m of trackline at AC
sigma = 200
N <- 181
beta1 <- -6e-6
fixed_beta1 <- beta1
beta2 <- 0
calcDv <- function(xs, 
                   ys, 
                   beta1_,
                   beta2_, 
                   N_,
                   meshspacing) {
  eta = beta1_*((xs + beta2_)^2 )#+ (ys + beta2_)^2)
  Z = sum(exp(eta)) * meshspacing^2
  D = N * exp(eta) / Z
    return(D)
}


nsims = 30
set.seed(1994)


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


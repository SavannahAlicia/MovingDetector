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
lambda0 = .08 #expected number of detections per m of trackline at AC
sigma = 300
beta1 <- -(1/40000)
beta2 <- -2500
beta3 <- 0 #THIS WILL CHANGE in 01_data_setup to scale total abundance 
flatD = 12/1000000 #per sqr m

nsims = 2#50
set.seed(1994)



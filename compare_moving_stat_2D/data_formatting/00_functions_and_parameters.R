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
setwd("~/Documents/UniStAndrews/MovingDetector")
Rcpp::sourceCpp("approx_movingdetectorlikelihood.cpp")


#--------------------------------------true parameters -------------------------
lambda0 = 1/25 #expected number of detections per m of trackline at AC
sigma = 500
beta1 <- -(1/40000)
beta2 <- -2500
beta3 <- log(0.082944/102.0933) #ensures that sum of Dmesh = sum Dmeshq
flatD = 36/1000000 #per sqr m

nsims = 50
set.seed(1994)



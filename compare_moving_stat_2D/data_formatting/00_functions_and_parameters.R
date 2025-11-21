#compare moving and stationary detector fits
library(secr)
library(lubridate)
library(parallel)
library(ggplot2)
library(sp)
library(sf)
library(gridExtra)
library(dplyr)
setwd("~/Documents/UniStAndrews/MovingDetector")
Rcpp::sourceCpp("approx_movingdetectorlikelihood.cpp")


#--------------------------------------true parameters -------------------------
lambda0 = .001
sigma = 300
beta1 <- -(1/40000)
beta2 <- -2500
beta3 <- 4.843736 #ensures that sum of Dmesh = sum Dmeshq
flatD = 12 #per sqr km

nsims = 50
set.seed(1994)



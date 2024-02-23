#likelihood formulation
library(lubridate)

#distance
distkxt <- function(k, x, t, dist_dat, timesnap = 10*60){
  #referencing a distance matrix data object
    distmat <- dist_dat$distmat
    timediffs <- abs(difftime(dist_dat$times, t, units = "secs"))
    closesttime <- which(timediffs == min(timediffs, na.rm = T))[1]
    if(timediffs[closesttime] <= timesnap){
      tindex <- closesttime
    } else {
      simpleError(paste("Closest time is more than", timesnap, "seconds away"))
    }
    distout <- distmat[k, x, tindex]
  return(distout)
}
#hazard function
haz <- function(t, x, k, lambda0, sigma, dist_dat){
  distkxt. <- distkxt(k = k, x = x, t = t, dist_dat)
  hazout <- lambda0 * exp(-(distkxt.^2/(2 * sigma^2)))
  return(hazout)
}
#survival function
surv <- function(timestart, timeend, timeincr, x, k, lambda0, sigma, dist_dat){ 
  #note time incr in seconds
  if(timeend < timestart) {
    survout <- NA
  } else {
    times <- as.array(seq(from = timestart, to = timeend, by = timeincr))
    hazs <- apply(as.array(1:length(times)), 1, FUN = function(tt){
      timet <- times[tt]
      haz(timet, x = x, k = k, lambda0 = lambda0, sigma = sigma, dist_dat = dist_dat)})
    integral <- sum(hazs * timeincr)
    survout <- exp(-integral)
  }
  return(survout)
}

#note D_mesh needs to be density at each mesh point (hrc) defined in dist_dat
lambdan <- function(dist_dat, D_mesh, timeincr, lambda0, sigma){
  meshx_array <- as.array(1:nrow(dist_dat$mesh))
  Dx_pdotxs <- apply(meshx_array, 1, FUN = function(meshx){
    Dx <- D_mesh[meshx]
    studystart <- min(dist_dat$times)
    studyend <- max(dist_dat$times)
    surv_eachtrap <- apply(as.array(1:length(dist_dat$traps)), 1, 
                                    FUN = function(trapk){
                                      surv(studystart, studyend, timeincr, 
                                           x = meshx, k = trapk, 
                                           lambda0 = lambda0, sigma = sigma,
                                           dist_dat = dist_dat)}
                          )
    surv_alltraps <- prod(surv_eachtrap)
    pdot <- 1 - surv_alltraps
    Dx_pdotx <- Dx * pdot
    return(Dx_pdotx)
  })
  integral <- sum(Dx_pdotxs)
  return(integral)
}

likelihood <- function(lambda0, sigma, D_mesh, timeincr, capthist, dist_dat){
  lambdan. <- lambdan(dist_dat, D_mesh, timeincr, lambda0, sigma)
  n <- nrow(capthist)
  firstterm <- exp(-lambdan.)/(factorial(n))
  integral_eachi <- apply(as.array(1:n), 1, FUN = function(i){
    DKprod_eachx <- apply(as.array(1:nrow(dist_dat$mesh)), 1, FUN = function(x){
      Sxhx_eachtrap <- 
        apply(as.array(1:length(dist_dat$traps)), 1,
              FUN = function(trapk){
                opentimeindx <- which(!is.na(colSums(dist_dat$distmat[trapk,,], na.rm = F)))
                starttime <- dist_dat$times[min(opentimeindx)]
                if(!is.na(capthist[i,trapk])){
                  dettime <- capthist[i, trapk]
                  Sx <- surv(starttime, dettime, timeincr, x, trapk, lambda0, sigma, dist_dat)
                  hx <- haz(dettime, x, trapk, lambda0, sigma, dist_dat)
                  Sxhxout <- Sx * hx
                } else { #individual never detected by this trap
                  endtime <- dist_dat$times[max(opentimeindx)]
                  Sx <- surv(starttime, endtime, timeincr, x, trapk, lambda0, sigma, dist_dat)
                  Sxhxout <- Sx
                }
              })
      Sxhx_alltraps <- prod(Sxhx_eachtrap)
      DKprod_out <- D_mesh[x] * Sxhx_alltraps 
      return(DKprod_out)
    })
    integral <- sum(DKprod_eachx)
    return(integral)
  })
  secondterm <- prod(integral_eachi)
  out <- firstterm * secondterm
}


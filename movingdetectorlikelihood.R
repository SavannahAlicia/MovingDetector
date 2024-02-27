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
    stop(paste("Closest time is more than", timesnap, "seconds away"))
  }
  distout <- distmat[k, x, tindex]
  return(distout)
}
#hazard function
hazdist <- function(lambda0, sigma, d){
  lambda0 * exp(-(d^2/(2 * sigma^2)))
}
haz <- function(t, x, k, lambda0, sigma, dist_dat){
  distkxt. <- distkxt(k = k, x = x, t = t, dist_dat)
  hazout <- hazdist(lambda0 = lambda0, sigma = sigma, d = distkxt.)
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
    integral <- sum(hazs) #* timeincr?
    survout <- exp(-integral)
  }
  return(survout)
}

#note D_mesh needs to be density at each mesh point (hrc) defined in dist_dat
lambdan <- function(dist_dat, D_mesh, timeincr, lambda0, sigma){
  meshx_array <- as.array(1:nrow(dist_dat$mesh))
  Dx_pdotxs <- apply(meshx_array, 1, FUN = function(meshx){
    Dx <- D_mesh[meshx]
    surv_eachtrap <- apply(as.array(1:nrow(dist_dat$traps)), 1, 
                                    FUN = function(trapk){
                                      opentimeindx <- which(!is.na(colSums(dist_dat$distmat[trapk,,], na.rm = F)))
                                      topentime <- dist_dat$times[min(opentimeindx)]
                                      tclosetime <- dist_dat$times[max(opentimeindx)]
                                      surv(topentime, tclosetime, timeincr, 
                                           x = meshx, k = trapk, 
                                           lambda0 = lambda0, sigma = sigma,
                                           dist_dat = dist_dat)}
                          )
    surv_alltraps <- prod(surv_eachtrap)
    pdot <- 1 - surv_alltraps
    Dx_pdotx <- Dx * pdot
    return(Dx_pdotx)
  })
  integral <- mean(Dx_pdotxs)
  return(integral)
}

negloglikelihood <- function(lambda0, sigma, D_mesh, timeincr, capthist, dist_dat){
  lambdan. <- lambdan(dist_dat, D_mesh, timeincr, lambda0, sigma)
  n <- nrow(capthist)
  integral_eachi <- apply(as.array(1:n), 1, FUN = function(i){
    DKprod_eachx <- apply(as.array(1:nrow(dist_dat$mesh)), 1, FUN = function(x){
      Sxhx_eachtrap <- 
        apply(as.array(1:nrow(dist_dat$traps)), 1,
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
    integral <- mean(DKprod_eachx) #sum * mesh area
    return(integral)
  })
  lognfact <- sum(apply(as.array(0:(n-1)), 1, FUN = function(n_){ log(n-n_)}))
  out <- -1 * (-lambdan. - lognfact + n * sum(log(integral_eachi)))
 }


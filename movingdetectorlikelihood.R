#likelihood formulation
library(lubridate)

#' Calculate distance using dist_dat
#'
#' @param k index for trap
#' @param x index for home range center (or mesh)
#' @param t time as POSIXct object
#' @param dist_dat list of "traps" (data frame of trap names), "mesh" (secr obj),
#' "times" (POSIXct), and "distmat" (trap by mesh by times array)
#' @param timesnap tolerance for difference between t and recorded times in 
#' dist_dat, defaults to 10 minutes
#'
#' @return distance from distmat corresponding to k,x, and time closest to t
#' @export
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

#' Hazard(d)
#'
#' @param lambda0 
#' @param sigma
#' @param d distance in m 
#'
#' @return detection rate
#' @export
hazdist <- function(lambda0, sigma, d){
  lambda0 * exp(-(d^2/(2 * sigma^2)))
}

#' Hazard(t)
#'
#' @param t time index for dist_dat
#' @param x home range center index for dist_dat
#' @param k trap index for dist_dat
#' @param lambda0 
#' @param sigma 
#' @param dist_dat list of "traps" (data frame of trap names), "mesh" (secr obj),
#' "times" (POSIXct), and "distmat" (trap by mesh by times array)
#'
#' @return detection rate
#' @export
haz <- function(t, x, k, lambda0, sigma, dist_dat){
  distkxt. <- distkxt(k = k, x = x, t = t, dist_dat)
  hazout <- hazdist(lambda0 = lambda0, sigma = sigma, d = distkxt.)
  return(hazout)
}

#' Survival Function
#' 
#' @param timestart start time
#' @param timeend end time
#' @param timeincr increment of time used for integral (should be consistent)
#' @param x home range center index for dist_dat
#' @param k trap index for dist_dat
#' @param lambda0 
#' @param sigma 
#' @param dist_dat list of "traps" (data frame of trap names), "mesh" (secr obj),
#' "times" (POSIXct), and "distmat" (trap by mesh by times array)
#' 
#' @return probability home range center x is not detected by trap k from 
#' timestart to timeend
#' @export
surv <- function(timestart, timeend, timeincr, x, k, lambda0, sigma, dist_dat){ 
  #note time incr in seconds
  if(timeend < timestart) {
    survout <- NA
  } else {
    timesopen <- as.array(seq(from = timestart, to = timeend, by = timeincr))
    hazs <- apply(as.array(1:length(timesopen)), 1, FUN = function(tt){
      timet <- timesopen[tt]
      haz(timet, x = x, k = k, lambda0 = lambda0, sigma = sigma, dist_dat = dist_dat)})
    integral <- sum(hazs) #* timeincr?
    survout <- exp(-integral)
  }
  return(survout)
}

#' Lambda_n (from likelihood)
#' 
#' @param dist_dat list of "traps" (data frame of trap names), "mesh" (secr obj),
#' "times" (POSIXct), and "distmat" (trap by mesh by times array)
#' @param D_mesh density at each mesh pt
#' @param timeincr time increment for integral
#' @param lambda0
#' @param sigma 
#' 
#' @return rate of detected home range centers
#' @export
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
  integral <- sum(Dx_pdotxs) * attr(mesh, "area")
  return(integral)
}

#' Negative log likelihood
#' 
#' @param lambda0 
#' @param sigma 
#' @param D_mesh density at each mesh pt
#' @param timeincr increment of time used for integral (should be consistent)
#' @param capthist capture history array of dimension individual by trap that 
#' contains NA if not captured or time of first capture
#' @param dist_dat list of "traps" (data frame of trap names), "mesh" (secr obj),
#' "times" (POSIXct), and "distmat" (trap by mesh by times array)
#' 
#' @return negative log-likelihood value
#' @export
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
    integral <- sum(DKprod_eachx) * attr(mesh, "area") #all mesh same size
    return(integral)
  })
  lognfact <- sum(apply(as.array(0:(n-1)), 1, FUN = function(n_){ log(n-n_)}))
  out <- -1 * (-lambdan. - lognfact + n * sum(log(integral_eachi)))
  return(out)
}

create_distdat <- function(traps, mesh){
  #distance data object (created from traps and mesh)
  traptomesh <- apply(as.array(1:nrow(mesh)), 1, FUN = 
                        function(meshrow){
                          apply(as.array(1:nrow(traps)), 
                                1, FUN = 
                                  function(traprow){
                                    dist(rbind(traps[traprow, c("x", "y")],
                                               mesh[meshrow,c("x","y")]), method = "euclidean")
                                  })})
  dist_dat <- list(traps = data.frame(traps = unique(traps$trapID)),
                   mesh = as.data.frame(mesh),
                   times = sort(unique(traps$time)))
  dist_dat$distmat <- sapply(as.array(dist_dat$times),  FUN = function(timet){ 
    apply(as.array(1:nrow(dist_dat$mesh)), 1, FUN = function(meshcol){
      apply(as.array(1:nrow(dist_dat$traps)), 1, FUN = function(trapid){
        dist <- traptomesh[traps$trapID == trapid & traps$time == timet,meshcol]
        if (length(dist) == 0){
          dist <- NA
        }
        return(dist)
      })
    })
  }, simplify = "array")
  return(dist_dat)
}#see about using this in sim_capthist

#' Simulate moving detector capture history
#' 
#' @param pop defaults to null, population locations from sim.popn
#' @param traps dataframe of trapID, x y trap coordinates, and time of location
#' @param timeincr increment of time used for integral (should be consistent)
#' @param lambda0 
#' @param sigma 
#' @param D_mesh density at each mesh pt (for simulating population if not specified)
#' 
#' @return capture history array of dimension individual by trap that 
#' contains NA if not captured or time of first capture
#' @export
sim_capthist <- function(pop = NULL, traps, timeincrement, lambda0, sigma, D_mesh){
  if(is.null(pop)){
    pop <- sim.popn(D = D_mesh, core = mesh, model2D = "IHP", 
                    Ndist = "poisson", buffertype = "rect")
    rownames(pop) <- NULL
  }
  
  dist_dat_pop <- create_distdat(traps, pop)

  capthist <- lapply(as.list(1:nrow(dist_dat_pop$mesh)), 
                     FUN = function(hrcx){
                       lapply(as.list(1:nrow(dist_dat_pop$traps)),
                              FUN = function(trapk){
                                #can do a dominating process of rate lambda0*T and thin
                                opentimeindx <- which(!is.na(colSums(dist_dat_pop$distmat[trapk,,], na.rm = F)))
                                topentime <- dist_dat_pop$times[min(opentimeindx)]
                                tclosetime <- dist_dat_pop$times[max(opentimeindx)]
                                probseenxk <- 1 - surv(topentime, tclosetime, timeincrement, hrcx, trapk, lambda0, sigma, dist_dat_pop)
                                seenxk_bool <- rbinom(1,1, probseenxk)
                                if (seenxk_bool){
                                  timesopen <- seq(topentime, tclosetime, timeincrement)
                                  detattime <- apply(as.array(1:length(timesopen)), 1, FUN = function(t){
                                    notseenuntil <- surv(topentime, timesopen[t], timeincrement, hrcx, trapk, lambda0, sigma, dist_dat_pop)
                                    seenat <- haz(timesopen[t], hrcx, trapk, lambda0, sigma, dist_dat_pop)
                                    probdettime <- notseenuntil * seenat
                                  })
                                  capik <- timesopen[sample(length(detattime), 1, prob = detattime/sum(detattime))]
                                } else {
                                  capik <- ymd_hms(NA)
                                }
                                return(ymd_hms(capik))
                              })
                     }) 
  
  capthist_array <- t(structure(array(unlist(capthist), dim = c(4,nrow(dist_dat_pop$mesh))), class = c("POSIXct", "POSIXt")))
  return(capthist_array)
}

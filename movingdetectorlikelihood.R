#likelihood formulation
library(lubridate)
library(secr)

#' #' Calculate distance using dist_dat
#' #'
#' #' @param k index for trap
#' #' @param x index for home range center (or mesh)
#' #' @param t time as POSIXct object
#' #' @param dist_dat list of "traps" (data frame of trap names), "mesh" (secr obj),
#' #' "times" (POSIXct), and "distmat" (trap by mesh by times array)
#' #' @param timesnap tolerance for difference between t and recorded times in
#' #' dist_dat, defaults to 10 minutes
#' #'
#' #' @return distance from distmat corresponding to k,x, and time closest to t
#' #' @export
#' distkxt <- function(k, x, t, dist_dat, timesnap = 10*60){
#'   #referencing a distance matrix data object
#'   distmat <- dist_dat$distmat
#'   timediffs <- abs(difftime(dist_dat$times, t, units = "secs"))
#'   closesttime <- which(timediffs == base::min(timediffs, na.rm = T))[1]
#'   if(timediffs[closesttime] <= timesnap){
#'     tindex <- closesttime
#'   } else {
#'     stop(paste("Closest time is more than", timesnap, "seconds away"))
#'   }
#'   distout <- distmat[k, x, tindex]
#'   return(distout)
#' }
#' 
#' #' Hazard(d)
#' #'
#' #' @param lambda0
#' #' @param sigma
#' #' @param d distance in m
#' #'
#' #' @return detection rate
#' #' @export
#' hazdist <- function(lambda0, sigma, d){
#'   lambda0 * exp(-(d^2/(2 * sigma^2)))
#' }
#' 
#' #' Hazard(t)
#' #'
#' #' @param t time index for dist_dat
#' #' @param x home range center index for dist_dat
#' #' @param k trap index for dist_dat
#' #' @param lambda0
#' #' @param sigma
#' #' @param dist_dat list of "traps" (data frame of trap names), "mesh" (secr obj),
#' #' "times" (POSIXct), and "distmat" (trap by mesh by times array)
#' #'
#' #' @return detection rate
#' #' @export
#' haz <- function(t, x, k, lambda0, sigma, dist_dat){
#'   distkxt. <- distkxt(k = k, x = x, t = t, dist_dat)
#'   hazout <- hazdist(lambda0 = lambda0, sigma = sigma, d = distkxt.)
#'   return(hazout)
#' }
#' 
#' #' Survival Function
#' #'
#' #' @param timestart start time
#' #' @param timeend end time
#' #' @param timeincr increment of time used for integral (should be consistent)
#' #' @param x home range center index for dist_dat
#' #' @param k trap index for dist_dat
#' #' @param lambda0
#' #' @param sigma
#' #' @param dist_dat list of "traps" (data frame of trap names), "mesh" (secr obj),
#' #' "times" (POSIXct), and "distmat" (trap by mesh by times array)
#' #'
#' #' @return probability home range center x is not detected by trap k from
#' #' timestart to timeend
#' #' @export
#' surv <- function(timestart, timeend, timeincr, x, k, lambda0, sigma, dist_dat){
#'   #note time incr in seconds
#'   if(timeend < timestart) {
#'     survout <- NA
#'   } else {
#'     timesopen <- as.array(seq(from = timestart, to = timeend, by = timeincr))
#'     hazs <- apply(as.array(1:length(timesopen)), 1, FUN = function(tt){
#'       timet <- timesopen[tt]
#'       haz(timet, x = x, k = k, lambda0 = lambda0, sigma = sigma, dist_dat = dist_dat)})
#'     integral <- sum(hazs) * timeincr
#'     survout <- exp(-integral)
#'   }
#'   return(survout)
#' }
#' 
#' #' Lambda_n (from likelihood)
#' #'
#' #' @param dist_dat list of "traps" (data frame of trap names), "mesh" (secr obj),
#' #' "times" (POSIXct), and "distmat" (trap by mesh by times array)
#' #' @param D_mesh density at each mesh pt
#' #' @param timeincr time increment for integral
#' #' @param lambda0
#' #' @param sigma
#' #'
#' #' @return rate of detected home range centers
#' #' @export
#' lambdan <- function(dist_dat, D_mesh, timeincr, lambda0, sigma){
#'   meshx_array <- as.array(1:nrow(dist_dat$mesh))
#'   Dx_pdotxs <- apply(meshx_array, 1, FUN = function(meshx){
#'     Dx <- D_mesh[meshx]
#'     surv_eachtrap <- apply(as.array(1:nrow(dist_dat$traps)), 1,
#'                                     FUN = function(trapk){
#'                                       opentimeindx <- which(!is.na(colSums(dist_dat$distmat[trapk,,], na.rm = F)))
#'                                       topentime <- dist_dat$times[min(opentimeindx)]
#'                                       tclosetime <- dist_dat$times[max(opentimeindx)]
#'                                       surv(topentime, tclosetime, timeincr,
#'                                            x = meshx, k = trapk,
#'                                            lambda0 = lambda0, sigma = sigma,
#'                                            dist_dat = dist_dat)}
#'                           )
#'     surv_alltraps <- prod(surv_eachtrap)
#'     pdot <- 1 - surv_alltraps
#'     Dx_pdotx <- Dx * pdot
#'     return(Dx_pdotx)
#'   })
#'   integral <- sum(Dx_pdotxs) * attr(mesh, "area")
#'   return(integral)
#' }
#' 
#' #' Negative log likelihood
#' #'
#' #' @param lambda0
#' #' @param sigma
#' #' @param D_mesh density at each mesh pt
#' #' @param timeincr increment of time used for integral (should be consistent)
#' #' @param capthist capture history array of dimension individual by trap that
#' #' contains NA if not captured or time of first capture
#' #' @param dist_dat list of "traps" (data frame of trap names), "mesh" (secr obj),
#' #' "times" (POSIXct), and "distmat" (trap by mesh by times array)
#' #'
#' #' @return negative log-likelihood value
#' #' @export
#' negloglikelihood <- function(lambda0, sigma, D_mesh, timeincr, capthist, dist_dat){
#'   lambdan. <- lambdan(dist_dat, D_mesh, timeincr, lambda0, sigma)
#'   n <- nrow(capthist)
#'   integral_eachi <- apply(as.array(1:n), 1, FUN = function(i){
#'     DKprod_eachx <- apply(as.array(1:nrow(dist_dat$mesh)), 1, FUN = function(x){
#'       Sxhx_eachtrap <-
#'         apply(as.array(1:nrow(dist_dat$traps)), 1,
#'               FUN = function(trapk){
#'                 opentimeindx <- which(!is.na(colSums(dist_dat$distmat[trapk,,], na.rm = F)))
#'                 starttime <- dist_dat$times[min(opentimeindx)]
#'                 if(!is.na(capthist[i,trapk])){
#'                   dettime <- capthist[i, trapk]
#'                   Sx <- surv(starttime, dettime, timeincr, x, trapk, lambda0, sigma, dist_dat)
#'                   hx <- haz(dettime, x, trapk, lambda0, sigma, dist_dat)
#'                   Sxhxout <- Sx * hx
#'                 } else { #individual never detected by this trap
#'                   endtime <- dist_dat$times[max(opentimeindx)]
#'                   Sx <- surv(starttime, endtime, timeincr, x, trapk, lambda0, sigma, dist_dat)
#'                   Sxhxout <- Sx
#'                 }
#'               })
#'       Sxhx_alltraps <- prod(Sxhx_eachtrap)
#'       DKprod_out <- D_mesh[x] * Sxhx_alltraps
#'       return(DKprod_out)
#'     })
#'     integral <- sum(DKprod_eachx) * attr(mesh, "area") #all mesh same size
#'     return(integral)
#'   })
#'   lognfact <- sum(apply(as.array(0:(n-1)), 1, FUN = function(n_){ log(n-n_)}))
#'   out <- -1 * (-lambdan. - lognfact + n * sum(log(integral_eachi)))
#'   return(out)
#' }

create_distdat <- function(traps, mesh){
  #if only one trap, create dummy trap to preserve matrix dimensions
  #note: need to do this for time and mesh too
  dummytrap = F
  if(length(unique(traps$trapID)) == 1){
    dummytrap <- traps[1,]
    dummytrap$trapID <- max(traps$trapID) + 1
    traps = rbind(traps,
                 dummytrap)
    dummytrap = T
  }
  dummymesh = F
  if(nrow(mesh) == 1){
    meshdf <- as.data.frame(mesh) 
    dummymesh <- rbind(meshdf[nrow(meshdf),] + c(-.5,-.5),
                       meshdf[nrow(meshdf),] + c(1.5,0.5)
    )
    mesh <- make.mask(dummymesh, buffer = 0, spacing = 1)
    dummymesh <- T
  }
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
  distmat <- sapply(as.array(dist_dat$times),  FUN = function(timet){ 
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
  if (dummytrap){
    distmat <- array(distmat[1,,], dim = c(1, dim(distmat)[2], dim(distmat)[3]))
  }
  if (dummymesh){
    dist_dat$mesh <- dist_dat$mesh[1,]
    distmat <- array(distmat[,1,], dim = c(dim(distmat)[1], 1, dim(distmat)[3]))
  }
  dist_dat$distmat <- distmat
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
sim_capthist <- function(pop = NULL, traps, timeincrement, lambda0, sigma, D_mesh,
                         report_probseenxk = FALSE){
  if(is.null(pop)){
    pop <- sim.popn(D = D_mesh, core = mesh, model2D = "IHP", 
                    Ndist = "poisson", buffertype = "rect")
    rownames(pop) <- NULL
  }
  
  dist_dat_pop <- create_distdat(traps, pop)#this is now distances between pop hrcs and traps

  capthist <- lapply(as.list(1:nrow(dist_dat_pop$mesh)), #for each individual
                     FUN = function(hrcx){
    #capiks <- data.frame(ymd_hms(capik), nrow(dist_dat_pop$mesh), nrow(dist_dat_pop$traps))
    #for(hrcx in 1:nrow(dist_dat_pop$mesh)){
      #print(paste("mesh", hrcx))
      lapply(as.list(1:nrow(dist_dat_pop$traps)),
                              FUN = function(trapk){
      #for (trapk in 1:nrow(dist_dat_pop$traps)){
        #print(paste("trap", trapk))
    
                                opentimeindx <- which(!is.na(colSums(dist_dat_pop$distmat[trapk,,], na.rm = F)))
                                topentime <- dist_dat_pop$times[min(opentimeindx)]
                                tclosetime <- dist_dat_pop$times[max(opentimeindx)]
                                #prob individual is seen by this trap during the total time it's open
                                probseenxk <- 1 - surv_cpp(topentime, tclosetime, timeincrement, (hrcx - 1), (trapk - 1), lambda0, sigma, dist_dat_pop) #any issue with indexing dist_dat_pop?
                                if (report_probseenxk) { 
                                  return(probseenxk)
                                } else {
                                  seenxk_bool <- rbinom(1,1, probseenxk)
                                  if (seenxk_bool){
                                    timesopen <- seq(topentime, tclosetime, timeincrement)
                                    numt <- length(timesopen)
                                    seenfirstat_t <- apply(as.array(2:(numt)), 1, FUN = function(t){
                                      seenfirstatt <- 
                                        
                                      #probability the first detection happens at t given at least one detection happens
                                      #(surv_cpp(topentime, timesopen[t], timeincrement, (hrcx-1), (trapk-1), lambda0, sigma, dist_dat_pop) + 1e-16) *
                                      #haz_cpp(timesopen[t], (hrcx-1), (trapk-1), lambda0, sigma, dist_dat_pop, timeincrement) * timeincrement
                                        
                                      #probability that at least one detection happens in increment ending at t given at least one detection happens
                                      (surv_cpp(topentime, timesopen[(t-1)], timeincrement, (hrcx-1), (trapk-1), lambda0, sigma, dist_dat_pop) + 1e-16) *
                                      (1 - surv_cpp(timesopen[(t-1)], timesopen[(t)], timeincrement, (hrcx-1), (trapk - 1), lambda0, sigma, dist_dat_pop))
                                      
                                      return(seenfirstatt)
                                      })
                                    seenfirstatt_givenseen <- c(0, seenfirstat_t/probseenxk)
                                    capik <- timesopen[sample(x = c(1:length(timesopen)), size = 1, replace = T,  prob = seenfirstatt_givenseen)]
                                    
                                  } else {
                                    capik <- ymd_hms(NA)
                                  }
                                  
                                  return(ymd_hms(capik)) 
                                  }
                                #return(probseenxk)
                                #capiks[hrcx, trapk] <- capik
                                #print(capik)
                              })
                     }) 
  if(report_probseenxk){
    capthist_array <- t(
        array(unlist(capthist), 
              dim = c(nrow(dist_dat_pop$traps),nrow(dist_dat_pop$mesh)))
      )
  } else {
  capthist_array <- t(
    structure(
      array(unlist(capthist), 
                                      dim = c(nrow(dist_dat_pop$traps),nrow(dist_dat_pop$mesh))), 
                                class = c("POSIXct", "POSIXt")
                                )
    )
  }
  return(capthist_array)
}

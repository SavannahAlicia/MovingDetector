#RTMB llk
#replace all apply with for loops



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
negloglikelihood <- function(pars){
  getAll(dist_dat)
  lambda0 <- exp(pars$loglambda0)
  sigma <- exp(pars$logsigma)
  D_mesh <- exp(pars$logD_mesh)
  out <- as.numeric(0)
  
  surv <- function(timestart, timeend, timeincr, x, k, lambda0, sigma, distmat){ 
    survout <- 0
    #note time incr in seconds
    if(timeend < timestart) {
      survout <- as.numeric(NA)
    } else {
      times <- as.array(seq(from = timestart, to = timeend, by = timeincr))
      hazs <- rep(0, length(times))
      for(tt in 1:length(times)){
        timet <- times[tt]
        hazs[tt] <- haz(timet, x = x, k = k, lambda0 = lambda0, sigma = sigma, distmat = distmat)
      }
      integral <- sum(hazs) #pay attention to time increment here
      survout <- as.numeric(exp(-integral))
    }
    return(survout)
  }
  
  haz <- function(t, x, k, lambda0, sigma, distmat){
    hazout <- 0
    distkxt_ <- distkxt(k = k, x = x, t = t, distmat)
    hazout <- hazdist(lambda0 = lambda0, sigma = sigma, d = distkxt_)
    return(hazout)
  }
  
  hazdist <- function(lambda0, sigma, d){
    h <- 0
    h <- lambda0 * exp(-(d^2/(2 * sigma^2)))
    return(h)
  }
  
  distkxt <- function(k, x, t, distmat, timesnap = 10*60){
    #referencing a distance matrix data object
    distout <- 0
    timediffs <- abs(difftime(times, t, units = "secs"))
    closesttime <- which(timediffs == min(timediffs, na.rm = T))[1]
    if(timediffs[closesttime] <= timesnap){
      tindex <- closesttime
    } else {
      stop(paste("Closest time is more than", timesnap, "seconds away"))
    }
    distout <- distmat[k, x, tindex]
    return(distout)
  }
  
  meshx_array <- as.array(1:nrow(mesh))
  Dx_pdotxs <- rep(0, length(meshx_array))
  for (meshx in meshx_array){
    Dx <- as.numeric(D_mesh[meshx])
    surv_eachtrap <- rep(0, nrow(traps))
    for(trapk in 1:length(traps)){
      opentimeindx <- which(!is.na(colSums(distmat[trapk,,], na.rm = F)))
      topentime <- times[min(opentimeindx)]
      tclosetime <- times[max(opentimeindx)]
      surv_eachtrap[trapk] <- surv(topentime, tclosetime, timeincr, 
                                   x = meshx, k = trapk, 
                                   lambda0 = lambda0, sigma = sigma,
                                   distmat = distmat)
    }
    surv_alltraps <- prod(surv_eachtrap)
    pdot <- 1 - surv_alltraps
    Dx_pdotx <- Dx * pdot
    Dx_pdotxs[meshx] <- (Dx_pdotx)
  }
  lambdan_ <- mean(Dx_pdotxs)
  n <- nrow(capthist)
  integral_eachi <- rep(0, n)
  for (i in 1:n){
    DKprod_eachx <- rep(0, nrow(mesh))
    for (x in 1:nrow(mesh)){
      Sxhx_eachtrap <- rep(0, nrow(traps))
      for(trapk in 1:nrow(traps)){
                opentimeindx <- which(!is.na(colSums(distmat[trapk,,], na.rm = F)))
                starttime <- times[min(opentimeindx)]
                if(!is.na(capthist[i,trapk])){
                  dettime <- capthist[i, trapk]
                  Sx <- as.numeric(surv(starttime, dettime, timeincr, x, trapk, lambda0, sigma, distmat))
                  hx <- as.numeric(haz(dettime, x, trapk, lambda0, sigma, distmat))
                  Sxhxout <- as.numeric(Sx * hx)
                } else { #individual never detected by this trap
                  endtime <- times[max(opentimeindx)]
                  Sx <- as.numeric(surv(starttime, endtime, timeincr, x, trapk, lambda0, sigma, distmat))
                  Sxhxout <- Sx
                }
              Sxhx_eachtrap[trapk] <- Sxhxout
              }
      Sxhx_alltraps <- prod(Sxhx_eachtrap)
      DKprod_out <- as.numeric(D_mesh[x]) * Sxhx_alltraps 
      DKprod_eachx[x] <- DKprod_out
    }
    integral_eachi[i] <- mean(DKprod_eachx) #sum * mesh area
  }
  ns <- 1:n
  logns <- log(ns)
  lognfact <- sum(logns)
  out <- (-lambdan_ - lognfact + n * sum(log(integral_eachi)))
  ADREPORT(sigma) ## If I want estimates on the real scale.
  ADREPORT(lambda0)
  ADREPORT(D_mesh)
  return(-out)
}

lambda0 = .9
sigma = 300
timeincrement = 6*60*10
meshgrid <- expand.grid(x = seq(000, 6000, 300), y = seq(000,6000, 300))
mesh <- make.mask(meshgrid, buffer = 0, spacing = 500)
D_mesh <- rep(.1, nrow(mesh))

#traps
#specify K moving traps
#traps should be dataframe with column for trapid, x, y, and time
traps <- rbind(
  data.frame(trapID = 1, x = seq(1500, 4500, 50), y = 1500, 
             time = seq(ymd_hms("2024-01-01 8:00:00"), ymd_hms("2024-01-01 14:00:00"), 
                        by = 360)),
  data.frame(trapID = 2, x = seq(1500, 4500, 50), y = 3500, 
             time = seq(ymd_hms("2024-01-02 8:00:00"), ymd_hms("2024-01-02 14:00:00"), 
                        by = 360)),
  data.frame(trapID = 3, x = seq(1500, 4500, 50), y = 2500, 
             time = seq(ymd_hms("2024-01-03 8:00:00"), ymd_hms("2024-01-03 14:00:00"), 
                        by = 360)),
  data.frame(trapID = 4, x = seq(1500, 4500, 50), y = 4500, 
             time = seq(ymd_hms("2024-01-01 9:00:00"), ymd_hms("2024-01-01 15:00:00"), 
                        by = 360))
)

dist_dat <- create_distdat(traps, mesh)
pop <- sim.popn(D = D_mesh, core = mesh, model2D = "IHP", 
                Ndist = "poisson", buffertype = "rect")
rownames(pop) <- NULL

capthist <- sim_capthist(pop, traps, timeincrement, lambda0, sigma, D_mesh)
dist_dat$capthist <- capthist
dist_dat$timeincr <- timeincrement

pars <- list(loglambda0 = as.numeric(0), logsigma = as.numeric(0), logD_mesh = as.numeric(rep(0, 144)))
myobj <- MakeADFun(func = negloglikelihood, parameters = pars)  ## First call is 100% R while turning into C++ and AD
fit <- nlminb(obj$par, obj$fn, obj$gr, control=list(iter.max=1000,eval.max=1000))

sdrep <- sdreport(obj)  ## Builds again to do ADREPORT values.
pl <- as.list(sdrep, "Est", report=TRUE)  ## Reported values.
plsd <- as.list(sdrep, "Std", report=TRUE)

pl$N
pl$sigma
pl$lambda


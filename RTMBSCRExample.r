## Build very basic SCR model in RTMB
library(RTMB)

dat <- list()
dat$traps <- expand.grid(x = 0:10, y = 0:10)
dat$range <- c(-3,13)
dat$mask <- expand.grid(x = seq(-dat$range[1], dat$range[2], by = 0.25),
  y = seq(-dat$range[1], dat$range[2], by = 0.25))
dat$dists2 <- t(apply(dat$mask, 1, FUN = function(x){(x[1]-dat$traps[,1])^2 + (x[2]-dat$traps[,2])^2}))

simSCR <- function(N = 50, sigma = 2, lambda = 5){
  dat$capt <<- NULL
  dat$ID <<- NULL
  S <- matrix(runif(N*2, dat$range[1], dat$range[2]), nrow = N, ncol = 2)
  J <- nrow(dat$traps)
  id <- 1
  for( i in 1:N ){
    d2 <- (S[i,1]- dat$traps[,1])^2 + (S[i,2] - dat$traps[,2])^2
    H <- lambda*exp(-0.5*d2/sigma^2)
    capts <- rpois(J, H)
    if(sum(capts) > 0){
      dat$capt <<- rbind(dat$capt, capts)
      dat$ID <<- c(dat$ID, id)
      id <- id + 1
    }
  }
  return(dat)
}

simSCR(50, 2, 5)

## Now build the model:

pars <- list(logLambda = 0, logSigma = 0)

SCRTMB <- function(pars){
  getAll(dat)
  lambda <- exp(pars$logLambda)
  sigma <- exp(pars$logSigma)
  m <- nrow(mask)
  J <- nrow(traps)
  K <- nrow(capt)
  ll <- 0
  
  Hij <- lambda*exp(-0.5*dists2/sigma^2)
  Hi <- rowSums(Hij)
  pDetect <- sum(1-exp(-Hi))/m
  logPDetect <- log(pDetect)
  for( k in 1:K){
    li <- 0
    for( i in 1:m ){
      li <- li + exp(sum(dpois(capt[k,], Hij[i,], log = TRUE)))
    }
    ll <- log(li) - logPDetect
  }
  ADREPORT(sigma) ## If I want estimates on the real scale.
  ADREPORT(lambda)
  N <- K/pDetect
  ADREPORT(N)
  return(-ll)
}

obj <- MakeADFun(SCRTMB, pars)  ## First call is 100% R while turning into C++ and AD

system.time(SCRTMB(pars))
system.time(obj$fn(pars)) ## WILD RIGHT????

fit <- nlminb(obj$par, obj$fn, obj$gr, control=list(iter.max=1000,eval.max=1000))

sdrep <- sdreport(obj)  ## Builds again to do ADREPORT values.
pl <- as.list(sdrep, "Est", report=TRUE)  ## Reported values.
plsd <- as.list(sdrep, "Std", report=TRUE)

## Got it right.
pl$N
pl$sigma
pl$lambda

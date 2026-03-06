
#compare fits between my script and secr
source("compare_moving_stat_2D/data_formatting/00_functions_and_parameters.R")
source("compare_moving_stat_2D/data_formatting/01_data_setup.R")
source("compare_moving_stat_2D/data_formatting/01.5_visualize.R")
source("compare_moving_stat_2D/data_formatting/02_simulate_and_fit.R")

#--------------Try calculating one llk value------------------------------------
#sim_fit <- function(sim){

simdch <- simulate_popandcapthist(tracksdf, D_mesh_v, lambda0, sigma,
                                  mesh, traps, trapspacing)
ch <- simdch$capthist
induse <- simdch$induse

mask <- mesh

covariates(mask) <- data.frame(cov = D_mesh_v,
                               x = mesh$x,
                               x2 = mesh$x^2)
trapscr <- read.traps(data = data.frame(TrapID = paste0("Trap", 1:nrow(traps)),
                                        x = traps$x,
                                        y = traps$y),
                      detector = "multi")
usage(trapscr) <- useall

rownames(trapscr) <- paste0("Trap", 1:nrow(traps))

eta = beta1*((mask$x/meshspacing + beta2)^2 )#+ (ys + beta2_)^2)
Z = sum(exp(eta)) * meshspacing^2/100^2
D = N * exp(eta) / Z

b0 = log(N/Z) 



#stepOrder is just the vector of last step where det happens
detinch <- which(ch == 1, arr.ind = T) #ikj
trap_n_vert <- ntrapsish/trap_n_horiz
CH <- make.capthist(captures = data.frame(
  Session = 1,
  ID = paste0("ind_", detinch[,1]),
  Occasion = detinch[,2],
  TrapID = paste0("Trap",detinch[,3])
),
traps = trapscr,
fmt = "trapID"
)

fit_secr <- secr.fit(CH, 
                     model = list(D ~ x + x2), 
                     detectfn = 'HHN',
                     mask = mask, 
                     start = c(b0,beta2,beta1,log(lambda0),log(sigma)))

start <- c(
  log(lambda0),
  log(sigma),
  beta1, 
  beta2, 
  log(N))

#quadratic density function
stat_nll <- function(v){
  lambda0_ <- exp(v[1])
  if (!is.finite(lambda0_)) return(1e12)
  
  sigma_ <- exp(v[2])
  if (!is.finite(sigma_)) return(1e12)
  
  D_mesh_  <- calcDv(scale(mesh[,1]),
                     mesh[,2],
                     v[3],
                     v[4],
                     exp(v[5]),
                     meshspacing)
  if (any(!is.finite(D_mesh_))) return(1e12)
  
  out <- negloglikelihood_stationary_cpp(lambda0_, sigma_,
                                         hazdenom, D_mesh_, 
                                         ch, useall,
                                         dist_trapmesh, as.matrix(mesh),
                                         mesharea = meshspacing^2)
  return(out)
}
wrapper <- function(v){
  out <- stat_nll(v)
  return(out$negloglik)
}

fitmy <- nlm(p = start,
              f = wrapper,
             print.level = 1,
              hessian = FALSE) 
diagmy <- statnll()

estmy <- fitmy$estimate
estsecr <- fit_secr$fit$par
mydf <- data.frame(lambda0 = exp(estmy[1]),
            sigma = exp(estmy[2]),
            beta1 = estmy[3],
            beta2 = estmy[4],
            N = exp(estmy[5])
)

secrdf <- data.frame(lambda0 = exp(estsecr[4]),
          sigma = exp(estsecr[5]),
          beta1 = estsecr[2],
          beta2 = estsecr[3],
          N = sum(exp(estsecr[1] + estsecr[2]*scale(mask$x) + estsecr[3]*scale(mask$x)^2)) * spacing(mask)^2/100^2)

truedf <- data.frame(lambda0 = lambda0,
                     sigma = sigma,
                     beta1 = beta1, 
                     beta2 = beta2,
                     N = N)
dat <-rbind(mydf, secrdf, truedf)
rownames(dat) <- c("my", "secr", "true")
dat



pdotsecr <- pdot(as.matrix(mesh), trapscr, detectfn = "HHN", 
     detectpar = list(lambda0 = lambda0,sigma =  sigma), 
     noccasions = dim(ch)[2])

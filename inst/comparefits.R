#compare fits between my script and Abinand's
source("compare_moving_stat_2D/data_formatting/00_functions_and_parameters.R")
source("compare_moving_stat_2D/data_formatting/01_data_setup.R")
#source("compare_moving_stat_2D/data_formatting/01.5_visualize.R")
source("compare_moving_stat_2D/data_formatting/02_simulate_and_fit.R")
source("inst/BananSim.R")
nsims = 30
trap_n_vert = round(ntrapsish/trap_n_horiz)
tracksteplength = round(trapspacing/nsteps_pertrap)

surv_obj <- setup_data(sigma,
                       N,
                       beta1,
                       beta2,
                       ntrapsish,
                       trackxmin,
                       trapspacing,
                       meshspacing,
                       trap_n_horiz,
                       nsteps_pertrap,
                       occreps
)
tracksdf <- surv_obj$tracksdf
traps <- surv_obj$traps
mesh <- surv_obj$mesh
D_mesh <- surv_obj$D_mesh_v
useall <- surv_obj$useall
dist_trapmesh <- surv_obj$dist_trapmesh
meanstepsize = mean(tracksdf$inc[tracksdf$inc !=0])

#create input objects for Abinand's
survObj <- simulateScrTrapsMask(
  nxTraps = trap_n_horiz,
  nyTraps = trap_n_vert,
  trapSpacing = trapspacing,
  nSteps = round(trapspacing/tracksteplength),
  maskSpacing = meshspacing,
  sigma = sigma, 
  N = N,
  b1 = beta1,
  b2 = beta2
)

#trapSteps <- survObj$trapSteps
trapSteps <- data.frame(matrix(nrow = nrow(tracksdf), ncol = 0))
trapSteps$x <- traps$x[tracksdf$trapno]
trapSteps$y <- traps$y[tracksdf$trapno]
trapSteps$xStep <- tracksdf$midx
trapSteps$yStep <- tracksdf$midy
trapSteps$effort <- tracksdf$inc/100
trapSteps$transect <- tracksdf$transect
trapSteps$TrapID <- paste0("Trap",tracksdf$trapno)
trapSteps <- trapSteps[which(tracksdf$rep == 1 & tracksdf$inc != 0),]
trapSteps$trapOrder <- rep((1:(trap_n_horiz * trapspacing/tracksteplength)),2)
trapSteps$stepOrder <- rep(c(1:10), (nrow(trapSteps)/10))
trapSteps$StepID <- paste0("Step", rownames(trapSteps))

mask <- mesh
covariates(mask) <- data.frame(cov = D_mesh,
                               x = mesh$x,
                               x2 = mesh$x^2/meshspacing^2)
saveRDS(mask, "inst/mask.RDS")

#do_a_sim_fit <- function(x){

pop <- sim_pop_C(D_mesh, 
                 as.matrix(mesh), 
                 meshspacing)

capthist_full <- sim_capthist_C(as.matrix(traps),
                                tracksdf, 
                                lambda0,
                                sigma,
                                D_mesh,
                                as.matrix(mesh),
                                meshspacing,
                                hazdenom=1,
                                pop,
                                dist_dat_pop = NULL,
                                report_probseenxk = F)

#standard scr likelihood
nocc <- length(unique(tracksdf$occ))
trapspacing <- sort(unique(traps$x))[2]- sort(unique(traps$x))[1]

if(length(which(apply((!is.na(capthist_full)), 1, sum)>0)) == 0){
  warning("Empty capture history.")
  capthist <- array(NA, dim = c(1, nocc, nrow(traps)))
} else {
  capthist <- capthist_full[which(apply((!is.na(capthist_full)), 1, sum)>0),,]
}

#use
induse <- create_ind_use_C(capthist,
                           as.matrix(traps),
                           trapspacing, 
                           tracksdf,
                           scenario = "everything")
saveRDS(induse, "inst/induse.Rds")
saveRDS(tracksdf, "inst/tracksdf.RDS")

#convert capthist to 1s and 0s
capthist[is.na(capthist)] <- 0
capthist[capthist!=0] <- 1

ch <- capthist
saveRDS(ch, "inst/ch.RDS")

trapscr <- read.traps(data = data.frame(TrapID = paste0("Trap", 1:nrow(traps)),
                        x = traps$x,
                        y = traps$y),
                        detector = "proximity")
rownames(trapscr) <- paste0("Trap", 1:nrow(traps))
saveRDS(trapscr, "inst/trapscr.RDS")


eta = beta1*((mask$x/meshspacing + beta2)^2 )#+ (ys + beta2_)^2)
Z = sum(exp(eta)) * meshspacing^2/100^2
D = N * exp(eta) / Z
  
b0 = log(N/Z) 
  
#stepOrder is just the vector of last step where det happens
detinch <- which(ch == 1, arr.ind = T) #ikj
CH <- make.capthist(captures = data.frame(
                          Session = 1,
                          ID = paste0("ind_", detinch[,1]),
                         Occasion = c(ifelse(detinch[,2] %% 2 ==0,1,2)), #so occasions 2 and 4 are 1, 1 and 3 are 2
                         TrapID = paste0("Trap",detinch[,3])
    ),
    traps = trapscr,
    fmt = "trapID"
  )
  
detinch <- data.frame(detinch) %>%
    arrange(as.character(dim1), dim3, dim2)
  
stepOrder <- apply(as.array(1:nrow(detinch)), 1, function(r){
    induse[detinch[r,1], detinch[r,3], detinch[r,2]]/10
  })
  
detinCH <- which(CH==1, arr.ind =T)
#reorder CH to put ind_10 at the end
CHorderi <- as.numeric(sapply(str_match_all(rownames(CH), "(ind_)(\\d+)"), function(x)x[,3]))
  
  
  # fit_abinand <- scrFitMov(capthist = CH, 
  #                          stepOrder = stepOrder,
  #                          mask = mask, 
  #                          trapSteps = trapSteps,
  #                          model = list(D ~ x + x2),
  #                          startparams = c(b0,beta2,beta1,log(lambda0),log(sigma)),
  #                          hessian = F)
  # 
set.seed(1)
simdch <- simulate_popandcapthist(tracksdf, D_mesh_v, lambda0, sigma,
                                  mesh, traps, trapspacing)
calc_nll(dist_trapmesh,
           useall,
           lambda0, 
           sigma, 
           D_mesh, 
           beta1, 
           beta2,
           N,
           hazdenom, 
           mesh, 
           capthistout = simdch,
           Dmod = "~x^2",
           meshspacing,
           meanstepsize)
set.seed(1)
  sim_fit(traps, tracksdf, mesh, meshspacing, dist_trapmesh, useall, lambda0, sigma, D_mesh_v, beta1, beta2, N, 1, Dmod = "~x^2")
  
  start0 <- c(
    log(lambda0),
    log(sigma),
    beta1, 
    beta2, 
    log(N))
  scaling_factors <- rep(1, length(start0)) #10^round(log10(abs(start0)))
  start <- start0/scaling_factors
  
  #moving detector likelihood
  nll <- function(v_scaled){
    v <- v_scaled * 1
    lambda0_ <- exp(v[1])
    if (!is.finite(lambda0_)) return(1e12)
    
    sigma_ <- exp(v[2])
    if (!is.finite(sigma_)) return(1e12)
    D_mesh_ <- calcDv(mesh[,1] ,
                      mesh[,2],
                      v[3],
                      v[4],
                      exp(v[5]),
                      meshspacing
    )
    
    out <- negloglikelihood_moving_cpp(lambda0 = lambda0_, 
                                       sigma = sigma_,  
                                       haz_denom = 1,
                                       D_mesh = D_mesh_,
                                       capthist = simdch$capthist, 
                                       usage = useall,
                                       indusage = simdch$induse, 
                                       distmat = dist_trapmesh,
                                       mesh = as.matrix(mesh),
                                       mesharea = meshspacing^2,
                                       meanstepsize = meanstepsize)
    
    return(out)
  }
  llkme2 <- nll(start)
  
  
  fit_me <- optim(par = start,
                  fn = nll,
                  hessian = F, method = "Nelder-Mead")
  est_abinand <- fit_abinand$estimate
  est_abinand <- c("b0" = est_abinand[1], 
                   "b2" = est_abinand[2],
                   "b1" = est_abinand[3],
                   "lambda0" = exp(est_abinand[4])/100, 
                   "sigma" = exp(est_abinand[5]))
  
  est_me <- fit_me$par
  est_me <- c("lambda0" = exp(est_me[1]), 
              "sigma" = exp(est_me[2]),
              "beta1" = est_me[3],
              "beta2" = est_me[4],
              "N" = exp(est_me[5]))
  
  true <- c(lambda0 = lambda0,
            sigma = sigma,
            beta1 = beta1,
            beta2 = beta2,
            N = N,
            beta0 = b0
  )
  
  out_ls <- list(est_me = est_me,
                 est_abinand = est_abinand,
                 sim = x)
#}
fit <- do_a_sim_fit()
fits <- lapply(as.list(1:nsims), do_a_sim_fit)
saveRDS(fits, "inst/comparefits_30.Rds")
fits <- readRDS("inst/comparefits_30.Rds")

beta1s <- data.frame(which = c(rep("my", length(fits)),
                                 rep("abi", length(fits))),
                     beta1est = c(sapply(fits, function(x){x$est_me["beta1"]}),
                               sapply(fits, function(x){x$est_abinand["b2"]})),
                     beta2est = c(sapply(fits, function(x){x$est_me["beta2"]}),
                               sapply(fits, function(x){x$est_abinand["b1"]})),
                     lambda0est = c(sapply(fits, function(x){x$est_me["lambda0"]}),
                               sapply(fits, function(x){x$est_abinand["lambda0"]})),
                     sigmaest = c(sapply(fits, function(x){x$est_me["sigma"]}),
                               sapply(fits, function(x){x$est_abinand["sigma"]})),
                     Nest = c(sapply(fits, function(x){x$est_me["N"]}),
                              sapply(fits, function(x){
                                sum(exp(x$est_abinand["b0"] + x$est_abinand["b2"]*(mask$x) + 
                                          x$est_abinand["b1"]*(mask$x/meshspacing)^2)) * spacing(mask)^2/100^2})
                              )
                     )



p1 <- ggplot() + 
  geom_boxplot(beta1s, 
               mapping = aes(y = (beta1est - beta1)/beta1*100,
                                     group = which,
                             fill = which))


p2 <- ggplot() + 
  geom_boxplot(beta1s,
               mapping = aes(y = beta2est,
                                     group = which,
                             fill = which))

p3 <- ggplot() + 
  geom_boxplot(beta1s,
               mapping = aes(y = (lambda0est-lambda0)/lambda0*100,
                             group = which,
                             fill = which))

p4 <- ggplot() + 
  geom_boxplot(beta1s,
               mapping = aes(y = (sigmaest-sigma)/sigma*100,
                             group = which,
                             fill = which))
p5 <- ggplot() + 
  geom_boxplot(beta1s,
               mapping = aes(y = (Nest-N)/N*100,
                             group = which,
                             fill = which)) 

grid.arrange(
  grobs = list(p1,p2,p3,p4,p5),
  widths = c(1,1),
  heights = c(1,1,1),
  layout_matrix = rbind(c(1,2),
                        c(3,4),
                        c(5,5)))

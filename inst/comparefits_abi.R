rm(list = ls())

source("compare_moving_stat_2D/data_formatting/00_functions_and_parameters.R")
source("compare_moving_stat_2D/data_formatting/01_data_setup.R")
source("compare_moving_stat_2D/data_formatting/02_simulate_and_fit.R")
source("inst/BananSim.R")

nsims = 30
ntrapsish = 150
trap_n_vert = round(ntrapsish/trap_n_horiz)
N = 69
meshspacing = 100
nsteps_pertrap = 5
trapspacing = 100
occreps = 2

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
meanstepsize = mean(tracksdf$inc[tracksdf$inc !=0])


survObj <- simulateScrTrapsMask(
  nxTraps = trap_n_horiz,
  nyTraps = trap_n_vert,
  trapSpacing = trapspacing,
  nSteps = nsteps_pertrap,
  maskSpacing = meshspacing,
  sigma = sigma, 
  N = N,
  b1 = beta1,
  b2 = beta2
)

trapSteps <- survObj$trapSteps
trapSteps <- data.frame(matrix(nrow = nrow(tracksdf), ncol = 0))
trapSteps$x <- traps$x[tracksdf$trapno]
trapSteps$y <- traps$y[tracksdf$trapno]
trapSteps$xStep <- tracksdf$midx
trapSteps$yStep <- tracksdf$midy
trapSteps$effort <- tracksdf$inc/100
trapSteps$transect = tracksdf$transect
trapSteps$TrapID <- paste0("Trap",tracksdf$trapno)
trapSteps <- trapSteps[which(tracksdf$rep == 1 & tracksdf$inc != 0),]
trapSteps$trapOrder <- rep((1:(trap_n_horiz * nsteps_pertrap)),nsteps_pertrap)
trapSteps$stepOrder <- rep(c(1:nsteps_pertrap), (nrow(trapSteps)/nsteps_pertrap))
trapSteps$StepID <- paste0("Step", rownames(trapSteps))

mask <- mesh

covariates(mask) <- data.frame(cov = D_mesh,
                               x = mesh$x,
                               x2 = mesh$x^2/meshspacing^2)


fit_secr <- fit_abi <- fit_sav <- list()

for (i in seq_len(nsims)){
  
  ch_out <- simulate_popandcapthist(tracksdf = tracksdf,
                                    D_mesh = D_mesh,
                                    lambda0 = lambda0, 
                                    sigma = sigma, 
                                    mesh = mesh,
                                    traps = traps,trapspacing = trapspacing)
  
  
  
  
  ch <- ch_out$capthist
  
  induse <- ch_out$induse
  
  
  trapscr <- read.traps(data = data.frame(TrapID = paste0("Trap", 1:nrow(traps)),
                                          x = traps$x,
                                          y = traps$y),
                        detector = "proximity")
  
  rownames(trapscr) <- paste0("Trap", 1:nrow(traps))
  
  
  
  eta = beta1*((mask$x/meshspacing + beta2)^2 )#+ (ys + beta2_)^2)
  Z = sum(exp(eta)) * meshspacing^2/100^2
  D = N * exp(eta) / Z
  
  b0 = log(N/Z) 
  
  
  
  #stepOrder is just the vector of last step where det happens
  detinch <- which(ch == 1, arr.ind = T) #ikj
  
  CH <- make.capthist(captures = data.frame(
    Session = 1,
    ID = paste0("ind_", detinch[,1]),
    Occasion = ((detinch[,2] - 1) %/% trap_n_vert) + 1,
    TrapID = paste0("Trap",detinch[,3])
  ),
  traps = trapscr,
  fmt = "trapID"
  )
  
  
  
   # fit_secr[[i]] <- secr.fit(CH,mask = mask, model = list(D ~ x + x2), detectfn = 'HHN',start = c(b0,beta2,beta1,log(lambda0*100),log(sigma)))
  
  
  detinch <- data.frame(detinch) %>% 
    arrange(as.character(dim1),dim2,dim3) 
  
  stepOrder <- apply(as.array(1:nrow(detinch)), 1, function(r){
    induse[detinch[r,1], detinch[r,3], detinch[r,2]]/(trapspacing/nsteps_pertrap)
  })

  fit_abi[[i]] <- scrFitMov(capthist = CH, 
                           stepOrder = stepOrder,
                           mask = mask, 
                           trapSteps = trapSteps,
                           model = list(D ~ x + x2),
                           startparams = c(b0,beta2,beta1,log(lambda0*100),log(sigma)),
                           hessian = F)
  
  
  
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
    sigma_ <- exp(v[2])
    D_mesh_ <- calcDv(mesh$x,
                      mesh$y,
                      v[3], 
                      v[4],
                      exp(v[5]),
                      meshspacing)
    
    out <- negloglikelihood_moving_cpp(lambda0 = lambda0_, 
                                       sigma = sigma_,  
                                       haz_denom = 1,
                                       D_mesh = D_mesh_,
                                       capthist = ch, 
                                       usage = surv_obj$useall,
                                       indusage = induse, 
                                       distmat = surv_obj$dist_trapmesh,
                                       mesh = as.matrix(mesh),
                                       mesharea = meshspacing^2,
                                       meanstepsize = meanstepsize)
    return(out)
  }
  


  fit_sav[[i]] <- nlm(f = nll,p = start, hessian = F,print.level = 1)
  
  
}

# save(fit_abi,fit_sav,fit_secr,file = 'comparefits_30.RData')


est_abi <- t(sapply(fit_abi, \(x) x$estimate)) %>% 
  data.frame() %>%
  setNames(c('b0','b1','b2','logLambda0','logSigma')) %>%
  rowwise() %>%
  mutate(lambda0 = exp(logLambda0),
         sigma = exp(logSigma),
         N = sum(exp(b0 + b1*covariates(mask)$x + b2*covariates(mask)$x2)) * spacing(mask)^2/100^2,
         model = 'Abi')


est_sav <- t(sapply(fit_sav, \(x) x$estimate)) %>% 
  data.frame() %>%
  setNames(c('logLambda0','logSigma','b1','b2','N')) %>%
  rowwise() %>%
  mutate(lambda0 = exp(logLambda0)*100,
         sigma = exp(logSigma),
         N = exp(N),
         model = 'Sav')

est_scr <- t(sapply(fit_secr, \(x) x$fit$estimate)) %>%
  data.frame() %>%
  setNames(c('b0','b1','b2','logLambda0','logSigma')) %>%
  rowwise() %>%
  mutate(lambda0 = exp(logLambda0)*4,
         sigma = exp(logSigma),
         N = sum(exp(b0 + b1*scale(mask$x) + b2*scale(mask$x)^2)) * spacing(mask)^2/100^2,
         model = 'Stationary')

expN <- N
trueSigma = 300
trueLambda0 = 0.8
allEst <- bind_rows(est_sav,est_abi,est_scr)

p1 <- allEst %>%
  ggplot()+
  geom_boxplot(aes(y = N, x = model, fill = model))+
  # geom_hline(yintercept = expN)+
  geom_hline(aes(yintercept = expN,col = model,group = model),alpha = 0.5) + 
  theme_bw()

p2 <- allEst %>%
  ggplot()+
  geom_boxplot(aes(y = sigma/trueSigma*100 - 100, x = model, fill = model))+
  geom_hline(yintercept = 0) + 
  ylab('sigma bias (%)')+
  theme_bw()

p3 <- allEst %>%
  ggplot()+
  geom_boxplot(aes(y = lambda0/trueLambda0*100 - 100, x = model, fill = model))+
  geom_hline(yintercept = trueLambda0)+
  ylab('lambda0 bias (%)')+
  theme_bw()


x1 <- unique(mask$x)

xlinD <- cbind(1,x1,(x1/spacing(mask))^2)

xlinD_sav <- cbind(1,x1,(x1/spacing(mask))^2)

xlinDs <- unique(fit_secr[[1]]$designD)

trueD = exp(b0 + beta2^2*beta1 + 2*beta1*beta2*x1/spacing(mask) + beta1*x1^2/spacing(mask)^2)

movD_abi <- exp(xlinD %*% t(as.matrix(est_abi[,1:3]))) %>%
  apply(.,1, quantile,c(0.025,0.5,0.975)) %>%
  t() %>%
  data.frame() %>%
  setNames(c('lcl','est','ucl')) %>%
  mutate(mod = 'abi',x = x1, trueD = trueD)




movD_sav <- sapply(1:nrow(est_sav), \(x) calcDv(x1,
                                                0,
                                                est_sav$b1[x], 
                                                est_sav$b2[x],
                                                est_sav$N[x],
                                                meshspacing)*300) %>% 
  apply(.,1, quantile,c(0.025,0.5,0.975)) %>%
  t() %>%
  data.frame() %>%
  setNames(c('lcl','est','ucl')) %>%
  mutate(mod = 'sav',x = x1, trueD = trueD)


statD <- exp(xlinDs %*% t(as.matrix(est_scr[,1:3]))) %>%
  apply(.,1, quantile,c(0.025,0.5,0.975)) %>%
  t() %>%
  data.frame() %>%
  setNames(c('lcl','est','ucl')) %>%
  mutate(mod = 'stationary',x = x1, trueD = trueD)

p4 <- bind_rows(movD_abi,movD_sav,statD) %>%
  ggplot(aes(x = x))+
  geom_line(aes(y = est, col = mod))+
  geom_line(aes(y = lcl, col = mod), linetype = 2)+
  geom_line(aes(y = ucl, col = mod), linetype = 2)+
  geom_line(aes(y = trueD))+
  theme_bw()

library(ggpubr)


ggarrange(ggarrange(p1,p2,p3,nrow = 1,common.legend = T),p4,
           nrow = 2)
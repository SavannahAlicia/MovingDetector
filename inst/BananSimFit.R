rm(list = ls())


source('inst/BananSim.R')


sigma = 300
trapSpacing = round(sigma/3)
maskSpacing = trapSpacing/2
nxTraps = 15
nyTraps = 13
nSteps = 10
nOccasionsTransect = 6
nsims = 30
trapxmin = -1000

N = 65
b1 = -0.015
b2 = 0
lambda0 = .008*100



survObj <- simulateScrTrapsMask(
  nxTraps = nxTraps,
  nyTraps = nyTraps,
  trapSpacing = trapSpacing,
  nSteps = nSteps,
  maskSpacing = maskSpacing,
  sigma = sigma, 
  N = N,
  b1 = b1,
  b2 = b2,
  trapxmin = trapxmin
)

trapSteps <- survObj$trapSteps

mask <- survObj$mask

#expN <- sum(exp(b0 + b1*covariates(mask)$cov))*spacing(mask)^2/100^2

eta = b1*((mask$x/spacing(mask) + b2)^2 )#+ (ys + beta2_)^2)
Z = sum(exp(eta)) * maskSpacing^2/100^2
D = N * exp(eta) / Z

b0 = log(N/Z) 

sum(exp(b0 + b2*mask$x + b1*mask$x^2/spacing(mask)^2))*spacing(mask)^2/100^2

simPop = replicate(nsims, sim.popn(D = D, core = mask,model2D = 'IHP'),simplify = F)

ggplot()+
  geom_tile(data = mask,aes(x = x, y = y, fill = covariates(mask)$cov))+
  geom_point(data = trapSteps, aes(x = x, y = y), col = 'firebrick')+
  geom_point(data = trapSteps, aes(x = xStep, y = yStep), col = 'grey80',shape = 3)+
  geom_point(data = simPop[[1]], aes(x = x, y = y))+
  coord_equal()+
  scale_fill_viridis_c(name = 'cov')

allCapthists <- lapply(simPop, \(p) simCapthist(pop = p,
                                                trapSteps = trapSteps, 
                                                mask = mask, 
                                                lambda0 = lambda0, 
                                                sigma = sigma,
                                                nOccasionsTransect = nOccasionsTransect))

capthist <- allCapthists[[1]]

sum(capthist)

traps <- traps(capthist)

usage(traps)

ggplot(traps)+
  geom_tile(aes(x = x, y = y, fill = apply(capthist,3,sum)))+
  scale_color_viridis_c()+
  coord_equal()


movFit <- lapply(allCapthists, \(x) scrFitMov(capthist = x, mask = mask, 
                                              trapSteps = trapSteps,
                                              model = list(D ~ x + x2),
                                              startparams = c(b0,b2,b1,log(lambda0),log(sigma)),
                                              hessian = F))
statFit <- lapply(allCapthists, \(x) secr.fit(capthist = x,
                                              mask = mask, 
                                              model = list(D ~ x + x2),
                                              detectfn = 'HHN',
                                              start = c(b0,b2,b1,log(lambda0),log(sigma))))

expN = N
trueSigma = sigma
trueLambda0 = lambda0
trueb0 = b0
trueb1 = b1
trueb2 = b2

movEstimates <- t(sapply(movFit, \(x) x$estimate)) %>% 
  data.frame() %>% 
  setNames(c('b0','b1','b2','logLambda0','logSigma')) %>% 
  rowwise() %>% 
  mutate(lambda0 = exp(logLambda0),
         sigma = exp(logSigma),
         N = sum(exp(b0 + b1*covariates(mask)$x + b2*covariates(mask)$x2)) * spacing(mask)^2/100^2,
         model = 'Moving') 

statEstimates <- t(sapply(statFit, \(x) x$fit$estimate)) %>% 
  data.frame() %>% 
  setNames(c('b0','b1','b2','logLambda0','logSigma')) %>% 
  rowwise() %>% 
  mutate(lambda0 = exp(logLambda0),
         sigma = exp(logSigma),
         N = sum(exp(b0 + b1*scale(mask$x) + b2*scale(mask$x)^2)) * spacing(mask)^2/100^2,
         model = 'Stationary') 


allEst <- bind_rows(movEstimates,statEstimates)

allEst %>% 
  filter(N<200) %>% 
  ggplot()+
  geom_boxplot(aes(y = N, x = model))+
  # geom_hline(yintercept = expN)+
  geom_hline(aes(yintercept = expN,col = model,group = model),alpha = 0.5)

allEst %>% group_by(model) %>% 
  summarise(mean(N))

allEst %>% 
  ggplot()+
  geom_boxplot(aes(y = sigma/trueSigma*100 - 100, x = model))+
  geom_hline(yintercept = 0)
  
allEst %>% 
  ggplot()+
  geom_boxplot(aes(y = logLambda0/log(trueLambda0)*100 - 100, x = model))+
  geom_hline(yintercept = 0)


x1 <- unique(mask$x)

xlinD <- cbind(1,x1,(x1/spacing(mask))^2)

xlinDs <- unique(statFit[[1]]$designD)

trueD = exp(trueb0 + trueb2*x1 + trueb1*(x1/spacing(mask))^2)

movD <- exp(xlinD %*% t(as.matrix(movEstimates[,1:3]))) %>% 
  apply(.,1, quantile,c(0.025,0.5,0.975)) %>% 
  t() %>% 
  data.frame() %>% 
  setNames(c('lcl','est','ucl')) %>% 
  mutate(mod = 'moving',x = x1, trueD = trueD)

statD <- exp(xlinDs %*% t(as.matrix(statEstimates[,1:3]))) %>% 
  apply(.,1, quantile,c(0.025,0.5,0.975)) %>% 
  t() %>% 
  data.frame() %>% 
  setNames(c('lcl','est','ucl')) %>% 
  mutate(mod = 'stationary',x = x1, trueD = trueD)

bind_rows(movD,statD) %>% 
  ggplot(aes(x = x))+
  geom_line(aes(y = est, col = mod))+
  geom_line(aes(y = lcl, col = mod), linetype = 2)+
  geom_line(aes(y = ucl, col = mod), linetype = 2)+
  geom_line(aes(y = trueD))+
  theme_bw()
  




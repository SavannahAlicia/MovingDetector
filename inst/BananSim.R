library(tidyverse)
library(secr)
library(rlang)


createSteps <- function(x,y,trapSpacing,nSteps){
  xStep = seq(x - trapSpacing/2+ trapSpacing/(2*nSteps), x + trapSpacing/2 - trapSpacing/(2*nSteps),length.out = nSteps)
  data.frame(x = x, y = y, xStep = xStep, yStep = y,effort = trapSpacing/nSteps/100)
}

simulateScrTrapsMask <- function(
    nxTraps = 21,
    nyTraps = 10,
    nSteps = 3,
    maskSpacing = 50,        #spacing between mask points
    sigma = 300,
    trapSpacing = NULL,   # Spacing between traps in grid
    N = 100,
    b1 = 0,
    b2 = 0){  # Spacial covariate smooth factor
  
  #Define spacing between traps as is recommended
  if(is.null(trapSpacing)) trapSpacing = sigma
  
  #Create trapping grid
  rawTrap = expand.grid(x = seq( - ((nxTraps-1)*trapSpacing/2), ((nxTraps-1)*trapSpacing/2),length.out = nxTraps),
                        y = seq(- ((nyTraps-1)*trapSpacing/2),((nyTraps-1)*trapSpacing/2),length.out = nyTraps)) %>%
    arrange(y,x) %>%
    mutate(transect = as.numeric(factor(y))) %>%
    mutate(TrapID = paste0('Trap',rownames(.)))
  
  trapSteps <- map(1:nrow(rawTrap), \(i) createSteps(rawTrap$x[i],rawTrap$y[i],trapSpacing,nSteps)) %>%
    bind_rows() %>%
    left_join(rawTrap) %>%
    arrange(yStep,xStep) %>%
    group_by(transect) %>%
    mutate(trapOrder = 1:n()) %>% 
    group_by(TrapID) %>% 
    mutate(stepOrder = 1:n()) %>%
    ungroup()  %>%
    mutate(StepID = paste0('Step',rownames(.)))
  
  
  ## Create trap object
  traps = read.traps(data = rawTrap %>% dplyr::select(x,y), detector = 'proximity')
  
  ## Create mask object at recommended buffer of 4*sigma
  mask = make.mask(traps, buffer = 4*sigma, spacing = maskSpacing)
  
  #
  # distmat = as.matrix(dist(mask))
  #
  # V = exp(-distmat/(maskSpacing*10))
  #
  # run = runif(nrow(mask),min = -1.5,max = 1.5)
  #
  # cov = t(chol(V))%*%run
  #
  # rm(distmat,V)
  #
  # covariates(mask)$cov = cov
  
  eta = b1*((mask$x/spacing(mask) + b2)^2)
  Z = sum(exp(eta)) * maskSpacing^2
  D = N * exp(eta) / Z
  
  ggplot(mask)+
    geom_tile(aes(x = x, y= y, fill= D))+
    scale_fill_viridis_c()+
    coord_equal()
  
  covariates(mask) = data.frame(cov = D, x= mask$x, x2 = mask$x^2/spacing(mask)^2)
  
  if(is.null(spacing)) spacing = spacing(mask)
  
  ## Create spatial covariate and smooth with a gaussian kernel
  
  return(list(mask = mask,trapSteps = trapSteps))
}




simCapthist = function(pop,trapSteps,mask,lambda0,sigma,nOccasionsTransect){
  
  ninds = nrow(pop)
  
  traps = trapSteps %>%
    dplyr::select(x,y,TrapID) %>%
    unique() %>%
    data.frame(.,row.names = .$TrapID) %>%
    read.traps(data = ., detector = 'proximity')
  
  
  effort <- trapSteps %>% group_by(TrapID) %>%
    summarise(usage = sum(effort)) %>% pull(usage)
  
  usage(traps) <- replicate(nOccasionsTransect,effort)
  
  transects = unique(trapSteps$transect)
  
  steps = read.traps(data = trapSteps %>% dplyr::select(xStep,yStep) %>%
                       setNames(c('x','y')) %>%
                       data.frame(.,row.names = trapSteps$StepID), detector = 'proximity')
  
  usage(steps) <- replicate(nOccasionsTransect,trapSteps$effort)
  
  capthistStat <- sim.capthist(traps = steps, popn = pop, detectfn = 'HHN', detectpar = list(lambda0 = lambda0, sigma = sigma), noccasions = nOccasionsTransect)
  
  dat <- data.frame(capthistStat) %>%
    rename(StepID = TrapID) %>%
    left_join(trapSteps) %>%
    ungroup() %>%
    group_by(ID,Occasion,transect) %>%
    arrange(trapOrder) %>%
    slice(1) %>%
    ungroup() %>%
    select(Session,ID,Occasion,TrapID,stepOrder) %>%
    mutate(ID = paste0('ind_',ID)) %>%
    data.frame()
  
  capthist = dat %>%
    make.capthist(.,traps = traps,fmt = 'trapID')
  
  stepOrder = dat %>% 
    pull(stepOrder)
  
  return(list(capthist = capthist, stepOrder = stepOrder))
  
}

##### Half normal encounter rate function

lambda_hhn <- function(d,lambda0,sigma){
  return(lambda0*exp(-d^2/(2*(sigma^2))))
}


modifymodel <- function(model, user_model){
  modifyList(model,user_model[intersect(names(model),names(user_model))])
}

lik_opt <- function(params,
                    desmat,
                    args.index,
                    distmat,
                    effort, # list of length nverttraps, each length n horiz traps
                    nOccasion, #occasion reps
                    dets,
                    inds,
                    n,
                    lambda_x,#hazard function 
                    noDets,
                    a){
  D = exp(desmat$D %*% params[args.index$D])
  lambda0 = exp(params[args.index$lambda0])
  sigma = exp(params[args.index$sigma])
  
  #hu list length ntransects, matrix of ntrap (per transect) by nmesh 
  lambda_x_mat <- lapply(1:length(distmat),\(x) t(t(lambda_x(distmat[[x]],lambda0,sigma)) %*% diag(effort[[x]])))
  
  ntotsurveys <- nOccasion * length(effort)
  nhoriz <- length(effort[[1]])
  nvert <- length(effort)
  
  #sum of hu nmesh by ntransect
  transectLambda <- sapply(lambda_x_mat,colSums)
  notseen_log <- t(do.call(rbind, replicate(5, t(-transectLambda), simplify = F)))
  

  prob_no_capt <-  exp(-rowSums(transectLambda*nOccasion)) ### Probability of 0 encounters
  
  prob_capt <- 1 - prob_no_capt ### probability of atleast 2 encounters
  Dx_pdotxs <- prob_capt * D
  
  llk_inner <- function(i){
    ind_dat <- dets[[i]]
    DKprod_eachx_log <- (
      rowSums(sapply(seq_len(nrow(ind_dat)), \(j){
        colSums(rbind(log(1 - dpois(0,lambda_x_mat[[ind_dat$transect[j]]][ind_dat$order[j],]/10)),
                      dpois(0,lambda_x_mat[[ind_dat$transect[j]]][ind_dat$order[j],]*ind_dat$last[j]/10,log = T),
                      ifelse(ind_dat$order[j]>1,1,0)*dpois(0,lambda_x_mat[[ind_dat$transect[j]]][c(1:(ind_dat$order[j]-1)),],log = T)))
      })) - rowSums(transectLambda %*% diag(noDets[i,]))
    )*D
     inner <- log(sum(exp(DKprod_eachx_log))*a)
    out_ls <- list(DKprod_eachx_log = DKprod_eachx_log,
                   inner = inner)
  }
  
  DKprod_eachx_log <- sapply(inds, function(i){llk_inner(i)$DKprod_eachx_log})
  lik_all_ind <- sapply(inds,function(i){llk_inner(i)$inner})
  
  ### the full likelihood
  loglik = - sum(D * a * prob_capt) - lfactorial(n) + sum(lik_all_ind)
  out_ls <- list(negloglik = -loglike,
                  notseen_log,
                  integral_eachi_log,
                  sumallhuexcept,
                  didntsurvivej_log,
                  DKprod_eachx_log,
                  Dx_pdotxs)
  return(-loglik)
}


scrFitMov <- function(capthist,
                                 stepOrder,
                                 mask,
                                 trapSteps,
                                 model = NULL, 
                                 startparams = NULL,
                                 hessian = F){
  
  ## Define the standard formula
  correct_model = list(D~1,lambda0~1,sigma~1)
  ##Name the list with the formulae
  names(correct_model) <- lapply(correct_model, f_lhs)
  
  ## Use initial parameters if provided
  if(!is.null(model)) {
    names(model) <- lapply(model, f_lhs)
    model = modifymodel(correct_model,model)
  }else{
    model = correct_model
  }
  
  ## Remove single encounters
  
  lambda_x <- lambda_hhn  ### Hazard half normal encounter rate
  
  ## Create design matrix
  desmat <- lapply(model,model.matrix,data = 
                     data.frame(cbind(D = 1,lambda0 = 1,sigma = 1, covariates(mask))))
  
  ## Calculate number of parameters
  npars = lapply(desmat,ncol)
  
  ### extract number of individuals detected
  n <- nrow(capthist)
  
  inds <- rownames(capthist)
  
  ### Area of each pixel
  a <- attr(mask,'spacing')^2/100^2
  
  trapOrder <- trapSteps %>%
    select(x,y,TrapID,transect,effort) %>%
    group_by(TrapID, transect) %>%
    mutate(effort = sum(effort)) %>%
    unique() %>%
    group_by(transect) %>%
    mutate(order = 1:n())
  
  nOccasion <- max(occasion(capthist))
  
  capthistMov <- data.frame(capthist) %>%
    arrange(ID) %>% 
    mutate(last = stepOrder-1) %>% 
    left_join(trapOrder) %>%
    dplyr::select(ID,transect,order,last)
  
  noDets <- nOccasion - expand.grid(ID = unique(capthistMov$ID), transect = unique(trapOrder$transect)) %>%
    left_join(capthistMov %>% group_by(ID,transect) %>% summarise(dets = n())) %>%
    arrange(ID) %>%
    mutate(dets = ifelse(is.na(dets),0,dets)) %>%
    pivot_wider(names_from = ID, values_from = dets) %>%
    dplyr::select(-1) %>%
    as.matrix() %>% t()
  
  dets <- capthistMov %>%
    split(.,.$ID,drop = T)
  
  traps <- trapOrder %>%
    dplyr::select(x,y,transect) %>%
    split(.,.$transect) %>%
    map(\(t) data.frame(t[,1:2]))
  
  effort <- trapOrder %>%
    dplyr::select(effort,transect) %>%
    split(.,.$transect) %>%
    map(\(t) t$effort)
  ## distance from traps to each mask point
  distmat = lapply(traps, \(x) edist(x[,1:2],mask))
  
  
  ## Set up the parameter vector
  if(!is.null(startparams)){
    params = startparams
  }else{
    autopars = autoini(capthist,mask)
    D = numeric(length = npars$D)
    D[1] = log(autopars$D)
    lambda0 = numeric(length = npars$lambda0)
    lambda0[1] = log(1 - exp(-autopars$g0))
    sigma = numeric(length = npars$sigma)
    sigma[1] = log(autopars$sigma)
    params = c(D,lambda0,sigma)
  }
  
  ## Vector of the index of each parameter
  args.index = split(1:length(params),rep(c('D','lambda0','sigma'),times = c(npars$D,npars$lambda0,npars$sigma)))
  
  ## Local function to optimise
  optimiser <- local ({
   })
  
  args = list(f = optimiser,p = params,hessian = hessian, print.level = 1)
  
  ## Run optimiser
  mod = do.call(nlm,args)
  
  ## Output optimised parameters, lambda matrix and p.
  
  return(mod)
  
}


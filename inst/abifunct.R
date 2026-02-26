library(rlang)
library(tidyverse)

lambda_hhn <- function(d,lambda0,sigma){
  return(lambda0*exp(-d^2/(2*(sigma^2))))
}




modifymodel <- function(model, user_model){
  modifyList(model,user_model[intersect(names(model),names(user_model))])
}


prep_dat_for_lik <- function(trapSteps,
                             capthistscr
){
  
  trapOrder <- trapSteps %>%
    select(x,y,TrapID,transect,effort) %>%
    group_by(TrapID, transect) %>%
    mutate(effort = sum(effort)) %>%
    unique() %>%
    group_by(transect) %>%
    mutate(order = 1:n())
  
  nOccasion <- max(occasion(capthistscr))
  
  capthistMov <- data.frame(capthistscr) %>%
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
  distmat = lapply(traps, \(x) edist(x[,1:2],mask))
  outs <- list(dets = dets,
               noDets = noDets,
               effort = effort,
               distmat = distmat)
  
  return(outs)
  
}

lik_opt <- function(startparams,
                    capthistscr,
                    mask,
                    distmat,
                    nOccasion,
                    dets,
                    noDets,
                    effort,
                    model = NULL
){
  n <- nrow(capthistscr)
  a <- attr(mask,'spacing')^2/100^2
  
  ## Define the standard formula
  correct_model = list(D~1, lambda0~1, sigma~1)
  ##Name the list with the formulae
  names(correct_model) <- lapply(correct_model, f_lhs)
  
  ## Use initial parameters if provided
  if(!is.null(model)) {
    names(model) <- lapply(model, f_lhs)
    model = modifymodel(correct_model,model)
  }else{
    model = correct_model
  }
  
  ## Create design matrix
  desmat <- lapply(model, model.matrix,
                   data = data.frame(cbind(D = 1, 
                                           lambda0 = 1, 
                                           sigma = 1, 
                                           covariates(mask))))
  
  ## Calculate number of parameters
  npars = lapply(desmat,ncol)
  ## Set up the parameter vector
  if(!is.null(startparams)){
    params = startparams
  }else{
    autopars = autoini(capthistscr,mask)
    D = numeric(length = npars$D)
    D[1] = log(autopars$D)
    lambda0 = numeric(length = npars$lambda0)
    lambda0[1] = log(1 - exp(-autopars$g0))
    sigma = numeric(length = npars$sigma)
    sigma[1] = log(autopars$sigma)
    params = c(D,lambda0,sigma)
  }
  
  ## Vector of the index of each parameter
  args.index = split(1:length(params),
                     rep(c('D','lambda0','sigma'),
                         times = c(npars$D,npars$lambda0,npars$sigma)))
  
  D = exp(desmat$D %*% params[args.index$D])
  lambda0 = exp(params[args.index$lambda0])
  sigma = exp(params[args.index$sigma])
  
  #hu list length ntransects, matrix of ntrap (per transect) by nmesh 
  lambda_x_mat <- lapply(1:length(distmat),function(x){
    t(t(lambda_hhn(distmat[[x]],lambda0,sigma)) %*% diag(effort[[x]]))
  } )
  
  #sum of hu nmesh by ntransect
  transectLambda <- sapply(lambda_x_mat,colSums)
  notseen_log <- t(do.call(rbind, 
                           replicate(nOccasion, t(-transectLambda), simplify = F)))
  
  
  prob_no_capt <-  exp(-rowSums(transectLambda*nOccasion)) ### Probability of 0 encounters
  
  prob_capt <- 1 - prob_no_capt ### probability of atleast 2 encounters
  Dx_pdotxs <- prob_capt * D
  
  llk_inner <- function(i){
    ind_dat <- dets[[i]]
    Js <- seq_len(nrow(ind_dat))
    DKprod_eachx_log <- (
      rowSums(sapply(Js, function(j){
        tr_j <- ind_dat$transect[j]
        o_j <- ind_dat$order[j]
        l_j <- ind_dat$last[j]
        colSums(rbind(
          log(1 - dpois(0, lambda_x_mat[[tr_j]][o_j,]/10)),
          dpois(0,lambda_x_mat[[tr_j]][o_j,] * l_j/10,log = T),
          ifelse(o_j>1,1,0)*dpois(0,lambda_x_mat[[tr_j]][c(1:(o_j-1)),]/10,log = T)
        ))
      })) - rowSums(transectLambda %*% diag(noDets[i,]))
    ) * D
    inner <- log(sum(exp(DKprod_eachx_log)) * a)
    out_ls <- list(DKprod_eachx_log = DKprod_eachx_log,
                   inner = inner)
  }
  inds <- rownames(capthistscr)
  DKprod_eachx_log <- sapply(inds, function(i){llk_inner(i)$DKprod_eachx_log})
  integral_eachi_log <- sapply(inds,function(i){llk_inner(i)$inner})
  
  ### the full likelihood
  loglik =  -sum(D * a * prob_capt) - lfactorial(n) + sum(integral_eachi_log)
  out_ls <- list(negloglik = -loglik,
                 notseen_log = notseen_log,
                 integral_eachi_log = integral_eachi_log,
                 #sumallhuexcept,
                 #didntsurvivej_log,
                 DKprod_eachx_log = DKprod_eachx_log,
                 pdotxs = prob_capt)
  return(out_ls)
}





## Local function to optimise
optimiser <- local ({
  wrap <- function(params){
    like_opt(params)
  }
})



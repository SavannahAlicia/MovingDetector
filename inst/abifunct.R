##### Half normal encounter rate function

lambda_hhn <- function(d,lambda0,sigma){
  return(lambda0*exp(-d^2/(2*(sigma^2))))
}


modifymodel <- function(model, user_model){
  modifyList(model,user_model[intersect(names(model),names(user_model))])
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
    lik_opt <- function(params){
      D = exp(desmat$D %*% params[args.index$D])
      lambda0 = exp(params[args.index$lambda0])
      sigma = exp(params[args.index$sigma])
      
      lambda_x_mat <- lapply(1:length(distmat),\(x) t(t(lambda_x(distmat[[x]],lambda0,sigma)) %*% diag(effort[[x]])))
      
      ### 21*1092 (traps x mask)
      
      transectLambda <- sapply(lambda_x_mat,colSums)
      
      prob_no_capt <-  exp(-rowSums(transectLambda*nOccasion)) ### Probability of 0 encounters
      
      prob_capt <- 1 - prob_no_capt ### probability of atleast 2 encounters
      
      lik_all_ind <- sapply(inds,\(i){
        ind_dat <- dets[[i]]
        log(sum(exp(
          rowSums(sapply(seq_len(nrow(ind_dat)), \(j){
            colSums(rbind(log(1 - dpois(0,lambda_x_mat[[ind_dat$transect[j]]][ind_dat$order[j],]/10)),
                          dpois(0,lambda_x_mat[[ind_dat$transect[j]]][ind_dat$order[j],]*ind_dat$last[j]/10,log = T),
                          ifelse(ind_dat$order[j]>1,1,0)*dpois(0,lambda_x_mat[[ind_dat$transect[j]]][c(1:(ind_dat$order[j]-1)),],log = T)))
          })) - rowSums(transectLambda %*% diag(noDets[i,]))
        )*D*a))
      })
      
      ### the full likelihood
      loglik = - sum(D * a * prob_capt) - lfactorial(n) + sum(lik_all_ind)
      
      return(-loglik)
    }
  })
  
  #return(optimiser(params))
  
  args = list(f = optimiser,p = params,hessian = hessian, print.level = 1)
  
  ## Run optimiser
  mod = do.call(nlm,args)
  
  ## Output optimised parameters, lambda matrix and p.
  
  return(mod)
  
}




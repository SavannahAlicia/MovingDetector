
#compare fits between my script and Abinand's
source("compare_moving_stat_2D/data_formatting/00_functions_and_parameters.R")
source("compare_moving_stat_2D/data_formatting/01_data_setup.R")
source("compare_moving_stat_2D/data_formatting/01.5_visualize.R")
source("compare_moving_stat_2D/data_formatting/02_simulate_and_fit.R")
source("inst/abifunct.R")

#--------------Try calculating one llk value------------------------------------
simdch <- simulate_popandcapthist(tracksdf, D_mesh_v, lambda0, sigma,
                                  mesh, traps, trapspacing)
ch <- simdch$capthist
induse <- simdch$induse

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


detinch <- data.frame(detinch) %>% 
  arrange(as.character(dim1),dim2,dim3) 

stepOrder <- apply(as.array(1:nrow(detinch)), 1, function(r){
  induse[detinch[r,1], detinch[r,3], detinch[r,2]]/10
})


fit_abi <- scrFitMov(capthist = CH, 
                     stepOrder = stepOrder,
                     mask = mask, 
                     trapSteps = trapSteps,
                     model = list(D ~ x + x2),
                     startparams = c(b0,beta2,beta1,log(lambda0*100),log(sigma)),
                     hessian = F)

fit_my <- fit_capthist(dist_trapmesh,
                       useall,
                       lambda0,
                       sigma,
                       D_mesh = D_mesh_v,
                       beta1,
                       beta2,
                       N,
                       hazdenom,
                       mesh,
                       capthistout = simdch,
                       Dmod = "~x^2",
                       meshspacing,
                       meanstepsize = mean(tracksdf$inc[tracksdf$inc != 0]),
                       fitstat = FALSE
)

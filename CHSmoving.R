library(ggplot2)
library(sf)
library(parallel)
library(sp)
library(openpopscr)
library(secr)

#-------------------------------functions---------------------------------------
#' Create polygons for each box centered at grid pt
#' 
#' @param grid dataframe with x and y columns specifying grid pts
#' @param spacing the spacing between grid pts, if not specified will calculate
#' @return SpatialPolygons object in list form 
#' @export
create_grid_polygons <- function(grid, spacing = NULL){
  #so we only have to create gridbboxes and lines list once
  if(is.null(spacing)){
    spacing <- min(dist(grid)[dist(grid$x)!=0])
  }
  #specify corners of box centered on grid pt
  grid_spacing <- spacing
  grid_bboxes <- data.frame(
    left_bound = grid$x - (grid_spacing/2),
    right_bound = grid$x + (grid_spacing/2), 
    lower_bound = grid$y - (grid_spacing/2),
    upper_bound = grid$y + (grid_spacing/2)
  )
  #create a polygon for box centered on grid pt
  create_polygon_for_grid_row <- function(grid_bboxes_row){
    upperleft <- as.numeric(c(grid_bboxes_row[1], grid_bboxes_row[4]))
    lowerleft <- as.numeric(c(grid_bboxes_row[1], grid_bboxes_row[3]))
    upperright <- as.numeric(c(grid_bboxes_row[2], grid_bboxes_row[4]))
    lowerright <- as.numeric(c(grid_bboxes_row[2], grid_bboxes_row[3]))
    boxpts <- rbind(upperleft, upperright, lowerright, lowerleft)
    colnames(boxpts) <- c("x", "y")
    rownames(boxpts) <- NULL
    p = sp::Polygon(boxpts)
    ps = sp::Polygons(list(p),1)
    sps = sp::SpatialPolygons(list(ps))
    proj4string(sps) <- sp::CRS(sf::st_crs(26916)$input)
    return(sps)
  }
  #create polygon box for all grid pts
  grid_bboxes$polygons  <- apply(X = grid_bboxes, MARGIN = 1, 
                                 FUN = create_polygon_for_grid_row)
  return(grid_bboxes)
}

#' Create named list of SpatialLines where name is track ID (occ)
#' @param tracksdf dataframe with xy, occ, trapno, and time
#' @return named list of SpatialLines objects
#' @export
create_line_spatlines <- function(tracksdf, scenario = "onison", tracksdfcolname = "ID",
                                  projto = sp::CRS(sf::st_crs(26916)$input)){
  #create a line for between each pt of track, put into named list by ID
  #check that it's not the last point in a track
  #(we don't want a line between tracks)
  check_pt <- function(row, tracksdf. = tracksdf){
    if (tracksdf.[row, tracksdfcolname] == tracksdf.[row+1, tracksdfcolname]){
      newline <- Line(tracksdf.[c(row, row+1),c("x", "y")])
    } else {
      newline <- NULL
    }
    return(newline)
  }
  linelist <- lapply(X = as.list(1:(nrow(tracksdf)-1)), FUN = check_pt)
  names(linelist) <- tracksdf[1:(nrow(tracksdf)-1), tracksdfcolname]
  linelist <- linelist[which(!lapply(linelist, is.null) == TRUE)]
  big_Lines_list <- list()
  nocc <- length(unique(tracksdf[, tracksdfcolname]))
  for (trackindex in 1:nocc){
    trackIDi <- unique(tracksdf[, tracksdfcolname])[trackindex]
    thetrack <- tracksdf[which(tracksdf[, tracksdfcolname] == trackIDi),]
    thelines <- linelist[which(names(linelist) == paste(trackIDi))]
    
    #check scenario
    if (scenario == "onison"){ #on effort is on effort
      newlines <- list(Lines(thelines[which(
        thetrack[1:(nrow(thetrack)-1),]$effort == "OnEffort")], 
        ID = paste(trackIDi)))
    } else if (scenario == "everything"){ #everything is on effort
      newlines <- list(Lines(thelines, ID = paste(trackIDi)))
    }
    #just creating desired class of object
    newlines_spat <- SpatialLines(newlines)
    #specify projection
    proj4string(newlines_spat) <- projto
    #adding new objects to list
    big_Lines_list <- append(big_Lines_list, newlines_spat)
    #naming list for track ID
    names(big_Lines_list)[trackindex] <- paste(trackIDi)
  }
  return(big_Lines_list)
}

#' Get length line that cross polygon
#' @param polygon a spatial polygon object
#' @param sp_line a Spatial Lines object
#' @return returns the length of the line within polygon
#' @export
get_length <- function(polygon, sp_line){
  boxlength = 0
  #attempts to chop line to just the bit within the polygon
  theintersection <- st_intersection(st_as_sf(polygon), st_as_sf(sp_line))
  #only calculate length of the line if there is a line that exists
  if(length(theintersection$geometry) > 0){
    boxlength <- as.numeric(st_length(theintersection))
  } #otherwise the length remains 0
  return(boxlength)
}

#'Calculate track lengths in grid boxes
#' function to find intersection between lines and grid box 
#' @param tracks_sp_lines list of lines labelled by track ID in name field
#' @param grid_polygons dataframe with a polygon in each row corresponding to 
#'                      that grid pt
#' @return column of summed intersections, each row is grid pt (trap)
#'         and sum is for all tracks specified in function
lengths_in_grid <- function(tracks_sp_lines, tracksinoc, grid_polygons){
  trackspl_to_use <- tracks_sp_lines[which(names(tracks_sp_lines) %in%
                                             as.character(tracksinoc))]
  
  if (length(trackspl_to_use) == 1){#just one track per occ now 
    inters <-  sapply(X = grid_polygons$polygons, 
                      FUN = function(x){get_length(x, trackspl_to_use[[1]])})
  } else {
    warning(paste("track length", paste0(tracksinoc, collapse = " "), "not equal to 1"))
    inters = rep(NA, length(grid_polygons$polygons))
  }
  return(inters)
}
#
#-------------Read files -------------------------------------------------------

lpoly <- readRDS("data/all_scenarios/larger_poly.Rds")
tracks <- readRDS("data/all_scenarios/lbl_trackdat.Rds")
olddat <- readRDS("data/onison/all_occasions/model_objs/chs_noneuclidean_extraprims2000.Rds")
sight <- readRDS("data/onison/all_occasions/sight_for_capthist_noopen_2000.Rds")
sightwn <-  readRDS("data/all_scenarios/sightings_data.Rds")
wpt_summaries <- readRDS("data/all_scenarios/wpt_summaries.Rds")
capthist <- readRDS("data/onison/all_occasions/capthistscr_noopen_2000.Rds")
prim <- readRDS("data/all_scenarios/all_occasions/primary.Rds")
traps <- readRDS("data/onison/all_occasions/trapscr_noopen_2000.Rds")
mesh <- readRDS("data/onison/all_occasions/meshscr_NSbuff_2000.Rds")
meshpoly <- readRDS("data/onison/all_occasions/meshpoly_2000.Rds")
distmat <- readRDS("data/onison/all_occasions/user_ne_dist_mat_2000.Rds")
meshdistmat <- readRDS("data/onison/all_occasions/mesh_dist_mat_2000.Rds")


#-------------Tracks dataframe -------------------------------------------------
#format one primary of data 
olddat$time()[11]+2004
surveys <-  unique(sight[which(sight$occ_key %in% which(prim$primary == 11)),"survey"])
oldoccs <- unique(sight[which(sight$occ_key %in% which(prim$primary == 11)),"occ_key"])
tracks <- tracks[tracks$ID %in% surveys,]
tracks <- tracks[order(tracks$ID, tracks$t),]

#sort by track ID, then t
tracksdf <- data.frame(occ = apply(as.array(1:nrow(tracks)), 1, 
                                   function(x){which(unique(tracks$ID) == tracks$ID[x])}),
                       x = tracks$x,
                       y = tracks$y,
                       effort = tracks$effort,
                       time = tracks$t) #time right now only determines order, not relative time
#note SIghting 1 in survey 710 (occasion 4) happens immediately, and the two trackpts aren't labelled on effort
tracksdf[which(tracksdf$occ == 4 & tracksdf$time < 101),"effort"] <- "OnEffort"
nocc <- length(unique(tracksdf$occ))
dx = 2000

#-------------Subset capthist -------------------------------------------------
capthist <- subset(capthist, occasions = oldoccs)
sight <- sight[sight$survey %in% surveys,]
sight$occ_key <- apply(as.array(sight$survey), 1, function(x){which(surveys == x)})
traps <- traps[which(rowSums(usage(traps)[,c(oldoccs)])>0),]
trapscr <- traps
usage(trapscr) <- usage(traps)[which(rowSums(usage(traps)[,c(oldoccs)])>0),c(oldoccs)]
distmatscr <- distmat[which(rowSums(usage(traps)[,c(oldoccs)])>0),]
#-------------Times of detections ----------------------------------------------
#need to identify time from survey start that each individual is detected in each survey for induse
sightwn <- sightwn[which(sightwn$SurveyNumber %in% surveys &
                           sightwn$CatalogID %in% unique(sight$ID) &
                           sightwn$OnEffort == "Yes"),]
sight <- left_join(sight, sightwn[,c("CatalogID", "SurveyNumber", "Sighting", "lat")], 
          by=c("ID"="CatalogID","survey" = "SurveyNumber", "lat" = "lat"))

sight$t <- NA

wpt_summaries <- wpt_summaries[which(as.numeric(wpt_summaries$ID) %in% surveys),]
for (i in 1:nrow(sight)){
  survi <- sight[i, "survey"]
  sightnumi <- sight[i, "Sighting"]
  wpt <- wpt_summaries[which(
    wpt_summaries$ID == survi & 
      wpt_summaries$SightingNum == sightnumi), 
    ]
  sight[i, "t"] <- wpt$t
}
#discard any second detections
pcomboIDsurv <- do.call(paste, sight[,c("ID", "survey")])
sight$freq <- apply(as.array(pcomboIDsurv), 1, function(x){length(which(pcomboIDsurv == x))})
repeaters <- unique(sight[sight$freq>1, c("ID", "survey")])
sight_single <- sight
for(row in 1:nrow(repeaters)){
  rID = repeaters[row, "ID"]
  rsurvey = repeaters[row, "survey"]
  sightrows = which(sight_single$ID == rID & sight_single$survey == rsurvey)
  print(length(sightrows))
  mults <- sight_single[sightrows,]
  sight_single <- sight_single[-sightrows,]
  sight_single <- rbind(sight_single, mults[order(mults$t),][1,])
}

sight_single$trap <- NA

#closest trap
for (tr in 1:nrow(sight_single)) {
  d <- sqrt((sight_single[tr, ]$x - traps[, 1])^2 + (sight_single[tr, ]$y - traps[, 2])^2)
  dmin <- min(d)
  if (dmin > 5000) {
    sight_single[tr, "trap"] <- NA
    next
  }
  sight_single[tr, "trap"] <- which.min(d)
}


#format this into a ch of times
cht <- apply(as.array(unique(sight_single$ID)), 1, 
      function(i){
        apply(as.array(surveys), 1, function(k){
            chik <- rep(NA, nrow(traps))
            s <- sight_single[which(sight_single$ID == i & 
                                      sight_single$survey ==k),]
            if(nrow(s) > 0){
              chik[s$trap] <- s$t
            }
            return(chik)
        })
      })

ch_t <- aperm(array(cht, dim = c(nrow(traps), nocc, length(unique(sight_single$ID)))),
               c(3,2,1))

trapcells <- create_grid_polygons(trapscr, spacing = 2000)
#use
#for each individual, create a use matrix
create_ind_use <- function(ch_t, trapcells, tracksdf){
  use <-mclapply(as.list(1:dim(ch_t)[1]), FUN = function(i){
    apply(as.array(1:dim(ch_t)[2]), 1, function(k){
      trackoccdf <- tracksdf[tracksdf$occ == k,]
      if(!all(is.na(ch_t[i,k,]))){
        #if i detected k, discard survey pts after detection in occasion k
        dettime <- ch_t[i,k,which(!is.na(ch_t[i,k,]))]
        trackoccdf <- trackoccdf[trackoccdf$time <= dettime,]
        #add up length of trackline in each grid cell
        useik <- lengths_in_grid(create_line_spatlines(trackoccdf, tracksdfcolname = "occ"), k, trapcells)
      } else {#if i wasn't detected in k, all traps used full amount, so skip subsetting
        useik <- usage(trapscr)[,k]
      }
      #(if there are covariates, perhaps could create matrix for them here, but for now ignore)
    })
    
  }, mc.cores = 3)
  return(use)
}

trapcells <- create_grid_polygons(trapscr, spacing = 2000)

# induse_ls <- create_ind_use(ch_t, trapcells, tracksdf)
# induse <- aperm(
#   array(unlist(induse_ls), 
#         dim =c(nrow(traps), nocc, dim(capthist)[1])), 
#   c(3, 1, 2))
# saveRDS(induse, "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS_results/induse.Rds")
induse <- readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS_results/induse.Rds")

#convert capthist back to 1s and 0s
ch_10 <- ch_t
ch_10[is.na(ch_10)] <- 0
ch_10[ch_10!=0] <- 1

#-------------------------fit normal secr model to data ------------------------
scrmesh <- make.mask(trapscr,
                     spacing = dx,
                     buffer = 0,
                     poly = meshpoly,
                     type = "polygon")
rawcap <- data.frame(session = 1,
                     ID = sight_single$ID,
                     occasion = sight_single$occ_key,
                     x = sight_single$x,
                     y = sight_single$y)

capthist <- make.capthist(captures = rawcap,
                          traps = trapscr,
                          fmt = "XY",
                          noccasions = length(unique(sight$survey)),
                          snapXY = TRUE,
                          tol = 5000)
args <- list(capthist = capthist, 
             mask = scrmesh, 
             detectfn = "HHN", 
            # method = "Nelder-Mead",
             details = list(userdist = distmatscr), 
             trace = FALSE)
models <- list(D ~ 1, 
               D ~ s(x, y, k = 5),
               D ~ s(x, y, k = 6),
               D ~ s(x, y, k = 7),
               D ~ s(x, y, k = 8),
               D ~ s(x, y, k = 9),
               D ~ s(x, y, k = 10),
               D ~ s(x, y, k = 11),
               D ~ s(x, y, k = 12),
               D ~ s(x, y, k = 13))
names <- c('null', #1
           'Dsmooth5', #2
           'Dsmooth6', #3
           'Dsmooth7', #4
           'Dsmooth8', #5
           'Dsmooth9', #6
           'Dsmooth10', #7
           'Dsmooth11', #8
           'Dsmooth12', #9
           'Dsmooth13') #10
fits <- list.secr.fit(model = models, constant = args, names = names)
AIC(fits[!apply(as.array(1:length(fits)), 1, function(x){any(is.na(attr(predictDsurface(fits[[x]], cl.D = T), "covariates")[,2]))})])
m0 <- fits[[2]]
lowerD <- attr(predictDsurface(m0, cl.D = T), "covariates")[,2]*100
upperD <- attr(predictDsurface(m0, cl.D = T), "covariates")[,3]*100
Dpar <- exp(m0$fit$par[m0$parindx$D]) #density per hectare
DdesignX <- m0$designD
denssurf <- exp(DdesignX %*% log(Dpar))*100
Ddiffco <- 100
spreadD <- function(mesh_dist_mat, lambda0, sigma, D){
  rowSums(apply(as.array(1:nrow(mesh_dist_mat)), 1, function(meshx){
    #this is the probability of detection at x 
    probdet = apply(as.array(1:nrow(mesh_dist_mat)), 1, function(x){
      lambda0 * exp(-(mesh_dist_mat[meshx,x]^2)/(2*sigma^2)) 
    })
    probdet = probdet/sum(probdet)
    D[meshx] * probdet #* area #keep density per square km, not per mesh
  }))
}
Dspreadsurf <- spreadD(meshdistmat, lambda0, sigma, denssurf)
Drowlimit <- (upperD-lowerD)<Ddiffco
ggplot() +
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 1) +
  geom_tile(data.frame(x = mesh$x, y = mesh$y, D = denssurf)[Drowlimit,], 
            mapping = aes(x = x, y =y, fill = D)) +
  geom_point(data = mesh, mapping = aes(x = x, y = y)) +
  geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "red", shape = "+") +
  geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
  scale_fill_viridis_c(name = "magma") +
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 0.3) +
  coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
           ylim = c(min(mesh$y), max(mesh$y))) +
  theme_bw()

Drowlimit <- apply(distmat, 2, min) < 1.5*sigma
ggplot() +
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 1) +
  geom_tile(data.frame(x = mesh$x, y = mesh$y, D = denssurf)[Drowlimit,], 
            mapping = aes(x = x, y =y, fill = D)) +
  geom_point(data = mesh, mapping = aes(x = x, y = y)) +
  geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "red", shape = "+") +
  geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
  scale_fill_viridis_c(name = "magma") +
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 0.3) +
  coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
           ylim = c(min(mesh$y), max(mesh$y))) +
  theme_bw()


sum(denssurf[(upperD-lowerD)<Ddiffco])*4
sum(denssurf)*4


#------------------------------Moving detector --------------------------------


move_fit <- function(capthist,
                     tracksdf, 
                     trapcells,
                     dist_trapmesh,
                     useall,
                     induse,
                     DdesignX, 
                     hazdenom, 
                     mesh, 
                     startpar0){
  
  mesh_mat <- as.matrix(mesh)
  par0 <- unlist(startpar0)
  #rescale for easy hessian
  scaling_factors <- 10^round(log10(abs(par0)))
  par <- par0/scaling_factors
  lambda0parindex <- which(names(par) == "lambda0")
  sigmaparindex <- which(names(par) == "sigma")
  Dparindex <- grep("D", names(par))
  
   #stationary detector likelihood
   stat_nll <- function(v_scaled){
     v <- v_scaled * scaling_factors 
    lambda0_ <- exp(v[lambda0parindex])
    sigma_ <- exp(v[sigmaparindex])
    D_mesh_ <- exp(DdesignX %*% v[Dparindex]) 
     out <- negloglikelihood_stationary_cpp(lambda0_, sigma_,
                                            hazdenom, D_mesh_, 
                                            capthist, useall,
                                            dist_trapmesh, mesh_mat)
     return(out)
   }
  #moving detector likelihood
  nll <- function(v_scaled){
    v <- v_scaled * scaling_factors 
    lambda0_ <- exp(v[lambda0parindex])
    sigma_ <- exp(v[sigmaparindex])
    D_mesh_ <- exp(DdesignX %*% v[Dparindex]) 
    out <- negloglikelihood_moving_cpp(lambda0_, sigma_,  
                                       hazdenom, D_mesh_,
                                       capthist, useall,
                                       induse, dist_trapmesh, mesh_mat)
    return(out)
  }
  
   start.time.sd <- Sys.time()
   fit_sd <- optim(par = par,
                   fn = stat_nll,
                   hessian = F, method = "Nelder-Mead")
   fit_sd$hessian <- numDeriv::hessian(stat_nll, x = fit_sd$par,
                                       method = "Richardson",
                                       method.args = list(eps = 1e-6, d = 1e-4, r = 4))
   fit.time.sd <- difftime(Sys.time(), start.time.sd, units = "secs")
  
  start.time.md <- Sys.time()
  fit_md <- optim(par = par,
                  fn = nll,
                  hessian = F, method = "Nelder-Mead")
  fit_md$hessian <- numDeriv::hessian(nll, x = fit_md$par,
                                      method = "Richardson",
                                      method.args = list(eps = 1e-6, d = 1e-4, r = 4))
  fit.time.md <- difftime(Sys.time(), start.time.md, units = "secs")
  

  
  assemble_CIs <- function(fit){
    fisher_info <- MASS::ginv(fit$hessian)
    prop_sigma <- sqrt(diag(fisher_info))
    prop_sigma <- diag(prop_sigma)
    upper <- fit$par+1.96*prop_sigma
    lower <- fit$par-1.96*prop_sigma
    interval <- data.frame(name = names(par),
                           value = fit$par * scaling_factors, 
                           sd = diag(prop_sigma),
                           upper = diag(upper) * scaling_factors, 
                           lower = diag(lower) * scaling_factors
    )
    
    return(interval)
  } 
  
  out <- list(statdet_est = assemble_CIs(fit_sd), 
              movdet_est = assemble_CIs(fit_md),
              statdet_time = fit.time.sd,
              movdet_time = fit.time.md)
  return(out)
}



setwd("~/Documents/UniStAndrews/MovingDetector")
Rcpp::sourceCpp("approx_movingdetectorlikelihood.cpp")
setwd("~/Documents/UniStAndrews/Dolphins/Charleston")

formulas <- list(D~s(x,y,k=5),
                 D~s(x,y,k=6),
                 D~s(x,y,k=7),
                 D~s(x,y,k=8),
                 D~s(x,y,k=9))

myfits <- lapply(as.list(1:5), function(f){
  formula = formulas[[f]]
  Dpar =  exp(fits[[f+1]]$fit$par[fits[[f+1]]$parindx$D]) 
  lambda0 = exp(fits[[f+1]]$fit$par[fits[[f+1]]$parindx$lambda0])
  sigma = exp(fits[[f+1]]$fit$par[fits[[f+1]]$parindx$sigma])
  
  startparf = list(D = log(Dpar), 
                   lambda0 = log(lambda0), 
                   sigma = log(sigma))
  
  datak <- data.frame(x = scrmesh$x, 
                      y = scrmesh$y)
  split <- interpret.gam(formula)
  sml =  mgcv::smoothCon(split$smooth.spec[[1]], data = datak, 
                         knots = NULL, absorb.cons = T, 
                         scale.penalty = T, 
                         null.space.penalty = F,
                         sparse.cons = 0, 
                         diagonal.penalty = F, 
                         apply.by = T, 
                         modCon = 0)
  DdesignX = cbind( rep(1, nrow(datak)), sml[[1]]$X)
  colnames(DdesignX) = c("(Intercept)", paste0("s(x,y).", 1:(ncol(DdesignX)-1)))
  
  m_move <- move_fit(capthist = ch_10,
                     tracksdf = tracksdf, 
                     trapcells = trapcells,
                     dist_trapmesh = distmatscr,
                     useall = usage(trapscr),
                     induse = induse,
                     DdesignX = DdesignX, 
                     hazdenom = 1, 
                     mesh = scrmesh, 
                     startpar0 = startparf)
  saveRDS(m_move, file = paste("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS_results/m_move", paste(formula)[3], ".Rds", sep = ""))
})


denssurf_stat <- exp(DdesignX %*% (m_move$statdet_est[m0$parindx$D,"value"]))*100
lambda0_stat <- exp(m_move$statdet_est[m0$parindx$lambda0, "value"])
sigma_stat <- exp(m_move$statdet_est[m0$parindx$sigma, "value"])
denssurf_move <- exp(DdesignX %*% (m_move$movdet_est[m0$parindx$D,"value"]))*100
lambda0_move <- exp(m_move$movdet_est[m0$parindx$lambda0, "value"])
sigma_move <- exp(m_move$movdet_est[m0$parindx$sigma, "value"])
diffdense <- denssurf_stat - denssurf_move


subarea <- which(apply(distmat, 2, min) < exp(m_move$movdet_est[10,2]))
Dspreadsurf_stat <- spreadD(meshdistmat, lambda0_stat, sigma_stat, denssurf_stat*as.numeric(1:length(denssurf_move) %in% subarea))
Dspreadsurf_move <- spreadD(meshdistmat, lambda0_move, sigma_move, denssurf_move*as.numeric(1:length(denssurf_move) %in% subarea))
Dspread_diff <- Dspreadsurf_stat- Dspreadsurf_move


lowcolor = "#000004FF" 
colorm1 = "#51127CFF" 
colorm2 = "#B63679FF"
colorm3 = "#FB8861FF"
highcolor = "#FCFDBFFF"

sharedmax = max(c(denssurf_stat[subarea,],
                denssurf_move[subarea,],
                denssurf[subarea,]
                ))
sharedmids = c(sharedmax*(1/5), sharedmax*(2/5), sharedmax*(3/5),  sharedmax*(4/5))

ggpubr::ggarrange(
# ggplot() +
#   geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
#           linewidth = .1, alpha = 1) +
#   geom_tile(data.frame(x = mesh$x, y = mesh$y, D = denssurf)[,], 
#             mapping = aes(x = x, y =y, fill = D)) +
#   #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
#   geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "red", shape = "+") +
#   geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
#   scale_fill_gradientn(colors = c(lowcolor, colorm1, colorm2, colorm3, highcolor),
#                        values = c(0, sharedmids, sharedmax)/sharedmax,
#                        limits = c(0, sharedmax),
#                        name = "AC Density", 
#                        breaks = c(0, sharedmids, sharedmax),
#                        labels = c(0, round(sharedmids,1),round(sharedmax, 2))) +  
#   geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
#           linewidth = .1, alpha = 0.3) +
#   coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
#            ylim = c(min(mesh$y), max(mesh$y))) +
#   ggtitle("secr") +
#   theme_bw(),
ggplot() +
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 1) +
  geom_tile(data.frame(x = mesh$x, y = mesh$y, D = denssurf_stat)[subarea,], 
            mapping = aes(x = x, y =y, fill = D)) +
  #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
  geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "red", shape = "+") +
  geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
  scale_fill_stepsn(
    colors = c(lowcolor, colorm1, colorm2, colorm3, highcolor),
    breaks = c(0,  sharedmids, sharedmax),
    values = (c(0,  sharedmids, sharedmax)/sharedmax),
    limits = c(0, sharedmax),
    name = "AC Density",
    labels = c(0, round(sharedmids,1),round(sharedmax, 2))
  ) + 
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 0.3) +
  coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
           ylim = c(min(mesh$y), max(mesh$y))) +
  ggtitle("Stationary Detector") +
  theme_bw(),
ggplot() +
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 1) +
  geom_tile(data.frame(x = mesh$x, y = mesh$y, D = denssurf_move)[subarea,], 
            mapping = aes(x = x, y =y, fill = D)) +
  #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
  geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "red", shape = "+") +
  geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
  scale_fill_stepsn(
    colors = c(lowcolor, colorm1, colorm2, colorm3, highcolor),
    breaks = c(0, sharedmids, sharedmax),
    values = (c(0,  sharedmids, sharedmax)/sharedmax),
    limits = c(0, sharedmax),
    name = "AC Density",
    labels = c(0, round(sharedmids,1),round(sharedmax, 2))
  ) + 
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 0.3) +
  coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
           ylim = c(min(mesh$y), max(mesh$y))) +
  ggtitle("Moving Detector") +
  theme_bw(),
ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

ggplot() +
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 1) +
  geom_tile(data.frame(x = mesh$x, y = mesh$y, D = diffdense)[subarea,], 
            mapping = aes(x = x, y =y, fill = D)) +
  #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
  geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "red", shape = "+") +
  geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
  scale_fill_viridis_c(
    #colors = c(lowcolor, colorm1, colorm2, colorm3, highcolor),
   # breaks = c(min(diffdense[subarea]),  sharedmids,100, sharedmax),
    #values = (c(min(diffdense[subarea]),  sharedmids,100, sharedmax)/sharedmax),
    #limits = c(min(diffdense[subarea]), max(diffdense[subarea])),
    name = "AC Density Difference"
  ) + 
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 0.3) +
  coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
           ylim = c(min(mesh$y), max(mesh$y))) +
  ggtitle("Stationary - Moving") +
  theme_bw()

sharedmax = max(c(Dspreadsurf_move,
                  Dspreadsurf_stat
))
sharedmids = c(sharedmax*(1/5), sharedmax*(2/5), sharedmax*(3/5),  sharedmax*(4/5))


ggpubr::ggarrange(
# ggplot() +
#   geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
#           linewidth = .1, alpha = 1) +
#   geom_tile(data.frame(x = mesh$x, y = mesh$y, D = Dspreadsurf)[,], 
#             mapping = aes(x = x, y =y, fill = D)) +
#   #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
#   geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "red", shape = "+") +
#   geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
#   scale_fill_gradientn(colors = c(lowcolor, colorm1, colorm2, colorm3, highcolor),
#                        values = c(0, sharedmids, sharedmax)/sharedmax,
#                        limits = c(0, sharedmax),
#                        name = "Animal D", 
#                        breaks = c(0, sharedmids, sharedmax),
#                        labels = c(0, round(sharedmids,1),round(sharedmax, 2))) +  
#   geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
#           linewidth = .1, alpha = 0.3) +
#   coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
#            ylim = c(min(mesh$y), max(mesh$y))) +
#   ggtitle("secr") +
#   theme_bw(),
ggplot() +
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 1) +
  geom_tile(data.frame(x = mesh$x, y = mesh$y, D = Dspreadsurf_stat)[,], 
            mapping = aes(x = x, y =y, fill = D)) +
  #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
  geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "red", shape = "+") +
  geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
  scale_fill_stepsn(
    colors = c(lowcolor, colorm1, colorm2, colorm3, highcolor),
    breaks = c(0,  sharedmids, sharedmax),
    values = (c(0,  sharedmids, sharedmax)/sharedmax),
    limits = c(0, sharedmax),
    name = "Animal D",
    labels = c(0, round(sharedmids,1),round(sharedmax, 2))
   ) +  
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 0.3) +
  coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
           ylim = c(min(mesh$y), max(mesh$y))) +
  ggtitle("Stationary Detector") +
  theme_bw() 
,
ggplot() +
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 1) +
  geom_tile(data.frame(x = mesh$x, y = mesh$y, D = Dspreadsurf_move)[,], 
            mapping = aes(x = x, y =y, fill = D)) +
  #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
  geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "red", shape = "+") +
  geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
  scale_fill_stepsn(
    colors = c(lowcolor, colorm1, colorm2, colorm3, highcolor),
    breaks = c(0,  sharedmids, sharedmax),
    values = (c(0,  sharedmids, sharedmax)/sharedmax),
    limits = c(0, sharedmax),
    name = "Animal D",
    labels = c(0, round(sharedmids,1),round(sharedmax, 2))
  ) + 
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 0.3) +
  coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
           ylim = c(min(mesh$y), max(mesh$y))) +
  ggtitle("Moving Detector") +
  theme_bw(),
ncol=2, nrow=1, common.legend = TRUE, legend="bottom")

ggplot() +
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 1) +
  geom_tile(data.frame(x = mesh$x, y = mesh$y, D = Dspread_diff), 
            mapping = aes(x = x, y =y, fill = D)) +
  #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
  geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "red", shape = "+") +
  geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
  scale_fill_viridis_c(
    #colors = c(lowcolor, colorm1, colorm2, colorm3, highcolor),
    # breaks = c(min(diffdense[subarea]),  sharedmids,100, sharedmax),
    #values = (c(min(diffdense[subarea]),  sharedmids,100, sharedmax)/sharedmax),
    #limits = c(min(diffdense[subarea]), max(diffdense[subarea])),
    name = "Animal Density Difference"
  ) + 
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 0.3) +
  coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
           ylim = c(min(mesh$y), max(mesh$y))) +
  ggtitle("Stationary - Moving") +
  theme_bw()

ggplot() +
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 1) +
  geom_tile(data.frame(x = mesh$x, y = mesh$y, D = Dspread_diff/Dspreadsurf_stat), 
            mapping = aes(x = x, y =y, fill = D)) +
  #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
  geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "red", shape = "+") +
  geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
  scale_fill_viridis_c(
    #colors = c(lowcolor, colorm1, colorm2, colorm3, highcolor),
    # breaks = c(min(diffdense[subarea]),  sharedmids,100, sharedmax),
    #values = (c(min(diffdense[subarea]),  sharedmids,100, sharedmax)/sharedmax),
    #limits = c(min(diffdense[subarea]), max(diffdense[subarea])),
    name = "Animal Density \n% Difference"
  ) + 
  geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
          linewidth = .1, alpha = 0.3) +
  coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
           ylim = c(min(mesh$y), max(mesh$y))) +
  ggtitle("Stationary - Moving") +
  theme_bw()


detpars <- data.frame(lambda0 = c(lambda0, lambda0_stat, lambda0_move),
                      sigma = c(sigma, sigma_stat, sigma_move),
                      name = c("secr", "stationary", "moving"))
detdat <- do.call(rbind, lapply(as.list(1:3), function(n){df=data.frame(x = seq(0,4*detpars$sigma[n], length.out = 20))
                                df$y = detpars$lambda0[n]*exp(-df$x^2/(2*detpars$sigma[n]^2))
                                df$name = detpars$name[n]
                                return(df)}))
ggplot() +
  geom_line(detdat[detdat$name %in% c("moving", "stationary"),],
            mapping = aes(x = x ,y = y, color = name, group = name),
            linewidth = 1.5) +
  geom_vline(data = detpars[detpars$name %in% c("moving","stationary"),],  
             mapping = aes(group = name, color = name, xintercept = sigma),
             linetype = "dashed",
             linewidth = 1.5) +
  scale_color_manual(name = "Model", 
                     values =  c("cornflowerblue", "goldenrod")) +
  theme_bw()





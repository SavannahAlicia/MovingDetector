library(ggplot2)
library(sf)
library(parallel)
library(sp)
library(openpopscr)
library(secr)
library(mgcv)
library(dplyr)
library(tidyr)

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
create_line_spatlines <- function(tracksdf, scenario = "onison", 
                                  tracksdfcolname = "ID",
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
all_poly_2000 <- readRDS("data/onison/all_occasions/all_poly_2000.Rds")


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
# AIC(fits[!apply(as.array(1:length(fits)), 1, function(x){any(is.na(attr(predictDsurface(fits[[x]], cl.D = T), "covariates")[,2]))})])
# m0 <- fits[[2]]
# lowerD <- attr(predictDsurface(m0, cl.D = T), "covariates")[,2]*100
# upperD <- attr(predictDsurface(m0, cl.D = T), "covariates")[,3]*100
# Dpar <- exp(m0$fit$par[m0$parindx$D]) #density per hectare
# DdesignX <- m0$designD
# denssurf <- exp(DdesignX %*% log(Dpar))*100
# Ddiffco <- 100
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
# Dspreadsurf <- spreadD(meshdistmat, lambda0, sigma, denssurf)
# Drowlimit <- (upperD-lowerD)<Ddiffco
# ggplot() +
#   geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
#           linewidth = .1, alpha = 1) +
#   geom_tile(data.frame(x = mesh$x, y = mesh$y, D = denssurf)[Drowlimit,], 
#             mapping = aes(x = x, y =y, fill = D)) +
#   geom_point(data = mesh, mapping = aes(x = x, y = y)) +
#   geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "red", shape = "+") +
#   geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
#   scale_fill_viridis_c(name = "magma") +
#   geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
#           linewidth = .1, alpha = 0.3) +
#   coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
#            ylim = c(min(mesh$y), max(mesh$y))) +
#   theme_bw()
# 
# Drowlimit <- apply(distmat, 2, min) < 1.5*sigma
# ggplot() +
#   geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
#           linewidth = .1, alpha = 1) +
#   geom_tile(data.frame(x = mesh$x, y = mesh$y, D = denssurf)[Drowlimit,], 
#             mapping = aes(x = x, y =y, fill = D)) +
#   geom_point(data = mesh, mapping = aes(x = x, y = y)) +
#   geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "red", shape = "+") +
#   geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
#   scale_fill_viridis_c(name = "magma") +
#   geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
#           linewidth = .1, alpha = 0.3) +
#   coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
#            ylim = c(min(mesh$y), max(mesh$y))) +
#   theme_bw()
# 
# 
# sum(denssurf[(upperD-lowerD)<Ddiffco])*4
# sum(denssurf)*4


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

#soap film smoother
cut <- list(matrix(c(1143961, bbox(all_poly_2000)[2,1], #anticlockwise, start bottom left
                     1156161,bbox(all_poly_2000)[2,1],
                     # bbox(all_poly_2000)[1,2], bbox(all_poly_2000)[2,1],
                     bbox(all_poly_2000)[1,2], 3669697,
                     bbox(all_poly_2000)[1,2], bbox(all_poly_2000)[2,2],
                     bbox(all_poly_2000)[1,1],  bbox(all_poly_2000)[2,2],
                     1150161, 3638997,
                     1143961, bbox(all_poly_2000)[2,1]), 
                   ncol = 2, 
                   byrow = TRUE))
box <- st_polygon(cut)
box <- st_geometry(box)
box <- st_set_crs(box, st_crs(all_poly_2000))
epoly <- st_intersection(st_as_sf(all_poly_2000), box)
boundline <- st_union(st_buffer(st_geometry(st_as_sf(epoly)),800))
crds <- coordinates(as_Spatial(st_cast(boundline, "LINESTRING")))[[1]][[1]]
bound <- list(list(x = crds[,1], y = crds[,2], f = rep(0, nrow(crds))))

knots_soap <- data.frame( x = c(#1162593, 
  1160659, 1158595, 
  1168527, 
  #1145050,
  #1176396, 
  #1147759,
  1167237, 1162593, 1176396),
  y = c(#3671596, 
    3659084, 3651086,
    3659986, 
    #3660760,
    #3666436,
    #3648506, 
    3647732, 3638444, 3653279))

knots_soap2 <- data.frame( x = c(
  1162593, 
  1160659, 1158595, 
  1168527, 
  1145050,
  1176396, 
  1147759,
  1167237, 1162593, 1176396),
  y = c(
    3671596, 
    3659084, 3651086,
    3659986, 
    3660760,
    3666436,
    3648506, 
    3647732, 3638444, 3653279))

formulas <- list(D~s(x,y,k=5),
                 D~s(x,y,k=6),
                 D~s(x,y,k=7),
                 D~s(x,y,k=8),
                 D~s(x,y,k=9),
                 D~s(x, y, bs = "so", xt = list(bnd = bound)),
                 D~s(x, y, bs = "so", xt = list(bnd = bound)))

get_X_mat <- function(f, knots = NULL){
  formula <- formulas[[f]]
  split <- interpret.gam(formula)
  sml =  mgcv::smoothCon(split$smooth.spec[[1]], data = scrmesh, 
                         knots = knots, absorb.cons = T, 
                         scale.penalty = T, 
                         null.space.penalty = F,
                         sparse.cons = 0, 
                         diagonal.penalty = F, 
                         apply.by = T, 
                         modCon = 0)
  DdesignX = cbind( rep(1, nrow(scrmesh)), sml[[1]]$X)
  colnames(DdesignX) = c("(Intercept)", paste0("s(x,y).", 1:(ncol(DdesignX)-1)))
  return(DdesignX)
  }

fit_smooth <- function(f, startother = NULL, addtl_name = ""){
  formula = formulas[[f]]
  if(is.null(startother)){
    Dpar =  exp(fits[[f+1]]$fit$par[fits[[f+1]]$parindx$D]) 
    lambda0 = exp(fits[[f+1]]$fit$par[fits[[f+1]]$parindx$lambda0])
    sigma = exp(fits[[f+1]]$fit$par[fits[[f+1]]$parindx$sigma])
  } else {
    Dpar =  exp(fits[[startother]]$fit$par[fits[[startother]]$parindx$D]) 
    lambda0 = exp(fits[[startother]]$fit$par[fits[[startother]]$parindx$lambda0])
    sigma = exp(fits[[startother]]$fit$par[fits[[startother]]$parindx$sigma])
  }
  
  startparf = list(D = log(Dpar), 
                   lambda0 = log(lambda0), 
                   sigma = log(sigma))
  DdesignX <- Xmats[[f]]
  
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
  saveRDS(m_move, file = paste("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS_results/m_move", paste(formula)[3], addtl_name, ".Rds", sep = ""))
  return(m_move)
}

Xmats <- lapply(as.list(1:5), get_X_mat)
Xmats[[6]] <- get_X_mat(f = 6, knots = knots_soap)
Xmats[[7]] <- get_X_mat(f = 6, knots = knots_soap2)

#myfits <- lapply(as.list(1:5), fit_smooth)
#myfits[[6]] <- fit_smooth(6, startother = 4)
#myfits[[7]] <- fit_smooth(7, startother = 8, addtl_name = "moreknots")
myfits <- list(readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS_results/m_moves(x, y, k = 5).Rds"),
               readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS_results/m_moves(x, y, k = 6).Rds"),
               readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS_results/m_moves(x, y, k = 7).Rds"),
               readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS_results/m_moves(x, y, k = 8).Rds"),
               readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS_results/m_moves(x, y, k = 9).Rds"),
               readRDS('~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS_results/m_moves(x, y, bs = "so", xt = list(bnd = bound)).Rds'),
               readRDS('~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS_results/m_moves(x, y, bs = "so", xt = list(bnd = bound))moreknots.Rds'))

AICs <- apply(as.array(1:length(myfits)), 1, function(i){
  AIC_fn <- function(n,L){2*n + 2*L}
  parloc <- list(lambda0 = which(myfits[[i]]$movdet_est$name == "lambda0"),
                 sigma = which(myfits[[i]]$movdet_est$name == "sigma"),
                 D = which(grepl("D", myfits[[i]]$movdet_est$name)))

  movL <- negloglikelihood_moving_cpp(lambda0 = exp(myfits[[i]]$movdet_est$value[parloc$lambda0]), 
                                      sigma = exp(myfits[[i]]$movdet_est$value[parloc$sigma]),  
                              timeincr = 1, 
                              D_mesh = exp(Xmats[[i]] %*% myfits[[i]]$movdet_est$value[parloc$D]),
                              capthist= ch_10, 
                              usage = usage(trapscr),
                              indusage = induse, 
                              distmat =  distmatscr,
                              mesh = as.matrix(scrmesh))
  statL <- negloglikelihood_stationary_cpp(lambda0 = exp(myfits[[i]]$movdet_est$value[parloc$lambda0]), 
                                           sigma = exp(myfits[[i]]$movdet_est$value[parloc$sigma]),  
                                           timeincr = 1, 
                                           D_mesh = exp(Xmats[[i]] %*% myfits[[i]]$movdet_est$value[parloc$D]),
                                           capthist= ch_10, 
                                           usage = usage(trapscr),
                                           distmat =  distmatscr,
                                           mesh = as.matrix(scrmesh))
  statAIC <- AIC_fn(n = length(unlist(parloc)),
                    L = statL)
  movAIC <- AIC_fn(n = length(unlist(parloc)),
                   L = movL)
  return(
    data.frame(stat = statAIC, 
                    mov = movAIC, 
                    mod = paste(formulas[[i]])[c(3)])
         )
})
myAICs <- do.call(rbind, AICs)


create_plots <- function(m_move, DdesignX, m0, label, subarea_sigmamult = NULL){
  
  denssurf_stat <- exp(DdesignX %*% (m_move$statdet_est[m0$parindx$D,"value"]))#*100 returns in km^2  now
  lambda0_stat <- exp(m_move$statdet_est[m0$parindx$lambda0, c("value", "lower", "upper")])
  sigma_stat <- exp(m_move$statdet_est[m0$parindx$sigma, c("value", "lower", "upper")])
  denssurf_move <- exp(DdesignX %*% (m_move$movdet_est[m0$parindx$D,"value"]))#*100
  lambda0_move <- exp(m_move$movdet_est[m0$parindx$lambda0, c("value", "lower", "upper")])
  sigma_move <- exp(m_move$movdet_est[m0$parindx$sigma, c("value", "lower", "upper")])
  diffdense <- denssurf_stat - denssurf_move
  
  if(is.null(subarea_sigmamult)){
    subarea <- 1:ncol(distmat)
  } else {
    subarea <- which(apply(distmat, 2, min) < exp(m_move$movdet_est[m0$parindx$sigma,2])) * subarea_sigmamult
  }
  #
  Dspreadsurf_stat <- spreadD(meshdistmat, lambda0_stat$value, sigma_stat$value, denssurf_stat*as.numeric(1:length(denssurf_move) %in% subarea))
  Dspreadsurf_move <- spreadD(meshdistmat, lambda0_move$value, sigma_move$value, denssurf_move*as.numeric(1:length(denssurf_move) %in% subarea))
  Dspread_diff <- Dspreadsurf_stat- Dspreadsurf_move
  
  
  lowcolor = "#000004FF" 
  colorm1 = "#51127CFF" 
  colorm2 = "#B63679FF"
  colorm3 = "#FB8861FF"
  highcolor = "#FCFDBFFF"
  
  sharedmax = max(c(denssurf_stat[subarea,],
                    denssurf_move[subarea,]
  ))
  sharedmids = c(sharedmax*(1/5), sharedmax*(2/5), sharedmax*(3/5),  sharedmax*(4/5))
  
  
  ACDstatplot <- ggplot() +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 1) +
    geom_tile(data.frame(x = mesh$x, y = mesh$y, D = denssurf_stat)[subarea,], 
              mapping = aes(x = x, y =y, fill = D)) +
    #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
    geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "white", shape = "+") +
    #geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
    scale_fill_viridis_c(limits = c(0, sharedmax),
                         name = "AC Density",
                         option = "magma") +
    # scale_fill_stepsn(
    #   colors = c(lowcolor, colorm1, colorm2, colorm3, highcolor),
    #   breaks = c(0,  sharedmids, sharedmax),
    #   values = (c(0,  sharedmids, sharedmax)/sharedmax),
    #   limits = c(0, sharedmax),
    #   name = "AC Density",
    #   labels = c(0, round(sharedmids,1),round(sharedmax, 2))
    # ) + 
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 0.3) +
    coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
             ylim = c(min(mesh$y), max(mesh$y))) +
    ggtitle("Stationary Detector") +
    theme_bw()+
    theme(legend.position = "none")
  ACDmovplot <- ggplot() +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 1) +
    geom_tile(data.frame(x = mesh$x, y = mesh$y, D = denssurf_move)[subarea,], 
              mapping = aes(x = x, y =y, fill = D)) +
    #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
    geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "white", shape = "+") +
    #geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") + 
    scale_fill_viridis_c(limits = c(0, sharedmax),
                         name = "AC Density",
                         option = "magma") +
    # scale_fill_stepsn(
    #   colors = c(lowcolor, colorm1, colorm2, colorm3, highcolor),
    #   breaks = c(0, sharedmids, sharedmax),
    #   values = (c(0,  sharedmids, sharedmax)/sharedmax),
    #   limits = c(0, sharedmax),
    #   name = "AC Density",
    #   labels = c(0, round(sharedmids,1),round(sharedmax, 2))
    # ) + 
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 0.3) +
    coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
             ylim = c(min(mesh$y), max(mesh$y))) +
    ggtitle("Moving Detector") +
    theme_bw()+
    theme(legend.position = "none")
  
  ACDlegendplot <- ggplot(data = data.frame(x = mesh$x, y = mesh$y, D = denssurf_move)[subarea,], 
                          mapping = aes(x = x, y =y, fill = D)) +
    scale_fill_viridis_c(limits = c(0, sharedmax),
                         name = "AC Density",
                         option = "magma") +
    geom_point(alpha = 0, shape = 0) +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "inside",
          legend.direction = "horizontal",
          legend.position.inside = c(0.5, 0.5), # move the legend to the center
          legend.key = element_rect(fill='NA'),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank())
  
  ACDdiffplot<- ggplot() +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 1) +
    geom_tile(data.frame(x = mesh$x, y = mesh$y, D = diffdense)[subarea,], 
              mapping = aes(x = x, y =y, fill = D)) +
    #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
    geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "white", shape = "+") +
    # geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
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
    theme_bw() +
    theme(legend.position = "none")
  
  ACDdifflegendplot <- ggplot(data.frame(x = mesh$x, y = mesh$y, D = diffdense)[subarea,], 
                              mapping = aes(x = x, y =y, fill = D)) +
    scale_fill_viridis_c(
      name = "AC Density\nDifference"
    ) + 
    geom_point(alpha = 0, shape = 0) +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "inside",
          legend.direction = "horizontal",
          legend.position.inside = c(0.5, 0.5), # move the legend to the center
          legend.key = element_rect(fill='NA'),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank())
  
  
  ggsave(file = paste("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS_results/", label, "ACDplot.png", sep = ""),
         plot = grid.arrange(
           grobs = list(ACDstatplot, ACDmovplot, ACDdiffplot,
                        ACDlegendplot, ACDdifflegendplot),
           widths = c((1), (1), (1)),
           heights = c(1,.2),
           layout_matrix = rbind(c(1,  2, 3),
                                 c(4, 4, 5))),
         width = 230,
         height = 230*.5,
         units = c("mm"),
         dpi = 300)
  
  sharedmax = max(c(Dspreadsurf_move,
                    Dspreadsurf_stat
  ))
  sharedmids = c(sharedmax*(1/5), sharedmax*(2/5), sharedmax*(3/5),  sharedmax*(4/5))
  
  
  animDstat <- ggplot() +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 1) +
    geom_tile(data.frame(x = mesh$x, y = mesh$y, D = Dspreadsurf_stat)[,], 
              mapping = aes(x = x, y =y, fill = D)) +
    #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
    geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "white", shape = "+") +
    # geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
    scale_fill_viridis_c(option = "magma",
                         name = "Animal Density",
                         limits = c(0, sharedmax)
    ) +  
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 0.3) +
    coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
             ylim = c(min(mesh$y), max(mesh$y))) +
    ggtitle("Stationary Detector") +
    theme_bw() +
    theme(legend.position = "none")
  
  animDmov <- ggplot() +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 1) +
    geom_tile(data.frame(x = mesh$x, y = mesh$y, D = Dspreadsurf_move)[,], 
              mapping = aes(x = x, y =y, fill = D)) +
    #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
    geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "white", shape = "+") +
    #geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
    scale_fill_viridis_c(limits = c(0, sharedmax),
                         name = "Animal Density",
                         option = "magma") +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 0.3) +
    coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
             ylim = c(min(mesh$y), max(mesh$y))) +
    ggtitle("Moving Detector") +
    theme_bw() +
    theme(legend.position = "none")
  
  animDlegend <- ggplot(data.frame(x = mesh$x, y = mesh$y, D = Dspreadsurf_move)[,], 
                        mapping = aes(x = x, y =y, fill = D)) +
    scale_fill_viridis_c(limits = c(0, sharedmax),
                         name = "Animal Density",
                         option = "magma") +
    geom_point(alpha = 0, shape = 0) +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "inside",
          legend.direction = "horizontal",
          legend.position.inside = c(0.5, 0.5), # move the legend to the center
          legend.key = element_rect(fill='NA'),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank())
  
  animDdiff <- ggplot() +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 1) +
    geom_tile(data.frame(x = mesh$x, y = mesh$y, D = Dspread_diff), 
              mapping = aes(x = x, y =y, fill = D)) +
    #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
    geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "white", shape = "+") +
    #geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
    scale_fill_viridis_c(
      name = "Animal Density Difference"
    ) + 
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 0.3) +
    coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
             ylim = c(min(mesh$y), max(mesh$y))) +
    ggtitle("Stationary - Moving") +
    theme_bw() + 
    theme(legend.position = "none")
  
  animDpercdiff <- ggplot() +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 1) +
    geom_tile(data.frame(x = mesh$x, y = mesh$y, D = Dspread_diff/Dspreadsurf_stat), 
              mapping = aes(x = x, y =y, fill = D)) +
    #geom_point(data = mesh, mapping = aes(x = x, y = y)) +
    geom_point(data = trapscr, mapping = aes(x = x, y =y), color = "white", shape = "+") +
    #geom_point(data = sight_single, mapping = aes(x = x, y =y), color = "green") +
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
  
  animDdifflegend <- 
    ggplot(data.frame(x = mesh$x, y = mesh$y, D = Dspread_diff), 
           mapping = aes(x = x, y =y, fill = D)) +
    scale_fill_viridis_c(limits = c(0, sharedmax),
                         name = "Animal Density\nDifference") +
    geom_point(alpha = 0, shape = 0) +
    scale_y_continuous(limits = c(0,1)) +
    scale_x_continuous(limits = c(0,1)) +
    theme_classic() +
    theme(axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          legend.position = "inside",
          legend.direction = "horizontal",
          legend.position.inside = c(0.5, 0.5), # move the legend to the center
          legend.key = element_rect(fill='NA'),
          panel.grid = element_blank(),
          axis.line = element_blank(),
          panel.border = element_blank())
  
  ggsave(file = paste("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS_results/", label, "animDplot.png", sep = ""),
         plot = grid.arrange(
           grobs = list(animDstat, animDmov, animDdiff,
                        animDlegend, animDdifflegend),
           widths = c((1), (1), (1)),
           heights = c(1,.2),
           layout_matrix = rbind(c(1,  2, 3),
                                 c(4, 4, 5))),
         width = 230,
         height = 230*.5,
         units = c("mm"),
         dpi = 300)
  
  
  
  detpars <- data.frame(lambda0 = unlist(c(lambda0_stat, lambda0_move)),
                        sigma = unlist(c(sigma_stat, sigma_move)),
                        name = names(c(lambda0_stat, lambda0_move)),
                        model = c(rep("stationary", 3), rep("moving", 3)))
  detdat <- do.call(rbind, lapply(as.list(1:6), function(n){
    df=
      data.frame(x = seq(0,4*max(detpars$sigma), length.out = 20))
    
    df$y = detpars$lambda0[n]*exp(-df$x^2/(2*detpars$sigma[n]^2))
    
    df$name = detpars$name[n]
    df$model = detpars$model[n]
    return(df)}))
  detdatwide <-  pivot_wider(detdat, values_from = y, names_from = name)
  
  detfctplot <- ggplot() +
    geom_line(detdat[detdat$name %in% c("value"),],
              mapping = aes(x = x, y = y, color = model, group = model),
              linewidth = 1.5) +
    geom_ribbon(detdatwide,
                mapping = aes(x = x, ymin = lower, ymax = upper, 
                              fill = model, group = model),
                color = "transparent",
                alpha = .5) +
    geom_point(data = data.frame(x = detpars[detpars$name %in% c("value"),"sigma"], y = 0, 
                                 model = detpars[detpars$name %in% c("value"),"model"]),  
               mapping = aes(group = model, color = model, x = x, y = y),
               size = 3, shape = 4, stroke = 1.5) +
    ylab("Detection rate") +
    xlab("Distance (m)")+
    scale_color_manual(name = "", 
                       values =  c("cornflowerblue", "goldenrod"),
                       labels = c("Moving", "Stationary")) +
    scale_fill_manual(name = "", 
                      values =  c("cornflowerblue", "goldenrod"),
                      labels = c("Moving", "Stationary")) +
    theme_bw()
  
  ggsave(file = paste("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS_results/", label, "detplot.png", sep = ""),
         plot = detfctplot,
         width = 169,
         height = 169*.7,
         units = c("mm"),
         dpi = 300)
  
  tracksdf$xnext <- c(tracksdf$x[2:(nrow(tracksdf))], NA)
  tracksdf$ynext <- c(tracksdf$y[2:(nrow(tracksdf))], NA)
  for(occi in unique(tracksdf$occ)){
    tracksdf[tracksdf$occ == occi,]$xnext[nrow(tracksdf[tracksdf$occ == occi,])] <- 
      tracksdf[tracksdf$occ == occi,]$ynext[nrow(tracksdf[tracksdf$occ == occi,])] <-
      NA
  }
  tracksdf$dist_betw <- c(apply(as.array(1:(nrow(tracksdf)-1)), 1, function(x){
    dist(rbind(c(tracksdf[x,c("x","y")]), c(tracksdf[(x),c("xnext", "ynext")])))
  }), NA)
  
  tracksdf$cumdist <- 0
  for(occi in unique(tracksdf$occ)){
    tracksdf[tracksdf$occ == occi,]$cumdist <- cumsum(tracksdf[tracksdf$occ == occi,]$dist_betw)
  }
  
  #
  
  tracksdf$arrowbool <- F
  cumdist_n <- plyr::round_any(tracksdf$cumdist, 2000)
  tracksdf$arrowbool[which(cumdist_n[-1] != cumdist_n[-nrow(tracksdf)])] <- T
  
  #visualize track direction
  tracklines_ls <- create_line_spatlines(tracksdf, tracksdfcolname = "occ")
  
  plotsurv <- function(occi){
    survplot <- ggplot() +
      geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
              linewidth = .1, alpha = 0.3) +
      geom_sf(st_as_sf(tracklines_ls[[occi]]), mapping = aes()
      ) +
      geom_segment(data = tracksdf[which(tracksdf$occ == occi & 
                                           tracksdf$effort == "OnEffort" &
                                           tracksdf$arrowbool),],
                   mapping = aes(x = x, y = y,
                                 xend = xnext, yend = ynext
                   ),
                   color = "red",
                   arrow = arrow(angle = 45, 
                                 ends = "last", 
                                 type = "open", 
                                 length = unit(0.1, "cm"))) +
      #scale_color_viridis_c(name = "Hours") +
      coord_sf(xlim = c(min(tracksdf[tracksdf$occ == occi,]$x), max(tracksdf[tracksdf$occ == occi,]$x)), 
               ylim = c(min(tracksdf[tracksdf$occ == occi,]$y), max(tracksdf[tracksdf$occ == occi,]$y))) +
      theme_bw() +
      annotate("text", x = -Inf, y = Inf,
               label = paste("Survey", occi), vjust = 1.2, hjust = -0.1) +
      theme(legend.position = "none",
            axis.text = element_blank(),
            panel.grid = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank())
    return(survplot)
  }
  surv1 <- plotsurv(1)
  surv2 <- plotsurv(2)
  surv3 <- plotsurv(3)
  surv4 <- plotsurv(4)
  surv5 <- plotsurv(5)
  surv6 <- plotsurv(6)
  
  ggsave(file = paste("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS_results/", label, "surveydirection.png", sep = ""),
         plot = 
           grid.arrange(
             grobs = list(surv1, surv2, surv3,
                          surv4, surv5, surv6),
             widths = c((1), (1)),
             heights = c(1, 1, 1),
             layout_matrix = rbind(c(1,  4),
                                   c(2, 5),
                                   c(3, 6)))
         ,
         width = 150,
         height = 250,
         units = c("mm"),
         dpi = 300)
  return(list(statN = sum(Dspreadsurf_stat),
              movN = sum(Dspreadsurf_move)))
}
create_plots(m_move = myfits[[6]],
             DdesignX = Xmats[[6]],
             m0 = fits[[4]],
             label = paste(formulas[[6]])[3])
create_plots(m_move = myfits[[3]],
             DdesignX = Xmats[[3]],
             m0 = fits[[4]],
             label = paste(formulas[[3]])[3])
create_plots(m_move = myfits[[3]],
             DdesignX = Xmats[[3]],
             m0 = fits[[4]],
             label = paste(paste(formulas[[3]])[3], "subarea1", sep = ""),
             subarea_sigmamult = 1)



colSums(usage(trapscr))

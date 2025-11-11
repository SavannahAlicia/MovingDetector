#compare moving and stationary detector fits
library(secr)
library(lubridate)
library(parallel)
library(ggplot2)
library(sp)
library(sf)
library(gridExtra)
library(dplyr)
setwd("~/Documents/UniStAndrews/MovingDetector")
Rcpp::sourceCpp("approx_movingdetectorlikelihood.cpp")

#-------------------------------functions---------------------------------------
#' Create polygons for each box centered at grid pt
#' 
#' @param grid dataframe with x and y columns specifying grid pts
#' @param spacing the spacing between grid pts, if not specified will calculate
#' @return SpatialPolygons object in list form 
#' @export
create_grid_polygons <- function(grid, spacing = NULL){
  #so we only have to create gridbboxes and lines list once
  
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
create_line_spatlines <- function(tracksdf, scenario = "everything",
                                  projto = sp::CRS(sf::st_crs(26916)$input)){
  #create a line for between each pt of track, put into named list by ID
  #check that it's not the last point in a track
  #(we don't want a line between tracks)
  check_pt <- function(row, tracksdf. = tracksdf){
    if (tracksdf.$occ[row] == tracksdf.$occ[row+1]){
      newline <- Line(tracksdf.[c(row, row+1),c("x", "y")])
    } else {
      newline <- NULL
    }
    return(newline)
  }
  linelist <- lapply(X = as.list(1:(nrow(tracksdf)-1)), FUN = check_pt)
  names(linelist) <- tracksdf$occ[1:(length(tracksdf$occ)-1)]
  linelist <- linelist[which(!lapply(linelist, is.null) == TRUE)]
  #Testing: check there are lines for each track segment, but not between tracksdf
  #dim(tracksdf)[1] - length(unique(tracksdf$occ)) == length(linelist)
  #for each track ID, create a Lines object by scenario
  #append that Lines object to a list of them
  big_Lines_list <- list()
  nocc <- length(unique(tracksdf$occ))
  for (trackindex in 1:nocc){
    #print(trackindex)
    trackIDi <- unique(tracksdf$occ)[trackindex]
    thetrack <- tracksdf[which(tracksdf$occ == trackIDi),]
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
    warning(paste("number of tracks per occasion", paste0(tracksinoc, collapse = " "), "not equal to 1"))
    inters = rep(NA, length(grid_polygons$polygons))
  }
  return(inters)
}

#' Simulate moving detector capture history
#' 
#' @param pop defaults to null, population locations from sim.popn
#' @param traps dataframe of x y trap coordinates for trap grid
#' @param tracksdf dataframe with track locations (xy), occ, time, and trap 
#'                  number for corresponding trap grid row
#' @param lambda0 
#' @param sigma 
#' @param D_mesh density at each mesh pt (for simulating population if not specified)
#' 
#' @return capture history array of dimension individual by trap by occ that 
#' contains NA if not captured or time of first capture (seconds since the first
#' time logged for that occasion) if captured
#' @export
sim_capthist <- function(pop = NULL, 
                         dist_dat_pop = NULL,
                         traps, 
                         tracksdf,
                         lambda0, 
                         sigma, 
                         D_mesh,
                         hazdenom, #for hazard rate
                         report_probseenxk = FALSE){
  if(is.null(pop)){
    pop <- sim.popn(D = D_mesh, core = mesh, model2D = "IHP", #D in sim.popn is inds per hectare
                    Ndist = "poisson", buffertype = "rect")
    rownames(pop) <- NULL
  }
  if(is.null(dist_dat_pop)){
  #this is now distances between pop hrcs and track locations
  dist_dat_pop <- calc_dist_matC(as.matrix(tracksdf[,c("x","y")]), 
                                 as.matrix(pop)) #ind by track location
  }
  capthist <- lapply(as.list(1:nrow(pop)), #for each individual
                     FUN = function(i){
                       lapply(as.list(unique(tracksdf$occ)),
                              FUN = function(occk){
                                #prob individual is seen this occasion (survey)
                                trackoccdf <- tracksdf[tracksdf$occ == occk,]
                                begintimek <-min(trackoccdf$time)
                                #NOTE if you want a distance denominated hazard, use length diffs for hazdenom
                                #if you want rate per unit time, use time increment in seconds
                                increments <- #c(0, 
                                  # difftime(trackoccdf$time[-1], 
                                  #          trackoccdf$time[-nrow(trackoccdf)], 
                                  #          unit = "secs"))
                                  c(0, apply(as.array(1:(nrow(trackoccdf)-1)), 1, function(x){
                                    sqrt((tracksdf[x, c("x")]-tracksdf[x+1, c("x")])^2 + (tracksdf[x, c("y")]-tracksdf[x+1, c("y")])^2)
                                    }))
                                hazs <- apply(as.array(dist_dat_pop[i,tracksdf$occ == occk]),
                                              1, FUN = function(d){
                                                hazdist_cpp(lambda0, sigma, d, hazdenom)
                                              })
                                integ <- sum(hazs * (increments/hazdenom)) #hazard is per time incr, need how many increments (or if dist, per dist incr)
                                survk <- exp(-integ)
                                probseenxk <- 1 - survk 
                                if (report_probseenxk) { 
                                  return(probseenxk)
                                } else {
                                  #seenxk_bool <- rbinom(1,1, probseenxk)
                                  #if (seenxk_bool){
                                  survive_until_t <- exp(-1 * (cumsum(hazs * (increments/hazdenom))))
                                  survive_t_inc <- exp(-1 * hazs * (increments/hazdenom))
                                  seenfirstat_t <- survive_until_t[-length(survive_until_t)] * (1 - survive_t_inc[-1]) #really seen first in time increment that ends in t and begins at t-1
                                  seenfirstatt <- c(0, seenfirstat_t)
                                  #need to record time of detection at traps as defined in grid
                                  det <- sample(x = c(1:(nrow(trackoccdf)+1)), size = 1, replace = T,  
                                                prob = c(seenfirstatt, survk))
                                  capik <- rep((NA), nrow(traps))
                                  if(det != (length(seenfirstatt) + 1)){ #if didn't survive detection
                                    trapdet <- trackoccdf[det,"trapno"]
                                    timedet <- difftime(trackoccdf[det, "time"], begintimek, units = "secs")
                                    capik[trapdet] <- timedet
                                  }
                                  return(capik) 
                                }
                                
                              })
                     }) 
  if(report_probseenxk){
    capthist_array <- t(
      array(unlist(capthist), 
            dim = c(nocc,nrow(pop))) 
    ) #ind by occasion
  } else {
    capthist_array <- aperm(
      array(unlist(capthist), 
            dim =c(nrow(traps), nocc,  nrow(pop))), 
      c(3, 2, 1)) #ind x occ x traps
    
  }
  return(capthist_array)
}

#' Create individual use matrix
#' 
#' @param ch capture history matrix i x k x j of times of detection
#' @param trapcells polygons of grid cells of each trap
#' @param tracksdf dataframe of track locations, occ, and times
#' 
#' @return list of use matrices for each individual
create_ind_use <- function(ch, trapcells, tracksdf){
  use <- mclapply(as.list(1:dim(ch)[1]), FUN = function(i){
    apply(as.array(1:dim(ch)[2]), 1, function(k){
      trackoccdf <- tracksdf[tracksdf$occ == k,]
      if(!all(is.na(ch[i,k,]))){
        #if i detected k, discard survey pts after detection in occasion k
        dettime <- min(trackoccdf$time) + ch[i,k,which(!is.na(ch[i,k,]))]
        trackoccdf <- trackoccdf[trackoccdf$time <= dettime,]
        #add up length of trackline in each grid cell
        useik <- lengths_in_grid(create_line_spatlines(trackoccdf), k, trap_cells)
      } else {#if i wasn't detected in k, all traps used full amount, so skip subsetting
        useik <- useall[,k]
      }
      #(if there are covariates, perhaps could create matrix for them here, but for now ignore)
    })
  }, mc.cores = 3)
  return(use)
}

#--------------------------------------true parameters -------------------------
lambda0 = .001
sigma = 300
beta1 <- -(1/40000)
beta2 <- -2500
flatD = 0.4

nsims = 50
set.seed(1994)




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


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

ggsave(file = paste("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/soapfilmsetup.png", sep = ""),
       plot = grid.arrange(
         grobs = list(  ggplot() +
                          geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
                                  linewidth = .1, alpha = 1) + 
                          geom_sf(boundline, mapping = aes()) +
                          geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
                                  linewidth = .1, alpha = .5) +
                          geom_point(data = knots_soap, mapping = aes(x = x, y = y)) +
                          # geom_point(data.frame(x = bound[[1]]$x, y = bound[[1]]$y), mapping = aes(x = x,y = y)) +
                          coord_sf(xlim = c(min(bound[[1]]$x), max(bound[[1]]$x)), 
                                   ylim = c(min(bound[[1]]$y), max(bound[[1]]$y))) +
                          theme_bw() +
                          theme(axis.title = element_blank()),
                        ggplot() +
                          geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
                                  linewidth = .1, alpha = 1) + 
                          geom_sf(boundline, mapping = aes()) +
                          geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
                                  linewidth = .1, alpha = .5) +
                          geom_point(data = knots_soap2, mapping = aes(x = x, y = y)) +
                          # geom_point(data.frame(x = bound[[1]]$x, y = bound[[1]]$y), mapping = aes(x = x,y = y)) +
                          coord_sf(xlim = c(min(bound[[1]]$x), max(bound[[1]]$x)), 
                                   ylim = c(min(bound[[1]]$y), max(bound[[1]]$y))) +
                          theme_bw() +
                          theme(axis.title = element_blank(),
                                axis.text.y = element_text(color = "transparent"))),
         widths = c((1), (1)),
         heights = c(1),
         layout_matrix = rbind(c(1,  2)))
       ,
       width = 169,
       height = 169/2,
       units = c("mm"),
       dpi = 300)


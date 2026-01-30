
create_plots <- function(m_move, DdesignX, label, subareacutoff = NULL){
  Dind <- grep("D", m_move$statdet_est$name)
  lind <- grep("lambda", m_move$statdet_est$name)
  sind <- grep("sigma", m_move$statdet_est$name)
  denssurf_stat <- exp(DdesignX %*% (m_move$statdet_est[Dind,"value"]))*1000000 #returns in m^2  now
  abund_stat <- sum(denssurf_stat) * 4
  lambda0_stat <- exp(m_move$statdet_est[lind, c("value", "lower", "upper")])
  sigma_stat <- exp(m_move$statdet_est[sind, c("value", "lower", "upper")])
  denssurf_move <- exp(DdesignX %*% (m_move$movdet_est[Dind,"value"]))*1000000
  abund_move <- sum(denssurf_move) * 4
  lambda0_move <- exp(m_move$movdet_est[lind, c("value", "lower", "upper")])
  sigma_move <- exp(m_move$movdet_est[sind, c("value", "lower", "upper")])
  diffdense <- denssurf_stat - denssurf_move
  difabund <- sum(diffdense) *4
  
  #------------------------determine pdet for cutoff -----------------------------
  calc_pdet_x <- function(lambda0, sigma){
    #probability of not being detected during k 
    #nmesh rows, nocc cols
    notdetk <- apply(as.array(1:dim(useall)[2]), 1, function(k){
      enc <- lambda0 * 
        exp(- distmatscr^2 / (2 * sigma^2)) * 
        useall[,k]
      x = exp(-colSums(enc)) 
    }
    )
    
    #product over occ
    notdet <- apply(notdetk, 1, prod)
    #1 minus quantity
    pdet <- 1-notdet
    return(pdet) #vector length x in mesh
  }
  
  pdet_df <- data.frame(x = mesh$x,
                        y = mesh$y,
                        pdetx_stat = calc_pdet_x(lambda0_stat$value, sigma_stat$value),
                        pdetx_mov = calc_pdet_x(lambda0_move$value, sigma_move$value))
  pdet_df$pdet_stat <- sum(pdet_df$pdetx_stat * (denssurf_stat/sum(denssurf_stat)))
  pdet_df$pdet_mov <- sum(pdet_df$pdetx_mov * (denssurf_move/sum(denssurf_move)))
  pdet_df$relpdetx_stat <- pdet_df$pdetx_stat/pdet_df$pdet_stat
  pdet_df$relpdetx_mov <- pdet_df$pdetx_mov/pdet_df$pdet_mov
  pdet_df$relpdetdiff <- pdet_df$relpdetx_stat - pdet_df$relpdetx_mov
  
  if(is.null(subareacutoff)){
    subarea_stat <- subarea_mov <- subarea_both <- 1:ncol(distmat)
  } else {
    subarea_stat <- which(pdet_df$relpdetx_stat > subareacutoff)
    subarea_mov <- which(pdet_df$relpdetx_mov > subareacutoff)
    subarea_both <- unique(c(subarea_stat, subarea_mov))
  }
  
  #--------------------pdet plot -------------------------------------------------
  sharedmax <- max(c(pdet_df$relpdetx_mov, pdet_df$relpdetx_stat))
  pdetstatplot <- ggplot() +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 1) +
    geom_tile(pdet_df, 
              mapping = aes(x = x, y =y, fill = relpdetx_stat)) +
    geom_point(data = pdet_df[pdet_df$relpdetx_stat < 0.5,], 
               mapping = aes(x = x, y =y), color = "red", shape = 4) +
    scale_fill_viridis_c(limits = c(0, sharedmax),
                         name = "p.(x)/p.") +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 0.3) +
    coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
             ylim = c(min(mesh$y), max(mesh$y))) +
    ggtitle("Stationary Detector") +
    theme_bw()+
    theme(legend.position = "none",
          axis.title = element_blank())
  
  pdetmovplot <- ggplot() +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 1) +
    geom_tile(pdet_df, 
              mapping = aes(x = x, y =y, fill = relpdetx_mov)) +
    geom_point(data = pdet_df[pdet_df$relpdetx_mov < 0.5,],
               mapping = aes(x = x, y =y), color = "red", shape = 4) +
    scale_fill_viridis_c(limits = c(0, sharedmax),
                         name = "p.(x)/p.") +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 0.3) +
    coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
             ylim = c(min(mesh$y), max(mesh$y))) +
    ggtitle("Moving Detector") +
    theme_bw()+
    theme(legend.position = "none",
          axis.title = element_blank())
  
  pdetdiffplot <- ggplot() +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 1) +
    geom_tile(pdet_df, 
              mapping = aes(x = x, y =y, fill = relpdetdiff)) +
    scale_fill_viridis_c(limits = c(min(pdet_df$relpdetdiff), max(pdet_df$relpdetdiff)),
                         name = "p.(x)/p.\ndifference",
                         option = "magma") +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 0.3) +
    coord_sf(xlim = c(min(mesh$x), max(mesh$x)), 
             ylim = c(min(mesh$y), max(mesh$y))) +
    ggtitle("Stationary - Moving") +
    theme_bw()+
    theme(legend.position = "none",
          axis.title = element_blank())
  
  pdetlegendplot <- ggplot(pdet_df, 
                           mapping = aes(x = x, y = y, fill = relpdetx_stat)) +
    scale_fill_viridis_c(limits = c(0, sharedmax),
                         name = "p.(x)/p.") +
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
  
  pdetdifflegendplot <- ggplot(pdet_df, 
                               mapping = aes(x = x, y =y, fill = relpdetdiff)) +
    scale_fill_viridis_c(limits = c(min(pdet_df$relpdetdiff), max(pdet_df$relpdetdiff)),
                         name = "p.(x)/p.\nDifference",
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
  
  ggsave(file = paste("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/", label, "relpdet.png", sep = ""),
         plot = grid.arrange(
           grobs = list(pdetstatplot, pdetmovplot, pdetdiffplot,
                        pdetlegendplot, pdetdifflegendplot),
           widths = c((1), (1), (1)),
           heights = c(1,.2),
           layout_matrix = rbind(c(1,  2, 3),
                                 c(4, 4, 5))),
         width = 230,
         height = 230*.5,
         units = c("mm"),
         dpi = 300)
  
  #-----------------------------densities---------------------------------------
  denssurf_stat_forspreading <- denssurf_stat
  denssurf_stat_forspreading[-subarea_stat,] <- 0
  abund_stat_inpdot <- sum(denssurf_stat_forspreading) * 4
  denssurf_mov_forspreading <- denssurf_move
  denssurf_mov_forspreading[-subarea_mov,] <- 0
  abund_move_inpdot <- sum(denssurf_mov_forspreading) * 4
  Dspreadsurf_stat <- spreadD(meshdistmat, lambda0_stat$value, sigma_stat$value, denssurf_stat_forspreading)
  Dspreadsurf_move <- spreadD(meshdistmat, lambda0_move$value, sigma_move$value, denssurf_mov_forspreading)
  Dspread_diff <- Dspreadsurf_stat - Dspreadsurf_move
  
  
  lowcolor = "#440154FF"#"#000004FF" 
  colorm1 =  "#3B528BFF" #"#51127CFF" 
  colorm2 = "#21908CFF"#"#B63679FF"
  colorm3 = "#5DC863FF" #"#FB8861FF"
  highcolor = "#FDE725FF"#"#FCFDBFFF"
  #-----------------------AC density---------------------------------------------  
  
  sharedmax = max(c(denssurf_stat[subarea_both],
                    denssurf_move[subarea_both]
  ))
  sharedmids = c(sharedmax*(1/5), sharedmax*(2/5), sharedmax*(3/5),  sharedmax*(4/5))
  
  ACDstatplot <- ggplot() +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 1) +
    geom_tile(data.frame(x = mesh$x, y = mesh$y, D = denssurf_stat)[subarea_stat,]
              , 
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
    theme(legend.position = "none",
          axis.title = element_blank())
  
  ACDmovplot <- ggplot() +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 1) +
    geom_tile(data.frame(x = mesh$x, y = mesh$y, D = denssurf_move)[subarea_mov,], 
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
    theme(legend.position = "none",
          axis.title = element_blank())
  
  ACDlegendplot <- ggplot(data = data.frame(x = mesh$x, y = mesh$y, D = denssurf_move)[subarea_both,], 
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
    geom_tile(data.frame(x = mesh$x, y = mesh$y, D = diffdense)[subarea_both,], 
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
    theme(legend.position = "none",
          axis.title = element_blank())
  
  ACDdifflegendplot <- ggplot(data.frame(x = mesh$x, y = mesh$y, D = diffdense)[subarea_both,], 
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
  
  
  ggsave(file = paste("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/", label, "ACDplot.png", sep = ""),
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
  #----------------------------- animal D ---------------------------------------  
  
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
    theme(legend.position = "none",
          text = element_text(size = 42),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
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
    theme(legend.position = "none",
          text = element_text(size = 42),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
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
          text = element_text(size = 42),
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
    theme(legend.position = "none",
          text = element_text(size = 42),
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank())
  
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
    theme_bw() + 
    theme(axis.title = element_blank())
  
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
          text = element_text(size = 42),
          panel.border = element_blank())
  
  ggsave(file = paste("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/", label, "animDplot.png", sep = ""),
         plot = grid.arrange(
           grobs = list(animDstat, animDmov, animDdiff,
                        animDlegend, animDdifflegend),
           widths = c((1), (1), (1)),
           heights = c(1,.2),
           layout_matrix = rbind(c(1,  2, 3),
                                 c(4, 4, 5))),
         width = 700,
         height = 700*.5,
         units = c("mm"),
         dpi = 300)
  
  #-----------------detection parameters ---------------------------------------
  
  detpars <- data.frame(lambda0 = unlist(c(lambda0_stat, lambda0_move)),
                        sigma = unlist(c(sigma_stat, sigma_move)),
                        name = names(c(lambda0_stat, lambda0_move)),
                        model = c(rep("stationary", 3), rep("moving", 3)))
  detdat <- do.call(rbind, lapply(as.list(1:6), function(n){
    df=
      data.frame(x = seq(0,4*max(detpars$sigma), length.out = 40))
    
    df$y = detpars$lambda0[n]*exp(-df$x^2/(2*detpars$sigma[n]^2))
    
    df$name = detpars$name[n]
    df$model = detpars$model[n]
    return(df)}))
  detdatwide <-  pivot_wider(detdat, values_from = y, names_from = name)
  
  detfctplot <- ggplot() +
    geom_line(detdat[detdat$name %in% c("value"),],
              mapping = aes(x = x/1000, y = y*1000, color = model, group = model),
              linewidth = 1.5) +
    geom_ribbon(detdatwide,
                mapping = aes(x = x/1000, ymin = lower*1000, ymax = upper*1000, 
                              fill = model, group = model),
                color = "transparent",
                alpha = .5) +
    geom_point(data = data.frame(x = detpars[detpars$name %in% c("value"),"sigma"]/1000,
                                 y = 0, 
                                 model = detpars[detpars$name %in% c("value"),"model"]),  
               mapping = aes(group = model, color = model, x = x, y = y),
               size = 3, shape = 4, stroke = 1.5) +
    ylab("Detection rate\n(detections per km)") +
    xlab("Distance (km)")+
    scale_color_manual(name = "", 
                       values =  c("#5F187FFF", "#F8765CFF"),
                       labels = c("Moving", "Stationary")) +
    scale_fill_manual(name = "", 
                      values =  c("#5F187FFF", "#F8765CFF"),
                      labels = c("Moving", "Stationary")) +
    theme_bw() +
    theme(text = element_text(size = 32))
  
  ggsave(file = paste("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/", label, "detplot.png", sep = ""),
         plot = detfctplot,
         width = 169,
         height = 169*.7,
         units = c("mm"),
         dpi = 300)
  
  return(list(statN = abund_stat,
              movN = abund_move,
              statNpdet = abund_stat_inpdot,
              movNpdet = abund_move_inpdot))
}  

#-------visualize track direction ----------------------------------------------
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
tracklines_ls <- create_line_list_C(tracksdf, scenario = "onison")

plotsurv <- function(occi){
  survplot <- ggplot() +
    geom_sf(data = st_as_sf(lpoly), mapping = aes(), fill = "#93c0d3", col = "#93c0d3",
            linewidth = .1, alpha = 0.3) +
    geom_segment(tracklines_ls[[occi]], mapping = aes(x = x1, y = y1,
                                                      xend = x2, yend = y2)
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

ggsave(file = paste("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/surveydirection.png", sep = ""),
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

create_plots(m_move = myfits[[6]],
             DdesignX = Xmats[[6]],
             label = paste(formulas[[6]])[3],
             subareacutoff = 0.5)
create_plots(m_move = myfits[[6]],
             DdesignX = Xmats[[6]],
             m0 = fits[[4]],
             label = paste(paste(formulas[[6]])[3], "nocutoff"),
             subareacutoff = 0)
create_plots(m_move = myfits[[2]],
             DdesignX = Xmats[[2]],
             label = paste(paste(formulas[[2]])[3]),
             subareacutoff = 0.5)
create_plots(m_move = myfits[[2]],
             DdesignX = Xmats[[2]],
             m0 = fits[[3]],
             label = paste(paste(formulas[[2]])[3], "nocutoff"),
             subareacutoff = 0)




colSums(usage(trapscr))

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", 
                "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
baseplot <-  ggplot() +
  geom_sf(data = st_as_sf(lpoly), fill = "#93c0d3", col = "#9aaac4") +
  coord_sf(xlim = c(min(traps$x), max(traps$x)), 
           ylim = c(min(traps$y), max(traps$y))) +
  theme_bw() +
  theme(text=element_text(size=2),
        plot.title=element_blank(), 
        axis.text.x = element_text(angle = 90),
        axis.ticks.length=unit(-0.1, "cm")
  )
#avg and max detections per detected individual,

detsperind <- apply(capthist, 1, sum)
dpi <- ggplot() +
  geom_bar(data.frame(x = detsperind),
           mapping =  aes(x = x),
           fill = cbbPalette[2], color = cbbPalette[2],
           alpha = .8) +
  ylab("Count") +
  ggtitle("B") +
  xlab("Detections per individual") +
  scale_y_continuous(limits = c(0, 250)) + 
  scale_x_continuous(breaks = c(0, 1, 2)) +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        plot.title=element_text(margin = margin(l = 4, b = -15)))

#number of different traps/detectors at which individuals were detected, 
trapsperind <- rowSums(apply(capthist, c(1,3), max))
tpi <- ggplot() +
  geom_bar(data.frame(x = trapsperind),mapping =  aes(x = x),
           fill = cbbPalette[2], color = cbbPalette[2],
           alpha = .8) +
  scale_y_continuous(limits = c(0, 250)) + 
  scale_x_continuous(breaks = c(0, 1, 2)) +
  xlab("Traps per individual") +
  ggtitle("C") +
  theme_bw() +
  theme(axis.title.y = element_blank(),
        plot.title=element_text(margin = margin(l = 4,b = -15))) 

totaleffort <- rowSums(usage(trapscr))/1000
#and number of individuals detected at a detector
indspertrap <- apply(capthist, 3, sum)/totaleffort

indspertrap[(indspertrap==0)]<- NA
trapdat <- cbind(traps, indspertrap)
cuts <- c(.01,  .5,  1, 1.5,  2,  7, 27)


ipt <- baseplot +
  #coord_sf(xlim = c(min(trapdat$x), max(trapdat$x)), 
  #         ylim = c(min(trapdat$y), max(trapdat$y))) +
  geom_point(trapdat, mapping = aes(x = x, y = y, fill = indspertrap), 
             size = 2, shape = 21) +
  scale_fill_stepsn(colours = viridisLite::magma(length(cuts)-1),
                    breaks = cuts,
                    values = scales::rescale(cuts, to = c(0,1), from = range(cuts)), 
                    limits = c(min(cuts), max(cuts)),
                    labels = c(cuts),
                    name = "",
                    na.value = "grey") +
  ggtitle("A") +
  theme_void() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(.865,.16),
        #axis.text.x = element_text(angle = 90),
        legend.background = element_rect(fill= "transparent", color = NA),
        plot.title=element_text(margin = margin(l = 4, b = -15)))


ggsave(file = paste("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/CHS/CHS_results/capthistsummary.png", sep = ""),
       plot = grid.arrange(
         grobs = list(dpi, tpi, ipt),
         widths = c(1,1),
         heights = c(1,.9),
         layout_matrix = rbind(c(3,1),c(3,2))
       ),
       width = 169*.8,
       height = 169*.8,
       units = c("mm"),
       dpi = 300)


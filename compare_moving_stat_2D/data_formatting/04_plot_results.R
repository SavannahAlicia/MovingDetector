library(dplyr)


create_plots <- function(sim_fits_out, 
                         Dmodel = "variable", 
                         plotcols = c("#5F187FFF", "#F8765CFF", "black"),#c("cornflowerblue", "goldenrod", "black"),
                         linesize = 1,
                         pointsize = 2,
                         fontsize = 24,
                         output){
  ##---converged models
  statconv <- unlist(lapply(sim_fits_out, function(x){
    if(is.null(x$stat_conv)){
      out = NA
    } else {
      out = x$stat_conv
    }
    })) == 0
  movconv <- unlist(lapply(sim_fits_out, function(x){
    if(is.null(x$mov_conv)){
      out = NA
    } else {
      out = x$mov_conv
    }
    })) == 0
  bothconv <- (statconv + movconv) == 2
  bothconv[is.na(bothconv)] <- FALSE
  
  sim_fits_out <- sim_fits_out[bothconv]
  
  ###------------------------compare computation time-----------------------------
  
  times <- do.call(rbind, lapply(as.list(1:length(sim_fits_out)), FUN = function(x){
    stattime = as.numeric(sim_fits_out[[x]]$statdet_time)
    movetime =as.numeric(sim_fits_out[[x]]$movdet_time)
    #should add something to exclude times from failed Hessians
    return(data.frame(stat = stattime, move = movetime))
  }))
  
  timeplot <- ggplot() +
    geom_boxplot(data = tidyr::pivot_longer(times, cols = c("stat", "move")),
                 mapping = aes(y = log(value), group = name, color = name),
                 linewidth = linesize) +
    scale_color_manual(name = "Model", labels = c("Moving", "Stationary"),
                                                  values = plotcols[1:2]) +
    ylab("Log seconds") +
    xlab("Model") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
      text = element_text(size = fontsize),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "none")
  ###--------------------------compare estimates ---------------------------------
  make_plot_dat<- function(sim_fits_out){
    
    stat_outs <- do.call(rbind,lapply(as.list(1:length(sim_fits_out)), FUN = function(x){
      df <- sim_fits_out[[x]]$statdet_est
      trueN = sim_fits_out[[x]]$n #just detected inds though

      # if(any(is.na(df))) {
      #   df <- df %>% 
      #     mutate(across(-name, ~NA)) #don't include results where hessian failed
      # }
      df$sim = rep(x,nrow(df))
      df$value = as.numeric(df$value)
      df$upper = as.numeric(df$upper)
      df$lower = as.numeric(df$lower)
      return(df)
    }))
    move_outs <-  do.call(rbind,lapply(as.list(1:length(sim_fits_out)), FUN = function(x){
      df <- sim_fits_out[[x]]$movdet_est
      trueN = sim_fits_out[[x]]$n
  
      # if(any(is.na(df))) {
      #   df <- df %>% 
      #     mutate(across(-name, ~NA)) #don't include results where hessian failed
      # }
      df$sim = rep(x,nrow(df))
      
      df$value = as.numeric(df$value)
      df$upper = as.numeric(df$upper)
      df$lower = as.numeric(df$lower)
      return(df)
    }))
    
    all_outs <- rbind(cbind(stat_outs, data.frame(model = rep("stationary", 
                                                              nrow(stat_outs)))),
                      cbind(move_outs, data.frame(model = rep("moving",nrow(move_outs)))))
    
    if(Dmodel == "flat"){
      Ndat <- all_outs[all_outs$name == "D",]
      Ndat$value <- log(exp(Ndat$value) * meshspacing^2 * nrow(mesh))
      Ndat$upper <- log(exp(Ndat$upper) * meshspacing^2 * nrow(mesh))
      Ndat$lower <- log(exp(Ndat$lower) * meshspacing^2 * nrow(mesh))
      Ndat$name = "N"
      all_outs <- rbind(all_outs, Ndat)
    }
    
    all_outs2 <- all_outs %>%
      group_by(name, model) %>%
      summarize(mean = mean(value, na.rm = T), 
                median = median(value, na.rm = T),
                meanupper = quantile(value, probs = .975, na.rm = T), 
                meanlower = quantile(value, probs = .025, na.rm = T))
    
    if (Dmodel == "variable"){

      meshstep = meshspacing/3
      newmeshxys = expand.grid(seq(min(mesh$x), max(mesh$x), meshstep),
                                seq(min(mesh$x), max(mesh$x), meshstep))

      
      #calculate D for meshx (row) and estimated betas (column)
      enames = sim_fits_out[[1]]$statdet_est$name
      
      statDdat <- apply(as.array(1:length(sim_fits_out)), 1, function(sim){
        calcDv(newmeshxys[,1],
               newmeshxys[,2],
        beta1_ = beta1, #sim_fits_out[[sim]]$statdet_est$value[enames == "beta1"],
        beta2_ = sim_fits_out[[sim]]$statdet_est$value[enames == "beta2"],
        N_ = exp(sim_fits_out[[sim]]$statdet_est$value[enames == "N"]),
        meshspacing)
      })
      moveDdat <- apply(as.array(1:length(sim_fits_out)), 1, function(sim){
          calcDv(newmeshxys[,1],
                 newmeshxys[,2],
          beta1_ = beta1, #sim_fits_out[[sim]]$movdet_est$value[enames == "beta1"],
          beta2_ = sim_fits_out[[sim]]$movdet_est$value[enames == "beta2"],
          N_ = exp(sim_fits_out[[sim]]$movdet_est$value[enames == "N"]),
          meshspacing)
      }) 

      #dataframe with x, y, trueD, stationarydets, movingdets, quantile columns
      D_plotdat <- data.frame(x = rep(newmeshxys[,1],3),
                              trueD = rep(
                                          calcDv(newmeshxys[,1], 
                                                 newmeshxys[,2], 
                                                 beta1, beta2, N,
                                                 meshspacing), 
                                          3),
                              stationarydets = c(rowMeans(statDdat),
                                                 apply(statDdat, 1, function(x) quantile(x, probs = .025, na.rm = T)),
                                                 apply(statDdat, 1, function(x) quantile(x, probs = .975, na.rm = T))
                                                 ),
                              movingdets = c(rowMeans(moveDdat),
                                             apply(moveDdat, 1, function(x) quantile(x, probs = .025, na.rm = T)),
                                             apply(moveDdat, 1, function(x) quantile(x, probs = .975, na.rm = T))
                                             ),
                              quantile = c(rep("mean", nrow(newmeshxys)),
                                           rep("2.5%",nrow(newmeshxys)),
                                           rep("97.5%",nrow(newmeshxys)))
      )
      #mean and quantiles for sum D
      moveNdat <- colSums(moveDdat) * meshspacing^2
     
      statNdat <- colSums(statDdat) * meshspacing^2
      
      Nplotdat <- data.frame(model = rep(c("True", "Moving", "Stationary"), 3),
                             quantile = c("mean", "mean", "mean",
                                          "2.5%", "2.5%", "2.5%",
                                          "97.5%", "97.5%", "97.5%"),
                             value = c(N,
                                       mean(moveNdat), 
                                       mean(statNdat),
                                       NA,
                                       quantile(moveNdat, probs = 0.025, na.rm = T),
                                       quantile(statNdat, probs = 0.025, na.rm = T),
                                       NA,
                                       quantile(moveNdat, probs = 0.975, na.rm = T),
                                       quantile(statNdat, probs = 0.975, na.rm = T))
                             )
    } else if(Dmodel == "flat"){
      D_plotdat <- data.frame(x = rep(mesh$x, 3),
                              y = rep(mesh$y, 3),
                              trueD = rep(D_mesh_f, 3),
                              stationarydets = c(rep(exp(as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                                all_outs2$name == "D", "mean"])), nrow(mesh)),
                                                 rep(exp(as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                                all_outs2$name == "D", "meanupper"])), nrow(mesh)),
                                                 rep(exp(as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                                all_outs2$name == "D", "meanlower"])), nrow(mesh))),
                              movingdets = c(rep(exp(as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                            all_outs2$name == "D", "mean"])), nrow(mesh)),
                                             rep(exp(as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                            all_outs2$name == "D", "meanupper"])), nrow(mesh)),
                                             rep(exp(as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                            all_outs2$name == "D", "meanlower"])), nrow(mesh)) 
                              ),
                              quantile = c(rep("mean", nrow(mesh)), rep("2.5%", nrow(mesh)), rep("97.5%",nrow(mesh)))
      )
      
      Nplotdat <- data.frame(model = rep(c("True", "Moving", "Stationary"), 3),
                             quantile = c("mean", "mean", "mean",
                                          "2.5%", "2.5%", "2.5%",
                                          "97.5%", "97.5%", "97.5%"),
                             value = c(flatD * nrow(mesh) * meshspacing^2,
                             exp(as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                        all_outs2$name == "D", "mean"]))* nrow(mesh) * meshspacing^2, 
                             exp(as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                        all_outs2$name == "D", "mean"]))* nrow(mesh) * meshspacing^2,
                                       NA,
                             exp(as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                        all_outs2$name == "D", "meanlower"]))* nrow(mesh) * meshspacing^2, 
                             exp(as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                        all_outs2$name == "D", "meanlower"]))* nrow(mesh) * meshspacing^2,
                                       NA,
                             exp(as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                        all_outs2$name == "D", "meanupper"]))* nrow(mesh) * meshspacing^2, 
                             exp(as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                        all_outs2$name == "D", "meanupper"]))* nrow(mesh) * meshspacing^2)
                             )
    }


    out <- list(all_outs = all_outs, 
                all_outs2 = all_outs2,
                D_plotdat = D_plotdat,
                Nplotdat = Nplotdat)
  }
  plotdat <- make_plot_dat(sim_fits_out)
  all_outs <- plotdat$all_outs
  all_outs2 <- plotdat$all_outs2
  D_plotdat <- plotdat$D_plotdat
  Nplotdat <- plotdat$Nplotdat
  
  ###--------------------------compare precision -------------------------------
  
  ###--------------------------create plots ------------------------------------
  
  lambda0plot <- 
    ggplot() +
    geom_pointrange(all_outs2[all_outs2$name == "lambda0",],
                    mapping = aes(ymin = exp(meanlower)*1000, 
                        ymax = exp(meanupper)*1000,
                        y = exp(mean)*1000,
                        col = model,
                        x = model),
                    size = pointsize/2,
                    linewidth = linesize) +
    # geom_density(all_outs[all_outs$name == "lambda0",], 
    #              mapping = aes(x = exp(value)*1000, #per km instead of m
    #                            col = model, fill = model),
    #              alpha = 0.5, size = linesize) +
    geom_hline(data = rbind( all_outs2[all_outs2$name == "lambda0", ], 
                               data.frame(name = "lambda0", model = "true", 
                                          mean = log(lambda0))), 
                 aes(yintercept = exp(c(mean))*1000,
                     col = model), size = linesize,
               linetype = "dashed") +
    # geom_vline(data = all_outs2[all_outs2$name == "lambda0",], 
    #            aes(xintercept = exp(c(meanlower))*1000, 
    #                col = model), 
    #            linetype = "dashed", size = linesize) +
    # geom_vline(data = all_outs2[all_outs2$name == "lambda0",], 
    #            aes(xintercept = exp(c(meanupper))*1000, 
    #                col = model),
    #            linetype = "dashed", size = linesize) +
    # #geom_vline(xintercept = lambda0*1000, size = linesize, col = "black") +
    xlab("Model") +
     ylab(expression(paste("\u03bb"[0], " (dets/km)"))) + 
    scale_x_discrete(labels = c("moving" = "Moving",
                                "stationary" = "Stationary")) +
  
    # xlim(0,lambda0*4*1000) +
    ylim(0,lambda0*10*1000) +
    scale_color_manual(name = "",
                       values = plotcols, 
                       labels = c("Moving", "Stationary", 
                                  expression("True \u03bb"[0]))) +
    scale_fill_manual(name = "",
                       values = plotcols, 
                       labels = c("Moving", "Stationary", 
                                  expression("True \u03bb"[0]))) +
    guides(fill = "none") +
    theme_classic() +
    theme(#axis.title = element_text(size = 10),
          legend.position = "none",
          #axis.text.y = element_blank(),
          text = element_text(size = fontsize),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank())
  #legend.title = element_text(size = 20),
  #legend.text = element_text(size = 20))
  

  sigmaplot <- 
    ggplot() +
    geom_pointrange(all_outs2[all_outs2$name == "sigma",],
                    mapping = aes(ymin = exp(meanlower)/1000, 
                                  ymax = exp(meanupper)/1000,
                                  y = exp(mean)/1000,
                                  col = model,
                                  x = model),
                    size = pointsize/2,
                    linewidth = linesize) +
    geom_hline(data = rbind(all_outs2[all_outs2$name == "sigma",],
                            data.frame(name = "sigma", model = "true",
                                       mean = log(sigma))),
               aes(yintercept = exp(mean)/1000,
                   col = model), 
               linetype = "dashed", size = linesize) +
    scale_color_manual(name = "",
                       labels = c("Moving", "Stationary", "True \u03C3"),
                       values = plotcols) +
    scale_y_continuous(limits = c(0*sigma/1000, 1.5*sigma/1000)) +
    scale_x_discrete(labels = c("moving" = "Moving",
                                "stationary" = "Stationary")) +
    ylab("\u03C3 (km)") +
    xlab("Model") +
    theme_classic() +
    theme(#axis.title = element_text(size = 10),
          legend.position = "none",
          text = element_text(size = fontsize),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank()
          #axis.text.y = element_blank()
          )
  #legend.title = element_text(size = 20),
  #legend.text = element_text(size = 20))
  
  # beta1plot <- ggplot() +
  #   geom_density(all_outs[all_outs$name == "beta1",], 
  #                mapping = aes(x = value, col = model), size = linesize) +
  #   geom_vline(data = all_outs2[all_outs2$name == "beta1",], 
  #              aes(xintercept = c(mean), col = model), size = linesize) +
  #   geom_vline(data = all_outs2[all_outs2$name == "beta1",], 
  #              aes(xintercept = c(meanlower), col = model), linetype = "dashed", size = linesize) +
  #   geom_vline(data = all_outs2[all_outs2$name == "beta1",], 
  #              aes(xintercept = c(meanupper), col = model), linetype = "dashed", size = linesize) +
  #   geom_vline(xintercept = beta1) +
  #   ylab("Frequency") +
  #   scale_color_manual(values = plotcols) +
  #   xlab("beta1") +
  #   scale_x_continuous(limits = c(beta1*4/3, beta1*1/3)) +
  #   theme_classic() +
  #   theme(axis.title = element_text(size = 10),
  #         axis.text.y = element_blank(),
  #         axis.title.y = element_blank(),
  #         legend.position = "none",
  #         legend.title = element_text(size = 20),
  #         legend.text = element_text(size = 20))
  
  beta2plot <- ggplot() +
    geom_density(all_outs[all_outs$name == "beta2",],
                 mapping = aes(x = value/1000, col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta2",], 
               aes(xintercept = c(mean)/1000, col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta2",], 
               aes(xintercept = c(meanlower)/1000, col = model), linetype = "dashed", size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta2",],
               aes(xintercept = c(meanupper)/1000, col = model), linetype = "dashed", size = linesize) +
    geom_vline(xintercept = beta2/1000) +
    scale_color_manual(values = plotcols) +
    xlab("beta2") +
    ylab("Frequency") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.position = "none",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20))
  
  
  Nplot <- ggplot() +
    geom_pointrange(all_outs2[all_outs2$name == "N",],
                    mapping = aes(ymin = exp(meanlower), 
                                  ymax = exp(meanupper),
                                  y = exp(mean),
                                  col = model,
                                  x = model),
                    size = pointsize/2,
                    linewidth = linesize) +
     geom_hline(data = rbind(all_outs2[all_outs2$name == "N",],
                             data.frame(name = "N", 
                                        model = "true",
                                        mean = log(N))),
               aes(yintercept = exp(mean), col = model), 
               linetype = "dashed",
               size = linesize) +
    scale_y_continuous(limits = c(0*N, 1.1*N)) +
    scale_x_discrete(labels = c("moving" = "Moving", "stationary" = "Stationary")) +
    scale_color_manual(name = "",
                      values = plotcols, 
                      labels = c("Moving", "Stationary", 
                                 expression("True"))) +
    xlab("Model") +
    ylab("N") +
    theme_classic() +
    theme(#axis.title = element_text(size = 10),
          #axis.text.y = element_blank(),
          #axis.title.y = element_blank(),
          #legend.position = "none",
          #legend.title = element_text(size = 20),
          #legend.text = element_text(size = 20),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      axis.ticks.x = element_blank(),
      text = element_text(size = fontsize)
      )
  

  # 
  # ggplot() + 
  #   geom_tile(data = D_plotdat[,c(1,2,3)], aes(x = x, y= y, fill = trueD)) +
  #   scale_color_manual(values = "black",  name = "D") +
  #   theme_classic() +
  #   labs(title = "True density") +
  #   theme(axis.text.y = element_blank()
  #         ,
  #         legend.position = "none") 
  # ggplot() + 
  #   geom_tile(data = D_plotdat[,c(1,2,4)], aes(x = x, y= y, fill = stationarydets)) +
  #   scale_fill_viridis_c(name = "D") +
  #   theme_classic() +
  #   labs(title = "Stationary detectors") +
  #   theme(axis.text.y = element_blank(),
  #         legend.position = "none") 
  # ggplot() + 
  #   geom_tile(data = D_plotdat[,c(1,2,5)], aes(x = x, y= y, fill = movingdets)) +
  #   scale_fill_viridis_c(name = "D") +
  #   theme_classic() +
  #   labs(title = "Moving detectors") +
  #   theme(axis.text.y = element_blank(),
  #         legend.position = "none") 
  
  
  D_plotdatlong <- tidyr::pivot_longer(D_plotdat, 
                                       cols = c("trueD", "stationarydets", "movingdets"))
  

  Dplot <- 
    ggplot() + 
    geom_line(data = D_plotdatlong[D_plotdatlong$quantile == "mean",], 
              mapping = aes(x = x/1000, y = value*1000000, col = name
                            ),linewidth = linesize/2) +
    geom_point(data = D_plotdatlong[which(D_plotdatlong$quantile == "2.5%" &
                                            D_plotdatlong$name != "trueD"),], 
              mapping = aes(x = x/1000, y = value*1000000, 
                            col = name), 
              #linetype = "dashed", 
              size = pointsize*2,#/2,
              shape = 18
              ) +
    geom_point(data = D_plotdatlong[which(D_plotdatlong$quantile == "97.5%" & 
                                            D_plotdatlong$name != "trueD"),], 
              mapping = aes(x = x/1000, y = value*1000000, col = name),
              #linetype = "dashed", 
              size = pointsize*2,#/2,
              shape = 16
              ) +
    scale_color_manual(values = c(plotcols, "black"), 
                       labels = c("Moving", "Stationary", "True"),
                       name = "") +
    #ylim(0,.5) +
    xlim(c(-beta2 - 1000, - beta2 + 1000)/1000)+
    #ylim(c(0, 400)) +
    ylab("D") +
    xlab("x") +
    theme_classic() +
    guides(
      color = guide_legend(
        override.aes = list(
          shape = NA,      # remove points
          linewidth = linesize
        )
      )
    ) +
    theme(#axis.title = element_text(size = 10),
          text = element_text(size = fontsize),
          #axis.text = element_text(size = 10),
          #legend.title = element_text(size = 10),
          #legend.text = element_text(size = 10)
          )
  
  if (Dmodel == "variable"){
    plotlist = list(lambda0plot, sigmaplot, 
                    #beta1plot,
                    beta2plot,
                    Dplot,
                    Nplot, timeplot)
    out = grid.arrange(
      grobs = plotlist,
      widths = c(1,1),
      heights = c(1,1,1,1),
      layout_matrix = rbind(c(1,6),
                            c(2,4),
                            c(7, 7),
                            5))
  }
  if(Dmodel == "flat"){
    plotlist <- list(lambda0plot, 
                          sigmaplot, 
                          Nplot,#Dplot,
                          timeplot)
    out = grid.arrange(
      grobs = plotlist,
      widths = c(1,1),
      heights = c(1,1,1),
      layout_matrix = rbind(c(1,3),
                            c(2,4)))
  }

  if(output == "plots"){
    return(out)
  } else if(output == "plotdat"){
    return(plotdat)
  } else if(output == "plotlist"){
    return(plotlist)
  }
  
}

# Check and create the directory
if (!dir.exists(paste(dirstart, "/plots", sep = ""))) {
  dir.create(paste(dirstart, "/plots", sep = ""), recursive = TRUE)
}

all_sim_fits_q <- readRDS(paste(dirstart, "variable_dens.Rds", sep = ""))
all_sim_fits <- readRDS(paste(dirstart, "flat_dens.Rds", sep = ""))

vpl <- create_plots(all_sim_fits_q, Dmodel = "variable", 
                    pointsize = .5,#1.5, 
                    fontsize = 8,
                    linesize = .5,#2,
                    output = "plotdat")
fpl <- create_plots(all_sim_fits, Dmodel = "flat", 
                    pointsize = .5,#1.5, 
                    linesize = 2, output = "plotdat")
ggsave(file = paste(dirstart, "plots/variable_moving_2D.png", sep = ""),
       plot = vpl,
       width = 169,
       height = 169*(1/2)*(4/3),
       units = c("mm"),
       dpi = 300)
ggsave(file = paste(dirstart, "plots/flat_moving_2D.png", sep = ""),
       plot = fpl,
       width = 169,
       height = 169*(1/2),
       units = c("mm"),
       dpi = 300)

ggsave(file = paste(dirstart, "plots/setup.png", sep = ""),
       plot = ggplot() +
         geom_raster(data.frame(x = mesh$x, y = mesh$y, D = D_mesh_v), 
                     mapping = aes(x = x, y = y, fill = D)) +
         geom_point(data = tracksdf, mapping = aes(x = x, y = y, group = occ),
                    size = 1.5,shape = "+", color= "white") +
         scale_fill_viridis_c() +
         scale_x_continuous(expand = c(0,0),
                            breaks = c(0, 1000, 2000, 3000, 4000, 5000),
                            labels = c("0","1","2",3,4,5)) +
         scale_y_continuous(expand = c(0,0),
                            breaks = c(0, 1000, 2000, 3000)) +
         geom_segment(aes(x = 1500, y = 2500, xend = 3500, yend = 2500),
                      arrow = arrow(length = unit(0.3, "cm")),
                      color = "grey") +
         #xlim(1000, 4000) +
         geom_point(data.frame(x = traps$x,
                               y = traps$y),
                    mapping = aes(x = x, y = y), shape = 21, , color = "pink", fill = "red", alpha = .7, size = 1.4) +
         theme_bw() +
         theme(axis.title = element_blank(),
               axis.text.y = element_blank()),
       width = 169,
       height = 169*(1/2),
       units = c("mm"),
       dpi = 300)

# plotdat_q <-create_plots(all_sim_fits_q, Dmodel = "variable", output = "plotdat")
# all_outs_q <- plotdat_q$all_outs
# c(t.test(all_outs_q[all_outs_q$name == "lambda0" & all_outs_q$model == "moving",]$sd, 
#        all_outs_q[all_outs_q$name == "lambda0" & all_outs_q$model == "stationary",]$sd,
#        paired = T, alternative = "less")$p.value,
# t.test(all_outs_q[all_outs_q$name == "sigma" & all_outs_q$model == "moving",]$sd, 
#        all_outs_q[all_outs_q$name == "sigma" & all_outs_q$model == "stationary",]$sd,
#        paired = T, alternative = "less")$p.value,
# t.test(all_outs_q[all_outs_q$name == "beta1" & all_outs_q$model == "moving",]$sd, 
#        all_outs_q[all_outs_q$name == "beta1" & all_outs_q$model == "stationary",]$sd, 
#        paired = T, alternative = "less")$p.value,
# t.test(all_outs_q[all_outs_q$name == "beta2" & all_outs_q$model == "moving",]$sd, 
#        all_outs_q[all_outs_q$name == "beta2" & all_outs_q$model == "stationary",]$sd, 
#        paired = T, alternative = "less")$p.value
# )
# c(
#   mean((all_outs_q[all_outs_q$name == "lambda0" & all_outs_q$model == "moving",]$sd^2-
#           all_outs_q[all_outs_q$name == "lambda0" & all_outs_q$model == "stationary",]$sd^2)/all_outs_q[all_outs_q$name == "lambda0" & all_outs_q$model == "stationary",]$sd^2),
#   mean((all_outs_q[all_outs_q$name == "sigma" & all_outs_q$model == "moving",]$sd^2-
#           all_outs_q[all_outs_q$name == "sigma" & all_outs_q$model == "stationary",]$sd^2)/all_outs_q[all_outs_q$name == "sigma" & all_outs_q$model == "stationary",]$sd^2),
#   mean((all_outs_q[all_outs_q$name == "beta1" & all_outs_q$model == "moving",]$sd^2-
#           all_outs_q[all_outs_q$name == "beta1" & all_outs_q$model == "stationary",]$sd^2)/all_outs_q[all_outs_q$name == "beta1" & all_outs_q$model == "stationary",]$sd^2),
#   mean((all_outs_q[all_outs_q$name == "beta2" & all_outs_q$model == "moving",]$sd^2-
#           all_outs_q[all_outs_q$name == "beta2" & all_outs_q$model == "stationary",]$sd^2)/all_outs_q[all_outs_q$name == "beta2" & all_outs_q$model == "stationary",]$sd^2)
#   
# )
# 
# 
# 
# plotdat_1 <-create_plots(all_sim_fits, Dmodel = "flat", output = "plotdat")
# all_outs_1 <- plotdat_1$all_outs
# c(t.test(all_outs_1[all_outs_1$name == "lambda0" & all_outs_1$model == "moving",]$sd, 
#        all_outs_1[all_outs_1$name == "lambda0" & all_outs_1$model == "stationary",]$sd,
#        paired = T, alternative = "less")$p.value,
# t.test(all_outs_1[all_outs_1$name == "sigma" & all_outs_1$model == "moving",]$sd, 
#        all_outs_1[all_outs_1$name == "sigma" & all_outs_1$model == "stationary",]$sd,
#        paired = T, alternative = "less")$p.value,
# t.test(all_outs_1[all_outs_1$name == "D" & all_outs_1$model == "moving",]$sd, 
#        all_outs_1[all_outs_1$name == "D" & all_outs_1$model == "stationary",]$sd, 
#        paired = T, alternative = "less")$p.value
# )
# 
# c(
#   mean((all_outs_1[all_outs_1$name == "lambda0" & all_outs_1$model == "moving",]$sd^2-
#           all_outs_1[all_outs_1$name == "lambda0" & all_outs_1$model == "stationary",]$sd^2)/all_outs_1[all_outs_1$name == "lambda0" & all_outs_1$model == "stationary",]$sd^2),
#   mean((all_outs_1[all_outs_1$name == "sigma" & all_outs_1$model == "moving",]$sd^2-
#           all_outs_1[all_outs_1$name == "sigma" & all_outs_1$model == "stationary",]$sd^2)/all_outs_1[all_outs_1$name == "sigma" & all_outs_1$model == "stationary",]$sd^2),
#   mean((all_outs_1[all_outs_1$name == "D" & all_outs_1$model == "moving",]$sd^2-
#           all_outs_1[all_outs_1$name == "D" & all_outs_1$model == "stationary",]$sd^2)/all_outs_1[all_outs_1$name == "D" & all_outs_1$model == "stationary",]$sd^2)
# )
# 

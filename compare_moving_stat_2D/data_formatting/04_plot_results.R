library(dplyr)


create_plots <- function(sim_fits_out, Dmodel = "variable",
                         plotcols = c("cornflowerblue", "goldenrod", "black"),
                         linesize = .3, output = "plots"){
  ###------------------------compare computation time-----------------------------
  
  times <- do.call(rbind, lapply(as.list(1:length(sim_fits_out)), FUN = function(x){
    stattime = as.numeric(sim_fits_out[[x]]$statdet_time)
    movetime = as.numeric(sim_fits_out[[x]]$movdet_time)
    return(data.frame(stat = stattime, move = movetime))
  }))
  
  timeplot <- ggplot() +
    geom_boxplot(data = tidyr::pivot_longer(times, cols = c("stat", "move")),
                 mapping = aes(y = log(value), group = name, color = name)) +
    scale_color_manual(name = "Model", labels = c("Moving", "Stationary"),
                                                  values = plotcols[1:2]) +
    ylab("Log seconds") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "none")
  ###--------------------------compare estimates ---------------------------------
  make_plot_dat<- function(sim_fits_out){
    
    stat_outs <- do.call(rbind,lapply(as.list(1:length(sim_fits_out)), FUN = function(x){
      df <- sim_fits_out[[x]]$statdet_est
      df$sim = rep(x,nrow(sim_fits_out[[1]]$statdet_est))
      return(df)
    }))
    move_outs <-  do.call(rbind,lapply(as.list(1:length(sim_fits_out)), FUN = function(x){
      df <- sim_fits_out[[x]]$movdet_est
      df$sim = rep(x,nrow(sim_fits_out[[1]]$statdet_est))
      return(df)
    }))
    
    all_outs <- rbind(cbind(stat_outs, data.frame(model = rep("stationary", 
                                                              nrow(stat_outs)))),
                      cbind(move_outs, data.frame(model = rep("moving",nrow(move_outs)))))
    all_outs$sd <- (((all_outs$upper)- (all_outs$value))/1.96)
    
    all_outs2 <- all_outs %>%
      group_by(name, model) %>%
      summarize(mean = mean(value), 
                median = median(value),
                meanupper = quantile(value, probs = .975), 
                meanlower = quantile(value, probs = .025),
                meansd = mean(sd))
    
    out <- list(all_outs= all_outs, all_outs2 =all_outs2)
  }
  plotdat <- make_plot_dat(sim_fits_out)
  all_outs <- plotdat$all_outs
  all_outs2 <- plotdat$all_outs2
  
  ###--------------------------compare precision -------------------------------
  
  ###--------------------------create plots ------------------------------------
  
  lambda0plot <- 
    ggplot() +
    geom_density(all_outs[all_outs$name == "lambda0",], 
                 mapping = aes(x = exp(value)*1000, #per km instead of m
                               col = model), size = linesize) +
    geom_vline(data = rbind( all_outs2[all_outs2$name == "lambda0", ], 
                             data.frame(name = "lambda0", model = "true", 
                                        mean = log(lambda0))), 
               aes(xintercept = exp(c(mean))*1000,
                   col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "lambda0",], 
               aes(xintercept = exp(c(meanlower))*1000, 
                   col = model), 
               linetype = "dashed", size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "lambda0",], 
               aes(xintercept = exp(c(meanupper))*1000, 
                   col = model),
               linetype = "dashed", size = linesize) +
    geom_vline(xintercept = lambda0*1000, size = linesize, col = "black") +
    xlab(expression("\u03bb"[0])) +
    ylab("Frequency") + 
    #xlim(.004,.006) +
    scale_color_manual(name = "",
                       values = plotcols, 
                       labels = c("Moving", "Stationary", 
                                  expression("True \u03bb"[0]))) +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          legend.position = "none",
          axis.text.y = element_blank())
  #legend.title = element_text(size = 20),
  #legend.text = element_text(size = 20))
  
  lambda0precisionplot <- 
    ggplot() +
    geom_density(all_outs[all_outs$name == "lambda0",], 
                 mapping = aes(x = invlogit(upper)- invlogit(lower), col = model),
                 size = linesize) +
    geom_vline(data = rbind( all_outs2[all_outs2$name == "lambda0", ], 
                             data.frame(name = "lambda0", model = "true", 
                                        mean = logit(lambda0))), 
               aes(xintercept = invlogit(c(mean)), col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "lambda0",], 
               aes(xintercept = invlogit(c(meanlower)), col = model), 
               linetype = "dashed", size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "lambda0",], 
               aes(xintercept = invlogit(c(meanupper)), col = model),
               linetype = "dashed", size = linesize) +
    geom_vline(xintercept = lambda0, size = linesize, col = "black") +
    xlab(expression("\u03bb"[0])) +
    ylab("Frequency") + 
    #xlim(.004,.006) +
    scale_color_manual(name = "",
                       values = plotcols, 
                       labels = c("Moving", "Stationary", 
                                  expression("True \u03bb"[0]))) +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          legend.position = "none",
          axis.text.y = element_blank())
  
  
  sigmaplot <- 
    ggplot() +
    geom_density(all_outs[all_outs$name == "sigma",], 
                 mapping = aes(x = exp(value)/1000, #in km instead of m
                               col = model), size = linesize) +
    geom_vline(data = rbind(all_outs2[all_outs2$name == "sigma",],
                            data.frame(name = "sigma", model = "true",
                                       mean = log(sigma))),
               aes(xintercept = exp(mean)/1000,
                   col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "sigma",],
               aes(xintercept = exp(meanlower)/1000,
                   col = model), 
               linetype = "dashed", size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "sigma",], 
               aes(xintercept = exp(meanupper)/1000,
                   col = model), 
               linetype = "dashed", size = linesize) +
    geom_vline(xintercept = sigma/1000, size = linesize, col = "black")+
    scale_color_manual(name = "",
                       labels = c("Moving", "Stationary", "True \u03C3"),
                       values = plotcols) +
    xlab("\u03C3") +
    ylab("Frequency") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          legend.position = "none",
          axis.title.y = element_blank(),
          axis.text.y = element_blank())
  #legend.title = element_text(size = 20),
  #legend.text = element_text(size = 20))
  
  beta1plot <- ggplot() +
    geom_density(all_outs[all_outs$name == "beta1",], 
                 mapping = aes(x = value, col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta1",], 
               aes(xintercept = c(mean), col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta1",], 
               aes(xintercept = c(meanlower), col = model), linetype = "dashed", size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta1",], 
               aes(xintercept = c(meanupper), col = model), linetype = "dashed", size = linesize) +
    geom_vline(xintercept = beta1) +
    ylab("Frequency") +
    scale_color_manual(values = plotcols) +
    xlab("beta1") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          axis.text.y = element_blank(),
          legend.position = "none",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20))
  
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
          legend.position = "none",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20))
  
  meshstep <- meshspacing/10
  if(Dmodel == "variable"){
    D_plotdat <- data.frame(x = rep(seq(min(mesh$x), max(mesh$x), meshstep),3),
                            y = rep(rep(mesh$y[1], length(seq(min(mesh$x), max(mesh$x), meshstep))),3),
                            trueD = rep(exp(beta1*(seq(min(mesh$x), max(mesh$x), meshstep) + beta2)^2), 3),
                            stationarydets = c(exp(as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                          all_outs2$name == "beta1", "mean"]) * 
                                                     (seq(min(mesh$x), max(mesh$x), meshstep) + 
                                                        as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                               all_outs2$name == "beta2", "mean"]))^2),
                                               exp(as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                          all_outs2$name == "beta1", "meanupper"]) * 
                                                     (seq(min(mesh$x), max(mesh$x), meshstep) + 
                                                        as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                               all_outs2$name == "beta2", "meanupper"]))^2),
                                               exp(as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                          all_outs2$name == "beta1", "meanlower"]) * 
                                                     (seq(min(mesh$x), max(mesh$x), meshstep) + 
                                                        as.numeric(all_outs2[all_outs2$model == "stationary" & 
                                                                               all_outs2$name == "beta2", "meanlower"]))^2)),
                            movingdets = c(exp(as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                      all_outs2$name == "beta1", "mean"]) * 
                                                 (seq(min(mesh$x), max(mesh$x), meshstep) + 
                                                    as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                           all_outs2$name == "beta2", "mean"]))^2),
                                           exp(as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                      all_outs2$name == "beta1", "meanupper"]) * 
                                                 (seq(min(mesh$x), max(mesh$x), meshstep) + 
                                                    as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                           all_outs2$name == "beta2", "meanupper"]))^2),
                                           exp(as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                      all_outs2$name == "beta1", "meanlower"]) * 
                                                 (seq(min(mesh$x), max(mesh$x), meshstep) + 
                                                    as.numeric(all_outs2[all_outs2$model == "moving" & 
                                                                           all_outs2$name == "beta2", "meanlower"]))^2)),
                            quantile = c(rep("mean", length(seq(min(mesh$x), max(mesh$x), meshstep))),
                                         rep("2.5%",length(seq(min(mesh$x), max(mesh$x), meshstep))),
                                         rep("97.5%",length(seq(min(mesh$x), max(mesh$x), meshstep))))
    )
    
  } else if(Dmodel == "flat"){
    D_plotdat <- data.frame(x = rep(mesh$x, 3),
                            y = rep(mesh$y, 3),
                            trueD = rep(D_mesh, 3),
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
    
  } else {
    stop("Please specify 'flat' or 'variable' for Dmodel")
  }
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
  
  
  D_plotdatlong <- tidyr::pivot_longer(D_plotdat, cols = c("trueD", "stationarydets", "movingdets"))
  
  Dplot <- 
    ggplot() + 
    geom_line(data = D_plotdatlong[D_plotdatlong$quantile == "mean",], mapping = aes(x = x/1000, y = value, col = name,
                                                                                     linewidth = name)) +
    geom_line(data = D_plotdatlong[D_plotdatlong$quantile == "2.5%",], mapping = aes(x = x/1000, y = value, col = name), linetype = "dashed", size = linesize) +
    geom_line(data = D_plotdatlong[D_plotdatlong$quantile == "97.5%",], mapping = aes(x = x/1000, y = value, col = name), linetype = "dashed", size = linesize) +
    scale_color_manual(values = c(plotcols, "black"), labels = c("Moving", "Stationary", "True"),
                       name = "") +
    scale_linewidth_manual(values = c(linesize*3, linesize*3, linesize), 
                           labels = c("Moving", "Stationary", "True"), 
                           name = "") +
    #ylim(0,.5) +
    xlim(c(-beta2 - 1000, - beta2 + 1000)/1000)+
    ylab("AC density") +
    xlab("x") +
    theme_classic() +
    guides(linewidth = "none") +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10))
  
  if (Dmodel == "variable"){
    out = grid.arrange(
      grobs = list(lambda0plot, sigmaplot, 
                   beta1plot, beta2plot,
                   Dplot,
                   timeplot),
      widths = c(1,1),
      heights = c(1,1,1,1),
      layout_matrix = rbind(c(1,2),
                            c(3,4),
                            5,
                            6))
  }
  if(Dmodel == "flat"){
    out = grid.arrange(
      grobs = list(lambda0plot, sigmaplot, 
                   
                   Dplot,
                   timeplot),
      widths = c(1,1),
      heights = c(1,1,1),
      layout_matrix = rbind(c(1,2),
                            3,
                            4))
  }
  if(output == "plots"){
    return(out)
  } else if(output == "plotdat")
    return(plotdat)
  
}

all_sim_fits_q <- readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/simulation_results/variable_dens.Rds")
all_sim_fits <- readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/simulation_results/flat_dens.Rds")

vpl <- create_plots(all_sim_fits_q, Dmodel = "variable")
fpl <- create_plots(all_sim_fits, Dmodel = "flat")
ggsave(file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/simulation_results/plots/variable_moving_2D.png",
       plot = vpl,
       width = 169,
       height = 169*(1/2)*(4/3),
       units = c("mm"),
       dpi = 300)
ggsave(file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/simulation_results/plots/flat_moving_2D.png",
       plot = fpl,
       width = 169,
       height = 169*(1/2),
       units = c("mm"),
       dpi = 300)

ggsave(file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_2D/simulation_results/plots/setup.png",
       plot = ggplot() +
         geom_raster(data.frame(x = mesh$x, y = mesh$y, D = D_mesh_q), 
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

plotdat_q <-create_plots(all_sim_fits_q, Dmodel = "variable", output = "plotdat")
all_outs_q <- plotdat_q$all_outs
c(t.test(all_outs_q[all_outs_q$name == "lambda0" & all_outs_q$model == "moving",]$sd, 
       all_outs_q[all_outs_q$name == "lambda0" & all_outs_q$model == "stationary",]$sd,
       paired = T, alternative = "less")$p.value,
t.test(all_outs_q[all_outs_q$name == "sigma" & all_outs_q$model == "moving",]$sd, 
       all_outs_q[all_outs_q$name == "sigma" & all_outs_q$model == "stationary",]$sd,
       paired = T, alternative = "less")$p.value,
t.test(all_outs_q[all_outs_q$name == "beta1" & all_outs_q$model == "moving",]$sd, 
       all_outs_q[all_outs_q$name == "beta1" & all_outs_q$model == "stationary",]$sd, 
       paired = T, alternative = "less")$p.value,
t.test(all_outs_q[all_outs_q$name == "beta2" & all_outs_q$model == "moving",]$sd, 
       all_outs_q[all_outs_q$name == "beta2" & all_outs_q$model == "stationary",]$sd, 
       paired = T, alternative = "less")$p.value
)
c(
  mean((all_outs_q[all_outs_q$name == "lambda0" & all_outs_q$model == "moving",]$sd^2-
          all_outs_q[all_outs_q$name == "lambda0" & all_outs_q$model == "stationary",]$sd^2)/all_outs_q[all_outs_q$name == "lambda0" & all_outs_q$model == "stationary",]$sd^2),
  mean((all_outs_q[all_outs_q$name == "sigma" & all_outs_q$model == "moving",]$sd^2-
          all_outs_q[all_outs_q$name == "sigma" & all_outs_q$model == "stationary",]$sd^2)/all_outs_q[all_outs_q$name == "sigma" & all_outs_q$model == "stationary",]$sd^2),
  mean((all_outs_q[all_outs_q$name == "beta1" & all_outs_q$model == "moving",]$sd^2-
          all_outs_q[all_outs_q$name == "beta1" & all_outs_q$model == "stationary",]$sd^2)/all_outs_q[all_outs_q$name == "beta1" & all_outs_q$model == "stationary",]$sd^2),
  mean((all_outs_q[all_outs_q$name == "beta2" & all_outs_q$model == "moving",]$sd^2-
          all_outs_q[all_outs_q$name == "beta2" & all_outs_q$model == "stationary",]$sd^2)/all_outs_q[all_outs_q$name == "beta2" & all_outs_q$model == "stationary",]$sd^2)
  
)



plotdat_1 <-create_plots(all_sim_fits, Dmodel = "flat", output = "plotdat")
all_outs_1 <- plotdat_1$all_outs
c(t.test(all_outs_1[all_outs_1$name == "lambda0" & all_outs_1$model == "moving",]$sd, 
       all_outs_1[all_outs_1$name == "lambda0" & all_outs_1$model == "stationary",]$sd,
       paired = T, alternative = "less")$p.value,
t.test(all_outs_1[all_outs_1$name == "sigma" & all_outs_1$model == "moving",]$sd, 
       all_outs_1[all_outs_1$name == "sigma" & all_outs_1$model == "stationary",]$sd,
       paired = T, alternative = "less")$p.value,
t.test(all_outs_1[all_outs_1$name == "D" & all_outs_1$model == "moving",]$sd, 
       all_outs_1[all_outs_1$name == "D" & all_outs_1$model == "stationary",]$sd, 
       paired = T, alternative = "less")$p.value
)

c(
  mean((all_outs_1[all_outs_1$name == "lambda0" & all_outs_1$model == "moving",]$sd^2-
          all_outs_1[all_outs_1$name == "lambda0" & all_outs_1$model == "stationary",]$sd^2)/all_outs_1[all_outs_1$name == "lambda0" & all_outs_1$model == "stationary",]$sd^2),
  mean((all_outs_1[all_outs_1$name == "sigma" & all_outs_1$model == "moving",]$sd^2-
          all_outs_1[all_outs_1$name == "sigma" & all_outs_1$model == "stationary",]$sd^2)/all_outs_1[all_outs_1$name == "sigma" & all_outs_1$model == "stationary",]$sd^2),
  mean((all_outs_1[all_outs_1$name == "D" & all_outs_1$model == "moving",]$sd^2-
          all_outs_1[all_outs_1$name == "D" & all_outs_1$model == "stationary",]$sd^2)/all_outs_1[all_outs_1$name == "D" & all_outs_1$model == "stationary",]$sd^2)
)


## 1 dimension
library(dplyr)

create_plots <- function(sim_fits_out, mesh, D_mesh,
                         sim_fits_outlin, meshlin, D_meshlin, 
                         Dmodel, 
                         plotcols = c("#721F81FF","#F1605DFF","black" ),
                         linesize = .3, output = "plots"){
  ###------------------------compare computation time-----------------------------
  
  times <- do.call(rbind, lapply(as.list(1:length(sim_fits_out)), FUN = function(x){
    lintime = as.numeric(sim_fits_outlin[[x]]$movdet_time)
    twoDtime = as.numeric(sim_fits_out[[x]]$movdet_time)
    return(data.frame(twoD = twoDtime, lin = lintime))
  }))
  
  timeplot <- ggplot() +
    geom_boxplot(data = tidyr::pivot_longer(times, cols = c("twoD", "lin")),
                 mapping = aes(y = log(value), group = name, color = name)) +
    scale_color_manual(name = "Model", labels = c("1D", "2D"),
                       values = plotcols[1:2]) +
    ylab("Log seconds") +
    xlab("") +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          legend.position = "none"
          )
  
  ###--------------------------compare estimates ---------------------------------
  make_plot_dat<- function(sim_fits_out, sim_fits_outlin){
    
    lin_outs <- do.call(rbind,lapply(as.list(1:length(sim_fits_outlin)), FUN = function(x){
      df <- sim_fits_outlin[[x]]$movdet_est
      df$sim = rep(x,nrow(sim_fits_outlin[[1]]$movdet_est))
      return(df)
    }))
    twoD_outs <-  do.call(rbind,lapply(as.list(1:length(sim_fits_out)), FUN = function(x){
      df <- sim_fits_out[[x]]$movdet_est
      df$sim = rep(x,nrow(sim_fits_out[[1]]$movdet_est))
      return(df)
    }))
    
    all_outs <- rbind(cbind(lin_outs, data.frame(model = rep("1D", 
                                                              nrow(lin_outs)))),
                      cbind(twoD_outs, data.frame(model = rep("2D",nrow(twoD_outs)))))
    all_outs$sd <- (((all_outs$upper)- (all_outs$value))/1.96)
    
    all_outs2 <- all_outs %>%
      dplyr::group_by(name, model) %>%
      dplyr::summarize(mean = mean(value), 
                median = median(value),
                meanupper = quantile(value, probs = .975), 
                meanlower = quantile(value, probs = .025),
                meansd = mean(sd, na.rm = T))
    
    out <- list(all_outs= all_outs, all_outs2 =all_outs2)
  }
  plotdat <- make_plot_dat(sim_fits_out, sim_fits_outlin)
  all_outs <- plotdat$all_outs
  all_outs2 <- plotdat$all_outs2
  
  
  ###--------------------------create plots ------------------------------------
  
  lambda0plot <- 
    ggplot() +
    geom_density(all_outs[all_outs$name == "lambda0",], 
                 mapping = aes(x = exp(value), col = model), size = linesize) +
    geom_vline(data = rbind( all_outs2[all_outs2$name == "lambda0", ], 
                             data.frame(name = "lambda0", model = "true", 
                                        mean = logit(lambda0))), 
               aes(xintercept = exp(c(mean)), col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "lambda0",], 
               aes(xintercept = exp(c(meanlower)), col = model), 
               linetype = "dashed", size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "lambda0",], 
               aes(xintercept = exp(c(meanupper)), col = model),
               linetype = "dashed", size = linesize) +
    geom_vline(xintercept = lambda0, size = linesize, col = "black") +
    xlab(expression("\u03bb"[0])) +
    ylab("Frequency") + 
    #xlim(.004,.006) +
    scale_color_manual(name = "",
                       values = plotcols, 
                       labels = c("1D", "2D", 
                                  expression("True \u03bb"[0]))
                       ) +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          legend.position = "none",
          axis.text.y = element_blank())
  #legend.title = element_text(size = 20),
  #legend.text = element_text(size = 20))
  
  lambda0precisionplot <- 
    ggplot() +
    geom_density(all_outs[all_outs$name == "lambda0",], 
                 mapping = aes(x = sd, col = model),
                 size = linesize) +
    xlab(expression("\u03bb"[0])) +
    ylab("Frequency") + 
    #xlim(.004,.006) +
    scale_color_manual(name = "",
                       values = plotcols, 
                       labels = c("1D", "2D", 
                                  expression("True \u03bb"[0]))) +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          legend.position = "none",
          axis.text.y = element_blank())
  
  
  sigmaplot <- 
    ggplot() +
    geom_density(all_outs[all_outs$name == "sigma",], 
                 mapping = aes(x = exp(value), col = model), size = linesize) +
    geom_vline(data = rbind(all_outs2[all_outs2$name == "sigma",],
                            data.frame(name = "sigma", model = "true",
                                       mean = log(sigma))),
               aes(xintercept = exp(mean), col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "sigma",],
               aes(xintercept = exp(meanlower), col = model), 
               linetype = "dashed", size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "sigma",], 
               aes(xintercept = exp(meanupper), col = model), 
               linetype = "dashed", size = linesize) +
    geom_vline(xintercept = sigma, size = linesize, col = "black")+
    scale_color_manual(name = "",
                       labels = c("1D", "2D", "True \u03C3"),
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
    geom_density(all_outs[all_outs$name == "beta1",], mapping = aes(x = value, col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta1",], aes(xintercept = c(mean), col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta1",], aes(xintercept = c(meanlower), col = model), linetype = "dashed", size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta1",], aes(xintercept = c(meanupper), col = model), linetype = "dashed", size = linesize) +
    geom_vline(xintercept = beta1) +
    ylab("Frequency") +
    scale_color_manual(values = plotcols) +
    xlab("beta1") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          legend.position = "none",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20))
  beta2plot <- ggplot() +
    geom_density(all_outs[all_outs$name == "beta2",], mapping = aes(x = value, col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta2",], aes(xintercept = c(mean), col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta2",], aes(xintercept = c(meanlower), col = model), linetype = "dashed", size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta2",], aes(xintercept = c(meanupper), col = model), linetype = "dashed", size = linesize) +
    geom_vline(xintercept = beta2) +
    scale_color_manual(values = plotcols) +
    xlab("beta2") +
    ylab("Frequency") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          legend.position = "none",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20))
  
  
  meshstep <- meshspacing
  if(Dmodel == "variable"){
  D_plotdat <- data.frame(x = c(rep(mesh$x,3)),
                          y = c(rep(mesh$y, 3)),
                          trueD = rep(D_mesh, 3),
                          twoDdets = c(
                                               exp(as.numeric(all_outs2[all_outs2$model == "2D" & 
                                                                          all_outs2$name == "beta1", "mean"]) * 
                                                     (mesh$x + 
                                                        as.numeric(all_outs2[all_outs2$model == "2D" & 
                                                                               all_outs2$name == "beta2", "mean"]))^2),
                                          
                                               exp(as.numeric(all_outs2[all_outs2$model == "2D" & 
                                                                          all_outs2$name == "beta1", "meanupper"]) * 
                                                     (mesh$x + 
                                                        as.numeric(all_outs2[all_outs2$model == "2D" & 
                                                                               all_outs2$name == "beta2", "meanupper"]))^2),
                                             
                                               exp(as.numeric(all_outs2[all_outs2$model == "2D" & 
                                                                          all_outs2$name == "beta1", "meanlower"]) * 
                                                     (mesh$x + 
                                                        as.numeric(all_outs2[all_outs2$model == "2D" & 
                                                                               all_outs2$name == "beta2", "meanlower"]))^2)),
                          lindets = c(
                                           exp(as.numeric(all_outs2[all_outs2$model == "1D" & 
                                                                      all_outs2$name == "beta1", "mean"]) * 
                                                 (mesh$x + 
                                                    as.numeric(all_outs2[all_outs2$model == "1D" & 
                                                                           all_outs2$name == "beta2", "mean"]))^2),
                                    
                                           exp(as.numeric(all_outs2[all_outs2$model == "1D" & 
                                                                      all_outs2$name == "beta1", "meanupper"]) * 
                                                 (mesh$x + 
                                                    as.numeric(all_outs2[all_outs2$model == "1D" & 
                                                                           all_outs2$name == "beta2", "meanupper"]))^2),
                                        
                                           exp(as.numeric(all_outs2[all_outs2$model == "1D" & 
                                                                      all_outs2$name == "beta1", "meanlower"]) * 
                                                 (mesh$x + 
                                                    as.numeric(all_outs2[all_outs2$model == "1D" & 
                                                                           all_outs2$name == "beta2", "meanlower"]))^2)),
                          quantile = c(rep("mean", length(mesh$x)),
                                       rep("2.5%",length(mesh$x)),
                                       rep("97.5%",length(mesh$x)))
  )
  
  D_plotdatlin <- data.frame(x = c(rep(meshlin$x,3)),
                              y = c(rep(meshlin$y, 3)),
                              trueD = rep(D_meshlin, 3),
                              twoDdets = c(
                                exp(as.numeric(all_outs2[all_outs2$model == "2D" & 
                                                           all_outs2$name == "beta1", "mean"]) * 
                                      (meshlin$x + 
                                         as.numeric(all_outs2[all_outs2$model == "2D" & 
                                                                all_outs2$name == "beta2", "mean"]))^2 +
                                      as.numeric(all_outs2[all_outs2$model == "2D" & 
                                                             all_outs2$name == "beta3", "mean"])),
                                
                                exp(as.numeric(all_outs2[all_outs2$model == "2D" & 
                                                           all_outs2$name == "beta1", "meanupper"]) * 
                                      (meshlin$x + 
                                         as.numeric(all_outs2[all_outs2$model == "2D" & 
                                                                all_outs2$name == "beta2", "meanupper"]))^2 +
                                      as.numeric(all_outs2[all_outs2$model == "2D" & 
                                                             all_outs2$name == "beta3", "meanupper"])),
                                
                                exp(as.numeric(all_outs2[all_outs2$model == "2D" & 
                                                           all_outs2$name == "beta1", "meanlower"]) * 
                                      (meshlin$x + 
                                         as.numeric(all_outs2[all_outs2$model == "2D" & 
                                                                all_outs2$name == "beta2", "meanlower"]))^2 +
                                      as.numeric(all_outs2[all_outs2$model == "2D" & 
                                                             all_outs2$name == "beta3", "meanlower"]))),
                              lindets = c(
                                exp(as.numeric(all_outs2[all_outs2$model == "1D" & 
                                                           all_outs2$name == "beta1", "mean"]) * 
                                      (meshlin$x + 
                                         as.numeric(all_outs2[all_outs2$model == "1D" & 
                                                                all_outs2$name == "beta2", "mean"]))^2 +
                                      as.numeric(all_outs2[all_outs2$model == "1D" & 
                                                             all_outs2$name == "beta3", "mean"])),
                                
                                exp(as.numeric(all_outs2[all_outs2$model == "1D" & 
                                                           all_outs2$name == "beta1", "meanupper"]) * 
                                      (meshlin$x + 
                                         as.numeric(all_outs2[all_outs2$model == "1D" & 
                                                                all_outs2$name == "beta2", "meanupper"]))^2 +
                                      as.numeric(all_outs2[all_outs2$model == "1D" & 
                                                             all_outs2$name == "beta3", "meanupper"])),
                                
                                exp(as.numeric(all_outs2[all_outs2$model == "1D" & 
                                                           all_outs2$name == "beta1", "meanlower"]) * 
                                      (meshlin$x + 
                                         as.numeric(all_outs2[all_outs2$model == "1D" & 
                                                                all_outs2$name == "beta2", "meanlower"]))^2 +
                                      as.numeric(all_outs2[all_outs2$model == "1D" & 
                                                             all_outs2$name == "beta3", "meanlower"]))),
                              quantile = c(rep("mean", length(meshlin$x)),
                                           rep("2.5%",length(meshlin$x)),
                                           rep("97.5%",length(meshlin$x)))
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
  )}
  
  ggplot() + 
    geom_tile(data = D_plotdatlin[D_plotdatlin$quantile == "mean",c(1,2,3)], aes(x = x, y= y, fill = trueD)) +
    scale_fill_viridis_c(name = "D") +
    theme_classic() +
    labs(title = "True density") +
    theme(axis.text.y = element_blank()
          ,
         # legend.position = "none"
          ) 
  ggplot() + 
    geom_tile(data = D_plotdatlin[D_plotdatlin$quantile == "mean",c(1,2,4)], aes(x = x, y= y, fill = twoDdets*streamwidth/1000)) +
    scale_fill_viridis_c(name = "D") +
    theme_classic() +
    labs(title = "2D") +
    theme(axis.text.y = element_blank(),
          #legend.position = "none"
          ) 
  ggplot() + 
    geom_tile(data = D_plotdatlin[D_plotdatlin$quantile == "mean",c(1,2,5)],
              aes(x = x, y= y, fill = lindets)) +
    scale_fill_viridis_c(name = "D") +
    theme_classic() +
    labs(title = "1D") +
    theme(axis.text.y = element_blank(),
          #legend.position = "none"
          ) 
  
  
  D_plotdatlong <- tidyr::pivot_longer(D_plotdat, cols = c("trueD", "twoDdets", "lindets"))
  
  Dplot <- 
    ggplot() + 
    geom_line(data = D_plotdatlong[D_plotdatlong$quantile == "mean",], mapping = aes(x = x, y = value, col = name), size = linesize) +
    geom_line(data = D_plotdatlong[D_plotdatlong$quantile == "2.5%",], mapping = aes(x = x, y = value, col = name), linetype = "dashed", size = linesize) +
    geom_line(data = D_plotdatlong[D_plotdatlong$quantile == "97.5%",], mapping = aes(x = x, y = value, col = name), linetype = "dashed", size = linesize) +
    scale_color_manual(values = c(plotcols, "black"), labels = c("1D", "2D", "True"),
                       name = "") +
    #ylim(0,.5) +
   # xlim(1500,3000)+
    ylab("AC density") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          axis.text = element_text(size = 10),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10))
  
  Dbyx <- unique(D_plotdat[D_plotdat$quantile == "mean",][order(D_plotdat[which(D_plotdat$quantile == "mean"), "x"]),][,-which(colnames(D_plotdat) == "y")])
  ggplot() +
    geom_line(Dbyx, mapping = aes(x = x, y = twoDdets)) +
    geom_line(Dbyx, mapping = aes(x = x, y = lindets))
  
  #ggsave(file = paste("writing_up/flatdens.png", sep = ""),
  # plot =
  grid.arrange(
    grobs = list(lambda0plot, sigmaplot, Dplot),
    widths = c(1,1),
    heights = c(1,1),
    layout_matrix = rbind(c(1,2), 3)
  )
  #,
  #  width = 169,
  #  height = 169*(1/2),
  #  units = c("mm"),
  #  dpi = 300)
  
  grid.arrange(
    grobs = list(lambda0plot, sigmaplot, beta1plot, beta2plot,Dplot),
    widths = c(1,1),
    heights = c(1,1,1),
    layout_matrix = rbind(c(1,2),c(3,4),5)
  )
}

sim_fits_out <- readRDS("~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/simulation_results/all_sim_fits1.Rds")

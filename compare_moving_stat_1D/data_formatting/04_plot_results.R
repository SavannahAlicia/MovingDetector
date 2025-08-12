## 1 dimension
library(dplyr)
library(boot)
library(ggplot2)
library(ggborderline)
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
                 mapping = aes(y = log(value), 
                               group = name, color = name),
                 linewidth = linesize) +
    scale_color_manual(name = "Model", labels = c("1D", "2D"),
                       values = plotcols[1:2]) +
    ylab("Log seconds") +
    xlab("") +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(color = "darkgrey", linewidth = linesize/2),
          panel.grid.minor.y = element_line(color = "darkgrey", linewidth = linesize/4),
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
                      cbind(twoD_outs, data.frame(model = rep("2D",
                                                              nrow(twoD_outs)))))
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
                                        mean = log(lambda0))), 
               aes(xintercept = exp(c(mean)), col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "lambda0",], 
               aes(xintercept = exp(c(meanlower)), col = model), 
               linetype = "dashed", size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "lambda0",], 
               aes(xintercept = exp(c(meanupper)), col = model),
               linetype = "dashed", size = linesize) +
    geom_vline(xintercept = lambda0, linetype = "dotted",
               size = linesize, col = "black") +
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
    geom_vline(xintercept = beta1, size = linesize) +
    ylab("Frequency") +
    scale_color_manual(values = plotcols) +
    xlab("beta1") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          legend.position = "none",
          axis.text.y = element_blank(),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20))
  beta2plot <- ggplot() +
    geom_density(all_outs[all_outs$name == "beta2",], mapping = aes(x = value, col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta2",], aes(xintercept = c(mean), col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta2",], aes(xintercept = c(meanlower), col = model), linetype = "dashed", size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta2",], aes(xintercept = c(meanupper), col = model), linetype = "dashed", size = linesize) +
    geom_vline(xintercept = beta2, size = linesize) +
    scale_color_manual(values = plotcols) +
    xlab("beta2") +
    ylab("Frequency") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          legend.position = "none",
          legend.title = element_text(size = 20),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          legend.text = element_text(size = 20))
  beta3plot <- ggplot() +
    geom_density(all_outs[all_outs$name == "beta3",], mapping = aes(x = value, col = model), size = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta3",], aes(xintercept = c(mean), col = model), linewidth = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta3",], aes(xintercept = c(meanlower), col = model), linetype = "dashed", linewidth = linesize) +
    geom_vline(data = all_outs2[all_outs2$name == "beta3",], aes(xintercept = c(meanupper), col = model), linetype = "dashed", linewidth = linesize) +
    geom_vline(xintercept = beta3, size = linesize) +
    scale_color_manual(values = plotcols) +
    xlab("beta3") +
    ylab("Frequency") +
    theme_classic() +
    theme(axis.title = element_text(size = 10),
          legend.position = "none",
          axis.text.y = element_blank(),
          #axis.title.y = element_blank(),
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 20))
  
  
  meshstep <- meshspacing
  
  if(Dmodel == "variable"){
    
  gv <- function(model, name, quantile)  {
    as.numeric(all_outs2[all_outs2$model == model & 
                           all_outs2$name == name, quantile])
  }
  pv <- function(meshx, model, quantile){
    exp(gv(model, "beta1", quantile) * (meshx + gv(model, "beta2", quantile))^2 + gv(model, "beta3", quantile))
  }
  pv_b <- function(meshx, beta1, beta2, beta3){
    exp(beta1 * (meshx + beta2)^2 + beta3)
  }
  D_plotdat <- data.frame(x = c(rep(mesh$x,3)),
                          y = c(rep(mesh$y, 3)),
                          trueD = rep(pv_b(mesh$x, beta1, beta2, beta3), 3),
                          trueDlin = rep(pv_b(mesh$x, beta1, beta2, beta3+log(streamwidth/1000)), 3),
                          twoDdets = c(pv(mesh$x, "2D", "mean"),
                                       pv(mesh$x, "2D", "meanupper"),
                                       pv(mesh$x, "2D", "meanlower")),
                          lindets = c(pv(mesh$x, "1D", "mean"),
                                      pv(mesh$x, "1D", "meanupper"),
                                      pv(mesh$x, "1D",  "meanlower")),
                          quantile = c(rep("mean", length(mesh$x)),
                                       rep("2.5%",length(mesh$x)),
                                       rep("97.5%",length(mesh$x)))
  )
  
  D_plotdatlin <- data.frame(x = c(rep(meshlin$x,3)),
                              y = c(rep(meshlin$y, 3)),
                             trueD = rep(pv_b(meshlin$x, beta1, beta2, beta3), 3),
                              trueDlin = rep(pv_b(meshlin$x, beta1, beta2, beta3+log(streamwidth/1000)), 3),
                             twoDdets = c(pv(meshlin$x, "2D", "mean"),
                                          pv(meshlin$x, "2D", "meanupper"),
                                          pv(meshlin$x, "2D", "meanlower")),
                             lindets = c(pv(meshlin$x, "1D", "mean"),
                                         pv(meshlin$x, "1D", "meanupper"),
                                         pv(meshlin$x, "1D",  "meanlower")),
                              quantile = c(rep("mean", length(meshlin$x)),
                                           rep("2.5%",length(meshlin$x)),
                                           rep("97.5%",length(meshlin$x)))
                              )
  
  } else if(Dmodel == "flat"){
  gvf <- function(model, quantile){
    exp(as.numeric(all_outs2[all_outs2$model == model & 
                               all_outs2$name == "D", quantile]))
  }  
  D_plotdat <- data.frame(x = rep(mesh$x, 3),
                          y = rep(mesh$y, 3),
                          trueD = rep(flatD, 3),
                          trueDlin = rep(flatD*(streamwidth/1000), 3),
                          twoDdets = c(rep(gvf("2D", "mean"), length(mesh$x)),
                                       rep(gvf("2D", "meanupper"), length(mesh$x)),
                                       rep(gvf("2D", "meanlower"), length(mesh$x))),
                          lindets = c(rep(gvf("1D", "mean"), length(mesh$x)),
                                      rep(gvf("1D", "meanupper"), length(mesh$x)),
                                      rep(gvf("1D", "meanlower"), length(mesh$x))),
                          quantile = c(rep("mean", nrow(mesh)), 
                                       rep("2.5%", nrow(mesh)), 
                                       rep("97.5%",nrow(mesh)))
  )
  
  D_plotdatlin <- data.frame(x = rep(meshlin$x, 3),
                          y = rep(meshlin$y, 3),
                          trueD = rep(flatD, 3),
                          trueDlin = rep(flatD*(streamwidth/1000), 3),
                          twoDdets = c(rep(gvf("2D", "mean"), length(meshlin$x)),
                                       rep(gvf("2D", "meanupper"), length(meshlin$x)),
                                       rep(gvf("2D", "meanlower"), length(meshlin$x))),
                          lindets = c(rep(gvf("1D", "mean"), length(meshlin$x)),
                                      rep(gvf("1D", "meanupper"), length(meshlin$x)),
                                      rep(gvf("1D", "meanlower"), length(meshlin$x))),
                          quantile = c(rep("mean", nrow(meshlin)), 
                                       rep("2.5%", nrow(meshlin)), 
                                       rep("97.5%",nrow(meshlin)))
  )
  
  }
  
  ggplot() + 
    geom_tile(data = D_plotdatlin[D_plotdatlin$quantile == "mean",c(1,2,4)], 
              aes(x = x, y= y, fill = trueDlin)) +
    scale_fill_viridis_c(name = "D", limits = c(0,1.6)) +
    theme_classic() +
    labs(title = "True density") +
    theme(axis.text.y = element_blank()
          ,
         # legend.position = "none"
          ) 
  ggplot() + 
    geom_tile(data = D_plotdatlin[D_plotdatlin$quantile == "mean",c(1,2,5)],
              aes(x = x, y= y, fill = twoDdets)) +
    scale_fill_viridis_c(name = "D", limits = c(0, 1.6)) +
    theme_classic() +
    labs(title = "2D") +
    theme(axis.text.y = element_blank(),
          #legend.position = "none"
          ) 
  ggplot() + 
    geom_tile(data = D_plotdatlin[D_plotdatlin$quantile == "mean",c(1,2,6)],
              aes(x = x, y= y, fill = lindets)) +
    scale_fill_viridis_c(name = "D", limits = c(0, 1.6)) +
    theme_classic() +
    labs(title = "1D") +
    theme(axis.text.y = element_blank(),
          #legend.position = "none"
          ) 
  
  ggplot() + 
    geom_tile(data = D_plotdat[D_plotdat$quantile == "mean",c(1,2,3)],
              aes(x = x, y= y, fill = trueD)) +
    scale_fill_viridis_c(name = "D", limits = c(0,1.2)) +
    theme_classic() +
    labs(title = "True density") +
    theme(axis.text.y = element_blank()
          ,
          # legend.position = "none"
    ) 
  ggplot() + 
    geom_tile(data = D_plotdat[D_plotdat$quantile == "mean",c(1,2,5)],
              aes(x = x, y= y, fill = twoDdets)) +
    scale_fill_viridis_c(name = "D", limits = c(0,1.2)) +
    theme_classic() +
    labs(title = "2D") +
    theme(axis.text.y = element_blank()
          ,
          # legend.position = "none"
    ) 
  ggplot() + 
    geom_tile(data = D_plotdat[D_plotdat$quantile == "mean",c(1,2,6)],
              aes(x = x, y= y, fill = lindets)) +
    scale_fill_viridis_c(name = "D", limits = c(0,1.2)) +
    theme_classic() +
    labs(title = "1D") +
    theme(axis.text.y = element_blank()
          ,
          # legend.position = "none"
    ) 
  
  
  D_plotdatlong <- tidyr::pivot_longer(D_plotdat, cols = c("trueD", "trueDlin", "twoDdets", "lindets"))
  D_plotdatlong$true <- D_plotdatlong$name
  D_plotdatlong$true[which(D_plotdatlong$true %in% c("trueD", "trueDlin"))] <- "true"
   
  Dplot <- 
    ggplot() + 
    #add color outline to true 1d and 2d lines
    geom_line(data = D_plotdatlong[which(D_plotdatlong$quantile == "mean" & D_plotdatlong$name %in% c("trueD")),], 
              mapping = aes(x = x, y = value, group = name), color = plotcols[2],
              size = (linesize * 2)) +
    geom_line(data = D_plotdatlong[which(D_plotdatlong$quantile == "mean" & D_plotdatlong$name %in% c("trueDlin")),], 
              mapping = aes(x = x, y = value, group = name), color = plotcols[1],
              size = (linesize * 2)) +
    #now actually do colored lines
    geom_line(data = D_plotdatlong[D_plotdatlong$quantile == "mean",], 
              mapping = aes(x = x, y = value, 
                            col = as.factor(true),
                            group = name),
              size = linesize) +
    geom_line(data = D_plotdatlong[which(D_plotdatlong$quantile == "2.5%" & D_plotdatlong$true != "true") ,], 
              mapping = aes(x = x, y = value, 
                            col = as.factor(true),
                            group = name), linetype = "dashed", size = linesize) +
    geom_line(data = D_plotdatlong[which(D_plotdatlong$quantile == "97.5%" & D_plotdatlong$true != "true"),], 
              mapping = aes(x = x, y = value, 
                            col = as.factor(true),
                            group = name), linetype = "dashed", size = linesize) +

    scale_color_manual(values = plotcols[c(1,3,2)], labels = c("1D", "True", "2D"),
                       name = "") +
    #ylim(0,10) +
    xlim(15000,30000)+
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
  if(Dmodel == "flat"){
    out =   grid.arrange(
      grobs = list(lambda0plot, sigmaplot, Dplot, timeplot),
      widths = c(1,1),
      heights = c(1,1,1),
      layout_matrix = rbind(c(1,2), 
                            3,
                            4)
    )
  } else if( Dmodel == "variable"){
    out =   grid.arrange(
      grobs = list(lambda0plot, sigmaplot, beta1plot, beta2plot, beta3plot, timeplot, Dplot),
      widths = c(1,1,1,1,1,1),
      heights = c(1,
                  1,
                  1,
                  1.2),
      layout_matrix = rbind(c(1,1,1,2,2,2),
                            c(3,3,3,4,4,4),
                            c(5,5,5,6,6,6),
                            c(7,7,7,7,7,7))
    )
  }
return(out)
}

fitsflat2D <- readRDS(file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/simulation_results/2D/all_sim_fits1.Rds")
fitsflat1D <- readRDS(file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/simulation_results/1D/all_sim_fits1.Rds")

create_plots(sim_fits_out = fitsflat2D, mesh = mesh2D, D_mesh = D_mesh2D,
             sim_fits_outlin= fitsflat1D, meshlin = meshlin,D_meshlin = D_meshlin, 
             Dmodel = "flat", 
             plotcols = c("aquamarine4","#F1605DFF","black" ),
             linesize = 1, output = "plots")
ggsave(file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/plots/flatDcomparison.png",
       plot = create_plots(sim_fits_out = fitsflat2D, mesh = mesh2D, D_mesh = D_mesh2D,
                           sim_fits_outlin= fitsflat1D, meshlin = meshlin,D_meshlin = D_meshlin, 
                           Dmodel = "flat", 
                           plotcols = c("aquamarine4","#F1605DFF","black" ),
                           linesize = 1, output = "plots"),
       width = 169,
       height = 169,
       units = c("mm"),
       dpi = 300)

fitsvar2D <- readRDS(file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/simulation_results/2D/all_sim_fitsq.Rds")
fitsvar1D <- readRDS(file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/simulation_results/1D/all_sim_fitsq.Rds")

create_plots(sim_fits_out = fitsvar2D, mesh = mesh2D, D_mesh = D_mesh2D_q,
             sim_fits_outlin = fitsvar1D, meshlin = meshlin,D_meshlin = D_meshlin_q, 
             Dmodel = "variable", 
             plotcols = c("aquamarine4","#F1605DFF","black" ),
             linesize = 1, output = "plots")
ggsave(file = "~/Documents/UniStAndrews/MovingDetector/compare_moving_stat_1D/plots/varDcomparison.png",
       plot = create_plots(sim_fits_out = fitsvar2D, mesh = mesh2D, D_mesh = D_mesh2D_q,
                           sim_fits_outlin = fitsvar1D, meshlin = meshlin,D_meshlin = D_meshlin_q, 
                           Dmodel = "variable", 
                           plotcols = c("aquamarine4","#F1605DFF","black" ),
                           linesize = 1, output = "plots"),
       width = 169,
       height = 169*4/3,
       units = c("mm"),
       dpi = 300)

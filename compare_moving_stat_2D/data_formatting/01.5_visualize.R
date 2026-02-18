
setup <- setup_data(
  sigma,
  N,
  beta1,
  beta2,
  ntrapsish,
  trackxmin,
  trapspacing,
  meshspacing,
  trap_n_horiz,
  nsteps_pertrap,
  occreps)

traps <- setup$traps
tracksdf <- setup$tracksdf
mesh <- setup$mesh
meshspacing <- setup$meshspacing
dist_trapmesh <- setup$dist_trapmesh
useall <- setup$useall
D_mesh_f <- setup$D_mesh_f
D_mesh_v <- setup$D_mesh_v

rm(setup)

################################################################################
#### ----------------------------visualize -----------------------------------
################################################################################
# pick which mesh to check

if(printplots){
  testmesh <- D_mesh_v
  
  layoutplot <- ggplot() +
    geom_raster(data.frame(x = mesh$x, y = mesh$y, D = testmesh),
                mapping = aes(x = x, y = y, fill = D)) +
    geom_point(data = tracksdf, mapping = aes(x = x, y = y, group = occ),
               size = 3,shape = "+", color= "white"
    ) +
    scale_fill_viridis_c()
  seetrap_plot <- layoutplot +
    #xlim(1000, 4000) +
    geom_point(data.frame(x = traps$x,
                          y = traps$y),
               mapping = aes(x = x, y = y),
               shape = 21, color = "pink",
               fill = "red",
               alpha = .7, size = 2) +
    scale_y_continuous(expand = c(0,0)) +
    scale_x_continuous(expand = c(0,0)) +
    theme_bw() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          legend.position = "none",
          text = element_text(size = fontsize))
  
  #getavgx <- function(x){
  # visualize a capture history
  testpop <- sim_pop_C(testmesh,
                       as.matrix(mesh),
                       meshspacing)
  
  testcapthist_full <- sim_capthist_C(as.matrix(traps),
                                      tracksdf,
                                      lambda0,
                                      sigma,
                                      testmesh,
                                      as.matrix(mesh),
                                      meshspacing,
                                      hazdenom,
                                      testpop,
                                      dist_dat_pop = NULL,
                                      report_hus = F,
                                      report_probseenxk = F)
  #average x value of trap of first det
  #sum((traps$x * apply(testcapthist_full, 3, function(j){sum(!is.na(j))})))/sum(!is.na(testcapthist_full))
  #}
  
  #mean(apply(as.array(1:500),1,getavgx))
  
  tocck = sample(nocc,1)
  tindi = sample(nrow(testpop),1)
  
  if(length(dim(testcapthist_full)) == 3){
    
    testcapthist <- testcapthist_full[which(apply((!is.na(testcapthist_full)), 1, sum)>0),,]
    testinduse <- create_ind_use_C(testcapthist_full,
                                   as.matrix(traps),
                                   trapspacing,
                                   tracksdf,
                                   scenario = "everything")
    
    example_ch_plot <- ggplot() +
      scale_fill_viridis_c(guide = "none")+
      geom_point(data.frame(x = traps$x,
                            y = traps$y,
                            induse1 = testinduse[tindi,,tocck]),
                 mapping = aes(x = x ,y = y, fill = induse1),
                 shape = 21, size = 2) +
      new_scale_fill() +
      geom_point(data.frame(x = testpop[,1],
                            y = testpop[,2],
                            detocci = apply(testcapthist_full[,tocck,], 1, sum, na.rm = T)>0,
                            interest = seq(1:nrow(testpop)) == tindi
      ), mapping = aes(x = x, y = y, fill = detocci, color = interest), shape = 21, stroke = 1.5) +
      scale_color_manual(values = c("transparent", "red"),
                         labels = NULL,
                         name = "Focus\nindividual\nand trap") +
      scale_fill_discrete(name = "Detected\nthis survey?") +
      geom_point(data.frame(x = traps$x[which(!is.na(testcapthist_full[tindi,tocck,]))],
                            y = traps$y[which(!is.na(testcapthist_full[tindi,tocck,]))]
      ), mapping = aes(x = x, y = y), color = "red", shape = 1,
      stroke = 1.5) +
      guides(
        fill = guide_legend(
          order = 1,
          override.aes = list(
            shape = 21,
            colour = "transparent",
            stroke = 1.5
          )
        ),
        color = "none"#guide_legend(
        #order = 2,
        #override.aes = list(
        #  shape = 21,
        # fill = "white",
        #  stroke = 1.5
        #)
        # )
      ) +
      theme_bw() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(size = fontsize))
    
    detmax <- max(c(apply(testcapthist_full, 1, function(i){sum(!is.na(i))}),
                    apply(testcapthist_full, 3, function(j){sum(!is.na(j))})))
    
    popdf <- data.frame(x = testpop[,1],
                        y = testpop[,2],
                        totdets = apply(testcapthist_full, 1, function(i){sum(!is.na(i))})
    )
    
    popdet_plot <- layoutplot +
      geom_point(popdf, mapping = aes(x = x, y = y,
                                      color = totdets),
                 size = 3) +
      scale_color_viridis_c(option = "magma", name = "dets",
                            limits = c(0, detmax)) +
      guides(fill = "none") +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      theme_bw() +
      theme(axis.text.y = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            legend.position = "none",
            text = element_text(size = fontsize))
    
    trapdet_plot <- layoutplot +
      geom_point(data.frame(x = traps$x,
                            y = traps$y,
                            dets = apply(testcapthist_full, 3, function(j){sum(!is.na(j))})),
                 mapping = aes(x = x, y = y, color = dets),
                 size = 3) +
      geom_vline(xintercept = -beta2 - 3*sigma, color = "white", linetype = "dashed") +
      scale_color_viridis_c(option = "magma", name = "Detections"#,
                            # limits = c(0, detmax)
      ) +
      guides(fill = "none") +
      scale_y_continuous(expand = c(0,0)) +
      scale_x_continuous(expand = c(0,0)) +
      theme_bw() +
      theme(axis.text.y = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(size = fontsize))
    
    grid.arrange(seetrap_plot,
                 example_ch_plot,
                 popdet_plot,
                 trapdet_plot)
    dim(testcapthist)
    summary(apply(testcapthist_full, 1, function(x){sum(!is.na(x))}))     # detections per individual
    
    
    
  } else if(length(dim(testcapthist_full)) == 2){ #just testing
    tracksdfmidpt <- calc_trackmidpts(tracksdf)
    middistsi <- calc_dist_matC(as.matrix(testpop[c(tindi, tindi+1),]),
                                tracksdfmidpt)[1,]
    distdf <- data.frame(x = tracksdf$x[tracksdf$occ == tocck],
                         y = tracksdf$y[tracksdf$occ == tocck],
                         dist = middistsi[tracksdf$occ == tocck],
                         hu = testcapthist_full[tindi, tracksdf$occ == tocck])
    
    ggplot() +
      scale_fill_viridis_c(guide = "none")+
      # geom_point(distdf,
      #           mapping = aes(x = x,
      #                         y = y,
      #                         fill = dist),
      #           size = 8, shape = 21) +
      geom_point(distdf,
                 mapping = aes(x = x ,
                               y = y,
                               fill = hu),
                 shape = 21, size = 4,
                 color = "white") +
      geom_point(data.frame(x = testpop[tindi,1],
                            y = testpop[tindi,2]
      ),
      mapping = aes(x = x,
                    y = y),
      shape = 21, stroke = 1.5) +
      
      theme_bw() +
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            axis.ticks = element_blank(),
            text = element_text(size = fontsize))
    
    
    
  }
  
}



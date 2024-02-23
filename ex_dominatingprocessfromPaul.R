simSTCR <- function(N = 50, sigma = 0.5, lambda = 0.5, beta = 1, 
                    StudyPeriod = 25, traps, xlim, ylim, keep0 = FALSE)
{
  J <- nrow(traps)
  
  detectFunc <- function(tdiff, x0, x, sigma, beta, s, lambda)
  {
    d2 <- (x[1] - (s[1] + (x0[1] - s[1])*exp(-beta*tdiff)))^2 + (x[2] - (s[2] + (x0[2] - s[2])*exp(-beta*tdiff)))^2
    g <- exp(-d2/(2*(1-exp(-2*beta*tdiff))*sigma^2))
    return(g)
  }
  
  
  S <- cbind(x = runif(N, xlim[1], xlim[2]), 
             y = runif(N, ylim[1], ylim[2]))
  capt.hist <- data.frame()
  ID <- 0
  capti <- data.frame()
  for(i in 1:N)
  {
    sk <- S[i,]
    mu0 <- rnorm(2, mean = sk, sd = sigma)
    mu0_s <- mu0
    tLast <- 0
    captured <- 0
    ## Animal attempts to teleport J*lambda*StudyPeriod times
    times <- sort(runif(rpois(1, J*lambda*StudyPeriod), 0, StudyPeriod))
    ## To detectors omegas.
    omegas <- sample(J, length(times), replace = TRUE)
    ## Now thin those jump attempts spatially and temporally to accept them according to STCR.
    for(i in 1:length(times))
    {
      ti <- times[i]
      pkj <- detectFunc(tdiff = ti - tLast, x0 = mu0, x = as.numeric(traps[omegas[i],]), 
                        sigma = sigma, beta = beta, s = sk, lambda = lambda)				
      
      ## If we accept jump, then the temporal correlation is reset at 
      ## ti, mu0 = trap[wi,]
      if(runif(1, 0, 1) < pkj){
        wi <- omegas[i]
        mu0 <- as.numeric(traps[wi,])
        capti <- rbind(capti, data.frame(site = wi, times = ti, ID = ID+1, x01 = mu0_s[1], x02 = mu0_s[2], sk1 = sk[1], sk2 = sk[2]))
        captured <- 1
        tLast <- ti
      }
    }
    ID <- ID + captured
    if(keep0 == TRUE & captured == 0)
    {
      capti <- rbind(capti, data.frame(site = wi, times = StudyPeriod + 10, ID = -1, x01 = mu0_s[1], x02 = mu0_s[2], sk1 = sk[1], sk2 = sk[2]))
    }
    # print(i)
  }
  capti
}

## Lewis, Peter A., and Gerald S. Shedler... by Van Dam-Bates, Paul (he, him / il, lui) (DFO/MPO)
Van Dam-Bates, Paul (he, him / il, lui) (DFO/MPO)
3:35 PM

## Lewis, Peter A., and Gerald S. Shedler. 
## "Simulation of nonhomogeneous Poisson processes by thinning." 
## Naval Research Logistics (NRL) 26.3 (1979): 403-413
simSTCRThinning <- function(N = 50, sigma = 0.5, lambda = 0.5, beta = 1, 
                            StudyPeriod = 25, traps, xlim, ylim, keep0 = FALSE)
{
  J <- nrow(traps)
  lamstar <- J*lambda + 0.1
  
  hazardFunc <- function(tdiff, x0, x, sigma, beta, s, lambda)
  {
    d2 <- (x[,1] - (s[1] + (x0[1] - s[1])*exp(-beta*tdiff)))^2 + (x[,2] - (s[2] + (x0[2] - s[2])*exp(-beta*tdiff)))^2
    g <- exp(-d2/(2*(1-exp(-2*beta*tdiff))*sigma^2))
    return(lambda*g)
  }
  
  
  S <- cbind(x = runif(N, xlim[1], xlim[2]), 
             y = runif(N, ylim[1], ylim[2]))
  capt.hist <- data.frame()
  ID <- 0
  capti <- data.frame()
  for(i in 1:N)
  {
    sk <- S[i,]
    mu0 <- rnorm(2, mean = sk, sd = sigma)
    mu0_s <- mu0
    ti <- 0
    tLast <- 0
    captured <- 0
    while(ti < StudyPeriod)
    {
      u1 <- runif(1)
      t1 <- - 1/lamstar * log(u1)
      ti <- ti + t1
      
      u2 <- runif(1)
      hkj <- hazardFunc(tdiff = ti - tLast, x0 = mu0, x = as.matrix(traps), 
                        sigma = sigma, beta = beta, s = sk, lambda = lambda)				
      hk <- sum(hkj)
      
      if(u2 <= hk/lamstar){
        wi <- sample(1:J, 1, prob = hkj)
        mu0 <- as.numeric(traps[wi,])
        capti <- rbind(capti, data.frame(site = wi, times = ti, ID = ID+1, x01 = mu0_s[1], x02 = mu0_s[2], sk1 = sk[1], sk2 = sk[2]))
        captured <- 1
        tLast <- ti
      }
    }
    ID <- ID + captured
    if(keep0 == TRUE & captured == 0)
    {
      capti <- rbind(capti, data.frame(site = wi, times = StudyPeriod + 10, ID = -1, x01 = mu0_s[1], x02 = mu0_s[2], sk1 = sk[1], sk2 = sk[2]))
    }
    # print(i)
  }
  capti
}

has context menu
## -----------------------------------------------------------------------------
##
## Script name: trial.R
##
## Purpose of script: main trial simulation loop.
##
## Author: MAJ
##
## Date Created: 2021-03-01
##
## -----------------------------------------------------------------------------
##
## Each simulated trial calls `trial1` once to simulate data and complete all
## interim analyses and the final analysis and then returns the dataset and
## data to post process trial decisions.
##
## -----------------------------------------------------------------------------

# library(colormap)
# library(matrixStats)
# library(parallel)
# library(emmeans)
# library(randomizr)
# library(data.table)
# library(rjags)
# library(runjags)

# remove when done:
# source("./R/zzz.R")
# globals$glab_brand <- c("AZ", "PF")
# globals$glab_trt <- rbind(c("CVD+FLU", "CVD+PBO"),
#                            c("PBO", "FLU"))
# source("./R/utils.R")
# source("./R/data.R")
# library(coadministr.stanc)
# looks = c(200, 300, 400, 500)
# n_per_yr = 2000
# # ni margin (on probability scale)
# delta_ni = 0.2
# beta = pr_adverse_evt(pae_coad = 0.1)
# sd_id = 0.5
# p_miss = 0.0
# Ncoad = 2
# Nbrand = 2
# Ncohort = 2
# a0 = 1; b0 = 1; eta0 = 0.01 


trial <- function(
  scenario = 0,
  looks = c(200, 300, 400, 500),
  n_per_yr = 2000,
  # ni margin (on probability scale)
  delta_ni = 0.2,
  beta = NULL,
  sd_id = 0.5,
  p_miss = 0.0,
  Ncoad = 2, #  flu vs pbo - estimate indep
  Nbrand = 2, # az, pf
  Ncohort = 2, # cvd dose 1, cvd dose 2
  a0 = 1, # hyperpar on beta_pdf for mu
  b0 = 1, # hyperpar on beta_pdf for mu
  eta0 = 0.01, # hyperpar on eta (a + b)
  keep_post_sampls = FALSE
){
  
  # Bundling up config
  lpar <- list(
    looks = looks,
    N_max = max(looks),
    n_per_yr = n_per_yr,
    delta_ni = delta_ni,
    beta = beta, sd_id = sd_id, p_miss = p_miss,
    Ncoad = Ncoad, Nbrand = Nbrand, Ncohort = Ncohort,
    a0 = a0, b0 = b0, eta0 = eta0
  )
  
  # logging
  set_hash()
  time0 <- Sys.time()
  message(get_hash(), " Started (scenario ", scenario, ")")
  
  # Data for this trial ----------------------------------
  
  # Ref grid
  drefgrid <- data.table(
    arm = rep(globals$glab_trt[1, ], each = 4),
    brand = rep(rep(globals$glab_brand, each = 2), len = 8),
    cohort = rep(1:2, len = 8)
  )
  
  # Set indexes, get dataset for this trial
  curr_look <- 1
  ixs <- ixe <- 0
  
  d <- get_data(J = lpar$N_max,
                beta = lpar$beta,
                sd_id = lpar$sd_id,
                n_per_yr = lpar$n_per_yr,
                p_miss = lpar$p_miss)
  
  # Results container
  lr <- lapply(1:length(lpar$looks), function(i) {
    ll <- list()
    ll$look      = i
    ll
  })
  # end Data for this trial -----------------------------
  
  lpost <- NULL
  if(keep_post_sampls){
    lpost <- list(); ip <- 1
  }

  while(lpar$looks[curr_look] < lpar$N_max){
    
    # Setup for this interim data range indexes for this interim ---------------
    l1 <- NULL; fitcount <- 0
    ixs <- ixe + 1
    ixe <- lpar$looks[curr_look]
    # current time is time of enrollment of last participant in this interim
    curr_time <- d[id == ixe & timept == 1, visit_time]
    
    if(is.na(curr_time)){
      message("Current time is NA")
    }
    # end Setup for this interim data range indexes for this interim ---------------
    
    # Data in binomial format (events out of trials) --------------
    dtbl <- d[id <= ixe & futime <= curr_time & timept == 1 & miss == 0,
              .(y=sum(y), n=.N, observed = mean(y)),
              keyby = .(arm, brand, cohort)]
    dtbl[, arm := globals$glab_trt[1, arm]]
    dtbl[, brand := globals$glab_brand[brand]]
    dtbl <- merge(drefgrid, dtbl, by = c("arm", "brand", "cohort"),
                  all.x = TRUE)[order(arm, brand, cohort)]
    dtbl[is.na(y), `:=`(y=0, n=0, observed=0)]
    
    # form of matrix is:
    # AZ1 PF1 
    # AZ2 PF2
    yflu <- array(dtbl[arm == "CVD+FLU", y], dim = c(2, 2))
    ypbo <- array(dtbl[arm == "CVD+PBO", y], dim = c(2, 2))
    nflu <- array(dtbl[arm == "CVD+FLU", n], dim = c(2, 2))
    npbo <- array(dtbl[arm == "CVD+PBO", n], dim = c(2, 2))
    
    ld = list(yflu = yflu, ypbo = ypbo, nflu = nflu, npbo = npbo,
              Ncoad = lpar$Ncoad, Nbrand = lpar$Nbrand, Ncohort = lpar$Ncohort,
              a0 = lpar$a0, b0 = lpar$b0, eta0 = lpar$eta0
    )
    # end Data in binomial --------------------------------------------
    
    # Fit model ---------------------------------------
    
    # pars are indexed cohort rows x brand cols
    while(is.null(l1) & fitcount < 10){
      l1 <- tryCatch({
        rstan::sampling(coadministr.stanc::stanc_betabinhier(), data = ld,
                        chains = 1, thin = 1,
                        iter = 2000, warmup = 1000, refresh = 0)
      }, error=function(e) { message(" Error " , e); return(NULL) })
      fitcount <- fitcount + 1
    }
    
    if(keep_post_sampls){
      lpost[[ip]] <- rstan::extract(l1)
      ip <- ip + 1
    }
    
    m <- rstan::extract(l1, pars = c("pflu", "ppbo"))
    post <- data.table(cbind(m$pflu[, , 1], m$pflu[, , 2],
                             m$ppbo[, , 1], m$ppbo[, , 2])
                       )
    colnames(post) <- c("p_cvdflu_az_1", "p_cvdflu_az_2",
                        "p_cvdflu_pf_1", "p_cvdflu_pf_2",
                        "p_cvdpbo_az_1", "p_cvdpbo_az_2",
                        "p_cvdpbo_pf_1", "p_cvdpbo_pf_2")
    # bind together az cohort and brand for flu vs pbo
    post1 <- post[, .(p_cvdflu_az = c(p_cvdflu_az_1, p_cvdflu_az_2),
                      p_cvdpbo_az = c(p_cvdpbo_az_1, p_cvdpbo_az_2))]
    post2 <- post[, .(p_cvdflu_pf = c(p_cvdflu_pf_1, p_cvdflu_pf_2),
                      p_cvdpbo_pf = c(p_cvdpbo_pf_1, p_cvdpbo_pf_2))]
    # end Fit model ---------------------------------------
    
    
    # Results -------------------------------------------------------
    
    if(nrow(dtbl) != 8) stop()
    
    # data available at this interim
    lr[[curr_look]]$dtbl      <- dtbl
    # posterior means cohort and brand
    lr[[curr_look]]$p         <- colMeans(post)
    # averages out the cohort variability
    lr[[curr_look]]$p_az      <- colMeans(post1)
    lr[[curr_look]]$p_pf      <- colMeans(post2)
    
    # differences in prob of ae for cmb vs mono across brand and cohort
    lr[[curr_look]]$delta <- c(mean(post[[1]] - post[[5]]),
                               mean(post[[2]] - post[[6]]),
                               mean(post[[3]] - post[[7]]),
                               mean(post[[4]] - post[[8]]))
    # difference in prob of ae for cmb vs mono across brand
    lr[[curr_look]]$delta_az <- mean(post1[[1]] - post1[[2]])
    lr[[curr_look]]$delta_pf <- mean(post2[[1]] - post2[[2]])
    
    # pr ae rate in cmb is higher than mono across brand and cohort
    lr[[curr_look]]$pr_inf <- c(mean(post[[1]] - post[[5]] > 0),
                                mean(post[[2]] - post[[6]] > 0),
                                mean(post[[3]] - post[[7]] > 0),
                                mean(post[[4]] - post[[8]] > 0))
    
    # pr ae rate in cmb is higher than mono across brand
    lr[[curr_look]]$pr_inf_az <- mean(post1[[1]] - post1[[2]] > 0)
    lr[[curr_look]]$pr_inf_pf <- mean(post2[[1]] - post2[[2]] > 0)
    
    # pr that ae rate in cmb is no more than delta_ni higher than to mono
    # across brand and cohort
    lr[[curr_look]]$pr_ninf   <- c(mean(post[[1]] - post[[5]] < lpar$delta_ni),
                                   mean(post[[2]] - post[[6]] < lpar$delta_ni),
                                   mean(post[[3]] - post[[7]] < lpar$delta_ni),
                                   mean(post[[4]] - post[[8]] < lpar$delta_ni))
    
    lr[[curr_look]]$pr_ninf_az <- mean(post1[[1]] - post1[[2]] < lpar$delta_ni)
    lr[[curr_look]]$pr_ninf_pf <- mean(post2[[1]] - post2[[2]] < lpar$delta_ni)
    
    lr[[curr_look]]$curr_time <- curr_time
    # end Results -------------------------------------------------------
    
    
    # Update trackers -----------------------------------------
    ixs <- ixe
    curr_look <- curr_look + 1
    # end Update trackers -----------------------------------------
    
    
  }# look
  
  
  # Final analysis when all participants reach end of followup
  
  
  # Data in binomial format, events out of trials --------------
  dtbl <- d[miss == 0,
            .(y=sum(y), n=.N, observed = mean(y)),
            keyby = .(arm, brand, cohort)]
  dtbl[, arm := globals$glab_trt[1, arm]]
  dtbl[, brand := globals$glab_brand[brand]]
  dtbl <- merge(drefgrid, dtbl, by = c("arm", "brand", "cohort"),
                all.x = TRUE)[order(arm, brand, cohort)]
  dtbl[is.na(y), `:=`(y=0, n=0, observed=0)]
  # put into a list for stan
  yflu <- array(dtbl[arm == "CVD+FLU", y], dim = c(2, 2))
  ypbo <- array(dtbl[arm == "CVD+PBO", y], dim = c(2, 2))
  nflu <- array(dtbl[arm == "CVD+FLU", n], dim = c(2, 2))
  npbo <- array(dtbl[arm == "CVD+PBO", n], dim = c(2, 2))
  
  ld = list(yflu = yflu, ypbo = ypbo, nflu = nflu, npbo = npbo,
            Ncoad = lpar$Ncoad, Nbrand = lpar$Nbrand, Ncohort = lpar$Ncohort,
            a0 = lpar$a0, b0 = lpar$b0, eta0 = lpar$eta0
  )
  # end Data in binomial --------------------------------------------
  
  
  # Fit model ---------------------------------------
  
  # pars are indexed cohort rows x brand cols
  l1 <- NULL; fitcount <- 0
  while(is.null(l1) & fitcount < 10){
    l1 <- tryCatch({
      rstan::sampling(coadministr.stanc::stanc_betabinhier(), data = ld,
                      chains = 1, thin = 1,
                      iter = 2000, warmup = 1000, refresh = 0)
    }, error=function(e) { message(" Error " , e); return(NULL) })
    fitcount <- fitcount + 1
  }
  
  if(keep_post_sampls){
    lpost[[ip]] <- rstan::extract(l1)
    ip <- ip + 1
  }
  
  m <- rstan::extract(l1, pars = c("pflu", "ppbo"))
  post <- data.table(cbind(m$pflu[, , 1], m$pflu[, , 2],
                           m$ppbo[, , 1], m$ppbo[, , 2])
  )
  colnames(post) <- c("p_cvdflu_az_1", "p_cvdflu_az_2",
                      "p_cvdflu_pf_1", "p_cvdflu_pf_2",
                      "p_cvdpbo_az_1", "p_cvdpbo_az_2",
                      "p_cvdpbo_pf_1", "p_cvdpbo_pf_2")
  # bind together az cohort and brand for flu vs pbo
  post1 <- post[, .(p_cvdflu_az = c(p_cvdflu_az_1, p_cvdflu_az_2),
                    p_cvdpbo_az = c(p_cvdpbo_az_1, p_cvdpbo_az_2))]
  post2 <- post[, .(p_cvdflu_pf = c(p_cvdflu_pf_1, p_cvdflu_pf_2),
                    p_cvdpbo_pf = c(p_cvdpbo_pf_1, p_cvdpbo_pf_2))]
  # end Fit model ---------------------------------------
  
  
  # Results -------------------------------------------------------
  
  if(nrow(dtbl) != 8) stop()
  
  # data available at this interim
  lr[[curr_look]]$dtbl      <- dtbl
  # posterior means cohort and brand
  lr[[curr_look]]$p         <- colMeans(post)
  # averages out the cohort variability
  lr[[curr_look]]$p_az      <- colMeans(post1)
  lr[[curr_look]]$p_pf      <- colMeans(post2)
  
  # differences in prob of ae for cmb vs mono across brand and cohort
  lr[[curr_look]]$delta <- c(mean(post[[1]] - post[[5]]),
                             mean(post[[2]] - post[[6]]),
                             mean(post[[3]] - post[[7]]),
                             mean(post[[4]] - post[[8]]))
  # difference in prob of ae for cmb vs mono across brand
  lr[[curr_look]]$delta_az <- mean(post1[[1]] - post1[[2]])
  lr[[curr_look]]$delta_pf <- mean(post2[[1]] - post2[[2]])
  
  # pr ae rate in cmb is higher than mono across brand and cohort
  lr[[curr_look]]$pr_inf <- c(mean(post[[1]] - post[[5]] > 0),
                              mean(post[[2]] - post[[6]] > 0),
                              mean(post[[3]] - post[[7]] > 0),
                              mean(post[[4]] - post[[8]] > 0))
  
  # pr ae rate in cmb is higher than mono across brand
  lr[[curr_look]]$pr_inf_az <- mean(post1[[1]] - post1[[2]] > 0)
  lr[[curr_look]]$pr_inf_pf <- mean(post2[[1]] - post2[[2]] > 0)
  
  # pr that ae rate in cmb is no more than delta_ni higher than to mono
  # across brand and cohort
  lr[[curr_look]]$pr_ninf   <- c(mean(post[[1]] - post[[5]] < lpar$delta_ni),
                                 mean(post[[2]] - post[[6]] < lpar$delta_ni),
                                 mean(post[[3]] - post[[7]] < lpar$delta_ni),
                                 mean(post[[4]] - post[[8]] < lpar$delta_ni))
  
  lr[[curr_look]]$pr_ninf_az <- mean(post1[[1]] - post1[[2]] < lpar$delta_ni)
  lr[[curr_look]]$pr_ninf_pf <- mean(post2[[1]] - post2[[2]] < lpar$delta_ni)
  
  lr[[curr_look]]$curr_time <- curr_time
  # end Results -------------------------------------------------------

  duration <- Sys.time() - time0
  message(get_hash(), " Completed, duration ", 
          round(as.numeric(duration, units = "secs"), 1), " secs")
  
  list(d = d,
       lr = lr,
       lpar = lpar,
       lpost = lpost)

}



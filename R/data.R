## -----------------------------------------------------------------------------
##
## Script name: data.R
##
## Purpose of script: creates dataset for trial
##
## Author: MAJ
##
## Date Created: 2021-03-01
##
## -----------------------------------------------------------------------------
##
##
## -----------------------------------------------------------------------------




pr_adverse_evt <- function(pae_pbo = 0.02, pae_flu = 0.05, 
                           pae_coad_az = 0, pae_coad_pf = 0,
                           pae_az1 = 0.5, pae_az2 = 0.4,
                           pae_pf1 = 0.3, pae_pf2 = 0.5
){
  
  # Probability of AE is dependent on covid dose (1 or 2) and brand.
  # Participants are enrolled either on covid dose 1 or dose 2. If they were
  # enrolled at dose 1, they are inelligible for randomisation at dose 2.
  # AZ has prob of AE 0.6, 0.4 at dose 1 and dose 2
  # PF has prob of AE 0.3, 0.6 at dose 1 and dose 2
  
  m <- expand.grid(brand = 1:2, cohort = 1:2, timept = 1:2, arm = 1:2)
  pr_ae <- data.table(m)[order(brand, cohort)]
  pr_ae[, pid := 1:nrow(m)]
  # add labels
  pr_ae[, lab := sapply(1:nrow(pr_ae), function(i){
    paste0(globals$glab_trt[pr_ae$timept[i], pr_ae$arm[i]], " (", globals$glab_brand[pr_ae$brand[i]], ")")
  }) ]
  # common pr of AE
  pr_ae[, pae0 := c(pae_az1, pae_pbo, pae_az1, pae_flu, # AZ @ dose 1 high reactivity
                    pae_az2, pae_pbo, pae_az2, pae_flu, # AZ @ dose 2 reactivity drops
                    pae_pf1, pae_pbo, pae_pf1, pae_flu, # PF @ dose 1 low reactivity
                    pae_pf2, pae_pbo, pae_pf2, pae_flu  # AZ @ dose 2 reactivity increases
  )]
  
  # brand level effects of coadministration with flu
  pr_ae[, es := 0]
  pr_ae[, pae := pae0]
  pr_ae[timept == 1 & arm == 1 & brand == 1, es := pae_coad_az]
  pr_ae[timept == 1 & arm == 1 & brand == 2, es := pae_coad_pf]
  pr_ae[timept == 1 & arm == 1, pae := pae0 + es]
  pr_ae[, lo := qlogis(pae)]
  
  pr_ae
}



get_data <- function(J = 100,
                    beta = NULL,
                    sd_id = 0.5,
                    n_per_yr = 100,
                    p_miss = 0){
  
  stopifnot(!is.null(beta))
  
  # cohort is covid dose 1 or dose 2
  n_cohorts <- 2
  # timepoint 1 gets cvd+flu or cvd+pbo
  # timepoint 2 gets pbo or flu
  n_timepoints <- 2
  n_arm <- 2
  n_brand <- length(globals$glab_brand)
  
  # participant only ever in cohort 1 OR cohort 2 no overlap
  id <- rep(1:J, each = n_timepoints)
  
  time1 <- cumsum(c(0, rexp(J-1, n_per_yr/365))) # co-administered shot
  time2 <- time1 + runif(J, 8, 14)               # single flu or pbo shot
  #
  fu1time <- time1 + runif(J, 6.5, 7.5) # self report on coad
  fu2time <- time2 + runif(J, 6.5, 7.5)  # self report on single shot
  
  # visit time is the time participant gets vaccination
  visit_time <- numeric(length(id))
  idx1 <- seq(1, length(id), by = 2)
  idx2 <- seq(2, length(id), by = 2)
  visit_time[idx1] <- time1
  visit_time[idx2] <- time2
  # fu time is the time at which we evalute the endpoint
  futime <- numeric(length(id))
  futime[idx1] <- fu1time
  futime[idx2] <- fu2time
  
  # 50% chance enrolled at dose 1 and 50% chance enrolled at dose 2
  pr_cohort <- 0.5
  cohort <- rbinom(J, 1, pr_cohort) + 1
  idx1 <- which(cohort == 1)
  idx2 <- which(cohort == 2)
  
  # 70% chance of being 1 = AZ
  # 30% chance of being 2 = PF
  pr_brand <- 0.3
  brand <- rbinom(J, 1, pr_brand) + 1
  
  # want 1:1 within cohort by brand
  block <- factor(paste0(cohort, "_", globals$glab_brand[brand]))
  
  # treatment allocation blocked on brand
  arm <- randomizr::block_ra(blocks = block,
                             conditions = 1:n_arm)
  # table(block, arm)
  block <- rep(block, each = n_timepoints)
  timept <- rep(1:n_timepoints, len = length(id))
  # Participants enrolled either at covid dose 1 or covid dose 2
  cohort <- rep(cohort, each = n_timepoints)
  # Participants will receive AZ or PF
  brand <- rep(brand, each = n_timepoints)
  # Treatment id either:
  # 1. CVD+FLU followed by PBO or
  # 2. CVD+PBO followed by FLU
  arm <- rep(arm, each = n_timepoints)
  
  lab <- sapply(seq_along(arm), function(i){
    paste0(globals$glab_trt[timept[i], arm[i]], " (", globals$glab_brand[brand[i]], ")")
  })
  
  d <- data.table(id, brand, cohort, block, timept, arm, visit_time, futime, lab)
  
  # pick up log odds of reaction
  d <- merge(d, beta[, .(brand, cohort, timept, arm, lo)],
             by = c("brand", "cohort", "timept", "arm"))
  setkey(d, id)
  
  # Within subject tendency for reacting to any vac
  lu_id <- rnorm(J, 0, sd_id)
  
  d[, u_id := lu_id[id]]
  d[, pae := plogis(lo + u_id)]
  d[, y := rbinom(.N, 1, prob = pae)]
  
  d[, miss := NA_integer_]
  d[, miss := rbinom(.N, 1, prob = p_miss)]
  
  d
}





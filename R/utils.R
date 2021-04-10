set_hash <- function(){
  globals$MYHASH <- digest::digest(as.character(runif(1)), algo = "xxhash32")
}
get_hash <- function(){
  globals$MYHASH
}

#' JAGS model string for trial inference
#'
#' Implementation of a hierarchical beta-binomial
#'
#' Implementation of JAGS code to independently model the adverse event rate
#' by treatment group (cvd+flu vs cvd+pbo) with partial pooling over
#' covid vaccine brand and cohort (first or second covid dose).
#' Participants are enrolled at either first or second covid dose and only
#' contribute in one or the other.
#'
#' @return string of JAGS model code
#' @export
#'
model_jags <- function(){
  
  # Partial pool across brand and cohort. Cohort membership
  # is determined by whether you were enrolled at cvd dose 1
  # or cvd dose 2. We are assuming that within a treatment
  # arm (CMB vs MONO), the probability of ae arises from some
  # common beta distribution.
  
  # https://stackoverflow.com/questions/47135726/error-slicer-stuck-at-value-with-infinite-density-running-binomial-beta-model
  # https://xianblog.wordpress.com/2011/07/28/11433/
  
  m1 <- "
model {
  ## likelihood
  for (i in 1:N){
     y_1[i] ~ dbin(p_cvdflu[i], n_1[i])
     y_2[i] ~ dbin(p_cvdpbo[i], n_2[i])
  }
  ## priors
  for (i in 1:N){
     # note the truncated distribution to prevent
     # Slicer stuck at value with infinite density
     p_cvdflu[i] ~ dbeta(a[1], b[1]) T(0.001,0.999)
     p_cvdpbo[i] ~ dbeta(a[2], b[2]) T(0.001,0.999)
  }
  ## hyperpriors for CMB vs MONO
  for(i in 1:2){
    mu[i] ~ dbeta(mua, mub) T(0.001,0.999)
    logeta[i] ~ dlogis(logn, 1)
    a[i] <- mu[i] * eta[i]
    b[i] <- (1-mu[i]) * eta[i]
    eta[i] <- exp(logeta[i])
  }}"
  return(m1)
}



post_to_data_table <- function(pp){
  
  nn <- names(pp)
  
  post <- do.call(cbind, lapply(1:length(pp), function(i){
    d <- switch(nn[i],
                pflu = data.table(cbind(pp$pflu[, , 1], pp$pflu[, , 2])),
                ppbo = data.table(cbind(pp$ppbo[, , 1], pp$ppbo[, , 2])),
                mu = data.table(pp$mu),
                eta = data.table(pp$eta),
                a = data.table(pp$a),
                b = data.table(pp$b),
                lp__ = data.table(pp$lp__)
    )
    names(d) <- paste0(nn[i], 1:ncol(d))
    d
  }))
  
  post
}
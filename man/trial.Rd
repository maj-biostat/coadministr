\name{trial}
\alias{trial}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Run a simulated trial}
\description{
Runs the all analyses for a coadministration trial
}
\usage{
trial <- function(
    scenario = 0,
    looks = c(200, 300, 400, 500), n_per_yr = 2000, delta_ni = 0.2, 
    beta = NULL, sd_id = 0.5, p_miss = 0.0,
    Ncoad = 2, Nbrand = 2, Ncohort = 2, 
    a0 = 1, b0 = 1, eta0 = 0.01, keep_post_sampls = FALSE
    )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{scenario}{an id for designating the current scenario}
  \item{looks}{sample size to run interim and final analyses}
  \item{n_per_yr}{enrollment per year}
  \item{delta_ni}{non-inferiority margin on probability scale}
  \item{beta}{simulation parameters, see \code{\link{pr_adverse_evt}}}
  \item{sd_id}{subject level variation}
  \item{p_miss}{probability of missingness}
  \item{Ncoad}{number of treatments - presently fixed to 2}
  \item{Nbrand}{number of brands - presently fixed to 2}
  \item{Ncohort}{number of cohorts - presently fixed to 2}
  \item{a0}{hyperparameter for beta prior}
  \item{b0}{hyperparameter for beta prior}
  \item{eta0}{hyperparameter for exponential prior}
  \item{keep_post_sampls}{returns posterior samples (or not)}
}
\details{
After setting up some preliminary data structures, this function simulates the data for this trial then proceeds through the interim analyses.
At each interim, the data accrued to date is used to model the adverse events for the coadministered covid and flu and the coadministered covid and placebo for each brand and cohort (where cohort denotes being enrolled at either dose 1 of the covid vaccine or dose 2 of the covid vaccine).

After all the interim analyses have been completed, a final analysis is run on all the data after followup is completed then the trial dataset, analysis summaries and config are returned to the user.

Best used in parallel per the example.
}
\value{
List of data and parameters used in the trial and trial analysis results
\item{d}{Simulated dataset used for trial}
\item{lpar}{Simulation parameters}
\item{lr}{Summaries of analysis results}
\item{lpost}{Posterior samples (if requested) otherwise NULL}
}
\references{
}
\author{
Mark
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
\dontrun{
library(data.table)
library(coadministr)
# Same incr in Pr AE in CVD+FLU regardless of whether cvd dose 1 and dose 2
# cvd+flu combination increases pr AE by 0.1 for AZ
# cvd+flu combination increases pr AE by 0.15 for PF
b <- pr_adverse_evt(pae_pbo = 0.02, pae_flu = 0.05, 
                    pae_coad_az = 0.1, pae_coad_pf = 0.15,
                    pae_az1 = 0.5, pae_az2 = 0.4,
                    pae_pf1 = 0.3, pae_pf2 = 0.5)
d <- get_data(J = 1000,
              beta = b,
              sd_id = 0.5,
              n_per_yr = 1000,
              p_miss = 0.1)
results <- coadministr::trial(
   scenario = 1,
   looks = seq(500, 1000, by = 100),
   n_per_yr = 2000,
   delta_ni = 0.2, 
   beta = b,
   sd_id = 0.5,
   p_miss = 0.1,
   Ncoad = 2,
   Nbrand = 2,
   Ncohort = 2,
   a0 = 1,
   b0 = 1,
   eta0 = 0.01,
   keep_post_sampls = TRUE
 )
 # Note that results$lpost[[1]]$pflu[, , 1] contains draws for AZ coad CVD+FLU 
 # at covid dose 1 and dose 2 in cols 1 and 2 respectively.
 # ditto for results$lpost[[1]]$pflu[, , 2] but for PF
 # ditto for results$lpost[[1]]$ppbo[, , 1] but for AZ CVD+PBO
 # ditto for results$lpost[[1]]$ppbo[, , 2] but for PF CVD+PBO
 str(results$lpost[[1]])
 post_to_data_table(results$lpost[[1]])
}
\dontrun{
library(parallel)
library(data.table)
ncores <- 10
results <- parallel::mclapply(X=1:10, mc.cores = ncores, FUN=function(x) {
 # should use tryCatch block
 b <- pr_adverse_evt(pae_pbo = 0.02, pae_flu = 0.05, pae_coad = 0.1)
 coadministr::trial(
   scenario = 1,
   looks = seq(500, 1000, by = 100),
   n_per_yr = 2000,
   delta_ni = 0.2, 
   beta = b,
   sd_id = 0.5,
   p_miss = 0.1,
   Ncoad = 2,
   Nbrand = 2,
   Ncohort = 2,
   a0 = 1,
   b0 = 1,
   eta0 = 0.01
 )
})
}
}


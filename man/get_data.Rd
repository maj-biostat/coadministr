\name{get_data}
\alias{get_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Simulate trial data}
\description{
Simulates data for the co-administration trial, design 2.
}
\usage{
get_data(J = 100, beta = NULL, sd_id = 0.5, 
    n_per_yr = 100, p_miss = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{J}{number of unique individuals}
  \item{beta}{required Pr AE by time, brand, combination see \code{\link{pr_adverse_evt}}}
  \item{sd_id}{subject level variability for repeat measure}
  \item{n_per_yr}{enrol per year}
  \item{p_miss}{Pr observation missing}
}
\details{
Individual level probabilities of occurrence of AE. For each participant there is a fixed 50\% that they are enrolled at the first or second CVD dose and there is a 70\% chance that they get AZ brand.

All participants have Pr AE conditional on which cohort they were enrolled in (covid dose 1 or 2), the brand they received and the timepoint. The first timepoint is that of the initial coadministered dose which is either with a coadministered flu vaccination or coad with placebo. The second timepoint is that of the vaccine given at seven days post the coadministered shot. There is random variability added to both the time of the 7 day follow up and the subsequent 7 day follow up to assess whether any AE occurred. As individuals have repeat measure, they are given a subject level variability.
}
\value{
\code{data.table} containing participant level data
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Mark
}
\note{
%%  ~~further notes~~
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
View(d)              
}
}

\name{pr_adverse_evt}
\alias{pr_adverse_evt}
\title{
Simulate trial data
}
\description{
Creates the trial data for rct comparing proportion of participants having moderate or adverse events associated with coadminstration of covid19 and
flu relative to coadministration of covid19 and placebo.
}
\usage{
pr_adverse_evt(pae_pbo = 0.02, pae_flu = 0.05, 
    pae_coad_az = 0, pae_coad_pf = 0,
    pae_az1 = 0.5, pae_az2 = 0.4,
    pae_pf1 = 0.3, pae_pf2 = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pae_pbo}{Probability of 7 day AE PBO group (single dose)}
  \item{pae_flu}{Probability of 7 day AE FLU group (single dose)}
  \item{pae_coad_az}{Increase in Pr 7 day AE in CVD+FLU for AZ}
  \item{pae_coad_pf}{Increase in Pr 7 day AE in CVD+FLU for PF}
  \item{pae_az1}{Pr AE for first cvd dose cohort in AZ brand}
  \item{pae_az2}{Pr AE for second cvd dose cohort in AZ brand}
  \item{pae_pf1}{Pr AE for first cvd dose cohort in PF brand}
  \item{pae_pf2}{Pr AE for second cvd dose cohort in PF brand}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\code{data.table} containing timepoint, arm, cohort, probability of AE
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
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.

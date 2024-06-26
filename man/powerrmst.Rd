% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmstnormal.R
\name{powerrmst}
\alias{powerrmst}
\title{Power calculation for testing restricted mean survival time difference adjusted for a normal covariate}
\usage{
powerrmst(alpha0, alpha1, gamma0, beta, crate, tau, maxE, n, effect, NN, alpha)
}
\arguments{
\item{alpha0}{Parameter to specify in Weibull model. See Details for more information.}

\item{alpha1}{Parameter to specify in Weibull model. See Details for more information. \eqn{\alpha_1  = 0 } means there are proportional hazards; \eqn{\alpha_1 \neq 0 }  means the proportional hazards assumption is violated}

\item{gamma0}{Parameter to specify in Weibull model. See Details for more information.}

\item{beta}{Coefficient of normal covariate}

\item{crate}{Censoring rate, assumes an exponential distribution}

\item{tau}{Pre-specified survival horizon time for adjusted restricted mean survival times}

\item{maxE}{Maximum enrollment time. Assumes uniform enrollment between [0, maxE]}

\item{n}{Sample size per group}

\item{effect}{Targeted effect size}

\item{NN}{Number of iterations used to calculate the maximum information}

\item{alpha}{Targeted type I error rate

 @inherit Imaxrmst references}
}
\value{
Power for a trial with the given parameters.
}
\description{
Power calculation for testing restricted mean survival time difference adjusted for a normal covariate
}
\details{
Calculates the power for testing survival probability difference given the sample size and effect size.
See Details section in \code{\link{Imaxasp}} on how trial data are simulated.
By using MC Monte Carlo replicates of datasets of total sample size N,
\eqn{I_{0}^N} and \eqn{I_{effect}^N} can be calculated,
which are the  information for a trial with total sample size N
 under the null hypothesis with effect size 0 and under a targeted alternative hypothesis
 where the effect size is \eqn{effect}, respectively.
This allows us to calculate \eqn{V_{0}^N = 1/I_{0}^N} and
  \eqn{V_{effect}^N=1/I_{effect}^N}, which are the variance of the estimated difference
   under the respective null and targeted alternative hypotheses for a fixed sample trial
    with total sample size N.
Using the large sample normal distribution of the effect size estimators, the power
 \eqn{\pi} can be expressed as
\deqn{\pi = \Phi \left(\frac{effect \sqrt{N} - z_{1-\alpha/2} \sqrt{V_0^N}}{\sqrt{V_{effect}^N}}\right),}
  where \eqn{\alpha} is the targeted type I error rate and \eqn{z_q=\Phi^{-1}(q)}
   for any \eqn{q \in (0,1)}.
}
\examples{
\dontrun{
set.seed(30)
#1397.7830
powerrmst(alpha0=1.5, alpha1=-0.3, gamma0=-log(0.4), beta=0, crate=0, tau=1,
maxE=2, n=199, effect=0.088, NN=200, alpha=0.05)
}


}

#' @title Calculates \eqn{I_{max}} for adjusted survival probability difference for a normal covariate
#'
#'
#' @description Calculates \eqn{I_{max}}, the maximum information for the trial, via Monte Carlo simulation.
#' This function simulates NN number of trials, calculates the
#' adjusted survival probability (SP) difference between the treatment and control group for each of the NN trials,
#'  and then takes the inverse of the variance of these differences to obtain \eqn{I_{max}}.
#'
#' @return Maximum information for a trial with the given parameters.
#'
#' @details Trial data are simulated using the following assumptions.
#' Assume the event time for a subject follows a Weibull distribution with survival function
#' \deqn{S(t|Z_W, \mathbf{Z}) = \exp(-\gamma t^{\alpha}),}
#' where the shape and rate parameters
#' \eqn{\gamma} and \eqn{\alpha}  depend on \eqn{Z_W} (treatment group indicator)
#' and \eqn{\mathbf{Z}} (vector of baseline covariates):
#'  \deqn{\alpha = \alpha_0 + \alpha_1Z_W,}
#' \deqn{\gamma = \gamma_0 \exp(\beta_WZ_W + \boldsymbol{\beta}^T \mathbf{Z}).}
#' \eqn{\alpha_0},  \eqn{\alpha_1}, and \eqn{\gamma_0} are parameters that are chosen to
#' simulate trial data under a proportional hazards model, delayed effect setting, or crossing curves setting.
#'
#' @param alpha0 Parameter to specify in Weibull model. See Details for more information.
#' @param alpha1 Parameter to specify in Weibull model. See Details for more information. \eqn{\alpha_1  = 0 } means there are proportional hazards; \eqn{\alpha_1 \neq 0 }  means the proportional hazards assumption is violated
#' @param gamma0 Parameter to specify in Weibull model. See Details for more information.
#' @param beta Coefficient of normal covariate
#' @param crate Censoring rate, assumes an exponential distribution
#' @param t0 Pre-specified time at which adjusted SPs for each group are calculated
#' @param maxE Maximum enrollment time. Assumes uniform enrollment between [0, maxE]
#' @param n Sample size per group
#' @param effect Targeted effect size
#' @param NN Number of iterations used to calculate the maximum information
#'
#'
#' @references Zhang, P.K., Logan, B.L., and Martens, M.J. (2024). Covariate-adjusted Group Sequential Comparisons of Survival Probabilities. \emph{arXiv}. https://arxiv.org/abs/2403.17117.
#' @references Zhang, X., Loberiza, F. R., Klein, J. P., and Zhang, M.-J. (2007). A SAS Macro for
#' Estimation of Direct Adjusted Survival Curves Based on a Stratified Cox Regression
#' Model. \emph{Comput Methods Programs Biomed} \strong{88(2)}, 95–101.
#' @references Zucker, D.M. (1998) Restricted Mean Life with Covariates: Modification and Extension
#' of a Useful Survival Analysis Method. \emph{J Am Stat Assoc} \strong{93(442)}, 702-709
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(567)
#' #467.13
#' Imaxasp(alpha0=1.5, alpha1=-1, gamma0=-log(0.4), beta=log(1.5), crate=-log(0.95), t0=1,
#' maxE=2, n=200, effect=0, NN=10000)
#' }
Imaxasp <- function(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n, effect, NN) {
  asp.diff.est = NULL
  asp.diff.se = NULL
  minE = 0
  beta1 = rootasp(alpha0, alpha1, gamma0, beta, t0, effect)
  umax = t0 + maxE
  for (i in 1:NN)
  {

    E = runif(2*n, min=minE, max=maxE)          # enrollment times
    Z1 = (c(rep(0,n),rep(1,n)))                   # treatment indicator
    Z2 = as.matrix(rnorm(2*n))                             # covariates
    alpha = alpha0+alpha1*Z1
    gamma1 = gamma0*exp(beta1*Z1+beta*Z2)
    FT = rweibull(2*n,shape=alpha,scale=gamma1**(-1/alpha))
    CT = NULL
    if (crate == 0)
    {CT = Inf }
    else {CT = rexp(2*n, rate=crate)}
    Data = cbind(E,FT,CT,Z1,Z2)   #create dataset of enrollment times, failure times, censorimg times, treatment, covariates
    #Data
    #Data = Data[order(Data[,1],decreasing=FALSE),]

    # only use data up to max calendar time u
    Data = Data[Data[,1]<umax,]
    X = pmin(Data[,2], Data[,3], umax-Data[,1], t0)                    #calculate observed time X
    delta = as.numeric(Data[,2] < pmin(Data[,3],umax-Data[,1],t0))   #calculate censoring indicator
    #table(delta)
    #cbind(X,delta)

    # do asp calculations

    asp = asp(t0 = t0, Time = X, Status = delta, Z = as.matrix(Data[, 5:dim(Data)[2]]), TRT= Data[, 4])
    asp.diff.est[i] = as.numeric(asp$SPD)
    # asp.diff.se[i] = as.numeric(aspest$V2)
    #adjsp.diff.ci = adjsp.diff.est + c(-1,1)*qnorm(1-type1/2)*adjsp.diff.se
    #adjsp.Z = adjsp.diff.est / adjsp.diff.se
    #cbind(adjsp.diff.est, adjsp.diff.se, adjsp.Z, adjsp.diff.ci[1], adjsp.diff.ci[2])
    # calculate max information for fixed design
  }

  Imax = 1/var(asp.diff.est) #1047.361
  #Imax =1047.361
  return(Imax)
}





#' @title Power calculation for testing survival probability difference adjusted for a normal covariate
#'
#' @return Power for a trial with the given parameters.
#'
#' @details Calculates the power for testing survival probability difference given the sample size and effect size.
#' See Details section in \code{\link{Imaxasp}} on how trial data are simulated.
#' By using MC Monte Carlo replicates of datasets of total sample size N,
#' \eqn{I_{0}^N} and \eqn{I_{effect}^N} can be calculated,
#' which are the  information for a trial with total sample size N
#'  under the null hypothesis with effect size 0 and under a targeted alternative hypothesis
#'  where the effect size is \eqn{effect}, respectively.
#' This allows us to calculate \eqn{V_{0}^N = 1/I_{0}^N} and
#'   \eqn{V_{effect}^N=1/I_{effect}^N}, which are the variance of the estimated difference
#'    under the respective null and targeted alternative hypotheses for a fixed sample trial
#'     with total sample size N.
#' Using the large sample normal distribution of the effect size estimators, the power
#'  \eqn{\pi} can be expressed as
#' \deqn{\pi = \Phi \left(\frac{effect \sqrt{N} - z_{1-\alpha/2} \sqrt{V_0^N}}{\sqrt{V_{effect}^N}}\right),}
#'   where \eqn{\alpha} is the targeted type I error rate and \eqn{z_q=\Phi^{-1}(q)}
#'    for any \eqn{q \in (0,1)}.
#'
#'
#'
#'
#'
#' @inherit Imaxasp details
#' @inherit Imaxasp references
#'
#' @inheritParams Imaxasp
#' @param alpha Targeted type I error rate
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(14)
#'  #0.8012848
#'  powerasp(alpha0 = 1.5, alpha1= -1, gamma0=-log(0.4), beta=log(1.5), crate=-log(0.95), t0=1,
#'  maxE=2, n = 190, effect=0.135, NN=10000, alpha=0.05) #0.797
#' }
powerasp <- function(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n, effect, NN, alpha) {
  Veffect = 1/Imaxasp(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n, effect, NN)
  V0= 1/Imaxasp(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n, effect = 0, NN)

  zbeta = (sqrt(n)*effect - qnorm(1-alpha/2)*sqrt(V0*n))/sqrt(Veffect*n)
  if (effect ==0)
    return(pnorm(zbeta)*2)
  else
    return(pnorm(zbeta))
}





#' @title Sample size calculation for testing survival probability difference adjusted for a normal covariate
#'
#' @return Sample size for a trial with the given parameters
#'
#' @details Calculates the sample size for testing survival probability difference given power and effect size via Monte Carlo simulation.
#'  See Details section in \code{\link{Imaxasp}} on how trial data are simulated.
#'  By using MC replicate datasets with an initial sample size M, we can calculate
#'  \eqn{V_{0}^M=1/I_{0}^M} and \eqn{V_{effect}^M=1/I_{effect}^M}, which represent the variances
#'  of the estimated differences in M patient trials.
#'  Under the large sample normal distribution for the effect size estimator,
#'   the required sample size depends on the unit variance,
#'   which is the reciprocal of the amount of information contributed by a single patient.
#'   To obtain the unit variances under the null and targeted alternative hypotheses,
#'    \eqn{V_{0}^1} and \eqn{V_{effect}^1}, we rescale \eqn{V_{0}^M} and \eqn{V_{effect}^M} by calculating
#' \eqn{V_{0}^1 = V_{0}^M \cdot M } and \eqn{V_{effect}^1 = V_{effect}^M \cdot M }.
#' With the individual variances per patient, we can calculate the actual required total sample size N for the trial using the equation
#' \deqn{N = \frac{\left(z_{1-\alpha/2}\sqrt{V_0^1} + z_{\pi} \sqrt{V_{effect}^1}\right)^2}{effect^2},}
#'
#'
#'
#' @inherit powerasp details
#' @inherit Imaxasp references
#'
#' @inheritParams powerasp
#' @param m Sample size used to calculate the maximum information, Imax
#' @param pi Targeted power
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  set.seed(10)
#'  Nasp(alpha0 = 1.5, alpha1=-1, gamma0=-log(0.4), beta=log(1.5), crate=-log(0.95),
#'   t0=1, maxE=2, m=500,
#'   effect=0.135, NN=10000, alpha=0.05, pi = 0.8) #189.78*2
#' }
Nasp <- function(alpha0, alpha1, gamma0, beta, crate, t0, maxE, m, effect, NN, alpha,pi=0.8) {
  Veffect = 1/Imaxasp(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n=m, effect, NN)
  V0= 1/Imaxasp(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n=m, effect = 0, NN)
  M = m
  N = (qnorm(1-alpha/2)*sqrt(V0*M) + qnorm(pi)*sqrt(Veffect*M))^2/effect^2

  return(N)
}








#' @title Effect size calculation for testing survival probability difference adjusted for a normal covariate
#'
#' @return Effect size for a trial given the parameters
#'
#' @details Calculates the effect size for testing survival probability difference given sample size and power via Monte Carlo simulation.
#'  See Details section in \code{\link{Imaxasp}} on how trial data are simulated.
#'   Because the targeted effect size impacts the sample size formula through the two terms
#'    \eqn{effect} and \eqn{V_{effect}^N},
#'     we use an iterative procedure to obtain the required \eqn{effect}.
#' We first calculate an initial variance in the treatment group, \eqn{V_{10}^N},
#' by calculating \eqn{I_{0}^N} and initializing \eqn{V_{10}^N = 1/I_{0}^N}.
#' Then, the initial effect size is calculated as
#' \deqn{effect_0 = z_{1-\alpha/2} \sqrt{V_0^N} + z_{\pi} \sqrt{V_{10}^N}.}
#' Using this effect size, \eqn{V_{effect_0}^N} is calculated by taking the reciprocal
#'  of the information for a trial with effect size \eqn{effect_0}.
#'  \eqn{V_{10}^N} is then updated to \eqn{V_{effect_0}^N}.
#'  This updating of \eqn{effect_0} and \eqn{V_{10}^N} is repeated
#  until convergence occurs.
#'
#'
#'
#'
#'
#'
#'
#'
#'
#' @inherit Imaxasp references
#'
#' @inheritParams Imaxasp
#'
#' @param alpha Targeted type I error rate
#' @param pi Targeted power
#' @param max.iter Maximum number of iterations to calculate the effect size
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1234)
#' ESasp(alpha0 = 1.5, alpha1=-1, gamma0=-log(0.4), beta=log(1.5), crate=-log(0.95), t0=1,
#' maxE=2, n=190, NN = 10000, alpha=0.05, pi=0.8, max.iter=10) #0.1340863
#' }
ESasp <- function(alpha0, alpha1, gamma0, beta, crate, t0,
                  maxE, n, NN, alpha=0.05, pi = 0.8, max.iter){
  zalpha = qnorm(1-alpha/2)
  zbeta = qnorm(pi)
  beta1 = NULL
  asp.diff.est = NULL
  #N = n
  Imax = Imaxasp(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n, effect = 0, NN)
  Vnull = 1/Imax
  V10 = 1/Imax
  V11 = 1/Imax
  effect = NULL
  for (j in 1:max.iter)
  {
    effect = (zalpha * sqrt(Vnull*n) + zbeta * sqrt(V10*n)) / sqrt(n)      # calculate delta

    Imax1 = Imaxasp(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n, effect, NN)
    V11 = 1 / Imax1
    diff = abs(V11-V10)
    #diff
    if (diff/V10 < 0.005)
    {
      if (V11 > V10)
      {
        effect = (zalpha*sqrt(Vnull*n) + zbeta*sqrt(V11*n)) / sqrt(n)    #calculate delta
        beta1 = rootasp(alpha0, alpha1, gamma0, beta, t0, effect)
      }    else {
        effect = (zalpha*sqrt(Vnull*n) + zbeta*sqrt(V10*n)) / sqrt(n)    #calculate delta
        beta1 =  rootasp(alpha0, alpha1, gamma0, beta, t0, effect)
      }
      break
    }
    V10 = V11
    cat("beta1: ",beta1, "\n")
  }
  return(list(beta1 = beta1, effect = effect))
}

#' @title Calculates \eqn{I_{max}} for adjusted SPs for a normal covariate
#'
#'
#' @description Calculates \eqn{I_{max}}, the maximum information for the trial, via Monte Carlo simulation.
#' This function simulates NN number of trials, calculates NN number of
#' adjusted SP differences between the treatment and control group, and
#' then takes the inverse of the variances of these differences to obtain \eqn{I_{max}}.
#'
#' @details Trial data are simulated using the following assumptions.
#' Assume the event time for a subject follows a Weibull distribution with survival function
#' \deqn{S(t|Z_W, \mathbf{Z}) = \exp(-\gamma t^{\alpha}),}
#' where the shape and rate parameters
#' \eqn{\gamma} and \eqn{\alpha}  depend on \eqn{Z_W} (treatment group indicator)
#' and \eqn{\mathbf{Z}} (vector of baseline covariates):
#'  \deqn{\alpha = \alpha_0 + \alpha_1Z_W,}
#' \deqn{\gamma = \gamma_0 \exp(\beta_WZ_W + \mathbf{\beta}^T \mathbf{Z}).}
#' \eqn{\alpha_0},  \eqn{\alpha_1}, and \eqn{\gamma_0} are parameters that are chosen to
#' simulate trial data under a proportional hazards model, delayed effect setting, or crossing curves setting.
#'
#' @param alpha0 - Parameter to specify in Weibull model. See Details for more information.
#' @param alpha1 - Parameter to specify in Weibull model. See Details for more information. \eqn{\alpha_1  = 0 } means there are proportional hazards; \eqn{\alpha_1 \neq 0 }  means the proportional hazards assumption is violated
#' @param gamma0 - Parameter to specify in Weibull model. See Details for more information.
#' @param beta2 - Coefficient of normal covariate
#' @param crate - Censoring rate, assumes an exponential distribution
#' @param t0 - Pre-specified time at which adjusted SPs for each group are calculated
#' @param maxE - Maximum enrollment time. Assumes uniform enrollment between [0, E]
#' @param n - Sample size per group
#' @param effect - Targeted effect size
#' @param NN - Number of iterations used to calculate the maximum information
#'
#'
#' @references Zhang, P.K., Logan, B.L., and Martens, M.J. (2024). Covariate-adjusted Group Sequential Comparisons of Survival Probabilities. \emph{arXiv}
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
#' Imaxasp(alpha0=1.5, alpha1=-1, gamma0=-log(0.4), beta2=0, crate=0, t0=1,
#' maxE=2, n=200, effect=0, NN=10000)
#' }
Imaxasp <- function(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect, NN) {
  asp.diff.est = NULL
  asp.diff.se = NULL
  minE = 0
  beta1 = rootasp(alpha0, alpha1, gamma0, beta2, t0, effect)
  umax = t0 + maxE
  for (i in 1:NN)
  {

    E = runif(2*n, min=minE, max=maxE)          # enrollment times
    Z1 = (c(rep(0,n),rep(1,n)))                   # treatment indicator
    Z2 = as.matrix(rnorm(2*n))                             # covariates
    alpha = alpha0+alpha1*Z1
    gamma1 = gamma0*exp(beta1*Z1+beta2*Z2)
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

    # do rmst calculations

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





#' @title Power calculation for a normal covariate
#'
#' @description Calculates the power given sample size and effect size via Monte Carlo simulation.
#' This function first calculates \eqn{V_{0}} and  \eqn{V_{effect}},
#' the variance in the control group and treatment group, respectively,
#' by calculating the maximum information for both groups.
#' The power is then calculated using the equation
#' \deqn{(\sqrt(N)*effect - z_{1-alpha/2}*\sqrt(V_0*N))/sqrt(V_{effect}*N)}
#'
#' \eqn{I_{max}}, the maximum information for the trial, via Monte Carlo simulation.
#' This function simulates NN number of trials, calculates NN number of
#' adjusted SP differences between the treatment and control group, and
#' then takes the inverse of the variances of these differences to obtain \eqn{I_{max}}.
#'
#' @details See Details section in \code{\link{Imaxasp}} on how trial data are simulated.
#'
#' @inherit Imaxasp details
#' @inherit Imaxasp references
#'
#' @inheritParams Imaxasp
#' @param alpha - Targeted type I error rate
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  powerasp(alpha0 = 1.5, alpha1= -1, gamma0=-log(0.4), beta2=0, crate=0, t0=1,
#'  maxE=2, n = 191, effect=0.14, NN=10000, alpha=0.05)
#' }
powerasp <- function(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect, NN, alpha) {
  Veffect = 1/Imaxasp(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect, NN)
  V0= 1/Imaxasp(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect = 0, NN)
  N = n
  zbeta = (sqrt(N)*effect - qnorm(1-alpha/2)*sqrt(V0*N))/sqrt(Veffect*N)
  if (effect ==0)
    return(pnorm(zbeta)*2)
  else
    return(pnorm(zbeta))
}





#' @title Sample size calculation for a normal covariate
#'
#' @description Calculates the sample size given power and effect size via Monte Carlo simulation.
#' This function first calculates \eqn{V_{0}} and  \eqn{V_{effect}},
#' the variance in the control group and treatment group, respectively,
#' by calculating the maximum information for both groups using an intial sample size of m .
#' The sample size N is then calculated using the equation
#' \deqn{N = (z_{1-\alpha/2}*\sqrt(V_0*M) + z_{1-\beta}*\sqrt(V_{effect}*M))^2/effect^2}
#'
#' @inherit powerasp details
#' @inherit Imaxasp references
#'
#' @inheritParams powerasp
#' @param m - Sample size used to calculate the maximum information, Imax
#' @param beta - Targeted type II error rate
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#'  set.seed(1234)
#'  Nasp(alpha0 = 1.5, alpha1=-1, gamma0=-log(0.4), beta2=0, crate=0, t0=1, maxE=2, m=400,
#'   effect=0.14, NN=10000, alpha=0.05, beta = 0.2)
#' }
Nasp <- function(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, m, effect, NN, alpha,beta) {
  Veffect = 1/Imaxasp(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n=m, effect, NN)
  V0= 1/Imaxasp(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n=m, effect = 0, NN)
  M = m
  N = (qnorm(1-alpha/2)*sqrt(V0*M) + qnorm(1-beta)*sqrt(Veffect*M))^2/effect^2

  return(N)
}








#' @title Effect size calculation for a normal covariate
#'
#' @description Calculates the effect size given sample size and power via Monte Carlo simulation
#' This function first calculates an initial variance in the treatment group, \eqn{V_{10}},
#'  by calculating the maximum information
#' for an effect size of 0.
#' Then the initial effect size is calculated as #'
#' \deqn{effect = (zalpha * sqrt(V0*N) + zbeta * sqrt(V10*N)) / sqrt(N)}
#' Using this effect size, \eqn{V_{11}} is calculated and
#'  set as the updated variance of the treatment group.
#'  \eqn{V_{10}} is set as \eqn{V_{11}} and this process is repeated
#'  until convergence occurs.
#'
#'
#'
#'
#' @inherit powerasp details
#'
#' @description Calculate effect size given power and sample size
#'
#'
#' @inherit Imaxasp references
#'
#' @inheritParams Imaxasp
#'
#' @param alpha - Targeted type I error rate
#' @param beta - Targeted type II error rate
#' @param max.iter - Maximum number of iterations to calculate the effect size
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ESasp(alpha0 = 1.5, alpha1=-1, gamma0=-log(0.4), beta2=0, crate=0, t0=1,
#' maxE=2, n=191, NN = 10000, alpha=0.05, beta = 0.2, max.iter=10)
#' }
ESasp <- function(alpha0, alpha1, gamma0, beta2, crate, t0,
                  maxE, n, NN, alpha=0.05, beta = 0.2, max.iter){
  zalpha = qnorm(1-alpha/2)
  zbeta = qnorm(1-beta)
  beta1 = NULL
  asp.diff.est = NULL
  N = n
  Imax = Imaxasp(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect = 0, NN)
  Vnull = 1/Imax
  V10 = 1/Imax
  V11 = 1/Imax
  effect = NULL
  for (j in 1:max.iter)
  {
    effect = (zalpha * sqrt(Vnull*N) + zbeta * sqrt(V10*N)) / sqrt(N)      # calculate delta

    Imax1 = Imaxasp(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect, NN)
    V11 = 1 / Imax1
    diff = abs(V11-V10)
    #diff
    if (diff/V10 < 0.005)
    {
      if (V11 > V10)
      {
        effect = (zalpha*sqrt(Vnull*N) + zbeta*sqrt(V11*N)) / sqrt(N)    #calculate delta
        beta1 = rootasp(alpha0, alpha1, gamma0, beta2, t0, effect)
      }    else {
        effect = (zalpha*sqrt(Vnull*N) + zbeta*sqrt(V10*N)) / sqrt(N)    #calculate delta
        beta1 =  rootasp(alpha0, alpha1, gamma0, beta2, t0, effect)
      }
      break
    }
    V10 = V11
    cat("beta1: ",beta1, "\n")
  }
  return(list(beta1 = beta1, effect = effect))
}

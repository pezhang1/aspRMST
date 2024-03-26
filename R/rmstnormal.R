
#' Title
#'
#' @param alpha0 - parameter to specify in Weibull model
#' @param alpha1 - parameter to specify in Weibull model. alpha1 = 0 means there are proportional hazards; alpha1 != 0 means the proportional hazards assumption is violated
#' @param gamma0 - parameter to specify in Weibull model
#' @param beta2 -vector of coefficients for non-treatment group binary variables
#' @param crate - censoring rate, assumes an exponential distribution
#' @param t0 - pre-specified time at which adjusted restricted mean survival times for each group are calculated
#' @param maxE - maximum enrollment time. Assumes uniform enrollment between [0,E]
#' @param n - sample size per group
#' @param effect - targeted effect size
#' @param NN - number of iterations
#'
#' @return
#' @export
#'
#' @examples
Imaxrmst <- function(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect, NN) {
  minE = 0
  beta1 = root(alpha0, alpha1, gamma0, beta2, t0, effect)
  umax = t0 + maxE
  rmst.diff.est = NULL
  rmst.diff.se = NULL
  for (i in 1:NN)
  {
    E = runif(2*n, min=minE, max=maxE)          # enrollment times
    Z1 = (c(rep(0,n),rep(1,n)))                   # treatment indicator
    Z2 = as.matrix(rnorm(2*n))                             # covariates
    alpha = alpha0+alpha1*Z1
    gamma1 = gamma0*exp(beta1*Z1+beta2*Z2)
    FT = rweib(2*n, alpha, gamma1)
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

    rmstest = rmst(t0 = t0, Time = X, Status = delta, Z = as.matrix(Data[, 5:dim(Data)[2]]), TRT= Data[, 4])
    rmst.diff.est[i] = rmstest$Delta
    rmst.diff.se[i] = as.numeric(rmstest$DSE2)
    #adjsp.diff.ci = adjsp.diff.est + c(-1,1)*qnorm(1-type1/2)*adjsp.diff.se
    #adjsp.Z = adjsp.diff.est / adjsp.diff.se
    #cbind(adjsp.diff.est, adjsp.diff.se, adjsp.Z, adjsp.diff.ci[1], adjsp.diff.ci[2])
    # calculate max information for fixed design
  }

  Imax = 1/var(rmst.diff.est) #1047.361
  #Imax =1047.361
  return(Imax)
}









#' Title
#'
#' @param alpha0 - parameter to specify in Weibull model
#' @param alpha1 - parameter to specify in Weibull model. alpha1 = 0 means there are proportional hazards; alpha1 != 0 means the proportional hazards assumption is violated
#' @param gamma0 - parameter to specify in Weibull model
#' @param beta2 -vector of coefficients for non-treatment group binary variables
#' @param crate - censoring rate, assumes an exponential distribution
#' @param t0 - pre-specified time at which adjusted restricted mean survival times for each group are calculated
#' @param maxE - maximum enrollment time. Assumes uniform enrollment between [0,E]
#' @param n - sample size per group
#' @param effect - targeted effect size
#' @param NN - number of iterations
#' @param alpha - targeted type I error rate
#'
#' @return
#' @export
#'
#' @examples
powerrmst <- function(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect, NN, alpha) {
  Veffect = 1/Imaxrmst(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect, NN)
  V0= 1/Imaxrmst(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect = 0, NN)
  N = n
  zbeta = (sqrt(N)*effect - qnorm(1-alpha/2)*sqrt(V0*N))/sqrt(Veffect*N)
  if (effect ==0)
    return(pnorm(zbeta)*2)
  else
    return(pnorm(zbeta))
}









#' Title
#'
#' @param alpha0 - parameter to specify in Weibull model
#' @param alpha1 - parameter to specify in Weibull model. alpha1 = 0 means there are proportional hazards; alpha1 != 0 means the proportional hazards assumption is violated
#' @param gamma0 - parameter to specify in Weibull model
#' @param beta2 -vector of coefficients for non-treatment group binary variables
#' @param crate - censoring rate, assumes an exponential distribution
#' @param t0 - pre-specified time at which adjusted restricted mean survival times for each group are calculated
#' @param maxE - maximum enrollment time. Assumes uniform enrollment between [0,E]
#' @param m - sample size used to calculate the maximum information, Imax
#' @param effect - targeted effect size
#' @param NN - number of iterations
#' @param alpha - targeted type I error rate
#' @param beta- targeted type II error rate
#'
#' @return
#' @export
#'
#' @examples
Nrmst <- function(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, m, effect, NN, alpha,beta) {
  Veffect = 1/Imaxrmst(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n=m, effect, NN)
  V0= 1/Imaxrmst(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n=m, effect = 0, NN)
  M = m
  N = (qnorm(1-alpha/2)*sqrt(V0*M) + qnorm(1-beta)*sqrt(Veffect*M))^2/effect^2

  return(N)
}








#' Title
#'
#' @param alpha0 - parameter to specify in Weibull model
#' @param alpha1 - parameter to specify in Weibull model. alpha1 = 0 means there are proportional hazards; alpha1 != 0 means the proportional hazards assumption is violated
#' @param gamma0 - parameter to specify in Weibull model
#' @param beta2 -vector of coefficients for non-treatment group binary variables
#' @param crate - censoring rate, assumes an exponential distribution
#' @param t0 - pre-specified time at which adjusted restricted mean survival times for each group are calculated
#' @param maxE - maximum enrollment time. Assumes uniform enrollment between [0,E]
#' @param n - sample size
#' @param NN - number of iterations
#' @param alpha - targeted type I error rate
#' @param beta - targeted type II error rate
#' @param max.iter - maximum number of iterations to calculate the effect size
#'
#' @return
#' @export
#'
#' @examples
ESrmst <- function(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, NN, alpha=0.05, beta = 0.2, max.iter){
  zalpha = qnorm(1-alpha/2)
  zbeta = qnorm(1-beta)
  beta1 = NULL
  asp.diff.est = NULL
  N = n
  Imax = Imaxrmst(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect = 0, NN)
  Vnull = 1/Imax
  V10 = 1/Imax
  V11 = 1/Imax
  effect = NULL
  for (j in 1:max.iter)
  {
    effect = (zalpha * sqrt(Vnull*N) + zbeta * sqrt(V10*N)) / sqrt(n)      # calculate delta

    Imax1 = Imaxrmst(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect, NN)
    V11 = 1 / Imax1
    diff = abs(V11-V10)
    #diff
    if (diff/V10 < 0.005)
    {
      if (V11 > V10)
      {
        effect = (zalpha*sqrt(Vnull*N) + zbeta*sqrt(V11*N)) / sqrt(N)    #calculate delta
        beta1 = root(alpha0, alpha1, gamma0, beta2, t0, effect)
      }    else {
        effect = (zalpha*sqrt(Vnull*N) + zbeta*sqrt(V10*N)) / sqrt(N)    #calculate delta
        beta1 =  root(alpha0, alpha1, gamma0, beta2, t0, effect)
      }
      break
    }
    V10 = V11
    #  cat("beta1final: ",beta1, "\n")
  }
  return(list(beta1 = beta1, effect = effect))
}
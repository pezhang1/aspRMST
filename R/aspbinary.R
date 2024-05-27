#' @title Calculates \eqn{I_{max}} for adjusted survival probability difference for binary covariates
#'
#' @inherit Imaxasp description
#' @inherit Imaxasp return
#' @inherit Imaxasp details
#'
#' @inherit Imaxasp references
#'
#' @inheritParams Imaxasp
#' @param beta Vector of coefficients for binary covariates
#' @param p Vector of probabilities for binary covariates
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(888)
#' Imaxaspbinary(alpha0=1.5, alpha1=-1, gamma0=-log(0.4), beta=log(2)/sqrt(2), crate=0, t0=1,
#' maxE=2, n=200, effect= 0, NN = 10000,p=c(0.5, 0.3)) #572.23
#' }
Imaxaspbinary <- function(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n, effect, NN,p) {
  asp.diff.est = NULL
  asp.diff.se = NULL
  minE = 0
  beta1 = rootaspbinary(alpha0, alpha1, gamma0, beta, t0, effect,p)
  umax = t0 + maxE
  for (i in 1:NN)
  {

    E = runif(2*n, min=minE, max=maxE)          # enrollment times
    Z1 = c(rep(0,n),rep(1,n))                   # treatment indicator
    Z0s = matrix(nrow = 2*n, ncol = length(p))
    for (j in 1:length(p))
    { Z0s[,j ] = (rbinom(2*n , size = 1 ,prob = p[j]))} # - p1) / sqrt(p1*(1-p1))                             # covariates
    # Z3 = (rbinom(2*n , size = 1 ,prob = p2) -p2) / sqrt(p2*(1-p2))
    alpha = alpha0+alpha1*Z1
    gamma1 = gamma0*exp(beta1*Z1 + rowSums(beta * Z0s))
    FT = rweibull(2*n,shape=alpha,scale=gamma1**(-1/alpha))
    CT = NULL
    if (crate == 0)
    {CT = Inf }
    else {CT = rexp(2*n, rate=crate)}
    Data = cbind(E,FT,CT,Z1,Z0s)   #create dataset of enrollment times, failure times, censorimg times, treatment, covariates
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


#' @title Power calculation for testing survival probability difference adjusted for binary covariates
#' @inherit powerasp return
#'
#' @inherit powerasp details
#'
#' @inheritParams powerasp
#' @param beta Vector of coefficients for binary covariates
#' @param p Vector of probabilities for binary covariates
#'
#'
#' @inherit Imaxasp references
#'@return
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(16)
#' poweraspbinary(alpha0 = 1.5, alpha1=-1, gamma0=-log(0.4), beta=log(2)/sqrt(2), crate=0, t0=1,
#' maxE=2, n=193, effect=0.121, NN=10000, alpha = 0.05, p =c(0.5, 0.3)) #0.800
#' }
poweraspbinary <- function(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n, effect, NN, alpha=0.05,p) {
  Veffect = 1/Imaxaspbinary(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n, effect, NN,p)
  V0= 1/Imaxaspbinary(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n, effect = 0, NN,p)
  N = n
  zbeta = (sqrt(n)*effect - qnorm(1-alpha/2)*sqrt(V0*n))/sqrt(Veffect*n)
  if (effect ==0)
    return(pnorm(zbeta)*2)
  else
    return(pnorm(zbeta))
}






#' @title Sample size calculation for testing survival probability difference adjusted for binary covariates
#' @inherit Nasp return
#'
#' @inherit Nasp details
#'
#' @inherit Imaxasp references
#'
#' @inheritParams Nasp
#' @param beta Vector of coefficients for binary covariates
#'
#'
#' @param p Vector of probabilities for binary covariates
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(8)
#' Naspbinary(alpha0 = 1.5, alpha1=-1, gamma0=-log(0.4), beta=log(2)/sqrt(2), crate=0, t0=1,
#'  maxE=2, m=500, effect=0.121, NN=10000, alpha=0.05, pi = 0.8, p =c(0.5, 0.3)) #192.27*2

#' }
#'
Naspbinary <- function(alpha0, alpha1, gamma0, beta, crate, t0, maxE, m, effect, NN, alpha=0.05,pi = 0.8, p) {
  Veffect = 1/Imaxaspbinary(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n=m, effect, NN, p)
  V0= 1/Imaxaspbinary(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n=m, effect = 0, NN, p)
  M = m
  N = (qnorm(1-alpha/2)*sqrt(V0*M) + qnorm(pi)*sqrt(Veffect*M))^2/effect^2

  return(N)
}
















#' @title Effect size calculation for testing survival probability difference adjusted for binary covariates
#'
#' @inherit ESasp return
#'
#'
#' @inheritParams ESasp
#'
#' @inherit ESasp details
#'
#' @inherit Imaxasp references
#' @param beta Vector of coefficients for binary covariates
#' @param p Vector of probabilities for binary covariates
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(23)
#' ESaspbinary(alpha0 = 1.5, alpha1 = -1, gamma0 = -log(0.4), beta=log(2)/sqrt(2), crate=0, t0=1,
#'  maxE=2, n=193, NN = 10000, alpha=0.05, pi = 0.8, max.iter=10, p =c(0.5,0.3)) #0.1218265
#' }
ESaspbinary <- function(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n, NN, alpha=0.05, pi=0.8, max.iter, p){
  zalpha = qnorm(1-alpha/2)
  zbeta = qnorm(pi)
  beta1 = NULL
  asp.diff.est = NULL
  #N = n
  Imax = Imaxaspbinary(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n, effect = 0, NN,p)
  Vnull = 1/Imax
  V10 = 1/Imax
  V11 = 1/Imax
  effect = NULL
  for (j in 1:max.iter)
  {
    effect = (zalpha * sqrt(Vnull*n) + zbeta * sqrt(V10*n)) / sqrt(n)      # calculate delta
    beta1 = rootaspbinary(alpha0, alpha1, gamma0, beta, t0, effect,p)
    Imax1 = Imaxaspbinary(alpha0, alpha1, gamma0, beta, crate, t0, maxE, n, effect, NN,p)
    V11 = 1 / Imax1
    diff = abs(V11-V10)
    #diff
    if (diff/V10 < 0.005)
    {
      if (V11 > V10)
      {
        effect = (zalpha*sqrt(Vnull*n) + zbeta*sqrt(V11*n)) / sqrt(n)    #calculate delta
        beta1 = rootaspbinary(alpha0, alpha1, gamma0, beta, t0, effect,p)
      }    else {
        effect = (zalpha*sqrt(Vnull*n) + zbeta*sqrt(V10*n)) / sqrt(n)    #calculate delta
        beta1 =  rootaspbinary(alpha0, alpha1, gamma0, beta, t0, effect,p)
      }
      break
    }
    V10 = V11
    cat("beta1: ",beta1, "\n")
  }
  return(list(beta1 = beta1, effect = effect))
}

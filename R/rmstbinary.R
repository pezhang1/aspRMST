

#' Title
#'
#' @inheritParams Imaxrmst
#' @param p - vector of probabilities for non-treatment group binary variables
#'
#' @inherit Imaxrmst references
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{}
Imaxrmstbinary <- function(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect, NN,p ) {
  minE = 0
  beta1 = rootbinary(alpha0, alpha1, gamma0, beta2, t0, effect,p)
  umax = t0 + maxE
  rmst.diff.est = NULL
  rmst.diff.se = NULL
  for (i in 1:NN)
  {
    E = runif(2*n, min=minE, max=maxE)          # enrollment times
    Z1 = c(rep(0,n),rep(1,n))                   # treatment indicator
    Z0s = matrix(nrow = 2*n, ncol = length(p))
    for (j in 1:length(p))
    { Z0s[,j ] = (rbinom(2*n , size = 1 ,prob = p[j]))} # - p1) / sqrt(p1*(1-p1))                             # covariates
    # Z3 = (rbinom(2*n , size = 1 ,prob = p2) -p2) / sqrt(p2*(1-p2))
    alpha = alpha0+alpha1*Z1
    gamma1 = gamma0*exp(beta1*Z1 + rowSums(beta2 * Z0s))
    FT = rweib(2*n, alpha, gamma1)
    CT = NULL
    if (crate == 0)
    {CT = Inf }
    else {CT = rexp(2*n, rate=crate)}
    Data = cbind(E,FT,CT,Z1,Z0s)   #create dataset of enrollment times, failure times, censorimg times, treatment, covariates
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

  Imax = (1/var(rmst.diff.est)) #1047.361
  #Imax =1047.361
  return(Imax)
}



#' Title
#'
#' @inheritParams powerrmst
#' @param p - vector of probabilities for non-treatment group binary variables
#'
#' @inherit Imaxrmst references
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{}
powerrmstbinary <- function(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect, NN, alpha,p) {
  Veffect = 1/Imaxrmstbinary(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect, NN,p)
  V0= 1/Imaxrmstbinary(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect = 0, NN,p)
  N = n
  zbeta = (sqrt(N)*effect - qnorm(1-alpha/2)*sqrt(V0*N))/sqrt(Veffect*N)
  if (effect ==0)
    return(pnorm(zbeta)*2)
  else
    return(pnorm(zbeta))
}



#' Title
#'
#' @inheritParams Nrmst
#' @param p - vector of probabilities for non-treatment group binary variables
#'
#' @inherit Imaxrmst references
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{}
Nrmstbinary <- function(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, m, effect, NN, alpha,beta,p) {
  Veffect = 1/Imaxrmstbinary(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n=m, effect, NN,p)
  V0= 1/Imaxrmstbinary(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n=m, effect = 0, NN,p)
  M = m
  N = (qnorm(1-alpha/2)*sqrt(V0*M) + qnorm(1-beta)*sqrt(Veffect*M))^2/effect^2

  return(N)
}



#' Title
#'
#' @inheritParams ESrmst
#' @param p - vector of probabilities for non-treatment group binary variables
#'
#' @inherit Imaxrmst references
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{}
ESrmstbinary <- function(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, NN, alpha=0.05, beta = 0.2, max.iter=10,p){
  zalpha = qnorm(1-alpha/2)
  zbeta = qnorm(1-beta)
  beta1 = NULL
  asp.diff.est = NULL
  N = n
  Imax = Imaxrmstbinary(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect = 0, NN,p)
  Vnull = 1/Imax
  V10 = 1/Imax
  V11 = 1/Imax
  effect = NULL
  for (j in 1:max.iter)
  {
    effect = (zalpha * sqrt(Vnull*N) + zbeta * sqrt(V10*N)) / sqrt(n)      # calculate delta

    Imax1 = Imaxrmstbinary(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n, effect, NN,p)
    V11 = 1 / Imax1
    diff = abs(V11-V10)
    #diff
    if (diff/V10 < 0.005)
    {
      if (V11 > V10)
      {
        effect = (zalpha*sqrt(Vnull*N) + zbeta*sqrt(V11*N)) / sqrt(N)    #calculate delta
        beta1 = rootbinary(alpha0, alpha1, gamma0, beta2, t0, effect,p)
      }    else {
        effect = (zalpha*sqrt(Vnull*N) + zbeta*sqrt(V10*N)) / sqrt(N)    #calculate delta
        beta1 =  rootbinary(alpha0, alpha1, gamma0, beta2, t0, effect,p)
      }
      break
    }
    V10 = V11
    #  cat("beta1final: ",beta1, "\n")
  }
  return(list(beta1 = beta1, effect = effect))
}

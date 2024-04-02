
#' Title
#'
#' @inheritParams Imaxasp
#' @param t0 - pre-specified time at which adjusted restricted mean survival times for each group are calculated
#' @return
#' @export
#'
#'
#'
#' @references Zhang, P.K., Logan, B.L., and Martens, M.J. (2024). Covariate-adjusted Group Sequential Comparisons of Survival Probabilities. \emph{arXiv}
#' @references Zhang, X., Loberiza, F. R., Klein, J. P., and Zhang, M.-J. (2007). A SAS macro for
#' estimation of direct adjusted survival curves based on a stratified Cox regression
#' model. \emph{Comput Methods Programs Biomed} \strong{88(2)}, 95–101.
#' @references Zhang, X. (2013). omparison of restricted mean survival times between treatments
#' based on a stratified cox model. \emph{Bio-Algorithms and Med-Systems} \strong{9(4)}, 183-189
#' @references Zucker, D.M. (1998) Restricted mean life with covariates: modification and extension
#' of a useful survival analysis method. \emph{J Am Stat Assoc} \strong{93(442)}, 702-709
#'
#'
#'
#' @examples
#' \dontrun{
#' set.seed(1234)
#' Imaxrmst(alpha0 = 1.5, alpha1 = -0.3, gamma0=-log(0.4), beta2=0, crate=0, t0=1, maxE=2, n=200, effect=0.08788154, NN = 100)
#' }
#'
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
#'
#' @inheritParams powerasp
#'
#'
#' @param t0 - pre-specified time at which adjusted restricted mean survival times for each group are calculated
#'
#'  @inherit Imaxrmst references
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1234)
#' powerrmst(alpha0=1.5, alpha1=-0.3, gamma0=-log(0.4), beta2=0, crate=0, t0=1, maxE=2, n=199, effect=0.08788154, NN=200, alpha=0.05)
#' }
#'
#'
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
#' @inheritParams Nasp
#' @param t0 - pre-specified time at which adjusted restricted mean survival times for each group are calculated
#'
#' @inherit Imaxrmst references
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1234)
#' Nrmst(alpha0 = 1.5, alpha1 = -0.3, gamma0 = -log(0.4), crate = 0, t0=1, maxE=2, m = 400, beta2=0, effect = 0.08788154, NN = 100, alpha = 0.05, beta=0.2)
#' }
Nrmst <- function(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, m, effect, NN, alpha,beta) {
  Veffect = 1/Imaxrmst(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n=m, effect, NN)
  V0= 1/Imaxrmst(alpha0, alpha1, gamma0, beta2, crate, t0, maxE, n=m, effect = 0, NN)
  M = m
  N = (qnorm(1-alpha/2)*sqrt(V0*M) + qnorm(1-beta)*sqrt(Veffect*M))^2/effect^2

  return(N)
}








#' Title
#'
#' @inheritParams powerasp
#' @param t0 - pre-specified time at which adjusted restricted mean survival times for each group are calculated
#' @inherit Imaxrmst references
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' set.seed(1234)
#' ESrmst(alpha0=1.5, alpha1=-0.3, gamma0=-log(0.4), beta2=0, crate=0, t0=1, maxE=2, n=199, NN=100, alpha=0.05, beta = 0.2, max.iter=10)
#' }
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

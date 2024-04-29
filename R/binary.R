#' @title Calculates \eqn{\beta_W} given a RMST effect size
#'
#'
#' @description Calculates \eqn{\beta_W}, the coefficient for the treatment group variable, that corresponds to a given RMST effect size
#'
#'
#'
#' @inheritParams  rootrmst
#' @param p - Vector of probabilities for binary covariates
#'
#'
#'
#' @export
#'
#' @examples
#' rootrmstbinary(alpha0 = 1.5, alpha1=-0.3, gamma0=-log(0.4), beta2 =0, t0 = 1,
#'  effect=0.1, p =c(0.2,0.5))
rootrmstbinary = function(alpha0, alpha1, gamma0, beta2, t0, effect,p ) {


  Diff=function(beta1)
  {
    Alpha0 = alpha0
    Alpha1 = alpha0+alpha1
    values <- expand.grid(replicate(length(p), 0:1, simplify = FALSE))
    standard <- matrix(nrow = dim(values)[1], ncol =dim(values)[2] )
    probs <- matrix(nrow = dim(values)[1], ncol =dim(values)[2] )
    for (i in 1:length(p)) {
      probs[, i] = p[i]^values[,i]*(1-p[i])^(1-values[,i])
      #    standard[, i] = (values[, i] - p[i]) / sqrt(p[i]*(1-p[i]))
    }
    Lambda0s = NULL
    Lambda1s = NULL
    # values = standard
    f0s = f1s =  NULL
    # num = length(p)
    for (i in 1:dim(values)[1]) {

      Lambda0s[i] = gamma0*exp((sum(beta2*values[i, ])))
      Lambda1s[i] = gamma0*exp(beta1+(sum(beta2*values[i, ])))
      f0s[i] = gamma(1+1/Alpha0) * pgamma(t0**Alpha0,shape=1/Alpha0,rate=Lambda0s[i]) / Lambda0s[i]**(1/Alpha0) *prod(probs[i,])
      f1s[i] = gamma(1+1/Alpha1) * pgamma(t0**Alpha1,shape=1/Alpha1,rate=Lambda1s[i]) / Lambda1s[i]**(1/Alpha1) *prod(probs[i,])
    }
    f0 = sum(f0s)
    f1 = sum(f1s)
    temp = f1-f0
    temp2 =temp - effect
    return(temp2)
  }
  return(uniroot(Diff, c(-10,10))$root)
}

















#' @title Calculates \eqn{\beta_W} given a SP effect size
#'
#' @description Calculates \eqn{\beta_W}, the coefficient for the treatment group variable, that corresponds to a given SP effect size
#'
#'
#'
#' @inheritParams  rootasp
#'
#' @param p - Vector of probabilities for binary covariates
#'
#'
#' @export
#'
#' @examples
#' rootaspbinary(alpha0 = 1.5, alpha1=-1, gamma0=-log(0.4), beta2 =0, t0 = 1,
#' effect=0.1334313, p =c(0.2,0.5))
rootaspbinary = function(alpha0, alpha1, gamma0, beta2, t0, effect,p) {


  Diff=function(beta1)
  {
    Alpha0 = alpha0
    Alpha1 = alpha0+alpha1
    values <- expand.grid(replicate(length(p), 0:1, simplify = FALSE))
    standard <- matrix(nrow = dim(values)[1], ncol =dim(values)[2] )
    probs <- matrix(nrow = dim(values)[1], ncol =dim(values)[2] )
    for (i in 1:length(p)) {
      probs[, i] = p[i]^values[,i]*(1-p[i])^(1-values[,i])
      #   standard[, i] = (values[, i] - p[i]) / sqrt(p[i]*(1-p[i]))
    }
    Lambda0s = NULL
    Lambda1s = NULL
    #   values = standard
    S0s = S1s =  NULL
    # num = length(p)
    for (i in 1:dim(values)[1]) {
      S0s[i] = exp(-gamma0*exp((sum(beta2*values[i, ])))*t0**alpha0)*prod(probs[i,])
      S1s[i] = exp(-gamma0*exp(beta1+(sum(beta2*values[i, ])))*t0**(alpha0+alpha1))*prod(probs[i,])
    }
    S0 = sum(S0s)
    S1 = sum(S1s)
    temp = S1-S0
    temp2 =temp - effect
    return(temp2)
  }
  return(uniroot(Diff, c(-10,10))$root)
}


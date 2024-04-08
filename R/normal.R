#' Title
#'
#' @param alpha0 - parameter to specify in Weibull model
#' @param alpha1 - parameter to specify in Weibull model. alpha1 = 0 means there are proportional hazards; alpha1 != 0 means the proportional hazards assumption is violated
#' @param gamma0 - parameter to specify in Weibull model
#' @param beta2 - vector of coefficients for non-treatment group binary variables
#' @param t0 - pre-specified time at which adjusted restricted mean survival times for each group are calculated
#' @param effect - targeted effect size
#'
#'
#' @export
#'
#' @examples
#' root(alpha0 = 1.5, alpha1=-0.3, gamma0=-log(0.4), beta2 =0, t0=1, effect=0.1)
root = function(alpha0, alpha1, gamma0, beta2, t0, effect) {

  Diff=function(beta1)
  {
    grmst2=function(z)
    {
      Alpha0 = alpha0
      Alpha1 = alpha0+alpha1
      Lambda0 = gamma0*exp(beta2*z)
      Lambda1 = gamma0*exp(beta1+beta2*z)
      f0 = gamma(1+1/Alpha0) * pgamma(t0**Alpha0,shape=1/Alpha0,rate=Lambda0) / Lambda0**(1/Alpha0)
      f1 = gamma(1+1/Alpha1) * pgamma(t0**Alpha1,shape=1/Alpha1,rate=Lambda1) / Lambda1**(1/Alpha1)
      temp = (f1-f0)*dnorm(z,0,1)
      return(temp)
    }
    temp2 = integrate(grmst2, -20, 20)$value - effect
    return(temp2)
  }
  return(uniroot(Diff, c(-10,10))$root)
}


#' Title
#'
#' @param alpha0 - parameter to specify in Weibull model
#' @param alpha1 - parameter to specify in Weibull model. alpha1 = 0 means there are proportional hazards; alpha1 != 0 means the proportional hazards assumption is violated
#' @param gamma0 - parameter to specify in Weibull model
#' @param beta2 - vector of coefficients for non-treatment group binary variables
#' @param t0 - pre-specified time at which adjusted survival probabilities for each group are calculated
#' @param effect - targeted effect size
#'
#'
#' @export
#'
#' @examples
#' rootasp(alpha0 = 1.5, alpha1=-1, gamma0=-log(0.4), beta2 =0, t0=1, effect=0.140996938)
rootasp = function(alpha0, alpha1, gamma0, beta2, t0, effect) {

  G = function(beta1)
  {
    f1 = function(z)
    {
      s1 = exp(-gamma0*exp(beta1+beta2*z)*t0**(alpha0+alpha1))
      f1 = s1*dnorm(z,0,1)
    }
    f0 = function(z)
    {
      s0 = exp(-gamma0*exp(beta2*z)*t0**alpha0)
      f0 = s0*dnorm(z,0,1)
    }
    S1 = integrate(f1, -10, 10)
    S0 = integrate(f0, -10, 10)
    f = S1$value - S0$value - effect
    return(f)
  }
  return(uniroot(G, c(-10,10))$root)
}




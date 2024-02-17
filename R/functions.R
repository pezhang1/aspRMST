



#' Title
#'
#' @param n
#' @param shape
#' @param rate
#'
#' @return
#' @export
#'
#' @keywords internal
#'
#' @examples
rweib = function(n,shape,rate)
{
  rweibull(n,shape=shape,scale=rate**(-1/shape))
}






#' Title
#'
#' @param alpha0 -
#' @param alpha1 -
#' @param gamma0 -
#' @param beta2 -
#' @param t0 -
#' @param effect -
#'
#' @return
#' @export
#'
#'
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

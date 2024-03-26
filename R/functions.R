



#' Title
#'
#' @param n - sample size
#' @param shape - shape parameter for Weibull distribution
#' @param rate - rate parameter for Weibull distribution
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




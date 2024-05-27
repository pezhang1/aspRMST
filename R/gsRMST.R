


#' Group sequential test of adjusted restricted mean survival times
#'
#' @inherit Imaxrmst references
#'
#' @description Performs a group sequential test of adjusted restricted mean survival times (RMST)s between two treatment groups
#' when provided the time point of analysis, the observed times, event indicator,
#' covariates, treatment group indicator, enrollment times, and calendar times of analysis.
#'
#' @param tau Pre-specified survival horizon time for adjusted restricted mean survival times
#' @param Time Observed times
#' @param Status Event indicator (0 = censored, 1 = observed)
#' @param Imax Maximum information
#' @param Z Non-treatment group covariates
#' @param TRT Treatment group indicator (0 = Control, 1 = Treatment)
#' @param E Enrollment times
#' @param alpha Targeted type I error rate
#' @param u Calendar times of analysis
#'
#' @returns
#'  \itemize{
#'   \item Z Vector of test statistics
#'   \item Crit Vector of critical values
#' }
#' @export
#'
#'
#'
#'
#'
#' @examples
#' \dontrun{
#' set.seed(1000)
#' alpha0 = 1.5
#' alpha1 = -0.3
#' beta1 = -0.61335306
#' n=200
#' crate =0
#' beta2=0
#' tau=1
#' gamma0 = -log(0.4)
#' maxE = 2
#' u = c(1.3, 1.8, 3)
#' E = runif(2*n, min=0, max=maxE)          # enrollment times
#' Z1 = (c(rep(0,n),rep(1,n)))                   # treatment indicator
#' Z2 = as.matrix(rnorm(2*n))                             # covariates
#' alpha = alpha0+alpha1*Z1
#' gamma1 = gamma0*exp(beta1*Z1+beta2*Z2)
#' FT = rweibull(2*n,shape=alpha,scale=gamma1**(-1/alpha))
#' CT = NULL
#' if (crate == 0)
#' {CT = Inf } else
#' {CT = rexp(2*n, rate=crate)}
#' Data = cbind(E,FT,CT,Z1,Z2)
#' delta = as.numeric(FT < CT)
#' X = pmin(FT,CT)
#' test <- gsRMST(tau = tau,Time = X,Status = delta,Z = Z2,TRT = Z1, Imax=1100, E = E, alpha = 0.05, u =u)
#' test$Z #1.051475 1.990320 2.232713
#' test$Crit # 2.818662 2.383526 2.022924
#' }
#'
#'
#'
#'
gsRMST = function(tau,Time,Status,Z,TRT, E, Imax, alpha = 0.05, u) {
  data = as.data.frame(cbind(TRT, E, Time, Status, Z))
  rmstdiff = rmstdiffse = NULL
  for (i in 1:length(u)) {
    data1 = data[ data$E < u[i],]
    E1 = data1$E
    Z1 = as.matrix(data1[, 5:dim(data)[2]])
    TRT1 = data1$TRT
    Time1 = data1$Time
    Status1 =  (data1$Status != 0)
    Status1 = Status1 & (Time1 < pmin(u[i] - E1, tau))
    X = pmin( Time1, u[i]-E1,tau)
    output <- rmst(tau = tau , X, Status1, Z1 , TRT1)
    rmstdiff[i] = output$muD
    rmstdiffse[i] = output$SED
  }
 Imax = 1/rmstdiffse[length(u)]^2
 print("Imax: ", Imax)
  IF = (1/rmstdiffse^2)/Imax

  Z = rmstdiff/rmstdiffse
  Crit =gsDesign(k = length(u), test.type=2, timing = IF, sfu = sfPower, sfupar = 3)$upper$bound
  return(list(Z = Z, Crit = Crit))

}



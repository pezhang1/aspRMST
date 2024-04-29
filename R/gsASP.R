#' Group sequential test of adjusted SPs
#'
#' @description Performs a group sequential test of adjusted SPs between two treatment groups
#' when provided the time point of analysis, the observed times, event indicator,
#' covariates, treatment group indicator, enrollment times, and calendar times of analysis.
#'
#'
#' @inherit Imaxasp references
#'
#' @param t0 -  Pre-specified time point of analysis
#' @param Time - Observed times
#' @param Status - Event indicator (0 = censored, 1 = observed)
#' @param Z - Non-treatment group covariates
#' @param TRT - Treatment group indicator (0 = Control, 1 = Treatment)
#' @param E - Enrollment times
#' @param alpha - Targeted type I error rate
#' @param u - Calendar times of analysis
#'
#'
#'
#'
#'
#'  @returns
#'  \itemize{
#'   \item Z - Vector of test statistics
#'   \item Crit - Vector of tritical values
#'    }
#' @export
#'
#' @examples
#' \dontrun{
#' library(survival)
#'
#' alpha0 = 1.5
#' alpha1 = -1
#' beta1 = -0.399782432
#' n=200
#' crate =0
#' beta2=0
#' t0=1
#' gamma0 = -log(0.4)
#' maxE = 2
#' u = c(1.6, 2.1, 3)
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
#' test <- gsASP(t0 = t0,Time = X,Status = delta,Z = Z2,TRT = Z1, E = E, alpha = 0.05, u =u)
#' test$Crit #2.802284 2.283383 2.037149
#' test$Z #2.545480 2.792856 2.628482
#' }
#'
gsASP = function(t0,Time,Status,Z,TRT, E, alpha = 0.05, u) {


  data = as.data.frame(cbind(TRT, E, Time, Status, Z))
  spd = spdse = NULL
  for (i in 1:length(u)) {
    data1 = data[ data$E < u[i],]
    E1 = data1$E
    Z1 = as.matrix(data1[, 5:dim(data)[2]])
    TRT1 = data1$TRT
    Time1 = data1$Time
    Status1 =  (data1$Status != 0)
    Status1 = Status1 & (Time1 < pmin(u[i] - E1, t0))
    X = pmin( Time1, u[i]-E1,t0)
    output <- asp(t0 = t0 , X, Status1, Z1 , TRT1)
    #output <- spse(t = t0 , Time = X, Status = status, Z=as.matrix(data1[, vars]) , D=data1$treat - 1)
    spd[i] = output$SPD
    spdse[i] = output$SED
  }
  Imax = 1/spdse[length(u)]^2
  IF = (1/spdse^2)/Imax
  Z = spd/spdse
  Crit =gsDesign(k = length(u), test.type=2, timing = IF, sfu = sfPower, sfupar = 3)$upper$bound
  return(list(Z = Z, Crit = Crit))
}




#' Group sequential test of adjusted restricted mean survival times
#'
#' @inherit Imaxrmst references
#'
#' @description Performs a group sequential test of adjusted restricted mean survival times (RMST)s between two treatment groups
#' when provided the time point of analysis, the observed times, event indicator,
#' covariates, treatment group indicator, enrollment times, and calendar times of analysis.
#'
#' @inheritParams gsASP
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
#' alpha0 = 1.5
#' alpha1 = -0.3
#' beta1 = -0.61335306
#' n=200
#' crate =0
#' beta2=0
#' t0=1
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
#' test <- gsRMST(t0 = t0,Time = X,Status = delta,Z = Z2,TRT = Z1, E = E, alpha = 0.05, u =u)
#' test$Z #1.051475 1.990320 2.232713
#' test$Crit # 2.818662 2.383526 2.022924
#' }
#'
#'
#'
#'
gsRMST = function(t0,Time,Status,Z,TRT, E, alpha = 0.05, u) {
  data = as.data.frame(cbind(TRT, E, Time, Status, Z))
  rmstdiff = rmstdiffse = NULL
  for (i in 1:length(u)) {
    data1 = data[ data$E < u[i],]
    E1 = data1$E
    Z1 = as.matrix(data1[, 5:dim(data)[2]])
    TRT1 = data1$TRT
    Time1 = data1$Time
    Status1 =  (data1$Status != 0)
    Status1 = Status1 & (Time1 < pmin(u[i] - E1, t0))
    X = pmin( Time1, u[i]-E1,t0)
    output <- rmst(t0 = t0 , X, Status1, Z1 , TRT1)
    #output <- spse(t = t0 , Time = X, Status = status, Z=as.matrix(data1[, vars]) , D=data1$treat - 1)
    rmstdiff[i] = output$muD
    rmstdiffse[i] = output$SED
  }
  Imax = 1/rmstdiffse[length(u)]^2
  IF = (1/rmstdiffse^2)/Imax
  Z = rmstdiff/rmstdiffse
  Crit =gsDesign(k = length(u), test.type=2, timing = IF, sfu = sfPower, sfupar = 3)$upper$bound
  return(list(Z = Z, Crit = Crit))

}



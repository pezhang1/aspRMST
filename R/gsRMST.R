


#' Title
#'
#' @param t0
#' @param Time
#' @param Status
#' @param Z
#' @param TRT
#' @param E
#' @param alpha
#' @param u
#'
#' @return
#' @export
#'
#' @examples
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
    rmstdiff[i] = output$Delta
    rmstdiffse[i] = output$DSE2
  }
  Imax = 1/rmstdiffse[length(u)]^2
  IF = (1/rmstdiffse^2)/Imax
  Z = rmstdiff/rmstdiffse
  Crit =gsDesign(k = length(u), test.type=2, timing = IF, sfu = sfPower, sfupar = 3)$upper$bound
  return(list(Z = Z, Crit = Crit))

}



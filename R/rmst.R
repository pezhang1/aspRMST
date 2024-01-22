#' rmst
#'
#' @description Calculates restricted mean survival time difference between two treatment groups
#' at a pre-specified time point via a treatment-stratified Cox proportional hazards model.
#' Variance and standard error estimates for the difference are also calculated.
#' Restricted mean survival time estimates and standard error / variance estimates
#' are also provided for each group.
#'
#' @author Peter Zhang, Brent Logan, Michael Martens
#'
#' @details
#'
#' @returns
#'  \itemize{
#'   \item Delta - restricted mean survival time differene estimate
#'   \item DSE2 - standard error estimate of restricted mean survival time difference estimate
#'   \item mu0 - restricted mean survival time of treatment group 0
#'   \item mu1 - restricted mean survival time of treatment group 1
#'   \item se10 - standard error of restricted mean survival time of treatment group 0
#'   \item se11 - standard error of restricted mean survival time of treatment group 1
#' }
#'
#' @references
#'
#' @seealso survival
#'
#' @param t0 - pre-specified time point, rmst is calculated over(0, t0)
#' @param Time - Observed times
#' @param Status - Censoring indicator (0 = Censored, 1 = Observed)
#' @param Z - Non-treatment group covariates
#' @param TRT - Treatment group indicator (0 = Control, 1 = Treatment)
#'
#' @return
#' @export
#'
#' @examples
#' library(survival)
#' t0 = 1
#' TRT = c(rep(0,5),rep(1,5))
#' Z = cbind(rnorm(10))
#' FT = c(0.30,  0.58,  0.41,  0.0333,  0.58,  0.10,  0.83,  2.45,  8.7, 17.1)
#' CT = c(0.20, 0.8, 0.73, 0.24, 0.5, 0.25, 1.3, 3.6, 10, 20.5)
#' X = pmin(FT,CT)
#' Status = as.numeric(FT <= CT)     & (X <= t0)
#' Time = pmin(X,t0)
#' output = rmst(t0,Time,Status,Z,TRT)
#' output$Delta
#' output$DSE2
#'
rmst = function(t0,Time,Status,Z,TRT){
  # Calculate RMST difference and its standard error


  SPSE = asp(t0,Time,Status,Z,TRT)
  n0 = SPSE$n0
  n1 = SPSE$n1
  n  = SPSE$n
  SigmaInv = SPSE$SigmaInv
  x0 = SPSE$X0
  x1 = SPSE$X1
  dx0 = diff(c(x0,t0),lag=1)
  dx1 = diff(c(x1,t0),lag=1)



  # Calculate RMST difference
  # Calculate B1i(u,t) and B2i(u,t) at t0.
  # Calculate Psi(iu,t) at t0.
  # Calculate B3(u,t) at t0.
  # Calculate mui(u|Z) at t0.
  sp0 = C01 = Gamma0 = Lambda0 = NULL
  C02 = QQ0 = matrix(0,length(x0),dim(Z)[2])
  cs0 = matrix(0,dim(Z)[1],length(x0))
  for(k in 1:length(x0))
  {
    spse0      = asp(x0[k],Time,Status,Z,TRT)
    sp0[k]     = spse0$S0
    C01[k]     = spse0$c01
    Gamma0[k]  = spse0$gamma0
    Lambda0[k] = spse0$Lambda00
    C02[k,]    = spse0$c02
    QQ0[k,]    = spse0$Q0
    cs0[,k]    = spse0$CS0
  }
  Temp0 = matrix(0,length(x0),length(x0))
  diag(Temp0) = Gamma0
  for(j in 1:(length(x0)-1))
  {
    for(k in (j+1):length(x0))
    {
      Temp0[j,k] = Temp0[k,j] = Gamma0[j]
    }
  }
  mu0 = sum(sp0*dx0)
  B10 = (n/n0)* rbind(C01*dx0) %*% Temp0 %*% cbind(C01*dx0)
  Psi0 = t(QQ0)%*%cbind(C01*dx0) - t(C02)%*%cbind(Lambda0*dx0)
  dgamma0 = c(0,diff(Gamma0,lag=1))
  B30 = t(Psi0) %*% SigmaInv %*% Psi0
  V10 = as.numeric( (B10+B30)/n )
  SE10 = sqrt(V10)
  cmu0 = as.numeric(cs0 %*% cbind(dx0))

  sp1 = C11 = Gamma1 = Lambda1 = NULL
  C12 = QQ1 = matrix(0,length(x1),dim(Z)[2])
  cs1 = matrix(0,dim(Z)[1],length(x1))
  for(k in 1:length(x1))
  {
    spse1      = asp(x1[k],Time,Status,Z,TRT)
    sp1[k]     = spse1$S1
    C11[k]     = spse1$c11
    Gamma1[k]  = spse1$gamma1
    Lambda1[k] = spse1$Lambda01
    C12[k,]    = spse1$c12
    QQ1[k,]    = spse1$Q1
    cs1[,k]    = spse1$CS1
  }
  Temp1 = matrix(0,length(x1),length(x1))
  diag(Temp1) = Gamma1
  for(j in 1:(length(x1)-1))
  {
    for(k in (j+1):length(x1))
    {
      Temp1[j,k] = Temp1[k,j] = Gamma1[j]
    }
  }
  mu1 = sum(sp1*dx1)
  B11  = (n/n1)* rbind(C11*dx1) %*% Temp1 %*% cbind(C11*dx1)
  Psi1 = t(QQ1)%*%cbind(C11*dx1) - t(C12)%*%cbind(Lambda1*dx1)
  dgamma1 = c(0,diff(Gamma1,lag=1))
  B31 = t(Psi1) %*% SigmaInv %*% Psi1
  V11 = as.numeric( (B11+B31)/n )
  SE11 = sqrt(V11)
  cmu1 = as.numeric(cs1 %*% cbind(dx1))

  Delta = mu1 - mu0
  B3 = t(Psi1-Psi0) %*% SigmaInv %*% (Psi1-Psi0)
  DV1 = as.numeric( (B10+B11+B3)/n )
  DSE1 = sqrt(DV1)



  # Calculate mean and variance of mu1(u|Z)-mu0(t|Z) at t0.
  v10 = V10 + var(cmu0)/n
  se10 = sqrt(v10)

  v11 = V11 + var(cmu1)/n
  se11 = sqrt(v11)

  DeltaV = var(cmu1-cmu0)
  DV2 = DV1 + DeltaV/n
  DSE2 = sqrt(DV2)

  list(mu0=mu0,mu1=mu1,Delta=Delta,DSE2=DSE2,
       se10=se10, se11=se11)
}








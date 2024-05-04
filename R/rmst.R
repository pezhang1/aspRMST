#' Calculate adjusted restricted mean survival times
#'
#'
#' @description Calculates RMST difference between two treatment groups
#' at a pre-specified time point via a treatment-stratified Cox proportional hazards model.
#' Standard error estimates for the estimated differences are also calculated.
#' RMST estimates and standard error estimates
#' are also provided for each group.
#'
#' @inherit Imaxrmst references
#'
#' @import stats
#' @importFrom survival coxph
#'
#'
#'
#'
#'
#'
#'
#' @inherit Imaxrmst references
#'
#'
#'
#'
#'
#' @param t0 Pre-specified time point, RMST estimate is calculated over(0, t0)
#' @param Time Observed times
#' @param Status Event indicator (0 = Censored, 1 = Observed)
#' @param Z Non-treatment group covariates
#' @param TRT Treatment group indicator (0 = Control, 1 = Treatment)
#'
#' @returns
#'  \itemize{
#'   \item muD Adjusted SP difference estimate
#'   \item SED Standard error estimate of adjusted SP difference estimate
#'   \item mu0 Adjusted SP of treatment group 0
#'   \item mu1 Adjusted SP of treatment group 1
#'   \item SE0 Standard error of adjusted SP of treatment group 0
#'   \item SE1 Standard error of adjusted SP of treatment group 1
#' }
#'
#'
#' @export
#'
#'
#'
#' @examples
#' set.seed(1234)
#' t0 = 1
#' n0 = 400
#' n1 = 400
#' n = n0 + n1
#' alpha0 = 1.5
#' alpha1 = -0.3
#' gamma0 = -log(0.4)
#' beta1 = -0.5
#' beta2 = log(1.5)
#' crate = -log(0.95)
#' TRT = c(rep(0,n0),rep(1,n1))                        # treatment indicator
#' Z = cbind(rnorm(n))                               # covariates
#' alpha = alpha0+alpha1*TRT
#' gamma1 = gamma0*exp(beta1*TRT+beta2*Z)
#' FT = rweibull(n,shape=alpha,scale=gamma1**(-1/alpha))
#' CT = rexp(n, rate=crate)
#' X = pmin(FT,CT)
#' Status = as.numeric(FT <= CT)     & (X <= t0)
#' Time = pmin(X,t0)
#' rmst(t0,Time,Status,Z,TRT)$muD #0.01355406
#' rmst(t0,Time,Status,Z,TRT)$SED # 0.04315677
#'
rmst = function(t0,Time,Status,Z,TRT){
  # Calculate RMST difference and its standard error


  SPSE = aspinternal(t0,Time,Status,Z,TRT)
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
    spse0      = aspinternal(x0[k],Time,Status,Z,TRT)
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
    spse1      = aspinternal(x1[k],Time,Status,Z,TRT)
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

  list(mu0=mu0,mu1=mu1,muD=Delta,SED=DSE2,
       SE0=se10, SE1=se11)
}




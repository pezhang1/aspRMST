#' Adjusted SP
#'
#' @inherit Imaxasp references
#' @description Calculates the adjusted SP difference between two treatment groups
#' at a pre-specified time point via a treatment-stratified Cox proportional hazards model.
#' Standard error estimates for the estimated difference are also calculated.
#' Adjusted SP estimates and standard error estimates
#' are also provided for each group.
#'
#'
#'
#'
#'
#' @import stats
#' @importFrom survival coxph
#'
#'
#'
#'
#'
#'
#' @param t0 - Pre-specified time point of analysis
#' @param Time - Observed times
#' @param Status - Event indicator (0 = censored, 1 = observed)
#' @param Z - Non-treatment group covariates
#' @param TRT - Treatment group indicator (0 = control, 1 = treatment)
#'
#' @returns
#'  \itemize{
#'   \item SPD - Adjusted SP difference estimate
#'   \item SED - Standard error estimate of adjusted SP difference estimate
#'   \item S0 - Adjusted SP estimate for treatment group 0
#'   \item S1 - Adjusted SP estimate for treatment group 1
#'   \item SE0 - Standard error estimate of adjusted SP estimate for treatment group 0
#'   \item SE1 - Standard error estimate of adjusted SP estimate for treatment group 1
#' }
#' @export
#'
#' @examples
#' set.seed(1234)
#' t0 = 1
#' TRT = c(rep(0,5),rep(1,5))
#' Z = cbind(rnorm(10))
#' FT = c(0.30,  0.58,  0.41,  0.0333,  0.58,  0.10,  0.83,  2.45,  8.7, 17.1)
#' CT = c(0.20, 0.8, 0.73, 0.24, 0.5, 0.25, 1.3, 3.6, 10, 20.5)
#' X = pmin(FT,CT)
#' Status = as.numeric(FT <= CT)     & (X <= t0)
#' Time = pmin(X,t0)
#' output = asp(t0,Time,Status,Z,TRT)
#' output$SPD
#' output$SED
asp = function(t0,Time,Status,Z,TRT) {
  n= length(TRT)
  n1 = sum(TRT)
  n0 = n-n1

  #run stratified Cox model
  Cox = coxph(Surv(Time,Status) ~ Z + strata(TRT))
  beta = cbind(as.numeric(coef(Cox)))
  SigmaInv = as.matrix(vcov(Cox))

  Time0 = Time[TRT==0]
  Time1 = Time[TRT==1]
  Status0 = Status[TRT==0]
  Status1 = Status[TRT==1]
  Z0 = as.matrix(subset(Z, subset= (TRT==0), drop = TRUE))
  Z1 = as.matrix(subset(Z, subset= (TRT==1), drop = TRUE))

  ZBETA  = as.numeric(Z%*%beta)
  ZBETA0 = as.numeric(Z0%*%beta)
  ZBETA1 = as.numeric(Z1%*%beta)


  X0 = unique(sort(Time0*Status0*(Time0<=t0)))
  X1 = unique(sort(Time1*Status1*(Time1<=t0)))


  temp1 = temp2 = NULL
  temp3 = temp4 = matrix(0,length(X0),dim(Z0)[2])
  for(k in 1:length(X0)) {
    NR0 = as.numeric(Time0>=X0[k]) #Number at risk for group0
    temp1[k] = sum(Status0[Time0==X0[k]])/sum(NR0*exp(ZBETA0))
    temp2[k]  = n0*sum(Status0[Time0==X0[k]])/(sum(NR0*exp(ZBETA0))^2)
    for(j in 1:dim(Z0)[2])
    {
      temp4[k,j] = sum(NR0*exp(ZBETA0)*Z0[,j])
    }
    temp3[k,] = sum(Status0[Time0==X0[k]])*temp4[k,]/(sum(NR0*exp(ZBETA0))^2)
  }
  Lambda00 = sum(temp1)
  gamma0 = sum(temp2)
  Q0 = cbind(apply(temp3,2,sum))


  temp5 = temp6 = NULL
  temp7 = temp8 = matrix(0,length(X1),dim(Z1)[2])
  for(k in 1:length(X1))
  {
    NR1 = as.numeric(Time1>=X1[k])
    temp5[k] = sum(Status1[Time1==X1[k]])/sum(NR1*exp(ZBETA1))
    temp6[k]  = n1*sum(Status1[Time1==X1[k]])/(sum(NR1*exp(ZBETA1))^2)
    for(j in 1:dim(Z1)[2])
    {
      temp8[k,j] = sum(NR1*exp(ZBETA1)*Z1[,j])
    }
    temp7[k,] = sum(Status1[Time1==X1[k]])*temp8[k,]/(sum(NR1*exp(ZBETA1))^2)
  }
  Lambda01 = sum(temp5)
  gamma1 = sum(temp6)
  Q1 = cbind(apply(temp7,2,sum))

  # Calculate marginal survival function estimate at t0
  # Calculate ci1(t) at t0
  # Calculate ci2(t) at t0
  temp1 = temp5 = temp2 = temp6 = NULL
  temp3 = temp7 = matrix(0,dim(Z)[1],dim(Z)[2])
  for(i in 1:dim(Z)[1])
  {
    ZBETAi = as.numeric(Z[i,]%*%beta)
    temp1[i] = exp(-exp(ZBETAi)*Lambda00)
    temp5[i] = exp(-exp(ZBETAi)*Lambda01)
    temp2[i]  = exp(-exp(ZBETAi)*Lambda00)*exp(ZBETAi)
    temp6[i]  = exp(-exp(ZBETAi)*Lambda01)*exp(ZBETAi)
    for(j in 1:dim(Z)[2])
    {
      temp3[i,j] = exp(-exp(ZBETAi)*Lambda00)*exp(ZBETAi)*Z[i,j]
      temp7[i,j] = exp(-exp(ZBETAi)*Lambda01)*exp(ZBETAi)*Z[i,j]
    }
  }
  CS0 = temp1
  CS1 = temp5
  S0  = mean(CS0)
  S1  = mean(CS1)
  SPD = mean(CS1-CS0)
  SPV = var(CS1-CS0)
  c01 = mean(temp2)
  c11 = mean(temp6)
  c02 = cbind(apply(temp3,2,mean))
  c12 = cbind(apply(temp7,2,mean))

  D0 = c01*Q0 - Lambda00*c02
  D1 = c11*Q1 - Lambda01*c12
  DD = D1-D0


  # Calculate standard error of SE difference
  V1 = c01^2*gamma0/n0 + c11^2*gamma1/n1 + t(DD)%*%SigmaInv%*%DD

  SE0 = c01^2*gamma0/n0
  SE1 = c11^2*gamma1/n1

  V2 = V1 + SPV/n
 # SE1 = sqrt(V1)
 # SED = sqrt(as.numeric(V2))

  return(list(S0=S0,S1=S1,SPD=SPD,SE0 = SE0 ,SE1=SE1,SED=as.numeric(sqrt(V2))))
}


#' Title
#'
#' @inheritParams asp
#'
#' @returns
#'  \itemize{
#'   \item SPD - Adjusted SP difference estimate
#'   \item SED - Standard error estimate of adjusted SP difference estimate
#'   \item S0 - Adjusted SP of treatment group 0
#'   \item S1 - Adjusted SP of treatment group 1
#'   \item se10 - Standard error of adjusted SP of treatment group 0
#'   \item se11 - Standard error of adjusted SP of treatment group 1
#' }
#' @export
#'
#' @keywords internal
#'
#'
aspinternal = function(t0,Time,Status,Z,TRT) {
  n= length(TRT)
  n1 = sum(TRT)
  n0 = n-n1

  #run stratified Cox model
  Cox = coxph(Surv(Time,Status) ~ Z + strata(TRT))
  beta = cbind(as.numeric(coef(Cox)))
  SigmaInv = as.matrix(vcov(Cox))

  Time0 = Time[TRT==0]
  Time1 = Time[TRT==1]
  Status0 = Status[TRT==0]
  Status1 = Status[TRT==1]
  Z0 = as.matrix(subset(Z, subset= (TRT==0), drop = TRUE))
  Z1 = as.matrix(subset(Z, subset= (TRT==1), drop = TRUE))

  ZBETA  = as.numeric(Z%*%beta)
  ZBETA0 = as.numeric(Z0%*%beta)
  ZBETA1 = as.numeric(Z1%*%beta)


  X0 = unique(sort(Time0*Status0*(Time0<=t0)))
  X1 = unique(sort(Time1*Status1*(Time1<=t0)))


  temp1 = temp2 = NULL
  temp3 = temp4 = matrix(0,length(X0),dim(Z0)[2])
  for(k in 1:length(X0)) {
    NR0 = as.numeric(Time0>=X0[k]) #Number at risk for group0
    temp1[k] = sum(Status0[Time0==X0[k]])/sum(NR0*exp(ZBETA0))
    temp2[k]  = n0*sum(Status0[Time0==X0[k]])/(sum(NR0*exp(ZBETA0))^2)
    for(j in 1:dim(Z0)[2])
    {
      temp4[k,j] = sum(NR0*exp(ZBETA0)*Z0[,j])
    }
    temp3[k,] = sum(Status0[Time0==X0[k]])*temp4[k,]/(sum(NR0*exp(ZBETA0))^2)
  }
  Lambda00 = sum(temp1)
  gamma0 = sum(temp2)
  Q0 = cbind(apply(temp3,2,sum))


  temp5 = temp6 = NULL
  temp7 = temp8 = matrix(0,length(X1),dim(Z1)[2])
  for(k in 1:length(X1))
  {
    NR1 = as.numeric(Time1>=X1[k])
    temp5[k] = sum(Status1[Time1==X1[k]])/sum(NR1*exp(ZBETA1))
    temp6[k]  = n1*sum(Status1[Time1==X1[k]])/(sum(NR1*exp(ZBETA1))^2)
    for(j in 1:dim(Z1)[2])
    {
      temp8[k,j] = sum(NR1*exp(ZBETA1)*Z1[,j])
    }
    temp7[k,] = sum(Status1[Time1==X1[k]])*temp8[k,]/(sum(NR1*exp(ZBETA1))^2)
  }
  Lambda01 = sum(temp5)
  gamma1 = sum(temp6)
  Q1 = cbind(apply(temp7,2,sum))

  # Calculate marginal survival function estimate at t0
  # Calculate ci1(t) at t0
  # Calculate ci2(t) at t0
  temp1 = temp5 = temp2 = temp6 = NULL
  temp3 = temp7 = matrix(0,dim(Z)[1],dim(Z)[2])
  for(i in 1:dim(Z)[1])
  {
    ZBETAi = as.numeric(Z[i,]%*%beta)
    temp1[i] = exp(-exp(ZBETAi)*Lambda00)
    temp5[i] = exp(-exp(ZBETAi)*Lambda01)
    temp2[i]  = exp(-exp(ZBETAi)*Lambda00)*exp(ZBETAi)
    temp6[i]  = exp(-exp(ZBETAi)*Lambda01)*exp(ZBETAi)
    for(j in 1:dim(Z)[2])
    {
      temp3[i,j] = exp(-exp(ZBETAi)*Lambda00)*exp(ZBETAi)*Z[i,j]
      temp7[i,j] = exp(-exp(ZBETAi)*Lambda01)*exp(ZBETAi)*Z[i,j]
    }
  }
  CS0 = temp1
  CS1 = temp5
  S0  = mean(CS0)
  S1  = mean(CS1)
  SPD = mean(CS1-CS0)
  SPV = var(CS1-CS0)
  c01 = mean(temp2)
  c11 = mean(temp6)
  c02 = cbind(apply(temp3,2,mean))
  c12 = cbind(apply(temp7,2,mean))

  D0 = c01*Q0 - Lambda00*c02
  D1 = c11*Q1 - Lambda01*c12
  DD = D1-D0


  # Calculate standard error of SE difference
  V1 = c01^2*gamma0/n0 + c11^2*gamma1/n1 + t(DD)%*%SigmaInv%*%DD

  #SE1 = c01^2*gamma0/n0
  #SE2 = c11^2*gamma1/n1

  V2 = V1 + SPV/n
  # SE1 = sqrt(V1)
  SE2 = sqrt(V2)

  return(list(Lambda00=Lambda00,Lambda01=Lambda01,S0=S0,S1=S1,SPD=SPD,SE2=SE2,CS0=CS0,CS1=CS1,
              gamma0=gamma0,gamma1=gamma1,Q0=Q0,Q1=Q1,c01=c01,c11=c11,c02=c02,c12=c12,D0=D0,D1=D1,DD=DD,SigmaInv=SigmaInv,
              X0=X0,X1=X1,n=n,n0=n0,n1=n1))
}

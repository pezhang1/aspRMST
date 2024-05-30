
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Introduction

Clinical trials often involve the analysis of time-to-event data to
determine the efficacy of a treatment. In the presence of baseline
covariates, the hazard ratio (HR) under the Cox proportional hazards
model is often used to quantify treatment effects. However, when the HR
is nonproportional over time, a Cox PH model may lead to incorrect
conclusions \#because the resulting estimates of the HR depend on the
follow-up duration and censoring distribution; the interpretation of the
HR estimate becomes less useful clinically in this context as well As an
alternative summary measure to the HR for describing the magnitude of
the treatment effect, the difference in adjusted survival probabilities
({}) at a fixed follow-up time has been proposed to quantify the
difference in survival outcomes between treatment arms. The SP measure
remains clinically interpretable even when survival curves have NPH,
facilitating communications between physicians and patients..

An alternative yet useful summary measure is the restricted mean
survival time (RMST) with time horizon $\tau$, also called $\tau$-year
mean survival time, which is the mean survival time of all subjects in
the study population followed up to the specific time point $\tau$. The
restricted time is often considered as information on patient outcomes
may not be available beyond a certain time point due to limited study
follow-up. The RMST summarizes the temporal event time distribution
during the time interval \[0, $\tau$\] and corresponds to the area under
the survival curve up to the survival time $\tau$. The RMST is an
alternative summary measure to the survival function and provides
information beyond the survival probability only at a single time point.

<!-- badges: start -->
<!-- badges: end -->

The goal of aspRMST is to calculate adjusted survival probability (SP)
estimates, adjusted restricted mean survival time (RMST) estimates, and
the estimated variances of these estimated differences. Sample size,
power, and effect size calculations are included for designing studies
that will test adjusted SP differences and adjusted RMST differences.
Finally, the package features two data analysis functions that conduct
group sequential tests of the adjusted SP difference and adjusted RMST
difference between two groups using error spending functions. When
provided the time point of analysis, observed times, event indicator,
covariates, treatment group indicator, enrollment times, and calendar
times of analyses, the functions will return the critical values and
test statistics.

## Installation

You can install aspRMST by downloading it from Github:
<https://github.com/pezhang1/aspRMST>

## Methods

The package allows the user to calculate the required sample size given
a specified power and targeted effect size, the power given a specified
sample size and targeted effect size, and the detectable effect size
given a specified power and sample size for both testing the difference
in adjusted SPs and difference in adjusted RMSTs between two groups for
a fixed sample setting. Trial data are simulated assuming the event time
for each subject follows a Weibull distribution with survival function

$$S(t|Z_W, \mathbf{Z}) = \exp(-\gamma t^{\alpha})$$,

where $Z_W = I(\mbox{group = treatment})$ is the treatment indicator and
{} are the covariates. The shape and rate parameters $\alpha$ and
$\gamma$ depend on $Z_W$ and $\mathbf{Z}$ through the forms
$\alpha = \alpha_0 + \alpha_1Z_W$ and
$\gamma= \gamma_0 \exp(\beta_WZ_W + \boldsymbol{\beta}^T \mathbf{Z})$
Note that when $\alpha_1 =0$, the hazard ratio is $\exp{(\beta_W)}$,
meaning that there are PH in this scenario; when $\alpha_1 \not=0$, the
PH assumption no longer holds. Functions are provided to perform sample
size, power, and effect size calculations for multiple binary
covariates, or for a single normal covariate.

In order to calculate the sample size, power, or targeted effect size,
the information for a trial under a fixed sample setting needs to be
calculated first. The effect size is the difference in marginal SPs for
adjusted SP comparisons and the difference in marginal RMSTs for
adjusted RMST comparisons, and depends on the Weibull model parameters ,
, and . For a given effect size $effect$, $\beta_W$ can be obtained
numerically by solving for the value of $\beta_W$ such that
$S_1(t_0)-S_0(t_0)=effect$ or $\mu_1(\tau)-\mu_0(\tau)=effect$. The
asymptotic distributions of the adjusted SP and RMST differences will be
used to perform sample size / power calculations. The package includes
functions to calculate $I^N$, the information for a fixed sample trial
with total sample size N via Monte Carlo simulation. Equal allocation
between the treatment and control group is used in this version of the
package. MC number of trials are simulated and the adjusted SP
difference or adjusted RMST difference between the treatment and control
group for each of the MC trials is calculated. Then, the inverse of the
variance of these differences is calculated to obtain $I^N$.

We now look at estimating the power given total sample size N and effect
size via Monte Carlo simulation. By using MC Monte Carlo replicates of
datasets of total sample size N, $I_{0}^N$ and $I_{effect}^N$ can be
calculated, which are the information for a trial with total sample size
$N$ under the null hypothesis with effect size 0 and under a targeted
alternative hypothesis where the effect size is $effect$, respectively.
This allows us to calculate $V_{0}^N = 1/I_{0}^N$ and
$V_{effect}^N=1/I_{effect}^N$, which are the variance of the estimated
difference under the respective null and targeted alternative hypotheses
for a fixed sample trial with total sample size N. Using the large
sample normal distribution of the effect size estimators, the power
$\pi$ can be expressed as
$$\pi = \Phi \left(\frac{effect - z_{1-\alpha/2} \sqrt{V_0^N}}{\sqrt{V_{effect}^N}}\right),$$
where $\alpha$ is the targeted type I error rate and $z_q=\Phi^{-1}(q)$
for any $q \in (0,1)$.

Given the sample size and power, we can also estimate the effect size.
Because the targeted effect size impacts the sample size formula through
the two terms $effect$ and $V_{effect}^N$, we use an iterative procedure
to obtain the required $effect$. We first calculate an initial variance
in the treatment group, $V_{10}^N$, by calculating $I_{0}^N$ and
initializing $V_{10}^N = 1/I_{0}^N$. Then, the initial effect size is
calculated as
$$effect_0 = z_{1-\alpha/2} \sqrt{V_0^N} + z_{\pi} \sqrt{V_{10}^N}.$$
Using this effect size, $V_{effect_0}^N$ is calculated by taking the
reciprocal of the information for a trial with effect size $effect_0$.
$V_{10}^N$ is then updated to $V_{effect_0}^N$. This updating of
$effect_0$ and $V_{10}^N$ is repeated until convergence occurs.

Consider estimation of the sample size given power and effect size via
Monte Carlo simulation. By using MC replicate datasets with an initial
sample size M, we can calculate $V_{0}^M=1/I_{0}^M$ and
$V_{effect}^M=1/I_{effect}^M$, which represent the variances of the
estimated differences in M patient trials. Under the large sample normal
distribution for the effect size estimator, the required sample size
depends on the unit variance, which is the reciprocal of the amount of
information contributed by a single patient. To obtain the unit
variances under the null and targeted alternative hypotheses, $V_{0}^1$
and $V_{effect}^1$, we rescale $V_{0}^M$ and $V_{effect}^M$ by
calculating $V_{0}^1 = V_{0}^M *M$ and $V_{effect}^1 = V_{effect}^M *M$.
With the individual variances per patient, we can calculate the actual
required total sample size N for the trial using the equation
$$N = \frac{\left(z_{1-\alpha/2}\sqrt{V_0^1} + z_{\pi} \sqrt{V_{effect}^1}\right)^2}{effect^2},$$

## Example1

This is an example to obtain restricted mean survival time estimates for
a single analysis using simulated data. Data are simulated under a
delayed effect setting with alpha0=1.5, alpha1=-0.3, and
gamma0=-log(0.4). The coefficient for the treatment indicator variable
is set to betaW=-0.5, and the coefficient for the covariate Z, which was
simulated from a standard normal distribution, is set to beta=log(1.5).
The censoring rate was set to exp(-log(0.95)) distribution to yield 5%
censoring per year. The sample size was set to 400 per group, and the
restriction window was set to \[0, 1\]. The estimated restricted mean
survival time difference is 0.0743 and the estimated standard error of
the estimated difference is 0.0217.

``` r
library(aspRMST)
#> Loading required package: survival
#> Loading required package: gsDesign
set.seed(1)
tau = 1
n0 = 400
n1 = 400
n = n0 + n1
alpha0 = 1.5
alpha1 = -0.3
gamma0 = -log(0.4)
betaW = -0.5
beta = log(1.5)
crate = -log(0.95)
TRT = c(rep(0,n0),rep(1,n1))                        # treatment indicator
Z = cbind(rnorm(n))                               # covariates
alpha = alpha0+alpha1*TRT
gamma1 = gamma0*exp(betaW*TRT+beta*Z)
FT = rweibull(n,shape=alpha,scale=gamma1**(-1/alpha))
CT = rexp(n, rate=crate)
X = pmin(FT,CT)
Status = as.numeric(FT <= CT)     & (X <= tau)
Time = pmin(X,tau)
rmst(tau,Time,Status,Z,TRT)$muD 
#> [1] 0.07427065
rmst(tau,Time,Status,Z,TRT)$SED 
#> [1] 0.02166575
```

## Example2

This is an example that performs a group sequential test of the equality
of RMSTs using simulated data. Data are simulated under a delayed effect
setting with alpha0=1.5, alpha1=-0.3, and gamma0=-log(0.4). The
coefficient for the treatment indicator variable is set to betaW=-0.613,
and the coefficient for the covariate Z, which was simulated from a
standard normal distribution, is set to beta=log(1.2). The censoring
rate was set to 0, i.e, no censoring. The sample size was set to 200 per
group, and the restriction window was set to \[0, 1\]. Enrollment was
assumed to be uniform between \[0, maxE\] = \[0, 2\]. Interim analyses
were performed at 1.3, 1.8, and 3 years. The maximum information was set
to 1100. The test statistics were 0.80, 2.56, 2.76; and the critical
values were 2.86, 2.43, and 2.02. A significant difference would be
found at stage 2.

``` r
library(aspRMST)
set.seed(1000)
alpha0 = 1.5
alpha1 = -0.3
betaW = -0.613
n=200
crate =0
beta=log(1.2)
tau=1
gamma0 = -log(0.4)
maxE = 2
u = c(1.3, 1.8, 3)
E = runif(2*n, min=0, max=maxE)          # enrollment times
Z1 = (c(rep(0,n),rep(1,n)))                   # treatment indicator
Z2 = as.matrix(rnorm(2*n))                             # covariates
alpha = alpha0+alpha1*Z1
gamma1 = gamma0*exp(betaW*Z1+beta*Z2)
FT = rweibull(2*n,shape=alpha,scale=gamma1**(-1/alpha))
CT = NULL
if (crate == 0)
 {CT = Inf } else
{CT = rexp(2*n, rate=crate)}
Data = cbind(E,FT,CT,Z1,Z2)
delta = as.numeric(FT < CT)
X = pmin(FT,CT)
test <- gsRMST(tau = tau, Time = X,Status = delta,Z = Z2,TRT = Z1, Imax=1100, E = E, alpha = 0.05, u =u)
test$TStat 
#> [1] 0.7030141 2.5626664 2.7628206
test$Crit 
#> [1] 2.856434 2.429439 2.016581
```

## Example3

This is an example that performs a sample size calculation for testing
equality of adjusted SPs given power and effect size. Data are simulated
under a crossing curves scenario using alpha0 = 1.5, alpha1=-1, and
gamma0=-log(0.4). Two binary covariates are used and simulated under a
Bern(0.5) and Bern(0.3) distributions. The coefficients for both
covariates are set to log(2)/sqrt(2). The censoring rate is set to 0,
i.e., no censoring. The pre-specified time point of analysis is 1 year.
Enrollment was assumed to be uniform between \[0, maxE\] = \[0, 2\].
Effect size was set to 0.121 difference in survival probabilities
between the treatment and control groups. Targeted type I error rate is
5% and targeted power is 80%. 10000 Monter Carlo iterations were used to
calculate the information for the trial. M was set to 1000 to obtain the
variances of the estimated differences in M patients. This variance was
rescaled to the unit variance to calculate the required sample size. See
the Methods section for more details on M. The calculated total sample
size was 385.

``` r
library(aspRMST)
set.seed(8)
Naspbinary(alpha0 = 1.5, alpha1=-1, gamma0=-log(0.4), beta=log(2)/sqrt(2), crate=0, t0=1,
maxE=2, M=1000, effect=0.121, MC=10000, alpha=0.05, pi = 0.8, p =c(0.5, 0.3)) 
#> [1] 384.5381
```

## Example4

This is an example that performs a power calculation for testing
equality of adjusted SPs given sample size and effect size. Data are
simulated under a crossing curves scenario using alpha0 = 1.5,
alpha1=-1, and gamma0=-log(0.4). Two binary covariates are used and
simulated under a Bern(0.5) and Bern(0.3) distributions. The
coefficients for both covariates are set to log(2)/sqrt(2). The
censoring rate is set to 0, i.e., no censoring. The pre-specified time
point of analysis is 1 year. Enrollment was assumed to be uniform
between \[0, maxE\] = \[0, 2\]. Effect size was set to 0.121 difference
in survival probabilities between the treatment and control groups.
Targeted type I error rate is 5% and targeted power is 80%. 10000 Monter
Carlo iterations were used to calculate the information for the trial.
The total sample size N was set to 386. The calculated power was 0.796.

``` r
library(aspRMST)
set.seed(16)
poweraspbinary(alpha0 = 1.5, alpha1=-1, gamma0=-log(0.4), beta=log(2)/sqrt(2), crate=0, t0=1,
maxE=2, N=386, effect=0.121, MC=10000, alpha = 0.05, p =c(0.5, 0.3)) 
#> [1] 0.7960318
```

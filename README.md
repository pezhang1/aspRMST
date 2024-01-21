
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Introduction

Clinical trials often involve the analysis of time-to-event data to
determine the efficacy of a treatment. In the presence of baseline
covariates, the hazard ratio (HR) under the Cox proportional hazards
model is often used to quantify treatment effects. However, when the HR
is nonproportional over time, a Cox PH model may lead to incorrect
conclusions; \#because the resulting estimates of the HR depend on the
follow-up duration and censoring distribution; the interpretation of the
HR estimate becomes less useful clinically in this context as well As an
alternative summary measure to the HR for describing the magnitude of
the treatment effect, the difference in survival probabilities ({}) at a
fixed follow-up time has been proposed to quantify the difference in
survival outcomes between treatment arms. The SP measure remains
clinically interpretable even when survival curves have NPH,
facilitating communications between physicians and patients..

An alternative yet useful summary measure is the restricted mean
survival time (RMST) at time *τ*, also called *τ*-year mean survival
time, which is the mean survival time of all subjects in the study
population followed up to the specific time point *τ*. The restricted
time is often considered, as information on patient outcomes may not be
available beyond a certain time point due to limited study follow-up.
The RMST summarizes the temporal event time distribution during the time
interval \[0, *τ*\] and corresponds to the area under the survival curve
up to the survival time *τ*. The RMST is an alternative summary measure
to the survival function and provides information beyond the survival
probability only at a single time point \# aspRMST

<!-- badges: start -->
<!-- badges: end -->

The goal of aspRMST is to …

## Installation

You can install the development version of aspRMST like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(aspRMST)
#> Loading required package: survival
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.

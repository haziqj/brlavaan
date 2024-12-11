
<!-- README.md is generated from README.Rmd. Please edit that file -->

# brlavaan

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check.yaml](https://github.com/haziqj/sem-bias/actions/workflows/R-CMD-check.yaml/badge.svg?branch=package)](https://github.com/haziqj/sem-bias/actions/workflows/R-CMD-check.yaml)
[![codecov](https://codecov.io/gh/haziqj/sem-bias/graph/badge.svg?token=00UGXV3BMK)](https://codecov.io/gh/haziqj/sem-bias)
<!-- badges: end -->

Apply empirical bias reduced methods to fit a variety of latent variable
models, including confirmatory factor analysis, structural equation
modelling, and latent growth curve models.

## Installation

You can install the development version of `{brlavaan}` from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("haziqj/sem-bias")
```

## Example

Here’s the beloved classic example – Bollen (1989)’s political democracy
example.

![](https://lavaan.ugent.be/figures/sem.png)

``` r
library(brlavaan)
#> Loading required package: lavaan
#> This is lavaan 0.6-19
#> lavaan is FREE software! Please report any bugs.
data("PoliticalDemocracy", package = "lavaan")
mod <- "
  # latent variables 
  ind60 =~ x1 + x2 + x3 
  dem60 =~ y1 + y2 + y3 + y4 
  dem65 =~ y5 + y6 + y7 + y8 
   
  # regressions
  dem60 ~ ind60 
  dem65 ~ ind60 + dem60 
  
  # residual covariances 
  y1 ~~ y5
  y2 ~~ y4 + y6 
  y3 ~~ y7 
  y4 ~~ y8
  y6 ~~ y8
"
fit <- brsem(model = mod, data = PoliticalDemocracy) 
summary(fit)
#> brlavaan 0.0.2.9007 ended normally after 71 iterations
#> 
#>   Estimator                                         ML
#>   Bias reduction method                       IMPLICIT
#>   Plugin penalty                         Huber penalty
#>   Optimization method                           NLMINB
#>   Number of model parameters                        31
#> 
#>   Number of observations                            75
#> 
#> 
#> Parameter Estimates:
#> 
#>   Standard errors                             Standard
#>   Information                                 Expected
#>   Information saturated (h1) model          Structured
#> 
#> Latent Variables:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   ind60 =~                                            
#>     x1                1.000                           
#>     x2                2.183    0.136   16.019    0.000
#>     x3                1.816    0.149   12.174    0.000
#>   dem60 =~                                            
#>     y1                1.000                           
#>     y2                1.255    0.179    7.001    0.000
#>     y3                1.058    0.149    7.121    0.000
#>     y4                1.263    0.142    8.876    0.000
#>   dem65 =~                                            
#>     y5                1.000                           
#>     y6                1.197    0.168    7.129    0.000
#>     y7                1.292    0.159    8.111    0.000
#>     y8                1.273    0.158    8.085    0.000
#> 
#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   dem60 ~                                             
#>     ind60             1.530    0.400    3.828    0.000
#>   dem65 ~                                             
#>     ind60             0.571    0.223    2.561    0.010
#>     dem60             0.830    0.097    8.547    0.000
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>  .y1 ~~                                               
#>    .y5                0.678    0.382    1.777    0.076
#>  .y2 ~~                                               
#>    .y4                1.429    0.758    1.886    0.059
#>    .y6                2.191    0.778    2.816    0.005
#>  .y3 ~~                                               
#>    .y7                0.838    0.647    1.294    0.196
#>  .y4 ~~                                               
#>    .y8                0.416    0.475    0.877    0.381
#>  .y6 ~~                                               
#>    .y8                1.436    0.612    2.347    0.019
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .x1                0.085    0.020    4.186    0.000
#>    .x2                0.127    0.073    1.730    0.084
#>    .x3                0.485    0.094    5.170    0.000
#>    .y1                1.977    0.470    4.204    0.000
#>    .y2                7.888    1.473    5.356    0.000
#>    .y3                5.393    1.015    5.314    0.000
#>    .y4                3.400    0.794    4.281    0.000
#>    .y5                2.544    0.519    4.906    0.000
#>    .y6                5.251    0.977    5.377    0.000
#>    .y7                3.629    0.762    4.764    0.000
#>    .y8                3.488    0.747    4.672    0.000
#>     ind60             0.484    0.093    5.196    0.000
#>    .dem60             4.297    0.990    4.340    0.000
#>    .dem65             0.222    0.231    0.961    0.336
```

By default, the implicit reduced bias ML estimator is used. To switch to
the *explicit* RBM, or to add a plugin penalty term, specify these as a
list to the `estimator.args` argument.

``` r
# for explicit RBM
fit <- brsem(model = mod, data = PoliticalDemocracy, 
             estimator.args = list(rbm = "explicit"))  
# for implicit RBM with plugin penalty
fit <- brsem(model = mod, data = PoliticalDemocracy, 
             estimator.args = list(rbm = "implicit", 
                                   plugin_penalty = brlavaan:::pen_huber))
```

To switch off the bias reduction, set `rbm = "none"`:

``` r
fit_ML  <- brsem(mod, PoliticalDemocracy, estimator.args = list(rbm = "none"))
fit_lav <-   sem(mod, PoliticalDemocracy)
tinytest::expect_equal(
   coef(fit_ML),
   coef(fit_lav),
   tolerance = 1e-4
)
#> ----- PASSED      : <-->
#>  call| tinytest::expect_equal(coef(fit_ML), coef(fit_lav), tolerance = 0.0001)
```

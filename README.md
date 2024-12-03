
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

You can install the development version of brlavaan from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("haziqj/sem-bias@package")
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
#> brlavaan 0.0.2.9001 ended normally after 65 iterations
#> 
#>   Estimator                                         ML
#>   Bias reduction method                       IMPLICIT
#>   Plugin penalty                             pen_ridge
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
#>     x2                2.181    0.136   16.050    0.000
#>     x3                1.814    0.149   12.191    0.000
#>   dem60 =~                                            
#>     y1                1.000                           
#>     y2                1.271    0.179    7.104    0.000
#>     y3                1.064    0.149    7.123    0.000
#>     y4                1.275    0.144    8.871    0.000
#>   dem65 =~                                            
#>     y5                1.000                           
#>     y6                1.202    0.168    7.171    0.000
#>     y7                1.297    0.160    8.117    0.000
#>     y8                1.279    0.158    8.097    0.000
#> 
#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   dem60 ~                                             
#>     ind60             1.518    0.395    3.846    0.000
#>   dem65 ~                                             
#>     ind60             0.564    0.222    2.542    0.011
#>     dem60             0.833    0.098    8.496    0.000
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>  .y1 ~~                                               
#>    .y5                0.686    0.378    1.813    0.070
#>  .y2 ~~                                               
#>    .y4                1.305    0.736    1.773    0.076
#>    .y6                2.045    0.751    2.724    0.006
#>  .y3 ~~                                               
#>    .y7                0.814    0.636    1.279    0.201
#>  .y4 ~~                                               
#>    .y8                0.410    0.468    0.876    0.381
#>  .y6 ~~                                               
#>    .y8                1.375    0.600    2.293    0.022
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .x1                0.085    0.020    4.180    0.000
#>    .x2                0.127    0.073    1.735    0.083
#>    .x3                0.485    0.094    5.171    0.000
#>    .y1                1.989    0.466    4.270    0.000
#>    .y2                7.474    1.412    5.295    0.000
#>    .y3                5.290    0.996    5.310    0.000
#>    .y4                3.314    0.782    4.239    0.000
#>    .y5                2.541    0.515    4.932    0.000
#>    .y6                5.093    0.953    5.344    0.000
#>    .y7                3.569    0.751    4.753    0.000
#>    .y8                3.419    0.735    4.650    0.000
#>     ind60             0.486    0.093    5.200    0.000
#>    .dem60             4.191    0.971    4.318    0.000
#>    .dem65             0.232    0.229    1.015    0.310
```

By default, the implicit reduced bias ML estimator (`iRBM`) is used. To
switch to the *explicit* RBM, or to add a plugin penalty term, specify
these as a list to the `estimator.args` argument.

``` r
# for explicit RBM
fit <- brsem(model = mod, data = PoliticalDemocracy, 
             estimator.args = list(rbm = "eRBM"))  
# for implicit RBM with plugin penalty
fit <- brsem(model = mod, data = PoliticalDemocracy, 
             estimator.args = list(rbm = "iRBM", plugin_penalty = "pen_ridge"))
```

To switch off the bias reduction, set `rbm = FALSE`:

``` r
fit_ML  <- brsem(mod, PoliticalDemocracy, estimator.args = list(rbm = FALSE))
fit_lav <-   sem(mod, PoliticalDemocracy)
tinytest::expect_equal(
   coef(fit_ML),
   coef(fit_lav),
   tolerance = 1e-4
)
#> ----- PASSED      : <-->
#>  call| tinytest::expect_equal(coef(fit_ML), coef(fit_lav), tolerance = 0.0001)
```

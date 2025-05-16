
<!-- README.md is generated from README.Rmd. Please edit that file -->

# brlavaan

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R CMD
check](https://github.com/haziqj/sem-bias/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/haziqj/sem-bias/actions/workflows/R-CMD-check.yaml)
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

## In brief

`{brlavaan}` provides explicit and implicit bias-reduction maximum
likelihood estimation of latent variable models that are `{lavaan}`
compatible. The main functions are

1.  `brsem()` for structural equation models.
2.  `brcfa()` for confirmatory factor analysis models.
3.  `brgrowth()` for growth curve models.
4.  `brlavaan()` for a more general interface to fit any model supported
    by `{lavaan}`.

Note the convention of `br` prefix for all functions in this package.
This is to distinguish between the bias-reduced methods and the standard
maximum likelihood methods provided by `{lavaan}`.

> \[!WARNING\]  
> Plugin penalties are not yet supported.

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

# lavaan fit (ML)
fit_lav <- sem(model = mod, data = PoliticalDemocracy)

# Bias-reduced fit (by default, implicit method is performed)
fit_iRBM <- brsem(model = mod, data = PoliticalDemocracy) 
summary(fit_iRBM)
#> brlavaan 0.1.1.9005 ended normally after 72 iterations
#> 
#>   Estimator                                         ML
#>   Bias reduction method                       IMPLICIT
#>   Plugin penalty                                  NONE
#>   Optimization method                           NLMINB
#>   Number of model parameters                        31
#> 
#>   Number of observations                            75
#> 
#> 
#> Parameter Estimates:
#> 
#>   Standard errors                             Standard
#>   Information                                 Observed
#>   Observed information based on                Hessian
#> 
#> Latent Variables:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   ind60 =~                                            
#>     x1                1.000                           
#>     x2                2.176    0.139   15.607    0.000
#>     x3                1.814    0.153   11.886    0.000
#>   dem60 =~                                            
#>     y1                1.000                           
#>     y2                1.254    0.185    6.791    0.000
#>     y3                1.050    0.148    7.093    0.000
#>     y4                1.254    0.150    8.374    0.000
#>   dem65 =~                                            
#>     y5                1.000                           
#>     y6                1.185    0.171    6.932    0.000
#>     y7                1.268    0.159    7.970    0.000
#>     y8                1.257    0.162    7.750    0.000
#> 
#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   dem60 ~                                             
#>     ind60             1.483    0.400    3.712    0.000
#>   dem65 ~                                             
#>     ind60             0.570    0.236    2.417    0.016
#>     dem60             0.832    0.098    8.510    0.000
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>  .y1 ~~                                               
#>    .y5                0.646    0.382    1.693    0.090
#>  .y2 ~~                                               
#>    .y4                1.292    0.707    1.828    0.068
#>    .y6                2.140    0.732    2.923    0.003
#>  .y3 ~~                                               
#>    .y7                0.816    0.640    1.275    0.202
#>  .y4 ~~                                               
#>    .y8                0.407    0.476    0.856    0.392
#>  .y6 ~~                                               
#>    .y8                1.312    0.572    2.294    0.022
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .x1                0.083    0.020    4.071    0.000
#>    .x2                0.122    0.072    1.699    0.089
#>    .x3                0.472    0.091    5.169    0.000
#>    .y1                1.896    0.481    3.943    0.000
#>    .y2                7.362    1.353    5.443    0.000
#>    .y3                5.161    1.004    5.143    0.000
#>    .y4                3.198    0.779    4.107    0.000
#>    .y5                2.406    0.513    4.694    0.000
#>    .y6                4.919    0.896    5.492    0.000
#>    .y7                3.474    0.748    4.642    0.000
#>    .y8                3.271    0.719    4.547    0.000
#>     ind60             0.447    0.086    5.197    0.000
#>    .dem60             4.019    0.958    4.197    0.000
#>    .dem65             0.203    0.231    0.880    0.379

# Should be different
tinytest::expect_equal(coef(fit_iRBM), coef(fit_lav))
#> ----- FAILED[data]: <-->
#>  call| tinytest::expect_equal(coef(fit_iRBM), coef(fit_lav))
#>  diff| Mean relative difference: 0.01134476
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
tinytest::expect_equal(
   coef(fit_ML),
   coef(fit_lav),
   tolerance = 1e-4
)
#> ----- PASSED      : <-->
#>  call| tinytest::expect_equal(coef(fit_ML), coef(fit_lav), tolerance = 1e-04)
```

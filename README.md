
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
fit <- brsem(model = mod, data = PoliticalDemocracy, information = "expected") 
summary(fit)
#> brlavaan 0.1.0 did NOT end normally after 18 iterations
#> ** WARNING ** Estimates below are most likely unreliable
#> 
#>   Bias reduction method                       IMPLICIT
#>   Plugin penalty                                  NONE
#>   Estimator                                         ML
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
#>     x2                2.193    0.299    7.328    0.000
#>     x3                1.821    0.261    6.984    0.000
#>   dem60 =~                                            
#>     y1                1.000                           
#>     y2                1.321    3.758    0.351    0.725
#>     y3                1.057    2.316    0.456    0.648
#>     y4                1.278    3.489    0.366    0.714
#>   dem65 =~                                            
#>     y5                1.000                           
#>     y6                1.314    1.926    0.682    0.495
#>     y7                1.430    1.652    0.866    0.387
#>     y8                1.419    2.007    0.707    0.480
#> 
#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   dem60 ~                                             
#>     ind60            -0.073    0.251   -0.292    0.771
#>   dem65 ~                                             
#>     ind60             0.075    0.225    0.333    0.739
#>     dem60             0.047    0.670    0.070    0.944
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>  .y1 ~~                                               
#>    .y5                0.110    0.420    0.261    0.794
#>  .y2 ~~                                               
#>    .y4                0.036    1.260    0.028    0.978
#>    .y6                0.044    0.812    0.054    0.957
#>  .y3 ~~                                               
#>    .y7                0.036    0.670    0.053    0.958
#>  .y4 ~~                                               
#>    .y8                0.046    0.680    0.067    0.946
#>  .y6 ~~                                               
#>    .y8                0.026    1.069    0.024    0.981
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .x1                0.024    0.032    0.750    0.453
#>    .x2                1.104    0.237    4.663    0.000
#>    .x3                0.930    0.185    5.023    0.000
#>    .y1                3.360    0.818    4.106    0.000
#>    .y2                7.691    1.845    4.168    0.000
#>    .y3                5.349    1.116    4.792    0.000
#>    .y4                5.466    1.457    3.752    0.000
#>    .y5                3.326    0.719    4.625    0.000
#>    .y6                5.553    1.369    4.056    0.000
#>    .y7                5.329    1.295    4.116    0.000
#>    .y8                5.174    1.395    3.710    0.000
#>     ind60             0.334    0.066    5.021    0.000
#>    .dem60             0.210    0.633    0.332    0.740
#>    .dem65             0.317    0.527    0.602    0.547
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
#>  call| tinytest::expect_equal(coef(fit_ML), coef(fit_lav), tolerance = 1e-04)
```

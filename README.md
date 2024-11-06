
<!-- README.md is generated from README.Rmd. Please edit that file -->

# brlavaan

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check](https://github.com/haziqj/sem-bias/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/haziqj/sem-bias/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/haziqj/sem-bias/graph/badge.svg)](https://app.codecov.io/gh/haziqj/sem-bias)
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
#> This is lavaan 0.6-18
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

# fit the model using explicit bias-reduced estimation
fit <- brsem(model = mod, data = PoliticalDemocracy, estimator = "eBRM") 

summary(fit)
#> brlavaan 0.0.1.9002 ended normally after 69 iterations
#> 
#>   Estimator                                       EBRM
#>   Penalty                                        FALSE
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
#>     x2                2.175    0.140   15.526    0.000
#>     x3                1.813    0.153   11.835    0.000
#>   dem60 =~                                            
#>     y1                1.000                           
#>     y2                1.253    0.182    6.893    0.000
#>     y3                1.049    0.151    6.922    0.000
#>     y4                1.252    0.145    8.657    0.000
#>   dem65 =~                                            
#>     y5                1.000                           
#>     y6                1.184    0.169    6.994    0.000
#>     y7                1.266    0.160    7.887    0.000
#>     y8                1.255    0.159    7.911    0.000
#> 
#> Regressions:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>   dem60 ~                                             
#>     ind60             1.484    0.403    3.682    0.000
#>   dem65 ~                                             
#>     ind60             0.570    0.225    2.537    0.011
#>     dem60             0.832    0.099    8.441    0.000
#> 
#> Covariances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>  .y1 ~~                                               
#>    .y5                0.648    0.366    1.771    0.077
#>  .y2 ~~                                               
#>    .y4                1.291    0.706    1.828    0.067
#>    .y6                2.140    0.737    2.902    0.004
#>  .y3 ~~                                               
#>    .y7                0.815    0.617    1.321    0.186
#>  .y4 ~~                                               
#>    .y8                0.409    0.450    0.908    0.364
#>  .y6 ~~                                               
#>    .y8                1.312    0.569    2.306    0.021
#> 
#> Variances:
#>                    Estimate  Std.Err  z-value  P(>|z|)
#>    .x1                0.083    0.020    4.179    0.000
#>    .x2                0.122    0.071    1.728    0.084
#>    .x3                0.472    0.091    5.173    0.000
#>    .y1                1.898    0.452    4.200    0.000
#>    .y2                7.363    1.379    5.338    0.000
#>    .y3                5.169    0.970    5.328    0.000
#>    .y4                3.200    0.749    4.271    0.000
#>    .y5                2.409    0.494    4.879    0.000
#>    .y6                4.919    0.915    5.374    0.000
#>    .y7                3.477    0.723    4.810    0.000
#>    .y8                3.273    0.701    4.670    0.000
#>     ind60             0.447    0.087    5.156    0.000
#>    .dem60             4.022    0.936    4.298    0.000
#>    .dem65             0.207    0.221    0.934    0.350
```

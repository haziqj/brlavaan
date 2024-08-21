
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Empirical bias reducing methods for Structural Equation Models

<!-- badges: start -->
<!-- badges: end -->

Last modified: 2024-08-21

## Two factor SEM toy example

Consider the following two factor SEM.

``` r
n <- 1000
dat <- simulateData(
  model = "
    eta1 =~ 1*y1 + 0.8*y2 + 0.6*y3
    eta2 =~ 1*y4 + 0.8*y5 + 0.6*y6
    eta2 ~ 0.3*eta1
  ",
  sample.nobs = n
)
head(dat)
#>           y1          y2          y3         y4         y5           y6
#> 1  0.6813519  0.09573879  0.22580113 -1.1693866  0.3164890  0.353604230
#> 2  0.4179522  0.36136792  1.14284107  1.1741100  0.3797073  2.173895157
#> 3 -0.3504405 -0.74416172 -0.19416321 -0.1794421  1.2365683 -0.909332334
#> 4 -0.8222861 -0.72477549 -0.51935492 -1.0326628  0.1856429 -1.665372696
#> 5 -2.4349577 -0.57915228 -1.97972924 -4.1888494 -0.3345901  0.377470995
#> 6 -0.3704730  2.18032880 -0.07861724 -1.3124645 -3.2150772 -0.009681741
```

![](README_files/figure-gfm/sempath-1.png)<!-- -->

Experiment: For each sample size `n` in `c(15, 20, 50, 100, 1000)`, we
simulate `B = 100` datasets and estimate the model parameters in two
ways. The first is maximum likelihood using `{lavaan}`. The second is
the implicit RBM with a penalty term. We then analyse the bias produced
by both methods.

``` r
tab_mab
#> # A tibble: 5 × 4
#>   n            ml  iRBMp   nok
#>   <fct>     <dbl>  <dbl> <dbl>
#> 1 15     826.     3.10     100
#> 2 20    2336.     1.93     100
#> 3 50      98.3    0.504    100
#> 4 100      0.213  0.213    100
#> 5 1000     0.0626 0.0624   100
```

``` r
tab_mab_by_type
#> # A tibble: 15 × 5
#>    n     type               ml  iRBMp   nok
#>    <fct> <chr>           <dbl>  <dbl> <dbl>
#>  1 15    Loadings      73.9    0.947    100
#>  2 15    Regressions   95.7    0.806    100
#>  3 15    Variances   1293.     4.47     100
#>  4 20    Loadings      30.7    0.830    100
#>  5 20    Regressions    5.13   0.644    100
#>  6 20    Variances   3779.     2.65     100
#>  7 50    Loadings       0.842  0.293    100
#>  8 50    Regressions    0.302  0.297    100
#>  9 50    Variances    159.     0.635    100
#> 10 100   Loadings       0.182  0.179    100
#> 11 100   Regressions    0.137  0.138    100
#> 12 100   Variances      0.239  0.239    100
#> 13 1000  Loadings       0.0502 0.0500   100
#> 14 1000  Regressions    0.0340 0.0340   100
#> 15 1000  Variances      0.0723 0.0721   100
```

``` r
p_bias
```

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->

Notes:

1.  Most of the bias of the ML parameter estimates is due to things
    going in wrong in the optimiser. The implicit RBM method is much
    more stable. If we set an arbitrary cut off for the ML (say consider
    only parameters \< 10 in absolute value) then the bias difference
    would be much smaller.

2.  I coded the objective function in a lazy way. Currently only works
    for this specific example. Ideally we dive into the `{lavaan}`
    source code and create a new objective function there…

3.  Noticed lots of numerical issues in the `optim()` routine. When the
    optimiser takes a weird turn in the parameter space things go
    haywire – the `JJJ()` and `HHH()` returns crazy values and the
    Hessian is not invertible. For now I cheated and used L-BFGS-B
    bounds to keep the optimiser in check.

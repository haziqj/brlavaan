
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Empirical bias reducing methods for Structural Equation Models

<!-- badges: start -->
<!-- badges: end -->

Last modified: 2024-08-31

## Toy example: Two factor SEM

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
#>            y1          y2         y3         y4         y5          y6
#> 1 -2.32655973 -3.17504296 -1.4907335 -1.8155530 -1.9055517 -0.19450996
#> 2 -1.31702987 -0.93324702 -1.3247469  0.3539171 -0.6759197 -1.15713914
#> 3  0.73643565  1.56472795 -1.7258649 -0.1202840  1.0968165 -0.09233550
#> 4 -0.01455502  0.08529215  0.6883789  0.1118176 -1.3122460  0.20347977
#> 5  3.42542923  2.87883536  1.9634456  0.5696859 -0.4487400  0.06150547
#> 6  0.64577916  0.18000289  0.7398690 -2.0939880 -2.7324616 -0.74960168
```

![](README_files/figure-gfm/sempath-1.png)<!-- -->

Experiment: For each sample size `n` in `c(15, 20, 50, 100, 1000)`, we
simulate `B = 1000` datasets and estimate the model parameters in four
ways.

1.  Maximum likelihood
2.  Explicit RBM
3.  Implicit RBM
4.  Implicit RBM with plugin penalty

Tables and graph below show results (mean bias).

    #> # A tibble: 15 × 10
    #>    n     type   bias_ML bias_eRBM bias_iRBM bias_iRBMp  mse_ML mse_eRBM mse_iRBM
    #>    <fct> <chr>    <dbl>     <dbl>     <dbl>      <dbl>   <dbl>    <dbl>    <dbl>
    #>  1 15    Load…  1.89    -9.31      0.219      0.212    3.57e+0  8.67e+1  4.78e-2
    #>  2 15    Regr…  0.718   -0.165    -0.201     -0.198    5.15e-1  2.72e-2  4.03e-2
    #>  3 15    Vari… -0.708    0.908    -0.178     -0.183    5.01e-1  8.24e-1  3.16e-2
    #>  4 20    Load…  1.65    -3.28      0.206      0.203    2.73e+0  1.08e+1  4.26e-2
    #>  5 20    Regr…  0.390   -3.45     -0.178     -0.179    1.52e-1  1.19e+1  3.16e-2
    #>  6 20    Vari… -0.517   -0.156    -0.151     -0.153    2.67e-1  2.44e-2  2.29e-2
    #>  7 50    Load…  0.220   -0.659     0.0850     0.0799   4.85e-2  4.34e-1  7.22e-3
    #>  8 50    Regr…  0.0536  -0.558    -0.0692    -0.0681   2.87e-3  3.12e-1  4.79e-3
    #>  9 50    Vari… -0.0924   0.144    -0.0532    -0.0554   8.53e-3  2.08e-2  2.83e-3
    #> 10 100   Load…  0.0349   0.0139    0.0433     0.0423   1.22e-3  1.94e-4  1.87e-3
    #> 11 100   Regr…  0.0109  -0.0419   -0.0111    -0.0120   1.20e-4  1.76e-3  1.22e-4
    #> 12 100   Vari… -0.0290  -0.0281   -0.0228    -0.0226   8.42e-4  7.91e-4  5.20e-4
    #> 13 1000  Load…  0.00128 -0.000779  0.00148    0.00148  1.64e-6  6.07e-7  2.19e-6
    #> 14 1000  Regr…  0.00160  0.000855  0.00273    0.00273  2.57e-6  7.32e-7  7.46e-6
    #> 15 1000  Vari… -0.00162 -0.000286 -0.000332  -0.000337 2.63e-6  8.20e-8  1.10e-7
    #> # ℹ 1 more variable: mse_iRBMp <dbl>

![](README_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

# Two factor SEM (D&R 2022)

# Growth model (D&R 2022)

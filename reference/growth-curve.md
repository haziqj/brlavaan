# Model specification for latent growth curve models

Model specification for latent growth curve models

## Usage

``` r
txt_mod_growth_pop(rel)

txt_mod_growth(rel)
```

## Arguments

- rel:

  Reliability of the growth curve model. Either 0.8 or 0.5.

## Value

A character string with the model specification, used for either fitting
lavaan models or simulating data.

## Examples

``` r
txt_mod_growth_pop(0.8)
#> [1] "\n    # intercept with coefficients fixed to 1\n    i =~  1*Day0 + 1*Day1 + 1*Day2 + 1*Day3 + 1*Day4 +\n          1*Day5 + 1*Day6 + 1*Day7 + 1*Day8 + 1*Day9\n\n    # slope with coefficients fixed to 0:9 (number of days)\n    s =~  0*Day0 + 1*Day1 + 2*Day2 + 3*Day3 + 4*Day4 +\n          5*Day5 + 6*Day6 + 7*Day7 + 8*Day8 + 9*Day9\n\n    i ~~ 550*i\n    i ~ 0*1\n\n    s ~~ 100*s\n    s ~ 0*1\n\n    i ~~ 40*s\n\n    Day0 ~~ 500*Day0\n    Day1 ~~ 500*Day1\n    Day2 ~~ 500*Day2\n    Day3 ~~ 500*Day3\n    Day4 ~~ 500*Day4\n    Day5 ~~ 500*Day5\n    Day6 ~~ 500*Day6\n    Day7 ~~ 500*Day7\n    Day8 ~~ 500*Day8\n    Day9 ~~ 500*Day9\n    "
txt_mod_growth(0.5)
#> [1] "\n  # intercept with coefficients fixed to 1\n  i =~  1*Day0 + 1*Day1 + 1*Day2 + 1*Day3 + 1*Day4 +\n        1*Day5 + 1*Day6 + 1*Day7 + 1*Day8 + 1*Day9\n\n  # slope with coefficients fixed to 0:9 (number of days)\n  s =~  0*Day0 + 1*Day1 + 2*Day2 + 3*Day3 + 4*Day4 +\n        5*Day5 + 6*Day6 + 7*Day7 + 8*Day8 + 9*Day9\n\n  i ~~ i\n  i ~ 1\n\n  s ~~ s\n  s ~ 1\n\n  i ~~ s\n\n  # fix intercepts\n  Day0 ~ 0*1\n  Day1 ~ 0*1\n  Day2 ~ 0*1\n  Day3 ~ 0*1\n  Day4 ~ 0*1\n  Day5 ~ 0*1\n  Day6 ~ 0*1\n  Day7 ~ 0*1\n  Day8 ~ 0*1\n  Day9 ~ 0*1\n\n  # apply equality constraints\n  Day0 ~~ v*Day0\n  Day1 ~~ v*Day1\n  Day2 ~~ v*Day2\n  Day3 ~~ v*Day3\n  Day4 ~~ v*Day4\n  Day5 ~~ v*Day5\n  Day6 ~~ v*Day6\n  Day7 ~~ v*Day7\n  Day8 ~~ v*Day8\n  Day9 ~~ v*Day9\n  "
```

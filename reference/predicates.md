# Predicates for `brlavaan` objects

Predicates for `brlavaan` objects

## Usage

``` r
is_ML(x)

is_eRBM(x)

is_iRBM(x)

is_iRBMp(x, quietly = FALSE)
```

## Arguments

- x:

  An object of class `brlavaan`, or an output from
  [`fit_sem()`](https://haziqj.ml/brlavaan/reference/fit_sem.md).

- quietly:

  If `TRUE`, suppresses messages.

## Value

`TRUE` or `FALSE`, or if not suppressed, a message indicating which
plugin penalty was used.

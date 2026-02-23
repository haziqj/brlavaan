# Functions for plugin penalty

Functions for plugin penalty

## Usage

``` r
pen_ridge(x, target = rep(0, length(x)), call = FALSE, ...)

pen_ridge_bound(x, lb, ub, call = FALSE, ...)

pen_huber(x, lb, ub, thres = 1, call = FALSE, ...)
```

## Arguments

- x:

  A numeric vector.

- target:

  The target values for components of `x` in the ridge penalty. Default
  is a vector of zeros.

- call:

  Logical. If `TRUE`, return the name of the penalty function.

- ...:

  Additional arguments which may be called. See Details.

- lb:

  The lower bounds for components of `x`.

- ub:

  The upper bounds for components of `x`.

- thres:

  The threshold for the Huber penalty.

## Value

A single numeric value, which is the penalty value.

# Fit two factor SEM using lavaan syntax manually

Fit two factor SEM using lavaan syntax manually

## Usage

``` r
fit_twofac(model, data, rbm = c("none", "explicit", "implicit"), start = NULL)
```

## Arguments

- model:

  A description of the user-specified model. Typically, the model is
  described using the lavaan model syntax. See
  [lavaan::model.syntax](https://rdrr.io/pkg/lavaan/man/model.syntax.html)
  for more information. Alternatively, a parameter table (eg. the output
  of the
  [`lavaan::lavaanify()`](https://rdrr.io/pkg/lavaan/man/model.syntax.html)
  function) is also accepted.

- data:

  An optional data frame containing the observed variables used in the
  model. If some variables are declared as ordered factors, lavaan will
  treat them as ordinal variables.

- rbm:

  The type of RBM method to use. One of `"none"`, `"explicit"`, or
  `"implicit"` (although, short forms are accepted, e.g. `"exp"`,
  `"iRBM`, etc.)

- start:

  A numeric vector of starting values

## Value

A list with the following components:

- coefficients:

  The estimated coefficients

- stderr:

  The standard errors of the estimated coefficients

- bias:

  The bias of the estimated coefficients

- timing:

  The elapsed time

- converged:

  Whether the optimization converged

- optim:

  The optimization object

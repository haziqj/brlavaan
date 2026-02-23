# Fit growth curve models using empirical bias-reducing methods

Fit growth curve models using empirical bias-reducing methods

## Usage

``` r
brgrowth(
  model,
  data,
  estimator = "ML",
  estimator.args = list(rbm = "implicit", plugin_pen = NULL),
  information = "observed",
  ...
)
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

- estimator:

  The estimator to use. **Currently only "ML" is supported.**

- estimator.args:

  A list containing RBM arguments. Possible arguments are

  - `rbm`: The type of RBM method to use. One of `"none"`, `"explicit"`,
    or `"implicit"` (although, short forms are accepted, e.g. `"exp"`,
    `"iRBM"`, etc.)

  - `plugin_pen`: The type of penalty to use. One of `NULL`,
    `"pen_ridge"`, or `"pen_ridge_bound"`.

  - `info_pen` The type of information matrix to use for the penalty
    term.

  - `info_bias` The type of information matrix to use for the bias term
    of the explicit reduced bias method.

- information:

  The type of information matrix to use for calculation of standard
  errors. Defaults to `"observed"`, although `"expected"` and
  `"first.order"` is also permitted.

- ...:

  Additional arguments to pass to the
  [lavaan::lavaan](https://rdrr.io/pkg/lavaan/man/lavaan.html) function.

## Value

An object of class `brlavaan` which is a subclass of the
[lavaan::lavaan](https://rdrr.io/pkg/lavaan/man/lavaan-class.html)
class.

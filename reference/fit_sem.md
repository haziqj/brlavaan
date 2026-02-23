# Fit a structural equation model using bias-reducing methods

This is the workhorse function to fit a structural equation model using
bias-reducing methods which we use for our simulations. For endusers, we
recommend using the
[`brsem()`](https://haziqj.ml/brlavaan/reference/brsem.md),
[`brcfa()`](https://haziqj.ml/brlavaan/reference/brcfa.md), or
[`brgrowth()`](https://haziqj.ml/brlavaan/reference/brgrowth.md)
functions which is a wrapper around this function and outputs something
more user-friendly.

## Usage

``` r
fit_sem(
  model,
  data,
  rbm = c("implicit", "explicit", "none"),
  plugin_pen = NULL,
  debug = FALSE,
  lavfun = "sem",
  maxgrad = FALSE,
  fn.scale = 1,
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

- rbm:

  The type of RBM method to use. One of `"none"`, `"explicit"`, or
  `"implicit"` (although, short forms are accepted, e.g. `"exp"`,
  `"iRBM`, etc.)

- plugin_pen:

  The type of penalty to use. One of `NULL`, `"pen_ridge"`, or
  `"pen_ridge_bound"`. User-specified penalty functions are possible.
  See description for details.

- debug:

  If TRUE, the function will return a list of intermediate results for
  debugging purposes.

- lavfun:

  The lavaan function to use. Default is "sem".

- maxgrad:

  If TRUE, the function will return the maximum gradient value.

- fn.scale:

  **\[experimental\]** A scaling factor for the log-likelihood function.

- ...:

  Additional arguments to pass to the
  [lavaan::lavaan](https://rdrr.io/pkg/lavaan/man/lavaan.html) function.

## Value

An object of class `brlavaan` which is a subclass of the
[lavaan::lavaan](https://rdrr.io/pkg/lavaan/man/lavaan-class.html)
class.

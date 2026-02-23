# Fit a latent variable model using bias-reducing methods

Fit a latent variable model using bias-reducing methods

## Usage

``` r
brlavaan(
  model,
  data,
  estimator = "ML",
  estimator.args = list(rbm = "implicit", plugin_pen = NULL),
  information = "observed",
  lavfun = "sem",
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

- lavfun:

  The lavaan function to use. Default is "sem".

- ...:

  Additional arguments to pass to the
  [lavaan::lavaan](https://rdrr.io/pkg/lavaan/man/lavaan.html) function.

## Value

An object of class `brlavaan` which is a subclass of the
[lavaan::lavaan](https://rdrr.io/pkg/lavaan/man/lavaan-class.html)
class.

## Details

The
[`pen_ridge()`](https://haziqj.ml/brlavaan/reference/plugin-penalties.md)
function applies a ridge regression style penalty \$f(x) = \|\| x
\|\|^2\$ that shrinks the parameters to zero. The
[`pen_ridge_bound()`](https://haziqj.ml/brlavaan/reference/plugin-penalties.md)
function applies a penalty that shrinks the parameters to the bounds of
the parameter space. The bounds are calculated by `{lavaan}` – see the
paper for more details.

The `info_pen` and `info_bias` arguments default to `"observed"`. It is
not recommended to change this to say `"expected"`, especially for small
sample sizes as the bias-reducing properties of the estimators are not
guaranteed.

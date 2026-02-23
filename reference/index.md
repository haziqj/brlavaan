# Package index

## Bias-reduced estimation

Main lavaan-style functions for bias-reduced estimation in SEM.

- [`brcfa()`](https://haziqj.ml/brlavaan/reference/brcfa.md) : Fit
  confirmatory factor analysis model using empirical bias-reducing
  methods
- [`brsem()`](https://haziqj.ml/brlavaan/reference/brsem.md) : Fit a
  structural equation model using empirical bias-reducing methods
- [`brgrowth()`](https://haziqj.ml/brlavaan/reference/brgrowth.md) : Fit
  growth curve models using empirical bias-reducing methods
- [`brlavaan()`](https://haziqj.ml/brlavaan/reference/brlavaan.md) : Fit
  a latent variable model using bias-reducing methods
- [`brlavaan-class`](https://haziqj.ml/brlavaan/reference/brlavaan-class.md)
  : brlavaan Class

## Penalty functions

Penalty functions for iRBM \[*For future implementation*\].

- [`pen_ridge()`](https://haziqj.ml/brlavaan/reference/plugin-penalties.md)
  [`pen_ridge_bound()`](https://haziqj.ml/brlavaan/reference/plugin-penalties.md)
  [`pen_huber()`](https://haziqj.ml/brlavaan/reference/plugin-penalties.md)
  : Functions for plugin penalty

## Manual eBR and iBR fit functions

For internal testing purposes.

- [`fit_sem()`](https://haziqj.ml/brlavaan/reference/fit_sem.md) : Fit a
  structural equation model using bias-reducing methods
- [`fit_twofac()`](https://haziqj.ml/brlavaan/reference/fit_twofac.md) :
  Fit two factor SEM using lavaan syntax manually
- [`fit_growth()`](https://haziqj.ml/brlavaan/reference/fit_growth.md) :
  Fit growth curve models using lavaan syntax manually

## Helper functions for simulation studies

Functions to facilitate simulation study in our paper.

- [`gen_data_growth()`](https://haziqj.ml/brlavaan/reference/gen-data.md)
  [`gen_data_twofac()`](https://haziqj.ml/brlavaan/reference/gen-data.md)
  : Generate data for simulation studies

- [`truth_growth()`](https://haziqj.ml/brlavaan/reference/true-values.md)
  [`truth_twofac()`](https://haziqj.ml/brlavaan/reference/true-values.md)
  [`truth()`](https://haziqj.ml/brlavaan/reference/true-values.md) :
  Retrieve the true values of the models used for simulation

- [`txt_mod_growth_pop()`](https://haziqj.ml/brlavaan/reference/growth-curve.md)
  [`txt_mod_growth()`](https://haziqj.ml/brlavaan/reference/growth-curve.md)
  : Model specification for latent growth curve models

- [`txt_mod_twofac_pop()`](https://haziqj.ml/brlavaan/reference/two-fac.md)
  [`txt_mod_twofac()`](https://haziqj.ml/brlavaan/reference/two-fac.md)
  : Model specification for two-factor SEM models

- [`is_ML()`](https://haziqj.ml/brlavaan/reference/predicates.md)
  [`is_eRBM()`](https://haziqj.ml/brlavaan/reference/predicates.md)
  [`is_iRBM()`](https://haziqj.ml/brlavaan/reference/predicates.md)
  [`is_iRBMp()`](https://haziqj.ml/brlavaan/reference/predicates.md) :

  Predicates for `brlavaan` objects

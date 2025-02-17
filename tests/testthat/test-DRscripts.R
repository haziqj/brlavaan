skip()

## ----- Growth curve models ---------------------------------------------------
dat <- gen_data_growth(n = 15, rel = 0.8, dist = "Normal", scale = 1)
mod <- txt_mod_growth(0.8)
tru <- truth(dat)
fit <- growth(mod, dat, start = tru, bounds = "standard", meanstructure = TRUE)

# Jackknife Bounds: bias-corrected estimates (JB)
jack.out <- sd_jack_bc(fit = fit)

# Bootstrap Bounds: bias-corrected estimates (BB) - ci.norm
boot.out <- sd_boot_bc(fit = fit, ndraws = 10L)

# Ozenne bias-corrected estimates (Ozenne)
ozenne.out <- try(yr_ozenne_bcs_growth(fit, return.se = TRUE), silent = TRUE)

# REML
out.REML <- dr_reml_growth(dat)

## ----- Two-factor SEM --------------------------------------------------------
dat <- gen_data_twofac(n = 15, rel = 0.5, dist = "Normal", meanstructure = TRUE)
mod <- txt_mod_twofac(0.5)
tru <- truth(dat)
fit <- sem(mod, dat, start = tru, bounds = "standard", meanstructure = TRUE)

# Jackknife Bounds: bias-corrected estimates (JB)
jack.out <- sd_jack_bc(fit = fit)

# Bootstrap Bounds: bias-corrected estimates (BB)
boot.out <- sd_boot_bc(fit = fit, ndraws = 10L)

# Ozenne bias-corrected estimates (Ozenne)
ozenne.out <- try(yr_ozenne_bcs_twofac(fit, return.se = TRUE), silent = TRUE)

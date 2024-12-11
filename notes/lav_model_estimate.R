# lavaan:::lav_model_estimate
# This function has the pack and unpack parameters (look at the objective_function() and gradient_function() functions)

library(tidyverse)
library(lavaan)
library(tinytest)

HS.model <- "
  visual  =~ x1 + a*x2 + x3
  textual =~ x4 + a*x5 + x6
  speed   =~ x7 + a*x8 + x9
"

fit0 <- cfa(HS.model, data = HolzingerSwineford1939, do.fit = FALSE)
theta.init <- coef(fit0)

# extract some internal slots for housekeeping
lavsamplestats <- fit0@SampleStats
lavmodel       <- fit0@Model
lavdata        <- fit0@Data
lavoptions     <- fit0@Options
lavpartable    <- fit0@ParTable
lavcache       <- fit0@Cache

tmp <- lavaan:::lav_model_estimate(
    lavmodel = lavmodel,
    lavpartable = lavpartable,
    lavh1 = NULL,
    lavsamplestats = lavsamplestats,
    lavdata = lavdata,
    lavoptions = lavoptions,
    lavcache = lavcache,
    start = "simple",
    do.fit = TRUE
  )

# "unpack" to long version of the parameters
x.unpack <- lav_model_get_parameters(lavmodel)
parscale <- rep(1, length(x.unpack))
z.unpack <- x.unpack * parscale

# "pack" to make the short version of the parameters
if (lavmodel@eq.constraints) {
  z.pack <- as.numeric(
    (z.unpack - lavmodel@eq.constraints.k0) %*% lavmodel@eq.constraints.K
  )
} else {
  z.pack <- z.unpack
}
start.x <- z.pack

# and for gradients:
# if (lavmodel@eq.constraints) {
#   dx <- as.numeric(dx %*% lavmodel@eq.constraints.K)
# }
#
# and I guess for hessians:
# if (lavmodel@eq.constraints) {
#   H <- t(lavmodel@eq.constraints.K) %*% H %*% lavmodel@eq.constraints.K
# }






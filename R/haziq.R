# From Yves  Rosseel (17/11/2023)
library(lavaan)

HS.model <- ' visual  =~ x1 + x2 + x3
              textual =~ x4 + x5 + x6
              speed   =~ x7 + x8 + x9 '

# ad hoc values for the 21 parameter values
# (in the order of the parameter table; see parTable(fit))
custom.values <- c(0.6, 0.7, 1.1, 0.9, 1.2, 1.1, 0.5, 1.1, 0.8, 0.4, 0.4, 0.4, 0.8, 0.5, 0.6, 0.8, 1, 0.4, 0.4, 0.3, 0.2)

fit <- cfa(HS.model, data = HolzingerSwineford1939, start = custom.values,
           optim.method = "none", baseline = FALSE, se = "robust.huber.white")
J <- lavInspect(fit, "information.firstorder")
# or
J <- lavTech(fit, "information.firstorder")

# using the lower-level function lav_model_information_firstorder()
fit2 <- cfa(HS.model, data = HolzingerSwineford1939, do.fit = FALSE)
# extract slots
lavmodel <- fit2@Model
lavsamplestats <- fit2@SampleStats
lavdata <- fit2@Data
lavoptions <- fit2@Options

# update lavmodel with 'new' set of values
lavmodel2 <- lav_model_set_parameters(lavmodel, x = custom.values)

J2 <- lavaan:::lav_model_information_firstorder(lavmodel = lavmodel2,
         lavsamplestats = lavsamplestats, lavdata = lavdata,
         lavoptions = lavoptions)

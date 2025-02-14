library(lavaan)

model <- "
  i =~ 1*t1 + 1*t2 + 1*t3 + 1*t4
  s =~ 0*t1 + 1*t2 + 2*t3 + 3*t4

  i~~i
  s~~s
  i~~s
  i~1
  s~1

  t1 ~~ v*t1
  t2 ~~ v*t2
  t3 ~~ v*t3
  t4 ~~ v*t4
"

# Equality constraints (no bounds)
fit <- growth(model, data = Demo.growth[1:10, ])
coef(fit)
fit@Model@eq.constraints  # TRUE
fit@Model@eq.constraints.K
fit@Model@eq.constraints.k0

(x <- coef(fit))
x_pack <- as.numeric(
  (x - fit@Model@eq.constraints.k0) %*% fit@Model@eq.constraints.K
)
x_pack

# Equality constraints (with bounds)
fit <- growth(model, data = Demo.growth[1:10, ], bounds = "standard",
              ceq.simple = TRUE)
coef(fit)
fit@Model@eq.constraints  # FALSE
fit@Model@eq.constraints.K # empty
fit@Model@eq.constraints.k0 # empty

x <- coef(fit)
x_pack <- x



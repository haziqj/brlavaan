library(tidyverse)
library(brlavaan)
dat <- readxl::read_xlsx("inst/influencer/raw SMI data.xlsx")
mycols <- c(
  ML = "#E31A1C",
  lav = "#FB9A99",
  eRBM = "#A6CEE3",
  iRBM = "#1F78B4",
  Jackknife = "#B2DF8A",
  Bootstrap = "#33A02C",
  `Ozenne et al.` = "#FDBF6F",
  REML = "#FF7F00"
)

# Rename covariates
dat <-
  dat |>
  select(-1) |>
  rename(
    gender     = `Gender`,
    age        = `Age`,
    hrs_socmed = `How many hours do you spend on social media platform per day?`,
    platforms  = `Social media platforms that you use (you can choose more than one)`
  ) |>
  mutate(
    gender     = factor(gender),
    age        = factor(age),
    hrs_socmed = factor(hrs_socmed)
  )

# Rename items and codify 1--5 scale
item_names <- c(
  paste0("phy", 1:5),   # Physical
  paste0("tru", 1:5),   # Trustworthiness
  paste0("exp", 1:5),   # Expertise
  paste0("coo", 1:4)    # Coolness
)
dat <-
  dat |>
  rename_with(~ item_names, .cols = all_of(names(dat)[5:23])) |>
  mutate(across(phy1:coo4, function(x) {
    as.numeric(factor(x, levels = c(
      "Completely disagree",
      "Mostly disagree",
      "Slightly disagree",
      "Neither disagree nor agree",
      "Slightly agree",
      "Mostly agree",
      "Completely agree"
    )))})
  )

# Remove missing values
dat <- dat |>
  drop_na() |>
  filter(gender != "Prefer not to say") |>
  filter(age != "Less than 18 years old")

# Model fit
mod <- '
  # Measurement model
  exprt =~ exp1 + exp2 + exp3 + exp4 + exp5
  attrc =~ phy1 + phy2 + phy3 + phy4 + phy5
  trust =~ tru1 + tru2 + tru3 + tru4 + tru5
  cool  =~ coo1 + coo2 + coo3 + coo4

  # Structural model
  trust ~ upper(10)*exprt + attrc + cool
'

# First, the overall fit
fit1 <- sem(mod, dat, std.lv = TRUE, information = "observed")
fit2 <- fit_sem(mod, dat, std.lv = TRUE, rbm = "explicit", verbose = TRUE,
                start = coef(fit1))
fit3 <- fit_sem(mod, dat, std.lv = TRUE, rbm = "implicit", verbose = TRUE,
                start = coef(fit1))

# Next, the multigroup SEM
fit_lav  <-   sem(mod, dat, std.lv = TRUE, group = "gender", information = "observed")
fit_eRBM <- fit_sem(mod, dat, std.lv = TRUE, group = "gender", verbose = TRUE,
                    start = coef(fit_lav), rbm = "explicit")
fit_iRBM <- fit_sem(mod, dat, std.lv = TRUE, group = "gender", verbose = TRUE,
                    start = coef(fit_lav), rbm = "implicit",
                    control = list(trace = 1, iter.max = 20, eval.max = 40))

# See results
list(
  MLovr = fit1, eRBMovr = fit2, iRBMovr = fit3,
  ML = fit_lav, eRBM = fit_eRBM, iRBM = fit_iRBM
) |>
  map(\(x) {

    get_se <- function(z) {
      if (inherits(z, "lavaan")) {
        res <- z@ParTable$se
        res <- res[z@ParTable$free > 0]
      } else {
        res <- z$stderr
      }
      res
    }

    tibble(
      param = names(coef(x)),
      est = as.numeric(coef(x)),
      se = get_se(x),
    )
  }) |>
  bind_rows(.id = "method") |>
  filter(grepl("[a-z]~[a-z]", param)) |>
  mutate(
    group = case_when(
      grepl("g2", param) ~ "Male",
      grepl("ovr", method) ~ "Overall",
      TRUE ~ "Female"
    ),
    group = factor(group, levels = rev(c("Overall", "Female", "Male"))),
    param = gsub(".g2", "", param),
    method = gsub("ovr", "", method),
    signif = abs(est) > 1.96 * se
  ) |>
  ggplot(aes(x = est, y = group, col = method)) +
  geom_pointrange(aes(xmin = est - 1.96 * se,
                      xmax = est + 1.96 * se),
                  position = position_dodge(width = 0.3),
                  size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  guides(fill = guide_legend(rev = TRUE), col = guide_legend(rev = TRUE)) +
  coord_cartesian(xlim = c(-10, 15)) +
  facet_grid(. ~ param, scales = "free_x") +
  # scale_fill_manual(values = mycols) +
  scale_colour_manual(values = mycols) +
  # scale_alpha_manual(values = c(0.2, 1)) +
  theme_bw() +
  labs(col = NULL)


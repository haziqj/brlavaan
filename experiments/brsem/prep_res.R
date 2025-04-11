library(tidyverse)
library(brlavaan)
library(glue)
library(gt)
here::i_am("experiments/brsem/prep_res.R")
if (file.exists(here::here("experiments/brsem/results.RData")))
  stop("Results already exist. Delete the file if you want to rerun the analysis.")

twofacpars <- c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")
growthpars <- c("v", "i~~i", "s~~s", "i~~s")

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

## ----- Download --------------------------------------------------------------

# Data files
data_files <- paste0("experiments/brsem/", c(
  "simu_id.RData",
  "simu_res_twofac_mp1.RData",
  "simu_res_twofac_mp2.RData",
  "simu_res_twofac_mp3.RData",
  "simu_res_growth.RData"
))
if (any(!file.exists(here::here(data_files)))) {
  cat("Downloading data files...\n")
  zipfile <- here::here("experiments/brsem/data.zip")
  download.file(
    "https://files.au-1.osf.io/v1/resources/z8aky/providers/osfstorage/?zip=",
    destfile = zipfile
  )
  unzip(zipfile, exdir = here::here("experiments/brsem"))
  file.remove(zipfile)
}

load(here::here("experiments/brsem/simu_id.RData"))
simu_id <- bind_rows(
  mutate(simu_id, model = "twofac", .after = simid),
  mutate(simu_id, model = "growth", .after = simid)
)

load(here::here("experiments/brsem/simu_res_twofac_mp1.RData"))
simu_res_twofac1 <- simu_res_twofac  # 20-30
load(here::here("experiments/brsem/simu_res_twofac_mp2.RData"))
simu_res_twofac2 <- simu_res_twofac  # 11-19
load(here::here("experiments/brsem/simu_res_twofac_mp3.RData"))
simu_res_twofac[11:19] <- simu_res_twofac2[11:19]
simu_res_twofac[20:30] <- simu_res_twofac1[20:30]
load(here::here("experiments/brsem/simu_res_growth.RData"))
simu_res <- c(simu_res_twofac, simu_res_growth)

## ----- Convergence statistics ------------------------------------------------
res_ours_nested <-
  simu_res |>
  bind_rows() |>
  filter(method %in% c("ML", "eRBM", "iRBM")) |>
  mutate(
    super_converged = converged & Sigma_OK &
      mapply(function(se_vec, mod) {
        if(any(is.na(se_vec))) return(FALSE)
        if(mod == "twofac") {
         return(all(se_vec < 5))
        } else if(mod == "growth") {
         return(all(se_vec < 500))
        } else {
         return(TRUE)
        }
      }, se, model)
  )

res_conv <-
  res_ours_nested |>
  summarise(count = sum(super_converged, na.rm = TRUE), .by = dist:method) |>
  mutate(count = count / max(count), .by = dist:n) |>
  pivot_wider(names_from = c(model, method, rel), values_from = count) |>
  group_by(dist)

## ----- Timing statistics -----------------------------------------------------
res_timing_n1000 <-
  simu_res |>
  bind_rows() |>
  filter(method != "lav", n == 1000) |>
  summarise(
    mean = mean(timing),
    sd = sd(timing, na.rm = TRUE),
    .by = c(model, method)
  )

## ----- Big results data frame ------------------------------------------------
res <-
  simu_res |>
  bind_rows() |>
  mutate(param = lapply(truth, names), .before = est) |>
  unnest(param:truth) |>
  distinct(sim, dist, model, rel, n, method, param, .keep_all = TRUE) |>
  mutate(
    type = case_when(
      grepl("f[xy]=~[xy][0-9]", param) ~ "Lambda",
      grepl("fy~fx", param) ~ "beta",
      grepl("[xy][0-9]~~[xy][0-9]|v", param) ~ "Theta",
      grepl("f[xy]~~f[xy]|[is]~~[is]", param) ~ "Psi",
      grepl("[is]~1", param) ~ "alpha",
      grepl("[xy][0-9]~1", param) ~ "nu",
      TRUE ~ NA
    ),
    across(c(est, truth, se), ~if_else(model == "growth" & type != "alpha", .x * 100, .x)),  # rescale back
    dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
    rel = factor(rel, levels = c("0.8", "0.5"), labels = c("Rel = 0.8", "Rel = 0.5")),
    method = factor(
      method,
      levels = c("ML", "lav", "eRBM", "iRBM", "JB", "BB", "Ozenne", "REML"),
      labels = names(mycols)
    ),
    bias = est - truth,
    relbias = bias / truth,
    covered = truth <= est + qnorm(0.975) * se & truth >= est - qnorm(0.975) * se
  )

# Load D&R two factor sims
dr_data_file <- here::here("experiments/brsem/simu_res_twofacDR.RData")
if (!file.exists(dr_data_file)) {
  source(here::here("experiments/brsem/_prepDR.R"))
} else {
  load(dr_data_file)
}

# For two-factor model BB & JB, get from res_dr
res <-
  res |>
  filter(!(method %in% c("Bootstrap", "Jackknife") & model == "twofac")) |>
  bind_rows(
    filter(res_dr, model == "twofac", method %in%  c("Bootstrap", "Jackknife"))
  )

## ----- Prepare data for plots ------------------------------------------------

# For the plots showing performance of our methods
plot_df <-
  res |>
  filter(method %in% c("ML", "eRBM", "iRBM")) |>
  filter( converged, !is.na(se)) |>
  # for each kind of model, filter bad standard errors
  filter(!(model == "twofac" & abs(se) > 5)) |>
  filter(!(model == "growth" & abs(se) > 500)) |>
  # which distribution and reliability we want?
  filter(dist == "Normal", rel == "Rel = 0.8") |>
  mutate(
    param = factor(param, levels = c(rev(twofacpars), rev(growthpars))),
    n = factor(n, labels = paste0("n = ", c(15, 20, 50, 100, 1000))),
  )

plot_df50 <-
  res |>
  filter(method %in% c("ML", "eRBM", "iRBM")) |>
  filter( converged, !is.na(se)) |>
  # for each kind of model, filter bad standard errors
  filter(!(model == "twofac" & abs(se) > 5)) |>
  filter(!(model == "growth" & abs(se) > 500)) |>
  # which distribution and reliability we want?
  filter(dist == "Normal", rel == "Rel = 0.5") |>
  mutate(
    param = factor(param, levels = c(rev(twofacpars), rev(growthpars))),
    n = factor(n, labels = paste0("n = ", c(15, 20, 50, 100, 1000))),
  )

# For plots showing performance of D&R methods
plot_drcomp <-
  res |>
  filter(param %in% c(twofacpars, growthpars), dist != "Kurtosis") |>
  # filter(converged, !is.na(se)) |>
  # for each kind of model, filter bad standard errors
  filter(!(model == "twofac" & abs(se) > 5)) |>
  filter(!(model == "growth" & abs(se) > 500)) |>
  summarise(
    relbias = mean(relbias, na.rm = TRUE, trim = 0.05),
    rmse = sqrt(mean(bias ^ 2, na.rm = TRUE, trim = 0.05)),
    .by = c(dist:param)
  ) |>
  mutate(
    n = as.numeric(factor(n)),
    param = factor(param, levels = c(twofacpars, growthpars), labels = c(
      expression(Theta["1,1"]),
      expression(Psi["1,1"]),
      expression(Psi["2,2"]),
      expression(beta),
      expression(Lambda["2,1"]),
      expression(Theta["1,1"]),
      expression(Psi["1,1"]),
      expression(Psi["2,2"]),
      expression(Psi["1,2"])
    )),
    rel = factor(rel, labels = c(
      expression(Reliability~"="~0.8),
      expression(Reliability~"="~0.5)
    )),
    dist = factor(dist, labels = c(
      expression(Normal),
      # expression(Kurtosis),
      expression(plain("Non-normal"))
    )),
    method = factor(method, levels = rev(names(mycols)))
  )

## ----- Tables ----------------------------------------------------------------

create_summdf <- function(model, rel) {
  i <- res$rel == glue("Rel = {rel}") & res$model == model

  res |>
    filter(dist != "Kurtosis", param %in% c(twofacpars, growthpars),
           method != "lav", i) |>
    filter(!(model == "twofac" & abs(se) > 5)) |>
    filter(!(model == "growth" & abs(se) > 500)) |>
    summarise(
      mean_bias = mean(relbias, na.rm = TRUE, trim = 0.05),
      med_bias = median(relbias, na.rm = TRUE),
      rmse = sqrt(mean(relbias ^ 2, na.rm = TRUE, trim = 0.05)),
      .by = c(dist:method, param)
    ) |>
    arrange(dist, rel, n, param, desc(method)) |>
    select(-model, -rel) |>
    pivot_wider(
      names_from = c(dist, n),
      values_from = c(mean_bias, med_bias, rmse)
    ) |>
    group_by(param) |>
    arrange(param, method)
}

create_covrdf <- function(model) {
  i <- res$model == model
  if (model == "growth") i <- i & res$method != "REML"

  res |>
    filter(dist != "Kurtosis", param %in% c(twofacpars, growthpars), i) |>
    filter(!(model == "twofac" & abs(se) > 5)) |>
    filter(!(model == "growth" & abs(se) > 500)) |>
    summarise(
      covr = mean(covered, na.rm = TRUE),
      .by = c(dist:method, param)
    ) |>
    select(-model) |>
    pivot_wider(
      names_from = c(dist, rel, n),
      values_from = covr
    ) |>
    group_by(param)
}

tab_bias <- function(tab, model) {

  out <-
    tab |>
    gt(rowname_col = "method") |>
    fmt_number(decimals = 2) |>
    fmt_markdown(param) |>
    tab_spanner(
      columns = matches("mean_bias_Normal_"),
      label = "Normal data"
    ) |>
    tab_spanner(
      columns = matches("mean_bias_Non-Normal_"),
      label = "Non-normal data"
    ) |>
    tab_spanner(
      columns = starts_with("mean_bias"),
      label = "Relative mean bias"
    ) |>
    tab_spanner(
      columns = matches("med_bias_Normal_"),
      label = "Normal data",
      id = "med_bias_normal"
    ) |>
    tab_spanner(
      columns = matches("med_bias_Non-Normal_"),
      label = "Non-normal data",
      id = "med_bias_nonnormal"
    ) |>
    tab_spanner(
      columns = starts_with("med_bias"),
      label = "Relative median bias"
    ) |>
    tab_spanner(
      columns = matches("rmse_Normal_"),
      label = "Normal data",
      id = "rmse_normal"
    ) |>
    tab_spanner(
      columns = matches("rmse_Non-Normal_"),
      label = "Non-normal data",
      id = "rmse_nonnormal"
    ) |>
    tab_spanner(
      columns = starts_with("rmse"),
      label = "Relative RMSE"
    ) |>
    cols_label(
      ends_with("_15") ~ "15",
      ends_with("_20") ~ "20",
      ends_with("_50") ~ "50",
      ends_with("_100") ~ "100",
      ends_with("_1000") ~ "1000"
    ) |>
    tab_options(table.font.size = "10px")

  if (model == "twofac") {
    out <-
      out |>
      tab_row_group(label = md("$\\lambda_{2,1}$"), rows = param == "fx=~x2") |>
      tab_row_group(label = md("$\\beta$"), rows = param == "fy~fx") |>
      tab_row_group(label = md("$\\mathbf\\Psi_{2,2}$"), rows = param == "fy~~fy") |>
      tab_row_group(label = md("$\\mathbf\\Psi_{1,1}$"), rows = param == "fx~~fx") |>
      tab_row_group(label = md("$\\theta_{1,1}$"), rows = param == "y1~~y1")
  } else if (model == "growth") {
    out <-
      out |>
      tab_row_group(label = md("$\\mathbf\\Psi_{1,2}$"), rows = param == "i~~s") |>
      tab_row_group(label = md("$\\mathbf\\Psi_{2,2}$"), rows = param == "s~~s") |>
      tab_row_group(label = md("$\\mathbf\\Psi_{1,1}$"), rows = param == "i~~i") |>
      tab_row_group(label = md("$\\theta_{1}$"), rows = param == "v")
  }

  out
}

tab_covr <- function(tab, model) {

  out <-
    tab |>
    gt(rowname_col = "method") |>
    fmt_number() |>
    tab_spanner(
      columns = matches("Normal_Rel = 0.8", ignore.case = FALSE),
      label = "Reliability = 0.8"
    ) |>
    tab_spanner(
      columns = matches("Normal_Rel = 0.5", ignore.case = FALSE),
      label = "Reliability = 0.5"
    ) |>
    tab_spanner(
      columns = matches("Normal_", ignore.case = FALSE),
      label = "Normal data"
    ) |>
    tab_spanner(
      columns = matches("Non-normal_Rel = 0.8", ignore.case = FALSE),
      label = "Reliability = 0.8",
      id = "nonnormal80"
    ) |>
    tab_spanner(
      columns = matches("Non-normal_Rel = 0.5", ignore.case = FALSE),
      label = "Reliability = 0.5",
      id = "nonnormal50"
    ) |>
    tab_spanner(
      columns = matches("Non-normal_", ignore.case = FALSE),
      label = "Non-Normal data"
    ) |>
    cols_label(
      ends_with("_15") ~ "15",
      ends_with("_20") ~ "20",
      ends_with("_50") ~ "50",
      ends_with("_100") ~ "100",
      ends_with("_1000") ~ "1000"
    )

  if (model == "twofac") {
    out <-
      out |>
      tab_row_group(label = md("$\\lambda_{2,1}$"), rows = param == "fx=~x2") |>
      tab_row_group(label = md("$\\beta$"), rows = param == "fy~fx") |>
      tab_row_group(label = md("$\\mathbf\\Psi_{2,2}$"), rows = param == "fy~~fy") |>
      tab_row_group(label = md("$\\mathbf\\Psi_{1,1}$"), rows = param == "fx~~fx") |>
      tab_row_group(label = md("$\\theta_{1,1}$"), rows = param == "y1~~y1")
  } else if (model == "growth") {
    out <-
      out |>
      tab_row_group(label = md("$\\mathbf\\Psi_{1,2}$"), rows = param == "i~~s") |>
      tab_row_group(label = md("$\\mathbf\\Psi_{2,2}$"), rows = param == "s~~s") |>
      tab_row_group(label = md("$\\mathbf\\Psi_{1,1}$"), rows = param == "i~~i") |>
      tab_row_group(label = md("$\\theta_{1}$"), rows = param == "v")
  }

  out
}

bias_twofac_80_df <- create_summdf("twofac", 0.8)
bias_twofac_50_df <- create_summdf("twofac", 0.5)
covr_twofac_df    <- create_covrdf("twofac")
bias_growth_80_df <- create_summdf("growth", 0.8)
bias_growth_50_df <- create_summdf("growth", 0.5)
covr_growth_df    <- create_covrdf("growth")

## ----- SAVE RESULTS ----------------------------------------------------------
save(twofacpars, growthpars, mycols, simu_id,
     res_conv, res_timing_n1000, plot_df, plot_df50, plot_drcomp,
     bias_twofac_80_df, bias_twofac_50_df, covr_twofac_df,
     bias_growth_80_df, bias_growth_50_df, covr_growth_df,
     tab_bias, tab_covr,
     file = here::here("experiments/brsem/results.RData"))

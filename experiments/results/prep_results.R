library(tidyverse)
library(brlavaan)
library(glue)
library(gt)

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

## ----- Part 1: Our methods ---------------------------------------------------

# To create the following data frames:
#
# 1. simu_id: A data frame with the simulation ids. (in Part 2 add seeds)
#
# 2. res_ours_nested: A nested data frame with the results of each simulation, used
# for counting convergences etc.
#
# 3. res_ours: The full (raw) results data frame.
#
# 4. res_filtered: The filtered results data frame, with bad estimates removed.

load(here::here("experiments/simu_id.RData"))
simu_id <- bind_rows(
  mutate(simu_id, model = "twofac", .after = simid),
  mutate(simu_id, model = "growth", .after = simid)
)

load(here::here("experiments/simu_res_twofac_new.RData"))
load(here::here("experiments/simu_res_growth_new.RData"))
simu_res_ours <- c(simu_res_twofac, simu_res_growth)

res_ours_nested <-
  simu_res_ours |>
  imap(\(x, idx) bind_cols(simid = idx, x$simu_res)) |>
  bind_rows() |>
  mutate(
    super_converged = converged & #max_loglik < 0 &
      unlist(lapply(se, \(x) !any(is.na(x))))
  )

res_ours <-
  res_ours_nested |>
  mutate(param = lapply(truth, names), .before = est) |>
  unnest(param:truth) |>
  select(!starts_with("info")) |>
  distinct(sim, dist, model, rel, n, method, param, .keep_all = TRUE) |>
  mutate(
    type = case_when(
      grepl("f[xy]=~[xy][0-9]", param) ~ "Lambda",
      grepl("fy~fx", param) ~ "beta",
      grepl("[xy][0-9]~~[xy][0-9]|v", param) ~ "Theta",
      grepl("f[xy]~~f[xy]|[is]~~[is]", param) ~ "Psi",
      grepl("[is]~1", param) ~ "alpha",
      TRUE ~ NA
    ),
    across(c(est, truth, se), ~if_else(model == "growth" & type != "alpha", .x * 100, .x)),  # rescale back
    dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
    rel = factor(rel, levels = c("0.8", "0.5"), labels = c("Rel = 0.8", "Rel = 0.5")),
    method = factor(method, levels = rev(names(mycols))),
    bias = est - truth,
    relbias = bias / truth,
    covered = truth <= est + qnorm(0.975) * se & truth >= est - qnorm(0.975) * se
  )

## ----- Part 2: Comparison to D&R methods -------------------------------------

# To create the following data frames:
#
# 1. res_ours: The full (raw) results data frame.
#
# 2. truth: The truth data frame.
#
# 3. res_dr: The full (raw) results data frame from D&R.

truth <-
  expand_grid(
    model = c("twofac", "growth"),
    rel = factor(c("Rel = 0.8", "Rel = 0.5"), levels = c("Rel = 0.8", "Rel = 0.5"))
  ) |>
  mutate(out = map2(model, rel, \(x, y) {
    rel <- as.numeric(gsub("Rel = ", "", y))
    out <- if (x == "twofac") truth_twofac(rel)[twofacpars] else truth_growth(rel)[growthpars]
    list(param = names(out), truth = as.numeric(out))
  })) |>
  unnest_wider(out) |>
  unnest(c(param, truth))

# Others
dr_file1 <- here::here("experiments/GCM_est_combined_final.RData")
dr_file2 <- here::here("experiments/2FSEM_est_combined_final.RData")
if (!file.exists(dr_file1))
  download.file("https://osf.io/vjq5m/download", destfile = dr_file1)
if (!file.exists(dr_file2))
  download.file("https://osf.io/cw5b7/download", destfile = dr_file2)

load(here::here("experiments/2FSEM_est_combined_final.RData"))
res_dr_twofac <-
  as_tibble(Results) |>
  unite("simu", jobid, iteration, sep = ".") |>
  mutate(
    simu = as.integer(factor(simu)),
    convergence = convergence == 1
  ) |>
  select(
    seed,
    sim = simu,
    method,
    dist,
    n = nobs,
    converged = convergence,
    rel,
    starts_with("est"),
    starts_with("se")
  ) |>
  select(-ends_with("~1")) |>
  pivot_longer(
    cols = c(starts_with("est_"), starts_with("se_")),
    names_to = c(".value", "param"),
    names_sep = "_"
  ) |>
  filter(method %in% c("MLB", "JB", "BB", "Ozenne", "REML")) |>
  distinct(sim, seed, method, dist, n, rel, param, est, se, converged) |>
  mutate(model = "twofac")

load(here::here("experiments/GCM_est_combined_final.RData"))
res_dr_growth <-
  as_tibble(Results) |>
  unite("simu", jobid, iteration, sep = ".") |>
  mutate(
    simu = as.integer(factor(simu)),
    convergence = convergence == 1
  ) |>
  select(
    seed,
    sim = simu,
    method,
    dist,
    n = nobs,
    converged = convergence,
    rel,
    `est_v` = `est_Day0~~Day0`,
    `est_i~~i`,
    `est_i~1`,
    `est_s~~s`,
    `est_s~1`,
    `est_i~~s`,
    `se_v` = `se_Day0~~Day0`,
    `se_i~~i`,
    `se_i~1`,
    `se_s~~s`,
    `se_s~1`,
    `se_i~~s`
  ) |>
  pivot_longer(
    cols = c(starts_with("est_"), starts_with("se_")),
    names_to = c(".value", "param"),
    names_sep = "_"
  ) |>
  filter(method %in% c("MLB", "JB", "BB", "Ozenne", "REML")) |>
  distinct(sim, seed, method, dist, n, rel, param, est, se, converged) |>
  mutate(model = "growth")

res_dr <- bind_rows(res_dr_twofac, res_dr_growth)

# Clean up
res_dr$dist <- factor(
  res_dr$dist,
  levels = c("Normal", "Kurtosis", "NonNormal"),
  labels = c("Normal", "Kurtosis", "Non-normal")
)
res_dr$dist <- as.character(res_dr$dist)

res_dr$rel <- factor(
  res_dr$rel,
  levels = c("REL80", "REL50"),
  labels = c(0.8, 0.5)
)
res_dr$rel <- as.numeric(as.character(res_dr$rel))

res_dr$n <- as.numeric(as.character(res_dr$n))

tmp_lev <- levels(res_dr$method)
tmp_lev <- gsub("JB", "Jackknife", tmp_lev)
tmp_lev <- gsub("BB", "Bootstrap", tmp_lev)
tmp_lev <- gsub("Ozenne", "Ozenne et al.", tmp_lev)
tmp_lev <- gsub("REML", "REML", tmp_lev)
levels(res_dr$method) <- tmp_lev

# Add truth and simid
levels(truth$rel) <- c(0.8, 0.5)
truth$rel <- as.numeric(as.character(truth$rel))
res_dr <- left_join(res_dr, truth)
res_dr <- left_join(res_dr, select(simu_id, -seed))

# Rearrange columns
res_dr <-
  res_dr |>
  mutate(timing = NA, scaled_grad = NA, max_loglik = NA, optim_message = NA, Sigma_OK = NA) |>
  select(simid, seed, sim, dist, model, rel, n, method, param, est, se, truth,
         timing, converged, scaled_grad, max_loglik, Sigma_OK, optim_message) |>
  arrange(simid, sim, dist, model, rel, n, method, param) |>
  mutate(
    type = case_when(
      grepl("f[xy]=~[xy][0-9]", param) ~ "Lambda",
      grepl("fy~fx", param) ~ "beta",
      grepl("[xy][0-9]~~[xy][0-9]|v", param) ~ "Theta",
      grepl("f[xy]~~f[xy]|[is]~~[is]", param) ~ "Psi",
      grepl("[is]~1", param) ~ "alpha",
      TRUE ~ NA
    ),
    simid = as.integer(simid),
    dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
    rel = factor(rel, levels = c("0.8", "0.5"), labels = c("Rel = 0.8", "Rel = 0.5")),
    method = factor(method, levels = rev(names(mycols))),
    bias = est - truth,
    relbias = bias / truth,
    covered = truth <= est + qnorm(0.975) * se & truth >= est - qnorm(0.975) * se
  )

tinytest::expect_equal(colnames(res_ours), colnames(res_dr))

## ----- Part 3: Summarised data -----------------------------------------------

# For convergence statistics
res_conv <-
  res_ours_nested |>
  summarise(count = sum(super_converged, na.rm = TRUE), .by = dist:method) |>
  mutate(count = count / max(count), .by = dist:n) |>
  pivot_wider(names_from = c(model, method, rel), values_from = count) |>
  group_by(dist)

# For the plots showing performance of our methods
plot_df <-
  res_ours |>
  filter(converged, !is.na(se)) |>
  filter(dist == "Normal", rel == "Rel = 0.8",
         param %in% c(twofacpars, growthpars)) |>
  mutate(
    param = factor(param, levels = c(rev(twofacpars), rev(growthpars))),
    n = factor(n, labels = paste0("n = ", c(15, 20, 50, 100, 1000))),
  )

# For plots showing performance of D&R methods
summ_ours <-
  res_ours |>
  filter(param %in% c(twofacpars, growthpars), dist != "Kurtosis",
         converged, !is.na(se)) |>
  summarise(
    relbias = mean(relbias, na.rm = TRUE, trim = 0.05),
    rmse = sqrt(mean(bias ^ 2, na.rm = TRUE, trim = 0.05)),
    .by = c(dist:param)
  )
summ_dr <-
  res_dr |>
  filter(param %in% c(twofacpars, growthpars), dist != "Kurtosis") |>
  summarise(
    relbias = mean(relbias, na.rm = TRUE, trim = 0.05),
    rmse = sqrt(mean(bias ^ 2, na.rm = TRUE, trim = 0.05)),
    .by = c(dist:param)
  )

plot_drcomp <-
  bind_rows(summ_ours, summ_dr) |>
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

# bind_rows(
#   res_dr, res_ours
# ) |>
#   filter(model == "growth", dist == "Non-normal", n == 15, rel == "Rel = 0.5", param == "v") |>
#   ggplot(aes(est, fill = method)) +
#   geom_histogram() + facet_grid(method ~ .)

## ----- Part 4: Tables --------------------------------------------------------

create_summdf <- function(model, rel) {
  ress <- bind_rows(res_ours, res_dr)
  i <- ress$rel == glue("Rel = {rel}") & ress$model == model

  ress |>
    filter(dist != "Kurtosis", param %in% c(twofacpars, growthpars), i) |>
    filter(!method %in% c("eRBM") | abs(relbias) < 1) |>
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
    group_by(param)
}

create_covrdf <- function(model) {
  ress <- bind_rows(res_ours, res_dr)
  i <- ress$model == model
  if (model == "growth") i <- i & ress$method != "REML"

  ress |>
    filter(dist != "Kurtosis", param %in% c(twofacpars, growthpars), i) |>
    filter(!method %in% c("eRBM") | abs(relbias) < 1) |>
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
     res_conv, plot_df, plot_drcomp,
     bias_twofac_80_df, bias_twofac_50_df, covr_twofac_df,
     bias_growth_80_df, bias_growth_50_df, covr_growth_df,
     tab_bias, tab_covr,
     file = here::here("experiments", "results.RData"))

# res_dr |>
#   filter(param == "v", n == 15) |>
#   summarise(est = modeest::mfv1(est, na_rm = TRUE), .by = dist:method) |>
#   pivot_wider(names_from = rel, values_from = est) |> print(n = Inf)

# bind_rows(res_dr, res_ours) |>
#   filter(param == "fy~fx", n == 15) |>
#   summarise(est = mean(est, na.rm = TRUE), .by = dist:method) |>
#   pivot_wider(names_from = rel, values_from = est) |>
#   arrange(dist:method) |>
#   print(n = Inf)
#
#
# ## TRANSFORM???
#
# x <-
#   res_ours |>
#   filter(model == "growth", dist == "Normal", rel == "Rel = 0.8", n == 15,
#          param == "v")
# xx <- x$est[x$method == "ML"]
#
# y <-
#   res_dr |>
#   filter(model == "growth", dist == "Normal", rel == "Rel = 0.8", n == 15,
#          param == "v")
# yy <- y$est[y$method == "MLB"]
#
# plot(sort(xx), sort(yy))
# mod <- lm(sort(yy) ~ I(50 + sort(xx)))
# summary(mod)
# plot(50 + sort(xx), sort(yy))

library(brlavaan)
library(tidyverse)
library(furrr)
theme_set(theme_bw())
here::i_am("experiments/obs_vs_exp/10-run_sims.R")
conflicted::conflict_prefer("filter", "dplyr")
mycols <- c(
  ML = "#E31A1C",
  eRBM = "#A6CEE3",
  iRBM = "#1F78B4",
  Jackknife = "#B2DF8A",
  Bootstrap = "#33A02C",
  `Ozenne et al.` = "#FDBF6F",
  REML = "#FF7F00"
)

## ----- Run simulations -------------------------------------------------------
ncores <- future::availableCores() - 2
future::plan(multisession, workers = ncores)
B <- 1000  # Number of simulations

simu_id <-
  expand_grid(
    dist = "Normal",
    model = c("twofac", "growth"),
    rel = 0.8,
    n = c(15, 20, 50, 100, 1000),
    info_pen = c("observed", "expected"),
    info_bias = c("observed", "expected"),
    info_se = c("observed", "expected")
  ) |>
  rownames_to_column(var = "simid")

# simu_res <- vector("list", length = nrow(simu_id))
# for (i in seq_len(nrow(simu_id))) {
#   dist  <- simu_id$dist[i]
#   model <- simu_id$model[i]
#   rel   <- simu_id$rel[i]
#   n     <- simu_id$n[i]
#
#   cli::cli_inform(">>> {Sys.time()} <<<\n\n[{i} / {nrow(simu_id)}] Now running {model} models ({dist}) rel = {rel}, n = {n}\n")
#   simu_res[[i]] <- sim_fun(
#     dist = dist,
#     model = model,
#     rel = rel,
#     n = n,
#     lavsim = FALSE,
#     lavfun = ifelse(model == "twofac", "sem", "growth"),
#     nsimu = B,
#     info_pen = simu_id$info_pen[i],
#     info_bias = simu_id$info_bias[i],
#     info_se = simu_id$info_se[i]
#   )
#   cat("\n")
#   save(simu_res, file = "experiments/obs_vs_exp/ove.RData")
# }

## ----- Analyse ---------------------------------------------------------------
load("experiments/obs_vs_exp/ove_mp1.RData")
simu_res1 <- simu_res
load("experiments/obs_vs_exp/ove_mp2.RData")
simu_res2 <- simu_res

# Weave results
simu_res <-
  map2(simu_res1, simu_res2, \(x, y) {
    simu_res <- bind_rows(x$simu_res, y$simu_res)
    errors <- c(x$error, y$error)
    list(simu_res = simu_res, error = errors)
  })

# Create results data frame
res <-
  simu_res |>
  imap(\(x, idx) bind_cols(simid = idx, x$simu_res)) |>
  bind_rows() |>
  drop_na()

# Where problems?
simu_id |>
  mutate(errors = map(simu_res, \(x) x$error)) |>
  unnest(errors)

# How many converged?
res |>
  summarise(
    count = sum(converged),
    .by = dist:method
  ) |>
  pivot_wider(names_from = method, values_from = count) |>
  print(n = Inf)

# False convergences
res |>
  filter(method == "iRBM") |>
  count(model, converged, optim_message)

# This one for bias plots etc
fullres <-
  res |>
  mutate(converged = case_when(
    optim_message == "false convergence (8)" ~ TRUE,
    TRUE ~ converged
  )) |>
  mutate(converged = all(converged), .by = c(simid, sim)) |>
  filter(converged) |>
  mutate(param = lapply(truth, names), .before = est) |>
  unnest(param:truth) |>
  filter(truth != 0) |>  # otherwise, cannot get relative bias
  mutate(
    bias = est - truth,
    relbias = bias / truth,
    covered = truth <= est + qnorm(0.975) * se & truth >= est - qnorm(0.975) * se
  )

fullres |>
  mutate(
    info_pen = factor(info_pen, labels = c("Pen = exp", "Pen = obs")),
    info_bias = factor(info_bias, labels = c("Bias = exp", "Bias = obs")),
    info_se = factor(info_se, labels = c("SE = exp", "SE = obs")),
  ) |>
  summarise(
    mean_bias = mean(relbias, na.rm = TRUE, trim = 0.05),
    .by = dist:method
  ) |>
  ggplot(aes(as.numeric(as.factor(n)), mean_bias, col = method)) +
  geom_line() +
  geom_hline(yintercept = 0, linetype = 2) +
  facet_grid(info_bias * info_se ~ model * info_pen) +
  scale_colour_manual(values = mycols) +
  scale_y_continuous(labels = scales::percent, name = "Relative mean bias") +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000), name = "Sample size") +
  theme_bw()

# Look at distribution
fullres |>
  filter(info_pen == info_bias) |>
  ggplot(aes(relbias, param, fill = method)) +
  geom_boxplot(outliers = FALSE) +
  facet_grid(model ~ info_pen, scales = "free_y", space = "free") +
  scale_fill_manual(values = mycols) +
  theme_bw() +
  theme(legend.position = "top")

# Coverage
fullres |>
  mutate(
    info_pen = factor(info_pen, labels = c("Pen = exp", "Pen = obs")),
    info_bias = factor(info_bias, labels = c("Bias = exp", "Bias = obs")),
    info_se = factor(info_se, labels = c("SE = exp", "SE = obs")),
  ) |>
  filter (info_pen == "Pen = obs", as.numeric(info_pen) == as.numeric(info_bias)) |>
  summarise(
    coverage = mean(covered, na.rm = TRUE),
    .by = dist:method
  ) |>
  ggplot(aes(as.numeric(as.factor(n)), coverage, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 0.95, linetype = 2) +
  facet_grid(info_bias * info_se ~ model * info_pen) +
  scale_fill_manual(values = mycols) +
  scale_y_continuous(labels = scales::percent, name = "Coverage rate") +
  scale_x_continuous(breaks = 1:5, labels = c(15, 20, 50, 100, 1000), name = "Sample size") +
  theme_bw() +
  coord_cartesian(ylim = c(0.7, 1))

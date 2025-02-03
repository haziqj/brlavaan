## ----- Organise simulation results -------------------------------------------
load(here::here("experiments/simu_res_growth.RData"))
load(here::here("experiments/simu_res_twofac.RData"))
growthpars <- c("v", "i~~i", "s~~s", "i~~s")
twofacpars <- c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")

# Create a simid table
simu_res <- c(simu_res_growth, simu_res_twofac)
simu_id <-
  expand_grid(
    model = c("growth", "twofac"),
    dist = c("Normal", "Kurtosis", "Non-normal"),
    rel = as.character(c(0.8, 0.5)),
    n = c(15, 20, 50, 100, 1000)
  ) |>
  rownames_to_column(var = "simid") |>
  mutate(simid = as.integer(simid))

# Create (nested) results data frame, used for counting convergences etc.
res_nested <-
  simu_res |>
  imap(\(x, idx) bind_cols(simid = idx, x$simu_res)) |>
  bind_rows() |>
  drop_na() |>
  select(-starts_with("info"))

# This is the full (raw) results data frame
res <-
  res_nested |>
  mutate(param = lapply(truth, names), .before = est) |>
  unnest(param:truth) |>
  select(!starts_with("info")) |>
  # For the growth model, keep the first instance of param == "v"
  filter(
    row_number() == 1,
    .by = c(simid, sim, dist, model, rel, n, method, param)
  )

# Create truth data frame
truth <-
  res |>
  summarise(
    across(truth, first),
    .by = c(model, param, rel)
  )

res |>
  filter(truth != 0, param %in% growthpars, dist != "Kurtosis") |>
  filter(converged, Sigma_OK) |>
  summarise(
    bias = mean((est - truth) / truth, na.rm = TRUE, trim = 0.05),
    .by = dist:param
  ) |>
  filter(model == "growth") |>
  mutate(n = factor(n)) |>
  ggplot(aes(as.numeric(n), bias, col = method)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_colour_manual(values = mycols) +
  facet_grid(param ~ dist * rel) +
  theme_bw()











## ----- Organise D&R sims -----------------------------------------------------
dr_file1 <- here::here("experiments/GCM_est_combined_final.RData")
dr_file2 <- here::here("experiments/2FSEM_est_combined_final.RData")
if (!file.exists(dr_file1))
  download.file("https://osf.io/vjq5m/download", destfile = dr_file1)
if (!file.exists(dr_file2))
  download.file("https://osf.io/cw5b7/download", destfile = dr_file2)

# Growth model
load("experiments/GCM_est_combined_final.RData")
res_dr_growth <-
  as_tibble(Results) |>
  unite("simu", jobid, iteration, sep = ".") |>
  mutate(
    simu = as.integer(factor(simu)),
    convergence = convergence == 1
  ) |>
  select(
    sim = simu,
    method,
    dist,
    n = nobs,
    converged = convergence,
    rel,
    `est_Day0~~Day0`,
    `est_i~~i`,
    `est_i~1`,
    `est_s~~s`,
    `est_s~1`,
    `est_i~~s`,
    `se_Day0~~Day0`,
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
  filter(method %in% c("JB", "BB", "Ozenne", "REML")) |>
  distinct(sim, method, dist, n, rel, param, est, se, converged) |>
  mutate(model = "growth")

# Two-factor SEM
load("experiments/2FSEM_est_combined_final.RData")
res_dr_twofac <-
  as_tibble(Results) |>
  unite("simu", jobid, iteration, sep = ".") |>
  mutate(
    simu = as.integer(factor(simu)),
    convergence = convergence == 1
  ) |>
  select(
    sim = simu,
    method,
    dist,
    n = nobs,
    converged = convergence,
    rel,
    `est_y1~~y1`,
    `est_fx~~fx`,
    `est_fy~~fy`,
    `est_fy~fx`,
    `est_fx=~x2`,
    `se_y1~~y1`,
    `se_fx~~fx`,
    `se_fy~~fy`,
    `se_fy~fx`,
    `se_fx=~x2`
  ) |>
  pivot_longer(
    cols = c(starts_with("est_"), starts_with("se_")),
    names_to = c(".value", "param"),
    names_sep = "_"
  ) |>
  filter(method %in% c("JB", "BB", "Ozenne", "REML")) |>
  distinct(sim, method, dist, n, rel, param, est, se, converged) |>
  mutate(model = "twofac")

# Combine results
res_dr <- bind_rows(res_dr_growth, res_dr_twofac)

# Clean up
levels(res_dr$dist) <- c("Kurtosis", "Non-normal", "Normal")
res_dr$dist <- factor(res_dr$dist, levels = c("Normal", "Kurtosis", "Non-normal"))

res_dr$param <- gsub("Day0~~Day0", "v", res_dr$param)

levels(res_dr$rel) <- c("0.5", "0.8")
res_dr$rel <- as.character(res_dr$rel)

res_dr$n <- as.numeric(as.character(res_dr$n))

tmp_lev <- levels(res_dr$method)
tmp_lev <- gsub("JB", "Jackknife", tmp_lev)
tmp_lev <- gsub("BB", "Bootstrap", tmp_lev)
tmp_lev <- gsub("Ozenne", "Ozenne et al.", tmp_lev)
tmp_lev <- gsub("REML", "REML", tmp_lev)
levels(res_dr$method) <- tmp_lev

# Add truth and simid
res_dr <- left_join(res_dr, truth)
res_dr <- left_join(res_dr, simu_id)

# Rearrange columns
res_dr <-
  res_dr |>
  mutate(timing = NA, optim_message = NA) |>
  select(simid, sim, dist, model, rel, n, method, param, est, se, truth, timing,
         converged, optim_message) |>
  arrange(simid, sim, dist, model, rel, n, method, param)

## ----- Combine our simulation results with D&R and simulation results --------
res_full <-
  bind_rows(
    bind_cols(res, tibble(resfrom = "ours")),
    bind_cols(res_dr, tibble(resfrom = "dr"))
  ) |>
  filter(truth != 0) |>  # otherwise, cannot get relative bias
  mutate(
    bias = est - truth,
    relbias = bias / truth,
    covered = truth <= est + qnorm(0.975) * se & truth >= est - qnorm(0.975) * se
  )

growthpars <- c("v", "i~~i", "s~~s", "i~~s")
twofacpars <- c("y1~~y1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2")

res_summ <-
  res_full |>
  # Here's where the non-convergences are filtered out
  mutate(converged = case_when(
    optim_message == "false convergence (8)" ~ TRUE,
    TRUE ~ converged
  )) |>
  mutate(converged = all(converged), .by = c(simid, sim, resfrom)) |>
  filter(converged) |>
  # Also only keep parameters we are interested in
  filter(param %in% c(growthpars, twofacpars)) |>
  summarise(
    count = n(),
    mean_bias = mean(relbias, na.rm = TRUE, trim = 0.05),
    median_bias = median(relbias, na.rm = TRUE),
    rmse = sqrt(mean(bias ^ 2, na.rm = TRUE, trim = 0.05)),
    coverage = mean(covered, na.rm = TRUE),
    .by = dist:param
  )

# Cleanup
rm(Results)
rm(tmp_lev)

list(
  res = res,
  res_nested = res_nested,
  res_dr = res_dr,
  res_full = res_full,
  res_summ = res_summ
) |>
  map(\(df) {
    df |>
      mutate(
        dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
        rel = factor(rel, levels = c("0.8", "0.5"), labels = c("Rel = 0.8", "Rel = 0.5")),
        method = factor(method, levels = rev(names(mycols))),
        n = as.factor(n)
      )
  }) |>
  list2env(envir = .GlobalEnv)


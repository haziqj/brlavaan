# Purpose of this script is to prepare the raw simulation results from Dhaene &
# Rosseel (2022) for the two factor model only (the growth model was re-run by
# us). Need only run this script once.

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

# dr_file1 <- here::here("experiments/GCM_est_combined_final.RData")
dr_file2 <- here::here("experiments/brsem/2FSEM_est_combined_final.RData")
# if (!file.exists(dr_file1))
#   download.file("https://osf.io/vjq5m/download", destfile = dr_file1)
if (!file.exists(dr_file2))
  download.file("https://osf.io/cw5b7/download", destfile = dr_file2)

load(here::here("experiments/brsem/2FSEM_est_combined_final.RData"))
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
  filter(method %in% c("JB", "BB", "Ozenne", "REML")) |>
  distinct(sim, seed, method, dist, n, rel, param, est, se, converged) |>
  mutate(model = "twofac")

# load(here::here("experiments/GCM_est_combined_final.RData"))
# res_dr_growth <-
#   as_tibble(Results) |>
#   unite("simu", jobid, iteration, sep = ".") |>
#   mutate(
#     simu = as.integer(factor(simu)),
#     convergence = convergence == 1
#   ) |>
#   select(
#     seed,
#     sim = simu,
#     method,
#     dist,
#     n = nobs,
#     converged = convergence,
#     rel,
#     `est_v` = `est_Day0~~Day0`,
#     `est_i~~i`,
#     `est_i~1`,
#     `est_s~~s`,
#     `est_s~1`,
#     `est_i~~s`,
#     `se_v` = `se_Day0~~Day0`,
#     `se_i~~i`,
#     `se_i~1`,
#     `se_s~~s`,
#     `se_s~1`,
#     `se_i~~s`
#   ) |>
#   pivot_longer(
#     cols = c(starts_with("est_"), starts_with("se_")),
#     names_to = c(".value", "param"),
#     names_sep = "_"
#   ) |>
#   filter(method %in% c("MLB", "JB", "BB", "Ozenne", "REML")) |>
#   distinct(sim, seed, method, dist, n, rel, param, est, se, converged) |>
#   mutate(model = "growth")
res_dr_growth <- NULL

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
         timing, converged, Sigma_OK) |>
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
  ) |>
  select(-simid)

tinytest::expect_equal(colnames(res), colnames(res_dr))
save(res_dr, file = here::here("experiments/brsem/simu_res_twofacDR.RData"))

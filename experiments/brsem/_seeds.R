library(tidyverse)
here::i_am("experiments/obtain_seeds.R")

simu_id <-
  expand_grid(
    dist = c("Normal", "Kurtosis", "Non-normal"),
    rel = c(0.8, 0.5),
    n = c(15, 20, 50, 100, 1000)
  ) |>
  rownames_to_column(var = "simid")

# dr_file1 <- here::here("experiments/GCM_est_combined_final.RData")
dr_file2 <- here::here("experiments/2FSEM_est_combined_final.RData")
# if (!file.exists(dr_file1))
#   download.file("https://osf.io/vjq5m/download", destfile = dr_file1)
if (!file.exists(dr_file2))
  download.file("https://osf.io/cw5b7/download", destfile = dr_file2)

# Two-factor SEM
load(here::here("experiments/2FSEM_est_combined_final.RData"))
simu_id <-
  Results |>
  as_tibble() |>
  select(dist, rel, n = nobs, seed) |>
  mutate(
    dist = case_when(
      dist == "NonNormal" ~ "Non-normal",
      TRUE ~ as.character(dist)
    ),
    rel = as.numeric(gsub("REL", "", as.character(rel))) / 100,
    n = as.numeric(as.character(n))
  ) |>
  distinct(dist, rel, n, seed) |>
  nest(seed = seed) |>
  mutate(seed = map(seed, \(x) unlist(x, use.names = FALSE))) |>
  left_join(x = simu_id)
save(simu_id, file = here::here("experiments/simu_id.RData"))

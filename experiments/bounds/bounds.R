library(brlavaan)
library(tidyverse)
library(furrr)
theme_set(theme_bw())

bounds_sim_fn <- function(
    dist = "Normal",
    model = "growth",
    rel = 0.8,
    n = 15,
    nsimu = 1,
    data_scale = 1 / 10,
    seeds = NULL
) {

  # Check seeds
  if (!is.null(seeds)) {
    if (length(seeds) != nsimu) {
      cli::cli_abort("Length of seeds must be equal to nsimu")
    }
  }

  # Generate all data ----------------------------------------------------------
  if (model == "growth") {
    gen_data <- gen_data_growth
    txt_mod <- txt_mod_growth
    lavfun <- "growth"
  } else if (model == "twofac") {
    gen_data <- gen_data_twofac
    txt_mod <- txt_mod_twofac  # this has arg meanstructure = FALSE by default
    lavfun <- "sem"
  } else {
    cli::cli_abort("Unknown model: {model}")
  }
  datasets <- NULL
  if (!is.null(seeds)) {
    datasets <- purrr::map(
      .x = seeds,
      .f = \(s) gen_data(n = n, rel = rel, dist = dist, lavsim = FALSE,
                         scale = data_scale, seed = s)
    )
  }
  mod <- txt_mod(rel)

  # Single run function --------------------------------------------------------
  single_sim <- function(j) {
    if (is.null(datasets)) {
      dat <- gen_data(n = n, rel = rel, dist = dist, lavsim = FALSE,
                      scale = data_scale)
    } else {
      dat <- datasets[[j]]
    }
    true_vals <- truth(dat)
    fit_list <- list()

    # lavaan fit
    fit_list$lav <-
      growth(mod, dat, start = true_vals, bounds = "none")
    fit_list$lavb <-
      growth(mod, dat, start = true_vals, bounds = "standard")

    # fit_sem
    fit_list$ML <-
      fit_sem(mod, dat, rbm = "none", plugin_pen = NULL, maxgrad = FALSE,
              lavfun = lavfun, bounds = "none", start = true_vals)
    fit_list$MLb <-
      fit_sem(mod, dat, rbm = "none", plugin_pen = NULL, maxgrad = FALSE,
              lavfun = lavfun, bounds = "standard", start = true_vals)

    nsimtypes <- length(fit_list)

    tibble::tibble(
      seed = seeds[j],
      sim = j,
      dist = dist,
      model = model,
      rel = rel,
      n = n,
      method = names(fit_list),
      est = lapply(fit_list, purrr::possibly(coef, NA)),
      truth = rep(list(true_vals), nsimtypes)
    )
  }

  # Run simulation -------------------------------------------------------------
  simu_res <- furrr::future_map(
    seq_len(nsimu),
    purrr::possibly(single_sim),
    .progress = TRUE,
    .options = furrr::furrr_options(seed = TRUE)
  )
  do.call("rbind", simu_res)
}

# Run simulation
ncores <- future::availableCores() - 1
future::plan(multisession, workers = ncores)
B <- 1000  # no. of sims

simu_id <-
  expand_grid(
    dist = c("Normal", "Kurtosis", "Non-normal"),
    rel = c(0.8, 0.5),
    n = c(15, 20, 50, 100, 1000)
  ) |>
  rownames_to_column(var = "simid") |>
  mutate(seeds = map(simid, \(x) 1234 * as.numeric(x) + seq_len(B)))

simu_res <- vector("list", length = nrow(simu_id))
for (i in seq_len(nrow(simu_id))) {
  dist  <- simu_id$dist[i]
  rel   <- simu_id$rel[i]
  n     <- simu_id$n[i]
  seeds <- simu_id$seeds[[i]]
  model <- "growth"

  cli::cli_inform(">>> {Sys.time()} <<<\n\n[{i} / {nrow(simu_id)}] Now running {model} models ({dist}) rel = {rel}, n = {n}\n")
  simu_res[[i]] <- bounds_sim_fn(
    dist = dist,
    model = model,
    rel = rel,
    n = n,
    nsimu = B,
    data_scale = 1 / 10,
    seeds = seeds
  )
  cat("\n")
  save(simu_res, file = "experiments/bounds/simu_res_bounds.RData")
}

# Analyse results
growthpars <- c("v", "i~~i", "s~~s", "i~~s")
mycols <- c(
  ML = "#E31A1C",
  MLb = "#FB9A99",
  lav = "#1F78B4",
  lavb = "#A6CEE3"
)

simu_res |>
  bind_rows() |>
  filter(dist != "Kurtosis") |>
  mutate(
    est = map(est, as.numeric),
    param = map(truth, names)
  ) |>
  unnest(c(est, truth, param)) |>
  filter(param %in% growthpars) |>
  distinct(across(seed:method), param, .keep_all = TRUE) |>
  mutate(bias = (est - truth) / truth) |>
  summarise(bias = mean(bias, trim = 0.05), .by = c(dist, model, rel, n, method, param)) |>
  mutate(
    n = as.numeric(factor(n)),
    param = factor(param, levels = growthpars),
    dist = factor(dist, levels = c("Normal", "Non-normal")),
    rel = factor(rel, levels = c(0.8, 0.5), labels = paste0("Rel =", c(0.8, 0.5)))
  ) |>
  ggplot(aes(n, bias, col = method)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggh4x::facet_nested(param ~ rel + dist, scales = "free_y") +
  scale_color_manual(values = mycols) +
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(-0.4, 0.2))

# CONCLUSION:
# - Using bounds = "none" gets us the expected mean bias results.

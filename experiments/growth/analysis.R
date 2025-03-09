library(tidyverse)
library(brlavaan)
here::i_am("experiments/growth/analysis.R")
load(here::here("experiments/simu_res_growthmp1.RData"))
simu_res_MP1 <- simu_res_growth
load(here::here("experiments/simu_res_growthmp2.RData"))
simu_res_MP2 <- simu_res_growth
load(here::here("experiments/simu_res_growthmp3.RData"))
simu_res_MP3 <- simu_res_growth

simu_id <-
  expand_grid(
    dist = c("Normal", "Kurtosis", "Non-normal"),
    model = "growth",
    rel = c(0.8, 0.5),
    n = c(15, 20, 50, 100, 1000)
  ) |>
  rownames_to_column(var = "simid")

simu_res <- vector("list", length = nrow(simu_id))
simu_res[which(sapply(simu_res_MP1, \(x) !is.null(x)))] <-
  simu_res_MP1[which(sapply(simu_res_MP1, \(x) !is.null(x)))]
simu_res[which(sapply(simu_res_MP2, \(x) !is.null(x)))] <-
  simu_res_MP2[which(sapply(simu_res_MP2, \(x) !is.null(x)))]
simu_res[which(sapply(simu_res_MP3, \(x) !is.null(x)))] <-
  simu_res_MP3[which(sapply(simu_res_MP3, \(x) !is.null(x)))]

growthpars <- c("v", "i~~i", "s~~s", "i~~s")
mycols <- c(
  ML = "#E31A1C",
  eRBM = "#A6CEE3",
  iRBM = "#1F78B4",
  Jackknife = "#B2DF8A",
  Bootstrap = "#33A02C",
  `Ozenne et al.` = "#FDBF6F",
  REML = "#FF7F00"
)

## ----- Prep results ----------------------------------------------------------

# Create (nested) results data frame, used for counting convergences etc.
res_nested <-
  simu_res |>
  imap(\(x, idx) bind_cols(simid = idx, x$simu_res)) |>
  bind_rows() |>
  drop_na() |>
  select(-starts_with("info")) |>
  mutate(
    super_converged = converged & max_loglik < 0 & Sigma_OK &
      unlist(lapply(scaled_grad, \(x) abs(max(x)) < 1e-2))
  )

# This is the full (raw) results data frame
res <-
  res_nested |>
  mutate(param = lapply(truth, names), .before = est) |>
  unnest(param:truth) |>
  select(!starts_with("info")) |>
  # For the growth model, keep the first instance of param == "v"
  distinct(simid, sim, dist, model, rel, n, method, param, .keep_all = TRUE) |>
  mutate(type = case_when(
    grepl("~1", param) ~ "alpha",
    grepl("~~", param) ~ "Psi",
    param == "v" ~ "Theta",
    TRUE ~ NA
  ))

## ----- Convergence -----------------------------------------------------------
tab_conv <-
  res_nested |>
  summarise(
    count = sum(converged),
    .by = dist:method
  ) |>
  pivot_wider(names_from = method, values_from = count) |>
  select(-starts_with("info")); print(tab_conv, n = Inf)

res_nested |>
  summarise(
    count = sum(converged),
    .by = dist:method
  ) |>
  mutate(
    dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
    rel = factor(rel, levels = c("0.8", "0.5"), labels = c("Rel = 0.8", "Rel = 0.5")),
    method = factor(method, levels = rev(names(mycols))),
    n = factor(n)
  ) |>
  ggplot(aes(n, count, fill = method)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_hline(yintercept = 1000, linetype = "dashed") +
  facet_grid(rel ~ dist) +
  scale_fill_manual(values = mycols) +
  scale_y_log10() +
  theme_bw() +
  labs(x = "Simulations", y = NULL, fill = NULL) +
  theme(legend.position = "none")

# res_nested |>
#   summarise(
#     count = sum(converged),
#     .by = dist:method
#   ) |>
#   mutate(
#     count = count / max(count),
#     .by = dist:n
#   ) |>
#   pivot_wider(names_from = method, values_from = count) |>
#   select(-starts_with("info")) |>
#   print(n = Inf)

p_conv_succ <-
  res_nested |>
  summarise(
    count = sum(converged),
    .by = dist:method
  ) |>
  mutate(
    count = count / max(count),
    .by = dist:n
  ) |>
  mutate(
    dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
    rel = factor(rel, levels = c("0.8", "0.5"), labels = c("Rel = 0.8", "Rel = 0.5")),
    method = factor(method, levels = rev(names(mycols))),
    n = factor(n)
  ) |>
  ggplot(aes(n, count, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  facet_grid(rel ~ dist) +
  scale_fill_manual(values = mycols) +
  scale_y_continuous(labels = scales::percent) +
  theme_bw() +
  labs(x = "Sample size (n)", y = "Success prop. (rel. to ML)", fill = NULL) +
  theme(legend.position = "none"); p_conv_succ

## ----- Plots -----------------------------------------------------------------


# Completion plot
p1 <-
  res_nested |>
  summarise(
    count = sum(converged),
    .by = dist:method
  ) |>
  mutate(
    dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
    rel = factor(rel, levels = c("0.8", "0.5"), labels = c("Rel = 0.8", "Rel = 0.5")),
    method = factor(method, levels = rev(names(mycols)))
  ) |>
  filter(n == 50, dist != "Kurtosis") |>
  ggplot(aes(count, method, fill = method)) +
  geom_bar(stat = "identity", width = 0.75) +
  geom_vline(xintercept = 1000, linetype = "dashed") +
  ggh4x::facet_nested(dist + rel ~ .) +
  scale_fill_manual(values = mycols) +
  scale_x_continuous(labels = scales::label_number(scale = 1/1000, suffix = "K")) +
  theme_bw() +
  labs(x = "Simulations", y = NULL, fill = NULL, title = " ") +
  theme(legend.position = "none"); p1

res_filtered <-
  res |>
  mutate(converged = all(converged), .by = c(simid, sim)) |>
  filter(converged) |>
  mutate(
    dist = factor(dist, levels = c("Normal", "Kurtosis", "Non-normal")),
    rel = factor(rel, levels = c("0.8", "0.5"), labels = c("Rel = 0.8", "Rel = 0.5")),
    method = factor(method, levels = rev(names(mycols))),
    bias = est - truth,
    relbias = bias / truth,
    covered = truth <= est + qnorm(0.975) * se & truth >= est - qnorm(0.975) * se
  )

truth <-
  res_filtered |>
  summarise(
    across(truth, first),
    .by = c(model, param, rel)
  )

# Centred distributions
p2 <-
  res_filtered |>
  filter(n == 50, dist != "Kurtosis", abs(bias) < 7) |>
  ggplot(aes(bias, type, fill = method)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = mycols) +
  # scale_y_discrete(labels = rev(c(
  #   expression(Theta),
  #   expression(Psi),
  #   expression(Lambda),
  #   expression(beta)
  # ))) +
  ggh4x::facet_nested(dist + rel ~ .) +
  theme_bw() +
  # coord_cartesian(xlim = c(-4, 4)) +
  labs(x = "Bias", y = NULL, fill = NULL, title = "Sample size n = 50") +
  theme(legend.position = "none"); p2

cowplot::plot_grid(
  p2, p1,
  rel_widths = c(4, 1)
)

p_dist_n50_all <-
  res_filtered |>
  filter(n == 50, dist == "Normal", abs(bias) < 7) |>
  ggplot(aes(bias, param, fill = method)) +
  geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = mycols) +
  facet_grid(rel ~ .) +
  theme_bw() +
  labs(x = "Bias", y = NULL, fill = NULL, title = "Sample size n = 50, normal dist.") +
  theme(legend.position = "bottom"); p_dist_n50_all

# res_filtered |>
#   filter(n == 50, dist != "Kurtosis", abs(relbias) < 2) |>
#   filter(param %in% twofacpars) |>
#   ggplot(aes(bias, param, fill = method)) +
#   geom_boxplot(alpha = 0.8, outlier.size = 0.5) +
#   geom_vline(xintercept = 0, linetype = "dashed") +
#   scale_fill_manual(values = mycols) +
#   ggh4x::facet_nested(dist + rel ~ .) +
#   theme_bw() +
#   labs(x = "Bias", y = NULL, fill = NULL, title = "Sample size n = 50") +
#   theme(legend.position = "bottom")

# Performance plot
the_trim <- 0.05
p_perf_n50_all <-
  res_filtered |>
  filter(n == 50, dist == "Normal", abs(bias) < 2) |>
  summarise(
    rmse = sqrt(mean(bias ^ 2, na.rm = TRUE, trim = the_trim)),
    meanbias = mean(bias, trim = the_trim),
    pu = mean(bias < 0),
    coverage = mean(covered, na.rm = TRUE),
    .by = c(dist:param, type)
  ) |>
  pivot_longer(rmse:coverage, names_to = "metric", values_to = "value") |>
  mutate(
    metric = factor(
      metric,
      levels = c("meanbias", "rmse", "pu", "coverage"),
      labels = c("Mean~bias", "RMSE", "Prob.~underest.",  "Coverage")
    ),
    rel = factor(rel, labels = c("Rel == 0.8", "Rel == 0.5")),
  ) |>
  ggplot(aes(value, param, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_vline(
    data = tibble(
      metric = factor(c("Mean~bias", "RMSE", "Prob.~underest.",  "Coverage")),
      value = c(0, 0, 0.5, 0.95),
    ),
    aes(xintercept = value),
    linetype = "dashed"
  ) +
  scale_fill_manual(values = mycols) +
  ggh4x::facet_nested(rel * type ~ metric, scales = "free", space = "free_y",
                      labeller = label_parsed) +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_continuous(limits = c(-0.3, 0.3)),
      scale_x_continuous(),
      scale_x_continuous(limits = c(0.35, 0.65), labels = scales::percent),
      scale_x_continuous(limits = c(0.7, 1), labels = scales::percent)
    )
  ) +
  theme_bw() +
  labs(x = NULL, y = NULL, fill = NULL, title = "Normal dist., n = 50") +
  theme(legend.position = "bottom"); p_perf_n50_all

the_trim <- 0.05
p_perf_n50_type <-
  res_filtered |>
  filter(n == 50, abs(bias) < 2) |>
  summarise(
    rmse = sqrt(mean(bias ^ 2, na.rm = TRUE, trim = the_trim)),
    meanbias = mean(relbias, trim = the_trim),
    pu = mean(relbias < 0),
    coverage = mean(covered, na.rm = TRUE),
    .by = c(dist:method, type)
  ) |>
  pivot_longer(rmse:coverage, names_to = "metric", values_to = "value") |>
  mutate(
    metric = factor(
      metric,
      levels = c("meanbias", "rmse", "pu", "coverage"),
      labels = c("Mean bias", "RMSE", "Prob. underest.", "Coverage")
    )
  ) |>
  ggplot(aes(value, type, fill = method)) +
  geom_bar(stat = "identity", position = "dodge", width = 0.7) +
  geom_vline(
    data = tibble(
      metric = factor(c("Mean bias", "RMSE", "Prob. underest.", "Coverage")),
      value = c(0, 0, 0.5, 0.95),
    ),
    aes(xintercept = value),
    linetype = "dashed"
  ) +
  scale_fill_manual(values = mycols) +
  ggh4x::facet_nested(rel + dist ~ metric, scales = "free", space = "free_y") +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_continuous(),
      scale_x_continuous(),
      scale_x_continuous(limits = c(0.35, 0.65), labels = scales::percent),
      scale_x_continuous(limits = c(0.7, 1), labels = scales::percent)
    )
  ) +
  theme_bw() +
  labs(x = NULL, y = NULL, fill = NULL, title = "n = 50") +
  scale_y_discrete(labels = rev(c(
    expression(Theta),
    expression(Psi),
    expression(alpha)
  ))) +
  theme(legend.position = "bottom"); p_perf_n50_type


# Bias vs sample size
p_biassampsize_ours <-
  res_filtered |>
  filter(abs(bias) < 2) |>
  summarise(
    bias = mean(bias, trim = 0),
    .by = c(dist:method, type)
  ) |>
  mutate(
    n = as.numeric(factor(n)),
    rel = factor(rel, labels = c("Rel == 0.8", "Rel == 0.5"))
  ) |>
  ggplot(aes(n, bias, col = method)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggh4x::facet_nested(type ~ rel + dist, labeller = label_parsed, scales = "free_y") +
  scale_colour_manual(values = mycols) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Sample size (n)", y = "Rel. mean bias", col = NULL) +
  guides(colour = guide_legend(nrow = 1, reverse = TRUE, position = "top")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)); p_biassampsize_ours

res_filtered |>
  filter(abs(bias) < 2) |>
  summarise(
    bias = mean(bias, trim = 0),
    .by = c(dist:method, param)
  ) |>
  filter(param %in% growthpars) |>
  mutate(n = as.numeric(factor(n))) |>
  ggplot(aes(n, bias, col = method)) +
  geom_line(linewidth = 0.8) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  ggh4x::facet_nested(param ~ rel + dist) +
  scale_colour_manual(values = mycols) +
  scale_x_continuous(labels = c(15, 20, 50, 100, 1000)) +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Sample size (n)", y = "Rel. mean bias") +
  theme_bw()


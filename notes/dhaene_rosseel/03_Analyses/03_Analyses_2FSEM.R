##==============================================================================
## Project:     Resampling based bias correction for small sample SEM
##
## Script:      Create plots & tables 2FSEM
##==============================================================================

library(dplyr)
library(ggplot2)
library(ggh4x)

load("02_Results/2FSEM_est_combined_final.RData")

head(Results, n = 10)
str(Results)

# Point estimates
point.est <- Results %>%
  filter(ci.type %in% c("MLB", "JB", "BB-ci.norm", "Ozenne", "REML")) %>%
  select(rel, dist, seed, method, tryerror, convergence, nobs, 
         grep("est", colnames(Results))) %>%
  tidyr::gather(key = parameter, value = est, 8:ncol(.)) %>%
  mutate(parameter = gsub(pattern = "est_", x = parameter, "")) %>%
  mutate(rel = as.factor(rel),
         dist = as.factor(dist),
         method = as.factor(method),
         parameter = as.factor(parameter))

# Confidence intervals
cis.lower <- Results %>%
  select(rel, dist, seed, method, ci.type, tryerror, convergence, nobs,
         grep("ci.lower", colnames(Results))) %>%
  tidyr::gather(key = parameter, 
                value = ci.lower, 
                grep("ci.lower", colnames(.))) %>%
  mutate(parameter = gsub(pattern = "ci.lower_", x = parameter, ""))

cis.upper <- Results %>%
  select(rel, dist, seed, method, ci.type, tryerror, convergence, nobs, 
         grep("ci.upper", colnames(Results))) %>%
  tidyr::gather(key = parameter, 
                value = ci.upper, 
                grep("ci.upper", colnames(.))) %>%
  mutate(parameter = gsub(pattern = "ci.upper_", x = parameter, ""))

cis <- full_join(cis.lower, cis.upper) %>%
  mutate(rel = as.factor(rel),
         dist = as.factor(dist),
         method = as.factor(method),
         ci.type = as.factor(ci.type),
         parameter = as.factor(parameter))

# Combine point estimates and confidence intervals into one df
all <- full_join(point.est, cis)

# Population values
true.theta.rel50 <- c("fx=~x2" = 0.7,
                      "fx=~x3" = 0.6,
                      "fy=~y2" = 0.7,
                      "fy=~y3" = 0.6,
                      "fy~fx"  = 0.25,
                      "x1~~x1" = 1,
                      "x2~~x2" = 0.49,
                      "x3~~x3" = 0.36,
                      "y1~~y1" = 1,
                      "y2~~y2" = 0.49,
                      "y3~~y3" = 0.36,
                      "fx~~fx" = 1,
                      "fy~~fy" = 1,
                      "x1~1"   = 0,
                      "x2~1"   = 0,
                      "x3~1"   = 0,
                      "y1~1"   = 0,
                      "y2~1"   = 0,
                      "y3~1"   = 0)

true.theta.rel80 <- c("fx=~x2" = 0.7,
                      "fx=~x3" = 0.6,
                      "fy=~y2" = 0.7,
                      "fy=~y3" = 0.6,
                      "fy~fx"  = 0.25,
                      "x1~~x1" = 0.25,
                      "x2~~x2" = 0.1225,
                      "x3~~x3" = 0.09,
                      "y1~~y1" = 0.25,
                      "y2~~y2" = 0.1225,
                      "y3~~y3" = 0.09,
                      "fx~~fx" = 1,
                      "fy~~fy" = 1,
                      "x1~1"   = 0,
                      "x2~1"   = 0,
                      "x3~1"   = 0,
                      "y1~1"   = 0,
                      "y2~1"   = 0,
                      "y3~1"   = 0)

# Compute coverage
all.cov <- all %>%
  mutate(pop.value = ifelse(test = rel == "REL50", 
                            yes = true.theta.rel50[as.character(.$parameter)],
                            no = ifelse(test = rel == "REL80", 
                                        yes = true.theta.rel80[as.character(.$parameter)], 
                                        no = NA)),
         coverage = ifelse(test = pop.value >= ci.lower & pop.value <= ci.upper, 
                           yes = 1,
                           no = 0),
         est.bias = est - pop.value,
         est.rel.bias = ifelse(pop.value != 0, 
                               yes = (est-pop.value)/pop.value,
                               no = (est-pop.value)/1),
         parameter = factor(parameter,
                            levels = c("x1~~x1",
                                       "x2~~x2",
                                       "x3~~x3",
                                       "y1~~y1",
                                       "y2~~y2",
                                       "y3~~y3",
                                       "fx~~fx", 
                                       "fy~~fy", 
                                       "fy~fx",
                                       "fx=~x2",
                                       "fx=~x3",
                                       "fy=~y2",
                                       "fy=~y3",
                                       "x1~1",
                                       "x2~1",
                                       "x3~1",
                                       "y1~1",
                                       "y2~1",
                                       "y3~1"))) %>%
  arrange(rel, dist, ci.type, nobs, parameter, coverage) %>%
  mutate(it = rep(1:2000, times = length(unique(.$rel))*
                    length(unique(.$dist))*
                    length(unique(.$ci.type))*
                    length(unique(.$nobs))*
                    length(unique(.$parameter)))) %>%
  filter(convergence ==  1)

# Summaries
overview <- all.cov %>%
  filter(ci.type %in% c("MLB", "JB", "BB-ci.norm", "Ozenne")) %>%
  plyr::ddply(.data = .,
              .variables = c("method", 
                             "parameter", 
                             "rel", 
                             "dist",
                             "nobs"),
              .fun = summarize,
              nsim = length(unique(seed)),
              mean = mean(est, na.rm = TRUE),
              sd = sd(est, na.rm = TRUE),
              var = var(est, na.rm = TRUE),
              median.bias = median(est.bias, na.rm = TRUE),
              mean.bias = mean(est.bias, na.rm = TRUE),
              median.bias.rel = median(est.rel.bias, na.rm = TRUE),
              mean.bias.rel = mean(est.rel.bias, na.rm = TRUE),
              MSE = var + mean.bias^2,
              RMSE = sqrt(MSE),
              coverage.rate = mean(as.numeric(coverage) - 1, na.rm = TRUE))

#-------------------------------------------------------------------------------
# Trimmed version
#-------------------------------------------------------------------------------

splits <- all.cov %>%
  group_by(rel, dist, method, nobs, parameter) %>%
  group_split()

pct.trim <- 0.01

for(i in 1:length(splits)){
  
  quants <- quantile(splits[[i]]$est, 
                     probs = c(pct.trim/2, 1-pct.trim/2),
                     na.rm = TRUE)
  
  splits[[i]] <- splits[[i]] %>%
    filter(est > quants[1] & est < quants[2])
  
  print(paste0(i, " out of ", length(splits)))
  
}

trimmed <- do.call("rbind.data.frame", splits)

overview.trimmed <- trimmed %>%
  filter(ci.type %in% c("MLB", "JB", "BB-ci.norm", "Ozenne")) %>%
  plyr::ddply(.data = .,
              .variables = c("method", 
                             "parameter", 
                             "rel", 
                             "dist",
                             "nobs"),
              .fun = summarize,
              nsim = length(unique(seed)),
              mean = mean(est, na.rm = TRUE),
              sd = sd(est, na.rm = TRUE),
              var = var(est, na.rm = TRUE),
              median.bias = median(est.bias, na.rm = TRUE),
              mean.bias = mean(est.bias, na.rm = TRUE),
              median.bias.rel = median(est.rel.bias, na.rm = TRUE),
              mean.bias.rel = mean(est.rel.bias, na.rm = TRUE),
              MSE = var + mean.bias^2,
              RMSE = sqrt(MSE),
              coverage.rate = mean(as.numeric(coverage), na.rm = TRUE))

#-------------------------------------------------------------------------------
# Function to create table for median bias, mean bias, & RMSE
#-------------------------------------------------------------------------------

create.table <- function(results, reliability) {
  
  tmp <- results %>%
    filter(parameter %in% c("x1~~x1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2") & 
           dist %in% c("Normal", "NonNormal"))

  params <- as.character(unique(tmp$parameter))
  nobs <- c(15, 20, 50, 100, 1000)
  est.methods <- c("MLB", "JB", "BB", "Ozenne")
  distr <- c("Normal", "NonNormal")

  for(i in seq_along(params)){
    
    assign(params[i],
           matrix(NA,
                  nrow = length(est.methods),
                  ncol = length(nobs)*length(distr),
                  dimnames = list(est.methods, 
                                  rep(nobs, times = 2))))  
    
  }
  
  mean.bias.rel <- list(mget(params))[[1]]
  median.bias.rel <- list(mget(params))[[1]]
  RMSE <- list(mget(params))[[1]]
  
  # relative median bias
  for(i in seq_along(params)){
    
    for(j in seq_along(nobs)){
        
        for(k in seq_along(distr)){
          
          column <- tmp %>%
            filter(parameter == params[i] & nobs == nobs[j] & 
                     rel == reliability & dist == distr[k]) %>%
            pull(median.bias.rel)
        
          median.bias.rel[[params[i]]][ , (k-1)*length(nobs) + j] <- column
          
        }
      }
  }
    
  # relative mean bias
  for(i in seq_along(params)){
    
    for(j in seq_along(nobs)){
      
      for(k in seq_along(distr)){
        
        column <- tmp %>%
          filter(parameter == params[i] & nobs == nobs[j] & 
                   rel == reliability & dist == distr[k]) %>%
          pull(mean.bias.rel)
        
        mean.bias.rel[[params[i]]][ , (k-1)*length(nobs) + j] <- column
        
      }
    }
  }
  
  # RMSE
  for(i in seq_along(params)){
    
    for(j in seq_along(nobs)){
      
      for(k in seq_along(distr)){
        
        column <- tmp %>%
          filter(parameter == params[i] & nobs == nobs[j] & 
                   rel == reliability & dist == distr[k]) %>%
          pull(RMSE)
        
        RMSE[[params[i]]][ , (k-1)*length(nobs) + j] <- column
        
      }
    }
  }
  
  as.data.frame(cbind(do.call(rbind.data.frame, mean.bias.rel),
                      do.call(rbind.data.frame, median.bias.rel), 
                      do.call(rbind.data.frame, RMSE)))

}

#-------------------------------------------------------------------------------
# Create actual table
#-------------------------------------------------------------------------------

nontrimmed.rel80 <- create.table(overview, reliability = "REL80")
nontrimmed.rel50 <- create.table(overview, reliability = "REL50")

trimmed.rel80 <- create.table(overview.trimmed, reliability = "REL80")
trimmed.rel50 <- create.table(overview.trimmed, reliability = "REL50")

xtable::xtable(nontrimmed.rel80)
xtable::xtable(nontrimmed.rel50)
xtable::xtable(trimmed.rel80)
xtable::xtable(trimmed.rel50)

#-------------------------------------------------------------------------------
# Plots mean bias, median bias, RSME
#-------------------------------------------------------------------------------

sd.color.values <- c("MLB" = "#00306C",
                     "JB" = "#50a17d",
                     "BB" = "#67cfa2",
                     "Ozenne" = "#1E90FF")

sd.shape.values <- c("MLB" = 15, 
                     "JB" = 5, 
                     "BB" = 9, 
                     "Ozenne" = 20)

sd.linetype.values <- c("MLB" = 1,
                        "JB" = 12,
                        "BB" = 3, 
                        "Ozenne" = 7)

sd.legend.labels <- c("ML", 
                      "Jackknife", 
                      "Bootstrap", 
                      "Ozenne et al.")

plot.values <- c("mean.bias.rel", 
                 "median.bias.rel", 
                 "RMSE")

ylab.labels <- c("relative mean bias",
                 "relative median bias", 
                 "RMSE")

ylimits <- list(c(-0.25, 0.25),
                c(-0.40, 0.30), 
                c(-0.2, 2.5))

filenames <- c("2FSEM_mean_trimmed.png",
               "2FSEM_median_trimmed.png",
               "2FSEM_RMSE_trimmed.png")

for(i in seq_along(plot.values)) {
  
  p <- ggplot(data = overview.trimmed %>%
                mutate(parameter = factor(parameter,
                                          levels = levels(.$parameter),
                                          labels = c(expression(theta["1,1"]),
                                                     expression(theta["2,2"]),
                                                     expression(theta["3,3"]),
                                                     expression(theta["4,4"]),
                                                     expression(theta["5,5"]),
                                                     expression(theta["6,6"]),
                                                     expression(psi["1,1"]),
                                                     expression(psi["2,2"]),
                                                     expression(beta),
                                                     expression(lambda["2,1"]),
                                                     expression(lambda["3,1"]),
                                                     expression(lambda["5,2"]),
                                                     expression(lambda["6,2"]),
                                                     expression(nu["1"]),
                                                     expression(nu["2"]),
                                                     expression(nu["3"]),
                                                     expression(nu["4"]),
                                                     expression(nu["5"]),
                                                     expression(nu["6"]))),
                       dist = factor(dist,
                                     levels = c("Normal", 
                                                "Kurtosis", 
                                                "NonNormal"),
                                     labels = c("Normal", 
                                                "Kurtosis", 
                                                "Non-Normal")),
                       rel = factor(rel,
                                    levels = c("REL80", "REL50"),
                                    labels = c("Reliability = .80", 
                                               "Reliability = .50"))) %>%
                filter(parameter %in% levels(.$parameter)[c(1, 7, 8, 9, 10)] & 
                         dist %in% c("Normal", "Non-Normal")),
              aes(x = nobs, 
                  y = get(plot.values[i]), 
                  group = method, 
                  color = method)) +
    geom_hline(yintercept = 0, color = "snow3", linetype = "dashed") +
    geom_point(aes(color = method, shape = method), size = 1.3) +
    geom_line(aes(color = method, linetype = method), lwd = 0.5, alpha = 0.75) +
    ggh4x::facet_nested(parameter ~ rel + dist, scales = "fixed",
                        labeller = labeller(.rows = label_parsed)) +
    labs(title = "", x ="sample size", y = ylab.labels[i]) +
    theme_bw(16) +
    theme(text = element_text(size = 16), 
          panel.grid.major = element_line(size = 0.15, 
                                          linetype = 'solid',
                                          colour = "lightgray"), 
          panel.grid.minor = element_line(size = 0.10,
                                          linetype = 'solid',
                                          colour = "lightgray"),
          legend.position = "top", 
          legend.box = "horizontal", 
          legend.title = element_blank(),
          legend.key.width = unit(3, "line"),
          panel.spacing = unit(2, "mm"), 
          axis.text = element_text(size = 10, 
                                   angle = 0, 
                                   hjust = 0.5),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          strip.background = element_rect(color = "black",
                                          size = 0.5, 
                                          fill = "gray90")) +
    scale_linetype_manual(values = sd.linetype.values, 
                          labels = sd.legend.labels) +
    scale_color_manual(values = sd.color.values,
                       labels = sd.legend.labels) +
    scale_shape_manual(values = sd.shape.values,
                       labels = sd.legend.labels) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    coord_cartesian(ylim = ylimits[[i]])
  
  print(p)
  
  ggsave(filename = filenames[i], 
         plot = last_plot(),
         width = 30,
         height = 25,
         units = "cm",
         dpi = "retina")
  
}

# Extra plot: RMSE inflation relative to MLB
RMSE.inflation <- overview.trimmed %>%
  select(method, dist, rel, parameter, nobs, RMSE) %>%
  tidyr::spread(key = method, value = RMSE) %>%
  mutate(JB = JB/MLB,
         BB = BB/MLB,
         Ozenne = Ozenne/MLB, 
         MLB = 1) %>%
  tidyr::gather(key = method, value = RMSE.infl, 5:ncol(.))

p <- ggplot(data = RMSE.inflation %>%
              mutate(parameter = factor(parameter,
                                        levels = levels(.$parameter),
                                        labels = c(expression(theta["1,1"]),
                                                   expression(theta["2,2"]),
                                                   expression(theta["3,3"]),
                                                   expression(theta["4,4"]),
                                                   expression(theta["5,5"]),
                                                   expression(theta["6,6"]),
                                                   expression(psi["1,1"]),
                                                   expression(psi["2,2"]),
                                                   expression(beta),
                                                   expression(lambda["2,1"]),
                                                   expression(lambda["3,1"]),
                                                   expression(lambda["5,2"]),
                                                   expression(lambda["6,2"]),
                                                   expression(nu["1"]),
                                                   expression(nu["2"]),
                                                   expression(nu["3"]),
                                                   expression(nu["4"]),
                                                   expression(nu["5"]),
                                                   expression(nu["6"]))),
                     dist = factor(dist,
                                   levels = c("Normal", 
                                              "Kurtosis", 
                                              "NonNormal"),
                                   labels = c("Normal", 
                                              "Kurtosis", 
                                              "Non-Normal")),
                     rel = factor(rel,
                                  levels = c("REL80", "REL50"),
                                  labels = c("Reliability = .80", 
                                             "Reliability = .50"))) %>%
              filter(parameter %in% levels(.$parameter)[c(1, 7, 8, 9, 10)] & 
                       dist %in% c("Normal", "Non-Normal")),
            aes(x = nobs, 
                y = RMSE.infl, 
                group = method, 
                color = method)) +
  geom_hline(yintercept = 0, color = "snow3", linetype = "dashed") +
  geom_point(aes(color = method, shape = method), size = 1.3) +
  geom_line(aes(color = method, linetype = method), lwd = 0.5, alpha = 0.75) +
  ggh4x::facet_nested(parameter ~ rel + dist, scales = "fixed",
                      labeller = labeller(.rows = label_parsed)) +
  labs(title = "", x ="sample size", y = ylab.labels[i]) +
  theme_bw(16) +
  theme(text = element_text(size = 16), 
        panel.grid.major = element_line(size = 0.15, 
                                        linetype = 'solid',
                                        colour = "lightgray"), 
        panel.grid.minor = element_line(size = 0.10,
                                        linetype = 'solid',
                                        colour = "lightgray"),
        legend.position = "top", 
        legend.box = "horizontal", 
        legend.title = element_blank(),
        legend.key.width = unit(3, "line"),
        panel.spacing = unit(2, "mm"), 
        axis.text = element_text(size = 10, 
                                 angle = 0, 
                                 hjust = 0.5),
        panel.border = element_rect(color = "black", fill = NA, size = 0.5),
        strip.background = element_rect(color = "black",
                                        size = 0.5, 
                                        fill = "gray90")) +
  scale_linetype_manual(values = sd.linetype.values, 
                        labels = sd.legend.labels) +
  scale_color_manual(values = sd.color.values,
                     labels = sd.legend.labels) +
  scale_shape_manual(values = sd.shape.values,
                     labels = sd.legend.labels) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  scale_y_continuous(breaks = c(1, 2, 3))

print(p)

ggsave(filename = "2FSEM_inflation_trimmed.png", 
       plot = last_plot(),
       width = 30,
       height = 25,
       units = "cm",
       dpi = "retina")

#-------------------------------------------------------------------------------
# Coverage
#-------------------------------------------------------------------------------

distr <- c("Normal", "NonNormal")
rels <- c("REL80", "REL50")
params <- c("x1~~x1", "fx~~fx", "fy~~fy", "fy~fx", "fx=~x2") 
ci.types <- c("MLB", "JB", "BB-ci.norm", "BB-ci.basic", "BB-ci.stud",
              "BB-ci.perc", "BB-ci.bca", "Ozenne")
nobss <- c(15, 20, 50, 100, 1000)

for(a in seq_along(distr)){
  
  for(b in seq_along(rels)){
    
    tmp <- trimmed %>%
      group_by(ci.type, dist, rel, nobs, parameter) %>%
      filter(dist == distr[a] & rel == rels[b] & nobs %in% nobss &
               parameter %in% params & ci.type %in% ci.types) %>%
      mutate(parameter = factor(parameter, levels = params)) %>%
      summarise(coverage.rate = mean(coverage, na.rm = TRUE)) %>%
      tidyr::spread(key = nobs, value = coverage.rate) %>%
      arrange(parameter) 
    
    assign(paste0("Results_", distr[a], "_", rels[b]), tmp)
  }
  
}


df.full <- cbind(Results_Normal_REL80[ , -3], 
                 Results_Normal_REL50[ , -c(1:4)],
                 Results_NonNormal_REL80[ , -c(1:4)],
                 Results_NonNormal_REL50[ , -c(1:4)]) %>%
  relocate(ci.type, .after = parameter) %>%
  mutate_if(is.numeric, round, 2)

xtable::xtable(df.full)

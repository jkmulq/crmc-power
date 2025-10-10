# QUick sim
rm(list = ls())
set.seed(83141984)

# Directory
setwd("~/Documents/Dev/crmc-power/R")

# Packages
library(tidyverse)
library(extraDistr)
library(ggfixest)
source("0_functions.R") # Source functions

# Seed
set.seed(1235901350)

## 0 Preliminaries
# 0.1 simulation parameters
J <- 30 # No. agencies
t <- 12 * 4 # No. years
rho_vec <- seq(-0.025, 0, 0.002) # Treatment effect sizes
rho <- 0.3
M <- 100 # Simulation reps
Njmin <- 199
Njmax <- 200

# Detect available cores
n_cores <- parallel::detectCores() - 1
cat("Using", n_cores, "cores\n")

store <- vector("list", length = M)
store_se <- vector("list", length = M)

for (i in 1:M){
  
  # Generate agency data
  agency_data <- generate_agency_data(J = J, t = t,
                                      Njmin = Njmin, Njmax = Njmax,
                                      delta = rho)

  # Make data
  data <- generate_data_fast(agency_data, delta = rho)

  # Unique id
  data <- data %>% mutate(uid = paste(id, tj, sep = "_"))

  # Callaway-Sant'Anna
  twfe_cs <- staggered::staggered_cs(df = data,
                                      i = "uid", # Cross-sectional unit identifier
                                      t = "period", # Time period column
                                      g = "tj", # First period observation is treated,
                                      y = "y", # Outcome variable
                                      estimand = "simple")
  
  # Roth-Sant'Anna
  twfe_rs <- staggered::staggered(df = data,
                                  i = "uid", # Cross-sectional unit identifier
                                  t = "period", # Time period column
                                  g = "tj", # First period observation is treated,
                                  y = "y", # Outcome variable
                                  estimand = "simple")
  
  # Sun-Abraham
  twfe_sa <- staggered::staggered_sa(df = data,
                                  i = "uid", # Cross-sectional unit identifier
                                  t = "period", # Time period column
                                  g = "tj", # First period observation is treated,
                                  y = "y", # Outcome variable
                                  estimand = "simple")
  
  store[[i]] <- c(as.numeric(twfe_sa["estimate"]),
                  as.numeric(twfe_rs["estimate"]),
                  as.numeric(twfe_cs["estimate"]))
  
  store_se[[i]] <- as.numeric(twfe_rs["se_neyman"])
  
}

store %>% 
  do.call(rbind, . ) %>% 
  colMeans()

store %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  setNames(c("Sun-Abraham", 
             "Roth-Sant'Anna", 
             "Callaway-Sant'Anna")) %>% 
  pivot_longer(., cols = everything()) %>% 
  ggplot(data = ., aes(x = value, fill = name)) + 
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Treatment Effects") + 
  geom_vline(xintercept = rho)

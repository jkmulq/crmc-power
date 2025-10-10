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
M <- 250 # Simulation reps
Njmin <- 199
Njmax <- 200

# Detect available cores
n_cores <- parallel::detectCores() - 1
cat("Using", n_cores, "cores\n")

store <- vector("list", length = M)
store_se <- vector("list", length = M)

# Generate agency data
agency_data <- generate_agency_data(J = J, t = t,
                                    Njmin = Njmin, Njmax = Njmax,
                                    delta = rho)

# Remove unnecessary agency_data columns 
agency_data$tj <- NULL
agency_data$treated <- NULL
agency_data$pj_post <- NULL

# Make data
data <- generate_rs_data_fast(agency_data, delta = rho)

for (i in 1:M){
  
  # Draw treatment periods for agency
  tj_periods <- gtools::permute(13:(13 + J - 1))
  
  # Join treatment dates onto data
  tj_key <- data.frame(agency = 1:J, tj = tj_periods)
  data <- left_join(data, tj_key, by = "agency")  
  
  # Create true outcome
  data <- data %>% 
    mutate(y = ifelse(period >= tj, y_treated, y_untreat))
  
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
  
  store_se[[i]] <- as.numeric(twfe_rs["se"])
  
  # Remove treatment column
  data$tj <- NULL
  
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

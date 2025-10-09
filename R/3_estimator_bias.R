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
rho <- -0.05
M <- 100 # Simulation reps

# Average number of probationers over time for each agency
Njmin <- 20
Njmax <- 200
Nj <- extraDistr::rdunif(n = J,
                         min = Njmin, # Minimum number
                         max = Njmax) # Maximum number

# True recidivism rate across agencies
pj <- runif(n = J, 
            min = 0.2, # Minimum rate. Must be >= min treatment effect
            max = 0.5) # Maximum rate. 1 - this > max treatment

# 0.2 Simulation parameter summary
sim_params <- data.frame("j" = 1:J,
                         "Nj" = Nj,
                         "pj_pre" = pj)

# 0.2 Setup storage matrix for data
agency_data <- tibble(agency = rep(0, J * t))
agency_data$agency <- sort(rep(1:J, t)) # Setup agency index
agency_data$period <- rep(1:t, J) # Setup time index

# Join on simulation params
agency_data <- left_join(agency_data, sim_params, by = join_by("agency" == "j"))

# Create time-varying agency size
agency_data$Njt <- agency_data$Nj + extraDistr::rdunif(n = J * t, min = -15, max = 15)

# Detect available cores
n_cores <- parallel::detectCores() - 1
cat("Using", n_cores, "cores\n")

store <- vector("list", length = M)

# Create dummy treatment variable in agency_data. 
# Loop overwrites this with actual treatment period


for (i in 1:M){
  
  # Treatment period for each agency (add 12 to account for first year)
  tj <- sample(13:t, size = J, replace = FALSE)
  
  # Join treatment dates onto agency_data
  agency_data$tj <- rep(tj, each = t)
  
  # Create treatment period dummy
  agency_data <- agency_data %>% 
    mutate(treated = ifelse(period - tj >= 0, 1, 0))
  
  # Create post treatment recidivism probability
  agency_data <- agency_data %>% 
    mutate(pj_post = ifelse(treated == 1, pj - rho, pj))
  
  # Make data
  data <- generate_data_fast(agency_data, delta = rho)
  
  # Unique id
  data <- data %>% mutate(uid = paste(id, tj, agency, sep = "_"))
  
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

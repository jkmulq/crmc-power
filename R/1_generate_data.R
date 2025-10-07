### DGP
### Jack Mulqueeney
### Oct 7, 2025

# Clean environment
rm(list = ls())

# Packages
library(tidyverse)
library(extraDistr)

# Seed
set.seed(1235901350)

## 0 Preliminaries
# 0.1 simulation parameters
J <- 30 # No. agencies
t <- 4 * 12 # No. years
rho <- 0.1 # Treatment effect size
  
# Average number of probationers over time for each agency
Nj <- extraDistr::rdunif(n = J, 
                         min = 50, # Minimum number
                         max = 200) # Maximum number

# True recidivism rate across agencies
pj <- runif(n = J, 
            min = 0.2, # Minimum rate
            max = 0.5) # Maximum rate

# Treatment period for each agency (add 12 to account for first year)
tj <- sample(1:J, size = J, replace = FALSE) + 12

# 0.2 Simulation parameter summary
sim_params <- data.frame("j" = 1:J,
                         "Nj" = Nj,
                         "pj" = pj,
                         "tj" = tj)

# 0.3 Setup storage matrix for data
data <- tibble(y = rep(0, J * t))
data$agency <- sort(rep(1:J, t)) # Setup agency index
data$period <- rep(1:t, J) # Setup time index

# Join on simulation params
data <- left_join(data, sim_params, by = join_by("agency" == "j"))

# Create treatment period dummy
data <- data %>% 
  mutate(treated = ifelse(period - tj >= 0, 1, 0))

## 1 Simulate data ## 
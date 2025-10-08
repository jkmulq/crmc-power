### Power analysis
### Jack Mulqueeney
### Oct 8, 2025

# Clean environment
rm(list = ls())

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
t <- 4 * 12 # No. years
rho_vec <- seq(-0.02, 0, 0.0025) # Treatment effect sizes
M <- 200 # Simulation reps

# Average number of probationers over time for each agency
Nj <- extraDistr::rdunif(n = J, 
                         min = 50, # Minimum number
                         max = 200) # Maximum number

# True recidivism rate across agencies
pj <- runif(n = J, 
            min = 0.2, # Minimum rate. Must be >= min treatment effect
            max = 0.5) # Maximum rate. 1 - this > max treatment

# Treatment period for each agency (add 12 to account for first year)
tj <- sample(1:J, size = J, replace = FALSE) + 12

# 0.2 Simulation parameter summary
sim_params <- data.frame("j" = 1:J,
                         "Nj" = Nj,
                         "pj_pre" = pj,
                         "tj" = tj)

# 0.3 Setup storage matrix for data
agency_data <- tibble(agency = rep(0, J * t))
agency_data$agency <- sort(rep(1:J, t)) # Setup agency index
agency_data$period <- rep(1:t, J) # Setup time index

# Join on simulation params
agency_data <- left_join(agency_data, sim_params, by = join_by("agency" == "j"))

# Create treatment period dummy
agency_data <- agency_data %>% 
  mutate(treated = ifelse(period - tj >= 0, 1, 0))

# Create time-varying agency size
agency_data$Njt <- agency_data$Nj + extraDistr::rdunif(n = J * t, min = -15, max = 15)

# Detect available cores
n_cores <- parallel::detectCores() - 1
cat("Using", n_cores, "cores\n")

# Simulation result list 
results <- vector("list", length = length(rho_vec))

# Loop across rho values (serial outer loop, parallel inner loop)
for (e in seq_along(rho_vec)) {
  
  rho <- rho_vec[e]
  
  cat("Running rho =", rho, "(", e, "/", length(rho_vec), ")\n")
  
  # Run simulations and store results
  results[[e]] <- run_sim_for_rho(rho, M, agency_data, n_cores)
  
}


## Reshape data to nice format
results_df <- do.call(rbind, results) %>% 
  as_tibble() %>% 
  setNames(c("delta", "coef", "TWFE Clustered", "TWFE Robust")) %>% 
  pivot_longer(., cols = -c("delta", "coef"), values_to = "vcv")

## Generate rejection rates for H0: rho = delta0
delta0 <- 0
conf_lev <- 0.05
z <- qnorm(1 - conf_lev / 2)

# Rejection rate calculations
reject_rates <- results_df %>% 
  mutate(t = abs( (coef - delta0) / sqrt(vcv) ),
         reject = (t > z))

# Power
power_res <- reject_rates %>% 
  group_by(delta, name) %>% 
  summarise(power = mean(reject))

# Graph
ggplot(data = power_res, aes(x = delta, y = power, colour = name)) +
  geom_point() + 
  geom_smooth(aes(colour = name), 
              se = FALSE, span = 0.75) +
  theme_minimal() +
  labs(title = "Power Curve for TWFE", 
       x = expression(delta), 
       y = "Power",
       caption = str_wrap("Points represent frequency with which a two-tailed 5% t-test rejects the false H0: delta = 0. Smoothed line generated using LOESS."), 
       90)  +
  theme(
    plot.title = element_text(hjust = 0.5), 
    plot.caption = element_text(hjust = 0),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.direction = "horizontal"
  ) +
  ylim(-0.05, 1.05) +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = conf_lev, 
             linetype = "dashed", 
             colour = "firebrick") +
  geom_hline(yintercept = 0.8, 
             linetype = "dashed", 
             colour = "steelblue") +
  annotate(
    "text",
    x = mean(power_res$delta),       
    y = conf_lev + 0.05,            
    label = bquote(alpha == .(conf_lev)),
    colour = "firebrick",
    hjust = 1,                        
    size = 3.5
  ) +
  annotate(
    "text",
    x = max(power_res$delta),       
    y = 0.8 - 0.075,                  
    label = "80% power",
    colour = "steelblue",
    hjust = 1,
    size = 3.5
  )

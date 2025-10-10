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
t <- 12 * 4 # No. years
rho_vec <- seq(-0.025, 0, 0.0025) # Treatment effect sizes
M <- 100 # Simulation reps
Njmin <- 20 # Min average number of probationers in an agency
Njmax <- 200 # Min average number of probationers in an agency

# Detect available cores
n_cores <- parallel::detectCores() - 1
cat("Using", n_cores, "cores\n")

# Simulation result list 
results <- vector("list", length = length(rho_vec))

# Loop across rho values (serial outer loop, parallel inner loop)
for (e in seq_along(rho_vec)) {
  
  # Select treatment effect
  rho <- rho_vec[e]
  
  # Generate agency data
  agency_data <- generate_agency_data(J = J, t = t, 
                                      Njmin = Njmin, Njmax = Njmax,
                                      delta = rho)
  
  # Make data
  data <- generate_data_fast(agency_data, delta = rho)
  
  cat("Running rho =", rho, "(", e, "/", length(rho_vec), ")\n")
  
  # Run simulations and store results
  results[[e]] <- run_sim_for_rho(rho, M, agency_data, n_cores)
  
}

## Reshape data to nice format
results_df <- do.call(rbind, results) %>% 
  as.data.frame() %>% 
  setNames(c("delta", 
             "TWFE", "TWFE Clustered SE", "TWFE Robust SE", 
             "TWFE (drops)", "TWFE (drops) Clustered SE", "TWFE (drops) Robust SE",
             "Callaway-Sant'Anna", "Callaway-Sant'Anna SE")) %>% 
  pivot_longer(., cols = c("TWFE Clustered SE", 
                            "TWFE Robust SE", 
                            "TWFE (drops) Clustered SE",
                            "TWFE (drops) Robust SE",
                            "Callaway-Sant'Anna SE"), 
               values_to = "se") %>% 
  mutate(across(.cols = -name, .fns = as.numeric))

## Generate rejection rates for H0: rho = delta0
delta0 <- 0
conf_lev <- 0.05
z <- qnorm(1 - conf_lev / 2)

# Rejection rate calculations
reject_rates <- results_df %>% 
  mutate(across(.cols = c("TWFE", "TWFE (drops)", "Callaway-Sant'Anna"),
                .fns = ~abs( (.x - delta0) / se ) > z))

# Power
# TWFE
twfe_power <- reject_rates %>% 
  filter(str_detect(name, "TWFE"), !str_detect(name, "(drops)")) %>% 
  group_by(delta, name) %>% 
  summarise(power = mean(TWFE)) %>% 
  ungroup()

# TWFE drops
twfe_drops_power <- reject_rates %>% 
  filter(str_detect(name, "TWFE"), str_detect(name, "(drops)")) %>% 
  group_by(delta, name) %>% 
  summarise(power = mean(`TWFE (drops)`)) %>% 
  ungroup()

# C-S
cs_power <- reject_rates %>% 
  filter(str_detect(name, "Callaway")) %>% 
  group_by(delta, name) %>% 
  summarise(power = mean(`Callaway-Sant'Anna`)) %>% 
  ungroup()

# Power combined
power_res <- rbind(twfe_power, twfe_drops_power, cs_power)


# Graph
p <- ggplot(data = power_res, aes(x = abs(delta), y = power, colour = name)) +
  geom_point() + 
  geom_smooth(aes(colour = name), 
              se = FALSE, span = 0.75) +
  theme_minimal() +
  labs(title = "Power Curves", 
       x = expression("|" * delta * "|"), 
       y = "Power",
       caption = str_wrap("Points represent frequency with which a two-tailed 5% t-test rejects the false H0: delta = 0. Smoothed line generated using LOESS.", 90))  +
  theme(
    plot.title = element_text(hjust = 0.5), 
    plot.caption = element_text(hjust = 0),
    legend.position = "bottom",
    legend.title = element_blank()
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
    x = 0.015,       
    y = conf_lev + 0.05,            
    label = paste0("alpha == ", conf_lev),
    parse = TRUE,
    colour = "firebrick",
    hjust = 1,                        
    size = 3.5
  ) +
  annotate(
    "text",
    x = 0.005,       
    y = 0.8 - 0.075,                  
    label = "80% power",
    colour = "steelblue",
    hjust = 1,
    size = 3.5
  ) +
  guides(
    colour = guide_legend(
      nrow = 3,
      ncol = 2,
      byrow = TRUE
    )
  )

# Print graph
p

## Save results and graph
results_loc <- "/Users/jmulqueeney/Documents/Dev/crmc-power/results"
graphs_loc <- "/Users/jmulqueeney/Documents/Dev/crmc-power/graphs"

write.csv(results_df, paste0(results_loc, 
                             "/power results (M", M, "J", J, "T", t, "Njmax", 500, ").csv"))

ggsave(paste0(graphs_loc,
               "/power curves (M", M, "J", J, "T", t, "Njmax", Njmax, ").png"),
        p)
ggsave(paste0(graphs_loc,
              "/power curves (M", M, "J", J, "T", t, "Njmax", Njmax, ").svg"),
       p)
ggsave(paste0(graphs_loc,
              "/power curves (M", M, "J", J, "T", t, "Njmax", Njmax, ").pdf"),
       p)

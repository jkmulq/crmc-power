### DGP
### Jack Mulqueeney
### Oct 7, 2025

# Clean environment
rm(list = ls())

# Packages
library(tidyverse)
library(extraDistr)

## 0 Preliminaries
# 0.1 simulation parameters
J <- 30 # No. agencies
t <- 4 * 12 # No. years

# Average number of probationers over time for each agency
Nj <- extraDistr::rdunif(n = J, 
                         min = 50, # Minimum number
                         max = 200) # Maximum number


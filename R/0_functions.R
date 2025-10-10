## Functions script
devtools::install_github("jonathandroth/staggered")
library(staggered) # For Roth & Sant'Anna (2023) Estimator
library(fixest)
library(data.table)

generate_data <- function(agency_data, delta){
  # Function that generates the probation-agency-time data with treatment effect
  # Inputs are a special agency_data dataframe, which contains all information
  # about agency parameters over time.
  # Other is delta = treatment effect
  
  # Infer J (no. of agencies)
  J <- max(agency_data$agency)
  
  # Create post-treatment probabilities conditional on true treatment effect
  agency_data <- agency_data %>% 
    mutate(pj_post = ifelse(treated == 1, pj + delta, pj))
  
  # Create vector for storage
  data <- vector("list", length = t) 

  
  for (time in 1:t){
    
    # Filter for time period t
    agency_data_t <- agency_data %>% 
      filter(period == time)
    
    # For each agency, simulate based on the parameters in that time period
    outcome_data_t <- lapply(1:J, function(j){
      
      # Select an agency's data
      df <- agency_data_t %>% 
        filter(agency == j)
      
      # Simulate outcomes. 
      # This is the key simulation step, 
      # where we have different outcomes across simulations for the same people
      outcome_jt <- rbern(n = df$Njt, prob = df$pj_post)
      
      # Outcome DF
      tibble(y = outcome_jt, 
             id = 1:(df$Njt), 
             agency = j, 
             period = time, 
             treated = df$treated)
      
    }) %>% 
      do.call(rbind, .)
    
    # Append data to list
    data[[time]] <- outcome_data_t
    
  }
  
  # Combine into normal dataframe
  data <- do.call(rbind, data)
  
  # Join on simulation parameters
  data <- left_join(data, agency_data, by = c("agency", "period")) 
  
  # Check that treatment years align
  check <- all(data$treated.x == data$treated.y)
  print(ifelse(check, "Treatment years align", stop()))
  
  # Check that group years align
  data <- data %>% 
    group_by(agency, period) %>% 
    mutate(Njt_imputed = n()) %>% 
    ungroup()
  
  check <- all(data$Njt == data$Njt_imputed)
  print(ifelse(check, "Number of probationers per period align", stop()))
  
  # Create event year variable
  data <- data %>% 
    mutate(e_period = period - tj)
  
  # Return
  return(data)
  
}


generate_agency_data <- function(J, t, Njmin, Njmax, delta){

  # Average number of probationers over time for each agency
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
  
  # Treatment period for each agency
  # tj <- sample(13:t, size = J, replace = FALSE)
  tj <- gtools::permute(13:t)[1:J]
  
  # Join treatment dates onto agency_data
  tj_key <- data.frame(agency = 1:J, tj = tj)
  agency_data <- left_join(agency_data, tj_key, by = "agency")
  
  # Create treatment period dummy
  agency_data <- agency_data %>% 
    mutate(treated = ifelse(period - tj >= 0, 1, 0))
  
  # Create post treatment recidivism probability
  agency_data <- agency_data %>% 
    mutate(pj_post = ifelse(treated == 1, pj_pre + delta, pj_pre))
  
  # Return data frame
  return(agency_data)
  
}

generate_data_fast <- function(agency_data, delta) {
  
  # infer J and t
  J <- max(agency_data$agency)
  t <- max(agency_data$period)
  
  # update post-treatment probabilities once
  agency_data$pj_post <- ifelse(agency_data$treated == 1, 
                                agency_data$pj_pre + delta, 
                                agency_data$pj_pre)
  
  # preallocate a list for results
  data_list <- vector("list", length = nrow(agency_data))
  
  # simulate all outcomes at once using vectorized rbern
  # (one row = one agency-period cell)
  for (i in seq_len(nrow(agency_data))) {
    df <- agency_data[i, ]
    y <- rbern(df$Njt, df$pj_post)
    data_list[[i]] <- data.frame(
      y = y,
      id = seq_len(df$Njt),
      agency = df$agency,
      period = df$period,
      treated = df$treated
    )
  }
  
  # bind results
  data <- data.table::rbindlist(data_list, use.names = TRUE)
  
  # join parameters using a key merge (faster than left join)
  agency_data_dt <- data.table::as.data.table(agency_data)
  data <- merge(data, agency_data_dt, by = c("agency", "period"), all.x = TRUE)
  
  # checks (can be disabled for production)
  stopifnot(all(data$treated.x == data$treated.y))
  stopifnot(all(data$Njt == ave(data$agency, data$agency, data$period, FUN = length)))
  
  # Return data
  data[]
}


# Function to run simulations for one rho
run_sim_for_rho <- function(rho, M, agency_data, n_cores) {
  store <- mclapply(1:M, function(i) {
    
    # Generate data
    data <- generate_data_fast(agency_data = agency_data, delta = rho)
    
    # Estimate treatment effects
    # Base TWFE
    twfe_reg <- fixest::feols(y ~ i(treated.x) | agency + period,
                          data = data, lean = FALSE, mem.clean = TRUE)
    
    
    # Estimate clustered and HC robust standard errors
    # Note: vcov = "hetero" calculates HC1 errors, equivalent to STATA's vce(robust) option.
    # vcov() and vcov_cluster() estimate variances, so sqrt() to get SE.
    twfe_cl <- sqrt(fixest::vcov_cluster(twfe_reg, cluster = ~agency))
    twfe_rb <- sqrt(stats::vcov(twfe_reg, "hetero"))
    
    # TWFE, dropping 12 months after treatment
    data <- data %>% mutate(eperiod = period - tj)
    twfe_reg_drops <- fixest::feols(y ~ i(treated.x) | agency + period,
                              data = data %>% filter(!(eperiod %in% 0:11)), # Drop 12 months after initial treatment 
                              lean = FALSE, mem.clean = TRUE)
    
    twfe_cl_drops <- sqrt(fixest::vcov_cluster(twfe_reg_drops, cluster = ~agency))
    twfe_rb_drops <- sqrt(stats::vcov(twfe_reg_drops, "hetero"))
    
    # Other estimators
    # Require a unique panel id. 
    # To generate, concatenate 'id' and 'tj' (which is the treatment cohort variable)
    data <- data %>% 
      mutate(uid = paste(id, tj, sep = "_"))
    
    # Callaway-Sant'Anna
    twfe_cs <- staggered::staggered_cs(df = data,
                                      i = "uid", # Cross-sectional unit identifier
                                      t = "period", # Time period column
                                      g = "tj", # First period observation is treated,
                                      y = "y", # Outcome variable
                                      estimand = "simple")
    
    # Return coefficients and standard errors
    out <- c(twfe_reg$coefficients[1], 
             twfe_cl, 
             twfe_rb,
             twfe_reg_drops$coefficients[1],
             twfe_cl_drops,
             twfe_rb_drops,
             twfe_cs[1:2])
    
    # Return
    return(out)
    
  }, mc.cores = n_cores, mc.set.seed = TRUE)
  
  store <- do.call(rbind, store)
  colnames(store) <- NULL # Remove column names
  store <- cbind(rep(rho, M), store) # Append the true treatment effect
  store <- as.data.frame(store) %>% 
    setNames(c("delta", 
               "twfe_est", "twfe_cl_se", "twfe_hc1_se", 
               "twfe_drops_est", "twfe_drops_cl_se", "twfe_drops_hc1_se", 
               "cs_est", "cs_se"))
  
  return(store)
}


# sim_params: data.frame with columns j (agency), Nj (size), pj_pre, tj
generate_panel_data <- function(sim_params, t, delta, balanced = TRUE, attrition_sd = 0) {
  # sim_params: data.frame with columns j, Nj, pj_pre, tj
  sim_dt <- as.data.table(sim_params)
  J <- nrow(sim_dt)
  
  # Create individuals for each agency
  # Step 1: create rows (agency, id)
  indiv_dt <- sim_dt[, .(id = seq_len(Nj)), by = .(j, Nj, pj_pre, tj)]
  
  # Step 2: expand to all periods
  panel <- indiv_dt[, .(period = seq_len(t)), by = .(j, id, Nj, pj_pre, tj)]
  setnames(panel, "j", "agency")
  
  # treated indicator (unit is treated from period >= tj)
  panel[, treated := as.integer(period >= tj)]
  
  # probability for this row (bounded to [0,1])
  panel[, p := pmin(pmax(pj_pre + delta * treated, 0), 1)]
  
  # simulate outcomes (vectorised)
  panel[, y := rbinom(.N, size = 1, prob = p)]
  
  # unique unit id across agencies
  panel[, uid := paste0(agency, "_", id)]
  
  # if you want to allow attrition (unbalanced panel), remove random draws:
  if (!balanced) {
    # attrition_sd: fraction sd of Nj to remove each period (simple example)
    set.seed(NULL)
    panel <- panel[, {
      n_keep <- max(1, round(Nj + rnorm(1, mean = 0, sd = attrition_sd * Nj)))
      keep_ids <- sample(seq_len(Nj), size = min(n_keep, Nj))
      .SD[id %in% keep_ids]
    }, by = .(agency, period)]
  }
  
  # return data.frame for downstream packages
  return(as.data.frame(panel))
}



generate_rs_data_fast <- function(agency_data, delta) {
  
  # infer J and t
  J <- max(agency_data$agency)
  t <- max(agency_data$period)
  
  # preallocate a list for results
  data_list <- vector("list", length = nrow(agency_data))
  
  # simulate all outcomes at once using vectorized rbern
  # (one row = one agency-period cell)
  for (i in seq_len(nrow(agency_data))) {
    df <- agency_data[i, ]
    y_untreat <- rnorm(df$Njt, mean = 0, sd = 1)
    y_treated <- rnorm(df$Njt, mean = 0 + delta, sd = 1)
    # y_untreat <- rbern(df$Njt, df$pj_pre)
    # y_treated <- rbern(df$Njt, df$pj_pre + delta)
    
    data_list[[i]] <- data.frame(
      y_untreat = y_untreat,
      y_treated = y_treated,
      id = seq_len(df$Njt),
      agency = df$agency,
      period = df$period
    )
  }
  
  # bind results
  data <- data.table::rbindlist(data_list, use.names = TRUE)
  
  # join parameters using a key merge (faster than left join)
  agency_data_dt <- data.table::as.data.table(agency_data)
  data <- merge(data, agency_data_dt, by = c("agency", "period"), all.x = TRUE)
  
  # create treatment probability column
  data$pj_post <- data$pj_pre + delta
  
  # checks (can be disabled for production)
  stopifnot(all(data$treated.x == data$treated.y))
  stopifnot(all(data$Njt == ave(data$agency, data$agency, data$period, FUN = length)))
  
  # Return data
  data[]
}


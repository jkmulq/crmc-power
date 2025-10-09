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
    
    
    # Other estimators
    # Require a unique panel id. 
    # To generate, concatenate 'id' and 'tj' (which is the treatment cohort variable)
    data <- data %>% 
      mutate(uid = paste(id, tj, sep = "_"))
    
    # Roth-Sant'Anna
    twfe_rs <- staggered::staggered(df = data,
                                          i = "uid", # Cross-sectional unit identifier
                                          t = "period", # Time period column
                                          g = "tj", # First period observation is treated,
                                          y = "y", # Outcome variable
                                          estimand = "simple")
    
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
             twfe_rs[1:2],
             twfe_cs[1:2])
    
    # Return
    return(out)
    
  }, mc.cores = n_cores)
  
  store <- do.call(rbind, store)
  colnames(store) <- NULL # Remove column names
  store <- cbind(rep(rho, M), store) # Append the treatment effect
  store <- as.data.frame(store) %>% 
    setNames(c("delta", "twfe_est", "twfe_cl_se", "twfe_hc1_se", 
               "rs_est", "rs_se",
               "cs_est", "cs_se"))
  
  return(store)
}

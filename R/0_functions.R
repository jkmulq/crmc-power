## Functions

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
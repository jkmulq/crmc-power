# QUick sim
store <- vector("list", length = M)

set.seed(83141984)

for (i in 1:M){
  
  # Make data
  data <- generate_data_fast(agency_data, delta = 0)
  
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
  
  
  # sunab <- sunab_reg_bin_pre_post$coeftable[2,1]
  rthsnt <- twfe_rs[1] %>% as.numeric
  calsnt <- twfe_cs[1] %>% as.numeric
  sunab <- twfe_sa[1] %>% as.numeric
  store[[i]] <- c(sunab, rthsnt, calsnt)
  
}

store %>% 
  do.call(rbind, . ) %>% 
  colMeans()

store %>% 
  do.call(rbind, .) %>% 
  as.data.frame() %>% 
  setNames(c("Sun-Abraham", "Roth-Sant'Anna", "Callaway-Sant'Anna")) %>% 
  pivot_longer(., cols = everything()) %>% 
  ggplot(data = ., aes(x = value, fill = name)) + 
  geom_density(alpha = 0.5) +
  theme_minimal() +
  labs(title = "Treatment Effects") + 
  geom_vline(xintercept = 0.05)

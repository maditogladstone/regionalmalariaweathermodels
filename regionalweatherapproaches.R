# A comparative analysis of approaches to modelling weather impact on malaria incidence in different climate regions 
total_s_t <- Sys.time()

packs_s_t <- Sys.time()

library(readxl)
library(tidyverse)
library(deSolve)
library(parallel)
library(furrr)
library(dplyr)
library(purrr)
library(stringr)

theme_weather_plots <- function(){
  
  theme_minimal() +
    theme(
      plot.title = element_text(size = 20, hjust = 0.5, face = "plain", colour = "black"),
      panel.border = element_rect(linewidth = 1, fill = NA, colour = "black"),
      strip.text.y = element_text(size = 20, colour = "black"),
      strip.text.x = element_text(size = 20, colour = "black"),
      axis.ticks = element_line(),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 5),  # Add extra space on the left margin to reduce squashing
      axis.text.x = element_text(size = 20, colour = "black"), 
      axis.text.y = element_text(size = 20, colour = "black"),
      axis.title.x = element_text(size = 20, colour = "black"),
      axis.title.y = element_text(size = 20, colour = "black"),
      legend.title = element_text(size = 20, colour = 'black'),
      legend.text = element_text(size = 20, colour = 'black'),
      legend.position = "bottom",
      legend.key.height = unit(1, "cm"),
      coord_fixed(ratio = 2)
    )
}

packs_e_t <- Sys.time()

proc_s_t <- Sys.time()

initial_states <- read_excel("vector_host_states_parameters.xlsx", sheet = "states") 

wide_states <- initial_states %>% 
  select(!c(description, source)) %>% 
  pivot_wider(names_from = "state", names_vary = "slowest", values_from = starts_with("patch")) 

vec_states <- wide_states %>% 
  as.matrix() %>% 
  as.vector()

names(vec_states) <- names(wide_states)

# parameters
model_parameters <- read_excel("vector_host_states_parameters.xlsx", sheet = "parameters") 

wide_parameters <- model_parameters %>% 
  select(!c(description, source)) %>%
  pivot_wider(names_from = "parameter", names_vary = "slowest", values_from = starts_with("patch")) 

vec_parameters <- wide_parameters %>% 
  as.matrix() %>% 
  as.vector()

names(vec_parameters) <- names(wide_parameters)

# temperature
temps_values <- bind_rows(
  read_excel("era5-x0.25_timeseries_pr,tas_timeseries_monthly_1950-2022_mean_historical_era5_x0.25_kano.xlsx", sheet = "tas"), # Kano State, Nigeria
  read_excel("era5-x0.25_timeseries_pr,tas_timeseries_monthly_1950-2022_mean_historical_era5_x0.25_benisha.xlsx", sheet = "tas"), # Benishangul Gemz, Ethiopia
  read_excel("era5-x0.25_timeseries_pr,tas_timeseries_monthly_1950-2022_mean_historical_era5_x0.25_limpopo.xlsx", sheet = "tas") # Limpopo Province, South Africa
) %>% 
  select(!c("code")) %>% 
  pivot_longer(cols = !name, names_to = "time", values_to = "value") %>% 
  mutate(
    year = year(ym(time)),
    month = month(ym(time)),
    day = day(ym(time)),
    time = as.numeric(ym(time)-min(ym(time))),
    patch = case_when(
      name == "Kano" ~ "patch_1",
      name == "Benishangul Gumz" ~ "patch_2",
      name == "Limpopo" ~ "patch_3"
    ),
    name = case_when(
      name == "Kano" ~ "Sub-Tropical",
      name == "Benishangul Gumz" ~ "Tropical",
      name == "Limpopo" ~ "Semi-Arid"
    )
  ) 

wide_temps <- temps_values %>% 
  select(patch, time, value) %>% 
  pivot_wider(names_from = "patch", values_from = "value") 

approx_temps_values <- bind_cols(
  time = seq(0, 365*nrow(wide_temps)/12, 1), 
  apply(select(wide_temps, !time), 2, function(col){
    approx(x = unique(wide_temps$time), y = col, xout = seq(0, 365*nrow(wide_temps)/12, 1), method = "linear", rule = 2)$y
  }) # approx. values
)

tempfunc <- apply(select(approx_temps_values, !time), 2, function(column){ 
  approxfun (approx_temps_values$time, column, rule = 2) 
})

approx_temps_values_long <- approx_temps_values %>% 
  as.data.frame() %>% 
  pivot_longer(cols = !time, names_to = "patch", values_to = "value") %>%
  mutate(
    name = case_when(
      patch == "patch_1" ~ "Sub-Tropical",
      patch == "patch_2" ~ "Tropical",
      patch == "patch_3" ~ "Semi-Arid"
    )
  )

# rainfall
rains_values <- bind_rows(
  read_excel("era5-x0.25_timeseries_pr,tas_timeseries_monthly_1950-2022_mean_historical_era5_x0.25_kano.xlsx", sheet = "pr"), # Kano State, Nigeria
  read_excel("era5-x0.25_timeseries_pr,tas_timeseries_monthly_1950-2022_mean_historical_era5_x0.25_benisha.xlsx", sheet = "pr"), # Benishangul Gemz, Ethiopia
  read_excel("era5-x0.25_timeseries_pr,tas_timeseries_monthly_1950-2022_mean_historical_era5_x0.25_limpopo.xlsx", sheet = "pr") # Limpopo Province, South Africa
) %>% 
  select(!c("code")) %>% 
  pivot_longer(cols = !name, names_to = "time", values_to = "value") %>% 
  mutate(
    year = year(ym(time)),
    month = month(ym(time)),
    day = day(ym(time)),
    time = as.numeric(ym(time)-min(ym(time))),
    patch = case_when(
      name == "Kano" ~ "patch_1",
      name == "Benishangul Gumz" ~ "patch_2",
      name == "Limpopo" ~ "patch_3"
    ),
    name = case_when(
      name == "Kano" ~ "Sub-Tropical",
      name == "Benishangul Gumz" ~ "Tropical",
      name == "Limpopo" ~ "Semi-Arid"
    )
  )

wide_rains <- rains_values %>% 
  select(patch, time, value) %>% 
  pivot_wider(names_from = "patch", values_from = "value") 

approx_rains_values <- bind_cols(
  time = seq(0, 365*nrow(wide_rains)/12, 1), 
  apply(select(wide_rains, !time), 2, function(col){
    approx(x = unique(wide_rains$time), y = col, xout = seq(0, 365*nrow(wide_rains)/12, 1), method = "linear", rule = 2)$y
  }) # approx. values
)

rainfunc <- apply(select(approx_rains_values, !time), 2, function(column){ 
  approxfun (approx_rains_values$time, column, rule = 2) 
})

approx_rains_values_long <- approx_rains_values %>% 
  as.data.frame() %>% 
  pivot_longer(cols = !time, names_to = "patch", values_to = "value") %>%
  mutate(
    name = case_when(
      patch == "patch_1" ~ "Sub-Tropical",
      patch == "patch_2" ~ "Tropical",
      patch == "patch_3" ~ "Semi-Arid"
    )
  ) 

proc_e_t <- Sys.time()

vectorhostModel <- list(
  "Approach_A" = function(times, states, parameters, no_patches, sigma){
    with(as.list(c(states, parameters)), {
      
      # weather values
      temperature <- sapply(tempfunc, function(x){x(times)})
      rainfall <- sapply(rainfunc, function(x){x(times)})
      
      # states
      E_a = matrix(states[paste0("patch_",1:no_patches,"_","E_a")], ncol = 1)
      L_a = matrix(states[paste0("patch_",1:no_patches,"_","L_a")], ncol = 1)
      P_a = matrix(states[paste0("patch_",1:no_patches,"_","P_a")], ncol = 1)
      S_m = matrix(states[paste0("patch_",1:no_patches,"_","S_m")], ncol = 1)
      E_m = matrix(states[paste0("patch_",1:no_patches,"_","E_m")], ncol = 1)
      I_m = matrix(states[paste0("patch_",1:no_patches,"_","I_m")], ncol = 1)
      S_h = matrix(states[paste0("patch_",1:no_patches,"_","S_h")], ncol = 1)
      E_h = matrix(states[paste0("patch_",1:no_patches,"_","E_h")], ncol = 1)
      A_h = matrix(states[paste0("patch_",1:no_patches,"_","A_h")], ncol = 1)
      Iu_h = matrix(states[paste0("patch_",1:no_patches,"_","Iu_h")], ncol = 1)
      Is_h = matrix(states[paste0("patch_",1:no_patches,"_","Is_h")], ncol = 1)
      Tu_h = matrix(states[paste0("patch_",1:no_patches,"_","Tu_h")], ncol = 1)
      Ts_h = matrix(states[paste0("patch_",1:no_patches,"_","Ts_h")], ncol = 1)
      R_h = matrix(states[paste0("patch_",1:no_patches,"_","R_h")], ncol = 1)
      Icov = matrix(states[paste0("patch_",1:no_patches,"_","Icov")], ncol = 1)
      
      # parameters
      K_e <- matrix(parameters[paste0("patch_",1:no_patches,"_","K_e")], ncol = 1)
      n_e <- matrix(parameters[paste0("patch_",1:no_patches,"_","n_e")], ncol = 1)
      theta <- matrix(parameters[paste0("patch_",1:no_patches,"_","theta")], ncol = 1)
      kappa_e <- matrix(parameters[paste0("patch_",1:no_patches,"_","kappa_e")], ncol = 1)
      mu_e <- matrix(parameters[paste0("patch_",1:no_patches,"_","mu_e")], ncol = 1)
      kappa_l <- matrix(parameters[paste0("patch_",1:no_patches,"_","kappa_l")], ncol = 1)
      mu_l <- matrix(1/(-4.4 + 1.31*(temperature + 2) - 0.03*(temperature + 2)^2), ncol = 1)
      kappa_p <- matrix(parameters[paste0("patch_",1:no_patches,"_","kappa_p")], ncol = 1)
      mu_p <- matrix(parameters[paste0("patch_",1:no_patches,"_","mu_p")], ncol = 1)    
      mu_m <- matrix((3.04 + 29.564*exp(-(temperature + 273.15 - 278)/2.7035))/30.4, ncol = 1)
      gamma_m <- matrix(parameters[paste0("patch_",1:no_patches,"_","gamma_m")], ncol = 1)
      a <- 1.25*matrix(1/(107.204 - 13.3523*temperature + 0.677509*temperature^2 - 0.0159732*temperature^3 + 0.000144876*temperature^4), ncol = 1)
      b <- matrix(parameters[paste0("patch_",1:no_patches,"_","b")], ncol = 1)
      c <- matrix(parameters[paste0("patch_",1:no_patches,"_","c")], ncol = 1)
      mu_h <- matrix(parameters[paste0("patch_",1:no_patches,"_","mu_h")], ncol = 1)
      mu_s <- matrix(parameters[paste0("patch_",1:no_patches,"_","mu_s")], ncol = 1)
      rho <- matrix(parameters[paste0("patch_",1:no_patches,"_","rho")], ncol = 1)
      gamma_h <- matrix(parameters[paste0("patch_",1:no_patches,"_","gamma_h")], ncol = 1)
      pa <- matrix(parameters[paste0("patch_",1:no_patches,"_","pa")], ncol = 1)
      omega <- matrix(parameters[paste0("patch_",1:no_patches,"_","omega")], ncol = 1)
      delta_r <- matrix(parameters[paste0("patch_",1:no_patches,"_","delta_r")], ncol = 1)
      delta_u <- matrix(parameters[paste0("patch_",1:no_patches,"_","delta_u")], ncol = 1)
      delta_s <- matrix(parameters[paste0("patch_",1:no_patches,"_","delta_s")], ncol = 1)
      nu <- matrix(parameters[paste0("patch_",1:no_patches,"_","nu")], ncol = 1)
      tau_u <- matrix(parameters[paste0("patch_",1:no_patches,"_","tau_u")], ncol = 1)
      tau_s <- matrix(parameters[paste0("patch_",1:no_patches,"_","tau_s")], ncol = 1)
      alpha <- matrix(parameters[paste0("patch_",1:no_patches,"_","alpha")], ncol = 1)
      zeta_a <- matrix(parameters[paste0("patch_",1:no_patches,"_","zeta_a")], ncol = 1)
      zeta_u <- matrix(parameters[paste0("patch_",1:no_patches,"_","zeta_u")], ncol = 1)
      zeta_s <- matrix(parameters[paste0("patch_",1:no_patches,"_","zeta_s")], ncol = 1)
      
      # total populations
      tot_m = S_m + E_m + I_m
      tot_h = S_h + E_h + A_h + Iu_h + Is_h + Tu_h + Ts_h + R_h
      
      # infectious humans 
      Infectious <- zeta_a*A_h + Iu_h + Is_h + zeta_u*Tu_h + zeta_s*Ts_h
      
      # foi
      eff_m <- 1 # prevents max. 50% of bites on infectious humans
      eff_h <- 1 # prevents max. 50% of bites by infectious mosquitoes
      lambda_m <- (a*b/(1 + eff_m*Icov))*(Infectious/tot_h)
      lambda_h <- (a*c/(1 + eff_h*Icov))*(tot_m/tot_h)*(I_m/tot_m)
      
      # aquatic mosquitoes
      dE_a = n_e*theta*pmax(0, 1 - tot_m/K_e)*tot_m - kappa_e*E_a - mu_e*E_a
      dL_a = kappa_e*E_a - kappa_l*L_a - mu_l*L_a
      dP_a = kappa_l*L_a - kappa_p*P_a - mu_p*P_a
      
      # adult mosquitoes
      dS_m = kappa_p*P_a - lambda_m*S_m - mu_m*S_m
      dE_m = lambda_m*S_m - gamma_m*E_m - mu_m*E_m
      dI_m = gamma_m*E_m - mu_m*I_m
      
      # human population
      dS_h = mu_h*tot_h + rho*R_h - lambda_h*S_h - mu_h*S_h
      dE_h = lambda_h*S_h - gamma_h*E_h - mu_h*E_h
      dA_h = pa*gamma_h*E_h + omega*Iu_h - delta_r*A_h - mu_h*A_h
      dIu_h = (1 - pa)*gamma_h*E_h + alpha*Is_h - omega*Iu_h - nu*Iu_h - delta_r*Iu_h - tau_u*Iu_h - mu_h*Iu_h
      dIs_h = nu*Iu_h - alpha*Is_h - tau_s*Is_h - mu_h*Is_h - 0*mu_s*Is_h
      dTu_h = tau_u*Iu_h - delta_u*Tu_h - mu_h*Tu_h
      dTs_h = tau_s*Is_h - delta_s*Ts_h - mu_h*Ts_h
      dR_h = delta_u*Tu_h + delta_s*Ts_h + delta_r*A_h + delta_r*Iu_h - rho*R_h - mu_h*R_h
      
      # intervention
      xi <- 1/365
      dIcov = sigma - xi*Icov
      
      # counts
      dCInc = (lambda_h*S_h/tot_h)*1000
      dCPrv = A_h + Iu_h + Is_h + Tu_h + Ts_h
      
      # output 
      output <- c(
        dE_a, dL_a, dP_a,
        dS_m, dE_m, dI_m, 
        dS_h, dE_h, dA_h, dIu_h, dIs_h, dTu_h, dTs_h, dR_h, 
        dIcov,
        dCInc, dCPrv
      )
      
      return(list(output))
    })
  },
  "Approach_B" = function(times, states, parameters, no_patches, sigma){
    with(as.list(c(states, parameters)), {
      
      # weather values
      temperature <- sapply(tempfunc, function(x){x(times)})
      rainfall <- sapply(rainfunc, function(x){x(times)})
      
      # states
      E_a = matrix(states[paste0("patch_",1:no_patches,"_","E_a")], ncol = 1)
      L_a = matrix(states[paste0("patch_",1:no_patches,"_","L_a")], ncol = 1)
      P_a = matrix(states[paste0("patch_",1:no_patches,"_","P_a")], ncol = 1)
      S_m = matrix(states[paste0("patch_",1:no_patches,"_","S_m")], ncol = 1)
      E_m = matrix(states[paste0("patch_",1:no_patches,"_","E_m")], ncol = 1)
      I_m = matrix(states[paste0("patch_",1:no_patches,"_","I_m")], ncol = 1)
      S_h = matrix(states[paste0("patch_",1:no_patches,"_","S_h")], ncol = 1)
      E_h = matrix(states[paste0("patch_",1:no_patches,"_","E_h")], ncol = 1)
      A_h = matrix(states[paste0("patch_",1:no_patches,"_","A_h")], ncol = 1)
      Iu_h = matrix(states[paste0("patch_",1:no_patches,"_","Iu_h")], ncol = 1)
      Is_h = matrix(states[paste0("patch_",1:no_patches,"_","Is_h")], ncol = 1)
      Tu_h = matrix(states[paste0("patch_",1:no_patches,"_","Tu_h")], ncol = 1)
      Ts_h = matrix(states[paste0("patch_",1:no_patches,"_","Ts_h")], ncol = 1)
      R_h = matrix(states[paste0("patch_",1:no_patches,"_","R_h")], ncol = 1)
      Icov = matrix(states[paste0("patch_",1:no_patches,"_","Icov")], ncol = 1)
      
      # parameters
      K_e <- matrix(parameters[paste0("patch_",1:no_patches,"_","K_e")], ncol = 1)
      n_e <- matrix(pmax(0, -0.61411*(temperature + 2)^3 + 38.93*(temperature + 2)^2 - 801.27*(temperature + 2) + 5391.4), ncol = 1)
      theta <- 0.1*matrix(0.00054*temperature^3 - 0.038*temperature^2 + 0.88*temperature + 1, ncol = 1)
      kappa_e <- 0.0001*matrix(pmax(0, -0.012*(temperature + 2)^3 + 0.81*(temperature + 2)^2 + 18*(temperature + 2) - 135.93), ncol = 1)
      mu_e <- matrix(pmax(0, 0.0033*(temperature + 2)^3 - 0.23*(temperature + 2)^2 + 5.3*(temperature + 2) - 40), ncol = 1)
      kappa_l <- 0.1*matrix(-0.002*(temperature + 2)^3 + 0.14*(temperature + 2)^2 - 3*(temperature + 2) + 22, ncol = 1)
      mu_l <- 0.1*matrix(0.00081*(temperature + 2)^3 - 0.056*(temperature + 2)^2 + 1.3*(temperature + 2) - 8.6, ncol = 1)
      kappa_p <- 0.1*matrix(pmax(0, 0.0018*(temperature + 2)^3 - 0.12*(temperature + 2)^2 + 2.7*(temperature + 2) - 20), ncol = 1)
      mu_p <- 0.001*matrix(-0.0034*(temperature + 2)^3 + 0.22*(temperature + 2)^2 + 4.9*(temperature + 2) + 34, ncol = 1)
      mu_m <- 0.001*matrix(-0.000091*temperature^3 + 0.059*temperature^2 + 1.3*temperature + 9.9, ncol = 1)
      gamma_m <- matrix(parameters[paste0("patch_",1:no_patches,"_","gamma_m")], ncol = 1)
      a <- 5*matrix(0.000203*(temperature^2 - 11.7*temperature)*sqrt(42.3 - temperature), ncol = 1)
      b <- matrix(parameters[paste0("patch_",1:no_patches,"_","b")], ncol = 1)
      c <- matrix(parameters[paste0("patch_",1:no_patches,"_","c")], ncol = 1)
      mu_h <- matrix(parameters[paste0("patch_",1:no_patches,"_","mu_h")], ncol = 1)
      mu_s <- matrix(parameters[paste0("patch_",1:no_patches,"_","mu_s")], ncol = 1)
      rho <- matrix(parameters[paste0("patch_",1:no_patches,"_","rho")], ncol = 1)
      gamma_h <- matrix(parameters[paste0("patch_",1:no_patches,"_","gamma_h")], ncol = 1)
      pa <- matrix(parameters[paste0("patch_",1:no_patches,"_","pa")], ncol = 1)
      omega <- matrix(parameters[paste0("patch_",1:no_patches,"_","omega")], ncol = 1)
      delta_r <- matrix(parameters[paste0("patch_",1:no_patches,"_","delta_r")], ncol = 1)
      delta_u <- matrix(parameters[paste0("patch_",1:no_patches,"_","delta_u")], ncol = 1)
      delta_s <- matrix(parameters[paste0("patch_",1:no_patches,"_","delta_s")], ncol = 1)
      nu <- matrix(parameters[paste0("patch_",1:no_patches,"_","nu")], ncol = 1)
      tau_u <- matrix(parameters[paste0("patch_",1:no_patches,"_","tau_u")], ncol = 1)
      tau_s <- matrix(parameters[paste0("patch_",1:no_patches,"_","tau_s")], ncol = 1)
      alpha <- matrix(parameters[paste0("patch_",1:no_patches,"_","alpha")], ncol = 1)
      zeta_a <- matrix(parameters[paste0("patch_",1:no_patches,"_","zeta_a")], ncol = 1)
      zeta_u <- matrix(parameters[paste0("patch_",1:no_patches,"_","zeta_u")], ncol = 1)
      zeta_s <- matrix(parameters[paste0("patch_",1:no_patches,"_","zeta_s")], ncol = 1)
      
      # total populations
      tot_m = S_m + E_m + I_m
      tot_h = S_h + E_h + A_h + Iu_h + Is_h + Tu_h + Ts_h + R_h
      
      # infectious humans 
      Infectious <- zeta_a*A_h + Iu_h + Is_h + zeta_u*Tu_h + zeta_s*Ts_h
      
      # foi
      eff_m <- 1 # prevents max. 50% of bites on infectious humans
      eff_h <- 1 # prevents max. 50% of bites by infectious mosquitoes
      lambda_m <- (a*b/(1 + eff_m*Icov))*(Infectious/tot_h)
      lambda_h <- (a*c/(1 + eff_h*Icov))*(tot_m/tot_h)*(I_m/tot_m)
      
      # aquatic mosquitoes
      dE_a = n_e*theta*pmax(0, 1 - tot_m/K_e)*tot_m - kappa_e*E_a - mu_e*E_a
      dL_a = kappa_e*E_a - kappa_l*L_a - mu_l*L_a
      dP_a = kappa_l*L_a - kappa_p*P_a - mu_p*P_a
      
      # adult mosquitoes
      dS_m = kappa_p*P_a - lambda_m*S_m - mu_m*S_m
      dE_m = lambda_m*S_m - gamma_m*E_m - mu_m*E_m
      dI_m = gamma_m*E_m - mu_m*I_m
      
      # human population
      dS_h = mu_h*tot_h + rho*R_h - lambda_h*S_h - mu_h*S_h
      dE_h = lambda_h*S_h - gamma_h*E_h - mu_h*E_h
      dA_h = pa*gamma_h*E_h + omega*Iu_h - delta_r*A_h - mu_h*A_h
      dIu_h = (1 - pa)*gamma_h*E_h + alpha*Is_h - omega*Iu_h - nu*Iu_h - delta_r*Iu_h - tau_u*Iu_h - mu_h*Iu_h
      dIs_h = nu*Iu_h - alpha*Is_h - tau_s*Is_h - mu_h*Is_h - 0*mu_s*Is_h
      dTu_h = tau_u*Iu_h - delta_u*Tu_h - mu_h*Tu_h
      dTs_h = tau_s*Is_h - delta_s*Ts_h - mu_h*Ts_h
      dR_h = delta_u*Tu_h + delta_s*Ts_h + delta_r*A_h + delta_r*Iu_h - rho*R_h - mu_h*R_h
      
      # intervention
      xi <- 1/365
      dIcov = sigma - xi*Icov
      
      # counts
      dCInc = (lambda_h*S_h/tot_h)*1000
      dCPrv = A_h + Iu_h + Is_h + Tu_h + Ts_h
      
      # output 
      output <- c(
        dE_a, dL_a, dP_a,
        dS_m, dE_m, dI_m, 
        dS_h, dE_h, dA_h, dIu_h, dIs_h, dTu_h, dTs_h, dR_h, 
        dIcov,
        dCInc, dCPrv
      )
      
      return(list(output))
    })
  },
  "Approach_C" = function(times, states, parameters, no_patches, sigma){
    with(as.list(c(states, parameters)), {
      
      # weather values
      temperature <- sapply(tempfunc, function(x){x(times)})
      rainfall <- sapply(rainfunc, function(x){x(times)})
      
      # states
      E_a = matrix(states[paste0("patch_",1:no_patches,"_","E_a")], ncol = 1)
      L_a = matrix(states[paste0("patch_",1:no_patches,"_","L_a")], ncol = 1)
      P_a = matrix(states[paste0("patch_",1:no_patches,"_","P_a")], ncol = 1)
      S_m = matrix(states[paste0("patch_",1:no_patches,"_","S_m")], ncol = 1)
      E_m = matrix(states[paste0("patch_",1:no_patches,"_","E_m")], ncol = 1)
      I_m = matrix(states[paste0("patch_",1:no_patches,"_","I_m")], ncol = 1)
      S_h = matrix(states[paste0("patch_",1:no_patches,"_","S_h")], ncol = 1)
      E_h = matrix(states[paste0("patch_",1:no_patches,"_","E_h")], ncol = 1)
      A_h = matrix(states[paste0("patch_",1:no_patches,"_","A_h")], ncol = 1)
      Iu_h = matrix(states[paste0("patch_",1:no_patches,"_","Iu_h")], ncol = 1)
      Is_h = matrix(states[paste0("patch_",1:no_patches,"_","Is_h")], ncol = 1)
      Tu_h = matrix(states[paste0("patch_",1:no_patches,"_","Tu_h")], ncol = 1)
      Ts_h = matrix(states[paste0("patch_",1:no_patches,"_","Ts_h")], ncol = 1)
      R_h = matrix(states[paste0("patch_",1:no_patches,"_","R_h")], ncol = 1)
      Icov = matrix(states[paste0("patch_",1:no_patches,"_","Icov")], ncol = 1)
      
      # survival of aquatic mosquitoes
      max_p_e <- 0.9
      max_p_l <- 0.75
      max_p_p <- 0.25
      R_l <- 50
      
      # parameters
      K_e <- matrix((1e2/0.01)*rainfall, ncol = 1)
      n_e <- matrix(parameters[paste0("patch_",1:no_patches,"_","n_e")], ncol = 1)
      theta <- matrix(parameters[paste0("patch_",1:no_patches,"_","theta")], ncol = 1)
      kappa_e <- matrix(pmin(1, (4*max_p_e/R_l^2)*rainfall*pmax(0, R_l - rainfall)), ncol = 1)
      mu_e <- matrix(parameters[paste0("patch_",1:no_patches,"_","mu_e")], ncol = 1)
      kappa_l <- matrix(exp(-1/(0.0557*(temperature + 2) - 0.06737))*pmin(1, (4*max_p_l/R_l^2)*rainfall*pmax(0, R_l - rainfall)), ncol = 1)
      mu_l <- matrix(0.0025*(temperature + 2)^2 - 0.094*(temperature + 2) + 1.0257, ncol = 1)
      kappa_p <- matrix(pmin(1, (4*max_p_p/R_l^2)*rainfall*pmax(0, R_l - rainfall)), ncol = 1)
      mu_p <- matrix(parameters[paste0("patch_",1:no_patches,"_","mu_p")], ncol = 1)
      mu_m <- matrix(1/(-4.40 + 1.31*temperature - 0.03*temperature^2), ncol = 1)
      gamma_m <- matrix(pmax(0, (temperature - 16)/111), ncol = 1)
      a <- 7.5*matrix(0.000203*(temperature^2 - 11.7*temperature)*sqrt(42.3 - temperature), ncol = 1)
      b <- matrix(parameters[paste0("patch_",1:no_patches,"_","b")], ncol = 1)
      c <- matrix(parameters[paste0("patch_",1:no_patches,"_","c")], ncol = 1)
      mu_h <- matrix(parameters[paste0("patch_",1:no_patches,"_","mu_h")], ncol = 1)
      mu_s <- matrix(parameters[paste0("patch_",1:no_patches,"_","mu_s")], ncol = 1)
      rho <- matrix(parameters[paste0("patch_",1:no_patches,"_","rho")], ncol = 1)
      gamma_h <- matrix(parameters[paste0("patch_",1:no_patches,"_","gamma_h")], ncol = 1)
      pa <- matrix(parameters[paste0("patch_",1:no_patches,"_","pa")], ncol = 1)
      omega <- matrix(parameters[paste0("patch_",1:no_patches,"_","omega")], ncol = 1)
      delta_r <- matrix(parameters[paste0("patch_",1:no_patches,"_","delta_r")], ncol = 1)
      delta_u <- matrix(parameters[paste0("patch_",1:no_patches,"_","delta_u")], ncol = 1)
      delta_s <- matrix(parameters[paste0("patch_",1:no_patches,"_","delta_s")], ncol = 1)
      nu <- matrix(parameters[paste0("patch_",1:no_patches,"_","nu")], ncol = 1)
      tau_u <- matrix(parameters[paste0("patch_",1:no_patches,"_","tau_u")], ncol = 1)
      tau_s <- matrix(parameters[paste0("patch_",1:no_patches,"_","tau_s")], ncol = 1)
      alpha <- matrix(parameters[paste0("patch_",1:no_patches,"_","alpha")], ncol = 1)
      zeta_a <- matrix(parameters[paste0("patch_",1:no_patches,"_","zeta_a")], ncol = 1)
      zeta_u <- matrix(parameters[paste0("patch_",1:no_patches,"_","zeta_u")], ncol = 1)
      zeta_s <- matrix(parameters[paste0("patch_",1:no_patches,"_","zeta_s")], ncol = 1)
      
      # total populations
      tot_m = S_m + E_m + I_m
      tot_h = S_h + E_h + A_h + Iu_h + Is_h + Tu_h + Ts_h + R_h
      
      # infectious humans 
      Infectious <- zeta_a*A_h + Iu_h + Is_h + zeta_u*Tu_h + zeta_s*Ts_h
      
      # foi
      eff_m <- 1 # prevents max. 50% of bites on infectious humans
      eff_h <- 1 # prevents max. 50% of bites by infectious mosquitoes
      lambda_m <- (a*b/(1 + eff_m*Icov))*(Infectious/tot_h)
      lambda_h <- (a*c/(1 + eff_h*Icov))*(tot_m/tot_h)*(I_m/tot_m)
      
      # aquatic mosquitoes
      dE_a = n_e*theta*pmax(0, 1 - tot_m/K_e)*tot_m - kappa_e*E_a - mu_e*E_a
      dL_a = kappa_e*E_a - kappa_l*L_a - mu_l*L_a
      dP_a = kappa_l*L_a - kappa_p*P_a - mu_p*P_a
      
      # adult mosquitoes
      dS_m = kappa_p*P_a - lambda_m*S_m - mu_m*S_m
      dE_m = lambda_m*S_m - gamma_m*E_m - mu_m*E_m
      dI_m = gamma_m*E_m - mu_m*I_m
      
      # human population
      dS_h = mu_h*tot_h + rho*R_h - lambda_h*S_h - mu_h*S_h
      dE_h = lambda_h*S_h - gamma_h*E_h - mu_h*E_h
      dA_h = pa*gamma_h*E_h + omega*Iu_h - delta_r*A_h - mu_h*A_h
      dIu_h = (1 - pa)*gamma_h*E_h + alpha*Is_h - omega*Iu_h - nu*Iu_h - delta_r*Iu_h - tau_u*Iu_h - mu_h*Iu_h
      dIs_h = nu*Iu_h - alpha*Is_h - tau_s*Is_h - mu_h*Is_h - 0*mu_s*Is_h
      dTu_h = tau_u*Iu_h - delta_u*Tu_h - mu_h*Tu_h
      dTs_h = tau_s*Is_h - delta_s*Ts_h - mu_h*Ts_h
      dR_h = delta_u*Tu_h + delta_s*Ts_h + delta_r*A_h + delta_r*Iu_h - rho*R_h - mu_h*R_h
      
      # intervention
      xi <- 1/365
      dIcov = sigma - xi*Icov
      
      # counts
      dCInc = (lambda_h*S_h/tot_h)*1000
      dCPrv = A_h + Iu_h + Is_h + Tu_h + Ts_h
      
      # output 
      output <- c(
        dE_a, dL_a, dP_a,
        dS_m, dE_m, dI_m, 
        dS_h, dE_h, dA_h, dIu_h, dIs_h, dTu_h, dTs_h, dR_h, 
        dIcov,
        dCInc, dCPrv
      )
      
      return(list(output))
    })
  }
)

sim_model <- function(mod_intervention, mod_approach){
  
  model_output <- ode(
    times = seq(0, 365*no_yrs, 1),
    y = vec_states,
    parms = vec_parameters,
    no_patches = no_patches,
    sigma = c(0.25, 0.45, 0.60, 0.75, 0.9)[mod_intervention]/365,
    func = vectorhostModel[[mod_approach]]
  )
  
  return(model_output)
}

num_cores <- detectCores()
cl <- makeCluster(num_cores)

# simulation and plotting times
no_yrs <- 72
start_year <- 1950 
skip_yrs <- 66 # number of years to exclude from plots
no_scenarios <- 5
no_patches <- 3

clusterExport(cl, varlist = c(
  "sim_model", "vectorhostModel", # model functions
  "no_yrs", "no_patches", # duration & areas
  "vec_states", "vec_parameters", # states & parameters
  "tempfunc", "rainfunc" # weather functions
))

clusterEvalQ(cl, c(library(deSolve), library(dplyr))) # ode solver & data manipulation

sim_s_t <- Sys.time()
result <- list(
  Approach_A = parSapply(cl, 1:no_scenarios, function(i) sim_model(mod_intervention = i, mod_approach = "Approach_A"), simplify = "array"),
  Approach_B = parSapply(cl, 1:no_scenarios, function(i) sim_model(mod_intervention = i, mod_approach = "Approach_B"), simplify = "array"),
  Approach_C = parSapply(cl, 1:no_scenarios, function(i) sim_model(mod_intervention = i, mod_approach = "Approach_C"), simplify = "array")
)
sim_e_t <- Sys.time()

stopCluster(cl)

pproc_s_t <- Sys.time()

convert_to_longer <- function(result, mod_approach){
  
  patch_names <- function(variable, suffix, no_patches){
    patches <- as.data.frame(
      map(1:no_patches, function(y){
        case_when(
          str_split_i(variable, pattern = "_", i = 2) == as.character(y) ~ paste0(suffix, "_", y)
        ) 
      }),
      col.names = paste0(suffix, "_", 1:no_patches)
    ) 
    
    new_names <- patches %>%
      mutate(
        patch = reduce(select(., starts_with(paste0(suffix, "_"))), coalesce) %>% factor()
      ) %>% select(patch) 
    
    return(new_names)
  }
  total_pops <- function(df_result, species, no_patches, no_scenario){
    
    pop_sums <- map(1:no_scenario, function(z){
      map(1:no_patches, function(y){
        map(species, function(x){
          select(df_result, starts_with(paste0("patch_", y)) & ends_with(paste0("_", x, ".", z))) %>% apply(., 1, sum) %>% 
            as.data.frame() %>% set_names(paste0("patch_", y))
        }) %>% bind_cols()
      }) %>% bind_cols()
    }) %>% bind_cols()
    
    names(pop_sums) <- unlist(map(1:no_scenario, function(z){
      map(1:no_patches, function(y){
        map(species, function(x){paste0("patch_", y, "_tot_", x, ".", z)
        }) 
      })
    })
    )
    
    return(pop_sums)
  }  
  cum_counts <- function(df_result, counts, no_patches, no_scenarios){
    
    diff_counts <- map(1:no_scenarios, function(z){map(1:no_patches, function(y){map(counts, function(x){c(0, diff(df_result[[paste0("patch_", y, "_C", x, ".", z)]]))
    }) %>% bind_cols()
    }) %>% bind_cols()
    }) %>% bind_cols()
    
    names(diff_counts) <- unlist(map(1:no_scenarios, function(z){
      map(1:no_patches, function(y){
        map(counts, function(x){paste0("patch_", y, "_", x, ".", z)
        })
      })
    })
    )
    
    return(diff_counts)
  }
  
  pproc <- result[[mod_approach]] %>% 
    as.data.frame() %>%
    mutate(
      bind_cols(total_pops(., c("a", "m", "h"), no_patches, 5)),
      bind_cols(cum_counts(., c("Inc", "Prv"), no_patches, 5))
    ) %>%
    pivot_longer(., cols = starts_with("time"), values_to = "time", names_to = "scenario_time") %>%
    select(!scenario_time) %>%
    pivot_longer(!time, names_to = "variable", values_to = "value") %>%
    mutate(
      scenario = as.factor(c("25%", "45%", "60%", "75%", "90%")[as.numeric(substring(variable, nchar(variable), nchar(variable)))]),
      species = case_when(
        substring(variable, nchar(variable) - nchar(".Z"), nchar(variable) - nchar(".Z")) == "a" ~ "aquatic",
        substring(variable, nchar(variable) - nchar(".Z"), nchar(variable) - nchar(".Z")) == "m" ~ "adult",
        substring(variable, nchar(variable) - nchar(".Z"), nchar(variable) - nchar(".Z")) == "h" ~ "host",
        substring(variable, nchar(variable) - nchar("YYY.Z"), nchar(variable) - nchar(".Z")) == "Icov" ~ "intervention",
        .default = "count"
      ), 
      state_var = substring(variable, nchar("patch_j_X"), nchar(variable) - nchar(".Z")),
      patch = as.factor(c("Semi-Arid", "Tropical", "Sub-Tropical")[as.numeric(substring(variable, nchar("patch_X"), nchar("patch_X")))]),
      approach = as.factor(mod_approach)
    )
  
  return(pproc)
}

output <- bind_rows(
  convert_to_longer(result, "Approach_A"),
  convert_to_longer(result, "Approach_B"),
  convert_to_longer(result, "Approach_C")
) 

incidence_summary <- function(df_result, mod_approach){
  
  cum_counts <- function(df_result, counts, no_patches, no_scenarios){
    
    diff_counts <- map(1:no_scenarios, function(z){map(1:no_patches, function(y){map(counts, function(x){c(0, diff(df_result[[paste0("patch_", y, "_C", x, ".", z)]]))
    }) %>% bind_cols()
    }) %>% bind_cols()
    }) %>% bind_cols()
    
    names(diff_counts) <- unlist(map(1:no_scenarios, function(z){
      map(1:no_patches, function(y){
        map(counts, function(x){paste0("patch_", y, "_", x, ".", z)
        })
      })
    })
    )
    
    return(diff_counts)
  }
  
  annual_values <- df_result[[mod_approach]] %>% 
    as.data.frame() %>% 
    mutate(
      bind_cols(cum_counts(., "Inc", no_patches, 5))
    ) %>% 
    select(!contains("CInc")) %>%
    select(time.1, contains("Inc")) %>%
    pivot_longer(!time.1, names_to = "variable", values_to = "value") %>%
    mutate(
      time = ceiling(time.1/365),
      patch = substring(variable, 1, nchar("patch_j")),
      scenario = as.factor(c("25%", "45%", "60%", "75%", "90%")[as.numeric(substring(variable, nchar(variable), nchar(variable)))]),
      approach = mod_approach,
      name = case_when(
        patch == "patch_1" ~ "Semi-Arid",
        patch == "patch_2" ~ "Tropical",
        patch == "patch_3" ~ "Sub-Tropical"
      )
    ) %>%
    group_by(time, name, scenario, approach) %>%
    summarize(
      annual_value = sum(value),
      .groups = "keep"
    ) 
  
  return(annual_values)
}

annual_incidence <- bind_rows(
  incidence_summary(result, "Approach_A"),
  incidence_summary(result, "Approach_B"),
  incidence_summary(result, "Approach_C")
)

pproc_e_t <- Sys.time()

plot_s_t <- Sys.time()

# temperature
monthly_temperature <- ggplot() +
  geom_point(data = filter(temps_values, time >= 365*skip_yrs, time <= 365*no_yrs), aes(x = time, y = value, colour = name), size = 4) +
  geom_line(data = filter(approx_temps_values_long, time >= 365*skip_yrs, time <= 365*no_yrs), aes(x = time, y = value, colour = name), linewidth = 2) +
  scale_x_continuous(
    breaks = seq(0, max(output$time), by = 365),
    labels = start_year + seq(0, (max(output$time) %/% 365), 1)
  ) +
  labs(
    title = "Mean monthly temperatures",
    x = "Time (days)",
    y = expression("Value ("*degree*"C)"),
    colour = "Area"
  ) +
  theme(legend.position = "bottom") +
  theme_weather_plots()

# rainfall
monthly_rainfall <- ggplot() +
  geom_point(data = filter(rains_values, time >= 365*skip_yrs, time <= 365*no_yrs), aes(x = time, y = value, colour = name), size = 4) +
  geom_line(data = filter(approx_rains_values_long, time >= 365*skip_yrs, time <= 365*no_yrs), aes(x = time, y = value, colour = name), linewidth = 2) +
  scale_x_continuous(
    breaks = seq(0, max(output$time), by = 365),
    labels = start_year + seq(0, (max(output$time) %/% 365), 1)
  ) +
  labs(
    title = "Mean monthly rainfall",
    x = "Time (days)",
    y = "Value (mm)",
    colour = "Area"
  ) +
  theme(legend.position = "bottom") +
  theme_weather_plots()

# mosquito populations
aquatic_mosquitoes <- output %>% 
  filter(time >= 365*skip_yrs, state_var %in% c("E_a", "L_a", "P_a")) %>% 
  ggplot() +
  geom_line(aes(x = time, y = value, colour = as_factor(state_var), linetype = scenario), linewidth = 2) +
  facet_grid(cols = vars(approach), rows = vars(patch), scales = "free") +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_number(suffix = "M", scale = 1e-6)
  ) +
  scale_x_continuous(
    breaks = seq(0, max(output$time), by = 365), 
    labels = start_year + seq(0, (max(output$time) %/% 365), 1)
  ) +
  labs(
    title = "Aquatics",
    colour = "population",
    linetype = "coverage"
  ) +
  theme_weather_plots()

adult_mosquitoes <- output %>% 
  filter(time >= 365*skip_yrs, state_var %in% c("S_m", "E_m", "I_m")) %>% 
  ggplot() +
  geom_line(aes(x = time, y = value, colour = as_factor(state_var), linetype = scenario), linewidth = 2) +
  facet_grid(cols = vars(approach), rows = vars(patch), scales = "free") +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_number(suffix = "M", scale = 1e-6)
  ) +
  scale_x_continuous(
    breaks = seq(0, max(output$time), by = 365), 
    labels = start_year + seq(0, (max(output$time) %/% 365), 1)
  ) +
  labs(
    title = "Adults",
    colour = "population",
    linetype = "coverage"
  ) +
  theme_weather_plots()

# host populations
human_hosts <- output %>% 
  filter(time >= 365*skip_yrs, state_var %in% c("S_h", "E_h", "A_h", "Iu_h", "Is_h", "Tu_h", "Ts_h", "R_h")) %>% 
  ggplot() +
  geom_line(aes(x = time, y = value, colour = as_factor(state_var), linetype = scenario), linewidth = 2) +
  facet_grid(cols = vars(approach), rows = vars(patch), scales = "free") +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_number(suffix = "K", scale = 1e-3)
  ) +
  scale_x_continuous(
    breaks = seq(0, max(output$time), by = 365), 
    labels = start_year + seq(0, (max(output$time) %/% 365), 1)
  ) +
  labs(
    title = "Hosts",
    colour = "population",
    linetype = "coverage"
  ) +
  theme_weather_plots()

daily_incidence <- output %>%
  filter(time >= 365*skip_yrs, state_var == "Inc") %>%
  ggplot() +
  geom_line(aes(x = time, y = value, colour = as_factor(scenario)), linewidth = 2) +
  facet_grid(cols = vars(approach), rows = vars(patch), scales = "free") +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_number(suffix = "", scale = 1e0)
  ) +
  scale_x_continuous(
    breaks = seq(0, max(output$time), by = 365),
    labels = start_year + seq(0, (max(output$time) %/% 365), 1)
  ) +
  labs(
    title = "Incidence",
    colour = "coverage",
    linetype = "coverage",
    y = "cases per 1000" 
  ) +
  theme_weather_plots()

line_annual_incidence <- annual_incidence %>%
  filter(time > 0) %>%
  ggplot() +
  geom_line(aes(x = start_year+time, y = annual_value, colour = as_factor(scenario)), linewidth = 2) +
  facet_grid(cols = vars(approach), rows = vars(name), scales = "free") +
  scale_y_continuous(
    limits = c(0, NA),
    labels = scales::label_number(suffix = "", scale = 1e-3)
  ) +
  labs(
    title = "Annual incidence",
    x = "year",
    y = "cases per 100 000",
    colour = "coverage"
  ) +
  theme_weather_plots()

output_plots <- list(
  "Monthly_Temperature" = monthly_temperature,
  "Monthly_Rainfall" = monthly_rainfall,
  "Aquatic_Mosquitoes" = aquatic_mosquitoes,
  "Adult_Mosquitoes" = adult_mosquitoes,
  "Human_Hosts" = human_hosts,
  "Daily_Incidence" = daily_incidence,
  "Line_Annual_Incidence" = line_annual_incidence
)

plot_e_t <- Sys.time()

total_e_t <- Sys.time()

print(paste(num_cores, "cores used for parallel simulation"))
packs_e_t - packs_s_t
proc_e_t - proc_s_t
sim_e_t - sim_s_t
pproc_e_t - pproc_s_t
plot_e_t - plot_s_t
total_e_t - total_s_t

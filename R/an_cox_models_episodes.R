###############################################################################
#R script author: Kirsty Andresen
#Date: 16 September 2024
#Description: Creates crude + adj + time since diag strat Cox models for each cancer & outcome event (exacerbations) 
############# hosp/death episodes are recurrent events.
###############################################################################
# 
# # # ####################
# # # #TESTING#
# # # # ####################
# # # # outcomes_all <- "asthma"
# outcomes_chronic <- c("asthma")
# outcome <- "asthma"
# cancersites <- c("lun", "col", "thy")
# # 
# # site <- "lun"
# # # ####################

# Define outcome-specific variables

covariates_list <- c("age_up", "gender", "region", "imd5", "bmi_cat", "smokstatus", "b_asthma", "b_copd", "b_bronch", 
                "b_ckd", "b_cvdgrouped", "b_hypertension", "b_diab_cat", "b_cld", "b_multiple_sclerosis")


err <- data.frame(cancer = numeric(), 
                  outcome = numeric(),
                  invalid_pat = numeric())
  
results <- data.frame(
  events_exp = numeric(),
  events_unexp = numeric(),
  cancer = character(),
  outcome = character(),
  hr = numeric(),
  ci_lower = numeric(),
  ci_upper = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for(outcome in outcomes_chronic) {
  
  
  exc_num <- paste0(outcome, "_exc_num")
  exc_contact <- paste0(outcome, "_exc_contact")
  history_var <- paste0("b_", outcome)
  
  
  episode <- readRDS(paste0(path_datafiles_for_analysis, "cr_episodes_", outcome, "_new.rds"))
  
  episode <- episode %>%
    mutate(
      yob = as.Date(paste(yob, "01", "01", sep = "-")),  # Convert yob to July 1st
      age_up = as.numeric(d_excfup - yob) / 365.25,  # Calculate age at start of time at risk
      age_up_cat = as.factor(case_when(
        age_up < 40 ~ "18-39",
        age_up < 50 ~ "40-49",
        age_up < 60 ~ "50-59",
        age_up < 70 ~ "60-69",
        age_up < 80 ~ "70-79",
        TRUE ~ "80+")))
  

  for (site in cancersites) {
    
  print(paste0("cancer type: ", site))
  
  cr_episode <- episode %>%
    filter(cancer == site) %>% 
    group_by(setid, e_patid) %>%
    arrange(setid, e_patid, tstart) %>% 
    ungroup()
  
  factors <- c("age_up_cat", "imd5", "bmi_cat", "smokstatus", "gender", "region", "eth5_hes")
  cr_episode <- cr_episode %>%
    mutate(across(all_of(factors), as.factor))
  
  covariates <- setdiff(covariates_list, history_var)
  
  if(site == "pro" | site == "bre" | site == "ute" | site == "ova" | site == "cer") {
    covariates <- setdiff(covariates, "gender")
  }
  
  print(head(cr_episode))
  
  missing_summary <- colSums(is.na(cr_episode[ , all_of(covariates)]))
  print(missing_summary)
  cr_episode <- cr_episode %>%
    filter(if_all(all_of(covariates), ~ !is.na(.))) # complete case
  
  
  #Count number of outcome events in exposed and unexposed group and save in a variable
  n_exp <- sum(cr_episode$exposed == 1 & cr_episode$episode == 1)
  n_unexp <- sum(cr_episode$exposed == 0 & cr_episode$episode == 1) 
  N_tot <- sum(cr_episode$episode == 1)
  
  #complete case analysis
  cr_episode <- cr_episode %>% filter(!is.na(region))
  
  #Remove invalid rows and document
  invalid_pats <- cr_episode %>% filter(tstart >= time) %>% distinct(setid, e_patid)
  
  err <- rbind(err, data.frame(cancer = site, outcome = outcome, invalid_pat = nrow(invalid_pats)))
  
  cr_episode <- cr_episode %>% anti_join(invalid_pats, by="setid", "e_patid")

  stset <- as.formula(paste0("Surv(tstart, time, episode) ~ exposed + cluster(e_patid)"))
  stset
  
  cox_model_epi <- tryCatch({
    coxph(stset, data = cr_episode)
  }, warning = function(w) {
    cat("Warning in fitting Cox model: ", w$message, "\n")
    return(NULL)  # Return NULL so the script can continue
  }, error = function(e) {
    # Handle errors in case they occur
    cat("Error in fitting Cox model: ", e$message, "\n")
    return(NULL)
  })
  
  
  # Check if the model ran successfully or not
  if (is.null(cox_model_epi)) {
    cat("Cox model did not converge. Continuing with the rest of the script...\n")
    results <- rbind(results, data.frame(events_exp = n_exp,
                                         events_unexp = n_unexp,
                                         model_type = "Crude",
                                         cancer = site,
                                         outcome = outcome,
                                         hr = c(NA),
                                         ci_lower = c(NA),
                                         ci_upper = c(NA)))
  } else {
  
  sum(is.na(cox_model_epi$y))

  print(summary(cox_model_epi))
  
  # Extract the HR, CI, and p-value
  hr <- exp(coef(cox_model_epi))[1]  # Hazard ratio
  ci <- exp(confint(cox_model_epi))[1,] # Confidence intervals
  p_value <- summary(cox_model_epi)$coefficients[, "Pr(>|z|)"]  # p-value
  model_type <- "Crude"
  
  # Append the results to the results data frame
  results <- rbind(results, data.frame(events_exp = n_exp,
                                       events_unexp = n_unexp,                                    
                                       model_type = model_type,
                                       cancer = site,
                                       outcome = outcome,
                                       hr = hr,
                                       ci_lower = ci[1],
                                       ci_upper = ci[2]
                                       ))
  }  
  
################################################################################
#ADJUSTED MODEL
################################################################################

  #check missing in episode and total_time_exc
  missing_outcome <- cr_episode %>% filter(is.na(cr_episode$episode)) %>% pull(e_patid)
  missing_time <- cr_episode %>% filter(is.na(cr_episode$time)) %>% pull(e_patid)
  
stset_adj <- as.formula(
paste("Surv(tstart, time, episode) ~ exposed + cluster(e_patid) +", paste(covariates, collapse = " + ")))
  

print(stset_adj)
  cox_model_epi_adj <- tryCatch({
    coxph(data = cr_episode,
          formula = stset_adj)
  }, warning = function(w) {
    cat("Warning in fitting Cox model: ", w$message, "\n")
    return(NULL)  # Return NULL so the script can continue
  }, error = function(e) {
    # Handle errors in case they occur
    cat("Error in fitting Cox model: ", e$message, "\n")
    return(NULL)
  }
        )

if (is.null(cox_model_epi_adj)) {
  cat("Cox model did not converge. Continuing with the rest of the script...\n")
  results <- rbind(results, data.frame(events_exp = n_exp,
                                       events_unexp = n_unexp,
                                       model_type = "Adjusted",
                                       cancer = site,
                                       outcome = outcome,
                                       hr = c(NA),
                                       ci_lower = c(NA),
                                       ci_upper = c(NA)))
} else {

      adj_summary <- summary(cox_model_epi_adj)
        
      print(summary(cox_model_epi_adj))
      
      # Extract the HR, CI, and p-value for the first row

      # Hazard ratio for the first variable (exposed)
      hr <- exp(adj_summary$coefficients[1, "coef"])
      
      # Confidence intervals
      ci_lower <- adj_summary$conf.int[1, "lower .95"]
      ci_upper <- adj_summary$conf.int[1, "upper .95"]
      
      # p-value
      p_value <- adj_summary$coefficients[1, "Pr(>|z|)"]
      model_type <- "Adjusted"
      
      # Append the results to the results data frame
      results <- rbind(results, data.frame(events_exp = n_exp,
                                           events_unexp = n_unexp,
                                           model_type = model_type,
                                           cancer = site,
                                           outcome = outcome,
                                           hr = hr,
                                           ci_lower = ci_lower,
                                           ci_upper = ci_upper))
}      
###################################################
###COX MODEL INTERACTION - TIME SINCE DIAGNOSIS##
################################################

cr_episode <- cr_episode %>% mutate(time0 = 0)
  
##Split the data by time since diagnosis
split_survival <- survSplit(
  data = cr_episode,
  cut = c(365, 1095, 1825, 3650),  # Cut points (days)
  end = "time",  # follow-up time in days
  event = "episode",
  start = "time0",  # Start of follow-up
  episode = "survival" # Column with event status (1 = event, 0 = censored)
)

split_survival <- split_survival %>%
  mutate(
    age_up = age_up + (time0/365.25),  # Age at start of each interval
    age_up_cat = case_when(
      age_up < 40 ~ "0-39",
      age_up < 50 ~ "40-49",
      age_up < 60 ~ "50-59",
      age_up < 70 ~ "60-69",
      age_up < 80 ~ "70-79",
      TRUE ~ "80+"
    )
  )

##Convert to factor variables
split_survival$survival <- as.factor(split_survival$survival)

model_type <- c("0-1 years", "1-3 years", "3-5 years", "5-10 years", "10+ years")
out <- c(outcome, outcome, outcome, outcome, outcome)
csite <- c(site, site, site, site, site)
events_exp <- c(n_exp, n_exp, n_exp, n_exp, n_exp)
events_unexp <- c(n_unexp, n_unexp, n_unexp, n_unexp, n_unexp)


## : interactions
adjusted_params_split <- as.formula(
  paste("Surv(time0, time, episode) ~ exposed:survival + cluster(e_patid) +" , paste(covariates, collapse = "+"))
)

print(adjusted_params_split)


print(table(split_survival$exposed, split_survival$survival))

cox_model_split <- tryCatch(
  {
    cox_model_split <-coxph(formula = adjusted_params_split, data = split_survival)
    
  },
  
  warning = function(w) {
    cat("Warning in fitting Cox model: ", w$message, "\n")
    return(NULL)  # Return NULL so the script can continue
  },
  
  error = function(e) {
    # Handle errors in case they occur
    cat("Error in fitting Cox model: ", e$message, "\n")
    return(NULL)
  }
)


# Check if the model ran successfully or not
if (is.null(cox_model_split)) {
  cat("Cox model did not converge. Continuing with the rest of the script...\n")
  results <- rbind(results, data.frame(events_exp = c(NA, NA, NA, NA, NA),
                                       events_unexp = c(NA, NA, NA, NA, NA),
                                       model_type = model_type, 
                                       cancer = csite,
                                       outcome = out,
                                       hr = c(NA, NA, NA, NA, NA),
                                       ci_lower = c(NA, NA, NA, NA, NA),
                                       ci_upper = c(NA, NA, NA, NA, NA)))
} else {
  
  print(summary(cox_model_split))
  
  #Calculate linear combinations
  x <-cox_model_split %>% broom::tidy(conf.int = TRUE, exponentiate = TRUE) %>%   filter(grepl(":", term, fixed = TRUE)) %>% dplyr::select (term, estimate, conf.low, conf.high)
  
  hr <- x$estimate
  ci_lower <- x$conf.low
  ci_upper <- x$conf.high
  
  
  results <- rbind(results, data.frame(events_exp = events_exp,
                                       events_unexp = events_unexp,
                                       model_type = model_type,
                                       cancer = csite,
                                       outcome = out,
                                       hr = hr,
                                       ci_lower = ci_lower,
                                       ci_upper = ci_upper))
  
} #else





}
}

rownames(results) <- NULL
print(results)
saveRDS(results, file = paste0(path_datafiles_for_analysis, "results_exc.rds"))

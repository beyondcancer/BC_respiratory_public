# ####################
# #TESTING#
# # #####################
# # outcomes_all <- "asthma"
# outcomes_chronic <- c("asthma")
# # outcome <- "asthma"
# cancersites <- c("nhl")
# # site <- "bre"
# # ###################

library(lmtest)

## Create data frame to store results

# Define outcome-specific variables

covariates_list <- c("age_up_cat", "gender", "region", "imd5", "bmi_cat", "smokstatus", "b_asthma", "b_copd", "b_bronch",
                     "b_ckd", "b_cvdgrouped", "b_hypertension", "b_diab_cat")


err <- data.frame(cancer = numeric(), 
                  outcome = numeric(),
                  invalid_pat = numeric())

## Create data frame to store results
lrt_results <- data.frame(cancer = character(),
                          outcome = character(),
                          lrt_age = numeric(),
                          lrt_gender = numeric(),
                          lrt_smok = numeric(),
                          lrt_imd = numeric())

results <- data.frame(
  events_exp = numeric(),
  events_unexp = numeric(),
  cancer = character(),
  outcome = character(),
  hr = numeric(),
  ci_lower = numeric(),
  ci_upper = numeric(),
  p_value = numeric(),
  cat = character(),
  stringsAsFactors = FALSE
)


  
  for (outcome in outcomes_chronic) {

    print(paste0("Initiating loop for "))
    print(paste0("Outcome: ", outcome))
    
    exc_num <- paste0(outcome, "_exc_num")
    exc_contact <- paste0(outcome, "_exc_contact")
    history_var <- paste0("b_", outcome)
    
    episode <- readRDS(paste0(path_datafiles_for_analysis, "cr_episodes_", outcome, "_new.rds"))    
    
    episode <- episode %>%
      mutate(
        yob = as.Date(paste(yob, "01", "01", sep = "-")),  # Convert yob to July 1st
        age_up = as.numeric(d_excfup - yob) / 365.25,  # Calculate age at d_exc_fup
        age_up_cat = as.factor(case_when(
          age_up < 40 ~ 1,
          age_up < 60 ~ 2,
          TRUE ~ 3)))
    
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
      
      #Count number of outcome events in exposed and unexposed group and save in a variable
      n_exp <- sum(cr_episode$exposed == 1 & cr_episode$episode == 1)
      n_unexp <- sum(cr_episode$exposed == 0 & cr_episode$episode == 1) 
      N_tot <- sum(cr_episode$episode == 1)
      
      missing_summary <- colSums(is.na(cr_episode[ , all_of(covariates)]))
      print(missing_summary)
      cr_episode <- cr_episode %>%
        filter(if_all(all_of(covariates), ~ !is.na(.))) # complete case
      
      
      base_model <- tryCatch(
        coxph(as.formula(paste0("Surv(tstart, time, episode) ~ exposed + cluster(e_patid) + ", paste(covariates, collapse = "+"))), 
              data = cr_episode),
        error = function(e) {
          cat("Error in fitting Cox model:", e$message, "\n")
          return(NULL)
        },
        warning = function(w) {
          cat("Warning in fitting Cox model:", w$message, "\n")
          return(NULL)
        }
      )
      
      
      safe_lrtest <- function(base_model, new_model) {
        if (is.null(base_model) || is.null(new_model)) {
          return(NA)  # Return NA if either model is missing
        } else {
          return(lrtest(base_model, new_model)$`Pr(>Chisq)`[2])  # Extract LRT p-value
        }
      }
      
      
      
      
    ################################################################################
    # COX MODELS - STRATIFIED BY AGE, SEX, AND SMOKING STATUS
    ################################################################################
    
    ###AGE########
    
    age_strat_params <- as.formula(paste0("Surv(tstart, time, episode) ~ exposed:age_up_cat + cluster(e_patid) + ", paste(covariates, collapse = "+")))
    print(age_strat_params)
    
    
    
    #Define variables for results table
    model_type <- c("18-39", "40-59", ">=60")
    out <- c(outcome, outcome, outcome)
    csite <- c(site, site, site)
    events_exp <- c(n_exp, n_exp, n_exp)
    events_unexp <- c(n_unexp, n_unexp, n_unexp)
    
    
    age_stratified_model <- tryCatch(
      {
        coxph(formula = age_strat_params, data = cr_episode)
      },
      warning = function(w) {
        cat("Warning in fitting stratified Cox model: ", w$message, "\n")
        return(NULL)
      },
      error = function(e) {
        cat("Error in fitting stratified Cox model: ", e$message, "\n")
        return(NULL)
      }
    )
    
    if (is.null(age_stratified_model)) {
      results <- rbind(results, data.frame(events_exp = events_exp,
                                           events_unexp = events_unexp,
                                           model_type = model_type,
                                           cancer = csite,
                                           outcome = out,
                                           hr = c(NA, NA, NA),
                                           ci_lower = c(NA, NA, NA),
                                           ci_upper = c(NA, NA, NA),
                                           p_value = c(NA, NA, NA),
                                           lrtest = rep(NA, 3),
                                           cat = rep("Age", 3)))
      results <- unique(results)
      
    } else {
      
      # Extract HR, CI, and p-value
      print(summary(age_stratified_model))
      
      #Calculate linear combinations
      
      #Calculate linear combinations
      x <-age_stratified_model %>% broom::tidy(conf.int = TRUE, exponentiate = TRUE) %>%   filter(grepl(":", term, fixed = TRUE)) %>% dplyr::select (term, estimate, conf.low, conf.high, p.value)
      
      hr <- x$estimate
      ci_lower <- x$conf.low
      ci_upper <- x$conf.high
      p_value <- x$p.value
      lrt_age <- safe_lrtest(base_model, age_stratified_model)
      lrt_age <- rep(lrt_age, nrow(x))
      
      results <- rbind(results, data.frame(events_exp = events_exp,
                                           events_unexp = events_unexp,
                                           model_type = model_type,
                                           cancer = csite,
                                           outcome = out,
                                           hr = hr,
                                           ci_lower = ci_lower,
                                           ci_upper = ci_upper,
                                           p_value = p_value,
                                           lrtest = lrt_age,
                                           cat = rep("Age", 3)))
      results <- unique(results)
      
    } #else age
    
    ###GENDER########
    
    gender_strat_params <- as.formula(paste0("Surv(tstart, time, episode) ~ exposed:gender + cluster(e_patid) +", paste(covariates, collapse = "+")))
    print(gender_strat_params)
    
    #Define variables for results table
    model_type <- c("Male", "Female")
    out <- c(outcome, outcome)
    csite <- c(site, site)
    events_exp <- c(n_exp, n_exp)
    events_unexp <- c(n_unexp, n_unexp)     
    
    gender_stratified_model <- tryCatch(
      {
        coxph(formula = gender_strat_params, data = cr_episode)
      },
      warning = function(w) {
        cat("Warning in fitting stratified Cox model: ", w$message, "\n")
        return(NULL)
      },
      error = function(e) {
        cat("Error in fitting stratified Cox model: ", e$message, "\n")
        return(NULL)
      }
    )
    
    if (is.null(gender_stratified_model)) {
      results <- rbind(results, data.frame(events_exp = events_exp,
                                           events_unexp = events_unexp,
                                           model_type = model_type, 
                                           cancer = csite,
                                           outcome = out,
                                           hr = c(NA, NA),
                                           ci_lower = c(NA, NA),
                                           ci_upper = c(NA, NA),
                                           p_value = c(NA, NA),
                                           lrtest = rep(NA, 2),
                                           cat = rep("Gender", 2)))
      
      results <- unique(results)
      
    } else {
      
      # Extract HR, CI, and p-value
      print(summary(gender_stratified_model))
      
      #Calculate linear combinations
      
      x <-gender_stratified_model %>% broom::tidy(conf.int = TRUE, exponentiate = TRUE) %>% filter(grepl(":", term, fixed = TRUE)) %>% dplyr::select (term, estimate, conf.low, conf.high, p.value)
      
      hr <- x$estimate
      ci_lower <- x$conf.low
      ci_upper <- x$conf.high
      p_value <- x$p.value
      lrt_gender <- safe_lrtest(base_model, gender_stratified_model)
      lrt_gender <- rep(lrt_gender, nrow(x))
      
      results <- rbind(results, data.frame(events_exp = events_exp,
                                           events_unexp = events_unexp,
                                           model_type = model_type,
                                           cancer = csite,
                                           outcome = out,
                                           hr = hr,
                                           ci_lower = ci_lower,
                                           ci_upper = ci_upper,
                                           p_value = p_value, 
                                           lrtest = lrt_gender,
                                           cat = rep("Gender", 2)))
      results <- unique(results)
      
    } #else gender
    
    
    # ###SMOK########
    # 
    # 
    # smok_strat_params <- as.formula(paste0("Surv(tstart, time, episode) ~ exposed + exposed:smokstatus + cluster(e_patid) +", paste(covariates, collapse = "+")))
    # print(smok_strat_params)
    # 
    # model_type <- c("Never Smoker", "Smoker", "Ex-Smoker")
    # out <- c(outcome, outcome, outcome)
    # csite <- c(site, site, site)
    # events_exp <- c(n_exp, n_exp, n_exp)
    # events_unexp <- c(n_unexp, n_unexp, n_unexp)
    # 
    # smok_stratified_model <- tryCatch(
    #   {
    #     coxph(formula = smok_strat_params, data = cr_episode)
    #   },
    #   warning = function(w) {
    #     cat("Warning in fitting stratified Cox model: ", w$message, "\n")
    #     return(NULL)
    #   },
    #   error = function(e) {
    #     cat("Error in fitting stratified Cox model: ", e$message, "\n")
    #     return(NULL)
    #   }
    # )
    # 
    # if (is.null(smok_stratified_model)) {
    #   results <- rbind(results, data.frame(events_exp = events_exp, 
    #                                        events_unexp = events_unexp, 
    #                                        model_type = c("Never Smoker", "Smoker", "Ex-Smoker"), 
    #                                        cancer = csite,
    #                                        outcome = out,
    #                                        hr = c(NA, NA, NA), 
    #                                        ci_lower = c(NA, NA, NA), 
    #                                        ci_upper = c(NA, NA, NA),
    #                                        p_value = c(NA, NA, NA)))
    #   results <- unique(results)
    # } else {
    #   
    #   # Extract HR, CI, and p-value
    #   print(summary(smok_stratified_model))
    #   
    #   #Calculate linear combinations
    #   
    #   x <-smok_stratified_model %>% broom::tidy(conf.int = TRUE, exponentiate = TRUE) %>% filter(grepl(":", term, fixed = TRUE)) %>% dplyr::select (term, estimate, conf.low, conf.high, p.value)
    #   
    #   hr <- x$estimate
    #   ci_lower <- x$conf.low
    #   ci_upper <- x$conf.high
    #   p_value <- x$p.value

    #   
    #   results <- rbind(results, data.frame(events_exp = events_exp,
    #                                        events_unexp = events_unexp,
    #                                        model_type = model_type,
    #                                        cancer = csite,
    #                                        outcome = out,
    #                                        hr = hr,
    #                                        ci_lower = ci_lower,
    #                                        ci_upper = ci_upper,
    #                                        p_value = p_value,
    #                                        lrtest = lrt_smok))
    #   results <- unique(results)
    #   
    #   
    #   
    #   
    #   
    # } #else smok
    
    imd_strat_params <- as.formula(paste0("Surv(tstart, time, episode) ~ exposed:imd5 + cluster(e_patid) +", paste(covariates, collapse = "+")))
    imd_strat_params
    
    model_type <- c("1(Most deprived)", "2", "3", "4", "5(Least deprived)")
    out <- c(outcome, outcome, outcome, outcome, outcome)
    csite <- c(site, site, site, site, site)
    events_exp <- c(n_exp, n_exp, n_exp, n_exp, n_exp)
    events_unexp <- c(n_unexp, n_unexp, n_unexp, n_unexp, n_unexp)
    
    imd_stratified_model <- tryCatch(
      {
        coxph(formula = imd_strat_params, data = cr_episode)
      },
      warning = function(w) {
        cat("Warning in fitting stratified Cox model: ", w$message, "\n")
        return(NULL)
      },
      error = function(e) {
        cat("Error in fitting stratified Cox model: ", e$message, "\n")
        return(NULL)
      }
    )
    
    if (is.null(imd_stratified_model)) {
      results <- rbind(results, data.frame(events_exp = events_exp, 
                                           events_unexp = events_unexp, 
                                           model_type = model_type,
                                           cancer = csite,
                                           outcome = out,
                                           hr = c(NA, NA, NA, NA, NA), 
                                           ci_lower = c(NA, NA, NA, NA, NA), 
                                           ci_upper = c(NA, NA, NA, NA, NA),
                                           p_value = c(NA, NA, NA, NA, NA),
                                           lrtest= rep(NA, 5),
                                           cat = rep("Deprivation, 5")))
      results <- unique(results)
    } else {
      
      # Extract HR, CI, and p-value
      print(summary(imd_stratified_model))
      
      x <-imd_stratified_model %>% broom::tidy(conf.int = TRUE, exponentiate = TRUE) %>% filter(grepl(":", term, fixed = TRUE)) %>% dplyr::select (term, estimate, conf.low, conf.high, p.value)
      
      hr <- x$estimate
      ci_lower <- x$conf.low
      ci_upper <- x$conf.high
      p_value <- x$p.value
      lrt_imd <- safe_lrtest(base_model, imd_stratified_model)
      lrt_imd <- rep(lrt_imd, nrow(x))
      
      results <- rbind(results, data.frame(events_exp = events_exp,
                                           events_unexp = events_unexp,
                                           model_type = model_type,
                                           cancer = csite,
                                           outcome = out,
                                           hr = hr,
                                           ci_lower = ci_lower,
                                           ci_upper = ci_upper,
                                           p_value = p_value,
                                           lrtest = lrt_imd, 
                                           cat = rep("Deprivation, 5")))
      results <- unique(results)
    } #else imd  
   #  
   #  
   # cr_episode_eth <- cr_episode %>% filter(!is.na(eth5_hes)) %>% mutate(eth5_hes = as.factor(eth5_hes))
   # 
   # eth5_strat_params <- as.formula(paste0("Surv(tstart, time, episode) ~ exposed:eth5_hes + cluster(e_patid) +", paste(covariates, collapse = "+")))
   # print(eth5_strat_params)
   # 
   # model_type <- c("White", "Asian", "Black", "Other")
   # out <- c(outcome, outcome, outcome, outcome)
   # csite <- c(site, site, site, site)
   # events_exp <- c(n_exp, n_exp, n_exp, n_exp)
   # events_unexp <- c(n_unexp, n_unexp, n_unexp, n_unexp)
   # 
   # eth5_stratified_model <- tryCatch(
   #   {
   #     coxph(formula = eth5_strat_params, data = cr_episode_eth)
   #   },
   #   warning = function(w) {
   #     cat("Warning in fitting stratified Cox model: ", w$message, "\n")
   #     return(NULL)
   #   },
   #   error = function(e) {
   #     cat("Error in fitting stratified Cox model: ", e$message, "\n")
   #     return(NULL)
   #   }
   # )
   # 
   # if (is.null(eth5_stratified_model)) {
   #   results <- rbind(results, data.frame(events_exp = events_exp, 
   #                                        events_unexp = events_unexp, 
   #                                        model_type = model_type,
   #                                        cancer = csite,
   #                                        outcome = out,
   #                                        hr = c(NA, NA, NA, NA), 
   #                                        ci_lower = c(NA, NA, NA, NA), 
   #                                        ci_upper = c(NA, NA, NA, NA)))
   #   results <- unique(results)
   # } else {
   #   
   #   # Extract HR, CI, and p-value
   #   print(summary(eth5_stratified_model))
   #   
   #   x <-eth5_stratified_model %>% broom::tidy(conf.int = TRUE, exponentiate = TRUE) %>% filter(grepl(":", term, fixed = TRUE)) %>% dplyr::select (term, estimate, conf.low, conf.high, p.value)
   #   
   #   hr <- x$estimate
   #   ci_lower <- x$conf.low
   #   ci_upper <- x$conf.high
   #   
   #   results <- rbind(results, data.frame(events_exp = events_exp,
   #                                        events_unexp = events_unexp,
   #                                        model_type = model_type,
   #                                        cancer = csite,
   #                                        outcome = out,
   #                                        hr = hr,
   #                                        ci_lower = ci_lower,
   #                                        ci_upper = ci_upper))
   #   results <- unique(results)
   # } #else eth5
    
    # lrt_smok <- safe_lrtest(base_model, smok_stratified_model)
    
    # lrt_results <- rbind(lrt_results, data.frame(cancer = site, 
    #                                              outcome = outcome,
    #                                              lrt_age = lrt_age,
    #                                              lrt_gender = lrt_gender, 
    #                                              lrt_smok = lrt_smok,
    #                                              lrt_imd = lrt_imd))
    # 
    # 
  } #outcome loop
} #cancer loop


results <- unique(results)
rownames(results) <- NULL
print(results)

#save summary tables as csv

saveRDS(results, file = paste0(path_datafiles_for_analysis, "an_cox_models_episodes_strat.rds"))
write.csv(lrt_results, file = paste0(path_results, "lrtresults_exc.csv"))


# 

# # ########################
# # ##FOR TESTING PURPOSES##
# cancersites <- cancersites[1:4]
# outcomes_chronic <- c("asthma")
# # #########################

library(lmtest)
## Create data frame to store results
lrt_results <- data.frame(cancer = character(),
                          outcome = character(),
                          lrt_age = numeric(),
                          lrt_gender = numeric(),
                          lrt_smok = numeric(),
                          lrt_imd = numeric())

results <- data.frame(events_exp = numeric(),
                      events_unexp = numeric(),
                      cancer = character(),
                      outcome = character(),
                      hr = numeric(),
                      ci_lower = numeric(),
                      ci_upper = numeric(),
                      p_value = numeric(),
                      lrtest = numeric(),
                      cat = character(),
                      stringsAsFactors = FALSE)

## Initiate cancer type loop
print(paste0("Initiating loop for cancer"))

for(site in cancersites) {
  
  print(paste0("cancer type: ", site))
  
  #Initiate outcome loop  
  
  print(paste0("Initiating loop for "))
  
  for (outcome in outcomes_chronic) {
    
    cr_outcome_an <- read.csv(paste0(path_datafiles_for_analysis, site, "_", outcome, "_cr_outcome_an.csv"))
    
    print(paste0("Outcome: ", outcome))
    
    #Create dynamic variables  
    history_var <- paste0("b_", outcome)  # e.g., b_asthma
    event_var <- paste0("any_", outcome, "_inc_dx")  # e.g., asthma_inc_dx
  
    #table of counts of event_var by exposed
    table(cr_outcome_an$exposed, cr_outcome_an[[event_var]], useNA = "ifany")

    covariates <- c("imd5", "bmi_cat", "smokstatus","b_asthma", "b_copd", "b_bronch", "b_ild",
                    "b_ckd", "b_hypertension", "b_cvdgrouped", "b_diab_cat")
    
    factors <- c("age_cat", "imd5", "bmi_cat", "smokstatus", "eth5_hes")
    cr_outcome_an <- cr_outcome_an %>% mutate_at(vars(factors), as.factor)
    covariates <- setdiff(covariates, c(history_var))
    
    # if(site == "pro" | site == "bre" | site == "ute" | site == "ova" | site == "cer") {
    #   covariates <- setdiff(covariates, "gender")
    # }
    # 
    #Count number of outcome events in exposed and unexposed group and save in a variable 
    n_exp <- sum(cr_outcome_an$exposed == 1 & cr_outcome_an[[event_var]] == 1)
    n_unexp <- sum(cr_outcome_an$exposed == 0 & cr_outcome_an[[event_var]] == 1) 
    
    
    ################################################################################
    # COX MODELS - STRATIFIED BY AGE, SEX, AND SMOKING STATUS
    ################################################################################
   
    cr_outcome_an <- cr_outcome_an %>%
      mutate(age_cat2 = case_when(
        as.numeric(as.character(age_cat)) %in% c(1) ~ 1,
        as.numeric(as.character(age_cat)) %in% c(2) ~ 2,
        as.numeric(as.character(age_cat)) %in% c(3, 4) ~ 3
      )) %>%
      mutate(age_cat2 = factor(age_cat2, levels = c(1, 2, 3), labels = c("18-39", "40-59", ">=60")))
    
       
     cr_outcome_an <- cr_outcome_an %>%
      mutate(age_cat = as.factor(age_cat), 
             age_cat2 = as.factor(age_cat2),
             gender = as.factor(gender),
             smokstatus = as.factor(smokstatus))
    
     base_model <- tryCatch(
       coxph(as.formula(paste("Surv(total_fup_var,", event_var, ") ~ exposed +strata(setid) +", paste(covariates, collapse = "+"))), 
             data = cr_outcome_an),
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
     
     
    
    ###AGE########
    
    age_strat_params <- as.formula(paste("Surv(total_fup_var,", event_var, ") ~ exposed:age_cat2 + strata(setid) + ", paste(covariates, collapse = "+")))
    age_strat_params
    
    #Define variables for results table
    model_type <- c( "18-39", "40-59", ">=60")
    out <- c(outcome, outcome, outcome)
    csite <- c(site, site, site)
    events_exp <- c(n_exp, n_exp, n_exp)
    events_unexp <- c(n_unexp, n_unexp, n_unexp)
    
    
    age_stratified_model <- tryCatch(
      {
        coxph(formula = age_strat_params, data = cr_outcome_an)
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
                                           lrtest = c(NA, NA, NA),
                                           cat = rep("Age", 3)))
      results <- unique(results)
      
    } else {
      
      # Extract HR, CI, and p-value
      
      print(summary(age_stratified_model))
      
      #Calculate linear combinations
      x <-age_stratified_model %>% broom::tidy(conf.int = TRUE, exponentiate = TRUE) %>%   filter(grepl(":", term, fixed = TRUE)) %>% dplyr::select (term, estimate, conf.low, conf.high,  p.value)
      
      hr <- x$estimate
      ci_lower <- x$conf.low
      ci_upper <- x$conf.high
      p_value <- x$p.value
      lrt_age <- safe_lrtest(base_model, age_stratified_model)
      #repeat for the number of rows in hr
      lrt_test <- rep(lrt_age, nrow(x))

      results <- rbind(results, data.frame(events_exp = events_exp,
                                           events_unexp = events_unexp,
                                           model_type = model_type,
                                           cancer = csite,
                                           outcome = out,
                                           hr = hr,
                                           ci_lower = ci_lower,
                                           ci_upper = ci_upper,
                                           p_value = p_value,
                                           lrtest = lrt_test,
                                           cat = rep("Age", nrow(x))))
      results <- unique(results)
    
      } #else age
 

    ###GENDER########
    

    
    gender_strat_params <- as.formula(paste("Surv(total_fup_var,", event_var, ") ~ exposed:gender  + strata(setid) +", paste(covariates, collapse = "+")))
    gender_strat_params
    
    #Define variables for results table
    model_type <- c("Male", "Female")
    out <- c(outcome, outcome)
    csite <- c(site, site)
    events_exp <- c(n_exp, n_exp)
    events_unexp <- c(n_unexp, n_unexp)     
    
    gender_stratified_model <- tryCatch(
      {
        coxph(formula = gender_strat_params, data = cr_outcome_an)
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
                                           lrtest = c(NA, NA),
                                           cat = rep("Gender", 2)))
      results <- unique(results)
    
      } else {

      # Extract HR, CI, and p-value
      print(summary(gender_stratified_model))
      
      #Calculate linear combinations
        
        #Calculate linear combinations
        x <-gender_stratified_model %>% broom::tidy(conf.int = TRUE, exponentiate = TRUE) %>% filter(grepl(":", term, fixed = TRUE)) %>% dplyr::select (term, estimate, conf.low, conf.high,  p.value)
        
        hr <- x$estimate
        ci_lower <- x$conf.low
        ci_upper <- x$conf.high
        p_value <- x$p.value
        lrt_gender <- safe_lrtest(base_model, gender_stratified_model)
        lrt_test <- rep(lrt_gender, nrow(x))
      
      results <- rbind(results, data.frame(events_exp = events_exp,
                                           events_unexp = events_unexp,
                                           model_type = model_type,
                                           cancer = csite,
                                           outcome = out,
                                           hr = hr,
                                           ci_lower = ci_lower,
                                           ci_upper = ci_upper, 
                                           p_value = p_value, 
                                           lrtest = lrt_test,
                                          cat = rep("Gender", nrow(x))))
      results <- unique(results)
      
    } #else gender
    
    
    # ###SMOK########
    # 
    # 
    # smok_strat_params <- as.formula(paste("Surv(total_fup_var,", event_var, ") ~ exposed:smokstatus + ", paste(covariates, collapse = "+")))
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
    #     coxph(formula = smok_strat_params, data = cr_outcome_an)
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
    #   x <-smok_stratified_model %>% broom::tidy(conf.int = TRUE, exponentiate = TRUE) %>% filter(grepl(":", term, fixed = TRUE)) %>% dplyr::select (term, estimate, conf.low, conf.high,  p.value)
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
    #                                        p_value = p_value))
    #   results <- unique(results)
    # } #else smok
    
    
###########DEPRIVATION
    
    imd_strat_params <- as.formula(paste("Surv(total_fup_var,", event_var, ") ~ exposed:imd5 + strata(setid) +", paste(covariates, collapse = "+")))
    imd_strat_params
    
    model_type <- c("1(Most deprived)", "2", "3", "4", "5(Least deprived)")
    out <- c(outcome, outcome, outcome, outcome, outcome)
    csite <- c(site, site, site, site, site)
    events_exp <- c(n_exp, n_exp, n_exp, n_exp, n_exp)
    events_unexp <- c(n_unexp, n_unexp, n_unexp, n_unexp, n_unexp)
    
    imd_stratified_model <- tryCatch(
      {
        coxph(formula = imd_strat_params, data = cr_outcome_an)
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
                                           lrtest = c(NA, NA, NA, NA, NA),
                                           cat = rep("Deprivation", 5)))
      results <- unique(results)
    } else {
      
      # Extract HR, CI, and p-value
      print(summary(imd_stratified_model))
      
      x <-imd_stratified_model %>% broom::tidy(conf.int = TRUE, exponentiate = TRUE) %>% filter(grepl(":", term, fixed = TRUE)) %>% dplyr::select (term, estimate, conf.low, conf.high,  p.value)
      
      hr <- x$estimate
      ci_lower <- x$conf.low
      ci_upper <- x$conf.high
      p_value <- x$p.value
      lrt_imd <- safe_lrtest(base_model, imd_stratified_model)
      lrt_test <- rep(lrt_imd, nrow(x))
      results <- rbind(results, data.frame(events_exp = events_exp,
                                           events_unexp = events_unexp,
                                           model_type = model_type,
                                           cancer = csite,
                                           outcome = out,
                                           hr = hr,
                                           ci_lower = ci_lower,
                                           ci_upper = ci_upper,
                                           p_value = p_value, 
                                           lrtest = lrt_test,
                                           cat = rep("Deprivation", nrow(x))))
      results <- unique(results)
    } #else imd  
    
    # lrt_smok <- safe_lrtest(base_model, smok_stratified_model)
    # 
    # lrt_results <- rbind(lrt_results, data.frame(cancer = site, 
    #                           outcome = outcome,
    #                           lrt_age = lrt_age,
    #                           lrt_gender = lrt_gender, 
    #                           lrt_smok = lrt_smok,
    #                           lrt_imd = lrt_imd))
    # 
    # 
    # ########### ETHNICITY
    # 
    # cr_outcome_an_eth <- cr_outcome_an %>% filter(!is.na(eth5_hes)) %>% mutate(eth5_hes = as.factor(eth5_hes))
    # 
    # eth5_hes_strat_params <- as.formula(paste("Surv(total_fup_var,", event_var, ") ~ exposed:eth5_hes + eth5_hes + ", paste(covariates, collapse = "+")))
    # print(eth5_hes_strat_params)
    # 
    # model_type <- c("White", "South Asian", "Black", "Other")
    # out <- c(outcome, outcome, outcome, outcome)
    # csite <- c(site, site, site, site)
    # events_exp <- c(n_exp, n_exp, n_exp, n_exp)
    # events_unexp <- c(n_unexp, n_unexp, n_unexp, n_unexp)
    # 
    # eth_stratified_model <- tryCatch(
    #   {
    #     coxph(formula = eth5_hes_strat_params, data = cr_outcome_an_eth)
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
    # if (is.null(eth_stratified_model)) {
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
    #   print(summary(eth_stratified_model))
    #   
    #   x <-eth_stratified_model %>% broom::tidy(conf.int = TRUE, exponentiate = TRUE) %>% filter(grepl(":", term, fixed = TRUE)) %>% dplyr::select (term, estimate, conf.low, conf.high,  p.value)
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
    # } #else eth5_hes  
    
  } #outcome loop
} #cancer loop


results <- unique(results)
rownames(results) <- NULL
print(results)

#save summary tables as csv

saveRDS(results, file = paste0(path_datafiles_for_analysis, "an_cox_models_first_strat.rds"))
# write.csv(lrt_results, file = paste0(path_results, "/lrtest_results.csv"))


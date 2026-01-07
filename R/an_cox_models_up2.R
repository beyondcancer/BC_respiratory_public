###############################################################################
#R script author: Kirsty Andresen
#Date: 16 September 2024
#Description: Creates crude + adj Cox models for each cancer & outcome event (chronic diseases)
###############################################################################
# # ########################
# # ##FOR TESTING PURPOSES##
# cancersites <- top9cancers
# # #########################


# Load data


## Create data frame to store results

results <- data.frame(events_exp = numeric(),
                      events_unexp = numeric(),
                      cancer = character(),
                      outcome = character(),
                      hr = numeric(),
                      ci_lower = numeric(),
                      ci_upper = numeric(),
                      p_value = numeric(),
                      stringsAsFactors = FALSE)


## Initiate cancer type loop
print(paste0("Initiating loop for cancer"))

for(site in cancersites) {
  
  print(paste0("cancer type: ", site))
  
  #Initiate outcome loop  
  
  print(paste0("Initiating loop for "))
  
  for (outcome in outcomes_chronic) {
    
    print(paste0("Outcome: ", outcome))
    
    cr_outcome_an <- read.csv(paste0(path_datafiles_for_analysis, site, "_", outcome, "_cr_outcome_an.csv"))
    
    #Create dynamic variables  
    history_var <- paste0("b_", outcome)  # e.g., b_asthma
    event_var <- paste0("any_", outcome, "_inc_dx")  # e.g., asthma_inc_dx
    covariates <- c("imd5", "bmi", "smokstatus","b_asthma", "b_copd", "b_bronch", "b_ild", "b_ckd", "b_hypertension", "b_cvdgrouped", "b_diab_cat", "b_cld", "b_multiple_sclerosis")
    
    #Remove history of var
    covariates <- setdiff(covariates, c(history_var))
    
    #
    # #Complete case anlysis
    # missing_summary <- colSums(is.na(cr_outcome_an[ , all_of(covariates)]))
    # print(paste0("Number of missing values: "))
    # print(missing_summary)
    #
    # cr_outcome_an <- cr_outcome_an %>%
    #   filter(if_all(all_of(covariates), ~ !is.na(.))) # complete case
    # print(paste0("Number of complete cases ", count(cr_outcome_an)))
    #

    #Count number of outcome events in exposed and unexposed group and save in a variable
    n_exp <- sum(cr_outcome_an$exposed == 1 & cr_outcome_an[[event_var]] == 1)
    n_unexp <- sum(cr_outcome_an$exposed == 0 & cr_outcome_an[[event_var]] == 1)

    # ################################################################################
    # # PROPORTIONAL HAZARDS TEST
    # ################################################################################
    # print("Proportional Hazards Test")
    # # Test proportional hazards assumption
    # cox_model_test <- coxph(Surv(total_fup_var, cr_outcome_an[[event_var]]) ~ exposed, data = cr_outcome_an)
    #
    # ph_test <- cox.zph(cox_model_test)
    #
    # # Print the test results
    # print(ph_test)
    #

    ################################################################################
    # COX MODEL - UNADJUSTED
    ################################################################################

    cox_model <- try(coxph(Surv(total_fup_var, cr_outcome_an[[event_var]]) ~ exposed + strata(setid), data = cr_outcome_an))

    if (inherits(cox_model, "try-error")) {
      cat("Cox model did not converge. Continuing with the rest of the script...\n")

    } else {

      print(summary(cox_model))

      # Extract the HR, CI, and p-value
      hr <- exp(coef(cox_model))  # Hazard ratio
      ci <- exp(confint(cox_model))  # Confidence intervals
      p_value <- summary(cox_model)$coefficients[, "Pr(>|z|)"]  # p-value
      model_type <- "Crude"

      # Append the results to the results data frame
      results <- rbind(results, data.frame(events_exp = n_exp,
                                           events_unexp = n_unexp,
                                           model_type = model_type,
                                           cancer = site,
                                           outcome = outcome,
                                           hr = hr,
                                           ci_lower = ci[1],
                                           ci_upper = ci[2]))

    }
    ################################################################################
    # COX MODEL - ADJUSTED
    ################################################################################

    adjusted_params <- as.formula(
      paste("Surv(total_fup_var,", event_var, ") ~ exposed + strata(setid) + ", paste(covariates, collapse = "+"))
    )

    cox_model_adj <- tryCatch(

      {
        coxph(data = cr_outcome_an, formula = adjusted_params)
      },

      warning = function(w) {
        cat("Warning in fitting Cox model : ", w$message, "\n")
        return(NULL)  # Return NULL so the script can continue
      },
      error = function(e) {
        # Handle errors in case they occur
        cat("Error in fitting Cox model: ", e$message, "\n")
        return(NULL)
      }

    )

    # Check if the model ran successfully or not
    if (is.null(cox_model_adj)) {
      cat("Cox model did not converge. Continuing with the rest of the script...\n")
      results <- rbind(results, data.frame(events_exp = n_exp,
                                           events_unexp = n_unexp,
                                           model_type = "Adjusted",
                                           cancer = site,
                                           outcome = outcome,
                                           hr = NA,
                                           ci_lower = NA,
                                           ci_upper = NA))
    } else {

      print(summary(cox_model_adj))


      # Extract the HR, CI, and p-value
      hr <- exp(coef(cox_model_adj)[1])  # Hazard ratio
      ci <- exp(confint(cox_model_adj)[1,])  # Confidence intervals
      p_value <- summary(cox_model_adj)$coefficients[1, "Pr(>|z|)"]  # p-value
      model_type <- "Adjusted"



      # Append the results to the results data frame
      results <- rbind(results, data.frame(events_exp = n_exp,
                                           events_unexp = n_unexp,
                                           model_type = model_type,
                                           cancer = site,
                                           outcome = outcome,
                                           hr = hr,
                                           ci_lower = ci[1],
                                           ci_upper = ci[2]))

    }

      ###COX MODEL INTERACTION - Time UP TO DIAGNOSIS###

      cutoffs <- c(1, 3, 5, 10)
      hr_up2 <- list()
      ci_lower_up2 <- list()
      ci_upper_up2 <- list()
      model_typeup2 <- list()
      siteup2 <- list()
      outcomeup2 <- list()


      for (cutoff in cutoffs) {
        

        name <- as.character(cutoff)
        model_name <- paste("up to", cutoff, "years")

        print(model_name)
        
        siteup2[[name]] <- site
        outcomeup2[[name]] <- outcome

        cr_outcome_up2 <- cr_outcome_an %>%
        mutate(up2_fup = pmin(total_fup_var, cutoff),
               up2_event_var = ifelse(total_fup_var <= cutoff, !!sym(event_var), 0))

        up2_params <- as.formula(paste("Surv(up2_fup, up2_event_var) ~ exposed+ strata(setid) +", paste(covariates, collapse = "+")))


        cox_model_up2 <- try(coxph(formula = up2_params, data = cr_outcome_up2))

        if (inherits(cox_model_up2, "try-error")) {
        print(cat("Cox model did not converge for cutoff", cutoff, "\n"))

        hr_up2<- NA
        ci_lower_up2[[name]] <- NA
        ci_upper_up2[[name]] <- NA
        model_typeup2[[name]] <- model_name

        } else {

      print(summary(cox_model_up2))

      hr_up2[[name]]<- exp(coef(cox_model_up2)[1])
      ci <- exp(confint(cox_model_up2)[1,])  # Confidence intervals
      ci_lower_up2[[name]] <- ci[1]
      ci_upper_up2[[name]] <- ci[2]
      model_typeup2[[name]] <- model_name
      }
      }

      resultsup2 <- data.frame(events_exp = n_exp,
                               events_unexp = n_unexp,
                               model_type = unlist(model_typeup2),
                               cancer = unlist(siteup2),
                               outcome = unlist(outcomeup2),
                               hr = unlist(hr_up2),
                               ci_lower = unlist(ci_lower_up2),
                               ci_upper = unlist(ci_upper_up2))

      results <- rbind(results, resultsup2)

      # ###COX MODEL INTERACTION - TIME SINCE DIAGNOSIS###
      # ################################################
      # 
      # #create a new variable to split the data by time since diagnosis
      # cr_outcome_an <- cr_outcome_an %>%
      #   mutate(total_fup_check = total_fup_var)
      # 
      # ##Split the data by time since diagnosis
      # split_survival <- survSplit(
      #   data = cr_outcome_an,
      #   cut = c(1, 3, 5),  # Cut points (years)
      #   end = "total_fup_var",  # follow-up time in years
      #   event = event_var,
      #   start = "start_time",  # Start of follow-up
      #   episode = "survival" # Column with event status (1 = event, 0 = censored)
      # )
      # 
      # ##Convert to factor variables
      # split_survival$survival <- as.factor(split_survival$survival)
      # split_survival$exposed <- as.factor(split_survival$exposed)
      # 
      # 
      # ## : interactions
      # adjusted_params_split <- as.formula(
      #   paste("Surv(start_time, total_fup_var,", event_var, ") ~ exposed * survival + strata(setid) +", paste(covariates, collapse = "+"))
      # )
      # print(adjusted_params_split)
      # 
      # 
      # print(table(split_survival$exposed, split_survival$survival))
      # print(table(split_survival[[event_var]], split_survival$survival))
      # 
      # cox_model_split <- tryCatch(
      #   {
      #     coxph(formula = adjusted_params_split,  data = split_survival, id=e_patid)
      #   },
      # 
      #   warning = function(w) {
      #     cat("Warning in fitting Cox model: ", w$message, "\n")
      #     return(NULL)  # Return NULL so the script can continue
      #   },
      # 
      #   error = function(e) {
      #     # Handle errors in case they occur
      #     cat("Error in fitting Cox model: ", e$message, "\n")
      #     return(NULL)
      #   }
      # )
      # 
      # 
      # # Check if the model ran successfully or not
      # if (is.null(cox_model_split)) {
      #   cat("Cox model did not converge. Continuing with the rest of the script...\n")
      #   results <- rbind(results, data.frame(events_exp = c(n_exp, n_exp, n_exp, n_exp),
      #                                        events_unexp = c(n_unexp, n_unexp, n_unexp, n_unexp),
      #                                        model_type = c("0-1 years", "1-3 years", "3-5 years", "5+ years"),
      #                                        cancer = c(site, site, site, site),
      #                                        outcome = c(outcome, outcome, outcome, outcome),
      #                                        hr = c(NA, NA, NA, NA),
      #                                        ci_lower = c(NA, NA, NA, NA),
      #                                        ci_upper = c(NA, NA, NA, NA)))
      # } else {
      # 
      #   print(summary(cox_model_split))
      # 
      #   #Calculate linear combinations
      # 
      #   # Extract coefficients
      #   hr <- c(
      #     coef(cox_model_split)["exposed1"],  # Main effect for exposed
      #     coef(cox_model_split)["exposed1"] + coef(cox_model_split)["exposed1:survival2"],  # Interval 1-3 years
      #     coef(cox_model_split)["exposed1"] + coef(cox_model_split)["exposed1:survival3"],  # Interval 3-5 years
      #     coef(cox_model_split)["exposed1"] + coef(cox_model_split)["exposed1:survival4"]   # Interval 5+ years
      #   )
      # 
      #   # Convert to hazard ratios
      #   hr <- exp(hr)
      # 
      #   # Confidence intervals
      #   ci <- confint(cox_model_split)
      #   ci_lower <- exp(c(
      #     ci["exposed1", 1],
      #     ci["exposed1", 1] + ci["exposed1:survival2", 1],
      #     ci["exposed1", 1] + ci["exposed1:survival3", 1],
      #     ci["exposed1", 1] + ci["exposed1:survival4", 1]
      #   ))
      #   ci_upper <- exp(c(
      #     ci["exposed1", 2],
      #     ci["exposed1", 2] + ci["exposed1:survival2", 2],
      #     ci["exposed1", 2] + ci["exposed1:survival3", 2],
      #     ci["exposed1", 2] + ci["exposed1:survival4", 2]
      #   ))
      # 
      # 
      #   model_type <- c("0-1 years", "1-3 years", "3-5 years", "5+ years")
      #   out <- c(outcome, outcome, outcome, outcome)
      #   csite <- c(site, site, site, site)
      #   events_exp <- c(n_exp, n_exp, n_exp, n_exp)
      #   events_unexp <- c(n_unexp, n_unexp, n_unexp, n_unexp)
      # 
      # 
      #   results <- rbind(results, data.frame(events_exp = events_exp,
      #                                        events_unexp = events_unexp,
      #                                        model_type = model_type,
      #                                        cancer = csite,
      #                                        outcome = out,
      #                                        hr = hr,
      #                                        ci_lower = ci_lower,
      #                                        ci_upper = ci_upper))
      # 
      # } #else
      # 
      # 
      # 
      # ################################################################################
      

  } #outcome
} #cancer

results <- unique(results)


rownames(results) <- NULL
print(results)

#save summary tables as rds

saveRDS(results, file = paste0(path_datafiles_for_analysis, "sens3/", "an_cox_models_first_dx_crude_adj_up2.rds"))

# remove all except cr_finaldataforanalysis_respiratory
#rm(list=setdiff(ls(), "cr_finaldataforanalysis_respiratory"))
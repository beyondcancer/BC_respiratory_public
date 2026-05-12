###############################################################################
#R script author: Kirsty Andresen
#Date: 06 May 2025
# Description: Creates Cox proportional hazards models (adjusted, chemotherapy,
#              and radiotherapy) for each cancer type and chronic disease outcome.
#              Models compare cancer patients with/without treatment to matched
#              cancer-free controls. Likelihood ratio tests compare treatment
#              models to the base adjusted model. Results saved as RDS file.
##############################################################################################################################################################

# # # ########################
# # # ##FOR TESTING PURPOSES##
cancersites <- c("nhl")
# outcomes_chronic <- c("copd")
# # # #########################
# #
# site <- c("bre")
# outcome <- c("copd")

# Load data


## Create data frame to store results

results <- data.frame(
  events_exp = numeric(),
  events_unexp = numeric(),
  cancer = character(),
  outcome = character(),
  hr = numeric(),
  ci_lower = numeric(),
  ci_upper = numeric(),
  p_value = numeric(),
  p_value = numeric(),
  lrtest = numeric(),
  stringsAsFactors = FALSE
)



## Initiate cancer type loop
print(paste0("Initiating loop for cancer"))

for (site in cancersites) {
  print(paste0("cancer type: ", site))
  
  #Initiate outcome loop
  
  print(paste0("Initiating loop for "))
  
  for (outcome in outcomes_chronic) {
    print(paste0("Outcome: ", outcome))
    
    cr_outcome_an <- read.csv(
      paste0(
        path_datafiles_for_analysis,
        "trt/",
        site,
        "_",
        outcome,
        "_cr_outcome_an.csv"
      )
    )
    
    #Create dynamic variables
    history_var <- paste0("b_", outcome)  # e.g., b_asthma
    event_var <- paste0("any_", outcome, "_inc_dx")  # e.g., asthma_inc_dx
    
    covariates <- c(
      "imd5",
      "bmi",
      "smokstatus",
      "b_asthma",
      "b_copd",
      "b_bronch",
      "b_ild",
      "b_ckd",
      "b_hypertension",
      "b_cvdgrouped",
      "b_diab_cat",
      "b_cld",
      "b_multiple_sclerosis"
    )
    
    covariates <- setdiff(covariates, c(history_var))
    
    # --- Convert covariates to factors ---
    factors <- c(
      "exposed_radio",
      "exposed_chemo",
      "age_cat",
      "gender",
      "imd5",
      "bmi_cat",
      "eth5_hes",
      "eth5_cprd",
      "smokstatus"
    )
    
    cr_outcome_an <- cr_outcome_an %>%
      mutate(across(all_of(factors), as.factor))
    
    missing_summary <- colSums(is.na(cr_outcome_an[, all_of(covariates)]))
    
    print(paste0("Number of missing values: "))
    print(missing_summary)
    
    cr_outcome_an <- cr_outcome_an %>%
      filter(if_all(all_of(covariates), ~ !is.na(.))) # complete case
    
    cr_outcome_an <- cr_outcome_an %>%
      dplyr::select(
        setid,
        e_patid,
        exposed,
        exposed_chemo,
        exposed_radio,
        exposed_hsct,
        excl_radio,
        excl_chemo,
        excl_hsct,
        new_index,
        indexdate,
        doexit_var,
        total_fup_var,
        event_var,
        all_of(covariates),
        age_cat,
        gender,
        eth5_hes
      )
    
    
    print(paste0("Number of complete cases ", count(cr_outcome_an)))
    
    
    #
    # ################################################################################
    # # PROPORTIONAL HAZARDS TEST
    # ################################################################################
    # print("Proportional Hazards Test")
    # # Test proportional hazards assumption
    # cox_model_test <- coxph(Surv(total_fup_var, cr_outcome_an[[event_var]]) ~ exposed_chemo, data = cr_outcome_an)
    #
    # ph_test <- cox.zph(cox_model_test)
    #
    # # Print the test results
    # print(ph_test)
    #
    ################################################################################
    # COX MODEL - ADJUSTED
    ################################################################################
    
    
    
    #
    # cox_model_adj <- tryCatch(
    #
    #   {
    #     coxph(data = cr_outcome_an, formula = adjusted_params)
    #   },
    #
    #   warning = function(w) {
    #     cat("Warning in fitting Cox model : ", w$message, "\n")
    #     return(NULL)  # Return NULL so the script can continue
    #   },
    #   error = function(e) {
    #     # Handle errors in case they occur
    #     cat("Error in fitting Cox model: ", e$message, "\n")
    #     return(NULL)
    #   }
    #
    # )
    #
    # # Check if the model ran successfully or not
    # if (is.null(cox_model_adj)) {
    #   cat("Cox model did not converge. Continuing with the rest of the script...\n")
    #   results <- rbind(results, data.frame(model_type = "Adjusted",
    #                                        cancer = site,
    #                                        outcome = outcome,
    #                                        hr = NA,
    #                                        ci_lower = NA,
    #                                        ci_upper = NA,
    #                                        p_value = NA,
    #                                        lrtest = NA
    #
    # )
    # )
    #
    # } else {
    #
    #   print(summary(cox_model_adj))
    #
    #
    #   # Extract the HR, CI, and p-value
    #   hr <- exp(coef(cox_model_adj)[1])  # Hazard ratio
    #   ci <- exp(confint(cox_model_adj)[1,])  # Confidence intervals
    #   p_value <- summary(cox_model_adj)$coefficients[1, "Pr(>|z|)"]  # p-value
    #   model_type <- "Adjusted"
    #
    #
    #
    #   # Append the results to the results data frame
    #   results <- rbind(results, data.frame(model_type = model_type,
    #                                        cancer = site,
    #                                        outcome = outcome,
    #                                        hr = hr,
    #                                        ci_lower = ci[1],
    #                                        ci_upper = ci[2],
    #                                        p_value = p_value,
    #                                        lrtest = NA
    #   )
    #   )
    #
    # }
    #
    
    adjusted_params <- as.formula(paste(
      "Surv(total_fup_var,",
      event_var,
      ") ~ exposed + strata(setid) + ",
      paste(covariates, collapse = "+")
    ))
    # ################################################################################
    # # COX MODEL - CHEMO
    # ################################################################################
    # 
    # #Drop prior to 2014
    # 
    # cr_outcome_chemo <- cr_outcome_an %>% mutate(
    #   excl_chemo = if_else(indexdate < as.Date("2014-01-01"), 1, 0),
    # ) %>% filter(excl_chemo == 0)
    # print(num_excluded_chemo <- nrow(cr_outcome_an) - nrow(cr_outcome_chemo))
    # 
    # #adjusted for lrt_test
    # 
    # cox_model_adj_chemo <- tryCatch({
    #   coxph(data = cr_outcome_chemo, formula = adjusted_params)
    # }, warning = function(w) {
    #   cat("Warning in fitting Cox model : ", w$message, "\n")
    #   return(NULL)  # Return NULL so the script can continue
    # }, error = function(e) {
    #   # Handle errors in case they occur
    #   cat("Error in fitting Cox model: ", e$message, "\n")
    #   return(NULL)
    # })
    # print(summary(cox_model_adj_chemo))
    # 
    # 
    # 
    # 
    # chemo_params <- as.formula(
    #   paste(
    #     "Surv(total_fup_var,",
    #     event_var,
    #     ") ~ exposed_chemo + strata(setid) + ",
    #     paste(covariates, collapse = "+")
    #   )
    # )
    # 
    # cox_model_chemo <- tryCatch({
    #   coxph(data = cr_outcome_chemo, formula = chemo_params)
    # }, warning = function(w) {
    #   cat("Warning in fitting Cox model : ", w$message, "\n")
    #   return(NULL)  # Return NULL so the script can continue
    # }, error = function(e) {
    #   # Handle errors in case they occur
    #   cat("Error in fitting Cox model: ", e$message, "\n")
    #   return(NULL)
    # })
    # 
    # # Check if the model ran successfully or not
    # if (is.null(cox_model_chemo)) {
    #   cat("Cox model did not converge. Continuing with the rest of the script...\n")
    #   results <- rbind(
    #     results,
    #     data.frame(
    #       model_type = c("Chemo No Treatment", "Chemo Treatment"),
    #       cancer = c(site, site),
    #       outcome = c(outcome, outcome),
    #       hr = c(NA, NA),
    #       ci_lower = c(NA, NA),
    #       ci_upper = c(NA, NA),
    #       p_value = c(NA, NA),
    #       lrtest = c(NA, NA)
    #     )
    #   )
    #   
    # } else {
    #   print(summary(cox_model_chemo))
    #   
    #   
    #   # For exposed_chemo = 1 (No Treatment) vs 0 (Cancer-free - reference)
    #   hr_no_treatment <- exp(coef(cox_model_chemo)["exposed_chemoNo Treatment"])
    #   ci_no_treatment <- exp(confint(cox_model_chemo)["exposed_chemoNo Treatment", ])
    #   p_no_treatment <- summary(cox_model_chemo)$coefficients["exposed_chemoNo Treatment", "Pr(>|z|)"]
    #   
    #   # Or by row number if names don't work
    #   p_no_treatment <- summary(cox_model_chemo)$coefficients[1, 5]
    #   p_treatment <- summary(cox_model_chemo)$coefficients[2, 5]
    #   # For exposed_chemo = 2 (Treatment) vs 0 (Cancer-free - reference)
    #   hr_treatment <- exp(coef(cox_model_chemo)["exposed_chemoTreatment"])
    #   ci_treatment <- exp(confint(cox_model_chemo)["exposed_chemoTreatment", ])
    #   p_treatment <- summary(cox_model_chemo)$coefficients["exposed_chemoTreatment", "Pr(>|z|)"]
    #   
    #   ##  Likelihood ratio test between cox_model_adj and cox_model_chemo
    #   # CALCULATE LRT HERE - with error handling
    #   lr_test_chemo <- tryCatch({
    #     anova(cox_model_adj_chemo, cox_model_chemo, test = "LRT")
    #   }, warning = function(w) {
    #     cat("Warning in LRT for chemo: ", w$message, "\n")
    #     return(NULL)
    #   }, error = function(e) {
    #     cat("Error in LRT for chemo: ", e$message, "\n")
    #     return(NULL)
    #   })
    #   
    #   # Extract p-value safely
    #   if (!is.null(lr_test_chemo)) {
    #     print("Likelihood Ratio Test - Chemo vs Adjusted Model")
    #     print(lr_test_chemo)
    #     lrt_p_chemo <- lr_test_chemo$`Pr(>|Chi|)`[2]
    #   } else {
    #     cat("LRT could not be calculated for chemo model\n")
    #     lrt_p_chemo <- NA
    #   }
    #   
    #   
    #   results <- rbind(
    #     results,
    #     data.frame(
    #       model_type = c("Chemo No Treatment", "Chemo Treatment"),
    #       cancer = site,
    #       outcome = outcome,
    #       hr = c(hr_no_treatment, hr_treatment),
    #       ci_lower = c(ci_no_treatment[1], ci_treatment[1]),
    #       ci_upper = c(ci_no_treatment[2], ci_treatment[2]),
    #       p_value = c(p_no_treatment, p_treatment),
    #       lrtest = c(
    #         lr_test_chemo$`Pr(>|Chi|)`[2],
    #         lr_test_chemo$`Pr(>|Chi|)`[2]
    #       )
    #       
    #     )
    #   )
    # }
    # 
    # ################################################################################
    # # COX MODEL - RADIO
    # ################################################################################
    # #Drop prior to 2012
    # 
    # cr_outcome_radio <- cr_outcome_an %>% 
    #   mutate(
    #     excl_radio = if_else(indexdate < as.Date("2012-01-01"), 1, 0)
    #   ) %>%
    #   filter(excl_radio == 0)
    # 
    # #print number excluded:
    # print(num_excluded_radio <- nrow(cr_outcome_an) - nrow(cr_outcome_radio))
    # 
    # #adjusted for lrt_test
    # 
    # cox_model_adj_radio <- tryCatch({
    #   coxph(data = cr_outcome_radio, formula = adjusted_params)
    # }, warning = function(w) {
    #   cat("Warning in fitting Cox model : ", w$message, "\n")
    #   return(NULL)  # Return NULL so the script can continue
    # }, error = function(e) {
    #   # Handle errors in case they occur
    #   cat("Error in fitting Cox model: ", e$message, "\n")
    #   return(NULL)
    # })
    # 
    # print(summary(cox_model_adj_radio))
    # 
    # 
    # 
    # radio_params <- as.formula(
    #   paste(
    #     "Surv(total_fup_var,",
    #     event_var,
    #     ") ~ exposed_radio + strata(setid) + ",
    #     paste(covariates, collapse = "+")
    #   )
    # )
    # 
    # cox_model_radio <- tryCatch({
    #   coxph(data = cr_outcome_radio, formula = radio_params)
    # }, warning = function(w) {
    #   cat("Warning in fitting Cox model : ", w$message, "\n")
    #   return(NULL)  # Return NULL so the script can continue
    # }, error = function(e) {
    #   # Handle errors in case they occur
    #   cat("Error in fitting Cox model: ", e$message, "\n")
    #   return(NULL)
    # })
    # 
    # # Check if the model ran successfully or not
    # if (is.null(cox_model_radio)) {
    #   cat("Cox model did not converge. Continuing with the rest of the script...\n")
    #   results <- rbind(
    #     results,
    #     data.frame(
    #       model_type = c("Radio No Treatment", "Radio Treatment"),
    #       cancer = c(site, site),
    #       outcome = c(outcome, outcome),
    #       hr = c(NA, NA),
    #       ci_lower = c(NA, NA),
    #       ci_upper = c(NA, NA),
    #       p_value = c(NA, NA),
    #       lrtest = c(NA, NA)
    #     )
    #   )
    # } else {
    #   print(summary(cox_model_radio))
    #   
    #   # For exposed_radio = 1 (No Treatment) vs 0 (Cancer-free - reference)
    #   hr_no_treatment <- exp(coef(cox_model_radio)["exposed_radioNo Treatment"])
    #   ci_no_treatment <- exp(confint(cox_model_radio)["exposed_radioNo Treatment", ])
    #   p_no_treatment <- summary(cox_model_radio)$coefficients["exposed_radioNo Treatment", "Pr(>|z|)"]
    #   # For exposed_radio = 2 (Treatment) vs 0 (Cancer-free - reference)
    #   hr_treatment <- exp(coef(cox_model_radio)["exposed_radioTreatment"])
    #   ci_treatment <- exp(confint(cox_model_radio)["exposed_radioTreatment", ])
    #   p_treatment <- summary(cox_model_radio)$coefficients["exposed_radioTreatment", "Pr(>|z|)"]
    #   
    #   
    #   ###  Likelihood ratio test between cox_model_adj and cox_model_chemo
    #   lr_test_radio <- tryCatch({
    #     anova(cox_model_adj_radio, cox_model_radio, test = "LRT")
    #   }, warning = function(w) {
    #     cat("Warning in LRT for radio: ", w$message, "\n")
    #     return(NULL)
    #   }, error = function(e) {
    #     cat("Error in LRT for radio: ", e$message, "\n")
    #     return(NULL)
    #   })
    #   
    #   # Extract p-value safely
    #   if (!is.null(lr_test_radio)) {
    #     print("Likelihood Ratio Test - radio vs Adjusted Model")
    #     print(lr_test_radio)
    #     lrt_p_radio <- lr_test_radio$`Pr(>|Chi|)`[2]
    #   } else {
    #     cat("LRT could not be calculated for radio model\n")
    #     lrt_p_radio <- NA
    #   }
    #   
    #   
    #   
    #   
    #   results <- rbind(
    #     results,
    #     data.frame(
    #       model_type = c("Radio No Treatment", "Radio Treatment"),
    #       cancer = site,
    #       outcome = outcome,
    #       hr = c(hr_no_treatment, hr_treatment),
    #       ci_lower = c(ci_no_treatment[1], ci_treatment[1]),
    #       ci_upper = c(ci_no_treatment[2], ci_treatment[2]),
    #       p_value = c(p_no_treatment, p_treatment),
    #       lrtest = c(
    #         lr_test_radio$`Pr(>|Chi|)`[2],
    #         lr_test_radio$`Pr(>|Chi|)`[2]
    #       )
    #     )
    #   )
    #   
    # } #else
    # 
    # 
    
    if (site == "nhl") {

      
    ################################################################################
    # COX MODEL - HSCT
    ################################################################################
   
    
     #Drop prior to Unexposed
    
    cr_outcome_hsct <- cr_outcome_an %>% 
      filter(excl_hsct == 0)
    
    #print number excluded:
    print(num_excluded_hsct <- nrow(cr_outcome_an) - nrow(cr_outcome_hsct))
    
    #adjusted for lrt_test
    
    cox_model_adj_hsct <- tryCatch({
      coxph(data = cr_outcome_hsct, formula = adjusted_params)
    }, warning = function(w) {
      cat("Warning in fitting Cox model : ", w$message, "\n")
      return(NULL)  # Return NULL so the script can continue
    }, error = function(e) {
      # Handle errors in case they occur
      cat("Error in fitting Cox model: ", e$message, "\n")
      return(NULL)
    })
    
    print(summary(cox_model_adj_hsct))
    
    
    
    hsct_params <- as.formula(
      paste(
        "Surv(total_fup_var,",
        event_var,
        ") ~ exposed_hsct + strata(setid) + ",
        paste(covariates, collapse = "+")
      )
    )
    
    cox_model_hsct <- tryCatch({
      coxph(data = cr_outcome_hsct, formula = hsct_params)
    }, warning = function(w) {
      cat("Warning in fitting Cox model : ", w$message, "\n")
      return(NULL)  # Return NULL so the script can continue
    }, error = function(e) {
      # Handle errors in case they occur
      cat("Error in fitting Cox model: ", e$message, "\n")
      return(NULL)
    })
    
    # Check if the model ran successfully or not
    if (is.null(cox_model_hsct)) {
      cat("Cox model did not converge. Continuing with the rest of the script...\n")
      results <- rbind(
        results,
        data.frame(
          model_type = c("HSCT No Treatment", "HSCT Treatment"),
          cancer = c(site, site),
          outcome = c(outcome, outcome),
          hr = c(NA, NA),
          ci_lower = c(NA, NA),
          ci_upper = c(NA, NA),
          p_value = c(NA, NA),
          lrtest = c(NA, NA)
        )
      )
    } else {
      print(summary(cox_model_hsct))
      
      # For exposed_hsct = 1 (No Treatment) vs 0 (Cancer-free - reference)
      hr_no_treatment <- exp(coef(cox_model_hsct)["exposed_hsctNo Treatment"])
      ci_no_treatment <- exp(confint(cox_model_hsct)["exposed_hsctNo Treatment", ])
      p_no_treatment <- summary(cox_model_hsct)$coefficients["exposed_hsctNo Treatment", "Pr(>|z|)"]
      # For exposed_hsct = 2 (Treatment) vs 0 (Cancer-free - reference)
      hr_treatment <- exp(coef(cox_model_hsct)["exposed_hsctTreatment"])
      ci_treatment <- exp(confint(cox_model_hsct)["exposed_hsctTreatment", ])
      p_treatment <- summary(cox_model_hsct)$coefficients["exposed_hsctTreatment", "Pr(>|z|)"]
      
      
      ###  Likelihood ratio test between cox_model_adj and cox_model_chemo
      lr_test_hsct <- tryCatch({
        anova(cox_model_adj_hsct, cox_model_hsct, test = "LRT")
      }, warning = function(w) {
        cat("Warning in LRT for hsct: ", w$message, "\n")
        return(NULL)
      }, error = function(e) {
        cat("Error in LRT for hsct: ", e$message, "\n")
        return(NULL)
      })
      
      # Extract p-value safely
      if (!is.null(lr_test_hsct)) {
        print("Likelihood Ratio Test - hsct vs Adjusted Model")
        print(lr_test_hsct)
        lrt_p_hsct <- lr_test_hsct$`Pr(>|Chi|)`[2]
      } else {
        cat("LRT could not be calculated for hsct model\n")
        lrt_p_hsct <- NA
      }
      
      
      
      
      results <- rbind(
        results,
        data.frame(
          model_type = c("hsct No Treatment", "hsct Treatment"),
          cancer = site,
          outcome = outcome,
          hr = c(hr_no_treatment, hr_treatment),
          ci_lower = c(ci_no_treatment[1], ci_treatment[1]),
          ci_upper = c(ci_no_treatment[2], ci_treatment[2]),
          p_value = c(p_no_treatment, p_treatment),
          lrtest = c(
            lr_test_hsct$`Pr(>|Chi|)`[2],
            lr_test_hsct$`Pr(>|Chi|)`[2]
          )
        )
      )
      
    } #else
    }
    
  } #outcome
} #cancer

results <- unique(results)

rownames(results) <- NULL
print(results)

#save summary tables as rds

saveRDS(results,
        file = paste0(
          path_datafiles_for_analysis,
          "an_cox_models_first_dx_trt_nhl_hsct.rds"
        ))

# remove all except cr_finaldataforanalysis_respiratory
#rm(list=setdiff(ls(), "cr_finaldataforanalysis_respiratory"))
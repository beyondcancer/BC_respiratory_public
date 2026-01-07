###############################################################################
#3. FIRST_DX Object for survival analysis
##############################################################################

# # ########################
# # ##FOR TESTING PURPOSES##
# outcomes_chronic <- "asthma"
# outcome <- "asthma"
# cancersites <- "col"
# site <- "col"
# # # cancersites <- c("ora")
# # # outcomes_chronic <- c("asthma")

###############################################################################
#1. Load in data

#cr_finaldataforanalysis_respiratory 

###############################################################################

##Initiate dataset to save the exclusions

population <- data.frame(cancer = character(),
                         outcome = character(),
                         with_b_exp = numeric(),
                         with_b_unexp =numeric(), 
                         npop_exp = numeric(),
                         npop_unexp = numeric(),
                         inv_fup_exp = numeric(),
                         inv_fup_unexp = numeric(),
                         uncertain_dx = numeric(),
                         n_events_exp = numeric(), 
                         n_events_unexp = numeric(),
                         t_fup_exp = numeric(),
                         t_fup_unexp = numeric(),
                         t_fup_tot = numeric(),
                         IR_exp = numeric(),
                         IR_unexp = numeric(),
                         stringsAsFactors = FALSE)


## Initiate cancer type loop
print(paste0("Initiating loop for cancer"))

   for(site in cancersites) {
    
    print(paste0("cancer type: ", site))
    
    cr_cancer_data <- cr_finaldataforanalysis_respiratory %>%
      filter(cancer == site) 
    
   print(paste0("Number of  patients within the ", site, "cohort: ", count(cr_cancer_data)))
    

    #Initiate outcome loop  
    print(paste0("Initiating loop for "))
      
      for (outcome in outcomes_chronic) {
        
        print(paste0("Outcome: ", outcome))
        #Create dynamic variables  
        history_var <- paste0("b_", outcome)  # e.g., b_asthma
        type_var <- paste0("first_event_", outcome, "_type")
        date_event_var <- paste0("dof_any_", outcome, "_inc_dx")  # e.g., dof_asthma_inc_dx
        event_var <- paste0("any_", outcome, "_inc_dx")  # e.g., asthma_inc_dx
        
        print(event_var)
        print(table(cr_cancer_data$exposed, cr_cancer_data[[event_var]]))
        
    #Count for exclusions
        count_with_history <- sum(cr_cancer_data[[history_var]] == 1)
        count_with_history_unexp <- sum(cr_cancer_data[[history_var]] == 1 & cr_cancer_data$exposed == 0)
        count_with_history_exp <- sum(cr_cancer_data[[history_var]] == 1 & cr_cancer_data$exposed == 1)
        
      
        print(paste0("People with a history of outcome: ", count_with_history))
        print(paste0("People with a history of outcome in the unexposed group: ", count_with_history_unexp))
        print(paste0("People with a history of outcome in the exposed group: ", count_with_history_exp))
        print(paste0("percent of total that have b_asthma = 1: ", (count_with_history/nrow(cr_cancer_data) * 100)))
      
        
      # Apply outcome-specific exclusions
      cr_outcome_an <- cr_cancer_data %>%
      filter(.data[[history_var]] == 0)  # exclude those with a history of outcome
      
      table(cr_outcome_an[[type_var]], cr_outcome_an[[event_var]])
      
      print(table(cr_outcome_an$exposed, cr_outcome_an[[event_var]], useNA = "ifany"))
      n_events_b4prev <- sum(cr_outcome_an[[event_var]] == 1)
      
      cr_outcome_an <- cr_outcome_an %>% mutate(
        uncertain_dx = if_else(.data[[type_var]] == "prev" & .data[[event_var]] == 1, 1, 0),
        !!sym(event_var) := if_else(uncertain_dx == 1, 0, !!sym(event_var)))

      print(table(cr_outcome_an$exposed, cr_outcome_an[[event_var]], useNA = "ifany"))
      
    # Modify censor date to censor at date of first event diagnosis
    cr_outcome_an <- cr_outcome_an %>% 
      mutate(doexit_var = pmin(doexit, .data[[date_event_var]], na.rm= TRUE)) #select earliest date between doexit and event_var
    
    # Calculate total follow-up time (in years)
    cr_outcome_an <- cr_outcome_an %>% 
      mutate(total_fup_days = pmax((doexit_var - indexdate), 0.5)) %>%
      mutate(total_fup_var = as.numeric(total_fup_days/365.25)) #calculate total follow-up time (years)
    
    count_excluded_f_up_exp <- sum(cr_outcome_an$exposed == 1 & cr_outcome_an$total_fup_var < 0)
    count_excluded_f_up_unexp <- sum(cr_outcome_an$exposed == 0 & cr_outcome_an$total_fup_var < 0)
    count_excluded_f_up <- sum(cr_outcome_an$total_fup_var < 0)
    count_uncertain_dx <- sum(cr_outcome_an$uncertain_dx == 1)
    print(paste0("Number of people with invalid follow-up time: ", count_excluded_f_up))
    print(paste0("Percent of incident cases censored because of uncertain dx: ", 
                 round( (count_uncertain_dx / 
                          n_events_b4prev) * 100, 2), 
                 "%"))    
    print(table(cr_outcome_an$exposed, cr_outcome_an[[event_var]], useNA = "ifany"))
    
    #Recode event_var to 0 if doexit_var is < date_event_var and event_var is 1
    cr_outcome_an <- cr_outcome_an %>%
      mutate(!!event_var := if_else(doexit < doexit_var, 0, .data[[event_var]]))
    
    print(paste0("There are ", count(cr_outcome_an)," patients once outcome specific exclusions are applied"))

    covariates <- c( "age_cat", "gender", "imd5", "bmi", "bmi_cat", "smokstatus", "b_asthma", "b_copd", "b_bronch", "b_ild", "b_resp_inf",
                    "b_ckd", "b_hypertension", "b_cvdgrouped", "b_diab_cat", "b_cld", "b_multiple_sclerosis")

    #Complete case analysis only: Drop observations with missing data in the covariates included in the model
    missing_summary <- colSums(is.na(cr_outcome_an[ , all_of(covariates)]))
    print(missing_summary)
    count_missing <- sum(apply(cr_outcome_an[ , all_of(covariates)], 1, function(x) any(is.na(x))))
    
    #Exclusions numbers
    print(paste0("Number of missing values: ", count_missing))    
    print(paste0("Number of people with invalid follow-up time: ", count_excluded_f_up))
    
    #Keep only complete cases
    cr_outcome_an <- cr_outcome_an %>%
      filter(if_all(all_of(covariates), ~ !is.na(.))) # complete case
    
    #Inclusion numbers
    print(paste0("Total numbers for complete case analysis: ", count(cr_outcome_an)))    
    hist(cr_outcome_an[[date_event_var]], breaks = 20, freq = T, main = paste0("Histogram of ", date_event_var))
    print(event_var)
    table(cr_outcome_an[[event_var]])
   
    print(table(cr_outcome_an$exposed, cr_outcome_an[[event_var]], useNA = "ifany"))
    
     
    #Keep relevant variables
    cr_outcome_an <- cr_outcome_an %>%
      dplyr::select(setid, e_patid, exposed, indexdate, doexit_var, total_fup_var, event_var, type_var, date_event_var, eth5_hes, eth5_cprd, all_of(covariates), uncertain_dx) 
    
    
    
    #Convert covariates to factors
    factors <- c("age_cat", "gender", "imd5", "bmi", "bmi_cat", "eth5_hes", "eth5_cprd", "smokstatus")
    cr_outcome_an <- cr_outcome_an %>%
    mutate(across(all_of(factors), as.factor))
    
    #Count number of outcome events in exposed and unexposed group and save in a variable 
    n_exp <- cr_outcome_an %>% filter(exposed == 1 & .data[[event_var]] == 1) %>% nrow()
    n_unexp <- cr_outcome_an %>% filter(exposed == 0 & .data[[event_var]] == 1) %>% nrow()

    # Save the dataset for each site and outcome
    file_name <- paste0(path_datafiles_for_analysis, "sens4/", site, "_", outcome, "_cr_outcome_an.csv")
    write.csv(cr_outcome_an, file_name, row.names = FALSE)
    print(paste0("Dataset saved to: ", file_name))
    
    #table of counts of event_var by exposed
    print("2x2 table exposed / outcome")
    print(table(cr_outcome_an$exposed, cr_outcome_an[[event_var]], useNA = "ifany"))

    print(colnames(cr_outcome_an))
    
    fup <- cr_outcome_an %>%
      group_by(exposed) %>%
      summarise(mean_fup = mean(total_fup_var, na.rm = TRUE),
                median_fup = median(!!sym(date_event_var), na.rm = TRUE),
                min_fup = min(!!sym(date_event_var), na.rm = TRUE),
                max_fup = max(!!sym(date_event_var), na.rm = TRUE))
    print("Follow-up time statistics")
    print(fup)
    
    t_fup_exp <- cr_outcome_an %>%
      filter(exposed == 1) %>%  # Keep only exposed individuals
      summarise(total_fup_exp = sum(total_fup_var, na.rm = TRUE)) %>%  # Summarize with sum()
      pull(total_fup_exp)  # Extract as a single value
    
    t_fup_unexp <- cr_outcome_an %>%
      filter(exposed == 0) %>%  # Keep only exposed individuals
      summarise(total_fup_unexp = sum(total_fup_var, na.rm = TRUE)) %>%  # Summarize with sum()
      pull(total_fup_unexp)  # Extract as a single value
    
    
    ##Complete dataset:
    population <- rbind(population, data.frame(cancer = site,
                                               outcome = outcome,
                                               with_b_exp = count_with_history_exp,
                                               with_b_unexp = count_with_history_unexp,
                                               inv_fup_exp = count_excluded_f_up_exp,
                                               inv_fup_unexp = count_excluded_f_up_unexp,
                                               uncertain_dx = sum(cr_outcome_an$uncertain_dx == 1, na.rm = TRUE),
                                               n_pop_exp = sum(cr_outcome_an$exposed == 1),
                                               n_pop_unexp = sum(cr_outcome_an$exposed == 0),
                                               n_events_exp = n_exp, 
                                               n_events_unexp = n_unexp,
                                               t_fup_exp = t_fup_exp,
                                               t_fup_unexp = t_fup_unexp,
                                               t_fup_tot = sum(cr_outcome_an$total_fup_var, na.rm = TRUE),
                                               IR_exp = (n_exp/t_fup_exp) * 1000,
                                               IR_unexp = (n_unexp/t_fup_unexp) * 1000
                                              ))
    
  } #end of outcome loop
} #end of cancer loop


population <- unique(population)

#save summary tables as csv
write.csv(population, paste0(path_results, "/first_dx_analysis_numbers_.csv"), row.names = FALSE)


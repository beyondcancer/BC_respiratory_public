################################################################################
# Landmark analysis: Create outcome datasets
#
# Adapts the standard outcome dataset creation script for the landmark analysis.
#   - Time 0 is the landmark point (new_index = indexdate + 1 year)
#   - Follow-up is measured from new_index, not indexdate
#   - Exposure is treatment receipt (exposed_chemo / exposed_radio) rather
#     than cancer status
#   - Only patients who survived to the landmark are included (already applied
#     in cr_treatments via filter(total_fup_days >= 365.25))
#   - Events occurring before new_index are excluded
#   - Matched sets where the exposed individual or all unexposed individuals
#     have been lost are removed
################################################################################

# # ########################
# ##FOR TESTING PURPOSES##
# cancersites <- c("bre", "lun", "oes", "nhl")
cancersites <- c("nhl")

# outcomes_chronic <- c("asthma")

# site <- c("leu")
# outcome <- c("asthma")
# # ########################

################################################################################
# Read in file
################################################################################

cr_treatments <- readRDS(paste0(path_datafiles_for_analysis, "cr_treatmentanalysis.rds"))

################################################################################
# Initialise population summary table
################################################################################

population_landmark <- data.frame(
  cancer                        = character(),
  outcome                       = character(),
  npop_exposed                  = numeric(),
  npop_unexposed                = numeric(),
  n_events_exposed              = numeric(),
  n_events_unexposed            = numeric(),
  mean_fup_exposed              = numeric(),
  mean_fup_unexposed            = numeric(),
  median_fup_exposed            = numeric(),
  median_fup_unexposed          = numeric(),
  npop_chemo_treatment          = numeric(),
  npop_chemo_no_treatment       = numeric(),
  npop_chemo_cancer_free        = numeric(),
  n_events_chemo_treatment      = numeric(),
  n_events_chemo_no_treatment   = numeric(),
  n_events_chemo_cancer_free    = numeric(),
  mean_fup_chemo_treatment      = numeric(),
  mean_fup_chemo_no_treatment   = numeric(),
  mean_fup_chemo_cancer_free    = numeric(),
  median_fup_chemo_treatment    = numeric(),
  median_fup_chemo_no_treatment = numeric(),
  median_fup_chemo_cancer_free  = numeric(),
  npop_radio_treatment          = numeric(),
  npop_radio_no_treatment       = numeric(),
  npop_radio_cancer_free        = numeric(),
  n_events_radio_treatment      = numeric(),
  n_events_radio_no_treatment   = numeric(),
  n_events_radio_cancer_free    = numeric(),
  mean_fup_radio_treatment      = numeric(),
  mean_fup_radio_no_treatment   = numeric(),
  mean_fup_radio_cancer_free    = numeric(),
  median_fup_radio_treatment    = numeric(),
  median_fup_radio_no_treatment = numeric(),
  median_fup_radio_cancer_free  = numeric(),
  stringsAsFactors              = FALSE
)

################################################################################
# Loops
################################################################################

print(paste0("Initiating loop for cancer"))

for (site in cancersites) {
  
  print(paste0("Cancer type: ", site))
  
  cr_cancer_data <- cr_treatments %>%
    filter(cancer == site)
  
  print(paste0("Initiating loop for outcomes"))
  
  for (outcome in outcomes_chronic) {
    
    print(paste0("Outcome: ", outcome))
    
    # --- Define dynamic variable names ---
    history_var    <- paste0("b_", outcome)
    type_var       <- paste0("first_event_", outcome, "_type")
    date_event_var <- paste0("dof_any_", outcome, "_inc_dx")
    event_var      <- paste0("any_", outcome, "_inc_dx")
    
    # --- Exclude those with a history of the outcome ---
    cr_outcome_an <- cr_cancer_data %>%
      filter(.data[[history_var]] == 0)
    
    # --- Handle uncertain diagnoses ---
    cr_outcome_an <- cr_outcome_an %>%
      mutate(
        uncertain_dx   = if_else(.data[[type_var]] == "prev" & .data[[event_var]] == 1, 1, 0),
        !!sym(event_var) := if_else(uncertain_dx == 1, 0, !!sym(event_var))
      )
    
    # -----Remove incident events in the first - year ---------
     cr_outcome_an <- cr_outcome_an %>%
      mutate(
        event_before_landmark = if_else(!is.na(.data[[date_event_var]]) & .data[[date_event_var]] < new_index, 1, 0),
        !!sym(event_var) := if_else(event_before_landmark == 1, 0, !!sym(event_var))
      ) %>%
      filter(event_before_landmark == 0)
    
    # --- Censor date: earliest of doexit or event date ---
    cr_outcome_an <- cr_outcome_an %>%
      mutate(doexit_var = pmin(doexit, .data[[date_event_var]], na.rm = TRUE))
    
    # --- Follow-up time measured from new_index (landmark time point) ---
    cr_outcome_an <- cr_outcome_an %>%
      mutate(
        total_fup_days = pmax((doexit_var - new_index), 0.5),
        total_fup_var  = as.numeric(total_fup_days / 365.25)
      )
    
    # --- Recode event if doexit precedes doexit_var ---
    cr_outcome_an <- cr_outcome_an %>%
      mutate(!!event_var := if_else(doexit < doexit_var, 0, .data[[event_var]]))
    
    # --- Covariates ---
    # Remove outcome history variable 
    covariates <- c("age_cat", "gender", "imd5", "bmi", "bmi_cat", "smokstatus",
                    "b_asthma", "b_copd", "b_bronch", "b_ild", "b_resp_inf",
                    "b_ckd", "b_hypertension", "b_cvdgrouped", "b_diab_cat",
                    "b_cld", "b_multiple_sclerosis")
    covariates <- setdiff(covariates, history_var)
    
    # --- Complete case ---
    cr_outcome_an <- cr_outcome_an %>%
      filter(if_all(all_of(covariates), ~ !is.na(.)))
    
    # --- Select relevant variables ---
    cr_outcome_an <- cr_outcome_an %>%
      dplyr::select(setid, e_patid, exposed, exposed_chemo, exposed_radio, exposed_hsct, excl_radio, excl_chemo, excl_hsct,
                    new_index, indexdate, doexit_var, total_fup_var,
                    all_of(event_var), all_of(type_var), all_of(date_event_var),
                    eth5_hes, eth5_cprd, all_of(covariates), uncertain_dx)
    
    # --- Convert covariates to factors ---
    factors <- c("age_cat", "gender", "imd5", "bmi", "bmi_cat",
                 "eth5_hes", "eth5_cprd", "smokstatus")
    cr_outcome_an <- cr_outcome_an %>%
      mutate(across(all_of(factors), as.factor))

    
    # --- Event counts by exposure category ---
    n_chemo_exp_trt   <- sum(cr_outcome_an$exposed_chemo == "Treatment"    & cr_outcome_an[[event_var]] == 1, na.rm = TRUE)
    n_chemo_exp_notrt <- sum(cr_outcome_an$exposed_chemo == "No Treatment" & cr_outcome_an[[event_var]] == 1, na.rm = TRUE)
    n_chemo_unexp     <- sum(cr_outcome_an$exposed_chemo == "Cancer-free"  & cr_outcome_an[[event_var]] == 1, na.rm = TRUE)
    
    n_radio_exp_trt   <- sum(cr_outcome_an$exposed_radio == "Treatment"    & cr_outcome_an[[event_var]] == 1, na.rm = TRUE)
    n_radio_exp_notrt <- sum(cr_outcome_an$exposed_radio == "No Treatment" & cr_outcome_an[[event_var]] == 1, na.rm = TRUE)
    n_radio_unexp     <- sum(cr_outcome_an$exposed_radio == "Cancer-free"  & cr_outcome_an[[event_var]] == 1, na.rm = TRUE)
    
    n_hsct_exp_trt   <- sum(cr_outcome_an$exposed_hsct == "Treatment"    & cr_outcome_an[[event_var]] == 1, na.rm = TRUE)
    n_hsct_exp_notrt <- sum(cr_outcome_an$exposed_hsct == "No Treatment" & cr_outcome_an[[event_var]] == 1, na.rm = TRUE)
    n_hsct_unexp     <- sum(cr_outcome_an$exposed_hsct == "Cancer-free"  & cr_outcome_an[[event_var]] == 1, na.rm = TRUE)
    
    # --- Save dataset ---
    file_name <- paste0(path_datafiles_for_analysis, "trt/", site, "_", outcome, "_cr_outcome_an.csv")
    write.csv(cr_outcome_an, file_name, row.names = FALSE)
    print(paste0("Dataset saved to: ", file_name))
    
    # --- Follow-up time statistics ---
    print("Follow-up time statistics - Chemo")
    print(table(cr_outcome_an$exposed_chemo, cr_outcome_an[[event_var]]))
    
    print("Follow-up time statistics - Radio")
    print(table(cr_outcome_an$exposed_radio, cr_outcome_an[[event_var]]))
    
     print("Follow-up time statistics - HSCT")
    print(table(cr_outcome_an$exposed_hsct, cr_outcome_an[[event_var]]))
    
    
    # --- Population summary: numbers and follow-up by exposure category ---
    population_landmark <- rbind(population_landmark, data.frame(
      
      cancer  = site,
      outcome = outcome,
      
      # --- Cancer exposure (exposed vs unexposed) ---
      npop_exposed   = sum(cr_outcome_an$exposed == 1, na.rm = TRUE),
      npop_unexposed = sum(cr_outcome_an$exposed == 0, na.rm = TRUE),
      
      n_events_exposed   = sum(cr_outcome_an$exposed == 1 & cr_outcome_an[[event_var]] == 1, na.rm = TRUE),
      n_events_unexposed = sum(cr_outcome_an$exposed == 0 & cr_outcome_an[[event_var]] == 1, na.rm = TRUE),
      
      mean_fup_exposed   = round(mean(cr_outcome_an$total_fup_var[cr_outcome_an$exposed == 1],   na.rm = TRUE), 2),
      mean_fup_unexposed = round(mean(cr_outcome_an$total_fup_var[cr_outcome_an$exposed == 0],   na.rm = TRUE), 2),
      
      median_fup_exposed   = round(median(cr_outcome_an$total_fup_var[cr_outcome_an$exposed == 1], na.rm = TRUE), 2),
      median_fup_unexposed = round(median(cr_outcome_an$total_fup_var[cr_outcome_an$exposed == 0], na.rm = TRUE), 2),
      
      # --- Chemo exposure: population size ---
      npop_chemo_treatment    = sum(cr_outcome_an$exposed_chemo == "Treatment",    na.rm = TRUE),
      npop_chemo_no_treatment = sum(cr_outcome_an$exposed_chemo == "No Treatment", na.rm = TRUE),
      npop_chemo_cancer_free  = sum(cr_outcome_an$exposed_chemo == "Cancer-free",  na.rm = TRUE),
      
      # --- Chemo exposure: event counts ---
      n_events_chemo_treatment    = n_chemo_exp_trt,
      n_events_chemo_no_treatment = n_chemo_exp_notrt,
      n_events_chemo_cancer_free  = n_chemo_unexp,
      
      # --- Chemo exposure: mean and median follow-up ---
      mean_fup_chemo_treatment    = round(mean(cr_outcome_an$total_fup_var[cr_outcome_an$exposed_chemo == "Treatment"],    na.rm = TRUE), 2),
      mean_fup_chemo_no_treatment = round(mean(cr_outcome_an$total_fup_var[cr_outcome_an$exposed_chemo == "No Treatment"], na.rm = TRUE), 2),
      mean_fup_chemo_cancer_free  = round(mean(cr_outcome_an$total_fup_var[cr_outcome_an$exposed_chemo == "Cancer-free"],  na.rm = TRUE), 2),
      
      median_fup_chemo_treatment    = round(median(cr_outcome_an$total_fup_var[cr_outcome_an$exposed_chemo == "Treatment"],    na.rm = TRUE), 2),
      median_fup_chemo_no_treatment = round(median(cr_outcome_an$total_fup_var[cr_outcome_an$exposed_chemo == "No Treatment"], na.rm = TRUE), 2),
      median_fup_chemo_cancer_free  = round(median(cr_outcome_an$total_fup_var[cr_outcome_an$exposed_chemo == "Cancer-free"],  na.rm = TRUE), 2),
      
      # --- Radio exposure: population size ---
      npop_radio_treatment    = sum(cr_outcome_an$exposed_radio == "Treatment",    na.rm = TRUE),
      npop_radio_no_treatment = sum(cr_outcome_an$exposed_radio == "No Treatment", na.rm = TRUE),
      npop_radio_cancer_free  = sum(cr_outcome_an$exposed_radio == "Cancer-free",  na.rm = TRUE),
      
      # --- Radio exposure: event counts ---
      n_events_radio_treatment    = n_radio_exp_trt,
      n_events_radio_no_treatment = n_radio_exp_notrt,
      n_events_radio_cancer_free  = n_radio_unexp,
      
      # --- Radio exposure: mean and median follow-up ---
      mean_fup_radio_treatment    = round(mean(cr_outcome_an$total_fup_var[cr_outcome_an$exposed_radio == "Treatment"],    na.rm = TRUE), 2),
      mean_fup_radio_no_treatment = round(mean(cr_outcome_an$total_fup_var[cr_outcome_an$exposed_radio == "No Treatment"], na.rm = TRUE), 2),
      mean_fup_radio_cancer_free  = round(mean(cr_outcome_an$total_fup_var[cr_outcome_an$exposed_radio == "Cancer-free"],  na.rm = TRUE), 2),
      
      median_fup_radio_treatment    = round(median(cr_outcome_an$total_fup_var[cr_outcome_an$exposed_radio == "Treatment"],    na.rm = TRUE), 2),
      median_fup_radio_no_treatment = round(median(cr_outcome_an$total_fup_var[cr_outcome_an$exposed_radio == "No Treatment"], na.rm = TRUE), 2),
      median_fup_radio_cancer_free  = round(median(cr_outcome_an$total_fup_var[cr_outcome_an$exposed_radio == "Cancer-free"],  na.rm = TRUE), 2)
      
    ))
    
  } #end of outcome loop
} #end of cancer loop

################################################################################
# Save population summary
################################################################################

population_landmark <- unique(population_landmark)
rownames(population_landmark) <- NULL

print(population_landmark)
# 
# write.csv(population_landmark,
#           paste0(path_results, "/trt/landmark_analysis_numbers.csv"),
#           row.names = FALSE)

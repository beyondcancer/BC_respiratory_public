################################################################################
# Landmark analysis: treatment exposure variables for respiratory cancer cohort
#
# Landmark analysis approach:
#   - The landmark time point is set at one year post cancer diagnosis (index date)
#   - Patients who die or are lost to follow-up before the landmark are excluded
#   - Treatment exposures (chemotherapy, radiotherapy) are assessed in the
#     window between cancer diagnosis and the landmark (0 to 365.25 days)
#   - new_index marks the landmark time point, from which survival outcomes
#     are measured
#
# For each treatment modality, a composite exposure variable is created:
#   0 = Cancer-free (unexposed cohort)
#   1 = Cancer diagnosis, no treatment received before landmark
#   2 = Cancer diagnosis, treatment received before landmark
#
# Note: HSCT exposure is coded but commented out pending data availability
################################################################################

# --- Load data ---
cr_finaldataforanalysis_respiratory <- readRDS(paste0(path_datafiles_for_analysis, 
                                                      "cr_finaldataforanalysis_respiratory.rds"))


hsct <- read_dta(paste0(path_datafiles, "/cr_listpat_hsct_aurum.dta"))

cr_finaldataforanalysis_respiratory <- left_join(cr_finaldataforanalysis_respiratory, hsct) %>% 
  mutate(hsct = ifelse(is.na(hsct), 0, hsct))


# check <- cr_finaldataforanalysis_respiratory %>% filter(cancer == "nhl") %>% group_by(exposed) %>% summarise(n = n(),
#                                                         chemo = sum(!is.na(dof_chemo)),
#                                                         radio = sum(!is.na(dof_radio)),
#                                                         hsct = sum(!is.na(dof_hsct)),
#                                                         percent_chemo = round(chemo/n*100, 1),
#                                                         percent_radio = round(radio/n*100, 1),
#                                                         percent_hsct = round(hsct/n*100, 1))
# 
# 


cr_treatments <- cr_finaldataforanalysis_respiratory %>%
  
  # --- Flag treatment exposures within landmark window (0 to 1 year) ---
  # Patients with no recorded event (NA date) are coded as 0 (not exposed)
  mutate(
    chemo = as.integer(if_else(!is.na(dof_chemo), dof_chemo >= indexdate & dof_chemo < indexdate + 365.25, FALSE)),
    radio = as.integer(if_else(!is.na(dof_radio), dof_radio >= indexdate & dof_radio < indexdate + 365.25, FALSE)),
    hsct = as.integer(if_else(!is.na(dof_hsct), dof_hsct >= indexdate & dof_hsct < indexdate + 365.25, FALSE))
  ) %>%
  
  # --- Apply landmark exclusion criterion ---
  # Exclude patients who did not survive to the landmark time point (1 year),
  # as their treatment exposure window is incomplete
  filter(total_fup_days >= 365.25) %>%
  
  mutate(
    # --- Define landmark time point ---
    # new_index is the landmark date (1 year post diagnosis), from which
    # survival outcomes will be measured
    new_index = as.Date(indexdate + 365.25),
    
    # --- Composite exposure variables ---
    # Combines cancer exposure status with treatment receipt in landmark window:
    #   0 = Cancer-free (unexposed cohort)
    #   1 = Cancer diagnosis, no treatment received before landmark
    #   2 = Cancer diagnosis, treatment received before landmark
    exposed_chemo = case_when(
      exposed == 1 & chemo == 1 ~ 2,
      exposed == 1 & chemo == 0 ~ 1,
      exposed == 0              ~ 0
    ),
    
    exposed_radio = case_when(
      exposed == 1 & radio == 1 ~ 2,
      exposed == 1 & radio == 0 ~ 1,
      exposed == 0              ~ 0
    ),
    
    exposed_hsct = case_when(
      exposed == 1 & hsct == 1 ~ 2,
      exposed == 1 & hsct == 0 ~ 1,
      exposed == 0             ~ 0
    )
  ) %>%
  
  # --- Label composite exposure variables as factors ---
  mutate(
    exposed_chemo = factor(exposed_chemo, levels = c(0, 1, 2),
                           labels = c("Cancer-free", "No Treatment", "Treatment")),
    exposed_radio = factor(exposed_radio, levels = c(0, 1, 2),
                           labels = c("Cancer-free", "No Treatment", "Treatment")),
     exposed_hsct = factor(exposed_hsct, levels = c(0, 1, 2),
                           labels = c("Cancer-free", "No Treatment", "Treatment"))
  ) %>%


  #Create exclusion flags for chemo and radio analysis. excl_radio = 1 if indexdate < 01-01-2012 and excl_chemo = 1 if indexdate < 01-01-2014.
  mutate(
    excl_chemo = if_else(indexdate < as.Date("2014-01-01"), 1, 0),
    excl_radio = if_else(indexdate < as.Date("2012-01-01"), 1, 0),
    excl_hsct = if_else(exposed == 0 & hsct == 1, 1, 0)
  )



# --- Sense check ---
# Inspect key variables to verify treatment flagging and exposure coding
check <- cr_treatments %>% filter(cancer == "nhl") %>%
  dplyr::select(e_patid, cancer, indexdate, new_index, exposed, doexit,
                dof_chemo, dof_radio, dof_hsct,
                chemo, radio, hsct, 
                exposed_chemo, exposed_radio, exposed_hsct, exposed, excl_chemo, excl_radio, excl_hsct, doendcprdfup)



gc()
saveRDS(cr_treatments, file = paste0(path_datafiles_for_analysis, "cr_treatmentanalysis.rds"))

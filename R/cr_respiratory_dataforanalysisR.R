##RUN PATHS.R to load core data, libraries and setwd

required_vars <- c(
  "setid", "e_patid", "cancer","exposed", "indexdate", 
  "cprd_db", "e_pracid", "index_year", "index_year_gr",
  "dostartcprdfup", "doendcprdfup", "dod", 
  "gender", "yob", "age", "age_cat", "region", "eth5_cprd", "eth5_hes", "imd5", "bmi", "bmi_cat", "smokstatus",
  "b_ckd", "b_ckd_stage", "b_hypertension", "b_cvdgrouped", "b_diab_cat", "b_multiple_sclerosis",
  "stage_tnm", "stage_leu", "grade_cns","dof_radio", "dof_chemo")

core_dataset <- core_data %>% dplyr::select(all_of(required_vars)) %>% zap_labels()


# Show summary of 'exposed' variable
table(core_dataset$exposed)

# Merge with respiratory outcomes
for (outcome in outcomes_all) {
  outcome_file <- file.path(path_datafiles, paste0("cr_listpat_", outcome, "_aurum_firstdx.dta"))
  outcome_data <- haven::read_dta(outcome_file) %>% zap_labels()
  core_dataset <- core_dataset %>%
    left_join(outcome_data) 
    }


# Chronic liver disease
cld_data <- read_dta(file.path(path_datafiles, "cr_listpat_cld_aurum.dta"))

core_dataset <- core_dataset %>%
  left_join(cld_data) %>%
  rename(b_cld = cld_hx) %>%
  mutate(b_cld = ifelse(is.na(b_cld), 0, b_cld))

core_dataset <- core_dataset %>%
  mutate(
    across(starts_with("dof_"), ~ as.Date(., format = "%d%b%Y")), # Convert "dof_" columns to Date format
    across(starts_with("dol_"), ~ as.Date(., format = "%d%b%Y")), # Convert "dol_" columns to Date format
    across(starts_with("b_"), ~ coalesce(., 0)),  # Replace NA with 0 for "b_" columns
    across(ends_with("_inc_dx") & !starts_with("dof"), ~ coalesce(., 0)) # Replace NA with 0 for "inc_dx" columns
  )


base_cohort <- core_dataset 

print(#see data type of 
  str(dplyr::select(base_cohort, starts_with("dof_") | starts_with("dol_") | starts_with("b_") | ends_with("inc_dx")))
)

#rm(core_data)
rm(core_dataset)

# Create date of pandemic cut-off
base_cohort <- base_cohort %>%
  mutate(pandemic = as.Date("2019-12-31"))

# Calculate follow-up time
base_cohort <- base_cohort %>%
  mutate(doexit = pmin(doendcprdfup, dod, pandemic, na.rm = TRUE),
         time_b4_fup = indexdate - dostartcprdfup) %>%
  mutate(total_fup_days = doexit - indexdate)

#Remove rows with indexdate after doexit
base_cohort <- base_cohort %>% filter(indexdate <= doexit)

# #histogram total_fup days
# hist(as.numeric(base_cohort$total_fup_days), breaks = 100)

# # Recode respiratory variables as binary (1, 0)
# baseline_vars <- c("b_asthma", "b_copd", "b_bronch", "b_ild", "b_flu", "b_pneumo", "b_covid", "b_urti", "b_ckd", "b_ckd_stage", "b_hypertension", "b_cvdgrouped", "b_diab_cat", "b_cld", "b_multiple_sclerosis")
#                    
# prev_vars <- c("asthma_prev", "copd_prev", "bronch_prev", "ild_prev", "flu_prev", "pneumo_prev", "covid_prev", "urti_prev")
# 
# 
# dof_inc_dx_vars <- c("dof_asthma_inc_dx", "dof_copd_inc_dx", "dof_bronch_inc_dx", "dof_ild_inc_dx")
# 
# first_event_vars <- c("first_event_asthma_type", "first_event_copd_type", "first_event_bronch_type", "first_event_ild_type")
# 
# base_cohort <- base_cohort %>%
#   mutate(across(c(all_of(baseline_vars), all_of(prev_vars), all_of(inc_dx_vars)), 
#                 ~ if_else(is.na(.), 0, .)))

#Create variable for baseline history of respiratory infection and for baseline history of chronic respiratory disease
base_cohort <- base_cohort %>%
  mutate(
    b_resp_inf = if_else(b_flu == 1 | b_urti == 1 | b_pneumo == 1, 1, 0),
    b_chronic_resp = if_else(b_asthma == 1 | b_copd == 1 | b_bronch == 1 | b_ild == 1, 1, 0)
  )

# #Ensure that outcome flags are 0 if the outcome date is after end of follow-up

      for (outcome in outcomes_all) {
        date <- paste0("dof_any_", outcome, "_inc_dx")
        event <- paste0("any_", outcome, "_inc_dx")
        
        base_cohort <- base_cohort %>%
          mutate(!!sym(event) := coalesce(if_else(!!sym(date) > doexit, 0, !!sym(event)), 0),
                 !!sym(date) := if_else(!!sym(date) > doexit, as.Date(NA), !!sym(date)))
      }


print(#see data type of 
  str(dplyr::select(base_cohort, starts_with("dof_") | starts_with("dol_") | starts_with("b_") | ends_with("inc_dx")))
)

# inc_dx_vars <- c("asthma_inc_dx", "copd_inc_dx", "bronch_inc_dx", "ild_inc_dx")
# any_inc_dx_vars <- c("any_asthma_inc_dx", "any_copd_inc_dx", "any_bronch_inc_dx", "any_ild_inc_dx")
# hosp_vars <- c("asthma_hosp", "copd_hosp", "bronch_hosp", "ild_hosp")
# death_vars <- c("asthma_death", "copd_death", "bronch_death", "ild_death")
# 
# dof_inc_dx_vars <- c("dof_asthma_inc_dx", "dof_copd_inc_dx", "dof_bronch_inc_dx", "dof_ild_inc_dx")
# dof_inc_dx_vars <- c("dof_asthma_inc_dx", "dof_copd_inc_dx", "dof_bronch_inc_dx", "dof_ild_inc_dx")
# dof_hosp_vars <- c("dof_asthma_hosp", "dof_copd_hosp", "dof_bronch_hosp", "dof_ild_hosp")
# dof_death_vars <- c("dof_asthma_death", "dof_copd_death", "dof_bronch_death", "dof_ild_death")
# 
# for (i in seq_along(outcomes_chronic)) {
#   
#   base_cohort <- base_cohort %>%
#     mutate(!!sym(inc_dx_vars[i]) := if_else(!!sym(dof_inc_dx_vars[i]) > doexit, 0, !!sym(inc_dx_vars[i])),
#            !!sym(hosp_vars[i]) := if_else(!!sym(dof_hosp_vars[i]) > doexit, 0, !!sym(hosp_vars[i])),
#            !!sym(death_vars[i]) := if_else(!!sym(dof_death_vars[i]) > doexit, 0, !!sym(death_vars[i])))
#   
#   base_cohort <- base_cohort %>%
#     mutate(!!sym(inc_dx_vars[i]) := if_else(is.na(!!sym(inc_dx_vars[i])), 0, !!sym(inc_dx_vars[i])),
#            !!sym(hosp_vars[i]) := if_else(is.na(!!sym(hosp_vars[i])), 0, !!sym(hosp_vars[i])),
#            !!sym(death_vars[i]) := if_else(is.na(!!sym(death_vars[i])), 0, !!sym(death_vars[i])))
# }


print(table(base_cohort$any_asthma_inc_dx, base_cohort$exposed, useNA= "ifany"))

# print(table(base_cohort$asthma_hosp, base_cohort$exposed, useNA= "ifany"))
# print(table(base_cohort$asthma_death, base_cohort$exposed, useNA= "ifany"))

# #Check number of dof dates that are outside of the range between indexdate and doexit
# test <- base_cohort %>% 
#   filter(any_asthma_inc_dx == 1 & dof_any_asthma_inc_dx > doexit)
# 
# base_cohort <-  base_cohort %>% mutate(across(starts_with("dof_"), ~ if_else(. > doexit, as.Date(NA), .)))
# 
# test2 <- base_cohort %>%
#   filter(dof_any_asthma_inc_dx > doexit)
# 
# test3 <- base_cohort %>% 
#   filter(dof_any_asthma_inc_dx < indexdate)

#Check number dof dates that are outside of the range between indexdate and doexit
test <- base_cohort %>%
  summarise(across(
    starts_with("dof_"), ~ sum(. < indexdate | . > doexit, na.rm = TRUE)  # Summing occurrences outside indexdate and doexit
  ))

fup <- base_cohort %>%
  group_by(exposed) %>%
  summarise(mean_fup = mean(total_fup_days, na.rm = TRUE),
            median_fup = median(total_fup_days, na.rm = TRUE),
            min_fup = min(total_fup_days, na.rm = TRUE),
            max_fup = max(total_fup_days, na.rm = TRUE))

print(fup)
base_cohort <- base_cohort %>% mutate(flag = ifelse(doexit == pandemic, 0 , 1)) 

fit <- survfit(Surv(total_fup_days, flag) ~ exposed, data = base_cohort)

print(ggsurvplot(fit))


ggplot(base_cohort, aes(x = total_fup_days, fill = as.factor(exposed))) +
  geom_histogram(position = "dodge", bins = 50) +
  labs(x = "Follow-Up Days", y = "Count", fill = "Exposure Group")

fup <- base_cohort %>%
  group_by(exposed) %>%
  summarise(mean_fup = mean(total_fup_days, na.rm = TRUE),
            median_fup = median(total_fup_days, na.rm = TRUE),
            min_fup = min(dof_any_asthma_inc_dx, na.rm = TRUE),
            max_fup = max(dof_any_asthma_inc_dx, na.rm = TRUE))

print(fup)


colnames(base_cohort)



# Save the final dataset
saveRDS(base_cohort, file.path(path_datafiles_for_analysis, "cr_respiratory_dataforanalysis.rds"))


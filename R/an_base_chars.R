###############################################################################
#R script author: Kirsty Andresen
#Date: 23 OCtober 2024
#Description: Creates baseline tables for overall dataset and by cancer site
#######################################################################

###############################################################################
#1. Read in data
###############################################################################

#Create indicator variable called treatment in exposed == 1. the variable should be treatment = 1 if dof_chemo is between indexdate and doexit, and treatment = 2 if dof_radio is between indexdate and doexit, and 0 if neither
cr_finaldataforanalysis_respiratory <- cr_finaldataforanalysis_respiratory %>% 
  mutate(treatment = case_when(dof_chemo >= indexdate & dof_chemo <= doexit ~ 1,
                               dof_radio >= indexdate & dof_radio <= doexit ~ 2,
                               TRUE ~ 3))
systemic_treat <-cr_finaldataforanalysis_respiratory %>% filter(indexdate >= "2014-01-01") %>% filter(exposed == 1) %>%
  summarise(n = n(), n_chemo = sum(treatment == 1), percent = n_chemo/n*100)
radio_treat <-cr_finaldataforanalysis_respiratory %>% filter(indexdate >= "2012-01-01") %>% filter(exposed == 1) %>%
  summarise(n = n(), n_chemo = sum(treatment == 1), percent = n_chemo/n*100)

custom_freqs <- c(systemic_treat$n, radio_treat$n, nrow(cr_finaldataforanalysis_respiratory %>% filter(exposed == 1)) - systemic_treat$n - radio_treat$n)
custom_percents <- c(round(systemic_treat$percent, 1), round(radio_treat$percent, 1), round(100 - systemic_treat$percent - radio_treat$percent, 1))
custom_cumpercents <- cumsum(custom_percents)


# Variables of interest for the table
variables <- c("total_fup_yr", "age", "gender", "smokstatus", "bmi_cat", "eth5_hes", "imd5",
               "b_chronic_resp", "b_asthma", "b_copd", "b_bronch", "b_ild",
               "b_cvdgrouped", "b_diab_cat", "b_hypertension", "b_cld","b_multiple_sclerosis", "treatment")

# Categorical variables 
categorical <- c("gender", "smokstatus", "bmi_cat", "eth5_hes",  "imd5", 
                 "b_chronic_resp","b_asthma", "b_copd", "b_bronch", "b_ild",
                 "b_cvdgrouped", "b_diab_cat", "b_hypertension", "b_cld", "b_multiple_sclerosis", "treatment")

#Ensure the categorical variables are factors
cr_finaldataforanalysis_respiratory <- cr_finaldataforanalysis_respiratory %>%
mutate(across(all_of(categorical), as.factor),
       total_fup_yr = total_fup_days/365.25) %>%
  mutate(total_fup_yr= as.numeric(total_fup_yr)) 



###############################################################################
#2. OVERALL
###############################################################################
library(tableone)
# Create the baseline characteristics table
table1 <- CreateTableOne(vars = variables, 
                         data = cr_finaldataforanalysis_respiratory, 
                         factorVars = categorical,
                         includeNA = TRUE,
                         test = FALSE,
                         strata = "exposed")

#Manually change treatments
table1[["CatTable"]][["1"]][["treatment"]][, "freq"] <- custom_freqs
table1[["CatTable"]][["1"]][["treatment"]][, "percent"] <- custom_percents
table1[["CatTable"]][["1"]][["treatment"]][, "cum.percent"] <- custom_cumpercents

fup <- "total_fup_yr"

# Prepare the table for word
t1 <- print(table1, 
            quote = FALSE,
            
            nonnormal = fup)

write.csv(t1, paste0(path_results, "/table1_overall.csv"))

###############################################################################
#2. BY CANCER
###############################################################################

#Create empty list to store each table by cancer site
table_c_1 <- list()

# Create table for each cancersite

for (site in cancersites){
  
  cancer_data <- cr_finaldataforanalysis_respiratory %>%
    filter(cancer == site)
  
  
  systemic_treat <-cancer_data %>% filter(indexdate >= "2014-01-01") %>% filter(exposed == 1) %>%
    summarise(n = n(), n_chemo = sum(treatment == 1), percent = n_chemo/n*100)
  radio_treat <-cancer_data %>% filter(indexdate >= "2012-01-01") %>% filter(exposed == 1) %>%
    summarise(n = n(), n_chemo = sum(treatment == 1), percent = n_chemo/n*100)
  
  custom_freqs <- c(systemic_treat$n, radio_treat$n, nrow(cancer_data %>% filter(exposed == 1)) - systemic_treat$n - radio_treat$n)
  custom_percents <- c(round(systemic_treat$percent, 1), round(radio_treat$percent, 1), round(100 - systemic_treat$percent - radio_treat$percent, 1))
  custom_cumpercents <- cumsum(custom_percents)

  
  # Create the baseline characteristics table
  table_c_1[[site]] <- CreateTableOne(vars = variables, 
                                      data =  cancer_data, 
                                      factorVars = categorical,
                                      includeNA = TRUE,
                                      test = FALSE,
                                      strata = "exposed")
  #Manually change treatments
  table_c_1[[site]][["CatTable"]][["1"]][["treatment"]][, "freq"] <- custom_freqs
  table_c_1[[site]][["CatTable"]][["1"]][["treatment"]][, "percent"] <- custom_percents
  table_c_1[[site]][["CatTable"]][["1"]][["treatment"]][, "cum.percent"] <- custom_cumpercents
  

  # Prepare the table for word
  tc1 <- print(table_c_1[[site]],  quote = FALSE, nonnormal = fup)
  
  write.csv(tc1, paste0(path_results, "/table1_", site, ".csv"))
  
}

source(paste0(path_dofiles, "an_output_basechars.R"))

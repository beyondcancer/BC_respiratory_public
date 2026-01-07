rm(list=ls())
#### MASTER R SCRIPT ####
options(scipen=999) 
#set working directory

setwd("J:/EHR-Working/Krishnan/Kirsty/BC_respiratory_outcomes/paths") # Change to your own path
  
## READ PATHS AND DATA

source("pathsR.R")


# #########################

cancer_labels <- c(
  "Oral Cavity (C00-C06)",
  "Oesophageal (C15)",
  "Stomach (C16)",
  "Colorectal (C18-C21)",
  "Liver (C22)",
  "Pancreas (C25)",
  "Lung (C34)",
  "Melanoma (C43)",
  "Breast (C50)",
  "Cervix (C53)",
  "Uterus (C54-C55)",
  "Ovary (C56)",
  "Prostate (C61)",
  "Kidney (C64 - C65)",
  "Bladder (C67)",
  "Brain/CNS (C71-C72)",
  "Thyroid (C73)",
  "NHL (C82-C85)",
  "Myeloma (C90)",
  "Leukemia (C91-C95)"
)

#Set colours for future graphs
col1 <- "#073B4C"
col2 <- "#06D6A0"
col3  <-"#EF476F"
col4 <- "#FFD166"
col5 <- "#118AB2" 
col6 <- "#A2ACAB"
col7 <- "#621244"

# Create a named vector for mapping
cancer_label_map <- setNames(cancer_labels, cancersites)
outcome_labels <- c("Asthma", "COPD", "Bronchiectasis", "ILD")
outcome_label_map <- setNames(outcome_labels, outcomes_chronic)


objects_path <- ls()
###############################################################################
#1. DATA MANAGEMENT

##Read core_data from stata
core_data <- haven::read_dta(paste0(path_datafiles_for_analysis, "cr_coredataset_aurum.dta"))

##Create the base analysis dataset
 
source(paste0(path_dofiles, "cr_respiratory_dataforanalysisR_new.R"))
rm(list = (setdiff(ls(), c(objects_path, "objects_path", "core_data"))))

#read in dataset
cr_respiratory_dataforanalysis <- readRDS(paste0(path_datafiles_for_analysis, "cr_respiratory_dataforanalysis.rds"))

#Apply exclusions for first_ever DX analysis. 
######Produces flowchart####################
source(paste0(path_dofiles, "cr_apply_exclusions.R"))

rm(list = (setdiff(ls(), c(objects_path, "objects_path", "core_data"))))

###############################################################################
#2a. BASELINE CHARACTERISTICS 
###############################################################################


cr_finaldataforanalysis_respiratory <- readRDS(paste0(path_datafiles_for_analysis, "cr_finaldataforanalysis_respiratory.rds"))
objects_path <- ls()

##Baseline Characteristics

source(paste0(path_dofiles, "an_base_chars.R"))
rm(list = (setdiff(ls(), c(objects_path, "objects_path"))))

################################################################################ 
# #########FOR TESTING PURPOSES########
# outcomes_all <- c("asthma")
# outcomes_chronic <- c("asthma")
# outcome <- c("asthma")

# cancersites <- "col"
# site <- "col"

# ####################################
# # cr_finaldataforanalysis_respiratory <- readRDS(paste0(path_datafiles_for_analysis, "cr_finaldataforanalysis_respiratory.rds"))
# # #set seed
# # #number of unique setids
# # length(unique(cr_finaldataforanalysis_respiratory$setid))
# # set.seed(123)
# # random_setids <- sample(unique(cr_finaldataforanalysis_respiratory$setid), size = 100000, replace = FALSE)
# # sample_data <- cr_finaldataforanalysis_respiratory %>% dplyr::filter(setid %in% random_setids)
# # saveRDS(sample_data, paste0(path_datafiles_for_analysis, "sample_data.rds"))
#  cr_finaldataforanalysis_respiratory <- readRDS(paste0(path_datafiles_for_analysis, "sample_data.rds"))
# ###############################################################################################################3

###############################################################################
#2a. CREATE - FIRST RECORD ANALYSIS
###############################################################################
# 
# outcomes_all <- "asthma"
# outcomes_chronic <- "asthma"
# outcome <- "asthma"
# cancersites <- cancersites[1:10]


cr_finaldataforanalysis_respiratory <- readRDS(paste0(path_datafiles_for_analysis, "cr_finaldataforanalysis_respiratory.rds"))
objects_path <- ls()

#Create data for first record analysis
source(paste0(path_dofiles, "cr_outcome_an_first_dx.R"))

rm(list = (setdiff(ls(), c(objects_path, "objects_path"))))


#############################################################################
#2b CREATE EPISODES
#############################################################################

cr_finaldataforanalysis_respiratory <- readRDS(paste0(path_datafiles_for_analysis, "cr_finaldataforanalysis_respiratory.rds"))
objects_path <- ls()

source(paste0(path_dofiles, "cr_episodes_respiratory_new.R"))
rm(list = (setdiff(ls(), c(objects_path, "objects_path"))))
###############################################################################
#ANALYSIS - HRs
###############################################################################

## Analysis for first_ever Dx analysis
source(paste0(path_dofiles, "an_cox_models.R"))
rm(list = (setdiff(ls(), c(objects_path, "objects_path"))))

#Also run for other analysis 

## Analysis for exacerbation Dx analysis
source(paste0(path_dofiles, "an_cox_models_episodes.R"))
rm(list = (setdiff(ls(), c(objects_path, "objects_path"))))

source(paste0(path_dofiles, "an_cox_models_strat.R"))
rm(list = (setdiff(ls(), c(objects_path, "objects_path"))))

## Analysis for exacerbation Dx analysis
source(paste0(path_dofiles, "an_cox_models_episodes_strat.R"))
rm(list = (setdiff(ls(), c(objects_path, "objects_path"))))

###############################################################################
#GRAPHS - FOREST PLOTS
###############################################################################

results_crude_adj_t <- readRDS(paste0(path_datafiles_for_analysis, "an_cox_models_first_dx_crude_adj.rds"))
results_strat <- readRDS(paste0(path_datafiles_for_analysis, "an_cox_models_first_strat.rds"))

results <- bind_rows(results_crude_adj_t, results_strat)
results$cancer <- factor(results$cancer, levels = cancersites)
saveRDS(results, paste0(path_results,"/results_firstdx.rds"))
write.csv(results, paste0(path_results, "/results_firstdx.csv"))
###############################################################################

results_exc_crude_adj_t <- readRDS(paste0(path_datafiles_for_analysis, "results_exc.rds"))
results_exc_strat <- readRDS(paste0(path_datafiles_for_analysis, "an_cox_models_episodes_strat.rds"))

results <- bind_rows(results_exc_crude_adj_t, results_exc_strat)
results$cancer <- factor(results$cancer, levels = cancersites)
saveRDS(results, paste0(path_results,"/results_exc.rds"))
write.csv(results, paste0(path_results, "/results_exc.csv"))
## Graphical forest plots for first_ever Dx analysis

source(paste0(path_dofiles, "an_forest_plot.R"))
rm(list = (setdiff(ls(), c(objects_path, "objects_path"))))

###############################################################################
#ANALYSIS - INCIDENCE
###############################################################################

# Adjusted incidence graphs for first ever dx (Direct standardisation)
source(paste0(path_dofiles, "an_incidence_rp.R"))
rm(list = (setdiff(ls(), c(objects_path, "objects_path"))))

# Adjusted incidence graphs
source(paste0(path_dofiles, "an_mcc_episodes.R"))
rm(list = (setdiff(ls(), c(objects_path, "objects_path"))))

###############################################################################
#GRAPHS - RISK DIFFERENCE
###############################################################################

# # Incidence tables (Crude) for first ever Dx
# source(paste0(path_dofiles, "an_describe_incidence.R"))
# rm(list = (setdiff(ls(), c(objects_path, "objects_path"))))

## Calculate risk difference
source(paste0(path_dofiles, "an_risk_dif.R"))
rm(list = (setdiff(ls(), c(objects_path, "objects_path"))))

###############################################################################
#SENSITIVITY ANALYSES
###############################################################################

# outcomes_chronic <- c("asthma","copd")
# cancersites <- top9cancers

# Sensitivity analysis for first_ever Dx analysis
cr_finaldataforanalysis_respiratory <- readRDS(paste0(path_datafiles_for_analysis, "cr_finaldataforanalysis_respiratory.rds"))
objects_path <- ls()
source(paste0(path_dofiles, "cr_outcome_an_first_dx_allprev.R"))
source(paste0(path_dofiles, "cr_outcome_an_first_dx_b5.R"))
source(paste0(path_dofiles, "cr_outcome_an_first_dx_smok.R"))

       ##########
source(paste0(path_dofiles, "an_cox_models_allprev.R"))
source(paste0(path_dofiles, "an_cox_models_b5.R"))
source(paste0(path_dofiles, "an_cox_models_up2.R"))
source(paste0(path_dofiles, "an_cox_models_smok.R"))
#########
source(paste0(path_dofiles, "an_forest_plot_sens.R"))

#####################
#EXACERBATIONS
####################

cr_finaldataforanalysis_respiratory <- readRDS(paste0(path_datafiles_for_analysis, "cr_finaldataforanalysis_respiratory.rds"))
objects_path <- ls()
source(paste0(path_dofiles, "cr_respiratory_new_30.R"))
source(paste0(path_dofiles, "an_cox_models_episodes_30.R"))

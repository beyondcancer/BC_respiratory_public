##########
# MASTER: Response to Reviewers Analysis Code
##########
# This script sources all analysis and figure scripts added or updated in
# response to reviewer comments. Run scripts in order as some depend on
# outputs from earlier steps.
#
# Sections:
#   1. Smoking effect modification analysis (Figures A13-A20)
#   2. Post-hoc treatment analysis - radiotherapy & HSCT (Figures A27-A28)
#   3. Risk difference analysis & proportional hazards check (Figure 2)
#   4. Follow-up time graph (Figure A4)
#
# Author: Kirsty Andresen
# Date: May 2026
##########

##########
# 1. Smoking effect modification analysis (Figures A13-A20)
# Adds smoking interaction term to stratified Cox models for both
# new-onset and exacerbation outcomes, then generates forest plots.
##########

# New-onset outcomes
source(file.path(path_dofiles, "review", "an_cox_models_strat_rev.R"))

# Exacerbation outcomes
source(file.path(path_dofiles, "review", "an_cox_models_strat_rev.R"))

# Forest plots - Figures A13-A20
source(file.path(path_dofiles, "review", "gr_forest_plot_strat_rev.R"))

##########
# 2. Post-hoc treatment analysis (new-onset only - Figures A27-A28)
# Recategorises exposure by receipt of chest radiotherapy (breast, lung,
# NHL, oesophageal; 2012 onwards) and HSCT (NHL only) to examine whether
# treatment drives observed associations.
##########

# Data preparation: create treatment-recategorised exposure variable
source(file.path(path_dofiles, "review", "cr_treatmentanalysis.R"))

# Cox models with recategorised exposure
source(file.path(path_dofiles, "review", "cr_an_cox_models_trt.R"))

# Forest plot - radiotherapy analysis (Figure A27)
source(file.path(path_dofiles, "review", "gr_forest_plot_trt_radio.R"))

# Forest plot - HSCT analysis (Figure A28)
source(file.path(path_dofiles, "review", "gr_forest_plot_trt_hsct.R"))

# Events and follow-up time table for treatment analysis
source(file.path(path_dofiles, "review", "gr_trt_events_fup.R"))

##########
# 3. Risk difference analysis & proportional hazards check (Figure 2 - new Table 1)
# Updates main Cox models and derives absolute risk differences.
# Also includes proportional hazards assumption checks (cox.zph).
##########

# Updated Cox models
source(file.path(path_dofiles, "review", "an_cox_models_review.R"))

# Risk difference calculations
source(file.path(path_dofiles, "review", "an_risk_dif_review.R"))

##########
# 4. Follow-up time graph (Figure A4)
# Plots median follow-up time by cancer type.
##########

source(file.path(path_dofiles, "review", "gr_median_fup_by_cancer.R"))
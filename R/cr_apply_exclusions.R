
###############################################################################
#R script author: Kirsty Andresen
#Date: 11 September 2024
#Description: Creates flowchart with patient exclusions for BC_respiratory cancers
###############################################################################


#Apply general exclusions
flowchart <- cr_respiratory_dataforanalysis %>% 
  mutate(excluded = case_when(age < 18 ~ 1,
                              time_b4_fup < 365 ~ 2,
                              doexit <= indexdate ~ 3,
                              doexit <= dostartcprdfup ~ 4,
                              total_fup_days <0 ~ 5) 
                              )%>%
  mutate(excluded= factor(excluded,
                          levels = 1:5,
                          labels=c("Under 18 years of age", 
                                   "Follow-up time <12 months before study entry",
                                   "End of follow up before indexdate",
                                   "End of follow-up before CPRD start date",
                                   "Total follow-up time <0 days")

                          ))%>%
  mutate(exposed_pre_mtch = case_when(exposed == 1 & is.na(excluded) ~  "Cancer diagnosis",
                                 exposed == 0 & is.na(excluded) ~ "No history of cancer"))

print(paste0( "Number of patients excluded: ", nrow(flowchart[!is.na(flowchart$excluded), ])))
print(table(flowchart$excluded))

###############################################################################
#2. Remove cases without controls and controls without cases 
##  after exclusions are applied and include into flowchart dataframe
###############################################################################

flowchart2 <- flowchart %>% 
  filter(is.na(excluded)) %>%
  arrange(setid, exposed) %>%
  group_by(setid) %>%
  mutate(anyunexposed = min(exposed)) %>%
  arrange(setid, desc(exposed)) %>%
  mutate(anyexposed = max(exposed)) %>%
  ungroup()
 
  print("Number of patients with no assigned controls")
  print(table(flowchart2$anyunexposed))
  print("Number of patients with no assigned cases")
  print(table(flowchart2$anyexposed))
  
         
flowchart <- flowchart %>% 
  left_join(flowchart2[c("e_patid", "setid", "anyunexposed", "anyexposed")], by=c("setid", "e_patid"))


flowchart <- flowchart %>%
  mutate(excluded_mtch = case_when(anyunexposed == 1 ~ 1, 
                              anyexposed == 0 ~ 2
                              ))%>%
  mutate(excluded_mtch= factor(excluded_mtch,
                          levels = 1:2,
                          labels=c("Cases with no assigned controls",
                                   "Controls with no assigned cases")
                          )) %>%
  mutate(exposed_an = case_when(exposed == 1 & is.na(excluded) & is.na(excluded_mtch) ~  "Cancer diagnosis",
                              exposed == 0  & is.na(excluded) & is.na(excluded_mtch) ~ "No history of cancer"))

#Count  those with exposed_an = 1 "Cancer diagnosis" in each category of cancer

cancer_summary <- flowchart %>%
  filter(exposed_an == "Cancer diagnosis") %>%
  group_by(cancer) %>%
  summarise(n_cancer = n()) %>%
  arrange(desc(n_cancer)) %>%
  mutate(label_cancer = paste0(cancer, ": ", n_cancer))
  
cancer_summary_label <- paste(cancer_summary$label, collapse = "\n")

control_summary <- flowchart %>%
  filter(exposed_an == "No history of cancer") %>%
  group_by(cancer) %>%
  summarise(n_controls = n()) %>%
  arrange(desc(n_controls)) %>%
  mutate(label_control = paste0(cancer, " controls : ", n_controls))

control_summary_label <- paste(control_summary$label_control, collapse = "\n")


colnames(flowchart)
table(flowchart$anyunexposed)
table(flowchart$excluded)


table(flowchart$exposed_pre_mtch)
table(flowchart$exposed_an)

#drop patients where both excluded and excluded_mtch are not missing
count(flowchart, exposed_an)

cr_finaldataforanalysis_respiratory <- flowchart %>% 
  filter(!is.na(exposed_an)) %>%
  dplyr::select(-excluded, -exposed_pre_mtch, -excluded_mtch, -anyunexposed, -anyexposed, -exposed_an)

print(colnames(cr_finaldataforanalysis_respiratory))
#save data set with exclusions as rds file using path_datafiles_for_analysis


rm(flowchart2)
rm(cr_respiratory_dataforanalysis)
gc()
saveRDS(cr_finaldataforanalysis_respiratory, file = paste0(path_datafiles_for_analysis, "cr_finaldataforanalysis_respiratory.rds"))

###############################################################################
#3. Create flowchart
###############################################################################


flow <- consort_plot(
  data = flowchart,
  orders = c(
    exposed = "Patients in core matched dataset",
    excluded = "Excluded",
    exposed_pre_mtch = "Eligible Patients",
    excluded_mtch = "",
    exposed_an = ""
  ),
  side_box = c("excluded"),
  allocation = "exposed_pre_mtch"
)

# Display flowchart
plot(flow)

ggsave(paste0(path_results,"/flowchart.png"), plot = build_grid(flow), width = 8, height = 8, units = "in", dpi = 300)

cancer_labels <- c(
  "Oral Cavity (C00-C06)",
  "Oesophageal (C15)",
  "Stomach (C16)",
  "Colorectal (C18-C20)",
  "Liver (C22)",
  "Pancreas (C25)",
  "Lung (C34)",
  "Melanoma (C43)",
  "Breast (C50)",
  "Cervix (C53)",
  "Uterus (C54-C55)",
  "Ovary (C56)",
  "Prostate (C61)",
  "Kidney (C64)",
  "Bladder (C67)",
  "Brain/CNS (C71-C72)",
  "Thyroid (C73)",
  "NHL (C82-C85)",
  "Myeloma (C90)",
  "Leukemia (C91-C95)"
)

cancer_label_map <- setNames(cancer_labels, cancersites)


#bind summary tables
summary_tables <- cancer_summary %>% 
  left_join(control_summary, by = "cancer") %>%
  mutate(cancer= cancer_label_map[cancer])

#save summary tables as csv
write.csv(summary_tables, paste0(path_results, "/finaldataforanalysisnumbers_percancer.csv"), row.names = FALSE)


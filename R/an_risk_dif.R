###############################################################################
#R script author: Kirsty Andresen
#Date: 20 January 2025
#Description: Create risk differences by cancer
################################################################################

#Load in incidence table
incidence_table <- read.csv(paste0(path_results, "/first_dx_analysis_numbers_.csv"))%>%
  dplyr::select(cancer, outcome,n_events_exp, n_events_unexp, IR_exp, IR_unexp, t_fup_exp, t_fup_unexp) %>%
  mutate(events_combined = paste("[", n_events_exp,"/", n_events_unexp,"]", sep="")) %>%
  mutate(IR_combined = paste("[", sprintf("%.1f", IR_exp), "/", sprintf("%.1f", IR_unexp), "]", sep = ""))
  

  
#load in HR table (adjusted HR)
results <- readRDS(paste0(path_results, "/results_firstdx.rds")) %>%
           dplyr::filter(model_type == "Adjusted")



# Define the corresponding labels for cancer sites
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



#reshape results to columns 
results <- results %>%
  tidyr::pivot_wider(names_from = model_type, values_from = hr)


#join incidence_table and results by outcome and site
m <- results %>%
  dplyr::left_join(incidence_table, by = c("outcome" = "outcome", "cancer" = "cancer"))

# Create a named vector for mapping

# Use mutate to create a new column with the corresponding labels
m <- m %>%
  mutate(cancer_label = cancer_label_map[cancer],
         outcome_label = outcome_label_map[outcome]) %>%
  mutate(cancer_label = paste(cancer_label, events_combined, sep = "\n")) %>%
  mutate(outcome_label = paste(outcome_label, IR_combined, sep = " "))

#Combine the two tables
RD <- m %>% dplyr::mutate(rate_unexp_t= IR_unexp * (1/Adjusted))

RD <- RD %>% mutate(IR_expSE = IR_exp/sqrt(n_events_exp),
                    IR_unexpSE = IR_unexp/sqrt(n_events_unexp),
                    IR_expCIupper = IR_exp + 1.96 * IR_expSE,
                    IR_expCIlower = IR_exp - 1.96 * IR_expSE,
                    IR_unexpCIupper = IR_unexp + 1.96 * IR_unexpSE,
                    IR_unexpCIlower = IR_unexp - 1.96 * IR_unexpSE,
                    RD = (IR_exp - rate_unexp_t) ,
                    RD_uadj = IR_exp - IR_unexp ,
                    NNH = 1/RD_uadj, 
                    SE_RD = sqrt((IR_exp/n_events_exp) + (IR_unexp/n_events_unexp)), # Standard Error
                    RD_uadj_upper = RD_uadj + 1.96 * SE_RD,
                    RD_uadj_lower = RD_uadj + 1.96 * SE_RD
                    ) # 95% CI Upper Bound

RD <- RD %>% dplyr::select(cancer, outcome, events_combined, RD, RD_uadj, RD_uadj_upper, RD_uadj_lower, NNH, outcome_label)

write.csv(RD, paste0(path_results, "/RD.csv"))

# Define color palette for outcomes
outcome_colors <- c("#0D5257", "#00BF6F", "#FFB81C", "#A2ACAB")

grob_list <- list()

# # Create plots for each site
# for (site in cancersites) {
#   RD_can <- RD %>%
#     dplyr::filter(cancer == site) %>%
#     mutate(outcome = as.factor(outcome))  # Ensure 'outcome' is a factor
#   
#   RD_plot <- RD_can %>%
#     ggplot(aes(x = outcome_label, y = RD_uadj, fill = outcome)) +  # Map 'outcome' to fill
#     geom_bar(stat = "identity", position = "stack") +
#     # geom_text(aes(label = events_combined), size = 2, hjust = 0.1) +
#     # geom_text(aes(label = round(NNH, 2), hjust = -0.2),  # Add NNH labels
#     #           position = position_stack(vjust = 0.5),
#     #           size = 3) +  # Adjust label size
#     scale_fill_manual(values = outcome_colors) +  # Use custom colors for outcomes
#     coord_flip() +
#     theme_minimal() +
#     labs(title = cancer_label_map[site],
#          x = "" , y = "") +
#     scale_y_continuous(limits = c(-3, 40)) +
#     theme_minimal(base_size = 20)
#     theme(legend.position = "none",
#           plot.title = element_text(size=9)) 
# 
#   grob_list[[site]] <- grid.grabExpr(print(RD_plot))
# }
# 
# # Combine plots
# grid_plot <- do.call(gridExtra::arrangeGrob, c(grob_list, ncol = 4))
# 
# rd_gph <- gridExtra::grid.arrange(
#   grid_plot,
#   top = grid::textGrob("Crude rate difference per 1000 person-yrs",
#                        gp = gpar(fontsize = 16, fontface = "bold"))
# )
# 
# # Save the combined plot
# ggsave(filename = paste0(path_results, "/RD_plot_uadj.png"),
#        plot = rd_gph,
#        width = 24,  height = 22, units = "cm")


for (site in cancersites) {
  RD_can <- RD %>%
    dplyr::filter(cancer == site) %>%
    mutate(outcome = as.factor(outcome))  # Ensure 'outcome' is a factor
  
  RD_plot <- RD_can %>%
    ggplot(aes(x = outcome_label, y = RD, fill = outcome)) +  # Map 'outcome' to fill
    geom_bar(stat = "identity", position = "stack") +
    # geom_text(aes(label = round(NNH, 2), hjust = -0.2),  # Add NNH labels
    #           position = position_stack(vjust = 0.5),
    #           size = 3) +  # Adjust label size
    scale_fill_manual(values = outcome_colors) +  # Use custom colors for outcomes
    coord_flip() +
    theme_minimal(base_size = 14) +
    labs(title = cancer_label_map[site],
         x = "", y = "") +
    scale_y_continuous(limits = c(-3, 40)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 14, face = "bold", color= col1)) # Set consistent x-axis range
  
  grob_list[[site]] <- grid.grabExpr(print(RD_plot))
}

# Combine plots
grid_plot <- do.call(gridExtra::arrangeGrob, c(grob_list, ncol = 4))

rd_gph <- gridExtra::grid.arrange(
  grid_plot,
  top = grid::textGrob("Adjusted Rate Difference per 1000 person-years",
                       gp = gpar(fontsize = 16, fontface = "bold"))
)

# Save the combined plot
ggsave(
  filename = paste0(path_results,"/RD_adj.tiff"),
  plot = rd_gph,
  width = 16,         # fits within A4 with margins
  height = 11,       # portrait layout
  units = "in",
  dpi = 300,
  bg = "white",
  device = "tiff",
  compression = "lzw"  # optional: smaller file size, high quality
)

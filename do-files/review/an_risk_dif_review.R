###############################################################################
#R script author: Kirsty Andresen
#Date: 20 January 2025
#Description: Create risk differences by cancer
################################################################################
# 
# #Load in incidence table
# incidence_table <- read.csv(paste0(path_results, "/first_dx_analysis_numbers_.csv"))%>%
#   dplyr::select(cancer, outcome,n_events_exp, n_events_unexp, IR_exp, IR_unexp, t_fup_exp, t_fup_unexp) %>%
#   mutate(events_combined = paste("[", n_events_exp,"/", n_events_unexp,"]", sep="")) %>%
#   mutate(IR_combined = paste("[", sprintf("%.1f", IR_exp), "/", sprintf("%.1f", IR_unexp), "]", sep = ""))
#   

#   
# #load in HR table (adjusted HR)
# results <- readRDS(paste0(path_results, "/results_firstdx.rds")) %>%
#            dplyr::filter(model_type == "Adjusted")

hr <- readRDS(paste0(path_datafiles_for_analysis, "rev/an_cox_models_first_dx_crude_adj.rds"))

inc_rate <- read.csv(paste0(path_datafiles_for_analysis, "rev/prop_hazards_inc_rate.csv"))


#Join hr column from results to review vars

RD <- inc_rate %>% left_join(hr, by = c("cancer", "outcome", "model_type"))

RD <- RD %>%
  mutate(
    # Incidence rate per 1,000 person-years
    IR_exp = (events_exp_adj / total_fup_exp_adj) * 1000,
    IR_unexp = (events_unexp_adj / total_fup_unexp_adj) * 1000,
    
    # 95% CI for exposed rate using exact Poisson
    IR_exp_lower = (qgamma(0.025, shape = events_exp_adj, rate = 1) / total_fup_exp_adj) * 1000,
    IR_exp_upper = (qgamma(0.975, shape = events_exp_adj + 1, rate = 1) / total_fup_exp_adj) * 1000,
    
    # 95% CI for unexposed rate using exact Poisson
    IR_unexp_lower = (qgamma(0.025, shape = events_unexp_adj, rate = 1) / total_fup_unexp_adj) * 1000,
    IR_unexp_upper = (qgamma(0.975, shape = events_unexp_adj + 1, rate = 1) / total_fup_unexp_adj) * 1000,
    
    # Risk difference
    IR_counterfactual = IR_exp * (1 / hr),
    RD = IR_exp - IR_counterfactual
  ) %>%
  mutate(
    # Variance of the log HR (from HR and its CIs)
    log_hr_se = (log(ci_upper) - log(ci_lower)) / (2 * 1.96),
    
    # Delta method variance of RD
    RD_var = (IR_exp * (1/hr))^2 * log_hr_se^2 +   # uncertainty from HR
      (1 - 1/hr)^2 * (IR_exp / events_exp_adj), # uncertainty from IR
    
    RD_se = sqrt(RD_var),
    RD_lower = RD - 1.96 * RD_se,
    RD_upper = RD + 1.96 * RD_se) %>% 
  mutate(cancer_label = cancer_label_map[cancer],
         outcome_label = outcome_label_map[outcome]) 

RD <- RD %>%
  mutate(events_combined = paste("[", events_exp_adj,"/", events_unexp_adj,"]", sep=""),
         IR_combined = paste("[", sprintf("%.1f", IR_exp), "/", sprintf("%.1f", IR_unexp), "]", sep = "")) %>%
  mutate(cancer_label = paste(cancer_label, events_combined, sep = "\n")) %>%
  mutate(outcome_label = paste(outcome_label, IR_combined, sep = " ")) 

#Extract the IR adjusted only

RD_overall <- RD %>% filter(model_type %in% "Adjusted") %>%
  group_by(cancer, outcome) %>%
  summarise(
    cancer_label = first(cancer_label),
    outcome_label = first(outcome_label)
  )

RD_time <- RD %>% filter(model_type %in% c("0-1 years", "1-3 years", "3-5 years", "5+ years")) %>%
  group_by(cancer, outcome) %>%
  summarise(
    RD = mean(RD, na.rm = TRUE),
    RD_lower = mean(RD_lower, na.rm = TRUE),
    RD_upper = mean(RD_upper, na.rm = TRUE),
  ) %>% left_join(RD_overall, by = c("cancer", "outcome")) 



write.csv(RD_time, paste0(path_results, "/rev/RD.csv"))

RD <- RD_time

col1<- "#000"
# Define color palette for outcomes
outcome_colors <- c("#0D5257", "#621244", "#FFB81C", "#A2ACAB")
#outcome_colors <- c("#000", "#000", "#000", "#000")

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
  
  header_df <- RD_can %>%
    filter(toupper(as.character(outcome)) == "ILD") %>%
    distinct(outcome_label) %>%
    mutate(y = -3)
  
 
  RD_plot <- RD_can %>%
    ggplot(aes(x = outcome_label, y = RD, fill = outcome)) +  # Map 'outcome' to fill
    geom_bar(stat = "identity", position = "stack") +
    # geom_text(aes(label = round(NNH, 2), hjust = -0.2),  # Add NNH labels
    #           position = position_stack(vjust = 0.5),
    #           size = 3) +  # Adjust label size
    scale_fill_manual(values = outcome_colors) +  # Use custom colors for outcomes
    coord_flip(clip = "off") +
    geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "red") +
    geom_text(
      data = header_df,
      aes(x = outcome_label, y = y),
      label = "Outcome IR [Survivor/Control]",
      size = 3.5,
      hjust = 1,      # ← align with left edge of ILD label
      vjust = -1.3,   # ← lift it higher into header space
      inherit.aes = FALSE
    ) +
    theme_minimal(base_size = 14) +
    labs(title = cancer_label_map[site],
         x = "", y = "") +
    scale_y_continuous(limits = c(-3, 40)) +
    theme(legend.position = "none",
          plot.title = element_text(size = 14, face = "bold", color= col1),
          #axis.title.y = element_text(size = 3.5, face ="italic"),
          plot.margin = margin(t = 3, r = 5, b = 5, l = 0)) # Set consistent x-axis range
  
  grob_list[[site]] <- grid.grabExpr(print(RD_plot))
}

# Combine plots
# grid_plot <- do.call(gridExtra::arrangeGrob, c(grob_list, ncol = 4, ))

grid_plot <- gridExtra::arrangeGrob(
  grobs  = grob_list,
  ncol   = 4,
  widths = grid::unit(rep(1, 4), "null")
)

rd_gph <- gridExtra::grid.arrange(
  grid_plot,
  top = grid::textGrob("Adjusted Rate Difference per 1000 person-years",
                       gp = gpar(fontsize = 16, fontface = "bold"))
)

pdf(paste0(path_results,"/rev/RD_adj.pdf"), 
    width = 16, 
    height = 11)
grid::grid.draw(rd_gph)
dev.off()

# # Save the combined plot
# ggsave(
#   filename = paste0(path_results,"/RD_adj.tiff"),
#   plot = rd_gph,
#   width = 16,         # fits within A4 with margins
#   height = 11,       # portrait layout
#   units = "in",
#   dpi = 300,
#   bg = "white",
#   device = "tiff",
#   compression = "lzw"  # optional: smaller file size, high quality
# )

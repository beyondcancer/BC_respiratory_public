
#forest plot faceted by outcome


library(ggh4x)
library(forcats)
library(stringr)
# # 
# col1 <- "#0D5257"
# col2 <- "#06D6A0"
# col3  <-"#EF476F"
# col4 <- "#FFD166"
# col5 <- "#118AB2"
# col6 <- "#A2ACAB"
# col7 <- "#621244"

col1 <- "#000000"
col2 <- "#3A3A3A"
col3 <- "#5E5E5E"
col4 <- "#828282"
col5 <- "#A6A6A6"
col6 <- "#C8C8C8"
col7 <- "#E0E0E0"

results <- readRDS(paste0(
  path_datafiles_for_analysis,
  "an_cox_models_first_dx_trt_nhl_hsct.rds"
))


#create variable study with the first word of model type
results <- results %>% mutate(study = word(model_type, 1))
results <- results %>%
  mutate(model_type = case_when(
    model_type == "hsct Treatment" ~ "HSC therapy",
    model_type == "hsct No Treatment" ~ "No HSC therapy",
    TRUE ~ model_type
  ))
# 
# results_radio <-  results %>% filter(study ==  "Radio" & (cancer == "bre" | cancer == "oes"| cancer == "nhl"| cancer == "lun"))
# 
# name <- "radio"
# ft <- "radio"
#

name <- "hsct"
ft <- "hsct"

#Ensure ordering is correct
results <- results %>%
  mutate(
    model_type = fct_inorder(model_type),
    model_type = fct_rev(model_type))

# #Create all combinations of cancer, model_type, outcome to ensure missing combinations are shown as NA
# all_combinations <- expand.grid(
#   cancer = unique(results$cancer),
#   model_type = unique(results$model_type),
#   outcome = unique(results$outcome))
# 
# # Join with existing data to fill in missing combinations with NA
# results <- results %>%
#   right_join(all_combinations, by = c("cancer", "model_type", "outcome")) %>%
#   arrange(outcome, cancer, model_type)
# 
#Format the numbers and create labels
results <- results %>%
  mutate(
    cancer = gsub("\\s*\\(", "\n(", cancer),
    hr = as.numeric(sprintf("%.2f", hr)),
    ci_lower = as.numeric(sprintf("%.2f", ci_lower)),
    ci_upper = as.numeric(sprintf("%.2f", ci_upper)),
    cancer_model = paste0(cancer, " (", model_type, ")"),
    outcome = outcome_label_map[outcome],
    cancer = cancer_label_map[as.character(cancer)],
    hr_label = paste0("    ",hr,"(", ci_lower,"-", ci_upper,")")
  ) 

#Determine y positions for each row
results <- results %>%
  mutate(
    model_type = fct_inorder(model_type), 
    cancer = fct_inorder(cancer),
    cancer = fct_rev(cancer)
  ) %>%
  arrange(outcome, cancer, model_type) %>%
  group_by(outcome) %>%
  mutate(label_y = row_number()) %>%
  ungroup()

#Set out plot parameters
xmin <- 0.25 #CI limits
xmax <- 19.1 #CI limits
x_label_col <- 19 # x position for HR text labels

model_types <- unique(results$model_type)

#Main data for the plot (with clipped estimates)
plot_data <- results %>%
  mutate(
    ci_lower_capped = pmax(ci_lower, xmin),
    ci_upper_capped = pmin(ci_upper, xmax),
    clipped_left = ci_lower < xmin,
    clipped_right = ci_upper > xmax,
    label_col = hr_label,
  )


#Data for the header row
hr_header_df <- plot_data %>%
  group_by(outcome) %>%
  summarise(
    x = x_label_col,
    y = max(label_y) + 0.5,   # adjust spacing here
    .groups = "drop"
  )

## Set the colour map

custom_colors <- rev(c(col1, col2, col3))
color_map <- setNames(custom_colors, unique(results$model_type))

# Get ymin/ymax for each cancer based on label_y order
# Zebra striping directly from label_y
# Zebra striping every two rows using label_y
zebra <- results %>%
  group_by(cancer) %>%
  summarise(
    ymin = min(label_y) - 0.32,
    ymax = max(label_y) + 0.32,
    .groups = "drop"
  ) %>%
  arrange(ymin) %>%
  mutate(
    stripe_index = row_number(),
    fill = if_else(stripe_index %% 2 == 0, "#EAEAEA", "#FFFFFF")
  )

# Calculate midpoint y between the two model types per cancer per outcome
pval_df <- plot_data %>%
  group_by(cancer, outcome) %>%
  summarise(
    y_mid  = mean(label_y),
    lrtest = first(lrtest),
    .groups = "drop"
  ) %>%
  mutate(pval_label = case_when(
    is.na(lrtest)  ~ "NA",
    lrtest < 0.001 ~ "LRT pval <0.001",
    TRUE           ~ paste0("LRT pval=", round(lrtest, 3))
  ))

##

y_labels <- plot_data %>%
  group_by(cancer) %>%
  summarise(label_y_mid = mean(label_y), .groups = "drop")


#PLOT

p <- ggplot() +
  
  #Create zebra striping
  geom_rect(
    data = zebra,
    aes(ymin = ymin, ymax = ymax, fill = fill),
    xmin = -Inf, xmax = Inf,
    inherit.aes = FALSE
  ) +
  
  #Add x ticks
  geom_vline(
    xintercept = c(0.25, 0.5, 2, 5, 10),
    linewidth = 0.4,
    linetype = "dashed",
    colour = "grey80"
  ) +
  
  #Add null line
  geom_vline(xintercept = 1,
             linetype = "dashed",
             linewidth = 0.5,
             color = "red"
  ) +
  
  scale_fill_identity(guide = "none") +
  
  #Add 95% error bars for time stratifiied 
  geom_errorbarh( data = subset(plot_data, !is.na(hr)),
                  aes(xmin = ci_lower_capped, xmax = ci_upper_capped, y = label_y, color = model_type),
                  height = 0.2, linewidth = 0.8
  ) +
  #Add 95% error bars for adjusted → capped lines for those out of range
  geom_point(
    data = subset(plot_data, !is.na(hr) & clipped_left),
    aes(x = xmin, y = label_y, color = model_type),
    shape = "<", size = 8,
    show.legend = FALSE
  ) +
  geom_point(
    data = subset(plot_data, !is.na(hr) & clipped_right),
    aes(x = xmax, y = label_y, color = model_type),
    shape = ">", size = 8,
    show.legend = FALSE
  ) +
  # Non-adjusted → normal points
  geom_point(
    data = subset(plot_data, !is.na(hr) & !clipped_left & !clipped_right),
    aes(x = hr, y = label_y, color = model_type, shape = model_type), size = 3) +  
  
  # Customize the ticks x axis (log scale) - capped at 10
  scale_x_log10(
    limits = c(0.25, 20),
    breaks = c(0.2, 0.5, 1, 2, 5, 10),
    labels = c("0.2", "0.5", "1", "2", "5", "10")
  ) +
  
  #Customise the width of the x axis (need to leave space for labels)
  coord_cartesian(xlim = c(0.25, 200)) +
  
  #Customise y axis and add labels
  scale_y_continuous(
    breaks = y_labels$label_y_mid,
    labels = y_labels$cancer,
    expand = expansion(mult = c(0.005, 0.02))
  ) + 
  geom_text(
    data = pval_df,
    aes(x = x_label_col + 0.7,   # adjust x to sit right of HR labels
        y = y_mid,
        label = pval_label),
    hjust = 0,
    size  = 2,
    color = col1,
    inherit.aes = FALSE
  ) +
  #Customise colors of graph
  scale_color_manual(
    values = color_map,
    name   = "Treatment",
    labels = c("No HSC therapy" = "No HSC therapy" , "HSC therapy" = "HSC therapy"),
    guide  = guide_legend(reverse = TRUE)
  ) +
  scale_shape_manual(
    values = c("No HSC therapy" = 16, "HSC therapy"  = 15),
    name   = "Treatment",
    labels = c("No HSC therapy" = "No HSC therapy" , "HSC therapy" = "HSC therapy"),
    guide  = guide_legend(reverse = TRUE)
  ) +
  # Facet by outcome
  facet_grid(
    cols = vars(outcome),
    scales = "free_y"
  ) +
  
  # Add titles to x axis and to legend)
  labs(
    x = "Hazard Ratio (95% CI)",
    y = NULL,
    fill = "Category",
    color = "Category"
  ) +
  
  # Add adjusted HR header
  geom_text(
    data = hr_header_df,
    aes(x = x, y = y),
    label = "adjHR (95% CI)",
    hjust = 0.1,
    vjust = 0.8,
    fontface = "bold",
    size = 2.2,
    color = col1,
    inherit.aes = FALSE
  ) +
  
  # Add HR labels for adjusted rows
  geom_text(
    data = plot_data,
    aes(x = x_label_col, y = label_y, label = label_col),
    color = col1,
    hjust = 0,
    vjust = 0.9,
    size =2,
    lineheight = 0.90,
    inherit.aes = FALSE
  ) +
  
  # Customize overall theme
  theme(
    plot.title.position = "plot",
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5, color = col1),
    legend.position = "bottom",
    legend.justification = "center",
    legend.box.just = "center",
    legend.box.spacing = unit(0.3, "cm"),
    legend.box.background = element_blank(),
    legend.margin = margin(t = 0.2, b = 0.2),  # reduce top & bottom space
    # legend.spacing.y = unit(0.1, "cm"),
    # legend.key.height = unit(0.8, "cm"),
    # legend.key.width = unit(0.5, "cm"),
    legend.title = element_text(size = 8, face = "bold", color = col1),
    legend.text = element_text(size = 8, face = "bold"),
    legend.key = element_rect(fill = NA, colour = NA),
    strip.text = element_text(face = "bold", color = col1, size = 8),
    strip.background = element_rect(fill = "white", color = col6, linewidth = 0.5),
    axis.ticks.x = element_line(color = "grey"),
    axis.text.x = element_text(size = 8, color = col1),
    axis.text.y = element_text(size = 8, , color = col1),
    axis.title.x = element_text(size = 8, face = "bold", color = col1),      
    #panel.spacing.x = unit(0.3, "lines"),
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  ) +
  ggtitle(paste0("Risk of respiratory outcomes \nby hsct treatment status")) 

pdf(paste0(path_results,"/trt/forest_plot_trt_", ft,".pdf"), 
    width = 19.5, 
    height = 24)

print(p)    

dev.off()

ggsave(
  filename = paste0(path_results,"/trt/forest_plot_trt_", ft,".tiff"),
  plot = p,
  width = 8,         # fits within A4 with margins
  height = 5,       # portrait layout
  units = "in",
  dpi = 300,
  bg = "white",
  device = "tiff",
  compression = "lzw"  # optional: smaller file size, high quality
)

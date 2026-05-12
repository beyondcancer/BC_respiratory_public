f_up <- cr_finaldataforanalysis_respiratory %>%
  mutate(total_fup_yr = as.numeric(total_fup_days) / 365.25) %>%
  dplyr::select(setid, e_patid, exposed, cancer, total_fup_yr) %>%
  filter(cancer != "cns_b")

# Create named mapping (assumes cancer is coded 1:20)
cancer_label_map <- setNames(cancer_labels, cancersites)
# Calculate median follow-up by cancer type and exposure
f_up_summary <- f_up %>%
  mutate(
    exposed = factor(exposed, levels = c(0, 1), labels = c("Cancer-free", "Cancer survivor")),
    cancer_label = factor(cancer_label_map[as.character(cancer)], levels = cancer_labels)
  ) %>%
  group_by(cancer_label, exposed) %>%
  summarise(median_fup = median(total_fup_yr, na.rm = TRUE), .groups = "drop")

# Plot
ggplot(f_up_summary, aes(x = median_fup, y = fct_rev(cancer_label), fill = exposed)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(
    aes(label = sprintf("%.2f", median_fup)),
    position = position_dodge(width = 0.7),
    hjust = -0.2,
    size = 3.2
  ) +
  scale_fill_manual(values = c("Cancer-free" = col7, "Cancer survivor" = col4)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    x = "Median follow-up (years)",
    y = NULL,
    fill = NULL,
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.text = element_text(colour = col1, face = "bold"),
    axis.text.y = element_text(size = 11, colour = col1),
    axis.title.x = element_text(size = 12, colour = col1, face = "bold"))

    ggsave(
      filename = paste0(path_results,"/median_followup_by_cancer.png"),
      width = 210,
      height = 148,
      units = "mm",
      dpi = 500,
      bg = "white"
    )
  
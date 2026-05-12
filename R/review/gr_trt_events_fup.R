library(ggplot2)
library(tidyr)
library(dplyr)

################################################################################
# Plot 1: Number of events in treated population by cancer site
################################################################################

events_plot_data <- population_landmark %>%
  dplyr::select(cancer, outcome,
                n_events_chemo_treatment,
                n_events_radio_treatment) %>%
  pivot_longer(
    cols      = c(n_events_chemo_treatment, n_events_radio_treatment),
    names_to  = "treatment",
    values_to = "n_events"
  ) %>%
  mutate(
    treatment  = recode(treatment,
                        "n_events_chemo_treatment" = "Chemotherapy",
                        "n_events_radio_treatment" = "Radiotherapy"),
    cancer_label = factor(cancer_label_map[cancer],
                          levels = cancer_label_map[names(cancer_label_map) %in% unique(cancer)])
  )

plot_events <- ggplot(events_plot_data, aes(x = cancer_label, y = n_events, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = n_events),
            position = position_dodge(width = 0.9),
            vjust    = -0.5,
            size     = 3) +
  facet_wrap(~ outcome, scales = "free_y") +
  scale_fill_manual(values = c("Chemotherapy" = "#2c5f8a", "Radiotherapy" = "#e07b39")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title   = "Number of outcome events in treated population by cancer site",
    x       = "Cancer site",
    y       = "Number of events",
    fill    = "Treatment",
    caption = "Treated = received chemotherapy or radiotherapy within 1 year of cancer diagnosis"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#d6e4f0"),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "bottom"
  )

print(plot_events)

ggsave(paste0(path_results, "/trt/plot_events_treated_by_cancer_new.png"),
       plot   = plot_events,
       width  = 12,
       height = 8,
       dpi    = 300)

################################################################################
# Plot 2: Median follow-up time in treated population by cancer site
################################################################################

fup_plot_data <- population_landmark %>%
  dplyr::select(cancer, outcome,
                median_fup_chemo_treatment,
                median_fup_radio_treatment) %>%
  pivot_longer(
    cols      = c(median_fup_chemo_treatment, median_fup_radio_treatment),
    names_to  = "treatment",
    values_to = "median_fup"
  ) %>%
  mutate(
    treatment    = recode(treatment,
                          "median_fup_chemo_treatment" = "Chemotherapy",
                          "median_fup_radio_treatment" = "Radiotherapy"),
    cancer_label = factor(cancer_label_map[cancer],
                          levels = cancer_label_map[names(cancer_label_map) %in% unique(cancer)])
  )

plot_fup <- ggplot(fup_plot_data, aes(x = cancer_label, y = median_fup, fill = treatment)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_wrap(~ outcome, scales = "free_y") +
  scale_fill_manual(values = c("Chemotherapy" = "#2c5f8a", "Radiotherapy" = "#e07b39")) +
  labs(
    title   = "Median follow-up time in treated population by cancer site",
    x       = "Cancer site",
    y       = "Median follow-up (years)",
    fill    = "Treatment",
    caption = "Follow-up measured from landmark (1 year post cancer diagnosis)"
  ) +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "#d6e4f0"),
    strip.text       = element_text(face = "bold"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    legend.position  = "bottom"
  )

print(plot_fup)

ggsave(paste0(path_results, "/trt/plot_median_fup_treated_by_cancer_new.png"),
       plot   = plot_fup,
       width  = 12,
       height = 8,
       dpi    = 300)
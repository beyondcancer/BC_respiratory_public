library(ggplot2)
library(dplyr)
library(forcats)

smok <- paste0(path_datafiles_for_analysis, "/rev/an_cox_models_first_strat_smok.rds")
smok_ep <- paste0(path_datafiles_for_analysis, "/rev/an_cox_models_episodes_strat_smok.rds")
s<- readRDS(smok_ep)

datasets <- list("New-onset" = smok, "Exacerbation" = smok_ep)

for (name in names(datasets)) {
  
  for (outcome in outcomes_chronic) {
    
    col1 <- "#073B4C"
    col2 <- "#06D6A0"
    col3  <-"#EF476F"
    col4 <- "#FFD166"
    col5 <- "#118AB2" 
    col6 <- "#A2ACAB"
    col7 <- "#621244"
    
    i <- datasets[[name]]
    
    results <- readRDS(i)
    
    # format decimal places
    results <- results %>%
      mutate(
        hr = ifelse(is.na(hr), NA, as.numeric(sprintf("%.2f", hr))),
        ci_lower = as.numeric(sprintf("%.2f", ci_lower)),
        ci_upper = as.numeric(sprintf("%.2f", ci_upper)),
      ) 
    
    results <- results %>% group_by(cancer, outcome) %>% mutate(hr_label = paste0(hr[2],"(", ci_lower[2],"-", ci_upper[2],")")) %>% ungroup() 
    
    # Use mutate to create a new column with the corresponding labels
    results <- results %>%
      mutate(cancer_label_only = cancer_label_map[cancer]) %>%
     # mutate(events_combined = paste("[", events_exp,"/", events_unexp,"]", sep="")) %>%
     mutate(cancer_label = paste(cancer_label_only, hr_label, sep = "\n")) 
    
    
    #cat = deprivation,5 should be renamed to "deprivation
  
    results <- results %>% mutate(results, cat = if_else(cat == "Deprivation, 5", "Deprivation", cat))
    
    # rename variables
    results <- results %>%
      dplyr::rename(group = outcome,
                    labeltext = cancer_label,
                    mean = hr, 
                    lower = ci_lower,
                    upper = ci_upper) 
    
    
    
    label_gph <- c("Crude", "Adjusted")
    crude_adj <- results %>%
      filter(model_type == "Crude" | model_type == "Adjusted") %>%
      mutate(model_type= factor(model_type, levels = label_gph)) 
    
    time_diag <- results %>%
      filter(grepl("Adjusted|years", model_type))
    
    out_strat <- results %>%
      filter(group == outcome)
    
    mapped_outcome <- outcome_label_map[out_strat$group]
    mapped_cancer <- cancer_label_map[out_strat$cancer]
    
    out_strat <- out_strat %>% 
      filter(model_type != "Crude" & ! grepl("years", model_type) & !grepl("White|Black|Asian|Other", model_type)) %>%
      mutate( cat = if_else(model_type== "Adjusted", "Full population", cat))
    head(out_strat)  
    # # Create a header row for each category
    # header_rows <- out_strat %>%
    #   distinct(cat) %>%
    #   mutate(
    #     model_type = "",
    #     est = NA_real_,
    #     lower = NA_real_,
    #     upper = NA_real_,
    #     is_summary = TRUE  # Mark as a summary/header row
    #   )
    # 
    # # Combine header rows and data rows
    # cancer_strat_grouped <- bind_rows(
    #   header_rows,
    #   out_strat %>% mutate(is_summary = FALSE)
    # ) %>%
    #   arrange(cat, desc(is_summary))  # Ensures header comes before its group
    # 
    # 
    cancer_strat_grouped <- out_strat %>% dplyr::select(cancer, cat, model_type, mean, lower, upper)
    
    cancer_strat_grouped <- cancer_strat_grouped %>% mutate(cat = if_else(cat == "Deprivation, 5", "Deprivation", cat))
    
    cancer_strat_grouped <- out_strat %>%
      dplyr::select(cancer, cat, model_type, mean, lower, upper, lrtest) %>%
      filter(
        case_when(
          cancer %in% c("ova", "pro", "cer", "ute", "bre") ~ !model_type %in% c("Male", "Female"),
          TRUE ~ TRUE
        )
      )
  
    
    
    
    head(cancer_strat_grouped)
    
    cat_colors <- c(
      "Full population" = col1,
      "Age" = col5,      # Pantone 7476
      "Deprivation" = col2,  # Accent grey
      "Gender" = col3 ,
      "Smoking" = col4 # Soft fill
    )
    
    plot_data <- cancer_strat_grouped %>%
      mutate(
        cancer = fct_inorder(cancer),
        model_type = factor(model_type, levels = c(
          "Adjusted",
          "18-39", "40-59", ">=60",
          "Male", "Female",
          "1(Most deprived)", "2", "3", "4", "5(Least deprived)",
          "Never Smoker", "Ex-Smoker", "Smoker"
        )),
        model_type = fct_rev(model_type),  # reverse the correct custom order
        label_y = as.numeric(model_type)
      )
    
    plot_data <- plot_data %>%
      mutate(cancer = factor(cancer, levels = names(cancer_label_map)))
    
    # zebra <- plot_data %>%
    #   distinct(cancer, model_type, label_y) %>%
    #   arrange(cancer, desc(label_y)) %>%  # make sure rows are ordered correctly
    #   group_by(cancer) %>%
    #   mutate(fill = rep(c("#EAEAEA","#FFFFFF"), length.out = n())) %>%
    #   ungroup()
    # 
    # zebra <- plot_data %>%
    #   distinct(cancer, cat, model_type, label_y) %>%
    #   arrange(cancer, desc(label_y)) %>%
    #   group_by(cancer, cat) %>%
    #   mutate(cat_index = cur_group_id()) %>%
    #   ungroup() %>%
    #   mutate(fill = ifelse(cat_index %% 2 == 0, "#EAEAEA", "#FFFFFF"))
    # 
    adjusted_lines <- plot_data %>%
      filter(model_type == "Adjusted") %>%
      distinct(cancer, label_y)
    
    # Set CI limits
    xmin <- 0.25
    xmax <- 45
    
    plot_data <- plot_data %>%
      mutate(cat = factor(cat, levels = c("Full population", "Age", "Gender", "Deprivation", "Smoking")))
    
    # Mutate in capped limits and censoring flags
    plot_data <- plot_data %>%
      mutate(
        lower_capped = pmax(lower, xmin),
        upper_capped = pmin(upper, xmax),
        clipped_left = lower < xmin,
        clipped_right = upper > xmax
      )
    
    lrtest_labels <- plot_data %>%
      group_by(cancer, cat) %>%
      summarise(
        label_y  = median(label_y, na.rm = TRUE),
        lrtest   = first(lrtest),
        .groups  = "drop"
      ) %>%
      mutate(
        label_y = ifelse(
          cat == "Age" & cancer %in% c("ova", "pro", "cer", "ute", "bre"),
          10,
          label_y
        ),
        is_sig       = !is.na(lrtest) & lrtest < 0.05,
        lrtest_label = ifelse(
          is.na(lrtest),
          "",
          ifelse(
            lrtest < 0.001,
            paste0("<0.001", ifelse(is_sig, "*", "")),
            paste0("", formatC(lrtest, format = "f", digits = 3), ifelse(is_sig, "*", ""))
          )
        )
      )
    
    header_row <- lrtest_labels %>%
      group_by(cancer) %>%
      slice(1) %>%
      mutate(
        cat          = first(cat),
        lrtest_label = "Pval*",
        label_y      = ifelse(
          cancer %in% c("ova", "pro", "cer", "ute", "bre"),
          10 + 1.5,
          max(lrtest_labels$label_y, na.rm = TRUE) + 1.5
        )
      )
    lrtest_labels <- bind_rows(header_row, lrtest_labels)
    
    # Plot
    p <- ggplot() +
      # # Zebra background rows
      # geom_rect(
      #   data = zebra,
      #   aes(ymin = label_y - 0.5, ymax = label_y + 0.5),
      #   xmin = -Inf, xmax = Inf,
      #   fill = zebra$fill,
      #   inherit.aes = FALSE
      # ) +
      # Error bars
      geom_errorbarh( data = plot_data,
                      aes(xmin = lower_capped, xmax = upper_capped, y = model_type, color = cat),
                      height = 0.4, linewidth = 0.8
      ) +
      geom_point(
        data = subset(plot_data, clipped_left),
        aes(x = xmin, y = model_type, color = cat),
        shape = "<", size = 5,
        show.legend = FALSE
      ) +
      
      # Right clip: ">" shape, colored by category
      geom_point(
        data = subset(plot_data, clipped_right),
        aes(x = xmax, y = model_type, color = cat),
        shape = ">", size = 5,
        show.legend = FALSE
      ) +
      # Points, larger if summary
      geom_point(
        data = filter(plot_data, mean >= xmin & mean <= xmax),
        aes(x = mean, y = model_type, color = cat),
        shape = 16, size = 5
      ) +
      # Reference line
      geom_vline(
        xintercept = c(0.5, 2, 5, 10,40),
        linetype = "dashed",
        color = "grey80",
        linewidth = 0.3
      ) +
      geom_vline(
        xintercept = 1,
        linetype = "dashed",
        color = "red",
        linewidth = 0.3
      ) +
      scale_x_log10(
        limits = c(0.1, 45),
        breaks = c(0.5, 1, 2, 5, 10, 40),
        labels = scales::label_number(accuracy = 0.1)
      ) +
      scale_y_discrete(expand = expansion(add = c(0.3, 1))) +      # # Lines at each category boundary
      # geom_hline(
      #   data = plot_data %>%
      #     group_by(cancer, cat) %>%
      #     summarise(boundary = min(label_y) - 0.5, .groups = "drop"),
      #   aes(yintercept = boundary),
      #   color = "grey50", linetype = "solid", linewidth = 0.5,
      #   inherit.aes = FALSE
      # ) +
      scale_color_manual(values = c(
        "Age" = "#118AB2",
        "Gender" = "#EF476F",
        "Deprivation" = "#06D6A0",
        "Full population" = "#073B4C",
        "Smoking" = col4
      )) +
      scale_fill_manual(values = c(
        "Age" = "#118AB2",
        "Gender" = "#EF476F",
        "Deprivation" = "#06D6A0",
        "Full population" = "#073B4C",
        "Smoking" = col4
      )) +
      geom_label(
        data = lrtest_labels,
        aes(
          x     = 0.4,
          y     = label_y,
          label = lrtest_labels$lrtest_label
        ),
        hjust       = 1.1,
        size        = 3.5,
        color       = "#0D5257",
        fill        = alpha("white", 0.6),
        label.size  = 0,
        inherit.aes = FALSE
      ) +
      facet_wrap(
        ~ cancer,
        ncol = 4,
        nrow = 5, 
        scales = "free_y",
        labeller = labeller(cancer = cancer_label_map),
      ) +
      labs(
        x = "Hazard Ratio (95% CI)",
        y = NULL,
        fill = "Category",
        color = "Category"
      ) +
      theme_minimal(base_size = 16) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        strip.text = element_text(size = 14, face = "bold"),
        axis.text.y = element_text(size = 14, color = "#0D5257"),
        axis.text.x = element_text(size = 14, color = "#0D5257"),
        axis.title.x = element_text(size = 14, face = "bold", color = "#0D5257"),
        axis.ticks.x = element_line(color = "grey", size = 0.1),
        legend.position = "bottom",
        legend.key.size = unit(0.4, "cm"),   # smaller keys/icons
        legend.text = element_text(size = 7), # smaller label text
        legend.title = element_text(size = 8) # smaller title (or use element_blank() to remove)
      ) +
      ggtitle(paste0(mapped_outcome, " ", name)) +
      theme(
        plot.title = element_text(
          size = 10, face = "bold", hjust = 0.5, color = "#0D5257"
        )
      )
    
    print(p)
    
    
    ggsave(
      filename = paste0(path_results,"/forest_plot_strat_smok", outcome, "_", name, ".tiff"),
      plot = p,
      width = 17,         # fits within A4 with margins
      height = 20,       # portrait layout
      units = "in",
      dpi = 300,
      bg = "white",
      device = "tiff",
      compression = "lzw"  # optional: smaller file size, high quality
    )
    
  }
  
}


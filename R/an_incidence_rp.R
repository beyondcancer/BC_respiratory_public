###############################################################################
#R script author: Kirsty Andresen
#Date: 16 September 2024
#Description: Creates Adjusted incidence curves thought direct standarisation 
#for each cancer & outcome event (chronic diseases)
###############################################################################
# 
# cancersites <- c("mye", "leu")
# # # ########################
# # # ##FOR TESTING PURPOSES##
# cancersites <- "col"
# outcomes_chronic <- "asthma"
# # # #########################
# Load data

col1 <- "#073B4C"
col2 <- "#06D6A0"
col3  <-"#EF476F"
col4 <- "#FFD166"
col5 <- "#118AB2" 
col6 <- "#A2ACAB"
col7 <- "#621244"

  
grob_list <- list()
  # c_list <- list()
## Initiate cancer type loop
print(paste0("Initiating loop for cancer"))

for (site in cancersites) {
  print(paste0("cancer type: ", site))

  #Initiate outcome loop
  
  print(paste0("Initiating loop for "))
  
  for (outcome in outcomes_chronic) {
    print(paste0("Outcome: ", outcome))
    #cancer
    cr_outcome_an <- read.csv(paste0(path_datafiles_for_analysis, site, "_", outcome, "_cr_outcome_an.csv"))
    
    history_var <- paste0("b_", outcome)  # e.g., asthma_inc_dx
    event_var <- paste0("any_", outcome, "_inc_dx")  # e.g., asthma_inc_dx
    

    # covariates <- c("imd5", "bmi", "smokstatus","b_asthma", "b_copd", "b_bronch", "b_ild", "b_flu", "b_pneumo", "b_urti",
    #               "b_ckd", "b_hypertension", "b_cvdgrouped", "b_diab_cat", "b_cld", "b_multiple_sclerosis")
    # 
    covariates <- c("age_cat", "gender", "imd5", "bmi_cat", "smokstatus","b_asthma", "b_copd", "b_bronch",
                    "b_ckd", "b_hypertension", "b_cvdgrouped", "b_diab_cat", "b_cld", "b_multiple_sclerosis")
    
    covariates <- setdiff(covariates, history_var)
    
    if(site == "pro" | site == "bre" | site == "ute" | site == "ova" | site == "cer") {
      covariates <- setdiff(covariates, "gender")
    }
    
    factors <- c("age_cat", "imd5", "bmi_cat", "smokstatus")
    cr_outcome_an <- cr_outcome_an %>%
      mutate(across(all_of(factors), as.factor))
    
    stset <- as.formula(paste(
      "Surv(total_fup_var,",
      event_var,
      ") ~ exposed +",
      paste(covariates, collapse = "+")
    ))
    print(stset)
    
    #################################################################################
    #ROYSTON PAMAR MODEL
    #################################################################################
    
    ## Check numbers
    print(table(cr_outcome_an$exposed, cr_outcome_an[[event_var]]))
    
    ##Missing data summary
    
    missing_summary <- colSums(is.na(cr_outcome_an))
    missing_summary
    
    print(paste0("Fitting RP model"))
    rp_fit <- tryCatch ({
      stpm2(formula = stset,
            data = cr_outcome_an, 
            tvc = list(exposed = 4), 
            df = 4)
      
    }, error = function(e) {
      message("Error in fitting RP model: ", e$message)
      NULL
    }, warning = function(w) {
      message("Warning in fitting RP model: ", w$message)
      NULL  
    })
    
    if (is.null(rp_fit)) {
      next  # Skip the current loop iteration if the model fails
    }  

    sum_survivors <- cr_outcome_an %>% filter(exposed == 1) %>%
      dplyr::select(exposed,all_of(covariates))
    
    sum_controls <- sum_survivors %>% mutate(exposed = 0)
    
    max_time <- max(cr_outcome_an$total_fup_var[cr_outcome_an$exposed == 1], na.rm = TRUE)
    times <- seq(0, max_time, by = 0.25)
    
    
    # This function is now compatible with predictnl()
    predict_fun <- function(fit, newdata, time) {
      expanded <- tidyr::crossing(newdata, total_fup_var = time)
      preds <- predict(fit, newdata = expanded, type = "surv", se.fit = FALSE)
      avg_surv <- tapply(preds, expanded$total_fup_var, mean)
      names(avg_surv) <- time
      return(avg_surv)
    }
    
    # Function to compute standardized survival with CIs using predictnl
    standsurv_ci <- function(fit, newdata, times, exposed_value) {

      nd <- newdata
      nd$exposed <- exposed_value
      
      pn <- predictnl(fit, newdata = nd, fun = predict_fun, time = times)
      
      tibble(
        time = as.numeric(rownames(pn)),
        surv = pn$fit,
        lower = pn$Estimate - 1.96 * pn$SE,
        upper = pn$Estimate + 1.96 * pn$SE
      )
    }
    
    # Run for exposed and unexposed
    cancer_surv <- standsurv_ci(rp_fit, sum_survivors, times, exposed_value = 1)
    control_surv <- standsurv_ci(rp_fit, sum_survivors, times, exposed_value = 0)
    
    # Add group labels and cumulative incidence
    cancer_surv <- cancer_surv %>%
      mutate(group = "Cancer survivor",
             cum_inc = 1 - surv,
             cum_inc_l = 1 - upper,
             cum_inc_u = 1 - lower)
    
    control_surv <- control_surv %>%
      mutate(group = "Control",
             cum_inc = 1 - surv,
             cum_inc_l = 1 - upper,
             cum_inc_u = 1 - lower)
    
    # Combine for plotting
    plot_data <- bind_rows(cancer_surv, control_surv)
    # Combine both datasets
    mapped_site <- cancer_label_map[site]
    mapped_outcome <- outcome_label_map[outcome]
    name <- paste0(outcome, "_", site)
    
  p <-  ggplot(plot_data, aes(x = time, y = cum_inc, color = group, fill = group)) +
      geom_line(linewidth = 1.2) +  # Plot survival curves
      geom_ribbon(aes(ymin = cum_inc_l, ymax = cum_inc_u), alpha = 0.2, color = NA) +  # Add confidence intervals
      labs(
        title = paste0(mapped_site, "\n ", mapped_outcome),
        x = "Time (Years)",
        y = "Cumulative events 1-S(t)",
        color = "Group",
        fill = "Group"
      ) +
      scale_color_manual(values = c("Cancer survivor" = "#FFB81C", "Control" = "#621244")) +  # Custom line colors
      scale_fill_manual(values = c("Cancer survivor" = "#FFB81C", "Control" = "#621244")) +  # Custom fill colors
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 8, color = col1),
          legend.spacing = unit(5, "mm"),  # Reduce space between legend items
          plot.margin = margin(1, 1, 1, 1),  # Minimize extra space around the plot
          axis.title.x = element_text(margin = margin(t = 5), face = "bold", color = col1),  # Reduce x-axis title spacing
          axis.title.y = element_text(margin = margin(r = 5), face= "bold", color = col1),   # Reduce y-axis title spacing
          legend.position = "top")
    
    print(p)
    ggsave(paste0(path_results, "/cum_inc_rp_", site, "_", outcome, ".png"), p, width = 14, height = 10)
  
    # Convert to grob and store in the list
    grob_list[[name]] <- grid.grabExpr(print(p))

  }
}
  
####################

# ===== STEP 2: Split Grobs by Outcome and Site =====
grob_by_outcome <- setNames(lapply(unique(sub("_.*", "", names(grob_list))), function(outcome) {
  unname(grob_list[grep(paste0("^", outcome, "_"), names(grob_list))])
}), unique(sub("_.*", "", names(grob_list))))

grob_by_site <- setNames(lapply(unique(sub(".*_", "", names(grob_list))), function(site) {
  unname(grob_list[grep(paste0("_", site, "$"), names(grob_list))])
}), unique(sub(".*_", "", names(grob_list))))

# ===== STEP 3: Save Outcome-Based Panels =====
for (outcome in names(grob_by_outcome)) {
  png(file.path(path_results, paste0("cum_inc_rp_allsites_", outcome, ".png")), width = 2000, height = 1500, res = 150)
  grid.newpage()
  grid.text(mapped_site, gpar(fontsize = 12, fontface = "bold", col = col1))
  do.call(grid.arrange, c(grob_by_outcome[[outcome]], ncol = min(4, length(grob_by_outcome[[outcome]]))))
  dev.off()
}

# ===== STEP 4: Save Site-Based Panels  =====
for (site in names(grob_by_site)) {
  png(file.path(path_results, paste0("cum_inc_rp_alloutcomes_", site, ".png")), width = 2000, height = 1500, res = 150)
  grid.newpage()
  grid.text(mapped_outcome, gpar(fontsize = 12, fontface = "bold", col = col1))
  do.call(grid.arrange, c(grob_by_site[[site]], ncol = min(2, length(grob_by_site[[site]]))))
  dev.off()
}



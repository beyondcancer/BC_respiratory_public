###############################################################################
# CREATE TIME-TO-EVENT EPISODES FOR CR_FINALDATAFORANALYSIS_RESPIRATORY
###############################################################################

#Initialise empty dataframe to collect flowchart data
N_char <- data.frame(outcome = as.character(),
                     n_outcome_anytime = as.numeric(),
                     n_outcome_anytime_wexc = as.numeric(),
                     n_excluded_outfup = as.numeric(),
                     n_excluded_excb4inc = as.numeric(),
                     n_excluded_no_tar = as.numeric(),
                     n_outcome_valid_wexc = as.numeric(),
                     n_outcome_anytime_noexc = as.numeric(),
                     total_pats = as.numeric(),
                     total_events = as.numeric(),
                     error_rows = as.numeric())

#Set up outcome loop
 for (outcome in outcomes_chronic) {

  print(paste("Outcome: ", outcome, sep = " "))

  #Load in data  
  episode_data <- haven::read_dta(
  paste0(path_datafiles_for_episodes, "/cr_listpat_", outcome, "_exacerbations.dta"))

  # Define outcome-specific variables
  b_outcome <- paste0("b_", outcome)
  outcome_exc <- paste0(outcome, "_exc")
  exc_num <- paste0(outcome, "_exc_num")
  exc_contact <- paste0(outcome, "_exc_contact")
  dof_outcome_inc_dx <- paste0("dof_any_", outcome, "_inc_dx")
  outcome_inc_dx <- paste0("any_", outcome, "_inc_dx")

  #Extract patients who have had asthma at any time in their follow_up
  outcome_anytime <- cr_finaldataforanalysis_respiratory %>% 
                  mutate(
                    d_excfup = case_when(
                      !!sym(b_outcome) == 1 ~ as.Date(indexdate),  # If asthma at baseline, take index date
                      !!sym(b_outcome) == 0 & !!sym(outcome_inc_dx) == 0 ~ as.Date(NA),  # Ensure NA is a Date
                      TRUE ~ as.Date(!!sym(dof_outcome_inc_dx))    # Otherwise, take dof_outcome_inc_dx
                    )) %>%
                    filter(!is.na(d_excfup))

  #Select relevant variables used to filter invalid rows out of the exacerbation dataset
  doexit_data <-  outcome_anytime %>%  
                                   dplyr::select(setid, e_patid, exposed, indexdate, doexit, 
                                   !!sym(b_outcome), !!sym(outcome_inc_dx), !!sym(dof_outcome_inc_dx), d_excfup) 
  
  #Join to episode data
  exc_all <- episode_data %>%
    left_join(doexit_data) 
  
  #Patients with no exacerbations at any point
  outcome_anytime_noexc <- outcome_anytime %>% anti_join(episode_data, by = c("setid", "e_patid")) %>% distinct(setid, e_patid, .keep_all = TRUE) %>% dplyr::select(setid, e_patid, indexdate, doexit, d_excfup)
  
  #Collect counts 
  count_npat_outcome_anytime <- doexit_data %>% distinct(setid, e_patid) %>% nrow()
  count_npat_anytime_noexc_og <- outcome_anytime_noexc %>% distinct(setid, e_patid) %>% nrow()
  count_npat_anytime_noexc_og
  # Identify patients in the exacerbation dataset but not the cr_finaldataforanalisys dataset
  excluded <- exc_all %>% filter(is.na(doexit)) %>% distinct(setid, e_patid) 
  
  #remove exluded patients from episode dataset
  exc_anytime <- exc_all%>%
    filter(!is.na(doexit))

  #Filter based on date conditions (keep observations between indexdate and doexit) & censor maxdate at doexit    ## ADD BACK IN
  exc_anytime_wfup <- exc_anytime %>%
  filter((obsdate >= indexdate & obsdate <= doexit) & (maxdate >= indexdate)) %>%
    mutate(maxdate = pmin(doexit, maxdate))


 #create data frame of patients in data_int and not in data
  outfup <- exc_anytime %>% anti_join(exc_anytime_wfup, by =c("setid", "e_patid")) %>% distinct(setid, e_patid, .keep_all = TRUE) %>% dplyr::select(setid, e_patid, indexdate, doexit, d_excfup)
  
  # Add invalid episode patients to the outcome_anytime_noexc dataset
  outcome_anytime_noexc <- rbind(outcome_anytime_noexc, outfup) %>% distinct(setid, e_patid, .keep_all = TRUE) %>% dplyr::select(setid, e_patid, indexdate, doexit, d_excfup)
  
  count_npat_anytime_noexc <- outcome_anytime_noexc %>% distinct(setid, e_patid) %>% nrow()
  count_npat_anytime_noexc
  
  #Collect counts 
  count_excluded_events_notincohort <- excluded %>% nrow()
  
  count_events_exc <- exc_anytime %>%  nrow()
  count_events_exc_wfup <- exc_anytime_wfup %>% nrow()
  count_excluded_events_outfup <- count_events_exc - count_events_exc_wfup
  
  count_npats_exc <- exc_anytime %>% distinct(setid, e_patid) %>% nrow()
  count_npats_exc_fup <- exc_anytime_wfup %>% distinct(setid, e_patid) %>% nrow()
  count_excluded_npats_outfup <- count_npats_exc - count_npats_exc_fup
  

  #Create the episode indicator, flag and time at risk indicator. 
  
    #Tstart is 0 at:
        #indexdate if they had a history of outcome before indexdate, or
        #dof_asthma_inc_dx if they developed the outcome after indexdate
        #flag = 0 will eventually be excluded as it is time not at risk
 
  
#Remove EVENTS before first record of outcome
  exc_inc_wfup <- exc_anytime_wfup %>%
        filter(!(!!sym(b_outcome) == 0 & !is.na(!!sym(dof_outcome_inc_dx)) & (!!sym(dof_outcome_inc_dx) > obsdate)))
      
      b4inc <- exc_anytime_wfup %>% anti_join(exc_inc_wfup, by =c("setid", "e_patid")) %>% distinct(setid, e_patid, .keep_all = TRUE) %>% dplyr::select(setid, e_patid, indexdate, doexit, d_excfup)
      
      outcome_anytime_noexc <- rbind(outcome_anytime_noexc, b4inc) %>% distinct(setid, e_patid, .keep_all = TRUE) %>% dplyr::select(setid, e_patid, indexdate, doexit, d_excfup)
      
      count_npat_anytime_noexc <- outcome_anytime_noexc %>% distinct(setid, e_patid) %>% nrow()
      count_npat_anytime_noexc
      
      missing_summary <- colSums(is.na(exc_inc_wfup[ , ]))
      missing_summary
    
      count_events_excb4inc <-b4inc %>%  nrow()
      count_npats_b4bind = exc_inc_wfup %>% distinct(setid, e_patid) %>% nrow()
      
      exc_inc_wfup <- exc_inc_wfup %>%
        mutate(
          episode = 0, 
          flag = 0,
          tstart = as.numeric(obsdate - indexdate),
          time = as.numeric(maxdate - indexdate),
          type = "og"
        )
        
    #Step 1: Create rows with the episode indicator, flag indicator and time at risk for indexdate to first obsdate
      index_to_first_obs <- exc_inc_wfup %>%
        group_by(setid, e_patid) %>%
        filter(obsdate != d_excfup) %>%
        slice_head(n = 1) %>%
        mutate(
          episode = 1,
          flag = 1,
          tstart = as.numeric(d_excfup - indexdate),
          time = as.numeric(obsdate - indexdate),
          maxdate = obsdate
        ) %>%
        mutate(type = "first")
    
    #Step 2: Create the episode indicator, flag indicator and time at risk between episodes
      between_risk_periods <- exc_inc_wfup %>%
        arrange(e_patid, obsdate) %>%
        group_by(setid, e_patid) %>%
        filter(n() > 1) %>%
        mutate(
          tstart = lag(time, order_by = obsdate),
          time = as.numeric(obsdate - indexdate),
          episode = 1,
          flag = 1
        ) %>%
        filter(!is.na(tstart) & lag(maxdate) < doexit) %>%
        mutate(type = "between")
    
    #Step 3: Add row of data between the last episode and doexit if last episode maxdate < doexit
      #episode indicator = 0 (censored)
      last_maxdate_to_exit <- exc_inc_wfup %>%
        group_by(setid, e_patid) %>%
        filter(maxdate < doexit & doexit != d_excfup) %>%
        slice_tail(n = 1) %>%
        mutate(
          episode  = 0,
          flag = 1,
          !!sym(exc_num) := NA,
          !!sym(exc_contact) := NA,
          obsdate = maxdate,
          maxdate = doexit,
          tstart = as.numeric(obsdate - indexdate),
          time = as.numeric(maxdate - indexdate)
        ) %>%
        mutate(type = "last")

# # Check
# step1<- rbind(data, index_to_first_obs) %>%
#   arrange(setid, e_patid, tstart)
# step2 <- rbind(step1, between_risk_periods) %>%
#   arrange(setid, e_patid, tstart)
# step3 <- rbind(step2, last_maxdate_to_exit) %>%
#   arrange(setid, e_patid, tstart)
# step4 <- step3 %>% filter(flag == 1) %>%
#   arrange(setid, e_patid, tstart) %>%
#   ungroup()
# episode_dataset<- step4 

# pat1a <- step1 %>% filter(e_patid == )
# pat1b <- step2 %>% filter(e_patid == )
# pat1c <- step3 %>% filter(e_patid == )
# pat1d <- step4 %>%filter(e_patid == )
# 

  #Join rows of data 
    episode_dataset <- bind_rows(index_to_first_obs, between_risk_periods, last_maxdate_to_exit) %>%
      #filter(flag == 1) %>% #Remove time not at risk intervals
      #filter(time - tstart != 0) %>% #Remove intervals of time that are 1 day long (error rows that are created when obsdate = d_excfup)
      arrange(setid, e_patid, tstart) %>% ungroup()
    
    all <- episode_dataset %>% dplyr::select(setid, e_patid, indexdate, doexit, d_excfup, episode, tstart, time, type)
    
    #Count patient with no time at risk (all of their f-up is an exacerbation)
    invalid_tar_exc_pre <- exc_inc_wfup %>% anti_join(all, by = c("setid", "e_patid")) %>% distinct(setid, e_patid, .keep_all = TRUE) %>% dplyr::select(setid, e_patid, indexdate, obsdate, maxdate, doexit, d_excfup, type)
    
    add_tar <- invalid_tar_exc_pre %>%
      filter(maxdate - obsdate > 0 & obsdate == d_excfup) %>%
      mutate(episode = 1, 
             tstart = as.numeric(d_excfup - indexdate),
             time = as.numeric((d_excfup - indexdate)+ 0.5),
             type = "tar")
   
     all <- bind_rows(all, add_tar) %>% arrange(setid, e_patid, tstart)
    
    invalid_tar_exc <- invalid_tar_exc_pre %>% anti_join(add_tar, by = c("setid", "e_patid")) %>% distinct(setid, e_patid, .keep_all = TRUE) %>% dplyr::select(setid, e_patid, indexdate, obsdate, maxdate, doexit, d_excfup, type)
    
    count_invalid_tar_exc_pre <- invalid_tar_exc_pre %>% nrow()
    count_invalid_tar_exc <- invalid_tar_exc %>% nrow()
    
    #Calculate time for no exacerbation patients 
    outcome_anytime_noexc <- outcome_anytime_noexc %>%
      mutate(
        episode = 0,
        tstart = as.numeric(d_excfup - indexdate),
        time = as.numeric(doexit - indexdate),
        type = "noexc")
    
    outcome_anytime_noexc %>% summarise(
      missing_d_excfup = sum(is.na(d_excfup)),
      missing_doexit = sum(is.na(doexit))
    )
       
  #Join episode and noexacerbaion data    
  full_ep = bind_rows(all, outcome_anytime_noexc) %>% left_join(outcome_anytime) %>% arrange(setid, e_patid, tstart) 
  
  
  factors <- c("imd5", "bmi_cat", "smokstatus", "eth5_hes","age_cat")
  full_ep <- full_ep %>%
    mutate(across(all_of(factors), as.factor))

  #Count total with events and total without events
  count_npat <- full_ep %>% distinct(setid, e_patid) %>% nrow()
  count_events <- full_ep %>% filter(episode == 1) %>% nrow()
  count_events_valid_noexc <- full_ep %>% filter(is.na(episode)) %>% nrow()
    missing_summary
    
  count_events_valid_exc <- episode_dataset %>% filter(episode == 1) %>% nrow()
  count_npat_valid_noexc <- full_ep %>% filter(is.na(episode)) %>% distinct(setid, e_patid) %>% nrow()
  
  invalid_tar_noexc <-full_ep %>% filter(tstart >= time & type == "noexc") %>% distinct(setid, e_patid, .keep_all = TRUE) %>% dplyr::select(setid, e_patid, indexdate, doexit, d_excfup, type)
  all_invalid <- bind_rows(invalid_tar_exc, invalid_tar_noexc) %>% distinct(setid, e_patid, .keep_all = TRUE) %>% dplyr::select(setid, e_patid, indexdate, doexit, d_excfup, type)
  
  count_invalid_tar_noexc <- invalid_tar_noexc %>% nrow()
  count_invalid_tar = count_invalid_tar_exc + count_invalid_tar_noexc
  
  full_ep <- full_ep %>% anti_join(invalid_tar_noexc)
  
  error_events_df <-full_ep %>% filter(tstart >= time)
  error_events <- full_ep %>% filter(tstart >= time) %>% nrow()
  

#Collect flowchart data 
N_char <- rbind(N_char, data.frame(outcome = outcome,
                                   n_outcome_anytime = count_npat_outcome_anytime,
                                   n_outcome_no_exc = outcome_anytime_noexc %>% nrow(),
                                   n_episodes = all %>% distinct(e_patid, setid) %>% nrow(),
                                   n_excluded_outfup = count_excluded_npats_outfup,
                                   n_excluded_excb4inc = count_events_excb4inc,
                                   n_excluded_no_tar = count_invalid_tar,
                                   total_pats = (full_ep %>% distinct(setid, e_patid) %>% nrow()),
                                   total_events = (full_ep %>% filter(episode == 1) %>% nrow()),
                                   error_time = error_events)) 



#Save dataset 
print(full_ep)
print(N_char)

saveRDS(full_ep, paste0(path_datafiles_for_analysis, "cr_episodes_", outcome, "_new.rds"))

write.csv(N_char, paste0(path_results, "cr_episodes_respiratory_new_counts.csv"), row.names = FALSE)
}

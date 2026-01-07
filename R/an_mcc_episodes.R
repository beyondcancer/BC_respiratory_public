

#################################################################################
#                                                                               
#                        MEAN CUMULATIVE COUNT                                  
#                                                                               
# Last updated February 26, 2015.                                                  
#################################################################################
# R code to accompany:
#
#   Dong H, Robison LL, Leisenring WM, Martin LJ, Armstrong GT, Yasui Y. (in press).
#   Estimating the burden of recurrent events in the presence of competing risks: 
#   The method of mean cumulative count. American Journal of Epidemiology. 
#
# NOTE:
#   Our purpose of providing a simple hypothetical example and the computation code is 
#   that it would serve as a useful tutorial for researchers who want to learn how to 
#   apply the method of Mean Cumulative Count. Should you have any question, please conctact 
#   Huiru Dong by email: huiru@ualberta.ca.
#
# R version 3.1.1 (2014-07-10)
# Platform: x86_64-w64-mingw32/x64 (64-bit)


############ -------------  PARAMETERS -------------- #########
############----- id:    Participants's ID
############----- time:  Follow up time of participant to interest/competing-risk event/censoring
############----- cause: Event of interest is coded as 1; competing-risk event as 2; censoring as 0 
############----- type:  Character string specifying the type of calcultion. 
#				 Possible values are "MCC" (MCC calculation by equation) and "SCI" (sum of cumulative incidence). 



#################################################################################

col1 <- "#073B4C"
col2 <- "#06D6A0"
col3  <-"#EF476F"
col4 <- "#FFD166"
col5 <- "#118AB2" 
col6 <- "#A2ACAB"
col7 <- "#621244"

######### Run this only once to install the package #########
# install.packages("reshape")
# install.packages("etm")

library(reshape)
library(splines)
library(survival)
library(etm)

###################   This is the main function to run  ###################

MCC=function(id,time,cause,type)
{
  if(type=="MCC"){
    outdata=MCC_eq(id,time,cause)
  }
  if(type=="SCI"){
    outdata=SCI(id,time,cause)
  }
  return(outdata)
}

# ###################  Calculate the Sum of Cumulative Incidence  ########## 
# 
# SCI=function(id,time,cause)
# {
#   indata=data.frame(id,time,cause)
#   ii=order(indata$id,indata$time,-indata$cause)  ##if time is the same for 0/1 or 0/2 of the same person, make 0 last
#   sort_data=indata[ii,]
#   sort_id=sort_data$id
#   
#   #The time interval we are interested#
#   freq=1
#   time.interval=seq(min(time),max(time),freq)
# 
#   #Calculate the total number of participants: N#
#   Nodup_id=unique(id)
#   N=length(unique(id)) 
#   
#   
#   ######## Huiru's code has loop for each person and loop for each id, which was very slow (I actually cannot 
#   ### have it work when there in my CCSS --too much time) I will revise it to avoid looping pfr person 1: no. of person
#   
#   sort_data$first <- ave(sort_data$time, sort_data$id, FUN = seq_along)
#   M=max(sort_data$first) 		#maximum numer of events
#   
#   ### Take the last row so we know the maximum number of events per person
#   event.number=do.call(rbind, lapply(split(sort_data, sort_data$id), tail, 1))[,c("id","first")]
#   colnames(event.number)=c("id","maxE")
#   
#   alldata=merge(sort_data,event.number,by.x="id",by.y="id")
#   
#   #### make data for cumulative incidence, with dimension M*N
#   data.MCC=NULL
#   for(i in 1:M){
#     if(i==1){
#       data_temp=alldata[alldata$first==1,]
#       data_temp$m_event=i
#       data.MCC=rbind(data.MCC,data_temp)
#     }
#     if(i>1){
#       ## for those with i or more records, take the ith record
#       the_ith=alldata[alldata$first==i,]
#       ##for those with <i records, take the last record. If the last record is an event 1, then change it to 0;
#       ## because it contributed in data (i-1) already.
#       the_last=alldata[alldata$maxE<i & alldata$first==alldata$maxE,]
#       the_last$cause[the_last$cause==1]=0
#       CIdata_ith=rbind(the_ith,the_last)
#       CIdata_ith$m_event=i
#       data.MCC=rbind(data.MCC,CIdata_ith)
#     }
#   }
#   #Calculate cumulative incidence for M times event#
#   #Store the results in matrix MCC.base#
#   MCC.base=NULL
#   for (j in 1:M)
#   {
#     dataj=data.MCC[data.MCC$m_event==j,]
#     if(max(dataj$cause!=0)){  ##in example simple2, the 2nd time all cause=0
#       #Calculate cumulative incidence
#       P=summary(etmCIF(Surv(time, cause!=0)~ 1,data=dataj,etype=cause,failcode=1))[[1]]$'CIF 1'$P
#       #Calculate the change/increase
#       Deta_P=c(P[1],diff(P))  ## at time 1, the difference is P-0
#       
#       #Keep the time points
#       Time=summary(etmCIF(Surv(time, cause!=0)~ 1,data=dataj,etype=cause,failcode=1))[[1]]$'CIF 1'$time
#       #Combine all the information
#       if (is.null(P)==FALSE & is.null(Time)==FALSE)
#       { MCC.base=rbind(MCC.base,cbind(Time,P, Deta_P, CumI=j))}
#     }
#   }
#   
#   #Only keep the time points that affect MCC#
#   Nodup_MCC.base=MCC.base[!duplicated(MCC.base[,c(2,4)]),]
#   
#   #Sort by event dates#
#   jj=order(Nodup_MCC.base[,1])
#   sort_MCC.base=Nodup_MCC.base[jj,]
#   
#   #MCC is calculated by using the summation of cumulative incidences for all event of interest (first and recurrent)#
#   MCC=cumsum(sort_MCC.base[,3])
#   
#   combine_MCC=cbind(sort_MCC.base,MCC)
#   
#   #Only show the time points that have MCC change, and remove irrelevant information#
#   MCC.final=aggregate(combine_MCC[,5],list(combine_MCC[,1]),max)###
#   colnames(MCC.final)=c("Time","SumCIs")###
#   
#   return(MCC.final)
# }
# #-------------------------------------------------------------------------------------------------#



######### -----------------------  MCC function by equation ---------------##################
MCC_eq=function(id,time,cause)
{
  indata=data.frame(id,time,cause)
  ii=order(indata$id,indata$time,-indata$cause)  ##if time is the same for 0/1 or 0/2 of the same person, make 0 last
  sort_data=indata[ii,]
  #### take the last row of each person; if some people had the last row as event=1, need to add a row withe vent=0;
  lastr=aggregate(sort_data, list(sort_data$id), tail, 1)
  lastr=lastr[lastr$cause==1,]
  if(dim(lastr)[1]>0){ ## only if there are people whose last row is an event
    lastr$cause=0
    sort_data=rbind(sort_data,lastr[,c("id","time","cause")])
    ii=order(sort_data$id,sort_data$time,-sort_data$cause)  ##if time is the same for 0/1 or 0/2 of the same person, make 0 last
    sort_data=sort_data[ii,]
  }
  
  ntotal=length(unique(id)) 
  indata=sort_data
  time=sort_data$time
  cause=sort_data$cause
  id=sort_data$id
  count=rep(1,length(id))
  
  freq_cause=aggregate(count~time+cause, data=indata,sum)
  
  lifetable_1 <- cast(freq_cause, time~cause,value="count",fill=0)  
  colnames(lifetable_1)[colnames(lifetable_1)=="1"]="event"
  colnames(lifetable_1)[colnames(lifetable_1)=="0"]="censor"
  colnames(lifetable_1)[colnames(lifetable_1)=="2"]="cmprk"
  lifetable_1
  
  ### need to consider the situation that there is no censor 0, or event 1, competing risk 2 in the data at all.
  cause_in=unique(freq_cause$cause)  
  if(0 %in% cause_in==FALSE){
    censor=rep(0,dim(lifetable_1)[1])
    lifetable_1=data.frame(lifetable_1,censor)
  }
  if(1 %in% cause_in==FALSE){
    event=rep(0,dim(lifetable_1)[1])
    lifetable_1=data.frame(lifetable_1,event)
  }
  if(2 %in% cause_in==FALSE){
    cmprk=rep(0,dim(lifetable_1)[1])
    lifetable_1=data.frame(lifetable_1,cmprk)
  }
  
  
  ## n at risk at j = n at risk at j-1 -C-R, so get the running sum of C and R over time
  sum_censor=cumsum(lifetable_1[,"censor"])
  sum_cmprk=cumsum(lifetable_1[,"cmprk"])
  lifetable_2=cbind(lifetable_1,sum_censor,sum_cmprk)
  nrisk=ntotal-(sum_censor+sum_cmprk)
  nrisk_previous=c(ntotal,nrisk[1:(length(nrisk)-1)]) ## at the first time point, n at risk is the original number
  
  lifetable=data.frame(time=lifetable_1$time,nrisk=nrisk_previous,lifetable_1[,c("censor","event","cmprk")])
  lifetable
  
  surv_prob=1-lifetable$cmprk/lifetable$nrisk
  overall_surv=cumprod(surv_prob)
  overall_surv_previous=c(1,overall_surv[1:(length(overall_surv)-1)])   ###KM(Tj-1) is used in the MCC equation.
  Ave_events=overall_surv_previous*lifetable$event/lifetable$nrisk
  MCC=cumsum(Ave_events)
  
  MCCtable=data.frame(lifetable,MCC)
  MCC.final=do.call(rbind, lapply(split(MCCtable, MCCtable$MCC), head, 1))[,c("time","MCC")]
  rownames(MCC.final)=NULL
  MCCtable = data.frame(lifetable, MCC)
  MCC.final = do.call(rbind, lapply(split(MCCtable, MCCtable$MCC), head, 1))[,c("time", "MCC")]
  rownames(MCC.final) = NULL
  
  return(list(MCC_final = MCC.final, MCC_table = MCCtable))
}
#-------------------------------------------------------------------------------------------------#
# #### MASTER R SCRIPT ####
# options(scipen=999)
# #set working directory
# 
# setwd("J:/EHR-Working/Krishnan/Kirsty/BC_respiratory_outcomes/paths")
# 
# ## READ PATHS AND DATA
# 
# source("pathsR.R")


grob_mcc <- list()
grob_sci <- list()

# # ########TEST###########################
# # outcomes_chronic <- c("asthma", "copd")
# cancersites <- c("nhl")
# # ######################################

#cancersites <- c("lun", "mye", "leu", "nhl")


for (outcome in outcomes_chronic ) {

# ------------- Run the function with examples ------------#
#Run with data

episode <- readRDS(paste0(path_datafiles_for_analysis, "cr_episodes_", outcome, "_new.rds"))

for (site in cancersites) {

name <- paste0(outcome, "_", site)
  
cr_episode <- episode %>% dplyr::filter(cancer == site) 

cr_episode <- cr_episode %>% mutate(id = as.numeric(paste0(setid,e_patid)), 
                                    tstart_yrs = as.numeric(tstart/365.25),
                                    time_yrs = as.numeric(time/365.25)) 

data_mcc <- cr_episode %>% dplyr::select(setid, id, exposed, tstart_yrs, time_yrs, episode) %>% arrange(id, tstart_yrs) 

data_mcc_exp <- data_mcc %>% filter(exposed == 1)
data_mcc_unexp <- data_mcc %>% filter(exposed == 0)

mcc_exp <- MCC(id = data_mcc_exp$id, time = data_mcc_exp$time_yrs, cause = data_mcc_exp$episode, type = "MCC")
mcc_unexp <- MCC(id = data_mcc_unexp$id, time = data_mcc_unexp$time_yrs, cause = data_mcc_unexp$episode, type = "MCC")

mcc_table_exp <- mcc_exp$MCC_table
mcc_table_unexp <- mcc_unexp$MCC_table
mcc_exp <- mcc_exp$MCC_final
mcc_unexp <- mcc_unexp$MCC_final

# Define the target time points
target_times <- c(0, 5, 10, 15, 20)

risk_table_exp <- mcc_table_exp %>%
  mutate(closest_target = target_times[apply(abs(outer(time, target_times, "-")), 1, which.min)]) %>%
  group_by(closest_target) %>%
  slice_min(abs(time - closest_target), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(time, nrisk, MCC) %>%
  mutate(time = as.character(round(time, 0)),
         MCC = round(MCC, 1), 
         exposed =  1) 


# View results
print(risk_table_exp)

risk_table_unexp <- mcc_table_unexp %>%
  mutate(closest_target = target_times[apply(abs(outer(time, target_times, "-")), 1, which.min)]) %>%
  group_by(closest_target) %>%
  slice_min(abs(time - closest_target), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(time, nrisk, MCC) %>%
  mutate(time = as.character(round(time, 0)),
         MCC = round(MCC, 1), 
         exposed =  0)  # Round MCC values to 2 decimal places


risktable <- rbind(risk_table_exp, risk_table_unexp) 

risktable_plot <- risktable %>% 
  tidyr::pivot_wider(names_from = time, values_from = nrisk) %>%
  arrange(desc(exposed)) %>%
  rename_with(~gsub("nrisk_", "", .x), starts_with("nrisk_")) %>%  # Clean column names
  rename_with(~gsub("MCC_", "MCC_", .x), starts_with("MCC_")) %>% 
  group_by(exposed) %>%
  summarise(across(everything(), ~ na.omit(.)[1], .names = "{.col}"), .groups = "drop") %>%
  dplyr::select(-MCC) %>%
  mutate(exposed = ifelse(exposed == 1, "Cancer Survivor", "Control")) 

risktablelab <- risktable %>% dplyr::select(-MCC) %>%
  tidyr::pivot_wider(names_from = exposed, values_from = nrisk, names_prefix = "nrisk_") %>%
  mutate(nrisk = paste0(time,"\n", nrisk_1, "\n", nrisk_0))%>%  # Create formatted nrisk column
  dplyr::select(time, nrisk) 
# %>%
#   bind_rows(tibble(
#     time = "-1",  # Slightly before 0 to separate labels
#     nrisk = "  \nS:\nC:"  # Multi-line label
#   ))

# Keep only required columns
print(risktablelab)
custom_labels <- setNames(risktablelab$nrisk, risktablelab$time)
# Ensure Cancer Survivor is first
print(paste(site, outcome))
print(risktable)

risk_table_grob <- ggtexttable(
  risktable,
  rows = NULL,   # Remove row labels
  theme = ttheme("minimal")
)

# Add a group label to each MCC result
mcc_exp$Group <- "Cancer Survivor"
mcc_unexp$Group <- "Control"

# Combine the results into one dataframe
mcc_combined <- bind_rows(mcc_exp, mcc_unexp)

mapped_site <- cancer_label_map[site]
mapped_outcome <- outcome_label_map[outcome]
y_min <- min(mcc_combined$MCC, na.rm = TRUE)
y_max <- max(mcc_combined$MCC, na.rm = TRUE)
x_lab <- if (site == "thy" | site == "nhl"| site == "mye" | site == "leu") {
  "Time (Years)\nN at risk (Survivors/Control)"
} else {
  " "
}

# Plot MCC for all groups
mcc <- ggplot(mcc_combined, aes(x = time, y = MCC, color = Group)) +
  geom_step(size = 1.2) + # Step plot for cumulative MCC
  labs(
    title = paste0 (mapped_site, "\n", mapped_outcome),
    x = x_lab,
    y = " ",
    color = NULL
  ) +
  scale_color_manual(values = c("Cancer Survivor" = col4, "Control" = col7)) +  # Custom line colors
  ylim(-0, 25) +
  theme_minimal(base_size=12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12, color = col1),
      legend.spacing = unit(5, "mm"),  # Reduce space between legend items
        plot.margin = margin(1, 1, 1, 20),  # Minimize extra space around the plot
        axis.title.x = element_text( size= 10, face = "bold", color = col1),  # Reduce x-axis title spacing
      axis.text.x = element_text(size = 10) ,# Reduce y-axis title spacing
      legend.position = "off") + 
  scale_x_continuous(
    breaks = as.numeric(names(custom_labels)),  # Position of ticks
    labels = custom_labels  # Custom labels from risk table
  ) 

print(mcc)


grob_mcc[[name]] <- grid.grabExpr(print(mcc))

# ####################
# 
# sci_exp <- MCC(id = data_mcc_exp$id, time = data_mcc_exp$time_yrs, cause = data_mcc_exp$episode, type = "SCI")
# 
# sci_unexp <- MCC(id = data_mcc_unexp$id, time = data_mcc_unexp$time_yrs, cause = data_mcc_unexp$episode, type = "SCI")
# 
# # Add a group label to each sci result
# sci_exp$Group <- paste0(site, " survivor")
# sci_unexp$Group <- "control"
# 
# # Combine the results into one dataframe
# sci_combined <- bind_rows(sci_exp, sci_unexp)
# 
# # Plot sci for all groups
# sci <- ggplot(sci_combined, aes(x = Time, y = SumCIs, color = Group)) +
#   geom_step(size = 2) + # Step plot for cumulative sci
#   labs(
#     title = paste0(outcome, " exacerbations"),
#     x = "Time",
#     y = "Sum of cumulative incidence",
#     color = "Group"
#   ) +
#   theme_minimal() +
#   theme(plot.title = element_text(hjust = 0.5))
# 
# ggsave(paste0(path_results, "/sci_", outcome, site, ".png"), sci, width = 14, height = 10)
# 
# grob_sci[[name]] <- grid.grabExpr(print(sci))



} #cancer

} #outcome 


# ===== STEP 2: Split Grobs by Outcome and Site =====
grob_by_outcome <- setNames(lapply(unique(sub("_.*", "", names(grob_mcc))), function(outcome) {
  unname(grob_mcc[grep(paste0("^", outcome, "_"), names(grob_mcc))])
}), unique(sub("_.*", "", names(grob_mcc))))

grob_by_site <- setNames(lapply(unique(sub(".*_", "", names(grob_mcc))), function(site) {
  unname(grob_mcc[grep(paste0("_", site, "$"), names(grob_mcc))])
}), unique(sub(".*_", "", names(grob_mcc))))

# # ===== STEP 3: Save Outcome-Based Panels =====
# for (outcome in names(grob_by_outcome)) {
#     png(file.path(path_results, paste0("mcc_", outcome, ".png")), width = 2000, height = 1500, res = 150)
#     grid.newpage()
#     grid.text(mapped_site, gpar(fontsize = 12, fontface = "bold", col = col1))
#     do.call(grid.arrange, c(grob_by_outcome[[outcome]], ncol = min(4, length(grob_by_outcome[[outcome]]))))
#     dev.off()
#   }
# 
# # ===== STEP 4: Save Site-Based Panels  =====
# for (site in names(grob_by_site)) {
#     png(file.path(path_results, paste0("mcc_", site, ".png")), width = 2000, height = 1500, res = 150)
#     grid.newpage()
#     grid.text(mapped_outcome, gpar(fontsize = 12, fontface = "bold", col = col1))
#     do.call(grid.arrange, c(grob_by_site[[site]], ncol = min(2, length(grob_by_site[[site]]))))
#     dev.off()
# }

# all panels for lun, mye, leu and nhl

grob_list_panel <- c(grob_by_site[["lun"]], grob_by_site[["mye"]], grob_by_site[["leu"]], grob_by_site[["nhl"]])
png(file.path(path_results, paste0("mcc_panel.png")), width = 2000, height = 1500, res = 170)
grid.newpage()
grid.arrange(grobs = grob_list_panel, ncol = 4)
dev.off()


print(grob_list_panel)
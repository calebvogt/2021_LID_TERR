## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
## creates edge list of total time any two animals spent with each other per day (A-B and B-A are not repeated in the dataset, represented only once in random order)
## Testing: create edge list for every continuous interaction between mice A and B, regardless of whether mouse C enters or leaves. 
## Testing: Only break an A-B interaction bout if either A or B leaves the interaction. 
## Goal: Get time spent tolerating any other mouse. 

library(tidyverse)
library(data.table)
library(readxl)
library(asnipe)
library(igraph)

wd <- getwd()
dd <- paste(getwd(), "data", sep = "/")
output_fp <- paste(getwd(), "output", sep = "/")

filenames <- list.files(dd, pattern = "*MOVEBOUT_GBI.csv", full.names = T)
meta <- read.csv("data/metadata.csv")
social_data = lapply(filenames, fread) ## READ IN ALL FILES

# clean social data for triaged mice from social interaction bouts. 
aa = 1
for(aa in 1:length(social_data)){
  df <- social_data[[aa]] ## PULL OUT EACH TRIAL'S DATAFRAME
  
  df2 <- df %>% 
    select(!(V1)) %>% 
    filter(day %in% (1:20))
  
  print(aa)
  social_data[[aa]] <- df2
}



## CREATE LISTS OF NAMES FOR MATCHING COLUMNS
males <- meta %>% 
  filter(sex == "M") %>% 
  dplyr::select(name) %>% 
  filter(!is.na(name))
male_list <- dplyr::pull(males, name)

females <- meta %>% 
  filter(sex == "F", na.rm = TRUE) %>% 
  dplyr::select(name) %>% 
  filter(!is.na(name))
female_list <- dplyr::pull(females, name)

trial_stats <- list()
aa = 1
for(aa in 1:length(social_data)){
  df <- social_data[[aa]] ## PULL OUT EACH TRIAL'S DATAFRAME
  colnames(df)<-gsub("C57-M-","",colnames(df))
  colnames(df)<-gsub("C57-F-","",colnames(df))
  colnames(df)<-gsub("NYOB-M-","",colnames(df))
  colnames(df)<-gsub("NYOB-F-","",colnames(df))
  
  mice_names <- colnames(df[,10:ncol(df)]) ## get mouse column names starting at col 10
  focal = mice_names[1]
  # focal ="Amy"
  all_mouse_list <- list()
  done_mice <- c()
  flag1 = 1
  for(focal in mice_names[1:length(mice_names)-1]) { ## choose mouse focal column and compare to each successive mouse column. 
    print(focal)
    done_mice <- c(done_mice, focal) ## keep track of which mice have been used as a focal and do not use them as a comparison so as not to repeat observations. 
    df2 <- df %>% 
      filter((!!as.symbol(focal)) == 1) %>% # pull all rows where focal mouse is present.
      mutate(ID1 = focal) %>% 
      relocate(ID1)
    
    dyad_interactions <- list()
    partner_names <- mice_names[! mice_names %in% done_mice] ## remove columns of mice that have already been used as a focal. 
    partner <- partner_names[1] ## select current partner to look at. 
    partner = "Aphrodite"
    flag2 = 1
    for(partner in partner_names[1:length(partner_names)]) {
      print(paste(focal, partner))
      df3 <- df2 %>% 
        filter((!!as.symbol(partner)) == 1) %>% ## select rows where partner ==1. This line screws up field time order! fixed with arrange
        mutate(ID2 = partner) %>% 
        relocate(ID2) %>% 
        arrange(field_time_start) ## critical step
      
      ## troubleshooting tools
      # x <- df7 %>% 
      #   group_by(ID1,ID2,day) %>% 
      #   summarize(total = sum(duration_s))
      # write.table(df6, "clipboard-16384", sep="\t", row.names=FALSE, col.names = TRUE)
      ## troubleshooting tools
      
      df4 <- df3 %>% 
        mutate(start = as.POSIXct(field_time_start, origin="1970-01-01"), 
               stop = as.POSIXct(field_time_stop, origin="1970-01-01")) %>% 
        select(trial, day, zone, start, stop, ID1, ID2) %>% ## all good
        filter(!(start==stop))## remove rows with exact same start/stop time
        
      ## concatenate continuous flocking events for focal-partner pair into single row. i.e. the time those two mice ACTUALLY spent together, independent of other mouse behavior. 
      df5 <- df4 %>% 
        select(start,stop)
        
      temp <- unlist(df5) ## put all values into single vector
      myDupeVec <- unique(temp[duplicated(temp)]) ## get unique elements of vector
      noDupesList <- lapply(df5, function(i) i[!(i %in% myDupeVec)]) ## use the unique elements to remove non unique elements
       
      df6 <- do.call(data.frame, noDupesList) ## slap that baby back together. 
      
      df7 <- df4 %>% 
        merge(., df6, by="start") %>% 
        rename(stop=stop.y) %>% ## y because you want the stop from df6!!!
        select(trial, day, zone, ID1, ID2, start, stop) %>% 
        mutate(duration_s = as.numeric(stop)-as.numeric(start)) ## change to numeric so you dont get value as minutes
       
      dyad_interactions[[flag2]] <- df7
      flag2 = flag2 + 1
    }
    
    all_mouse_list[[flag1]] <- do.call("rbind", dyad_interactions)
    flag1 <- flag1 + 1
  }
  trial_stats[[aa]] <- do.call("rbind", all_mouse_list)

}
df8 <- do.call("rbind", trial_stats)
 
# test <- df8 %>% 
#   group_by(trial, ID1, ID2,day) %>% 
#   summarize(sum_duration_s = sum(duration_s))

write.csv(df8, "data/ALLTRIAL_MOVEBOUT_GBI_edgelist.csv")


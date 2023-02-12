## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(lubridate)
library(plotrix)

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
df <- rfid_data %>%
  dplyr::select(-c(V1)) %>%
  filter(noon_day %in% 1:20) %>% 
  filter(!(antenna %in% 13:16)) ## remove water towers

# options(scipen = 3)
df1 <- df %>% 
  filter(sex=="M") %>% 
  group_by(trial, antenna, noon_day) %>% 
  mutate(total_antenna_reads=n()) %>% 
  group_by(trial, strain, name, antenna, noon_day, total_antenna_reads) %>% 
  tally() %>% 
  rename(mus_antenna_reads = n) %>% 
  mutate(mus_perc_antenna_reads = (mus_antenna_reads/total_antenna_reads)*100) %>% 
  group_by(trial, antenna, noon_day) %>% 
  mutate(rank_order = rank(-mus_perc_antenna_reads)) %>% ## create rank order
  mutate(trial_antenna_day = paste(trial,antenna, noon_day, sep = "_")) %>% 
  arrange(trial_antenna_day, rank_order)

write.csv(df1, "data/ALLTRIAL_RFID_antenna_ownership.csv")


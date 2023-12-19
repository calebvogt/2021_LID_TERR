## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

# load dependencies -------------------------------------------------------
library(tidyverse)
library(data.table)
library(readxl)
library(lubridate)
library(funModeling)
library(ggplot2)
library(lme4)
library(lmerTest)
library(scales)
library(ggplot2)
library(viridis)
library(ggplot2)
library(viridis)
library(gtools)
library(svglite)

# Set working directory and load data -------------------------------------
wd <- setwd("C:/Users/caleb/Box/1_projects/2021_LID_TERR")
# wd <- setwd("C:/Users/ayalab/Box/LID_TERR_2021/2021_LID_TERR_CV_analysis")
meta <- read.csv("data/metadata.csv")
data <- as.data.frame(fread("data/ALLTRIAL_RFID.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

# RFID reads in zones over time line plot  ------------------------------
df <- data %>% 
  filter(noon_day %in% 1:20) %>% 
  filter(!(antenna %in% 13:16)) %>% 
  mutate(sex_treatment = paste(sex, treatment, sep = "-"))

## 1v2 zone males during first 5 days of the trial
custom_labels <- c("Cameron"="two-zone male 1","Carter"="two-zone male 2","Diego" = "two-zone male 3","Josiah" = "two-zone male 4",
                   "Isaiah"="one-zone male 1","Ibex"="one-zone male 2","Jeremy"="one-zone male 3","Dillard"="one-zone male 4")

df %>%  filter(sex_treatment=="M-early", noon_day %in% 1:5) %>% 
  mutate(name = factor(name,levels=c("Cameron","Carter","Diego","Josiah", ## set the order of the males
                                      "Isaiah","Ibex","Jeremy","Dillard"))) %>% 
ggplot(., aes(x=time_sec/86400+1,y=antenna)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  labs(x="night",y="zone") +
  # geom_vline(xintercept=6, linetype="solid", color="red", size=1) + ## add vertical line at invasion day 8 behind plotted data
  scale_x_continuous(breaks=seq(1,5,1)) +
  scale_y_continuous(breaks = seq(1,12,1), limits=c(1,12)) +
  theme_classic()+
  theme(strip.text = element_text(face = "bold", size = 10),
        panel.grid.major= element_line(color = "grey", size=0.5),
        panel.grid.major.x = element_blank(), ## remove vertical gridlines
        panel.background = element_blank()) +
  facet_wrap(~name,nrow=4,ncol=2,labeller=labeller(name=custom_labels))
# ggsave("output/plot.svg", device = "svg", bg = "white")


df %>%  filter(sex_treatment=="M-early", noon_day %in% 1:20) %>% 
  mutate(name = factor(name,levels=c("Cameron","Carter","Diego","Josiah", ## set the order of the males
                                     "Isaiah","Ibex","Jeremy","Dillard"))) %>% 
  ggplot(., aes(x=time_sec/86400+1,y=antenna)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  labs(x="night",y="zone") +
  geom_vline(xintercept=6, linetype="solid", color="red", size=1) + ## add vertical line at invasion day 
  scale_x_continuous(breaks=seq(0,20,5), limits=c(1,20.2)) +
  scale_y_continuous(breaks = seq(1,12,1), limits=c(1,12)) +
  theme_classic()+
  theme(strip.text = element_text(face = "bold", size = 10),
        panel.grid.major= element_line(color = "grey", size=0.5),
        panel.grid.major.x = element_blank(), ## remove vertical gridlines
        panel.background = element_blank()) +
  facet_wrap(~name,nrow=4,ncol=2,labeller=labeller(name=custom_labels))
# ggsave("output/plot.svg", device = "svg", bg = "white")


df %>% filter(trial == "T001", sex_treatment == "M-early") %>% 
  ggplot(., aes(x = time_sec/86400, y = antenna)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  xlab("Day") +
  scale_x_continuous(breaks=seq(1,20,1)) +
  scale_y_continuous(breaks = seq(1,12,1), limits=c(1,12)) +
  facet_wrap(~name) +
  theme(axis.text.x=element_text(angle = 65, hjust = 0.75))
# ggsave(filename = "output/rfid_T001_M_early_antenna_use.png", device = "png", bg = "white")


## T001_M_late
df %>% filter(trial == "T001", sex_treatment == "M-late") %>% 
  ggplot(., aes(x = time_sec/86400, y = antenna)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  labs(x="night",y="zone") +
  scale_x_continuous(breaks=seq(1,20,1)) +
  scale_y_continuous(breaks = seq(1,12,1), limits=c(1,12)) +
  theme_classic()+
  theme(strip.text = element_text(face = "bold", size = 10),
        panel.grid.major= element_line(color = "grey", size=0.5),
        panel.grid.major.x = element_blank(), ## remove vertical gridlines
        panel.background = element_blank(),
        axis.text.x=element_text(angle = 65, hjust = 0.75)) +
  facet_wrap(~code,nrow=3,ncol=4) 
# ggsave(filename = "output/rfid_T001_M_late_antenna_use.png", device = "png", bg = "white")


## T001_F_early
df %>% filter(trial == "T001", sex_treatment == "F-early") %>% 
  ggplot(., aes(x = time_sec/86400, y = antenna)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  xlab("Day") +
  scale_x_continuous(breaks=seq(1,20,1)) +
  scale_y_continuous(breaks = seq(1,12,1), limits=c(1,12)) +
  facet_wrap(~name) +
  theme(axis.text.x=element_text(angle = 65, hjust = 0.75))
ggsave(filename = "output/rfid_T001_F_early_antenna_use.png", device = "png", bg = "white")


## T001_F_late
df %>% filter(trial == "T001", sex_treatment == "F-late") %>% 
  ggplot(., aes(x = time_sec/86400, y = antenna)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  xlab("Day") +
  scale_x_continuous(breaks=seq(1,20,1)) +
  scale_y_continuous(breaks = seq(1,12,1), limits=c(1,12)) +
  facet_wrap(~name) +
  theme(axis.text.x=element_text(angle = 65, hjust = 0.75))
ggsave(filename = "output/rfid_T001_F_late_antenna_use.png", device = "png", bg = "white")



# Table S1: Numbers and Descriptive Stats  -------------------------------------
# Mean RFID reads per trial

df <- rfid_data %>% 
  group_by(trial) %>%  
  tally()
mean(df$n)
plotrix::std.error(df$n)

# Mean RFID reads per mouse per night 
df <- rfid_data %>% 
  group_by(name, noon_day) %>%  
  tally()
mean(df$n)
std.error(df$n)

## 
df <- rfid_data %>% 
  group_by(trial) %>% 
  tally()

##
df <- rfid_data %>% 
  filter(trial == "T001") %>% 
  group_by(noon_day, zone) %>% 
  tally()

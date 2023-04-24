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

# Set working directory and load data -------------------------------------
wd <- setwd("Z:/data/Caleb/2022_CV_LID_mong/rfid_data")
wd <- setwd("C:/Users/caleb/Box/1_CV_projects/2022_CV_LID_mong/rfid_data")

meta <- read.csv("metadata.csv")
data <- as.data.frame(fread("ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))


# RFID reads in zones over time line plot  ------------------------------
df <- data %>% 
  filter(noon_day %in% 1:4)

## individual
ggplot(df, aes(x = time_sec/86400+1, y = zone, group=1)) + ## x axis day number
  geom_point(na.rm=TRUE, size=0.5) +
  geom_line(linewidth=0.5) +
  xlab("Night") +
  scale_x_continuous(breaks=seq(1,4,1), limits=c(0.5,4.8)) +
  scale_y_continuous(breaks = seq(1,8,1), limits=c(0.5,8.5)) +
  theme_classic()+
  theme(axis.text.x=element_text(angle=60, hjust=1),
        strip.text = element_text(face = "bold", size = 10),
        panel.grid.major= element_line(color = "grey", size=0.5),
        panel.grid.major.x = element_blank(), ## remove vertical gridlines
        panel.background = element_blank()) +
  facet_wrap(~code)
# ggsave(filename = "output/rfid_eads_lineplot.png", device = "png", bg = "white")






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

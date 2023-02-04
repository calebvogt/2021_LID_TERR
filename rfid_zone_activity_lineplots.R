## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(readxl)
library(scales)

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
df <- rfid_data %>% 
  filter(noon_day %in% 1:20) %>% 
  filter(!(antenna %in% 13:16)) %>% 
  mutate(sex_treatment = paste(sex, treatment, sep = "-"))

## T001_M_early
df %>% filter(trial == "T001", sex_treatment == "M-early") %>% 
  ggplot(., aes(x = time_sec/86400, y = antenna)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  xlab("Day") +
  scale_x_continuous(breaks=seq(1,20,1)) +
  scale_y_continuous(breaks = seq(1,12,1), limits=c(1,12)) +
  facet_wrap(~name) +
  theme(axis.text.x=element_text(angle = 65, hjust = 0.75))
ggsave(filename = "output/rfid_T001_M_early_antenna_use.png", device = "png", bg = "white")


## T001_M_late
df %>% filter(trial == "T001", sex_treatment == "M-late") %>% 
  ggplot(., aes(x = time_sec/86400, y = antenna)) +
  geom_point(na.rm=TRUE, size=1) +
  geom_line() +
  xlab("Day") +
  scale_x_continuous(breaks=seq(1,20,1)) +
  scale_y_continuous(breaks = seq(1,12,1), limits=c(1,12)) +
  facet_wrap(~name) +
  theme(axis.text.x=element_text(angle = 65, hjust = 0.75))
ggsave(filename = "output/rfid_T001_M_late_antenna_use.png", device = "png", bg = "white")


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



# individual plots --------------------------------------------------------
ids <- unique(df$name)
i=ids[20]
for(i in ids[1:length(ids)]) {
  current_mouse <- print(i)
  move_df <- subset(df, name == current_mouse)
  current_mouse_info <- paste(unique(move_df$trial), unique(move_df$drop), unique(move_df$sex), sep = " ")

  
  p <- ggplot(move_df) +
    aes(x = field_time, y = antenna, color = reader.id) +
    geom_point(na.rm=TRUE, size=1) +
    ggtitle(paste(current_mouse_info, current_mouse, sep = " ")) +
    xlab("Date") + 
    ylab("antenna") +
    scale_x_datetime(breaks = "1 day", labels=date_format("%m-%d")) +
    scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(color = "black", size = 1, linetype = "solid"))
  
  p
  
  ggsave(filename=paste("rfid_antennareads",current_mouse_info, current_mouse, sep = "_"), 
         plot=p, 
         width = 5, 
         height = 4, 
         dpi = 300, 
         units = "in", 
         device='png', 
         path = output_fp)
} 


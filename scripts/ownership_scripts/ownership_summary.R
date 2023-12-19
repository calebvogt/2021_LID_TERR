## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggeffects)

data <- as.data.frame(fread("data/ALLTRIAL_RFID_antenna_ownership.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

df2 <- data %>% 
  filter(rank_order==1 & mus_perc_antenna_reads > 50) %>% 
  group_by(name, noon_day) %>% 
  tally() %>% 
  rename(zones_owned = n) %>% 
  arrange(noon_day)

(p <- ggplot(df2, aes(x = noon_day, y = zones_owned, color = name)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    # geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    # scale_x_continuous(limits = c(0.8,20.3), breaks = seq(1, 20, by = 1)) +
    # scale_y_continuous(limits = c(1,8), breaks = seq(1, 8, by = 1)) +
    theme_classic() +
    xlab("Day") +
    # ylab("Unique zones visited") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.21,0.8))
          legend.position = "right")
)
# ggsave(p, filename = "output/rfid_unique_zones_visited.png", device = "png", bg = "white")
# ggsave(p, filename = "output/rfid_unique_zone_antennas_visited.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


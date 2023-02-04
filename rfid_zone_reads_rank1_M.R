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

df0 <- rfid_data %>% 
  filter(!(antenna %in% 13:16)) ## remove water tower antennas from antenna calculation

# options(scipen = 3)
data <- df0 %>% 
  filter(sex=="M") %>% 
  group_by(trial, antenna, noon_day) %>% 
  mutate(total_antenna_reads=n()) %>% 
  group_by(trial, strain, name, antenna, noon_day, total_antenna_reads) %>% 
  tally() %>% 
  rename(mus_antenna_reads = n) %>% 
  mutate(mus_prop_totalantenna_reads = (mus_antenna_reads/total_antenna_reads)) %>% 
  group_by(trial, antenna, noon_day) %>% 
  mutate(rank_order = rank(-mus_prop_totalantenna_reads)) %>% ## create rank order
  mutate(trial_antenna_day = paste(trial,antenna, noon_day, sep = "_")) %>% 
  arrange(trial_antenna_day, rank_order)

# % male-sourced reads within each antenna belonging to top single male by strain --------
df <- data %>% 
  filter(rank_order==1) %>% ## grab top ranked male
  filter(noon_day %in% 1:5) %>% 
  mutate(trial_antenna = paste(trial, antenna, sep="_"))

## boxplot
(p <- ggplot(df, aes(x=as.factor(noon_day), y=mus_prop_totalantenna_reads*100, fill = as.factor(strain))) +
    geom_boxplot() +
    xlab("Day") +
    ylab("% male-sourced zone RFID reads from rank 1 male") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8),
          legend.position = "" )
)
ggsave(p, filename = "output/rfid_perc_reads_rank1_M_early_boxplot.png", device = "png", bg = "white")
# ggsave(p, filename = "output/rfid_perc_rank1_early_M_reads_boxplot.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


## lineplot of percentages
df2 <- df %>% 
  group_by(strain, noon_day) %>%
  summarise(mean = mean(mus_prop_totalantenna_reads*100), 
            sd = sd(mus_prop_totalantenna_reads*100),
            count = n(), 
            se = (sd/(sqrt(count))), 
            mode = mode(mus_prop_totalantenna_reads*100))

(p <- ggplot(df2, aes(x = noon_day, y = mean, color = strain)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    xlab("Day") +
    ylab("% male-sourced antenna RFID reads from rank 1 male") + 
    # scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) +
    scale_y_continuous(limits = c(70,100), breaks = seq(70, 100, by = 10)) +
    scale_color_manual(breaks = c("C57", "NYOB"),
                       values=c("sienna", "skyblue4")) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8),
          legend.position = "" )
)
ggsave(p, filename = "output/rfid_perc_rank1_M_reads.png", device = "png", bg = "white")
# ggsave(p, filename = "output/rfid_perc_rank1_M_reads.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")



## Stats
df$trial <- as.factor(df$trial)
df$strain <- as.factor(df$strain)
df$noon_day <- as.numeric(df$noon_day)

m1 = lmer(asin(sqrt(mus_prop_totalantenna_reads)) ~ strain*log(noon_day) + (1|name), data = df) 
AIC(m1)
qqnorm(resid(m1))
qqline(resid(m1))
summary(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


# % male-sourced reads within each antenna belonging to top single male on Day 3 --------
## pools C57 and outbred males. 
df <- data %>% 
  filter(rank_order==1, noon_day == 2) %>% ## 
  group_by(noon_day) %>% 
  summarize(mean(mus_prop_totalantenna_reads),
            sd = sd(mus_prop_totalantenna_reads),
            count = n(),
            se = (sd/sqrt(count)))

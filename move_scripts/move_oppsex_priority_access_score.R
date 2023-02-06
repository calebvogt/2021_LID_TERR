## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)

move_data <- as.data.frame(fread("data/ALLTRIAL_MOVEBOUT.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
meta <- read_excel("data/metadata.xlsx", sheet = 1, skip = 1)

df <- move_data %>% 
  filter(sex  == "M") %>%
  # filter(sex  == "F") %>%
  mutate(duration_min = duration_s / 60) %>% 
  group_by(noon_day, zone) %>% 
  tally(sum(duration_min)) %>% 
  mutate(day_zone = paste0(noon_day, "_", zone)) %>% 
  dplyr::rename(total_duration_min = n) %>% #specify dplyr due to conflict
  ungroup() %>% 
  dplyr::select(day_zone, total_duration_min)

df1 <- move_data %>% 
  filter(sex  == "M") %>%
  # filter(sex  == "F") %>%
  mutate(duration_min = duration_s / 60) %>% 
  mutate(day_zone = paste0(noon_day, "_", zone)) %>% 
  group_by(name, day_zone, noon_day, zone) %>% 
  tally(sum(duration_min)) %>%
  dplyr::rename(mus_duration_min = n)#specify dplyr due to conflict

df2 <- merge(df1, df, by = "day_zone", all = TRUE) # bring in rest of males that did not win any days. 
df3 <- merge(df2, meta, by = "name", all = FALSE) # bring in metadata

df4 <- df3 %>%
  mutate(name_phase = paste0(name,"_", phase)) %>% 
  # filter(phase == "early") %>% ## choose 
  filter(noon_day %in% (1:12)) # first 20 days

#summing daily adjusted scores and taking single number to avoid guessing of slice_max when scores vacillate in the negative range. 
## Need to change this such that the complete is done, and then the late mice day 1-5 is removed and not counted against them. 
## make sure if an animal is never seen again, their line dies. If mouse disappears then reappears, -1 penalty. 
df5 <- df4 %>% 
  mutate(mus_percent_capture = (mus_duration_min / total_duration_min)) %>% 
  group_by(name, noon_day) %>% 
  mutate(penalty = if(any(mus_percent_capture > 0.5)) 0 else -1) %>% # PENALTY #1: If on any day you dont capture greater than 50% for any zone, take off 1 point 
  group_by(name) %>%
  complete(noon_day = 1:20, fill = list(penalty = -1, mus_percent_capture = 0)) %>% # PENALTY #2: Not observed at all penalty. ## CHANGE DAY
  arrange(name, noon_day) %>% 
  group_by(name, noon_day) %>% 
  mutate(sum_daily_capture = sum(mus_percent_capture)) %>% # for each day sum percent capture pre-penalty application. 
  group_by(name,noon_day) %>% 
  mutate(disc_col = paste0(name, "_", noon_day, "_", sum_daily_capture)) %>% #Need to do this to drop repeated rows. 
  distinct(disc_col, .keep_all = TRUE) %>% # drop repeated daily sums on days with two zone rows. 
  mutate(sum_daily_capture_penalty = sum(sum_daily_capture+penalty)) %>% 
  group_by(name) %>% 
  mutate(csum_daily_capture_penalty = cumsum(sum_daily_capture_penalty)) %>% 
  select(!c(sex, code, phase)) # remove this information so you can bring it back in next. 

# bring in metadata
df5a <- merge(df5, meta_short, by = "name", all = FALSE)

df6 <- df5a %>% 
  # filter(!(phase == "early"))
  filter(!(phase == "late")) 
# filter(!(phase == "late" & noon_day %in% (1:5))) ## Remove days 1-5

# df7 <- merge(df6, meta, by = "name", all = FALSE) # bring in metadata

df8 <- df6 %>%
  mutate(name_phase = paste0(name, "_", phase)) %>%
  # dplyr::select(trial.x, strain_sex, strain,sex, name, noon_day, sum_daily_capture, penalty, sum_daily_capture_penalty, csum_daily_capture_penalty) %>% 
  # dplyr::rename(trial = trial.x) %>% 
  mutate(label = if_else(noon_day == max(noon_day), as.character(name_phase), NA_character_))

## create early_PAS days 1:20
# early_PAS <- df8

## create late_PAS scores days 6:20
# late_PAS <- df8

#merge early_PAS and late_PAS
df9 <- rbind(early_PAS, late_pas)


# output csv
# write.table(df8, "clipboard-16384", sep="\t", row.names=F, col.names=F) #col.names = false for second time. 
# df8 <- read_excel("Figure_Data.xlsx", sheet = "Fig2e-f") #note that this is combined male and female data. you will need to filter by sex to produce the graphs below. 

# line plot
library(ggrepel)
(p <- ggplot(late_pas, aes(x=noon_day, y=sum_daily_capture_penalty, group = name_phase, color = name_phase)) + #y=csum_adj_mus_percent_capture_score
    geom_line(size =1, alpha = 0.5) +
    scale_x_continuous(breaks = seq(1,20,by=1), limits = c(1,25)) + # Change
    # scale_y_continuous(breaks = seq(-5,15, by = 5), limits = c(-5,15)) +
    geom_label_repel(aes(label = label),
                     # xlim = c(-Inf, Inf),
                     # ylim = c(-Inf, Inf),
                     force = 1,
                     box.padding = 0.5,
                     segment.size = 1,
                     nudge_x = 0.5,
                     direction = "y",
                     hjust = "left",
                     segment.curvature = -0.5,
                     segment.ncp = 2,
                     na.rm = TRUE) +
    xlab("Day") +
    # ylab("Cumulative PA Score") +
    ylab("Daily Priority Access Score") +
    geom_hline(yintercept=0, linetype="dashed", color = "black", size = 1) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8), 
          legend.position = "none") 
)
# ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.75, height=2.15, bg = "transparent") #main figs

# density plot of day 10 priority access scores. 
(p <- ggplot(subset(df8, noon_day == 10), aes(x=csum_daily_capture_penalty, group = strain_sex, fill = strain_sex)) + 
    geom_density(adjust = 0.25,alpha = 0.8) +
    scale_x_continuous(breaks = seq(-10,20,by=5), limits = c(-10,21)) +
    scale_y_continuous(breaks = seq(0,0.2,by=0.05), limits = c(0, 0.18)) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    xlab("Day 10 Priority Access Score Distribution") +
    ylab("Density") +
    geom_vline(xintercept=0, linetype="dashed", color = "black", size = 0.75) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7), 
          legend.key.size = unit(0.25, 'cm'),
          legend.position = "top")
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.5, height=2.15, bg = "transparent") #main figs

library(diptest)
df9 <- df8 %>% 
  filter(noon_day == 10, strain == "NYOB") 
dip.test(df9$csum_daily_capture_penalty)
dS <- (dip(df9$csum_daily_capture_penalty, full.result = TRUE))
plot(dS)

#Combined male and female plot. 
mf <- read_excel("Data_Stats_v5.xlsx", sheet = "Fig2E-F")

# density plot of day 10 priority access scores. 
(p <- ggplot(subset(mf, noon_day == 10), aes(x=csum_daily_capture_penalty, group = sex, fill = sex)) + 
    geom_density(adjust = 0.25,alpha = 0.4) +
    scale_x_continuous(breaks = seq(-10,20,by=5), limits = c(-10,21)) +
    scale_y_continuous(breaks = seq(0,0.2,by=0.05), limits = c(0, 0.18)) +
    scale_fill_manual(breaks = c("F", "M"),
                      values=c("red", "blue")) +
    xlab("Day 10 Priority Access Score Distribution") +
    ylab("Density") +
    geom_vline(xintercept=0, linetype="dashed", color = "black", size = 0.75) +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8, face = "bold"),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8, face = "bold"),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6), 
          legend.key.size = unit(0.25, 'cm'),
          legend.position = "top")
)
ggsave(p, file = "Rplot.svg", device = "svg", output_fp, width=2.75, height=2.15, bg = "transparent") #main figs


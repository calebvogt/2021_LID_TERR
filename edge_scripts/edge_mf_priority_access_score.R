## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
## Using the edgelist, create an evolving metric for each male of how many females he captures. 

library(tidyverse)
library(ggrepel)
library(magrittr)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(diptest)

data <- read.csv("data/ALLTRIAL_MOVEBOUT_GBI_edgelist.csv")
metadata <- read.csv("data/metadata.csv")

df0 <- data %>%
  group_by(trial, ID1, ID2,day) %>%
  summarize(sum_duration_s = sum(duration_s))

df <- df0 %>% 
  merge(., metadata[,c("name","code", "sex", "strain")], by.x="ID1", by.y="name") %>% 
  rename(ID1_code=code, ID1_sex=sex) %>% 
  merge(., metadata[,c("name","code", "sex")], by.x="ID2", by.y="name") %>% 
  rename(ID2_code=code, ID2_sex=sex) %>% 
  select(trial, strain, day, sum_duration_s, ID1, ID1_code, ID1_sex,ID2, ID2_code, ID2_sex) %>% 
  mutate(grouptype=paste0(ID1_sex,ID2_sex)) %>% 
  mutate(strain=ifelse(strain=="NYOB", "Outbred","C57"))

x <- unique(df$ID1)
y <- unique(df$ID2)
mice <- c(x,y)
mice <- unique(mice)
mouse_list <- list()
aa = mice[1]
for(aa in mice[1:length(mice)]) { ## loop through mice
  df1 <- df %>% 
    filter(ID1 == aa | ID2 == aa) ## pull mouse if they are in either interaction
  days <- unique(df$day)
  day_list <- list()
  bb=1
  for(bb in days[1:length(days)]) {
    df2 <- df1 %>% 
      filter(day == bb) %>%
      mutate(focal = aa, ID1 = na_if(ID1, aa), ID2 = na_if(ID2, aa)) %>% 
      mutate(partner = coalesce(ID1,ID2)) %>% 
      select(trial, day, focal, partner, sum_duration_s) %>% 
      merge(., metadata, by.x = "focal", by.y = "name") %>%
      mutate(focal_sex = sex) %>% 
      select(trial, day, sum_duration_s, focal, focal_sex, partner) %>%  
      merge(., metadata, by.x = "partner", by.y = "name") %>%
      mutate(partner_sex = sex) %>% 
      select(trial, day, sum_duration_s, focal, focal_sex, partner, partner_sex) %>% 
      filter(!(focal_sex == partner_sex)) ## remove same sex observations, keep only oppsex interactions
    
    if(nrow(df2)>0) {
      df3 <- df2  %>% 
        arrange(desc(sum_duration_s)) %>% 
        mutate(rank_order = c(1:nrow(.))) %>% 
        select(trial, day, focal, partner, rank_order, sum_duration_s)
      
      day_list[[bb]] <- df3
      print(paste(aa, "day", bb, "finshed"))
    }
    ## note that nas are produced for animals 
  }
  mouse_list[[aa]] <- do.call("rbind", day_list)
}
df4 <- do.call("rbind", mouse_list)

df5 <- df4 %>% 
  group_by(focal) %>% 
  merge(., metadata, by.x = "focal", by.y = "name") %>%
  select(trial, day, focal, code, sex, partner, rank_order,sum_duration_s)
  
df5a <- df5 %>% ## total time females spend with oppsex males per day
  filter(sex  == "F") %>% ## selects focal females and their oppsex obs
  group_by(focal, day) %>% 
  tally(sum(sum_duration_s)) %>% 
  mutate(f_day = paste0(focal, "_", day)) %>% 
  dplyr::rename(total_male_s = n) %>% #specify dplyr due to conflict
  ungroup() %>% 
  select(f_day, total_male_s)

options(scipen=999)
df6 <- df5 %>% ## time males spent in a given zone per day
  filter(sex  == "M") %>%
  mutate(f_day = paste0(partner, "_", day)) %>% 
  merge(., df5a, by="f_day") %>% 
  mutate(f_prop_capture = sum_duration_s/total_male_s) %>% 
  arrange(day, focal) %>% 
  group_by(focal, day) %>% 
  mutate(total_mf_s = sum(sum_duration_s), 
         # captured = ifelse(f_prop_capture>0.5,1,0), ## determine which females were captured 
         num_f_met = n(), ## num females met per day
         num_f_captured = sum(f_prop_capture>0.50), ## IMPORTANT: set f capture threshold
         prop_num_f_captured_met=num_f_captured/num_f_met,
         mean_daily_capture = mean(f_prop_capture),# for each day sum proportion female capture for all females
         sum_daily_capture = sum(f_prop_capture)) %>% ## for each day sum proportion female capture for all females
  group_by(focal, day) %>% 
  mutate(disc_col = paste0(focal, "_", day, "_", sum_daily_capture)) %>% #Need to do this to drop repeated rows. 
  distinct(disc_col, .keep_all = TRUE) %>% # drop repeated daily sums on days with two zone rows. 
  arrange(focal,day) %>% 
  group_by(focal) %>% 
  mutate(csum_daily_capture = cumsum(sum_daily_capture)) %>% 
  merge(., metadata[,c("name","treatment", "num_terr")], by.x="focal", by.y="name") %>% 
  select(trial,focal,code,sex,treatment,num_terr,day,total_mf_s,num_f_met,num_f_captured,prop_num_f_captured_met,mean_daily_capture,sum_daily_capture,csum_daily_capture) %>% 
  arrange(treatment,focal,day)

# # of females met --------------------
df <- df6 %>% 
  filter(treatment=="early")
  
df2 <- df %>% 
  group_by(num_terr,day) %>%
  summarise(mean = mean(num_f_met), 
            sd = sd(num_f_met), 
            count = n(), 
            sem = (sd/(sqrt(count))))

ggplot(df2, aes(x=day, y=mean, color = as.factor(num_terr))) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin =mean-sem,ymax=mean+sem), width = 0.2) +
    geom_vline(xintercept=6,color="black",linetype="dashed") +
    # scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
    # scale_y_continuous(limits = c(1,20), breaks = seq(2, 20, by = 2)) +
    scale_color_manual(breaks = c("1","2"),
                       values=c("goldenrod4","darkorchid4")) +
    theme_classic() +
    xlab("Night") +
    ylab("# females met") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.8,0.8))
          # legend.position = "none")
ggsave(file="output/edge_mf_pas_num_f_met.svg",device="svg",unit="in",width=3,height=3,bg = "transparent") 

# total male time with females --------------------
df <- df6 %>% 
  filter(treatment=="early")

df2 <- df %>% 
  group_by(num_terr,day) %>%
  summarise(mean = mean(total_mf_s), 
            sd = sd(total_mf_s), 
            count = n(), 
            sem = (sd/(sqrt(count))))

ggplot(df2, aes(x=day, y=mean, color = as.factor(num_terr))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin =mean-sem,ymax=mean+sem), width = 0.2) +
  geom_vline(xintercept=6,color="black",linetype="dashed") +
  # scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
  # scale_y_continuous(limits = c(1,20), breaks = seq(2, 20, by = 2)) +
  scale_color_manual(breaks = c("1","2"),
                     values=c("goldenrod4","darkorchid4")) +
  theme_classic() +
  xlab("Night") +
  ylab("Time with females (s)") +
  theme(axis.text.x = element_text(color = "black", size = 8),
        axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
        axis.text.y = element_text(color = "black", size = 8),
        axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"), 
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.background = element_rect(fill='transparent'),
        legend.position = c(0.8,0.8))
# legend.position = "none")
ggsave(file="output/edge_mf_time.svg",device="svg",unit="in",width=3,height=3,bg = "transparent") 




# # of females captured --------------------
df <- df6 %>% 
  filter(treatment=="early")

df2 <- df %>% 
  group_by(num_terr,day) %>%
  summarise(mean = mean(num_f_captured), 
            sd = sd(num_f_captured), 
            count = n(), 
            sem = (sd/(sqrt(count))))

ggplot(df2, aes(x=day, y=mean, color = as.factor(num_terr))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin =mean-sem,ymax=mean+sem), width = 0.2) +
  geom_vline(xintercept=6,color="black",linetype="dashed") +
  # scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
  # scale_y_continuous(limits = c(1,20), breaks = seq(2, 20, by = 2)) +
  scale_color_manual(breaks = c("1","2"),
                     values=c("goldenrod4","darkorchid4")) +
  theme_classic() +
  xlab("Night") +
  ylab("# of females captured") +
  theme(axis.text.x = element_text(color = "black", size = 8),
        axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
        axis.text.y = element_text(color = "black", size = 8),
        axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"), 
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.background = element_rect(fill='transparent'),
        legend.position = c(0.8,0.8))
# legend.position = "none")
ggsave(file="output/edge_num_f_captured.svg",device="svg",unit="in",width=3,height=3,bg = "transparent") 



# proportion of females captured out of total met --------------------
df <- df6 %>% 
  filter(treatment=="early")

df2 <- df %>% 
  group_by(num_terr,day) %>%
  summarise(mean = mean(prop_num_f_captured_met), 
            sd = sd(prop_num_f_captured_met), 
            count = n(), 
            sem = (sd/(sqrt(count))))

ggplot(df2, aes(x=day, y=mean, color = as.factor(num_terr))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin =mean-sem,ymax=mean+sem), width = 0.2) +
  geom_vline(xintercept=6,color="black",linetype="dashed") +
  # scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
  # scale_y_continuous(limits = c(1,20), breaks = seq(2, 20, by = 2)) +
  scale_color_manual(breaks = c("1","2"),
                     values=c("goldenrod4","darkorchid4")) +
  theme_classic() +
  xlab("Night") +
  ylab("Proportion of met females captured") +
  theme(axis.text.x = element_text(color = "black", size = 8),
        axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
        axis.text.y = element_text(color = "black", size = 8),
        axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"), 
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.background = element_rect(fill='transparent'),
        legend.position = c(0.8,0.2))
# legend.position = "none")
ggsave(file="output/edge_prop_f_captured.svg",device="svg",unit="in",width=3,height=3,bg = "transparent") 



# mean priority female access score --------------------
df <- df6 %>% 
  filter(treatment=="early")

df2 <- df %>% 
  group_by(num_terr,day) %>%
  summarise(mean = mean(mean_daily_capture), 
            sd = sd(mean_daily_capture), 
            count = n(), 
            sem = (sd/(sqrt(count))))

ggplot(df2, aes(x=day, y=mean, color = as.factor(num_terr))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin =mean-sem,ymax=mean+sem), width = 0.2) +
  geom_vline(xintercept=6,color="black",linetype="dashed") +
  # scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
  # scale_y_continuous(limits = c(1,20), breaks = seq(2, 20, by = 2)) +
  scale_color_manual(breaks = c("1","2"),
                     values=c("goldenrod4","darkorchid4")) +
  theme_classic() +
  xlab("Night") +
  ylab("Mean priority female access score") +
  theme(axis.text.x = element_text(color = "black", size = 8),
        axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
        axis.text.y = element_text(color = "black", size = 8),
        axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"), 
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.background = element_rect(fill='transparent'),
        legend.position = c(0.8,0.2))
# legend.position = "none")
ggsave(file="output/edge_m_priority_f_access_score.svg",device="svg",unit="in",width=3,height=3,bg = "transparent") 




# cumulative sum daily capture score (Priority access score) --------------------
df <- df6 %>% 
  filter(treatment=="early")

df2 <- df %>% 
  group_by(num_terr,day) %>%
  summarise(mean = mean(csum_daily_capture), 
            sd = sd(csum_daily_capture), 
            count = n(), 
            sem = (sd/(sqrt(count))))

ggplot(df2, aes(x=day, y=mean, color = as.factor(num_terr))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) +
  geom_errorbar(aes(ymin =mean-sem,ymax=mean+sem), width = 0.2) +
  geom_vline(xintercept=6,color="black",linetype="dashed") +
  # scale_x_continuous(limits = c(0.8,10.3), breaks = seq(1, 10, by = 1)) + 
  # scale_y_continuous(limits = c(1,20), breaks = seq(2, 20, by = 2)) +
  scale_color_manual(breaks = c("1","2"),
                     values=c("goldenrod4","darkorchid4")) +
  theme_classic() +
  xlab("Night") +
  ylab("Cumulative priority female access score") +
  theme(axis.text.x = element_text(color = "black", size = 8),
        axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
        axis.text.y = element_text(color = "black", size = 8),
        axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"), 
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.background = element_rect(fill='transparent'),
        legend.position = c(0.8,0.2))
# legend.position = "none")
ggsave(file="output/edge_csum_priority_f_access_score.svg",device="svg",unit="in",width=3,height=3,bg = "transparent") 



# delete ------------------------------------------------------------------

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


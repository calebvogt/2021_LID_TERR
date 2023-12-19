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

# edge_num_mf_captured --------------------
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



# edge_prop_mf_captured --------------------
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




## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)
library(asnipe)
library(igraph)
library(lme4)
library(lmerTest)

# load  & clean data ---------------------------------------------------------------
filenames <- list.files("data/", pattern = "*MOVEBOUT_GBI.csv", full.names = T)
social_data = lapply(filenames, fread) ## READ IN ALL FILES
meta <- read.csv("data/metadata.csv")

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

## count num. males and females met for each mouse per day
trial_stats <- list()
aa = 1
for(aa in 1:length(social_data)){
  df <- social_data[[aa]] ## PULL OUT EACH TRIAL'S DATAFRAME
  
  col_ids <- colnames(df[,10:ncol(df)]) ## get mouse column names starting at col 10
  bb = col_ids[15]
  all_mouse_list <- list()
  first_flag = 1
  bb="C57-M-Josiah"
  for(bb in col_ids[1:length(col_ids)]) {
    print(paste(bb, first_flag))
    df2 <- df %>% 
      filter((!!as.symbol(bb)) == 1) %>% # pull all rows where bb mouse is present. 
      mutate(name = bb) %>% 
      relocate(name)
    
    df3 <- df2 %>% ## get counts of spatiotemporal overlapping movebout events by day. useful!
      select(day, !!col_ids) %>% 
      group_by(day) %>% 
      summarise_all(sum) 
    
    m_met <- df3 %>% 
      mutate_at(vars(col_ids), funs(ifelse(. !=0, 1,0))) %>% ## convert non-0 entries to 1s
      select(day, contains("-M-"), -!!bb) %>% 
      mutate(sum = rowSums(.[-1])) %>% 
      rename(m_met=sum) %>%   
      select(day, m_met)
   
    f_met <- df3 %>% 
      mutate_at(vars(col_ids), funs(ifelse(. !=0, 1,0))) %>% ## convert non-0 entries to 1s
      select(day, contains("-F-"), -!!bb) %>% 
      mutate(sum = rowSums(.[-1])) %>% 
      rename(f_met=sum) %>%   
      select(day, f_met)
    
    df4 <- merge(m_met,f_met, by="day")
    
    df4$name <- bb
    
    all_mouse_list[[first_flag]] <- df4 %>% 
      mutate(name = gsub("C57-M-", "", name)) %>% 
      mutate(name = gsub("C57-F-", "", name)) %>% 
      merge(., meta[,c("name", "sex", "treatment", "num_terr")], by="name") %>% 
      select(name, sex,day, treatment, num_terr, m_met, f_met)
    first_flag <- first_flag+1
  }
  df5 <- do.call("rbind", all_mouse_list)
  
  trial_stats[[aa]] <- df5
}
data <- do.call("rbind", trial_stats)


# gbi_number_mf_encountered --------------------
df1 <- data %>% 
  mutate(num_terr = as.factor(num_terr)) %>% 
  filter(sex == "M") %>% 
  filter(num_terr %in% 0:2)

df2 <- df1 %>% 
  group_by(num_terr,day) %>%
  summarise(mean = mean(f_met), 
            sd = sd(f_met), 
            count = n(), 
            sem = (sd/(sqrt(count))))

ggplot(df2, aes(x=day, y=mean, color = as.factor(num_terr))) + 
  geom_line(size = 0.75) + 
  geom_point(size = 1.5) +
  # geom_errorbar(aes(ymin =mean-sem,ymax=mean+sem), width = 0.2) +
  geom_vline(xintercept=6,color="black",linetype="dashed") +
  scale_x_continuous(limits = c(0.8,22), breaks = seq(0, 20, by = 5)) +
  # scale_y_continuous(limits = c(1,20), breaks = seq(2, 20, by = 2)) +
  scale_color_manual(breaks = c("1","2", "0"),
                     values=c("goldenrod4","darkorchid4", "black"), 
                     labels = c("one-zone males", "two-zone males", "males without territory")) +
  theme_classic() +
  xlab("night") +
  ylab("# females encountered") +
  theme(axis.text.x = element_text(color = "black", size = 8),
        axis.title.x = element_text(color = "black", size = 8), 
        axis.text.y = element_text(color = "black", size = 8),
        axis.title.y = element_text(color = "black", size = 8), 
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent"), 
        panel.border = element_rect(colour = "black", fill=NA, size=0.25),
        legend.title = element_blank(),
        legend.text = element_text(size=6),
        legend.background = element_rect(fill='transparent'),
        legend.position = c(0.7,0.8))
# legend.position = "none")
ggsave(file="output/gbi_number_mf_encountered.svg",device="svg",unit="in",width=3,height=3,bg = "transparent") 


levels(df1$num_terr)

m1 = lmer(f_met ~ num_terr*day + (1|name), data = filter(df1, day %in% 1:5, num_terr %in% 1:2)) 
summary(m1)
AIC(m1)

m2 = lmer(f_met ~ num_terr*day + (1|name), data = filter(df1, day %in% 6:20, num_terr %in% 1:2)) 
summary(m2)
AIC(m2)

m3 = lmer(f_met ~ num_terr*day + (1|name), data = filter(df1, day %in% 6:20, num_terr %in% 0:2)) 
summary(m3)
AIC(m3)

write.table(summary(m2)$coef, "clipboard-16384", sep="\t", row.names=TRUE, col.names = TRUE)


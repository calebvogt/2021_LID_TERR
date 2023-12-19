## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(ggpubr)
library(lme4) 
library(lmerTest)

wd <- setwd("C:/Users/caleb/Box/LID_TERR_2021/2021_MZ_CV_LID_TERR_rfid_analysis")
data <- read.csv("data/ALLTRIAL_SNA_node_stats.csv")

df1 <- data %>% 
  mutate(num_terr = as.factor(num_terr)) %>% 
  filter(sex == "M") %>% 
  filter(num_terr %in% 0:2)
  
df2 <- df1 %>% 
  group_by(num_terr, day) %>%
  summarise(mean = mean(node_oppsex_edge_strength), sd = sd(node_oppsex_edge_strength), count = n(), sem = (sd/(sqrt(count))))

ggplot(df2, aes(x=day, y=mean, color=as.factor(num_terr))) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    geom_vline(xintercept=6,color="black",linetype="dashed") +
    # scale_x_continuous(limits = c(0.8,5.3), breaks = seq(1, 10, by = 1)) +
    # scale_y_continuous(breaks = seq(1, 20, by = 1)) +
    scale_color_manual(breaks = c("1","2", "0"),
                       values=c("goldenrod4","darkorchid4", "black"), 
                       labels = c("one-zone male", "two-zone male", "zero-zone male")) +
    theme_classic() +
    xlab("night") +
    ylab("social network node edge strength with females") +
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
          # legend.position = c(0.8,0.8), 
          legend.position = "")
ggsave(file="output/sna_node_oppsex_edge_strength.svg",device="svg",unit="in",width=3,height=3,bg = "transparent") 


levels(df1$num_terr)

m1 = lmer(node_oppsex_edge_strength ~ num_terr*day + (1|name), data = filter(df1, day %in% 1:5, num_terr %in% 1:2)) 
m2 = lmer(node_oppsex_edge_strength ~ num_terr*day + (1|name), data = filter(df1, day %in% 6:20, num_terr %in% 1:2)) 
m3=lmer(node_oppsex_edge_strength ~ num_terr+day + (day|name), data = filter(df1, day %in% 6:20, num_terr %in% 1:2)) 
m4 = lmer(node_oppsex_edge_strength ~ num_terr*day + (1|name), data = filter(df1, day %in% 6:20, num_terr %in% 0:2)) 

AIC(m1,m2,m3,m4)
summary(m3)
write.table(summary(m3)$coef, "clipboard-16384", sep="\t", row.names=TRUE, col.names = TRUE)


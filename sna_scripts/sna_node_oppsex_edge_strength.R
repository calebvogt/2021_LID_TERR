## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(ggpubr)
library(lme4) 
library(lmerTest)

df <- read.csv("data/ALLTRIAL_SNA_node_stats.csv")

# 1v2 night 1-5 ----------------------------------------------------------
df1 <- df %>% 
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
    xlab("Night") +
    ylab("Node edge strength with females") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.8,0.8), 
          legend.spacing = unit(0.1, "cm"))
          # legend.position = "")
ggsave(file="output/sna_node_oppsex_edge_strength.svg",device="svg",unit="in",width=4,height=4,bg = "transparent") 

m1 = lmer(node_edge_strength ~ num_terr*day + (1|name), data = subset(df1, num_terr == 1|2)) 
summary(m1)
AIC(m1)

m2 = lmer(node_degree_centrality ~ genotype*sex*log(day) + (1|trial) + (log(day)|name), data = df) 
summary(m2)

m3 = lmer(node_degree_centrality ~ genotype*sex*day + (1|trial) + (day|name), data = df) 
summary(m3)

AIC(m1,m2,m3)
write.table(summary(m2)$coef, "clipboard-16384", sep="\t", row.names=TRUE, col.names = TRUE)

# # daily models: change for days 1 through 10. 
# # i = 4
# for(i in 1:10){
#   d = lmer(node_degree_centrality ~ genotype*sex + (1|trial), data = subset(df, day == i)) 
#   #get daily contrast estimates
#   contrasts <- emmeans(d, pairwise ~ genotype*sex)
#   print(i)
#   print(contrasts$contrasts)
#   # anova(d)
#   # summary(d)
# }



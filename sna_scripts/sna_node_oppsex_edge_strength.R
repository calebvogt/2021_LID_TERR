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
  filter(day %in% 1:5 & sex == "M") %>% 
  filter(num_terr %in% 1:2)

df2 <- df1 %>% 
  group_by(num_terr, day) %>%
  summarise(mean = mean(node_oppsex_edge_strength), sd = sd(node_oppsex_edge_strength), count = n(), sem = (sd/(sqrt(count))))

(p <- ggplot(df2, aes(x=day, y=mean, color=as.factor(num_terr))) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(0.8,5.3), breaks = seq(1, 10, by = 1)) +
    # scale_y_continuous(breaks = seq(1, 20, by = 1)) +
    scale_color_manual(breaks = c("0","1", "2"),
                       values=c("black", "red", "blue")) +
    theme_classic() +
    xlab("Night") +
    ylab("Node opposite sex edge strength") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.2,0.9))
          # legend.position = "")
)
# ggsave(p, filename = "output/sna_node_edge_strength.png", device = "png", bg = "white")
# ggsave(p, filename = "output/sna_node_edge_strength.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


# 0v1v2 night 6-20 ----------------------------------------------------------
df1 <- df %>% 
  filter(day %in% 6:20 & sex == "M")

df2 <- df1 %>% 
  group_by(num_terr, day) %>%
  summarise(mean = mean(node_oppsex_edge_strength), sd = sd(node_oppsex_edge_strength), count = n(), sem = (sd/(sqrt(count))))

(p <- ggplot(df2, aes(x=day, y=mean, color=as.factor(num_terr))) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(5.8,20.3), breaks = seq(6, 20, by = 1)) +
    # scale_y_continuous(breaks = seq(1, 20, by = 1)) +
    scale_color_manual(breaks = c("0","1", "2"),
        values=c("black", "red", "blue")) +
    theme_classic() +
    xlab("Night") +
    ylab("Node opposite sex edge strength") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.70,0.9))
  # legend.position = "")
)
# ggsave(p, filename = "output/sna_node_degree_centrality.png", device = "png", bg = "white")
# ggsave(p, filename = "output/sna_node_degree_centrality.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

# STATS
m1 = lmer(node_edge_strength ~ num_terr*day + (1|name), data = df1) 
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



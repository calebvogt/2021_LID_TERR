## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(ggpubr)
library(lme4)
library(lmerTest)
library(ggpubr)

# load data ---------------------------------------------------------------
meta <- read.csv("data/metadata.csv")
gbi <- read.csv("data/ALLTRIAL_MOVEBOUT_GBI_summary.csv")
gbi$name <- gsub("C57-M-","",gbi$name)
gbi$name <- gsub("C57-F-","",gbi$name)
data <- gbi %>% 
  merge(., meta[,c("name", "num_terr")], by="name")

# time with females, nights 1-5, 0v1v2 males --------------------------------------
df <- data %>% 
  group_by(name,sex, code, day) %>% 
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]),
            social_s=sum(duration_s[mf_sum>1]),
            mf_s=sum(duration_s[m_sum>0&f_sum>0]),
            mm_s=sum(duration_s[m_sum>1&f_sum==0]),
            ff_s=sum(duration_s[m_sum==0&f_sum>1])) %>% 
  mutate(samesex_s=ifelse(sex=="M", mm_s,ff_s), 
         oppsex_s=mf_s, .after=social_s) %>% 
  merge(., meta[,c("name", "num_terr")], by = "name") 

df1 <- df %>% 
  filter(day %in% 1:5 & sex == "M") %>% 
  filter(num_terr %in% 1:2) %>% 
  group_by(num_terr, day) %>%
  summarise(mean = mean(oppsex_s/3600), sd = sd(oppsex_s/3600), count = n(), sem = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day, y=mean, color=as.factor(num_terr))) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(0.8,5.3), breaks = seq(1, 5, by = 1)) +
    # scale_y_continuous(breaks = seq(1, 20, by = 1)) +
    scale_color_manual(breaks = c("0","1", "2"),
                       values=c("black", "red", "blue")) +
    theme_classic() +
    xlab("Night") +
    ylab("Time with females (h)") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.1,0.9))
  # legend.position = "")
)
# ggsave(p, filename = "output/sna_node_edge_strength.png", device = "png", bg = "white")
# ggsave(p, filename = "output/sna_node_edge_strength.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


# time with females, nights 6-20, 0v1v2 males --------------------------------------
df <- data %>% 
  group_by(name,sex, code, day) %>% 
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]),
            social_s=sum(duration_s[mf_sum>1]),
            mf_s=sum(duration_s[m_sum>0&f_sum>0]),
            mm_s=sum(duration_s[m_sum>1&f_sum==0]),
            ff_s=sum(duration_s[m_sum==0&f_sum>1])) %>% 
  mutate(samesex_s=ifelse(sex=="M", mm_s,ff_s), 
         oppsex_s=mf_s, .after=social_s) %>% 
  merge(., meta[,c("name", "num_terr")], by = "name") 

df1 <- df %>% 
  filter(day %in% 6:20 & sex == "M") %>% 
  group_by(num_terr, day) %>%
  summarise(mean = mean(oppsex_s/3600), sd = sd(oppsex_s/3600), count = n(), sem = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day, y=mean, color=as.factor(num_terr))) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(5.8,20.3), breaks = seq(6, 20, by = 1)) +
    # scale_y_continuous(breaks = seq(1, 20, by = 1)) +
    scale_color_manual(breaks = c("0","1", "2"),
                       values=c("black", "red", "blue")) +
    theme_classic() +
    xlab("Night") +
    ylab("Time with females (h)") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.50,0.9))
  # legend.position = "")
)
# ggsave(p, filename = "output/sna_node_edge_strength.png", device = "png", bg = "white")
# ggsave(p, filename = "output/sna_node_edge_strength.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")







# prop. obs.  time with females, nights 1-5, 0v1v2 males --------------------------------------
df <- data %>% 
  group_by(name,sex, code, day) %>% 
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]),
            social_s=sum(duration_s[mf_sum>1]),
            mf_s=sum(duration_s[m_sum>0&f_sum>0]),
            mm_s=sum(duration_s[m_sum>1&f_sum==0]),
            ff_s=sum(duration_s[m_sum==0&f_sum>1])) %>% 
  mutate(samesex_s=ifelse(sex=="M", mm_s,ff_s), 
         oppsex_s=mf_s, .after=social_s) %>% 
  merge(., meta[,c("name", "num_terr")], by = "name") 

df1 <- df %>% 
  filter(day %in% 1:5 & sex == "M") %>% 
  filter(num_terr %in% 1:2) %>% 
  group_by(num_terr, day) %>%
  summarise(mean = mean(oppsex_s/sum_s), sd = sd(oppsex_s/sum_s), count = n(), sem = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day, y=mean, color=as.factor(num_terr))) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(0.8,5.3), breaks = seq(1, 5, by = 1)) +
    # scale_y_continuous(breaks = seq(1, 20, by = 1)) +
    scale_color_manual(breaks = c("0","1", "2"),
                       values=c("black", "red", "blue")) +
    theme_classic() +
    xlab("Night") +
    ylab("Prop. observation time spent with females") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.15,0.9))
  # legend.position = "")
)
# ggsave(p, filename = "output/sna_node_edge_strength.png", device = "png", bg = "white")
# ggsave(p, filename = "output/sna_node_edge_strength.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


# prop. obs. time with females, nights 6-20, 0v1v2 males --------------------------------------
df <- data %>% 
  group_by(name,sex, code, day) %>% 
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]),
            social_s=sum(duration_s[mf_sum>1]),
            mf_s=sum(duration_s[m_sum>0&f_sum>0]),
            mm_s=sum(duration_s[m_sum>1&f_sum==0]),
            ff_s=sum(duration_s[m_sum==0&f_sum>1])) %>% 
  mutate(samesex_s=ifelse(sex=="M", mm_s,ff_s), 
         oppsex_s=mf_s, .after=social_s) %>% 
  merge(., meta[,c("name", "num_terr")], by = "name") 

df1 <- df %>% 
  filter(day %in% 6:20 & sex == "M") %>% 
  group_by(num_terr, day) %>%
  summarise(mean = mean(oppsex_s/sum_s), sd = sd(oppsex_s/sum_s), count = n(), sem = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day, y=mean, color=as.factor(num_terr))) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(5.8,20.3), breaks = seq(6, 20, by = 1)) +
    # scale_y_continuous(breaks = seq(1, 20, by = 1)) +
    scale_color_manual(breaks = c("0","1", "2"),
                       values=c("black", "red", "blue")) +
    theme_classic() +
    xlab("Night") +
    ylab("Prop. observation time spent with females") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.50,0.9))
  # legend.position = "")
)
# ggsave(p, filename = "output/sna_node_edge_strength.png", device = "png", bg = "white")
# ggsave(p, filename = "output/sna_node_edge_strength.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")








# prop. obs. time with males, nights 6-20, 0v1v2 males --------------------------------------
df <- data %>% 
  group_by(name,sex, code, day) %>% 
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]),
            social_s=sum(duration_s[mf_sum>1]),
            mf_s=sum(duration_s[m_sum>0&f_sum>0]),
            mm_s=sum(duration_s[m_sum>1&f_sum==0]),
            ff_s=sum(duration_s[m_sum==0&f_sum>1])) %>% 
  mutate(samesex_s=ifelse(sex=="M", mm_s,ff_s), 
         oppsex_s=mf_s, .after=social_s) %>% 
  merge(., meta[,c("name", "num_terr")], by = "name") 

df1 <- df %>% 
  filter(day %in% 6:20 & sex == "M") %>% 
  group_by(num_terr, day) %>%
  summarise(mean = mean(samesex_s/sum_s), sd = sd(samesex_s/sum_s), count = n(), sem = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day, y=mean, color=as.factor(num_terr))) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(5.8,20.3), breaks = seq(6, 20, by = 1)) +
    # scale_y_continuous(breaks = seq(1, 20, by = 1)) +
    scale_color_manual(breaks = c("0","1", "2"),
                       values=c("black", "red", "blue")) +
    theme_classic() +
    xlab("Night") +
    ylab("Prop. observation time spent with males") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.50,0.9))
  # legend.position = "")
)
# ggsave(p, filename = "output/sna_node_edge_strength.png", device = "png", bg = "white")
# ggsave(p, filename = "output/sna_node_edge_strength.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")











# prop. obs. time alone, nights 6-20, 0v1v2 males --------------------------------------
df <- data %>% 
  group_by(name,sex, code, day) %>% 
  summarize(sum_s=sum(duration_s),
            alone_s=sum(duration_s[mf_sum==1]),
            social_s=sum(duration_s[mf_sum>1]),
            mf_s=sum(duration_s[m_sum>0&f_sum>0]),
            mm_s=sum(duration_s[m_sum>1&f_sum==0]),
            ff_s=sum(duration_s[m_sum==0&f_sum>1])) %>% 
  mutate(samesex_s=ifelse(sex=="M", mm_s,ff_s), 
         oppsex_s=mf_s, .after=social_s) %>% 
  merge(., meta[,c("name", "num_terr")], by = "name") 

df1 <- df %>% 
  filter(day %in% 6:20 & sex == "M") %>% 
  group_by(num_terr, day) %>%
  summarise(mean = mean(alone_s/sum_s), sd = sd(alone_s/sum_s), count = n(), sem = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=day, y=mean, color=as.factor(num_terr))) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.2) +
    scale_x_continuous(limits = c(5.8,20.3), breaks = seq(6, 20, by = 1)) +
    # scale_y_continuous(breaks = seq(1, 20, by = 1)) +
    scale_color_manual(breaks = c("0","1", "2"),
                       values=c("black", "red", "blue")) +
    theme_classic() +
    xlab("Night") +
    ylab("Prop. observation time spent alone") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=6),
          legend.background = element_rect(fill='transparent'),
          legend.position = c(0.50,0.9))
  # legend.position = "")
)
# ggsave(p, filename = "output/sna_node_edge_strength.png", device = "png", bg = "white")
# ggsave(p, filename = "output/sna_node_edge_strength.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")










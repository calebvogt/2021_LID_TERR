## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
df0 <- rfid_data %>% 
  filter(noon_day %in% 1:20) %>% 
  filter(!(antenna %in% 13:16)) %>% 
  mutate(sex_treatment = paste(sex, treatment, sep = "-"))

df <- df0 %>% 
  select(trial, sex_treatment, sex,strain, name, noon_day, antenna) %>%
  group_by(trial, sex_treatment, sex,strain, name, antenna) %>% 
  distinct(antenna, .keep_all = TRUE) %>% # get distinct antennas ever visited, keep associated day it was visited
  group_by(trial, sex_treatment,sex,strain, name, noon_day) %>% #regroup by day
  tally() %>% #tally unique antennas visited per day
  mutate(csum_novel_antennas = cumsum(n)) %>%  
  complete(noon_day = full_seq(1:20, period = 1)) %>% #fill in missing day rows for each mouse
  arrange(name, noon_day) %>% 
  fill(csum_novel_antennas) %>% ## fill cumulative sum data from last observed day
  select(trial, sex_treatment, sex,strain, name, noon_day,csum_novel_antennas) %>% 
  na.omit(csum_novel_antennas)
  
df1 <- df %>% 
  group_by(sex_treatment, noon_day) %>% 
  summarise(mean = mean(csum_novel_antennas), sd = sd(csum_novel_antennas),count = n(), se = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x=noon_day, y=mean, color = sex_treatment)) + 
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,20.3), breaks = seq(1,21,by=1)) +
    scale_y_continuous(limits =c(1,12), breaks = seq(1,12,by=1)) +
    # scale_color_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                       # values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    xlab("Day") +
    ylab("Cumulative novel zones visited") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.position = "right") +
    guides(color=guide_legend(override.aes=list(fill=NA)))
)
ggsave(p, filename = "output/rfid_csum_novel_zones_visited.png", device = "png", bg = "white")
# ggsave(p, filename = "output/rfid_data_cumulative_novel_antennas_visited.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


# STATS
df[sapply(df, is.character)] <- lapply(df[sapply(df, is.character)],as.factor)
df$csum_novel_antennas <- as.numeric(df$csum_novel_antennas)

mod2 = lmer(csum_novel_antennas ~ sex*strain + (1|trial), data = subset(df, noon_day == 10))
anova(mod2)
summary(mod2)
write.table(summary(mod2)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

em.mod2 = emmeans(mod2, pairwise ~ sex*strain)
write.table(em.mod2$emmeans, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
write.table(em.mod2$contrasts, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)


#logistic regression model of differences in all antennas visited by strain*sex interaction
df2 <- df
df2$every_antenna <- ifelse(df2$csum_novel_antennas == 8, 1, 0)
df3 <- subset(df2, noon_day == 10)
mod3 = glmer(every_antenna ~ sex_treatment + (1|trial), data = subset(df2, noon_day == 10), family ="binomial") #family binomial because response variable is 0 or 1. 
mod3.1 = glmer(every_antenna ~ (1|trial), data = subset(df2, noon_day == 10), family ="binomial") #family binomial because response variable is 0 or 1. 
aov(mod3, mod3.1)
anova(mod3,mod3.1)
anova(mod3)
summary(mod3)
write.table(summary(mod3)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)




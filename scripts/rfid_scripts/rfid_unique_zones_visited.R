## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)
library(ggeffects)

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
df0 <- rfid_data %>% 
  filter(noon_day %in% 1:20) %>% 
  filter(!(antenna %in% 13:16)) %>% 
  mutate(sex_treatment = paste(sex, treatment, sep = "_"))
  
df <- df0 %>% # 
  mutate(sex_treatment = paste(sex, treatment, sep = "-")) %>% 
  select(trial, sex_treatment, sex, strain, name, noon_day, antenna) %>% 
  group_by(trial, sex_treatment, sex,strain, name, noon_day, antenna) %>% 
  tally() %>%  # get number of visits to each antenna
  group_by(trial, sex_treatment, sex, strain,name, noon_day, .drop = FALSE) %>%
  tally() %>% # get number of unique antenna visits
  complete(noon_day = 1:max(noon_day), fill = list(n = 0)) %>% # fill in days where mouse doesnt appear with 0s
  dplyr::rename(unique_antennas_visited = n) %>% ## add filtering after this if necessary if animals are known to have died. 
  filter(!(sex_treatment == "M-late" & noon_day %in% 1:5)) %>% 
  filter(!(sex_treatment == "F-late" & noon_day %in% 1:5))


# All sex and strain groups -----------------------------------------------
df1 <- df %>% 
  group_by(sex_treatment, noon_day) %>%
  summarise(mean = mean(unique_antennas_visited), sd = sd(unique_antennas_visited),count = n(), se = (sd/(sqrt(count))))

(p <- ggplot(df1, aes(x = noon_day, y = mean, color = sex_treatment)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.2) +
    scale_x_continuous(limits = c(0.8,20.3), breaks = seq(1, 20, by = 1)) +
    # scale_y_continuous(limits = c(1,8), breaks = seq(1, 8, by = 1)) +
    theme_classic() +
    xlab("Day") +
    ylab("Unique zones visited") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.background = element_rect(fill='transparent'),
          # legend.position = c(0.21,0.8))
          legend.position = "right")
)
ggsave(p, filename = "output/rfid_unique_zones_visited.png", device = "png", bg = "white")
# ggsave(p, filename = "output/rfid_unique_zone_antennas_visited.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")


# STATS
df$trial <- as.factor(df$trial)
df$sex <- as.factor(df$sex)
df$strain <- as.factor(df$strain)
df$noon_day <- as.numeric(df$noon_day)
df$c57F <- ifelse(df$strain=="C57"&df$sex=="F",1,0)

m1 = lmer(unique_antennas_visited ~ strain*sex*log(noon_day) + (1|trial) + (1+log(noon_day)|name), data = df) 
m2 = lmer(unique_antennas_visited ~ strain*sex+log(noon_day) + (1|trial) + (1+log(noon_day)|name), data = df)
# Removing C57s demonstrates they are solely responsible for effect. 
# m3 = lmer(unique_antennas_visited ~ sex_treatment*noon_day + (1|name), data = subset(df, df$sex_treatment != "C57-F"))
# m4 = lmer(unique_antennas_visited ~ c57F*log(noon_day) + (1|trial) + (1|name), data = df) # keeps C57 females in model, but evaluates them separately. 
AIC(m1,m2,m3,m4)
qqnorm(resid(m1))
qqline(resid(m1))
summary(m1)
anova(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
write.table(summary(m4)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

emmeans(m1, pairwise ~ strain*sex*log(noon_day)) #tukey

ggpredict(mod1, terms = c("strain", "sex", "noon_day"), type = "fe") %>% 
  plot()

# daily models: change for days 1 through 10. 
d = lmer(unique_antennas_visited ~ strain*sex + (1|trial), data = subset(df, noon_day == 10)) 
em.d <- emmeans(d, pairwise ~ strain*sex) # tukey adjusted
write.table(em.d$contrasts, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)



# summary stats of male only maximum average daily antennas visited --------------------------------------------------------------
df1 <- df %>% 
  filter(sex=="M")

df2 <- df1 %>% 
  group_by(name) %>%
  summarise(mean_antenna = mean(unique_antennas_visited), sd_antenna = sd(unique_antennas_visited),count = n(), se_antenna = (sd_antenna/(sqrt(count))))

median(df2$mean_antenna)
range(df2$mean_antenna)

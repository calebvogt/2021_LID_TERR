## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(lme4)
library(lmerTest)
library(emmeans)

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
df0 <- rfid_data %>% 
  filter(noon_day %in% 1:20) %>% 
  mutate(sex_treatment = paste(sex, treatment, sep = "-"))

ids <- unique(df$name)
data_list <- list()
aa = ids[60]
for(aa in ids[1:length(ids)]){
  df1 <- df0 %>% 
    filter(name == aa) 
  
  # delete consecutive repeat antenna hits, only keep rows where antennas change. 
  df2 <- as.data.table(df1)[, .SD[1], by = rleid(df1$antenna)] 
  df2$dist <- NA
  n <- nrow(df2)
  if(n==1){
    df2$dist[1] <- 0
    data_list[[aa]] <- df2
  } else{
    df2$dist[2:n] <- sqrt((df2$zone_x[2:n] - df2$zone_x[1:n-1])^2 + (df2$zone_y[2:n] - df2$zone_y[1:n-1])^2)
    df2$dist[1] <- 0
    data_list[[aa]] <- df2
  }
}
df3 <- do.call("rbind", data_list)

df4 <- df3 %>% 
  group_by(trial, sex_treatment,strain, sex, name, noon_day) %>% 
  tally(sum(dist)) %>% 
  complete(noon_day = 1:20, fill = list(n = 0)) %>%  # fill in days where mouse doesnt appear with 0s.
  dplyr::rename(dist = n) %>% 
  filter(!(sex_treatment == "M-late" & noon_day %in% 1:5)) %>% 
  filter(!(sex_treatment == "F-late" & noon_day %in% 1:5))
  

df5 <- df4 %>% 
  group_by(sex_treatment, noon_day) %>%
  summarise(mean_n = mean(dist), sd_n = sd(dist),count = n(), se_n = (sd_n/(sqrt(count))))
  
(p <- ggplot(df5, aes(x = noon_day, y = mean_n, color = sex_treatment)) +
    geom_line(size = 0.75) + 
    geom_point(size = 1.5) +
    geom_errorbar(aes(ymin = mean_n - se_n, ymax = mean_n + se_n), width = 0.2) +
    scale_x_continuous(limits = c(0.8,20.3), breaks = seq(1, 20, by = 1)) +
    # scale_y_continuous(breaks = seq(0, 150, by = 20)) +
    theme_classic() +
    xlab("Day") +
    ylab("Minimum distance travelled (m)") +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"), 
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"), 
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.title = element_blank(),
          legend.text = element_text(size=8),
          legend.position = "right" )
)
ggsave(p, filename = "output/rfid_distance_travelled.png", device = "png", bg = "white")
# ggsave(p, filename = "output/rfid_data_distance_travelled.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

# STATS  
df7 <- df4 %>% 
  group_by(trial, strain, sex, name, noon_day) %>% 
  tally(sum(dist)) %>% 
  group_by(trial, strain, sex, name, noon_day) %>% 
  dplyr::rename(sum_dist = n)

# write.table(df7, "clipboard-16384", sep="\t", row.names=FALSE, col.names = TRUE)
# df7 <- read_excel("Figure_Data.xlsx", sheet = "Fig2a")

df7$trial <- as.factor(df7$trial)
df7$strain <- as.factor(df7$strain)
df7$sex <- as.factor(df7$sex)
df7$name <- as.factor(df7$name)
df7$c57_F<-ifelse(df7$strain!="C57"|df7$sex!="F","other","pC57_F")

mod1 = lmer(sum_dist ~ strain*sex*log(noon_day) + (1|trial) + (log(noon_day)|name), data = df7) 
mod2 = lmer(sum_dist ~ strain*sex*log(noon_day) + (log(noon_day)|name), data = df7) 
mod3 = lmer(sum_dist ~ strain*sex*log(noon_day) + (1|trial) + (1|name), data = df7) 
mod4 = lmer(sum_dist ~ strain+sex+log(noon_day) + (1|trial) + (log(noon_day)|name), data = df7) 
AIC(mod1,mod2, mod3)
qqnorm(resid(mod1))
anova(mod1)
summary(mod1)
emmeans(mod1, pairwise ~ strain*sex*log(noon_day))
write.table(summary(mod1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

# daily models: change for days 1 through 10 and report
d = lmer(sum_dist ~ strain*sex + (1|trial), data = subset(df7, noon_day == 10)) 
em.d <- emmeans(d, pairwise ~ strain*sex) #tukey adjusted
write.table(em.d$contrasts, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)
  



## Total distance travelled
df6a <- df4 %>% 
  group_by(trial, sex_treatment, name) %>% 
  tally(sum(dist)) %>% 
  dplyr::rename(total_dist = n)

# write.table(df6a, "clipboard-16384", sep="\t", row.names=FALSE, col.names = TRUE)
df6a <- read_excel("Figure_Data.xlsx", sheet = "FigS2a")

(p <- ggplot(df6a, aes(x=sex_treatment, y=total_dist, fill = sex_treatment)) + 
    geom_violin(width=1, alpha = 1) +
    geom_boxplot(width=0.1, color = "black", alpha=0.5, size = 0.5) +
    scale_y_continuous(limits = c(0,2000)) +
    scale_fill_manual(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                      values=c("sienna1", "sienna", "skyblue", "skyblue4")) +
    # scale_fill_discrete(breaks = c("C57-F", "C57-M", "NYOB-F", "NYOB-M"),
                        # values=c("sienna1", "sienna", "skyblue", "skyblue4"),
                        # labels=c('C57 Female', 'C57 Male', 'Outbred Female', 'Outbred Male')) +
    xlab("") +
    ylab("Total estimated dist. travelled (m)") +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 8),
          axis.title.x = element_text(color = "black", size = 8, face = "bold"),
          axis.text.y = element_text(color = "black", size = 8),
          axis.title.y = element_text(color = "black", size = 8, face = "bold"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent"), 
          legend.position = "none")
)
ggsave(p, filename = "output/rfid_data_distance_travelled_total.png", device = "png", bg = "white")
ggsave(p, filename = "output/rfid_data_distance_travelled_total.svg", device = "svg", width=2.5, height=2.15, bg = "transparent")

m1 <- lmer(total_dist ~ sex_treatment + (1|trial), data = df6a)
summary(m1)
anova(m1)
write.table(summary(m1)$coef, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

em.d <- emmeans(m1, pairwise ~ sex_treatment, adjust = 'none') 
write.table(em.d$contrasts, "clipboard", sep="\t", row.names=TRUE, col.names = TRUE)

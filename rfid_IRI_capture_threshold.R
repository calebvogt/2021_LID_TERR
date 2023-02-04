## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

# for each mouse, get groups of consecutive rfid reads in the same zone. 
# Get the interread intervals for each group, then combine all of these inter-read intervals. 
# Note that this will group consecutive rfid reads seperated by many days. 
# Then get the 99% condfidence interval across every observed interread interval for every read ever, across all mice.
## improved the older method in that it counted interread intervals purely on a zone basis, in that it did not end bouts
## when the mice were known to have visited another zone. Thus, would have over inflated the interread interval distribution.

library(tidyverse)
library(data.table)
library(readxl)

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
df <- rfid_data %>% 
  filter(noon_day %in% 1:20) ## get IRI only for reads between days 1-20

ids <- unique(df$name)
id_list <- list()
flag <- 1
aa = ids[1]
for(aa in ids[1:length(ids)]){
  print(paste("Processing ", aa, " mouse ",flag, " out of ", length(ids), sep=''))
  df1 <- df %>% 
    filter(name == aa) %>% 
    arrange(field_time) %>%
    mutate(visit_group = NA)
  df1$visit_group[1] <- 1
  
  if(nrow(df1)>1){
    # label consecutive rfid read groups (antenna visit groups)
    bb = 2
    for(bb in 2:nrow(df1)) {
      if(df1$antenna[bb] == df1$antenna[bb-1]){
        df1$visit_group[bb] <- df1$visit_group[bb-1] 
      } else{
        df1$visit_group[bb] <- df1$visit_group[bb-1] + 1 
      }
    }
    
    # find interread interval for the visit group and put into list. 
    visit_group <- unique(df1$visit_group)
    visit_group_list <- list()
    cc = 3
    for(cc in visit_group[1:length(visit_group)]){
      df2 <- df1 %>% 
        filter(visit_group == cc) %>% 
        mutate(diff = field_time - lag(field_time), diff_secs = as.numeric(diff, units = 'secs')) %>% 
        filter(diff_secs > 0)
      
      if(nrow(df2) > 0){
        visit_group_list[[cc]] <- df2
      }
    }
  } else{}
  
  id_list[[aa]] <- do.call("rbind", visit_group_list)
  flag <- flag+1
}
df3 <- do.call("rbind", id_list)
df4 <- df3[!is.na(df3$diff_secs),]
df5 <- df4 %>% 
  filter(!(antenna %in% 13:16)) ## remove water tower visits and all diff secs for the water tower visit groups

summary(df5$diff_secs)
sort(df3$diff_secs)[0.95*length(df3$diff_secs)]  ## 163s
sort(df3$diff_secs)[0.99*length(df3$diff_secs)]  ## 1061s 

## Graph
png(filename = "output/rfid_IRI_interval_threshold.png")
sdat <- summary(df3$diff_secs)
summStr <- paste(names(sdat), format(sdat, digits = 3), collapse = "; ")
op <- par(mar = c(7,4,4,2) + 0.1)
hist(df3$diff_secs, 
     xlim = c(0, 500),
     # log = "y",
     breaks = 30000,
     main = "",
     xlab = "Within zone inter-read interval (s)"
)
abline(v=c(163), col=c("blue"), lty=c(1,1,1), lwd=c(3,3,3))
title(sub = summStr, line = 5.5)
dev.off()



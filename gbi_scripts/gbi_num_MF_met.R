## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)
library(readxl)
library(asnipe)
library(igraph)

# load  & clean data ---------------------------------------------------------------
dd <- paste(getwd(), "data/", sep = "/")
filenames <- list.files(dd, pattern = "*MOVEBOUT_GBI.csv", full.names = T)
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


# Number of females met, night 1:5, 1v2 males -----------------------------------------------
df <- data %>% 
  filter(day %in% 1:5 & sex == "M")

df1 <- df %>% 
  group_by(num_terr, day) %>%
  summarise(mean = mean(f_met), 
            sd = sd(f_met), 
            count = n(), 
            sem = (sd/(sqrt(count))))

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
    ylab("# of females met") +
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


# Number of females met, night 6:20, 0v1v2 males -----------------------------------------------
df <- data %>% 
  filter(day %in% 6:20 & sex == "M")

df1 <- df %>% 
  group_by(num_terr, day) %>%
  summarise(mean = mean(f_met), 
            sd = sd(f_met), 
            count = n(), 
            sem = (sd/(sqrt(count))))

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
    ylab("# of females met") +
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


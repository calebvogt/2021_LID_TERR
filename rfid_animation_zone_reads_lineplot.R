## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(readxl)
library(gganimate)
library(gifski)
library(av)

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))
df <- rfid_data %>% 
  filter(noon_day %in% 1:20) %>% 
  filter(!(antenna %in% 13:16)) %>% 
  mutate(sex_treatment = paste(sex, treatment, sep = "_"))

ids <- unique(df$name)
i=ids[1]
for(i in ids[1:length(ids)]) {
  # SAVE CURRENT MOUSE AS VARIABLE
  current_mouse <- print(i)
  move_df <- filter(df, name == current_mouse)
  current_mouse_info <- paste(unique(move_df$trial), unique(move_df$sex), unique(move_df$treatment),unique(move_df$name), sep = "_")
  
  (p2 <- ggplot(move_df, aes(field_time, antenna, color = factor(name))) + 
    geom_line(na.rm=TRUE, color="red", size=1) +
    geom_point() +
    ggtitle(paste(current_mouse_info)) +
    scale_color_viridis_d() +
    labs(x = "Date", y = "antenna") +
    scale_y_continuous(breaks = seq(1,8,1), limits=c(1,8)) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust = 1),
          axis.line = element_line(color = "black", size = 1, linetype = "solid"))
  )
  
  # RENDER THE ANIMATION
  p2a <- p2 +
    geom_point(aes(group=seq_along(field_time))) +
    transition_reveal(field_time)
  
  # SAVE THE ANIMATION
  animate(p2a, duration = 10, fps = 10, renderer = gifski_renderer())
  anim_save(paste("output/gifs/","rfid_",current_mouse_info,"_zone_reads_lineplot.gif", sep=""))
}

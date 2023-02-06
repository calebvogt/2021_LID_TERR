## Created by Caleb C. Vogt, PhD Candidate @ Cornell University

library(tidyverse)
library(data.table)
library(readxl)
library(gganimate)
library(gifski)
library(av)

rfid_data <- as.data.frame(fread("data/ALLTRIAL_RFID_DATA.csv", stringsAsFactors = FALSE, fill = TRUE, header = TRUE, check.names = TRUE))

#get data of xy boxes points
box <- unique(rfid_data[c("antenna","zone_x","zone_y")])
box <- box %>% 
  mutate(x1=zone_x+1, x2=zone_x-1, y1=zone_y+1, y2=zone_y-1) %>% 
  mutate(antenna_type = "zone") %>% 
  mutate(antenna_type=replace(antenna_type, antenna > 12 & antenna < 17, "water"))

df <- rfid_data %>% 
  filter(noon_day %in% (1:20)) %>% 
  mutate(ant_rleid = antenna)

ids <- unique(df$name)
ids
i=ids[1]
for(i in ids[1:length(ids)]) {
  # SAVE CURRENT MOUSE AS VARIABLE
  current_mouse <- print(i)
  df1 <- subset(df, name == current_mouse) # ONLY KEEP OBSERVATIONS FROM CURRENT_MOUSE, REMOVE OTHER MICE
  df2 <- as.data.table(df1)[, .SD[1], by = rleid(df1$ant_rleid)] ## take only the transitions between novel antennas 
  df3 <- df2 %>% 
    complete(noon_day = 1:20)
  
  move_df <- df3 %>% 
    mutate(time = 1:nrow(df3)) %>% 
    mutate(zone_x_j = jitter(df3$zone_x), 
           zone_y_j = jitter(df3$zone_y))
  
  # SAVE CURRENT SUBJECT AS VARIABLE
  current_mouse_info <- paste(unique(move_df$trial), unique(move_df$sex), unique(move_df$treatment),unique(move_df$name), sep = "_")
  current_mouse_info <- current_mouse_info[current_mouse_info != "NA_NA_NA_NA"]
  print(current_mouse_info)
  
  (p <- ggplot() +
      geom_rect(data=box, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=antenna_type), color="black", alpha=0.5) +
      # geom_text(data=box, aes(x=x1+(x2-x1)/2, y=y1+(y2-y1)/2, label=antenna), size=4, vjust = -2.2) +
      xlim(0,15.24) +
      ylim(0,38.1) +
      scale_color_viridis_d() + ## change this to scale color manual and change zones to green, water to blue circles
      labs(x = "none", y = "none") +
      theme(plot.background = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, size = 3),
            panel.background = element_blank(),
            axis.line = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            axis.ticks = element_blank(),
            axis.text = element_blank(),
            legend.title = element_blank(),
            legend.text = element_blank(),
            legend.position = "none",
            # legend.text = element_text(size = 15),
            plot.title = element_text(size = 20))
  )
  
  #add animation and titles. ## improvements. pause on days. 
  p2a <- p +
    geom_point(data=move_df, aes(x=zone_x_j, y=zone_y_j, group=seq_along(as.factor(time))), size = 4) + 
    geom_line(data=move_df, aes(x=zone_x_j, y=zone_y_j, group=1), color="red", size=1, alpha = 0.3) +
    labs(title = paste(current_mouse_info, "Day: {move_df$noon_day[which.min(abs(move_df$time-frame_along))]}"), 
         subtitle = paste("Total Zone Transitions: {frame_along}")) +
    transition_reveal(time)
  # 
  # ## testing. Make transitions between days even in time.... not quite sure yet how to do this. 
  # (p2a <- p +
  #   geom_point(data=move_df, aes(x=zone_x_j, y=zone_y_j, group=seq_along(as.factor(time))), size = 4) + 
  #   geom_line(data=move_df, aes(x=zone_x_j, y=zone_y_j, group=1), color="red", size=1, alpha = 0.3) +
  #   labs(title = paste(current_mouse_info, "Day: {move_df$noon_day[which.min(abs(move_df$time-frame_along))]}"), 
  #        subtitle = paste("Total Zone Transitions: {frame_along}")) +
  #   transition_states(noon_day, state_length = 1, transition_length = 1)
  # )
  
  # SAVE GIF
  animate(p2a, width = 400, height = 600, renderer = gifski_renderer())
  anim_save(paste("output/gifs/","rfid_",current_mouse_info, "_paddock_movement.gif", sep=""))
}

## 1_create_ALLTRIAL_RFID_DATA
## Caleb C. Vogt, Cornell University

library(tidyverse)
library(readxl)
library(data.table)
library(plyr)

setwd(paste(getwd(), "data", sep = "/"))
wd <- getwd()
metadata <- read.csv("LID_TERR_2021_metadata_ID.csv",header=T,stringsAsFactors = T)
metadata$tag_1 <- as.numeric(format(metadata$tag_1, digits = 15))
metadata$tag_2 <- as.numeric(format(metadata$tag_2, digits = 15))
weather <- read.csv("LID_TERR_2021_metadata_weather.csv",header=T)
weather$weather_time <- as.POSIXct(weather$date, format = "%m/%d/%Y %H:%M")
folders <- list.dirs(wd, recursive=FALSE, full.names = TRUE)
folders_short <- list.dirs(wd, recursive=FALSE, full.names = FALSE)

options(digits=15)
i=folders[1]
flag=1
all_trial_list <- list()
for (i in folders[1:length(folders)]) {
  setwd(i)
  trial_var <- folders_short[flag]
  txt <- list.files(pattern = "*.txt")
  dfall <- do.call(bind_rows, lapply(txt, function(x) read_table(file = x, col_names = c("Scan Date","Scan Time","Download Date","Download Time","Reader ID", "Antenna ID", "HEX Tag ID", "DEC Tag ID","Signal,mV", "Is Duplicate"),
                                                                 col_types = NULL,
                                                                 locale = default_locale(), na = "NA", skip = 6, n_max = Inf,
                                                                 progress = show_progress(), comment = )))
  dfall <- dfall[ , c("Scan Date", "Scan Time", "Reader ID", "Antenna ID", "DEC Tag ID")] 
  colnames(dfall) <- c("scan.date", "scan.time", "reader.id", "antenna.id", "dec.tag.id")
  dfall <- dfall[!is.na(dfall$dec.tag.id), ]
  dfall$dec.tag.id <- as.numeric(format(dfall$dec.tag.id, digits = 15))
  unique(dfall$dec.tag.id)
  dfall <- unique(dfall)
  dfall$read_tag <- dfall$dec.tag.id
  dfnew <- merge(dfall, metadata, by.x ="dec.tag.id", by.y = "tag_1")
  dfnew1 <- merge(dfall, metadata, by.x = "dec.tag.id", by.y = "tag_2")
  dfall <- rbind.fill(dfnew, dfnew1)
  dfall <- dfall[ , c("drop", "strain", "sex", "ID", "name", "code", 
                      "reader.id", "antenna.id", "scan.date", "scan.time", "read_tag","treatment")] 
  dfall$Field.Time <- as.POSIXct(paste(dfall$scan.date, dfall$scan.time), format="%m/%d/%Y %H:%M:%OS")
  
  # SORT DATAFRAME BY DATE AND TIME AND CREATE DATA
  data <- dfall[order(dfall$Field.Time), ]
  
  # CHANGE ANTENNA IDS/ZONES TO NUMERICS
  data$antenna.id <- as.numeric(data$antenna.id)
  data$ID <- as.numeric(data$ID)
  data$reader.id <- as.numeric(data$reader.id)
  
  # CREATE FULL_IDS 
  data$full_ids <- paste0(data$strain,"-", data$sex,"-", data$ID)
  
  # note that if more than 10 antennas are present, grepl will take the first digit and change it. requires anchors for exact match. \\b
  # CREATE X AND Y COORDINATES FOR THE ZONES FOR PLOTTING ON A GRID
  data$x <- ifelse(grepl("\\b1\\b", data$antenna.id), 4.12, 
                   ifelse(grepl("\\b2\\b", data$antenna.id), 11.12,
                          ifelse(grepl("\\b3\\b", data$antenna.id), 4.12,
                                 ifelse(grepl("\\b4\\b", data$antenna.id), 11.12,
                                        ifelse(grepl("\\b5\\b", data$antenna.id), 4.12,
                                               ifelse(grepl("\\b6\\b", data$antenna.id), 11.12,
                                                      ifelse(grepl("\\b7\\b", data$antenna.id), 4.12,
                                                             ifelse(grepl("\\b8\\b", data$antenna.id), 11.12,
                  ifelse(grepl("\\b9\\b",data$antenna.id), 7.62,
                         ifelse(grepl("\\b10\\b",data$antenna.id),7.62,
                                ifelse(grepl("\\b11\\b",data$antenna.id),7.62,
                                       ifelse(grepl("\\b12\\b",data$antenna.id),7.62,
                                              ifelse(grepl("\\b13\\b",data$antenna.id),5.87,
                                                     ifelse(grepl("\\b14\\b",data$antenna.id),9.37,
                                                            ifelse(grepl("\\b15\\b",data$antenna.id),5.87,
                                                                   ifelse(grepl("\\b16\\b",data$antenna.id),9.37,
                                                                    "none"))))))))))))))))
  
  
  data$y <- ifelse(grepl("\\b1\\b", data$antenna.id), 7.92, 
                   ifelse(grepl("\\b2\\b", data$antenna.id), 7.92,
                          ifelse(grepl("\\b3\\b", data$antenna.id), 15.24,
                                 ifelse(grepl("\\b4\\b", data$antenna.id), 15.24,
                                        ifelse(grepl("\\b5\\b", data$antenna.id), 22.86,
                                               ifelse(grepl("\\b6\\b", data$antenna.id), 22.86,
                                                      ifelse(grepl("\\b7\\b", data$antenna.id), 30.18,
                                                             ifelse(grepl("\\b8\\b", data$antenna.id), 30.18,
                    ifelse(grepl("\\b9\\b",data$antenna.id),7.92,
                           ifelse(grepl("\\b10\\b",data$antenna.id),15.24,
                                  ifelse(grepl("\\b11\\b",data$antenna.id),22.86,
                                         ifelse(grepl("\\b12\\b",data$antenna.id),30.18,
                                                ifelse(grepl("\\b13\\b",data$antenna.id),11.58,
                                                       ifelse(grepl("\\b14\\b",data$antenna.id),11.58,
                                                              ifelse(grepl("\\b15\\b",data$antenna.id),26.52,
                                                                     ifelse(grepl("\\b16\\b",data$antenna.id),26.52,
                                      "none"))))))))))))))))

  colnames(data)
  colnames(data) <- c("paddock", "strain", "sex", "ID", "name", "code", "reader_ID", 
                      "antenna", "date", "time", "read_tag", "treatment","field_time", "full_ID", "zone_x", "zone_y")
  rfid <- data
  class(rfid$field_time)
  rfid$time_sec <- as.numeric(difftime(rfid$field_time,min(rfid$field_time),units="secs")+1)
  origin <- as.POSIXct(paste(rfid$date[1], "12:00:00", sep =" "), format="%m/%d/%Y %H:%M:%OS")
  rfid$noon_day <- ceiling(difftime(rfid$field_time, origin,  units = "days"))
  rfid$days_antenna <- paste(rfid$noon_day,rfid$antenna,sep="_")
  
  ## MERGE WEATHER DATA WITH RFID DATA
  rfid$weather_time <- round(as.POSIXct(rfid$field_time), "hour")
  class(rfid$weather_time)
  rfid <- left_join(rfid,weather, by= "weather_time")
  rfid1 <- rfid %>% 
    mutate(date = date.x, trial = trial_var) %>%
    relocate(date, .before=time) %>% 
    relocate(trial,.before=paddock) %>% 
    select(!(c(date.x, date.y)))
  
  ## WRITE THE DATA TO CSV INTO EACH FOLDER
  write.csv(rfid1, file = paste0(folders_short[flag],'_RFID_FULL_DATA.csv'))
  all_trial_list[[i]] <- rfid1
  flag <- flag+1
}
all_trial_data = do.call(bind_rows, all_trial_list)
setwd(wd)
write.csv(all_trial_data, "ALLTRIAL_RFID_DATA.csv")

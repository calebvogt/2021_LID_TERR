## Created by Caleb C. Vogt, PhD Candidate @ Cornell University
library(tidyverse)
library(data.table)

meta <- read_excel("data/LID_TERR_2021_metadata.xlsx", sheet = 1, skip = 0)
meta$num_terr<-ifelse(meta$name=="Cameron"|meta$name=="Carter"|meta$name=="Diego"|meta$name=="Josiah",2,
                      ifelse(meta$name=="Isaiah"|meta$name=="Ibex"|meta$name=="Dillard"|meta$name=="Jeremy",1,0))
meta$drop <- ifelse(meta$treatment=="early",1,2)
meta$strain_sex_name <- paste(meta$strain,meta$sex,meta$name,sep="-")
write.csv(meta, "data/metadata.csv")

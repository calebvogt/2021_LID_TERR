## Created by Matthew Zipple
library(export)

social_summary_data <- as.data.frame(fread("LID_2021_TERR_GBI_SUMMARY.csv", stringsAsFactors = FALSE))

# social_summary_data: MZ_Prop_female time with 0,1,2 territory males ---------------------------------------------------

master_stats <- social_summary_data
agg<-aggregate(master_stats$duration_s,list(master_stats$name,master_stats$f_sum,master_stats$day,master_stats$sex),FUN=sum)
names(agg)<-c("name","females","day","sex","duration")
agg<-subset(agg,sex=="M")

agg$terr<-ifelse(agg$name=="Cameron"|agg$name=="Carter"|agg$name=="Diego"|agg$name=="Josiah",2,
                 ifelse(agg$name=="Isaiah"|agg$name=="Ibex"|agg$name=="Dillard"|agg$name=="Jeremy",1,0))

agg<-subset(agg,day<=20)
agg$total_female_time<-as.numeric(agg$females*agg$duration)

agg_2<-aggregate(agg$duration,list(agg$terr,agg$day,agg$females),FUN=mean)
names(agg_2)<-c("terr","day","females","duration")

plot(duration~day,subset(agg_2,terr==2&females==1),type="l",col="darkorchid4",ylim=c(0,10000))
points(duration~day,subset(agg_2,terr==1&females==1),type="l",col="goldenrod4")
points(duration~day,subset(agg_2,terr==0&females==1),type="l",col="gray40")

plot(duration~day,subset(agg_2,terr==2&females==2),type="l",col="darkorchid4",ylim=c(0,2000))
points(duration~day,subset(agg_2,terr==1&females==2),type="l",col="goldenrod4")
points(duration~day,subset(agg_2,terr==0&females==2),type="l",col="gray40")

plot(duration~day,subset(agg_2,terr==2&females==3),type="l",col="darkorchid4",ylim=c(0,2000))
points(duration~day,subset(agg_2,terr==1&females==3),type="l",col="goldenrod4")
points(duration~day,subset(agg_2,terr==0&females==3),type="l",col="gray40")


agg_3<-aggregate(agg$total_female_time,list(agg$terr,agg$day),FUN=sum)
names(agg_3)<-c("terr","day","total_time")

plot(total_time/4~day,subset(agg_3,terr==2),type="l",col="darkorchid4",ylim=c(0,10000))
points(total_time/4~day,subset(agg_3,terr==1),type="l",col="goldenrod4")

summary(lmer(total_female_time~terr+(1|name),data=subset(agg,terr>0&day<=5)))
summary(lmer(total_female_time~terr+(1|name),data=subset(agg,terr>0&day==6)))
summary(lmer(total_female_time~terr+(1|name),data=subset(agg,terr>0&day>=6)))

summary(lmer(as.numeric(duration)~terr+(1|name),data=subset(agg,terr>0&day<=5&females==1)))
summary(lmer(as.numeric(duration)~terr+(1|name),data=subset(agg,terr>0&day==6&females==1)))
summary(lmer(as.numeric(duration)~terr+(1|name),data=subset(agg,terr>0&day>=6&females==1)))

summary(lmer(as.numeric(duration)~terr+(1|name),data=subset(agg,terr>0&day<=5&females==2)))
summary(lmer(as.numeric(duration)~terr+(1|name),data=subset(agg,terr>0&day==6&females==2)))
summary(lmer(as.numeric(duration)~terr+(1|name),data=subset(agg,terr>0&day>=6&females==2)))


summary(subset(df2,df2$'C57-F-Nicole'==1&df2$m_sum>0&day<=20))
as.data.frame(summary(subset(df2,df2$'C57-F-Nadia'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                               136,142,160,166,184,
                                                                               190,202,226,232,238,268,
                                                                               274,286,292),]
as.data.frame(summary(subset(df2,df2$'C57-F-Nadia'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                               136,142,160,166,184,
                                                                               190,202,226,232,238,268,
                                                                               274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Nancy'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                               136,142,160,166,184,
                                                                               190,202,226,232,238,268,
                                                                               274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Lily'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                              136,142,160,166,184,
                                                                              190,202,226,232,238,268,
                                                                              274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Lucy'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                              136,142,160,166,184,
                                                                              190,202,226,232,238,268,
                                                                              274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Nina'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                              136,142,160,166,184,
                                                                              190,202,226,232,238,268,
                                                                              274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Lorena'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                                136,142,160,166,184,
                                                                                190,202,226,232,238,268,
                                                                                274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Luna'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                              136,142,160,166,184,
                                                                              190,202,226,232,238,268,
                                                                              274,286,292),]


as.data.frame(summary(subset(df2,df2$'C57-F-Oasis'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                               136,142,160,166,184,
                                                                               190,202,226,232,238,268,
                                                                               274,286,292),]


as.data.frame(summary(subset(df2,df2$'C57-F-Kathryn'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                                 136,142,160,166,184,
                                                                                 190,202,226,232,238,268,
                                                                                 274,286,292),]


as.data.frame(summary(subset(df2,df2$'C57-F-Kittie'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                                136,142,160,166,184,
                                                                                190,202,226,232,238,268,
                                                                                274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Kayla'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                               136,142,160,166,184,
                                                                               190,202,226,232,238,268,
                                                                               274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Melody'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                                136,142,160,166,184,
                                                                                190,202,226,232,238,268,
                                                                                274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Ophelia'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                                 136,142,160,166,184,
                                                                                 190,202,226,232,238,268,
                                                                                 274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Meredith'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                                  136,142,160,166,184,
                                                                                  190,202,226,232,238,268,
                                                                                  274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Meredith'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                                  136,142,160,166,184,
                                                                                  190,202,226,232,238,268,
                                                                                  274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Oakland'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                                 136,142,160,166,184,
                                                                                 190,202,226,232,238,268,
                                                                                 274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Olivia'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                                136,142,160,166,184,
                                                                                190,202,226,232,238,268,
                                                                                274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Kista'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                               136,142,160,166,184,
                                                                               190,202,226,232,238,268,
                                                                               274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Madison'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                                 136,142,160,166,184,
                                                                                 190,202,226,232,238,268,
                                                                                 274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Morgan'==1&df2$m_sum>0&day<=20)))[c(82,88,94,100,112,124,
                                                                                136,142,160,166,184,
                                                                                190,202,226,232,238,268,
                                                                                274,286,292),]





summary(subset(df2,df2$'C57-F-Nicole'==1&df2$m_sum>0&day<=20&day>=6))
as.data.frame(summary(subset(df2,df2$'C57-F-Nadia'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                      136,142,160,166,184,
                                                                                      190,202,226,232,238,268,
                                                                                      274,286,292),]
as.data.frame(summary(subset(df2,df2$'C57-F-Nadia'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                      136,142,160,166,184,
                                                                                      190,202,226,232,238,268,
                                                                                      274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Nancy'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                      136,142,160,166,184,
                                                                                      190,202,226,232,238,268,
                                                                                      274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Lily'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                     136,142,160,166,184,
                                                                                     190,202,226,232,238,268,
                                                                                     274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Lucy'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                     136,142,160,166,184,
                                                                                     190,202,226,232,238,268,
                                                                                     274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Nina'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                     136,142,160,166,184,
                                                                                     190,202,226,232,238,268,
                                                                                     274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Lorena'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                       136,142,160,166,184,
                                                                                       190,202,226,232,238,268,
                                                                                       274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Luna'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                     136,142,160,166,184,
                                                                                     190,202,226,232,238,268,
                                                                                     274,286,292),]


as.data.frame(summary(subset(df2,df2$'C57-F-Oasis'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                      136,142,160,166,184,
                                                                                      190,202,226,232,238,268,
                                                                                      274,286,292),]


as.data.frame(summary(subset(df2,df2$'C57-F-Kathryn'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                        136,142,160,166,184,
                                                                                        190,202,226,232,238,268,
                                                                                        274,286,292),]


as.data.frame(summary(subset(df2,df2$'C57-F-Kittie'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                       136,142,160,166,184,
                                                                                       190,202,226,232,238,268,
                                                                                       274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Kayla'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                      136,142,160,166,184,
                                                                                      190,202,226,232,238,268,
                                                                                      274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Melody'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                       136,142,160,166,184,
                                                                                       190,202,226,232,238,268,
                                                                                       274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Ophelia'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                        136,142,160,166,184,
                                                                                        190,202,226,232,238,268,
                                                                                        274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Meredith'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                         136,142,160,166,184,
                                                                                         190,202,226,232,238,268,
                                                                                         274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Meredith'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                         136,142,160,166,184,
                                                                                         190,202,226,232,238,268,
                                                                                         274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Oakland'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                        136,142,160,166,184,
                                                                                        190,202,226,232,238,268,
                                                                                        274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Olivia'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                       136,142,160,166,184,
                                                                                       190,202,226,232,238,268,
                                                                                       274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Kista'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                      136,142,160,166,184,
                                                                                      190,202,226,232,238,268,
                                                                                      274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Madison'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                        136,142,160,166,184,
                                                                                        190,202,226,232,238,268,
                                                                                        274,286,292),]

as.data.frame(summary(subset(df2,df2$'C57-F-Morgan'==1&df2$m_sum>0&day<=20&day>=6)))[c(82,88,94,100,112,124,
                                                                                       136,142,160,166,184,
                                                                                       190,202,226,232,238,268,
                                                                                       274,286,292),]
fem<-read.csv("female_prop_males.csv")
mod<-glmmTMB(log((top.per.6.20/100)/(1-top.per.6.20/100))~terr+(1|male),data=fem)
mod1<-glmmTMB(log((top.per.6.20/100)/(1-top.per.6.20/100))~terr,data=fem)
summary(mod)
qqnorm(resid(mod))        
qqline(resid(mod))

boxplot(top.per.6.20~terr,data=fem,cex.axis=1.4)
library(export)
#graph2ppt(file="female_monop",height=7,width=7)



df_duration<-cbind(df2[,1:9],df2[,10:49]*as.numeric(df2$duration_s))

sum(subset(df_duration,df_duration$'C57-F-Madison'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Madison'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Madison'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Madison'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Nicole'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Nicole'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Nicole'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Nicole'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Nadia'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Nadia'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Nadia'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Nadia'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Nancy'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Nancy'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Nancy'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Nancy'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Lily'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Lily'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Lily'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Lily'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Lucy'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Lucy'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Lucy'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Lucy'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Nina'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Nina'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Nina'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Nina'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Lorena'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Lorena'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Lorena'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Lorena'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Luna'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Luna'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Luna'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Luna'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Oasis'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Oasis'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Oasis'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Oasis'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Kathryn'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Kathryn'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Kathryn'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Kathryn'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Kittie'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Kittie'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Kittie'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Kittie'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Kayla'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Kayla'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Kayla'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Kayla'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Melody'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Melody'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Melody'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Melody'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Ophelia'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Ophelia'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Ophelia'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Ophelia'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Meredith'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Meredith'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Meredith'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Meredith'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Oakland'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Oakland'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Oakland'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Oakland'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Olivia'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Olivia'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Olivia'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Olivia'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Kista'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Kista'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Kista'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Kista'>=1,1,0)),FUN=sum)

sum(subset(df_duration,df_duration$'C57-F-Morgan'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)])
aggregate(cbind(subset(df_duration,df_duration$'C57-F-Morgan'>=1&m_sum>0&day<=20&day>=6)[,c(14,15,16,17,19,21,23,24,27,28,31,32,34,38,39,40,45,46,48,49)]),
          by=list(ifelse(subset(df_duration,df_duration$'C57-F-Morgan'>=1&m_sum>0&day<=20&day>=6)$'C57-F-Morgan'>=1,1,0)),FUN=sum)

fem<-read.csv("female_prop_males.csv")
mod<-glmmTMB(log((top_dur_prop)/(1-top_dur_prop))~dur_terr+(1|dur_male),data=fem)
summary(mod)
summary(mod1)

qqnorm(resid(mod))        
qqline(resid(mod))

boxplot(top_dur_prop~dur_terr,data=fem,cex.axis=1.4)
library(export)




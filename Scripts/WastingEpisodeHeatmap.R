
# load packages
rm(list=ls())
library(tidyverse)

setwd("U:/data/WastIncDatasets")


#-------------------------------------
# Counting birth incidence
#-------------------------------------
load("cmc_inc.Rdata")
load("eu_inc.Rdata")
load("irc_inc.Rdata")
load("tdc_inc.Rdata")
load("cort_inc.Rdata")
load("mled_inc.Rdata")
load("dvds_inc.Rdata")
load("vita_inc.Rdata")
load("vb12_inc.Rdata")
load("zlbw_inc.Rdata")
load("zsga_inc.Rdata")
load("zinf_inc.Rdata")
load("cmpf_inc.Rdata")
load("fspp_inc.Rdata")
load("zmrt_inc.Rdata")



#Define set of short ID's I care about
datasets <- list(
  cmc_inc, 
  eu_inc, 
  irc_inc,
  tdc_inc,
  cort_inc_india,            
  mled_inc_india,
  #dvds_inc,
  vita_inc,
  vb12_inc,
  #zlbw_inc,
  #zsga_inc,
  #zinf_inc,
  cmpf_inc,
  fspp_inc,
  zmrt_inc)


d<-NULL

for(i in 1:length(datasets)){
  temp<-datasets[[i]] %>% subset(., select=c(STUDYID, SUBJID, AGEDAYS, WHZ, wasting_episode, sevwasting_episode)) %>% 
    arrange(SUBJID, AGEDAYS) %>% mutate(agemonth=floor(AGEDAYS/30.25), meanWHZ=mean(WHZ, na.rm=T), n=rank(meanWHZ)) %>%
    group_by(SUBJID, agemonth) %>% slice(1) %>% filter(AGEDAYS < 24 * 30.25) %>% ungroup() 
  d<-rbind(d, temp)
}

d$individ <- paste0(d$STUDYID,"-",d$SUBJID)
#d$n <- as.numeric(factor(d$individ, levels=unique(d$individ)))
d$n <- rank(d$meanWHZ)
d$wasting_episode[d$sevwasting_episode=="Severe Wasted"] <- "Severe Wasted"

ggplot(data=d, aes(x=agemonth, y=n, group=n, color=wasting_episode, alpha=0.5)) + geom_line() + theme_void() #+ 
  #facet_wrap(~STUDYID, nrow=1, scales="free")




d <- d %>% group_by(STUDYID,SUBJID) %>% 
  mutate(meanWHZ=mean(WHZ, na.rm=T)) %>%
  ungroup() %>%
  filter(AGEDAYS==min(AGEDAYS) | AGEDAYS==max(AGEDAYS) | sevwast_inc==1 | wast_inc==1 | wast_rec==1) %>% 
  arrange(meanWHZ, STUDYID, SUBJID)

d$incident_age[is.na(d$incident_age)] <- d$AGEDAYS[is.na(d$incident_age)] 

#d$n <- (1:nrow(d))*2

d$wasting_episode[d$sevwasting_episode=="Severe Wasted"] <- "Severe Wasted"

d$individ <- paste0(d$STUDYID,"-",d$SUBJID)
d$n <- as.numeric(factor(d$individ, levels=unique(d$individ)))*100

ggplot(data=d, aes(x=incident_age, y=n, group=individ, color=wasting_episode)) + geom_line() + theme_bw()


ggplot(data=d[d$individ==unique(d$individ)[345],], aes(x=incident_age, y=n, group=individ, color=STUDYID)) + geom_line() + theme_bw()
length(unique(d$individ))


df <- mled_inc_india

#df$n <- as.numeric(factor(df$SUBJID, levels=unique(df$SUBJID)))
#ggplot(data=df, aes(x=AGEDAYS, y=n, group=SUBJID, color=wasting_episode)) + geom_line() + theme_bw()

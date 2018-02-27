





rm(list=ls())
library(tidyverse)
library(survminer)


setwd("U:/R scripts")
source("Stunt_incidence_functions.R")

load("U:/data/Compiled Datasets/6A_long.Rdata")
setwd("U:/data")



d <- d %>% subset(., select=c(STUDYID, SUBJID, AGEDAYS, WHZ, HAZ))


#Drop cross sectional studies
unique(d$STUDYID)
d <- d %>% filter(STUDYID!="ki1000302-TTBW" & STUDYID!="ki1000304-NEOVITA" & STUDYID!="ki1000303-BtS" & STUDYID!="ki1000305-CordBloodCytometry" & STUDYID!="ki1000303-NBSP" & STUDYID!="ki1000304c-IMNCI")

#Drop studies that enrol acutely ill children
d <- d %>% filter(STUDYID!="ki1000301-DIVIDS" & STUDYID!="ki1000304-LBW" & STUDYID!="ki1000306-ZincSGA" & STUDYID!="ki1000304b-ZincInf")


head(d)


# Calculate failure time for cumulative incidence plots

survwast <- d %>% group_by(STUDYID, SUBJID) %>% 
              mutate(wast = as.numeric(WHZ < (-2)),
                     wastmeasure = cumsum(wast),
                    numevents = max(wastmeasure),
                  maxage = max(AGEDAYS)) %>%
  filter(wastmeasure==1 | (numevents==0 & AGEDAYS==maxage)) %>%
  mutate(event=ifelse(wastmeasure==1 ,"wast","no")) %>%
  subset(., select = c(STUDYID, SUBJID, AGEDAYS, event))


survstunt <- d %>% group_by(STUDYID, SUBJID) %>% 
  mutate(stunt = as.numeric(HAZ < (-2)),
         stuntmeasure = cumsum(stunt),
         numevents = max(stuntmeasure),
         maxage = max(AGEDAYS)) %>%
  filter(stuntmeasure==1 | (numevents==0 & AGEDAYS==maxage)) %>%
  mutate(event=ifelse(stuntmeasure==1 ,"stunt","no")) %>%
  subset(., select = c(STUDYID, SUBJID, AGEDAYS, event))


summary(survwast$AGEDAYS)
summary(survstunt$AGEDAYS)


surv <- rbind(survwast, survstunt)



# handles cuminc objects
fit <- cmprsk::cuminc(surv$AGEDAYS, factor(surv$event), cencode="no")
ggcompetingrisks(fit, conf.int = TRUE) + coord_cartesian(xlim=c(0,1000))

fit <- cmprsk::cuminc(surv$AGEDAYS, factor(surv$event), factor(surv$STUDYID), cencode="no")
ggcompetingrisks(fit, conf.int = TRUE) + coord_cartesian(xlim=c(0,1000))


fit <- cmprsk::cuminc(survwast$AGEDAYS, factor(survwast$event), cencode="no")
ggcompetingrisks(fit, conf.int = TRUE) + coord_cartesian(xlim=c(0,720))

fit <- cmprsk::cuminc(survstunt$AGEDAYS, factor(survstunt$event), cencode="no")
ggcompetingrisks(fit, conf.int = TRUE) + coord_cartesian(xlim=c(0,720))

fit <- cmprsk::cuminc(survstunt$AGEDAYS, factor(survstunt$event), factor(survstunt$STUDYID), cencode="no")
ggcompetingrisks(fit, conf.int = TRUE)


# df <- d %>% group_by(STUDYID, SUBJID) %>% 
#   mutate(anystunt06=ifelse(any(stuntinc[AGEDAYS<6*30.25]==1),1,0),
#          anystunt012=ifelse(any(stuntinc[AGEDAYS<12*30.25]==1),1,0),
#          anystunt018=ifelse(any(stuntinc[AGEDAYS<18*30.25]==1),1,0),
#          anystunt024=ifelse(any(stuntinc[AGEDAYS<24*30.25]==1),1,0)) %>% 
#   summarize(Neps06=sum(anystunt06, na.rm=T),
#             Neps012=sum(anystunt012, na.rm=T),
#             Neps018=sum(anystunt018, na.rm=T),
#             Neps024=sum(anystunt024, na.rm=T))


#Filter out children older than 24 months
d <- d %>% filter(AGEDAYS < 24 *30.25)


ind.df <- d %>% group_by(STUDYID, SUBJID) %>% 
  mutate(anystunt06=ifelse(any(HAZ[AGEDAYS<6*30.25] < (-2)),1,0),
         anystunt012=ifelse(any(HAZ[AGEDAYS<12*30.25] < (-2)),1,0),
         anystunt018=ifelse(any(HAZ[AGEDAYS<18*30.25] < (-2)),1,0),
         anystunt024=ifelse(any(HAZ[AGEDAYS<24*30.25] < (-2)),1,0),
         lead.age=lead(AGEDAYS),
         lag.age=lag(AGEDAYS),
         for.period=(lead.age-AGEDAYS)/2,
         back.period=(AGEDAYS-lag.age)/2,
         period.length=sum(for.period,back.period, na.rm=T),
         meas06=sum(AGEDAYS<6*30.25) > 0,
         meas012=sum(AGEDAYS<12*30.25) > 0,
         meas018=sum(AGEDAYS<18*30.25) > 0,
         meas024=sum(AGEDAYS<24*30.25) > 0
         ) 

Ndf <- ind.df %>% 
  ungroup() %>% group_by(STUDYID) %>%
  summarize(Neps06=sum(anystunt06, na.rm=T),
            Neps012=sum(anystunt012, na.rm=T),
            Neps018=sum(anystunt018, na.rm=T),
            Neps024=sum(anystunt024, na.rm=T),
            pt06=sum(period.length[anystunt06==0], na.rm=T),
            pt012=sum(period.length[anystunt012==0], na.rm=T),
            pt018=sum(period.length[anystunt018==0], na.rm=T),
            pt024=sum(period.length[anystunt024==0], na.rm=T),
            N=n(),
            N06=sum(meas06),
            N012=sum(meas012),
            N018=sum(meas018),
            N024=sum(meas024)
            )


Ndf <- Ndf %>% filter(pt06!=0)




#Wasting
ind.df.wast <- d %>% group_by(STUDYID, SUBJID) %>% 
  mutate(anywast06=ifelse(any(WHZ[AGEDAYS<6*30.25] < (-2)),1,0),
         anywast012=ifelse(any(WHZ[AGEDAYS<12*30.25] < (-2)),1,0),
         anywast018=ifelse(any(WHZ[AGEDAYS<18*30.25] < (-2)),1,0),
         anywast024=ifelse(any(WHZ[AGEDAYS<24*30.25] < (-2)),1,0),
         lead.age=lead(AGEDAYS),
         lag.age=lag(AGEDAYS),
         for.period=(lead.age-AGEDAYS)/2,
         back.period=(AGEDAYS-lag.age)/2,
         period.length=sum(for.period,back.period, na.rm=T),
         meas06=sum(AGEDAYS<6*30.25) > 0,
         meas012=sum(AGEDAYS<12*30.25) > 0,
         meas018=sum(AGEDAYS<18*30.25) > 0,
         meas024=sum(AGEDAYS<24*30.25) > 0
  ) 

Ndf.wast <- ind.df.wast %>% 
  ungroup() %>% group_by(STUDYID) %>%
  summarize(Neps06=sum(anywast06, na.rm=T),
            Neps012=sum(anywast012, na.rm=T),
            Neps018=sum(anywast018, na.rm=T),
            Neps024=sum(anywast024, na.rm=T),
            pt06=sum(period.length[anywast06==0], na.rm=T),
            pt012=sum(period.length[anywast012==0], na.rm=T),
            pt018=sum(period.length[anywast018==0], na.rm=T),
            pt024=sum(period.length[anywast024==0], na.rm=T),
            N=n(),
            N06=sum(meas06),
            N012=sum(meas012),
            N018=sum(meas018),
            N024=sum(meas024)
  )


Ndf.wast<- Ndf.wast %>% filter(pt06!=0)


save(Ndf, Ndf.wast, file="U:/Rally 6 India/Results/CI_df_6A.Rdata")

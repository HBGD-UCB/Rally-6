



rm(list=ls())
library(caret)
library(tidyverse)
theme_set(theme_bw())

source("U:/R scripts/Wast_incidence_functions.R")


#hbgdki pallet
tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728",
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")



setwd("U:/data/WastIncDatasets")


load("U:/data/Compiled Datasets/6A_long.Rdata")
d<-d[,-c(1:2)]
colnames(d)

d <- d %>% filter(!is.na(HAZ)) %>% 
  filter(HAZ > (-5) & HAZ < 5) %>% 
  filter(AGEDAYS < 24*30.25 & !is.na(AGEDAYS)) %>%
  filter(!is.na(enrolHAZ) & enrolHAZ > (-5) & enrolHAZ < 5)

d$enrolStunt <- as.numeric(d$enrolHAZ < (-2))
d$enrolWast <- as.numeric(d$enrolWHZ < (-2))


#remove grant identifiers
d$STUDYID<- gsub("^k.*?-" , "", d$STUDYID)

#Drop cross sectional studies, ages over 2 years, and missing anthropometry
unique(d$STUDYID)
d <- d %>% filter(AGEDAYS < 24*30 & !is.na(WHZ)) %>%
  filter(STUDYID!="TTBW" & STUDYID!="CordBloodCytometry" & STUDYID!="IMNCI")



#Seperate out cohorts that enrolled low anthro or sick children
d_small <- d %>% filter(STUDYID=="ZincSGA" | STUDYID=="DIVIDS" | STUDYID=="LBW" | STUDYID=="ZincInf")
d <- d %>% filter(STUDYID!="ZincSGA" & STUDYID!="DIVIDS" & STUDYID!="LBW" & STUDYID!="ZincInf")



p1<-ggplot(aes(x=AGEDAYS, y=WHZ, colour=factor(enrolStunt)), data=d) +
  geom_smooth() + facet_wrap(~STUDYID)  + 
  scale_color_manual(values=tableau10) +
  ggtitle("Cohort-specific WHZ stratified by stunting status at birth or enrollment") + 
  theme(strip.background = element_blank())
p1


p2<-ggplot(aes(x=AGEDAYS, y=HAZ, colour=factor(enrolStunt)), data=d) +
  geom_smooth() + facet_wrap(~STUDYID)  + 
  scale_color_manual(values=tableau10) +
  ggtitle("Cohort-specific HAZ stratified by stunting status at birth or enrollment") + 
  theme(strip.background = element_blank())
p2


p3<-ggplot(aes(x=AGEDAYS, y=WHZ, colour=factor(enrolWast)), data=d) +
  geom_smooth() + facet_wrap(~STUDYID)  + 
  scale_color_manual(values=tableau10) +
  ggtitle("Cohort-specific WHZ stratified by wasting status at birth or enrollment") + 
  theme(strip.background = element_blank())
p3


p4<-ggplot(aes(x=AGEDAYS, y=HAZ, colour=factor(enrolWast)), data=d) +
  geom_smooth() + facet_wrap(~STUDYID)  + 
  scale_color_manual(values=tableau10) +
  ggtitle("Cohort-specific HAZ stratified by wasting status at birth or enrollment") + 
  theme(strip.background = element_blank())
p4



d %>% group_by(STUDYID, SUBJID) %>% arrange(AGEDAYS) %>% slice(1) %>% ungroup() %>%
  group_by(STUDYID, enrolStunt) %>% summarize(mean=mean(HAZ, na.rm=T), max=max(HAZ, na.rm=T)) 


df <- d %>%  group_by(STUDYID, SUBJID) %>% arrange(AGEDAYS) %>% slice(1) %>% subset(., select=c(HAZ, enrolHAZ, WHZ, enrolWHZ, enrolStunt, enrolWast))


df2 <- d %>%  group_by(STUDYID, SUBJID) %>% arrange(STUDYID, SUBJID, AGEDAYS) %>% subset(., select=c(STUDYID, SUBJID, AGEDAYS, HAZ, enrolHAZ, WHZ, enrolWHZ, enrolStunt, enrolWast))

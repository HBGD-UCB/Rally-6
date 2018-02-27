
rm(list=ls())
library(tidyverse)

#hbgdki pallet
tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728",
               "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")



load("U:/data/Compiled Datasets/6A_long.Rdata")

#Drop extra study variables
d<-d[,-c(1:2)]

#remove grant identifiers
d$STUDYID<- gsub("^k.*?-" , "", d$STUDYID)

#Drop cross sectional studies, ages over 2 years, and missing anthropometry
unique(d$STUDYID)
d <- d %>% filter(AGEDAYS < 24*30 & !is.na(WHZ)) %>%
           filter(STUDYID!="TTBW" & STUDYID!="CordBloodCytometry" & STUDYID!="IMNCI")




#Drop ages over 2 years
d <- d %>% filter(AGEDAYS < 24*30 & !is.na(WHZ))

whz <- d %>% select(STUDYID, AGEDAYS, BRTHWEEK, WHZ) %>% rename(zscore=WHZ) %>% mutate(measure="WHZ")
haz <- d %>% select(STUDYID, AGEDAYS, BRTHWEEK, HAZ) %>% rename(zscore=HAZ) %>% mutate(measure="HAZ")

anthro<-rbind(whz, haz)

#Seperate out cohorts that enrolled low anthro or sick children
anthro_small <- anthro %>% filter(STUDYID=="ZincSGA" | STUDYID=="DIVIDS" | STUDYID=="LBW" | STUDYID=="ZincInf")
anthro <- anthro %>% filter(STUDYID!="ZincSGA" & STUDYID!="DIVIDS" & STUDYID!="LBW" & STUDYID!="ZincInf")


#Plot WHZ and HAZ from each cohort
p1<-ggplot(aes(x=AGEDAYS, y=zscore, colour=measure), data=anthro) +
  geom_smooth() + facet_wrap(~STUDYID)  + 
  scale_color_manual(values=tableau10) +
  ggtitle("Cohort-specific WHZ and HAZ over child ages") + 
  theme(strip.background = element_blank())
p1


p2<-ggplot(aes(x=AGEDAYS, y=zscore, colour=measure), data=anthro_small) +
  geom_smooth() + facet_wrap(~STUDYID)  + 
  scale_color_manual(values=tableau10) +
  ggtitle("Cohort-specific WHZ and HAZ over child ages:\nstudies enrolling ill children") + 
  theme(strip.background = element_blank())
p2


#-------------------------------
#Plot anthro over season
#-------------------------------

d <- anthro %>% filter(!is.na(BRTHWEEK))


table(d$BRTHWEEK)
d$julianday <- as.numeric(d$BRTHWEEK)*7 + d$AGEDAYS
d$julianday[d$julianday>365] <- d$julianday[d$julianday>365] - 365
d$julianday[d$julianday>365] <- d$julianday[d$julianday>365] - 365
d$julianday[d$julianday>365] <- d$julianday[d$julianday>365] - 365
d$julianday[d$julianday>365] <- d$julianday[d$julianday>365] - 365
summary(d$julianday)

p3<-ggplot(aes(x=julianday, y=zscore, colour=measure), data=d) +
  geom_smooth() + facet_wrap(~STUDYID)  + 
  scale_color_manual(values=tableau10) +
  ggtitle("Cohort-specific WHZ and HAZ over day of the year") + 
  theme(strip.background = element_blank())
p3




d <- anthro_small %>% filter(!is.na(BRTHWEEK))


table(d$BRTHWEEK)
d$julianday <- as.numeric(d$BRTHWEEK)*7 + d$AGEDAYS
d$julianday[d$julianday>365] <- d$julianday[d$julianday>365] - 365
d$julianday[d$julianday>365] <- d$julianday[d$julianday>365] - 365
d$julianday[d$julianday>365] <- d$julianday[d$julianday>365] - 365
d$julianday[d$julianday>365] <- d$julianday[d$julianday>365] - 365
summary(d$julianday)

p4<-ggplot(aes(x=julianday, y=zscore, colour=measure), data=d) +
  geom_smooth() + facet_wrap(~STUDYID)  + 
  scale_color_manual(values=tableau10) +
  ggtitle("Cohort-specific WHZ and HAZ over day of the year:\nstudies enrolling ill children") + 
  theme(strip.background = element_blank())
p4


#Check for gaps in ages
m <- ggplot(d, aes(julianday))
m + geom_histogram() + facet_wrap(~STUDYID)



#----------------------
# Primary plots
#----------------------


#Set plot directory
setwd("U:/Rally 6 India/Figures")

#Set plot width and height
w <- 10
h <- 7

ggsave("6a_growth_trajectories.pdf", p1, width = w, height = h, units = "in")
ggsave("6a_growth_trajectories_illchildren.pdf", p1, width = w, height = h, units = "in")
ggsave("6a_seasonal_Zscores.pdf", p1, width = w, height = h, units = "in")
ggsave("6a_seasonal_Zscores_illchildren.pdf", p1, width = w, height = h, units = "in")

pdf("6A Zscores over age and season.pdf", width=w,height=h, paper="USr")
p1
p2
p3
p4
dev.off()


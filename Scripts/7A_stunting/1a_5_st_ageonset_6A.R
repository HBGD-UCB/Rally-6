#-----------------------------------
# Stunting analysis
# Objective 1a
# Mean age of stunting onset
# Calculate in monthly cohorts only
#-----------------------------------
rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(metafor)
theme_set(theme_bw())

# load random effects function
source("U:/Scripts/Stunting/2-analyses/0_randomeffects.R")

load("U:/Data/Stunting/stunting_data.RData")
d <- d%>% filter(country=="INDIA")

# check that measurements are approx monthly
keep=d %>%
  group_by(studyid,country,subjid) %>%
  mutate(agelag=lag(agedays),
         agediff=agedays-agelag) %>%
  group_by(studyid,country) %>%
  summarise(min=min(agediff,na.rm=TRUE),
            mean=mean(agediff,na.rm=TRUE),
            max=max(agediff,na.rm=TRUE)) %>%
  filter(mean<=40) %>%
  select(studyid,country) %>%
  mutate(keep=1)

dmonthly=inner_join(d,keep)

# ---------------------------------------
# identify incident cases
# ---------------------------------------
sem<-function(x){
  sd(x)/sqrt(length(x))
}


inc = dmonthly %>%
  # ungroup(studyid) %>%
  # mutate(studyid=as.factor(studyid)) %>%
  group_by(studyid,subjid) %>%
  arrange(studyid,subjid,agedays) %>%
  
  # create id for measurement within person
  mutate(measid=seq_along(subjid)) %>%

  # create indicator for whether haz at t < haz at t-1
  mutate(hazlag=lag(haz)) %>%
  mutate(newcase=ifelse(measid==1,ifelse(haz< -2, 1,0),
                        ifelse(hazlag>= -2 & haz< -2,1,0))) %>%
  mutate(newcaselag=lag(newcase))%>%
  mutate(newcaselag=ifelse(measid==1,0,newcaselag))%>%
  mutate(cnewcaselag=cumsum(newcaselag)) %>%
  mutate(inccase=ifelse(cnewcaselag>=1,0,newcase)) %>%
  filter(inccase==1) %>%
  select(studyid,country,subjid,agedays)
 
# cohort specific mean
age.onset.study = inc %>%
  # make region
  mutate(region = ifelse(country=="BANGLADESH" | country=="INDIA"|
                  country=="NEPAL" | country=="PAKISTAN"|
                  country=="PHILIPPINES" ,"Asia",
                ifelse(country=="BURKINA FASO"|
                         country=="GUINEA-BISSAU"|
                         country=="MALAWI"|
                         country=="SOUTH AFRICA"|
                         country=="TANZANIA, UNITED REPUBLIC OF"|
                         country=="ZIMBABWE"|
                         country=="GAMBIA","Africa",
                       ifelse(country=="BELARUS","Europe",
                              "Latin America")))) %>%
  # concatenate country and study
  mutate(study_country=paste0(studyid,"-",country)) %>%
  group_by(region,study_country) %>%
  summarise(mn=mean(agedays),
            se=sem(agedays),
            Nmeas=n(),
            Nchild=sum(length(unique(subjid)))) %>%
  mutate(lb=se-qnorm(0.975)*se,
         ub=se+qnorm(0.975)*se) 

# sort by mean age
age.onset.study$study_country=factor(age.onset.study$study_country, 
            levels = age.onset.study$study_country[order(age.onset.study$mn)])

# estimate random effects, format results
pool.fit=rma(yi=age.onset.study$mn, 
             sei=age.onset.study$se,
             method="REML")

# results:
c(est=pool.fit$beta, se=pool.fit$se, lb=pool.fit$ci.lb, ub=pool.fit$ci.ub)


ggplot(age.onset.study,aes(x=study_country,y=mn))+
  geom_point(aes(size=Nchild,col=region))+coord_flip()+
  geom_hline(yintercept=pool.fit$ci.lb,linetype="dashed")+
  geom_hline(yintercept=pool.fit$ci.ub,linetype="dashed")+
  geom_hline(yintercept=pool.fit$beta)+
  scale_y_continuous(breaks=seq(0,850,50),labels=seq(0,850,50))+
  xlab("Study & Country") + ylab("Age in days")





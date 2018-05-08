#-----------------------------------
# Stunting analysis
# Objective 1a
# Calculate cumulative incidence (ever severe stunted) at
# 6, 12, 18, and 24 mo of age

# Cumulative incidence pooled using random effects
#-----------------------------------
rm(list=ls())
library(tidyverse)
library(metafor)
theme_set(theme_bw())
cbPalette <- c( overall="#56B4E9", strata="#999999" , pooled="#f7a809", pooled_unstrat="#009E73") #, #f7a809, "#56B4E9",  "#E69F00",)

# load random effects function
source("U:/Scripts/Stunting/2-analyses/0_randomeffects.R")

load("U:/Data/Stunting/stunting_data_6A.RData") 
d <- d%>% filter(country=="INDIA")

# define age windows
d = d %>% 
  mutate(agecat=ifelse(agedays<=6*30.4167,"6 months",
                       ifelse(agedays>6*30.4167 & agedays<=12*30.4167,"12 months",
                              ifelse(agedays>12*30.4167 & agedays<=18*30.4167,"18 months",
                                     ifelse(agedays>12*30.4167& agedays<=24*30.4167,"24 months",""))))) %>%
  mutate(agecat=factor(agecat,levels=c("6 months","12 months","18 months","24 months")))

# check age categories
d %>%
  group_by(agecat) %>%
  summarise(n=sum(!is.na(agedays)),
            min=min(agedays/30.4167),
            mean=mean(agedays/30.4167),
            max=max(agedays/30.4167))

# identify ever stunted children
evs = d %>%
  filter(!is.na(agecat)) %>%
  group_by(studyid,subjid) %>%
  arrange(studyid,subjid) %>%
  #create variable with minhaz by age category, cumulatively
  mutate(minhaz=ifelse(agecat=="6 months",min(haz[agecat=="6 months"]),
                       ifelse(agecat=="12 months",min(haz[agecat=="6 months"|agecat=="12 months"]),
                              ifelse(agecat=="18 months",min(haz[agecat=="6 months"|agecat=="12 months"|agecat=="18 months"]),
                                     min(haz))))) %>%
  # create indicator for whether the child was ever stunted
  # by age category
  group_by(studyid,subjid,agecat) %>%
  mutate(evst=ifelse(minhaz< -3,1,0)) 

# count incident cases per study by age
# exclude time points if number of measurements per age
# in a study is <50  
cuminc.data= evs%>%
  group_by(studyid,agecat) %>%
  summarise(
    nchild=length(unique(subjid)),
    nstudy=length(unique(studyid)),
    ncases=sum(evst),
    N=sum(length(evst))) %>%
  filter(N>=50)

cuminc.data

# estimate random effects, format results
ci.res=lapply(list("6 months","12 months","18 months","24 months"),function(x)
  fit.rma(data=cuminc.data,ni="N", xi="ncases",age=x))
ci.res=as.data.frame(do.call(rbind, ci.res))
ci.res[,4]=as.numeric(ci.res[,4])
ci.res = ci.res %>%
  mutate(est=est*100, lb=lb*100, ub=ub*100)
ci.res$agecat=factor(ci.res$agecat,levels=c("6 months","12 months","18 months","24 months"))
ci.res$ptest.f=sprintf("%0.0f",ci.res$est)

ci.res$stratacol<-"pooled"

# plot pooled cumulative incidence
pdf("U:/Figures/sevstunting-cuminc-pool_6A.pdf",width=8,height=8,onefile=TRUE)
ggplot(ci.res,aes(y=est,x=agecat, fill=stratacol, color=stratacol))+
  geom_point(size=4)+
  geom_linerange(aes(ymin=lb,ymax=ub),width=0.05, alpha=0.5, size = 3) +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12))  +
  scale_y_continuous(limits=c(0,100))+
  xlab("Age category")+
  ylab("Cumulative incidence per 100 children (95% CI)")+
  annotate("text",x=ci.res$agecat,y=5,label=ci.res$nmeas.f,size=3)+
  annotate("text",x=ci.res$agecat,y=1,label=ci.res$nstudy.f,size=3)+
  annotate("text",label=ci.res$ptest.f,x=ci.res$agecat,
           y=ci.res$est,hjust=-0.75,size=3)+
  ggtitle("Pooled cumulative incidence of severe stunting")
dev.off()

#--------------------------------------
# Sensitivity analysis excluding birth from CI
# only among birth cohort studies
#--------------------------------------
# birth cohorts
bc= d %>% group_by(studyid,country) %>% 
  summarise(minage=min(agedays)) %>%
  # filter(minage==1) %>%
  mutate(birthcohort=ifelse(minage==1,1,0)) %>%
  select(-minage)

# identify ever stunted children
# excluding birth in birth cohort studies
evs.sens.nobirth = left_join(d, bc, by=c("studyid","country")) %>%
  filter(birthcohort==1 & !is.na(agecat) & agedays>1) %>%
  group_by(studyid,subjid) %>%
  arrange(studyid,subjid) %>%
  #create variable with minhaz by age category, cumulatively
  mutate(minhaz=ifelse(agecat=="6 months",min(haz[agecat=="6 months"]),
                       ifelse(agecat=="12 months",min(haz[agecat=="6 months"|agecat=="12 months"]),
                              ifelse(agecat=="18 months",min(haz[agecat=="6 months"|agecat=="12 months"|agecat=="18 months"]),
                                     min(haz))))) %>%
  # create indicator for whether the child was ever stunted
  # by age category
  group_by(studyid,subjid,agecat) %>%
  mutate(evst=ifelse(minhaz< -3,1,0)) 

# count incident cases per study by age
# exclude time points if number of measurements per age
# in a study is <50  
cuminc.data.nobirth= evs.sens.nobirth%>%
  group_by(studyid,agecat) %>%
  summarise(
    nchild=length(unique(subjid)),
    nstudy=length(unique(studyid)),
    ncases=sum(evst),
    N=sum(length(evst))) %>%
  filter(N>=50)

cuminc.data.nobirth

# estimate random effects, format results
ci.res.nobirth=lapply(list("6 months","12 months","18 months","24 months"),function(x)
  fit.rma(data=cuminc.data.nobirth,ni="N", xi="ncases",age=x))
ci.res.nobirth=as.data.frame(do.call(rbind, ci.res.nobirth))
ci.res.nobirth[,4]=as.numeric(ci.res.nobirth[,4])
ci.res.nobirth = ci.res.nobirth %>%
  mutate(est=est*100, lb=lb*100, ub=ub*100)
ci.res.nobirth$agecat=factor(ci.res.nobirth$agecat,levels=c("6 months","12 months","18 months","24 months"))
ci.res.nobirth$ptest.f=sprintf("%0.0f",ci.res.nobirth$est)

ci.res.nobirth$stratacol<-"pooled"



# plot pooled cumulative incidence
pdf("U:/Figures/sevstunting-cuminc-pool-bc-nobirth_6A.pdf",width=8,height=8,onefile=TRUE)

ggplot(ci.res.nobirth,aes(y=est,x=agecat, fill=stratacol, color=stratacol))+
  geom_point(size=4)+
  geom_linerange(aes(ymin=lb,ymax=ub),width=0.05, alpha=0.5, size = 3) +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12))  +
  scale_y_continuous(limits=c(0,100))+
  xlab("Age category")+
  ylab("Cumulative incidence per 100 children (95% CI)")+
  annotate("text",x=ci.res.nobirth$agecat,y=5,label=ci.res.nobirth$nmeas.f,size=3)+
  annotate("text",x=ci.res.nobirth$agecat,y=1,label=ci.res.nobirth$nstudy.f,size=3)+
  annotate("text",label=ci.res.nobirth$ptest.f,x=ci.res.nobirth$agecat,
           y=ci.res.nobirth$est,hjust=-0.75,size=3)+
  ggtitle("Pooled cumulative incidence of severe stunting - birth cohorts only - CI excludes birth") 

dev.off()


# identify ever stunted children
# including birth in birth cohort studies
evs.sens.birth = left_join(d, bc, by=c("studyid","country")) %>%
  filter(birthcohort==1 & !is.na(agecat)) %>%
  group_by(studyid,subjid) %>%
  arrange(studyid,subjid) %>%
  #create variable with minhaz by age category, cumulatively
  mutate(minhaz=ifelse(agecat=="6 months",min(haz[agecat=="6 months"]),
                       ifelse(agecat=="12 months",min(haz[agecat=="6 months"|agecat=="12 months"]),
                              ifelse(agecat=="18 months",min(haz[agecat=="6 months"|agecat=="12 months"|agecat=="18 months"]),
                                     min(haz))))) %>%
  # create indicator for whether the child was ever stunted
  # by age category
  group_by(studyid,subjid,agecat) %>%
  mutate(evst=ifelse(minhaz< -3,1,0)) 

# count incident cases per study by age
# exclude time points if number of measurements per age
# in a study is <50  
cuminc.data.birth= evs.sens.birth%>%
  group_by(studyid,agecat) %>%
  summarise(
    nchild=length(unique(subjid)),
    nstudy=length(unique(studyid)),
    ncases=sum(evst),
    N=sum(length(evst))) %>%
  filter(N>=50)

cuminc.data.birth

# estimate random effects, format results
ci.res.birth=lapply(list("6 months","12 months","18 months","24 months"),function(x)
  fit.rma(data=cuminc.data.birth,ni="N", xi="ncases",age=x))
ci.res.birth=as.data.frame(do.call(rbind, ci.res.birth))
ci.res.birth[,4]=as.numeric(ci.res.birth[,4])
ci.res.birth = ci.res.birth %>%
  mutate(est=est*100, lb=lb*100, ub=ub*100)
ci.res.birth$agecat=factor(ci.res.birth$agecat,levels=c("6 months","12 months","18 months","24 months"))
ci.res.birth$ptest.f=sprintf("%0.0f",ci.res.birth$est)

ci.res.birth
ci.res.birth$stratacol<-"pooled"

# plot pooled cumulative incidence
pdf("U:/Figures/sevstunting-cuminc-pool-bc-birth_6A.pdf",width=8,height=8,onefile=TRUE)
ggplot(ci.res.birth,aes(y=est,x=agecat, fill=stratacol, color=stratacol))+
  geom_point(size=4)+
  geom_linerange(aes(ymin=lb,ymax=ub),width=0.05, alpha=0.5, size = 3) +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12))  +
  scale_y_continuous(limits=c(0,100))+
  xlab("Age category")+
  ylab("Cumulative incidence per 100 children (95% CI)")+
  annotate("text",x=ci.res.birth$agecat,y=5,label=ci.res.birth$nmeas.f,size=3)+
  annotate("text",x=ci.res.birth$agecat,y=1,label=ci.res.birth$nstudy.f,size=3)+
  annotate("text",label=ci.res.birth$ptest.f,x=ci.res.birth$agecat,
           y=ci.res.birth$est,hjust=-0.75,size=3)+
  ggtitle("Pooled cumulative incidence of severe stunting - birth cohorts only - CI includes birth")
dev.off()



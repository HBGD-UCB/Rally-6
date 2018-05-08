#-----------------------------------
# Stunting analysis
# Objective 1a
# Calculate incidence at
# 6, 12, 18, and 24 mo of age

# Incidence rate pooled using random effects
#-----------------------------------
rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(metafor)
theme_set(theme_bw())
cbPalette <- c( overall="#56B4E9", strata="#999999" , pooled="#f7a809", pooled_unstrat="#009E73") #, #f7a809, "#56B4E9",  "#E69F00",)

# load random effects function
source("U:/Scripts/Stunting/2-analyses/0_randomeffects.R")

load("U:/Data/Stunting/stunting_data_6A.RData")
d <- d%>% filter(country=="INDIA")


# define age windows
d = d %>% 
  mutate(agecat=ifelse(agedays==1,"Birth",
    ifelse(agedays<=6*30.4167,"6 months",
                       ifelse(agedays>6*30.4167 & agedays<=12*30.4167,"12 months",
                              ifelse(agedays>12*30.4167& agedays<=18*30.4167,"18 months",
                                ifelse(agedays>18*30.4167& agedays<=24*30.4167,"24 months","")))))) %>%
  mutate(agecat=factor(agecat,levels=c("Birth","6 months","12 months","18 months","24 months")))

# check age categories
d %>%
  group_by(agecat) %>%
  summarise(n=sum(!is.na(agedays)),
            min=min(agedays/30.4167),
            mean=mean(agedays/30.4167),
            max=max(agedays/30.4167))

# ---------------------------------------
# flag incident cases and define risk set
# ---------------------------------------
inc.prep = d %>%
  filter(!is.na(agecat) & agecat!="Birth") %>%
  group_by(studyid,subjid) %>%
  arrange(studyid,subjid,agedays) %>%
  
  # create id for measurement within person
  mutate(measid=seq_along(subjid)) %>%
  # duration between measurements 
  mutate(agedayslag=lag(agedays)) %>%
  mutate(agedayslag=ifelse(is.na(agedayslag),0,agedayslag)) %>%
  mutate(deltat=ifelse(measid==1 & agecat=="Birth",0,agedays-agedayslag)) %>%

  # create indicator for whether haz at t < haz at t-1
  mutate(hazlag=lag(haz)) %>%
  mutate(newcase=ifelse(measid==1,ifelse(haz< -2, 1,0),
                    ifelse(hazlag>= -2 & haz< -2,1,0))) %>%
  mutate(newcaselag=lag(newcase))%>%
  mutate(newcaselag=ifelse(measid==1,0,newcaselag))%>%
  mutate(cnewcaselag=cumsum(newcaselag)) %>%
  
  # create at risk variable
  mutate(atrisk=ifelse(cnewcaselag>=1,0,1)) %>%
  # create inc case variable
  mutate(inccase=ifelse(cnewcaselag>=1,0,newcase)) %>%
  
  # create delta t with half interval for row
  # with incident case assuming it occurred halfway through
  # the follow-up period
  mutate(deltat_half=deltat/2) %>%
  mutate(deltat2=ifelse(inccase==0,deltat,deltat_half)) %>%

  # create person days
  mutate(pdays=atrisk*deltat2) %>%
    
  # clean up
  select(-c(hazlag,newcase,newcaselag, cnewcaselag,agedayslag,
            deltat,deltat_half))

# manually calculate incident cases, person-time at risk at each time point
inc.prep %>%
  group_by(agecat) %>%
  summarise(inc.case=sum(inccase),ptar=sum(pdays)) %>%
  mutate(cruderate=inc.case/ptar)


# count incident cases and sum person time at risk per study by age
# exclude time points if number of children per age
# in a study is <50  
inc.data = inc.prep %>%
  group_by(studyid,country,agecat) %>%
  summarise(ptar=sum(pdays),
            ncase=sum(inccase),
            nchild=length(unique(subjid)),
            nstudy=length(unique(studyid))) %>%
  filter(nchild>=50)
    
# estimate random effects, format results
ir.res=lapply(list("6 months","12 months","18 months","24 months"),function(x)
  fit.rma(data=inc.data,ni="ptar", xi="ncase",age=x))
ir.res=as.data.frame(do.call(rbind, ir.res))
ir.res[,4]=as.numeric(ir.res[,4])
ir.res$agecat=factor(ir.res$agecat,levels=
        c("6 months","12 months","18 months","24 months"))

ir.res
ir.res$pt.f=paste0("N=",format(ir.res$nmeas,big.mark=",",scientific=FALSE),
                  " person-days")
ir.res$ptest.f=sprintf("%0.02f",ir.res$est*1000)
ir.res$stratacol<-"pooled"

pdf("U:/Figures/stunting-inc-pool_6A.pdf",width=8,height=8,onefile=TRUE)
ggplot(ir.res,aes(y=est*1000,x=agecat, fill=stratacol, color=stratacol))+
  geom_point(size=4)+
  geom_linerange(aes(ymin=lb*1000,ymax=ub*1000),width=0.05, alpha=0.5, size = 3) +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12))  +
  xlab("Age category")+
  ylab("Incidence rate per 1,000 child-days (95% CI)")+
  scale_y_continuous(limits=c(0,5.5))+
  annotate("text",x=ir.res$agecat,y=0.2,label=ir.res$pt.f,size=3)+
  annotate("text",x=ir.res$agecat,y=0.01,label=ir.res$nstudy.f,size=3)+
  annotate("text",label=ir.res$ptest.f,x=ir.res$agecat,
           y=ir.res$est*1000,hjust=-0.3,size=3)+
  ggtitle("Pooled stunting incidence rate")
dev.off()




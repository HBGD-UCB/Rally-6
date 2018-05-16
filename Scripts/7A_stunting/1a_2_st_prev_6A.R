#-----------------------------------
# Stunting analysis
# Objective 1a
# Calculate point prevalence at
# Birth, 6, 12, 18, and 24 mo of age

# Prevalence pooled using random effects
#-----------------------------------
rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(binom)
library(metafor)
theme_set(theme_bw())
cbPalette <- c( overall="#56B4E9", strata="#999999" , pooled="#f7a809", pooled_unstrat="#009E73") #, #f7a809, "#56B4E9",  "#E69F00",)


# load meta-analysis functions
source("U:/Scripts/Stunting/2-analyses/0_randomeffects.R")

fit.rma <- function(data,age,ni,xi){
  data=filter(data,agecat==age)
  fit<-rma(ni=data[[ni]], xi=data[[xi]], 
           method="REML", measure="PR")
  out=data %>%
    ungroup() %>%
    summarise(nstudies=length(unique(studyid)),
              nmeas=sum(data[[ni]][agecat==age])) %>%
    mutate(agecat=age,est=fit$beta, se=fit$se, lb=fit$ci.lb, ub=fit$ci.ub,
           nmeas.f=paste0("N=",format(sum(data[[ni]]),big.mark=",",scientific=FALSE),
                          " obs"),
           nstudy.f=paste0("N=",nstudies," studies"))
  return(out)
}


load("U:/Data/Stunting/stunting_data_6A.RData")


# define age windows
d = d %>% 
  arrange(studyid,subjid,agedays) %>%
  mutate(agecat=ifelse(agedays==1,"Birth",
      ifelse(agedays>5*30.4167 & agedays<7*30.4167,"6 months",
       ifelse(agedays>11*30.4167 & agedays<13*30.4167,"12 months",
              ifelse(agedays>17*30.4167 & agedays<19*30.4167,"18 months",
                     ifelse(agedays>23*30.4167& agedays<25*30.4167,"24 months","")))))) %>%
    mutate(agecat=factor(agecat,levels=c("Birth","6 months",
                                         "12 months","18 months","24 months"))) %>%
    mutate(stunted=ifelse(haz< -2, 1,0),sstunted=ifelse(haz< -3, 1,0))

# check age categories
d %>%
  group_by(agecat) %>%
  summarise(n=sum(!is.na(agedays)),
            min=min(agedays/30.4167),
            mean=mean(agedays/30.4167),
            max=max(agedays/30.4167))

# count measurements per study by age
# exclude time points if number of measurements per age
# in a study is <50
prev.data = d %>%
  filter(!is.na(agecat)) %>%
  group_by(studyid,agecat) %>%
  summarise(nmeas=sum(!is.na(haz)),
            prev=mean(stunted),
            nxprev=sum(stunted==1)) %>%
  filter(nmeas>=50) 
  
#set width and height
w <- 6
h <- 4


# estimate random effects, format results
prev.res=lapply(list("Birth","6 months","12 months","18 months","24 months"),function(x) 
  fit.rma(prev.data,ni="nmeas", xi="nxprev",age=x))
prev.res=as.data.frame(do.call(rbind, prev.res))
prev.res[,4]=as.numeric(prev.res[,4])
                prev.res = prev.res %>%
  mutate(est=est*100,lb=lb*100,ub=ub*100)
prev.res$agecat=factor(prev.res$agecat,levels=c("Birth","6 months","12 months","18 months","24 months"))
prev.res$ptest.f=sprintf("%0.0f",prev.res$est)
prev.res$stratacol<-"pooled"
# plot prevalence
pdf("U:/Figures/stunting-ptprev-pool_6A.pdf",width=w,height=h,onefile=TRUE)
ggplot(prev.res,aes(y=est,x=agecat, fill=stratacol, color=stratacol))+
  geom_point(size=4)+
  geom_linerange(aes(ymin=lb,ymax=ub),width=0.05, alpha=0.5, size = 3) +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12))  +
  #scale_color_manual(values=tableau10)+
  xlab("Age category")+
  ylab("Point prevalence (95% CI)")+
  scale_y_continuous(limits=c(-4,70))+
  annotate("text",x=prev.res$agecat,y=0,label=prev.res$nmeas.f,size=3)+
  annotate("text",x=prev.res$agecat,y=-4,label=prev.res$nstudy.f,size=3)+
  annotate("text",label=prev.res$ptest.f,x=prev.res$agecat,
           y=prev.res$est,hjust=-0.75,size=3)+
  ggtitle("Pooled point prevalence of stunting")
dev.off()



#Severe stunting
# count measurements per study by age
# exclude time points if number of measurements per age
# in a study is <50
sevprev.data = d %>%
  filter(!is.na(agecat)) %>%
  group_by(studyid,agecat) %>%
  summarise(nmeas=sum(!is.na(haz)),
            prev=mean(sstunted),
            nxprev=sum(sstunted==1)) %>%
  filter(nmeas>=50) 


# estimate random effects, format results
prev.res=lapply(list("Birth","6 months","12 months","18 months","24 months"),function(x) 
  fit.rma(sevprev.data,ni="nmeas", xi="nxprev",age=x))
prev.res=as.data.frame(do.call(rbind, prev.res))
prev.res[,4]=as.numeric(prev.res[,4])
prev.res = prev.res %>%
  mutate(est=est*100,lb=lb*100,ub=ub*100)
prev.res$agecat=factor(prev.res$agecat,levels=c("Birth","6 months","12 months","18 months","24 months"))
prev.res$ptest.f=sprintf("%0.0f",prev.res$est)
prev.res$stratacol<-"pooled"
# plot prevalence
pdf("U:/Figures/sevstunting-ptprev-pool_6A.pdf",width=w,height=h,onefile=TRUE)
ggplot(prev.res,aes(y=est,x=agecat, fill=stratacol, color=stratacol))+
  geom_point(size=4)+
  geom_linerange(aes(ymin=lb,ymax=ub),width=0.05, alpha=0.5, size = 3) +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12))  +
  #scale_color_manual(values=tableau10)+
  xlab("Age category")+
  ylab("Point prevalence (95% CI)")+
  scale_y_continuous(limits=c(-4,70))+
  annotate("text",x=prev.res$agecat,y=0,label=prev.res$nmeas.f,size=3)+
  annotate("text",x=prev.res$agecat,y=-4,label=prev.res$nstudy.f,size=3)+
  annotate("text",label=prev.res$ptest.f,x=prev.res$agecat,
           y=prev.res$est,hjust=-0.75,size=3)+
  ggtitle("Pooled point prevalence of severe stunting")
dev.off()




#Plot prevalence for a single study


load("U:/Data/Stunting/stunting_data_6A.RData")


# define age windows
d = d %>% 
  arrange(studyid,subjid,agedays) %>%
  mutate(agecat=ifelse(agedays==1,"Birth",
                       ifelse(agedays>5*30.4167 & agedays<7*30.4167,"6 months",
                              ifelse(agedays>11*30.4167 & agedays<13*30.4167,"12 months",
                                     ifelse(agedays>17*30.4167 & agedays<19*30.4167,"18 months",
                                            ifelse(agedays>23*30.4167& agedays<25*30.4167,"24 months","")))))) %>%
  mutate(agecat=factor(agecat,levels=c("Birth","6 months",
                                       "12 months","18 months","24 months"))) %>%
  mutate(stunted=ifelse(haz< -2, 1,0),sstunted=ifelse(haz< -3, 1,0)) %>% 
  filter(studyid=="ki0047075b-MAL-ED")

# check age categories
d %>%
  group_by(agecat) %>%
  summarise(n=sum(!is.na(agedays)),
            min=min(agedays/30.4167),
            mean=mean(agedays/30.4167),
            max=max(agedays/30.4167))

# count measurements per study by age
# exclude time points if number of measurements per age
# in a study is <50
prev.data = d %>%
  filter(!is.na(agecat)) %>%
  group_by(studyid,agecat) %>%
  summarise(nmeas=sum(!is.na(haz)),
            prev=mean(stunted),
            nxprev=sum(stunted==1)) %>%
  filter(nmeas>=50) 

#set width and height
w <- 6
h <- 4


# estimate random effects, format results
prev.res=lapply(list("Birth","6 months","12 months","18 months","24 months"),function(x) 
  fit.rma(prev.data,ni="nmeas", xi="nxprev",age=x))
prev.res=as.data.frame(do.call(rbind, prev.res))
prev.res[,4]=as.numeric(prev.res[,4])
prev.res = prev.res %>%
  mutate(est=est*100,lb=lb*100,ub=ub*100)
prev.res$agecat=factor(prev.res$agecat,levels=c("Birth","6 months","12 months","18 months","24 months"))
prev.res$ptest.f=sprintf("%0.0f",prev.res$est)
prev.res$stratacol<-"strata"


# plot prevalence
pdf("U:/Figures/stunting-ptprev-mled_6A.pdf",width=w,height=h,onefile=TRUE)
ggplot(prev.res,aes(y=est,x=agecat, fill=stratacol, color=stratacol))+
  geom_point(size=4)+
  geom_linerange(aes(ymin=lb,ymax=ub),width=0.05, alpha=0.5, size = 3) +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12))  +
  #scale_color_manual(values=tableau10)+
  xlab("Age category")+
  ylab("Point prevalence (95% CI)")+
  scale_y_continuous(limits=c(-4,70))+
  annotate("text",x=prev.res$agecat,y=0,label=prev.res$nmeas.f,size=3)+
  annotate("text",x=prev.res$agecat,y=-4,label=prev.res$nstudy.f,size=3)+
  annotate("text",label=prev.res$ptest.f,x=prev.res$agecat,
           y=prev.res$est,hjust=-0.75,size=3)+
  ggtitle("Point prevalence of stunting - Mal-ED India Cohort")
dev.off()



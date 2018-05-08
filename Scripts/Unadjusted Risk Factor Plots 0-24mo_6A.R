




#----------------------------------------------------
# Setup
#----------------------------------------------------

rm(list=ls())
library(tidyverse)
library("ggplot2")
library("ggthemes")
library("metafor")

#Plot themes
theme_set(theme_bw())

#Forest plot Color palette
cbPalette <- c( "#56B4E9" , rep("#999999",40))

#hbgdki pallet
tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728",
  "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")


#Load data
load("C:/Users/andre/Documents/HBGDki/Results/Risk_Factor_Results_wastinc0-24_noBirthInc.Rdata") %>% filter(COUNTRY=="INDIA")
noBW<-wastinc_024_unadj
load("C:/Users/andre/Documents/HBGDki/Results/Risk_Factor_Results_wastinc0-24.Rdata") %>% filter(COUNTRY=="INDIA")

#Merge in the no-birth incidence risk factors
names(noBW)
wastinc_024_unadj$birthorder <- noBW$birthorder
wastinc_024_unadj$homedelivery <- noBW$homedelivery
wastinc_024_unadj$vagbirth <- noBW$vagbirth
wastinc_024_unadj$breastfeeding <- noBW$breastfeeding
wastinc_024_unadj$nchild5 <- noBW$nchild5
wastinc_024_unadj$ncomp <- noBW$ncomp
wastinc_024_unadj$nroom <- noBW$nroom
wastinc_024_unadj$chicken <- noBW$chicken
wastinc_024_unadj$cow <- noBW$cow
wastinc_024_unadj$improved.floor <- noBW$improved.floor
wastinc_024_unadj$improved.sanitation <- noBW$improved.sanitation
wastinc_024_unadj$safe.water <- noBW$safe.water
wastinc_024_unadj$treat.water <- noBW$treat.water
wastinc_024_unadj$SOAP <- noBW$SOAP
wastinc_024_unadj$cleancook <- noBW$cleancook
wastinc_024_unadj$enrolstunt <- noBW$enrolstunt
wastinc_024_unadj$BIRTHLEN <- noBW$BIRTHLEN
wastinc_024_unadj$BIRTHWT <- noBW$BIRTHWT
wastinc_024_unadj$GAGEBRTH <- noBW$GAGEBRTH


#Drop trials with only yearly measurements
for(i in 1:length(wastinc_024_unadj)){
  wastinc_024_unadj[[i]] <- wastinc_024_unadj[[i]] %>% filter(
    STUDYID != "ki1000110-WASH-Bangladesh" &
    STUDYID != "ki1000111-WASH-Kenya" &
    STUDYID != "kiGH5241-JiVitA-3" &
    STUDYID != "ki1148112-iLiNS-DOSE" &
    STUDYID != "ki1148112-iLiNS-DYAD-M" &
    STUDYID != "ki1000304-ZnMort" )
}





#----------------------------------------------------
# Function to create fixed or random effects pooled RRs
#----------------------------------------------------

#Fit RE and FE models
meta_fun <- function(res, method, reference=NULL){
  
  #ID reference level
  levels<-levels(res$level)
  res<-res[!is.na(res$logRR.var),]
  if(is.null(reference)){
    reference<- levels[!(levels %in% unique(res$level))]
  }
  require(metafor)
  RMAest<-data.frame(study="Pooled estimate", res$variable[1],reference ,b=NA, se=NA)
  colnames(RMAest)<-c("study","variable","level","logRR.psi","logSE")

  for(j in 1:length(unique(res$level))){
    temp<-res[res$level==unique(res$level)[j],]
    
    fit<-NULL
    tryCatch({fit<-rma(yi=logRR.psi, vi=logRR.var, data=temp, method=method, measure="RR")},
           error=function(e){cat("ERROR : REML did not converge, trying ML \n")})
    if(is.null(fit)){fit<-rma(yi=logRR.psi, vi=logRR.var, data=temp, method="ML", measure="RR")}
    
  
  est<-data.frame(study="Pooled estimate", temp$variable[1],temp$level[1] ,fit$b, fit$se)
  colnames(est)<-c("study","variable","level","logRR.psi","logSE")
  RMAest<-rbind(RMAest, est)
  }
  
RMAest$RR<-exp(RMAest$logRR)
RMAest$RR.CI1<-exp(RMAest$logRR - 1.96 * RMAest$logSE)
RMAest$RR.CI2<-exp(RMAest$logRR + 1.96 * RMAest$logSE)

#rename SE column to var just to allow binding to original results:
#not used again in calculations
colnames(RMAest)[5]<-"logRR.var"

RMAest$level.order <- 1
if("Q1" %in% RMAest$level){
  RMAest$level.order <- c(3, 1, 2, 4)
}

  return(RMAest)
}





#----------------------------------------------------
# Clean results dfs function
#----------------------------------------------------

cleandf <- function(d, meta_method="REML", RF_levels=c("<=2600","2600-3000","3000-3400",">3400")){
  
  
  d$study <- d$STUDYID
  d$country <- d$COUNTRY
  d$level <- factor(d$level)
  
  #remove grant identifier
  d$study<- gsub("^k.*?-" , "", d$study)
  
  #create labels for country-specific cohorts
  d$country <- as.character(d$country)
  d$country[d$country=="TANZANIA, UNITED REPUBLIC OF"] <- "TANZANIA"
  d$study <- paste0(d$study, " ", d$country)
  
  #drop studies with <5 cases in either reference or comparison level
  d <- d %>% filter(is.na(b) | !(b<5 | d<5))
  
  d_meantab <- d %>% subset(.,select=c(study,meanLevel:meanY, mean.CI1, mean.CI2))
  colnames(d_meantab) <- c("Study","Level", "Number in quartile", "Cumulative incidence of wasting", "95% CI-lower", "95% CI-upper")

  d <- d %>% subset(.,select=c(study, variable, level, logRR.psi, logRR.var,
                                               RR, RR.CI1, RR.CI2))
  d$level.order <- 1

  #Estimate and merge in pooled RR's
  d_RE<-meta_fun(d, method=meta_method)
  d<-rbind(d_RE,d) %>% subset(., select=-c(logRR.psi,  logRR.var))

  # Clean dataframes for plotting
  #Sort levels/studies for plot
  rankedRR <- d %>% group_by(study) %>% summarize(maxRR = max(RR, na.rm=T)) %>% arrange(maxRR)
  rankedRR$order <- 1:nrow(rankedRR)
  d <- left_join(d, rankedRR, by="study")
  d$order[d$study=="Pooled estimate"] <- 1
  
  #Drop studies with no estimated RR and order by size of RR and arrange the dataframe
  d$level <- factor(d$level, levels=RF_levels)
  d <- d %>% filter(!(maxRR=="-Inf")) %>% arrange(order, level)
  d$x<-1:nrow(d)
  d$study <- factor(d$study , as.character(unique(d$study )))

  return(d)
}




#----------------------------------------------------
# RR Plot function
#----------------------------------------------------
RRplot_fun <- function(d, reflevel, title, units="", levels=NULL, free_Y=T){
  
  #change names of risk factor levels
  if(!is.null(levels)){
    d$level <- factor(as.numeric(d$level))
    levels(d$level) <- levels
  }
  
  #plot
n <- nrow(d)

plot_df <- d %>% filter(level!=reflevel)
plot_df$size<-ifelse(plot_df$study=="Pooled estimate",1,0.5)

RRplot<-ggplot(data=plot_df) + 
  labs(title = d$variable, x = "Cohort", y = paste0("Risk ratio (reference = ",reflevel," ", units,")")) +
  geom_hline(yintercept = 1) +
  scale_y_continuous(breaks=c(0.125,0.25,0.5,1,2,4), trans='log10') +
  coord_flip() +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
    scale_size_continuous(range = c(0.5, 1))+
  geom_pointrange( mapping=aes(x=study, y=RR, ymin=RR.CI1, ymax=RR.CI2, colour=study, size=size)) +
    theme(panel.border = element_blank(), 
    strip.background = element_blank())

if(free_Y){
  RRplot<-RRplot + facet_wrap(~level, ncol=1, scales = "free_y") 
}else{
  RRplot<-RRplot + facet_wrap(~level, ncol=1) 
}

RRplot<-RRplot + geom_vline(xintercept=1.5, color="grey20", linetype=2) +
   #geom_label(aes(label="Overall RR", x=1, y=6, size=.1)) + #Need to make smaller
    ggtitle(title) +theme(legend.position="none")#+ guides(colour = guide_legend(reverse=T)) 
return(RRplot)
}




#----------------------------------------------------
# Clean and plot risk factor results, saving the pooled estimates
#----------------------------------------------------

REpooled_wastinc024_unadj <- NULL

setwd("C:/Users/andre/Documents/HBGDki/Figures/6A")

pdf("Unadjusted Risk Factor Plots Wasting CI 0-24 6A.pdf", height=10, width=8)

for(i in 1:length(wastinc_024_unadj)){
  if(nrow(wastinc_024_unadj[[i]])>0 & 
     sum(is.na(wastinc_024_unadj[[i]]$RR),na.rm=T) != nrow(wastinc_024_unadj[[i]])
    ){
    
    if(length(unique(wastinc_024_unadj[[i]]$meanLevel))>13){
    plotlevels <- c("Q3", "Q1", "Q2", "Q4")  
    wastinc_024_unadj[[i]]$meanLevel <- wastinc_024_unadj[[i]]$level <- rep(plotlevels, nrow(wastinc_024_unadj[[i]])/4)
    wastinc_024_unadj[[i]]$level.order <- rep(c(3,1,2,4), nrow(wastinc_024_unadj[[i]])/4)
    }else{
    plotlevels <- unique(wastinc_024_unadj[[i]]$meanLevel)  
    }
                 
  Anywast <- cleandf(wastinc_024_unadj[[i]], RF_levels=plotlevels)   #c("<=2700","2700-3000","3000-3400",">3400"))
  Anywast_plot <-RRplot_fun(Anywast, 
                               reflevel=wastinc_024_unadj[[i]]$level[1], 
                               title=paste0("Risk factor: ",wastinc_024_unadj[[i]]$variable[1],"\nCumulative incidence ratios:\nAny wasting from 0-24 months age"), 
                               units="",
                               levels = unique(wastinc_024_unadj[[i]]$meanLevel))
  print(Anywast_plot)
  REpooled_wastinc024_unadj<-rbind(REpooled_wastinc024_unadj, Anywast[Anywast$study=="Pooled estimate",])
  }
}

dev.off()


save(REpooled_wastinc024_unadj, file="C:/Users/andre/Documents/HBGDki/Results/REpooled_wastinc024_unadj_6A.Rdata")







load("C:/Users/andre/Documents/HBGDki/Results/REpooled_wastinc024_unadj_6A.Rdata")



#Drop variables that don't make sense to pool
d <- REpooled_wastinc024_unadj %>% filter(variable!="month" & variable!="birthmonth")

#filter out birthweight
#d <- d %>% filter(variable!="BIRTHWT")

d$reflabel <- ""
d$reflabel[is.na(d$RR)] <- "(ref.)"
d$RR[is.na(d$RR)] <- 1


#FIX to offset (ref.) by a relative amount
d <- d %>% group_by(variable) %>% mutate(RRrange=diff(range(RR), na.rm=T))


#Pooled plot function
scaleFUN <- function(x) sprintf("%.2f", x)

RF_metaplot <- function(d, title="", yticks=c(0.5, 0.6, 0.7,0.8,0.9,1, 1/0.9, 1/0.8, 1/0.7, 1/0.6, 2)){

p<-ggplot(d, aes(x=level)) + 
  geom_point(aes(y=RR, fill=variable, color=variable), size = 4) +
  geom_linerange(aes(ymin=RR.CI1, ymax=RR.CI2, color=variable),
                 alpha=0.5, size = 3) +
  geom_text(aes( y=0.95, label=reflabel, colour=variable)) +
  labs(x = "Risk factor level", y = "Cumulative Incidence Ratio") +
  geom_hline(yintercept = 1) +
  coord_cartesian(ylim=c(0.7, 1/0.6)) +
  scale_y_continuous(breaks=yticks, trans='log10', labels=scaleFUN) +
  scale_fill_manual(values=rep(tableau10,4)) +
  scale_colour_manual(values=rep(tableau10,4)) +
  scale_size_continuous(range = c(0.5, 1))+
  theme(strip.background = element_blank(),
    legend.position="none",
    strip.text.x = element_text(size=12),
    axis.text.x = element_text(size=12)) +
  facet_wrap(~variable,  scales = "free_x") +
  ggtitle(title)

return(p)
}





#Group risk factors by category
unique(d$variable)
d <- d %>% arrange(variable, level.order)

#Example birthweight plot
d0 <- d %>% filter(variable=="BIRTHWT")

#Child
d1 <- d %>% filter(variable=="SEX" | variable=="birthorder" | variable=="homedelivery" | variable=="vagbirth" | variable=="breastfeeding" | 
                     variable=="nchild5" | variable=="DIARFL" | variable=="enrolstunt" | variable=="GAGEBRTH" | variable=="BIRTHLEN" | variable=="DURBRST")
#Parental
d2 <- d %>% filter(variable=="single" | variable=="MHTCM" | variable=="MWTKG" | variable=="MBMI" | variable=="MAGE" | variable=="MEDUCYRS" | variable=="FHTCM" | variable=="FAGE" | variable=="FEDUCYRS")
#House
d3 <- d %>% filter(variable=="nroom" | variable=="chicken" | variable=="cow" | variable=="improved.floor" | variable=="improved.sanitation" | variable=="safe.water" | variable=="treat.water" | variable=="SOAP" | variable=="cleancook" | variable=="HHwealth_quart")

#Drop risk factors with insufficient study, or that need more cleaning
d1 <- d1 %>% filter(variable!="DIARFL" & variable!="breastfeeding" & variable!="DURBRST")                                                                       
d2 <- d2 %>% filter(variable!="")                                                                       
d3 <- d3 %>% filter(variable!="improved.sanitation" & variable!="safe.water" & 
                      variable!="treat.water" & variable!="SOAP")                                                                       




#Rename variables and levels for plot axes

recode_hbgdki <- function(df){
  
df$level <- as.character(df$level)  
df$variable <- as.character(df$variable)  

#Arrange by magnitude of largest RR
df <- df %>% group_by(variable) %>% 
             mutate(maxRR=max(abs(1-RR), na.rm=T)) %>% 
             ungroup() 

df$level.order[df$variable=="MBMI"] <- c(4,1,2,3)

df$level.order[df$variable=="nchild5" & df$level=="0"] <- 1
df$level.order[df$variable=="nchild5" & df$level=="1"] <- 2
df$level.order[df$variable=="nchild5" & df$level=="2+"] <- 3
  
df$level[df$variable=="nchild5" & df$level=="0"] <- "Only\nchild"
df$level[df$variable=="nchild5" & df$level=="1"] <- "Two\nchildren"
df$level[df$variable=="nchild5" & df$level=="2+"] <- "3+\nchildren"

df$level.order[df$variable=="nroom" & df$level=="1"] <- 1
df$level.order[df$variable=="nroom" & df$level=="2"] <- 2
df$level.order[df$variable=="nroom" & df$level=="3"] <- 3
df$level.order[df$variable=="nroom" & df$level=="4+"] <- 4

df$level[df$variable=="nroom" & df$level=="1"] <- "One"
df$level[df$variable=="nroom" & df$level=="2"] <- "Two"
df$level[df$variable=="nroom" & df$level=="3"] <- "Three"
df$level[df$variable=="nroom" & df$level=="4+"] <- "Four +"

df$level[df$level=="Wealth Q1"] <- "Wealth\nQ1"
df$level[df$level=="Wealth Q2"] <- "Wealth\nQ2"
df$level[df$level=="Wealth Q3"] <- "Wealth\nQ3"
df$level[df$level=="Wealth Q4"] <- "Wealth\nQ4"

df$level.order[df$variable=="MBMI" & df$level==">24.093"] <- 4
df$level.order[df$variable=="MBMI" & df$level=="<=18.904"] <- 1
df$level.order[df$variable=="MBMI" & df$level=="18.904-21.05"] <- 2
df$level.order[df$variable=="MBMI" & df$level=="21.05-24.093"] <- 3

df$level[df$variable=="MBMI" & df$level==">24.093"] <- ">24"
df$level[df$variable=="MBMI" & df$level=="<=18.904"] <- "<=19"
df$level[df$variable=="MBMI" & df$level=="18.904-21.05"] <- "(19- 21]"
df$level[df$variable=="MBMI" & df$level=="21.05-24.093"] <- "(21- 24]"

df$level[df$variable=="MHTCM" & df$level=="43-51.2"] <- "(43-\n51]"
df$level[df$variable=="MHTCM" & df$level=="51.2-63"] <- "(51-\n63]"

df$level[df$variable=="MAGE" & df$level=="<=20.769"] <- "<=21"
df$level[df$variable=="MAGE" & df$level=="20.769-24"] <- "(21-24]"
df$level[df$variable=="MAGE" & df$level=="24-29"] <- "(24-29]"
   
df$level[df$variable=="FAGE" & df$level=="24-28"] <- "(24-28]"
df$level[df$variable=="FAGE" & df$level=="28-33"] <- "(28-33]"

df$level[df$variable=="FHTCM" & df$level=="168-174"] <- "(168-\n174]"
df$level[df$variable=="FHTCM" & df$level=="174-178"] <- "(174-\n178]"

df$level[df$variable=="MWTKG" & df$level=="43-51.2"] <- "(43-51]"
df$level[df$variable=="MWTKG" & df$level=="51.2-63"] <- "(51-63]"

df$level[df$variable=="BIRTHLEN" & df$level=="<=47.5"] <- "<=48"
df$level[df$variable=="BIRTHLEN" & df$level=="47.5-49.5"] <- "(48-50]"
df$level[df$variable=="BIRTHLEN" & df$level=="49.5-51.2"] <- "(50-51]"
df$level[df$variable=="BIRTHLEN" & df$level==">51.2"] <- ">51"



df <- df %>% arrange(desc(maxRR), level.order)
df$level<-factor(df$level, levels = unique(df$level))
df$variable<-factor(df$variable, levels = unique(df$variable))

df$level<-recode(df$level,
  "0"= "No",
  "1"= "Yes",
  " 0"= "No",
  " 1"= "Yes")  
  
df$variable<-recode(df$variable,
  "SEX"= "Gender",
  "birthorder"= "Birth Order",
  "homedelivery"= "Home delivery",
  "vagbirth"= "Vaginal birth",
  "single"= "Single mother",
  "nchild5"= "Number of children\nunder 5",
  "nroom"= "Number\nof rooms",
  "improved.floor"= "Improved floor",
  "cleancook"= "Clean cooking\nfuel",
  "enrolstunt"= "Stunted at\nenrollment",
  "HHwealth_quart"= "Quartile of\nhousehold asset-\nbased wealth",
  "GAGEBRTH"= "Gestational age\nat birth",
  "BIRTHLEN"= "Birth length",
  "MHTCM"= "Maternal\nheight",
  "MWTKG"= "Maternal\nweight",
  "MBMI"= "Maternal\nBMI",
  "MAGE"= "Maternal\nage",
  "MEDUCYRS"= "Maternal\neducation",
  "FHTCM"= "Paternal\nheight",
  "FAGE"= "Paternal\nage",
  "FEDUCYRS"= "Paternal\neducation",
  "chicken"="Chicken ownership",
  "cow"="Cow ownership")

  return(df)
}

d0 <- recode_hbgdki(d0)
d1 <- recode_hbgdki(d1)
d2 <- recode_hbgdki(d2)
d3 <- recode_hbgdki(d3)




#Change order of risk factor levels

d0$level <- as.character(d0$level)
d0$level <-c("Q1","Q2","Q3","Q4")
d0$level <- factor(d0$level)

p0<-ggplot(d0, aes(x=level)) + 
  geom_point(aes(y=RR, fill=variable, color=variable), size = 4) +
  geom_linerange(aes(ymin=RR.CI1, ymax=RR.CI2, color=variable),
                 alpha=0.5, size = 3) +
  geom_text(aes( y=0.95, label=reflabel, colour=variable)) +
  labs(x = "Risk factor level", y = "Cumulative Incidence Ratio") +
  geom_hline(yintercept = 1) +
  scale_y_continuous(breaks=c(0.7,0.8,0.9,1, 1/0.9, 1/0.8, 1/0.7, 1/0.6, 2), trans='log10', labels=scaleFUN) +
  scale_fill_manual(values=rep(tableau10,4)) +
  scale_colour_manual(values=rep(tableau10,4)) +
  scale_size_continuous(range = c(0.5, 1))+
  theme(strip.background = element_blank(),
    legend.position="none",
    strip.text.x = element_text(size=12),
    axis.text.x = element_text(size=12)) +
  ggtitle("Birthweight: pooled CIR, 1-24 months")

#p0

p1 <- RF_metaplot(d1, title="Child characteristics")
#p1

p2 <- RF_metaplot(d2, title="Parental characteristics")
#p2 

p3 <- RF_metaplot(d3, title="Household characteristics")
#p3

pdf("Unadjusted Risk Factor Plots Wasting CI Pooled 0-24 6A.pdf", height=9, width=9)
p0
p1
p2
p3
dev.off()



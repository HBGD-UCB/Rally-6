
# load packages
rm(list=ls())
library(tidyverse)
library(metafor)

source("~/Rally-6/Scripts/HBGDki_plotting_functions_6A.R")
source("~/Rally-6/Scripts/Meta-analysis functions_6A.R")


setwd("~/Rally-6/Results")
load("descriptive_epi_mean_monthly_cohorts_6A.Rdata")

#Drop studies that enrol acutely ill children

unique(d$STUDYID)
d <- d %>% filter(STUDYID!="ki1000301-DIVIDS" & STUDYID!="ki1000304-LBW" & STUDYID!="ki1000306-ZincSGA" & STUDYID!="ki1000304b-ZincInf")
                     
            

d<-prep_desc_data(d)


#Drop less than 25 in an age category (after it has contributed to the pooled statistic)
d <- d[d$N > 24, ]




#Set theme and colors
theme_set(theme_bw())
d$stratacol <- "strata"
d$stratacol[d$strata=="Overall"] <- "overall"
d$stratacol[d$pooled==1] <- "pooled"
d$stratacol[d$strata=="Overall" & d$pooled==1] <- "pooled_unstrat"
cbPalette <- c( overall="#56B4E9", strata="#999999" , pooled="#f7a809", pooled_unstrat="#009E73")

#Set plot width and height
w <- 10
h <- 7

#----------------------
# Primary plots
#----------------------


#Set plot directory

setwd("C:/Users/andre/Documents/Rally-6/Figures/")

#Wasting prevalence
p1 <- desc_epi_metaplot(d, stat="Prevalence\nof\nwasting",
                     ylabel="Wasting longitudinal prevalence",
                     title="Wasting longitudinal prevalence")
ggsave("WastPrev_metaplot.pdf", p1, width = w, height = h, units = "in")

#Severe wasting prevalence
p2 <- desc_epi_metaplot(d, stat="Prevalence\nof\nsevere\nwasting",
                     ylabel="Severe wasting longitudinal prevalence",
                     title="Severe wasting longitudinal prevalence")
ggsave("SevWastPrev_metaplot.pdf", p2, width = w, height = h, units = "in")


#Wasting IR
p3 <- desc_epi_metaplot(d, stat="Wasting\nincidence\nrate",
                     ylabel="Wasting incidence rate per 1000 days",
                     title="Wasting incidence rate")
ggsave("WastInc_metaplot.pdf", p3, width = w, height = h, units = "in")

#Sev wasting IR

p4 <- desc_epi_metaplot(d, stat="Severe\nwasting\nrecovery\nincidence\nrate",
                     ylabel="Severe wasting incidence rate per 1000 days",
                     title="Severe wasting incidence rate")
ggsave("SevWastInc_metaplot.pdf", p4, width = w, height = h, units = "in")

#Duration
p5<- desc_epi_metaplot(d, stat="Average\nduration\nof\nwasting",
                     ylabel="Duration of wasting (days)",
                     title="Average duration of wasting")
ggsave("Durationc_metaplot.pdf", p5, width = w, height = h, units = "in")



# 60 day recovery
p6 <- desc_epi_metaplot(d, stat="Percent\nwasting\nrecovered\nin 60 days",
                     ylabel="Percent wasting recovered in 60 days",
                     title="Percent wasting recovered in 60 days")
ggsave("WastRec60_metaplot.pdf", p6, width = w, height = h, units = "in")


#Unstratified comparison of 30,60,90 day recovery
d_unstrat<-d[d$strata=="Overall",]
d_unstrat <- d_unstrat %>% filter(statistic=="Percent\nwasting\nrecovered\nin 30 days" |
                                  statistic=="Percent\nwasting\nrecovered\nin 60 days" |
                                  statistic=="Percent\nwasting\nrecovered\nin 90 days")
d_unstrat$statistic<-as.character(d_unstrat$statistic)
d_unstrat$statistic[d_unstrat$statistic=="Percent\nwasting\nrecovered\nin 30 days"]<-"30 days"
d_unstrat$statistic[d_unstrat$statistic=="Percent\nwasting\nrecovered\nin 60 days"]<-"60 days"
d_unstrat$statistic[d_unstrat$statistic=="Percent\nwasting\nrecovered\nin 90 days"]<-"90 days"
d_unstrat$statistic<-as.factor(d_unstrat$statistic)
d_unstrat$Lower.95.CI <- d_unstrat$Lower.95.CI  * 100
d_unstrat$Mean <- d_unstrat$Mean  * 100
d_unstrat$Upper.95.CI <- d_unstrat$Upper.95.CI  * 100

p7 <- desc_epi_metaplot(d_unstrat, stat=NULL,
                     ylabel="Percent wasting recovered",
                     title="Percent wasting recovered within 30, 60, and 90 days",
                     xlabel="Recovery time")
ggsave("WastRec_unstrat_metaplot.pdf", p7, width = w, height = h, units = "in")




#Unstratified comparison of faltering into severe wasting
d_unstrat<-d[d$strata=="Overall",]
d_unstrat <- d_unstrat %>% filter(statistic=="Percent\nfalter to\nsevere\nwasting\nin 30 days" |
                                  statistic=="Percent\nfalter to\nsevere\nwasting\nin 60 days" |
                                  statistic=="Percent\nfalter to\nsevere\nwasting\nin 90 days")
d_unstrat$statistic<-as.character(d_unstrat$statistic)
d_unstrat$statistic[d_unstrat$statistic=="Percent\nfalter to\nsevere\nwasting\nin 30 days"]<-"30 days"
d_unstrat$statistic[d_unstrat$statistic=="Percent\nfalter to\nsevere\nwasting\nin 60 days"]<-"60 days"
d_unstrat$statistic[d_unstrat$statistic=="Percent\nfalter to\nsevere\nwasting\nin 90 days"]<-"90 days"
d_unstrat$statistic<-as.factor(d_unstrat$statistic)
d_unstrat$Lower.95.CI <- d_unstrat$Lower.95.CI  * 100
d_unstrat$Mean <- d_unstrat$Mean  * 100
d_unstrat$Upper.95.CI <- d_unstrat$Upper.95.CI  * 100

d_unstrat$Mean[d_unstrat$country_cohort=="Content Peru" | d_unstrat$country_cohort=="Mal-ED Nepal" | d_unstrat$country_cohort=="Mal-ED Tanz." | d_unstrat$country_cohort=="Mal-ED Brazil" | d_unstrat$country_cohort=="Mal-ED Peru"]<-NA
d_unstrat$Lower.95.CI[d_unstrat$country_cohort=="Content Peru" | d_unstrat$country_cohort=="Mal-ED Nepal" | d_unstrat$country_cohort=="Mal-ED Tanz." | d_unstrat$country_cohort=="Mal-ED Brazil" | d_unstrat$country_cohort=="Mal-ED Peru"]<-NA
d_unstrat$Upper.95.CI[d_unstrat$country_cohort=="Content Peru" | d_unstrat$country_cohort=="Mal-ED Nepal" | d_unstrat$country_cohort=="Mal-ED Tanz." | d_unstrat$country_cohort=="Mal-ED Brazil" | d_unstrat$country_cohort=="Mal-ED Peru"]<-NA

p8 <- desc_epi_metaplot(d_unstrat, stat=NULL,
                     ylabel="Percent wasting faltered to severe wasting",
                     title="Percent wasting faltered to severe wasting within 30, 60, and 90 days",
                     xlabel="Faltering time")
ggsave("WastFalter_unstrat_metaplot.pdf", p8, width = w, height = h, units = "in")




#Print all descriptive plots together
pdf("6A_Descriptive_epi_plots.pdf", width=w,height=h, paper="USr")
p1
p2
p3
p4
p5
p6
p7
p8
dev.off()


#Save plots
save(p1,
      p2,
      p3,
      p4,
      p5,
      p6,
      p7,
      p8, 
file="C:/Users/andre/Documents/Rally-6/Results/6A_Descriptive_epi_plots.Rdata")



#-------------------------------------------
# Pooled results only for sprint report out
#-------------------------------------------


dpooled <- d %>% filter(country_cohort=="Pooled")

head(dpooled)





#Order the statistics
unique(dpooled$statistic)

#Pull out recovery/faltering plots:
d_rec<-dpooled[grepl("Percent", dpooled$statistic),]
dpooled<-dpooled[!grepl("Percent", dpooled$statistic),]

dpooled$order<-NA
dpooled$order[dpooled$statistic=="Prevalence\nof\nwasting"] <-1
dpooled$order[dpooled$statistic=="Prevalence\nof\nsevere\nwasting"] <-2
dpooled$order[dpooled$statistic=="Wasting\nincidence\nrate"] <-3
dpooled$order[dpooled$statistic=="Severe\nwasting\nincidence\nrate"] <-4
dpooled$order[dpooled$statistic=="Average\nduration\nof\nwasting"] <-5
dpooled$order[dpooled$statistic=="Wasting\nrecovery\nincidence\nrate"] <-6
dpooled$order[dpooled$statistic=="Severe\nwasting\nrecovery\nincidence\nrate"] <-7
dpooled <- dpooled %>% arrange(order)
dpooled$statistic <- factor(dpooled$statistic, levels=unique(dpooled$statistic))

#Convert prevalences to precent
dpooled$Mean[dpooled$order<3] <- dpooled$Mean[dpooled$order<3] * 100
dpooled$Lower.95.CI[dpooled$order<3] <- dpooled$Lower.95.CI[dpooled$order<3] * 100
dpooled$Upper.95.CI[dpooled$order<3] <- dpooled$Upper.95.CI[dpooled$order<3] * 100

      
#Drop recovery incidence rates
dpooled <- dpooled %>% filter(order < 6)

  p <- ggplot(dpooled) +
              geom_point(aes(x=strata, y=Mean, fill=stratacol, color=stratacol), size = 4) +
              geom_linerange(aes(x=strata, ymin = Lower.95.CI, ymax = Upper.95.CI, color=stratacol), 
                 alpha=0.5, size = 3)  + 
  facet_wrap(~statistic, nrow = 1, scales="free")  + 
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text.x = element_text(size=8),
        axis.text.x = element_text(size=8)) +
  ylab("")+
  ggtitle("Pooled, age-stratified descriptive statisics\nof child wasting from birth to 24 months in India") + 
  xlab("Child age stratification")

  print(p)
  
  ggsave("Pooled_wast.pdf", p, width = 7.5, height = 4, units = "in")

  
  
  d_rec$Mean <- d_rec$Mean * 100
  d_rec$Lower.95.CI <- d_rec$Lower.95.CI * 100
  d_rec$Upper.95.CI <- d_rec$Upper.95.CI * 100

    p_rec <- ggplot(d_rec) +
              geom_point(aes(x=strata, y=Mean, fill=stratacol, color=stratacol), size = 4) +
              geom_linerange(aes(x=strata, ymin = Lower.95.CI, ymax = Upper.95.CI, color=stratacol), 
                 alpha=0.5, size = 3)  + 
  facet_wrap(~statistic, nrow = 2)  + 
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text.x = element_text(size=8),
        axis.text.x = element_text(size=8)) +
  ylab("Percent")+
  ggtitle("Pooled, age-stratified faltering to severe wasting and recovery from wasting\namong moderated wasted children") + 
  xlab("Child age stratification")

  print(p_rec)

  ggsave("Pooled_wastrec.pdf", p_rec, width = 7.5, height = 7.5, units = "in")






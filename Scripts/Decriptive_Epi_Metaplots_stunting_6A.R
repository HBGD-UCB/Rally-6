
# load packages
rm(list=ls())
library(tidyverse)
library(metafor)

source("C:/Users/andre/Documents/HBGDki//Rally-6/Scripts/HBGDki_plotting_functions_6A.R")
source("C:/Users/andre/Documents/HBGDki//Rally-6/Scripts/Meta-analysis functions_6A.R")


setwd("C:/Users/andre/Documents/HBGDki//Rally-6/Results")
load("descriptive_epi_mean_monthly_cohorts_stunting_6A.Rdata")


#Drop out acutely ill cohorts
unique(d$STUDYID)

d <- d %>% filter(STUDYID!="ki1000301-DIVIDS" & STUDYID!="ki1000304b-ZincInf" & 
                  STUDYID!="ki1000306-ZincSGA" & STUDYID!="ki1000304-LBW")
d <- droplevels(d)
d<-prep_desc_data(d)


#Drop less than 25 in an age category (after it has contributed to the pooled statistic)
d <- d[d$N > 24, ]
#Drop na catergory
d <- d[!is.na(d$country_cohort),]


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

setwd("C:/Users/andre/Documents/HBGDki/Rally-6/Figures/")

 # d$Mean <- d$Mean * 100
 #  d$Lower.95.CI <- d$Lower.95.CI * 100
 #  d$Upper.95.CI <- d$Upper.95.CI * 100

#Stunting prevalence
st_p1 <- desc_epi_metaplot(d, stat="Prevalence\nof\nstunting",
                     ylabel="Stunting longitudinal prevalence",
                     title="Stunting longitudinal prevalence")
ggsave("StuntPrev_metaplot.png", st_p1, width = w, height = h, units = "in")
st_p1

#Severe stunting prevalence
st_p2 <- desc_epi_metaplot(d, stat="Prevalence\nof\nsevere\nstunting",
                     ylabel="Severe stunting longitudinal prevalence",
                     title="Severe stunting longitudinal prevalence")
ggsave("SevStuntPrev_metaplot.png", st_p2, width = w, height = h, units = "in")


#Stunting IR
st_p3 <- desc_epi_metaplot(d, stat="Stunting\nincidence\nrate",
                     ylabel="Stunting incidence rate per 1000 days",
                     title="Stunting incidence rate")
ggsave("StuntInc_metaplot.png", st_p3, width = w, height = h, units = "in")

#Sev stunting IR

st_p4 <- desc_epi_metaplot(d, stat="Severe\nstunting\nrecovery\nincidence\nrate",
                     ylabel="Severe stunting incidence rate per 1000 days",
                     title="Severe stunting incidence rate")
ggsave("SevStuntInc_metaplot.png", st_p4, width = w, height = h, units = "in")

#Duration
st_p5<- desc_epi_metaplot(d, stat="Average\nduration\nof\nstunting",
                     ylabel="Duration of stunting (days)",
                     title="Average duration of stunting")
ggsave("Durationc_metaplot.png", st_p5, width = w, height = h, units = "in")



# 60 day recovery
st_p6 <- desc_epi_metaplot(d, stat="Percent\nstunting\nrecovered\nin 60 days",
                     ylabel="Percent stunting recovered in 60 days",
                     title="Percent stunting recovered in 60 days")
ggsave("StuntRec60_metaplot.png", st_p6, width = w, height = h, units = "in")


#Unstratified comparison of 30,60,90 day recovery
d_unstrat<-d[d$strata=="Overall",]
d_unstrat <- d_unstrat %>% filter(statistic=="Percent\nstunting\nrecovered\nin 30 days" |
                                  statistic=="Percent\nstunting\nrecovered\nin 60 days" |
                                  statistic=="Percent\nstunting\nrecovered\nin 90 days")
d_unstrat$statistic<-as.character(d_unstrat$statistic)
d_unstrat$statistic[d_unstrat$statistic=="Percent\nstunting\nrecovered\nin 30 days"]<-"30 days"
d_unstrat$statistic[d_unstrat$statistic=="Percent\nstunting\nrecovered\nin 60 days"]<-"60 days"
d_unstrat$statistic[d_unstrat$statistic=="Percent\nstunting\nrecovered\nin 90 days"]<-"90 days"
d_unstrat$statistic<-as.factor(d_unstrat$statistic)
d_unstrat$Lower.95.CI <- d_unstrat$Lower.95.CI  * 100
d_unstrat$Mean <- d_unstrat$Mean  * 100
d_unstrat$Upper.95.CI <- d_unstrat$Upper.95.CI  * 100

st_p7 <- desc_epi_metaplot(d_unstrat, stat=NULL,
                     ylabel="Percent stunting recovered",
                     title="Percent stunting recovered within 30, 60, and 90 days",
                     xlabel="Recovery time")
ggsave("StuntRec_unstrat_metaplot.png", st_p7, width = w, height = h, units = "in")




#Unstratified comparison of faltering into severe stunting
d_unstrat<-d[d$strata=="Overall",]
d_unstrat <- d_unstrat %>% filter(statistic=="Percent\nfalter to\nsevere\nstunting\nin 30 days" |
                                  statistic=="Percent\nfalter to\nsevere\nstunting\nin 60 days" |
                                  statistic=="Percent\nfalter to\nsevere\nstunting\nin 90 days")
d_unstrat$statistic<-as.character(d_unstrat$statistic)
d_unstrat$statistic[d_unstrat$statistic=="Percent\nfalter to\nsevere\nstunting\nin 30 days"]<-"30 days"
d_unstrat$statistic[d_unstrat$statistic=="Percent\nfalter to\nsevere\nstunting\nin 60 days"]<-"60 days"
d_unstrat$statistic[d_unstrat$statistic=="Percent\nfalter to\nsevere\nstunting\nin 90 days"]<-"90 days"
d_unstrat$statistic<-as.factor(d_unstrat$statistic)
d_unstrat$Lower.95.CI <- d_unstrat$Lower.95.CI  * 100
d_unstrat$Mean <- d_unstrat$Mean  * 100
d_unstrat$Upper.95.CI <- d_unstrat$Upper.95.CI  * 100

d_unstrat$Mean[d_unstrat$country_cohort=="Content Peru" | d_unstrat$country_cohort=="Mal-ED Nepal" | d_unstrat$country_cohort=="Mal-ED Tanz." | d_unstrat$country_cohort=="Mal-ED Brazil" | d_unstrat$country_cohort=="Mal-ED Peru"]<-NA
d_unstrat$Lower.95.CI[d_unstrat$country_cohort=="Content Peru" | d_unstrat$country_cohort=="Mal-ED Nepal" | d_unstrat$country_cohort=="Mal-ED Tanz." | d_unstrat$country_cohort=="Mal-ED Brazil" | d_unstrat$country_cohort=="Mal-ED Peru"]<-NA
d_unstrat$Upper.95.CI[d_unstrat$country_cohort=="Content Peru" | d_unstrat$country_cohort=="Mal-ED Nepal" | d_unstrat$country_cohort=="Mal-ED Tanz." | d_unstrat$country_cohort=="Mal-ED Brazil" | d_unstrat$country_cohort=="Mal-ED Peru"]<-NA

st_p8 <- desc_epi_metaplot(d_unstrat, stat=NULL,
                     ylabel="Percent stunting faltered to severe stunting",
                     title="Percent stunting faltered to severe stunting within 30, 60, and 90 days",
                     xlabel="Faltering time")
ggsave("StuntFalter_unstrat_metaplot.png", st_p8, width = w, height = h, units = "in")




#Print all descriptive plots together
pdf("6A_Descriptive_epi_plots_stunting.png", width=w,height=h, paper="USr")
st_p1
st_p2
st_p3
st_p4
st_p5
st_p6
st_p7
st_p8
dev.off()


#Save plots
save(st_p1,
      st_p2,
      st_p3,
      st_p4,
      st_p5,
      st_p6,
      st_p7,
      st_p8, 
file="C:/Users/andre/Documents/Rally-6/Results/6A_Descriptive_epi_plots_stunting.Rdata")












#-------------------------------------------
# Pooled results only for sprint report out
#-------------------------------------------

dpooled <- d %>% filter(country_cohort=="Pooled")

head(dpooled)


#Order the statistics
unique(dpooled$statistic)

#Pull out recovery/faltering plots:
dpooled<-dpooled[grepl("Prev", dpooled$statistic),]

dpooled$statistic <- relevel(dpooled$statistic, ref="Prevalence\nof\nstunting")
  
  dpooled$Mean <- dpooled$Mean * 100
  dpooled$Lower.95.CI <- dpooled$Lower.95.CI * 100
  dpooled$Upper.95.CI <- dpooled$Upper.95.CI * 100
  st_pooled <- ggplot(dpooled) +
              geom_point(aes(x=strata, y=Mean, fill=stratacol, color=stratacol), size = 4) +
              geom_linerange(aes(x=strata, ymin = Lower.95.CI, ymax = Upper.95.CI, color=stratacol), 
                 alpha=0.5, size = 3)  + 
  facet_wrap(~statistic, nrow = 1)  + 
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text.x = element_text(size=8),
        axis.text.x = element_text(size=8)) +
  ylab("Prevalence")+
  ggtitle("Pooled, age-stratified prevalence\nof child stunting from birth to 24 months in Indian cohorts") + 
  xlab("Child age stratification")


  ggsave("Pooled_stunt_prev.png", st_pooled, width = 7.5, height = 5, units = "in")

  
  
  
  #Save plots
save(st_p1,
      st_p2,
      st_p3,
      st_p4,
      st_p5,
      st_p6,
      st_p7,
      st_p8, 
     st_pooled,
file="C:/Users/andre/Documents/HBGDki/Rally-6/Results/6A_Descriptive_epi_plots_stunting.Rdata")














#------------------------------------
# Secondary plots (for presentation)
#------------------------------------




#get legend
df <- d[d$statistic=="Stunting\nincidence\nrate",]
df <- df %>% rename(Legend = stratacol)
df$Legend <- factor(df$Legend, levels=unique(df$Legend))
df$Legend<-recode(df$Legend,
"overall"= "Unstratified",
"pooled"= "Pooled age stratified",
 "pooled_unstrat"= "Pooled unstratified",
 "strata"= "Age stratified"
)
col_legend <- c( `Unstratified`="#56B4E9", `Age stratified`="#999999" , `Pooled age stratified`="#f7a809", `Pooled unstratified`="#009E73") #, #f7a809, "#56B4E9",  "#E69F00",)
cbPalette <- c( overall="#56B4E9", strata="#999999" , pooled="#f7a809", pooled_unstrat="#009E73") #, #f7a809, "#56B4E9",  "#E69F00",)

ggplot(df, aes(`Child age stratification`)) +
  geom_point(aes(x=strata, y=Mean, fill=Legend, color=Legend), size = 4) +
  geom_linerange(aes(x=strata, ymin = Lower.95.CI, ymax = Upper.95.CI, color=Legend),
                 alpha=0.5, size = 3) +
  scale_fill_manual(values=col_legend) +
  scale_colour_manual(values=col_legend) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position="right",
        strip.text.x = element_text(size=8),
        axis.text.x = element_text(size=8)) +
  ylab("Stunting incidence rate per 1000 days")+
  ggtitle("Stunting incidence rate")

#Plot single pooled estimate
ggplot(d[d$statistic=="Prevalence\nof\nstunting" & d$country_cohort=="Pooled",], aes(`Child age stratification`)) +
  geom_point(aes(x=strata, y=Mean*100, fill=stratacol, color=stratacol), size = 4) +
  geom_linerange(aes(x=strata, ymin = Lower.95.CI*100, ymax = Upper.95.CI*100, color=stratacol),
                 alpha=0.5, size = 3) +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12)) +
  ylab("Stunting prevalence")+
  ggtitle("Pooled stunting longitudinal prevalence")

#  [1] Pooled         CMC-V-BCS-2002 EU             IRC            TDC           
#  [6] COHORTS        MAL-ED         VITAMIN-A      Vitamin-B12    SAS-CompFeed  
# [11] SAS-FoodSuppl  ZnMort  

ggplot(d[d$statistic=="Prevalence\nof\nstunting" & d$country_cohort=="MAL-ED",], aes(`Child age stratification`)) +
  geom_point(aes(x=strata, y=Mean*100, fill=stratacol, color=stratacol), size = 4) +
  geom_linerange(aes(x=strata, ymin = Lower.95.CI*100, ymax = Upper.95.CI*100, color=stratacol),
                 alpha=0.5, size = 3) +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text.x = element_text(size=12),
        axis.text.x = element_text(size=12)) +
  ylab("Stunting prevalence")+
  ggtitle("MAL-ED India stunting longitudinal prevalence")

#Stunting IR
desc_epi_metaplot(d[d$country_cohort=="",], stat="Stunting\nincidence\nrate",
                  ylabel="Stunting incidence rate per 1000 days",
                  title="Stunting incidence rate")










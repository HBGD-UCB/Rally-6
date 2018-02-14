
# load packages
rm(list=ls())
library(tidyverse)
library(metafor)

source("C:/Users/andre/Documents/HBGDki/Scripts/HBGDki_plotting_functions.R")
source("C:/Users/andre/Documents/HBGDki/Scripts/Meta-analysis functions.R")


setwd("C:/Users/andre/Documents/HBGDki/Results")
load("descriptive_epi_mean_monthly_cohorts.Rdata")


d<-prep_desc_data(d)


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

setwd("C:/Users/andre/Documents/HBGDki/Figures/Descriptive Epi/Primary pooled")

#Wasting prevalence
p1 <- desc_epi_metaplot(d, stat="Prevalence\nof\nwasting",
                     ylabel="Wasting longitudinal prevalence",
                     title="Wasting longitudinal prevalencee")
ggsave("WastPrev_metaplot.pdf", p1, width = w, height = h, units = "in")

#Severe wasting prevalence
p2 <- desc_epi_metaplot(d, stat="Prevalence\nof\nsevere\nwasting",
                     ylabel="Severe wasting longitudinal prevalence",
                     title="Severe wasting longitudinal prevalencee")
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
pdf("Descriptive_epi_plots.pdf", width=w,height=h, paper="USr")
p1
p2
p3
p4
p5
p6
p7
p8
dev.off()


#------------------------------------
# Secondary plots (for presentation)
#------------------------------------

# 
# #get legend
# df <- d[d$statistic=="Wasting\nincidence\nrate",]
# df <- df %>% rename(Legend = stratacol)
# df$Legend <- factor(df$Legend, levels=unique(df$Legend))
# df$Legend<-recode(df$Legend, 
# "overall"= "Unstratified",          
# "pooled"= "Pooled age stratified",    
#  "pooled_unstrat"= "Pooled unstratified",              
#  "strata"= "Age stratified"            
# )
# col_legend <- c( `Unstratified`="#56B4E9", `Age stratified`="#999999" , `Pooled age stratified`="#f7a809", `Pooled unstratified`="#009E73") #, #f7a809, "#56B4E9",  "#E69F00",)
# cbPalette <- c( overall="#56B4E9", strata="#999999" , pooled="#f7a809", pooled_unstrat="#009E73") #, #f7a809, "#56B4E9",  "#E69F00",)
# 
# ggplot(df, aes(`Child age stratification`)) + 
#   geom_point(aes(x=strata, y=Mean, fill=Legend, color=Legend), size = 4) +
#   geom_linerange(aes(x=strata, ymin = Lower.95.CI, ymax = Upper.95.CI, color=Legend), 
#                  alpha=0.5, size = 3) +
#   scale_fill_manual(values=col_legend) +
#   scale_colour_manual(values=col_legend) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme(legend.position="right",
#         strip.text.x = element_text(size=8),
#         axis.text.x = element_text(size=8)) +
#   ylab("Wasting incidence rate per 1000 days")+
#   ggtitle("Wasting incidence rate")
# 
# #Plot single pooled estimate
# ggplot(d[d$statistic=="Wasting\nincidence\nrate" & d$country_cohort=="Pooled - All",], aes(`Child age stratification`)) + 
#   geom_point(aes(x=strata, y=Mean, fill=stratacol, color=stratacol), size = 4) +
#   geom_linerange(aes(x=strata, ymin = Lower.95.CI, ymax = Upper.95.CI, color=stratacol), 
#                  alpha=0.5, size = 3) +
#   scale_fill_manual(values=cbPalette) +
#   scale_colour_manual(values=cbPalette) +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   theme(strip.background = element_blank(),
#         legend.position="none",
#         strip.text.x = element_text(size=8),
#         axis.text.x = element_text(size=8)) +
#   ylab("Wasting incidence rate per 1000 days")+
#   ggtitle("Pooled wasting incidence rate")
# 
# 
# 
# 
# #Compare wasting/severe wasting IR across regions
# 
# df <- d[grep("Pool",d$country_cohort),]
# df1<-df[df$statistic=="Wasting\nincidence\nrate",]
# df2<-df[df$statistic=="Severe\nwasting\nincidence\nrate",]
# 
# df <- data.frame(df1[,1:3],ratio=df2$Mean/df1$Mean)
# 
# ggplot(df, aes(strata)) +
#   geom_point(aes(x=strata, y=ratio), size = 4) +
#   facet_wrap(~country_cohort, nrow=2) 

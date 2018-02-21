
# load packages
rm(list=ls())
library(tidyverse)
library(metafor)

source("~/Rally-6/Scripts/HBGDki_plotting_functions_6A.R")
source("~/Rally-6/Scripts/Meta-analysis functions_6A.R")


setwd("~/Rally-6/Results")
load("descriptive_epi_mean_monthly_cohorts_6A.Rdata")
d60<-d %>% mutate(washout=60)
load("descriptive_epi_mean_monthly_cohorts_30d_6A.Rdata")
d30<-d %>% mutate(washout=30)
load("descriptive_epi_mean_monthly_cohorts_90d_6A.Rdata")
d90<-d %>% mutate(washout=90)

d <- rbind(d30,d60,d90)

unique(d$statistic)


d <- d %>% filter(statistic=="Wasting\nincidence\nrate" | 
                  statistic=="Percent\nwasting\nrecovered\nin 60 days" | 
                  statistic=="Average\nduration\nof\nwasting") %>%
          filter(strata=="Overall")

d$strata<-factor(d$washout)
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

setwd("C:/Users/andre/Documents/Rally-6/Figures/")

#Wasting IR
p1 <- desc_epi_metaplot(d, stat="Wasting\nincidence\nrate",
                     ylabel="Wasting incidence rate per 1000 days",
                     title="Wasting incidence rate")

             

#ggsave("Washout_sensitivity_incidence.pdf", p1, width = w, height = h, units = "in")

#Duration
p2 <- desc_epi_metaplot(d, stat="Average\nduration\nof\nwasting",
                     ylabel="Duration of wasting (days)",
                     title="Average duration of wasting")


#ggsave("Washout_sensitivity_duration.pdf", p1, width = w, height = h, units = "in")


#Wasting IR
p3 <- desc_epi_metaplot(d, stat="Percent\nwasting\nrecovered\nin 60 days",
                     ylabel="Percent wasting recovered in 60 days",
                     title="Percent wasting recovered in 60 days")


#ggsave("Washout_sensitivity_recovery.pdf", p1, width = w, height = h, units = "in")



#Print all descriptive plots together
pdf("6A_washout_sensitivity.pdf", width=w,height=h, paper="USr")
p1
p2
p3
dev.off()


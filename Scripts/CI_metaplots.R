

rm(list=ls())

source("C:/Users/andre/Documents/Rally-6/Scripts/HBGDki_plotting_functions_6A.R")
source("C:/Users/andre/Documents/Rally-6/Scripts/Meta-analysis functions_6A.R")

load("C:/Users/andre/Documents/Rally-6/Results/CI_df_6A.Rdata")

#Drop studies that enrol acutely ill children
unique(Ndf$STUDYID)
Ndf <- Ndf %>% filter(STUDYID!="ki1000301-DIVIDS" & STUDYID!="ki1000304-LBW" & STUDYID!="ki1000306-ZincSGA" & STUDYID!="ki1000304b-ZincInf")
     


Ndf <- as.data.frame(Ndf)

res<-rbind(
data.frame(STUDYID=Ndf$STUDYID, binom.confint(Ndf$Neps06, Ndf$N06, method="exact"), agecat="0-6"),
data.frame(STUDYID=Ndf$STUDYID, binom.confint(Ndf$Neps012, Ndf$N012, method="exact"), agecat="0-12"),
data.frame(STUDYID=Ndf$STUDYID, binom.confint(Ndf$Neps018, Ndf$N018, method="exact"), agecat="0-18"),
data.frame(STUDYID=Ndf$STUDYID, binom.confint(Ndf$Neps024, Ndf$N024, method="exact"), agecat="Overall")
)




#Drop when n is 0
res <- res[res$n!=0, ]


mean_ci <- data.frame(STUDYID=res$STUDYID, COUNTRY="India", strata=res$agecat, N=res$n, statistic="revalence",
                      N=res$n, Mean=res$mean, SD=NA, `Robust SE`=NA, `Lower 95%CI`=res$lower, `Upper 95%CI`=res$upper)
#colnames(mean_ci)[6:10] <- c("Mean", "SD", "Robust SE", "Lower 95%CI", "Upper 95%CI")  

    
    
    
d<-prep_desc_data(mean_ci)
d$Mean <- d$Mean * 100
d$Lower.95.CI <- d$Lower.95.CI * 100
d$Upper.95.CI <- d$Upper.95.CI * 100





#Drop studies that enrol acutely ill children
unique(Ndf.wast$STUDYID)
Ndf.wast <- Ndf.wast %>% filter(STUDYID!="ki1000301-DIVIDS" & STUDYID!="ki1000304-LBW" & STUDYID!="ki1000306-ZincSGA" & STUDYID!="ki1000304b-ZincInf")
     


Ndf.wast <- as.data.frame(Ndf.wast)

res<-rbind(
data.frame(STUDYID=Ndf.wast$STUDYID, binom.confint(Ndf.wast$Neps06, Ndf.wast$N06, method="exact"), agecat="0-6"),
data.frame(STUDYID=Ndf.wast$STUDYID, binom.confint(Ndf.wast$Neps012, Ndf.wast$N012, method="exact"), agecat="0-12"),
data.frame(STUDYID=Ndf.wast$STUDYID, binom.confint(Ndf.wast$Neps018, Ndf.wast$N018, method="exact"), agecat="0-18"),
data.frame(STUDYID=Ndf.wast$STUDYID, binom.confint(Ndf.wast$Neps024, Ndf.wast$N024, method="exact"), agecat="Overall")
)




#Drop when n is 0
res <- res[res$n!=0, ]


mean_ci_wast <- data.frame(STUDYID=res$STUDYID, COUNTRY="India", strata=res$agecat, N=res$n, statistic="revalence",
                      N=res$n, Mean=res$mean, SD=NA, `Robust SE`=NA, `Lower 95%CI`=res$lower, `Upper 95%CI`=res$upper)
#colnames(mean_ci)[6:10] <- c("Mean", "SD", "Robust SE", "Lower 95%CI", "Upper 95%CI")  

    
dwast<-prep_desc_data(mean_ci_wast)
dwast$Mean <- dwast$Mean * 100
dwast$Lower.95.CI <- dwast$Lower.95.CI * 100
dwast$Upper.95.CI <- dwast$Upper.95.CI * 100






#Set theme and colors
theme_set(theme_bw())
d$stratacol <- "strata"
d$stratacol[d$strata=="Overall"] <- "overall"
d$stratacol[d$pooled==1] <- "pooled"
d$stratacol[d$strata=="Overall" & d$pooled==1] <- "pooled_unstrat"
dwast$stratacol <- "strata"
dwast$stratacol[dwast$strata=="Overall"] <- "overall"
dwast$stratacol[dwast$pooled==1] <- "pooled"
dwast$stratacol[dwast$strata=="Overall" & dwast$pooled==1] <- "pooled_unstrat"





cbPalette <- c( overall="#56B4E9", strata="#999999" , pooled="#f7a809", pooled_unstrat="#009E73")

#Set plot width and height
w <- 10
h <- 7

#----------------------
# Primary plots
#----------------------



#----------------------
# Primary plots
#----------------------


#Set plot directory

setwd("C:/Users/andre/Documents/Rally-6/Figures/")

#Stunting CI
p9 <- desc_epi_metaplot(d, stat="revalence",
                     ylabel="Cumulative Incidence of Stunting",
                     title="Cumulative Incidence\nof Stunting")

p9
ggsave("StuntCI_metaplot.pdf", p9, width = w, height = h, units = "in")

    


#Wasting CI
p10 <- desc_epi_metaplot(dwast, stat="revalence",
                     ylabel="Cumulative Incidence of Wasting",
                     title="Cumulative Incidence\nof Wasting")
 
p10
ggsave("WastCI_metaplot.pdf", p10, width = w, height = h, units = "in")


save(p9, p10, file="C:/Users/andre/Documents/Rally-6/Results/6A_Descriptive_epi_CI_plots.Rdata")


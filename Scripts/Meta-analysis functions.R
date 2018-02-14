

#---------------------------------
# Meta-analysis wrapper functions
#---------------------------------

descr_epi_metafun <- function(df, mthd="REML"){
  require(metafor)
  fit<-NULL
  if(grepl("ncidence", df$statistic[1])){
  fit<-rma(xi=df$N, ti=(df$N/df$Mean), data=df, method=mthd, measure="IR")
  }
  if(grepl("revalence", df$statistic[1])){
  fit<-rma(ni=df$N, xi=(df$Mean * df$N), data=df, method=mthd, measure="PR")
  }
  if(grepl("ercent", df$statistic[1])){
  tryCatch({fit<-rma(ni=df$N, xi=(df$Mean * df$N), data=df, method=mthd, measure="PR")},
           error=function(e){cat("ERROR : REML did not converge, trying ML \n")})
    if(is.null(fit)){fit<-rma(ni=df$N, xi=(df$Mean * df$N), data=df, method="ML", measure="PR")}
  }
  if(grepl("verage", df$statistic[1])){
  fit<-rma(yi=Mean, sei=(df$Upper.95.CI - df$Mean)/1.96, data=df, method="REML", measure="MN")
  }

  return(fit)
}




#---------------------------------
#Clean descriptive epi means for meta-analysis function
#---------------------------------

prep_desc_data <- function(d){
  require(metafor)

#Drop CMIN GUINEA-BISSAU (not monthly)
d<-d[d$country_cohort!="CMIN GUINEA-BISSAU",]

#mark region of each study
unique(d$COUNTRY)
d$region <- "South Asia"
d$region[d$COUNTRY=="TANZANIA" | d$COUNTRY=="SOUTH AFRICA"] <- "Sub-Saharan Africa"
d$region[d$COUNTRY=="GUATEMALA" | d$COUNTRY=="PERU" | d$COUNTRY=="BRAZIL"] <- "Latin America"

d <- d %>% select(country_cohort, region, strata, statistic, N, Mean,Lower.95.CI, Upper.95.CI)



#Calculate pooled estimates
MetaEst<-MetaEst_Africa<-MetaEst_Asia<-MetaEst_America<-NULL
for(i in 1:length(unique(d$statistic))){
  for(j in 1:length(unique(d$strata))){
  temp<-d[d$statistic==unique(d$statistic)[i] & d$strata==unique(d$strata)[j],]
  fit<-descr_epi_metafun(temp)
  MetaEst <- rbind(MetaEst, est<-data.frame(country_cohort="Pooled estimate", region="All", strata=unique(d$strata)[j], statistic=unique(d$statistic)[i], N=sum(temp$N), Mean=fit$beta,  Lower.95.CI=NA, Upper.95.CI=NA, se=fit$se))
  
  #calculate pooled estimates by continent
  fit<-descr_epi_metafun(temp[temp$region=="South Asia",])
  MetaEst_Asia <- rbind(MetaEst_Asia, est<-data.frame(country_cohort="Pooled estimate", region="South Asia", strata=unique(d$strata)[j], statistic=unique(d$statistic)[i], N=sum(temp$N), Mean=fit$beta,  Lower.95.CI=NA, Upper.95.CI=NA, se=fit$se))
  fit<-descr_epi_metafun(temp[temp$region=="Sub-Saharan Africa",])
  MetaEst_Africa <- rbind(MetaEst_Africa, est<-data.frame(country_cohort="Pooled estimate", region="Sub-Saharan Africa", strata=unique(d$strata)[j], statistic=unique(d$statistic)[i], N=sum(temp$N), Mean=fit$beta,  Lower.95.CI=NA, Upper.95.CI=NA, se=fit$se))
  fit<-descr_epi_metafun(temp[temp$region=="Latin America",])
  MetaEst_America <- rbind(MetaEst_America, est<-data.frame(country_cohort="Pooled estimate", region="Latin America", strata=unique(d$strata)[j], statistic=unique(d$statistic)[i], N=sum(temp$N), Mean=fit$beta,  Lower.95.CI=NA, Upper.95.CI=NA, se=fit$se))
  }
}

#Calculate pooled estimate 95% CI
MetaEst$Lower.95.CI<-MetaEst$Mean - 1.96 * MetaEst$se
MetaEst$Upper.95.CI<-MetaEst$Mean + 1.96 * MetaEst$se
MetaEst$country_cohort <- "Pooled - All"

MetaEst_Asia$Lower.95.CI<-MetaEst_Asia$Mean - 1.96 * MetaEst_Asia$se
MetaEst_Asia$Upper.95.CI<-MetaEst_Asia$Mean + 1.96 * MetaEst_Asia$se
MetaEst_Asia$country_cohort <- "Pooled - Asia"

MetaEst_Africa$Lower.95.CI<-MetaEst_Africa$Mean - 1.96 * MetaEst_Africa$se
MetaEst_Africa$Upper.95.CI<-MetaEst_Africa$Mean + 1.96 * MetaEst_Africa$se
MetaEst_Africa$country_cohort <- "Pooled - Africa"

MetaEst_America$Lower.95.CI<-MetaEst_America$Mean - 1.96 * MetaEst_America$se
MetaEst_America$Upper.95.CI<-MetaEst_America$Mean + 1.96 * MetaEst_America$se
MetaEst_America$country_cohort <- "Pooled - America"

Pooled <- rbind(MetaEst, MetaEst_Asia, MetaEst_Africa, MetaEst_America)
Pooled <- subset(Pooled, select = -c(se))

d$pooled <- 0
Pooled$pooled <- 1
d <- rbind(d, Pooled)
d <- rbind(d[d$region=="All",], d[d$region!="All",])
d$region <- factor(d$region, levels=unique(d$region))
d <- d %>% arrange(region, desc(pooled))
d$country_cohort <- as.character(d$country_cohort)
d$country_cohort <- factor(d$country_cohort, levels=unique(d$country_cohort))


#Rename levels for plot axes
d$country_cohort<-recode(d$country_cohort, "GMS-Nepal NEPAL"= "GMS-Nepal",          
"CMIN BANGLADESH"= "CMIN Bang.",          
"CMC-V-BCS-2002 INDIA"= "CMC India",    
 "EE PAKISTAN"= "EE Pakistan",              
 "IRC INDIA"= "IRC India",                
 "TDC INDIA"= "TDC India",                
 "ResPak PAKISTAN"= "ResPak",         
 "MAL-ED BANGLADESH"= "Mal-ED Bang.",       
  "MAL-ED INDIA"= "Mal-ED India",             
  "MAL-ED NEPAL"= "Mal-ED Nepal",                    
 "TanzaniaChild2 TANZANIA"= "Tanz. Child",  
 "MAL-ED SOUTH AFRICA"= "Mal-ED S. Afr.",      
 "MAL-ED TANZANIA"= "Mal-ED Tanz.",               
 "CONTENT PERU"= "Content Peru",             
 "Guatemala BSC GUATEMALA"= "Guat. BSC",  
 "CMIN PERU"= "CMIN Peru",                
 "Peru Huascar PERU"= "Huascar- Peru",       
 "MAL-ED BRAZIL"= "Mal-ED Brazil",            
 "MAL-ED PERU"= "Mal-ED Peru",
 "Pooled - America" = "Pooled-Amer."
)

return(d)
}

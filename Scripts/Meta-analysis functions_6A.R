

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
  
  d$country_cohort<- gsub("^k.*?-" , "", d$STUDYID)

d <- d %>% select(country_cohort, strata, statistic, N, Mean,Lower.95.CI, Upper.95.CI)



#Calculate pooled estimates
MetaEst<-NULL
for(i in 1:length(unique(d$statistic))){
  for(j in 1:length(unique(d$strata))){
  temp<-d[d$statistic==unique(d$statistic)[i] & d$strata==unique(d$strata)[j],]
  fit<-descr_epi_metafun(temp)
  MetaEst <- rbind(MetaEst, est<-data.frame(country_cohort="Pooled estimate",  strata=unique(d$strata)[j], statistic=unique(d$statistic)[i], N=sum(temp$N), Mean=fit$beta,  Lower.95.CI=NA, Upper.95.CI=NA, se=fit$se))
  }
}

#Calculate pooled estimate 95% CI
MetaEst$Lower.95.CI<-MetaEst$Mean - 1.96 * MetaEst$se
MetaEst$Upper.95.CI<-MetaEst$Mean + 1.96 * MetaEst$se
MetaEst$country_cohort <- "Pooled"
Pooled <- subset(MetaEst, select = -c(se))

d$pooled <- 0
Pooled$pooled <- 1
d <- rbind(d, Pooled)
d <- d %>% arrange(-pooled)
d$country_cohort <- as.character(d$country_cohort)
d$country_cohort <- factor(d$country_cohort, levels=unique(d$country_cohort))


#Rename levels for plot axes
d$country_cohort<-recode(d$country_cohort, 
"CMC-V-BCS-2002 INDIA"= "CMC India"    
)

return(d)
}

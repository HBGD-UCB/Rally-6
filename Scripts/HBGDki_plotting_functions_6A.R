


#----------------------------------
# HBGDki results plotting functions
# author: Andrew Mertens
#----------------------------------



#----------------------------------------------------
# Functions to create fixed or random effects pooled RRs
#----------------------------------------------------

#Fit RE and FE models
meta_fun <- function(res, method, reference=NULL){
  
  #ID reference level
  levels<-unique(res$level)
  res<-res[!is.na(res$logRR.var),]
  if(is.null(reference)){
    reference<- levels[!(levels %in% unique(res$level))]
  }
  require(metafor)
  RMAest<-data.frame(study="Pooled estimate", res$variable[1],reference ,b=NA, se=NA, RFlevel=NA)
  colnames(RMAest)<-c("study","variable","level","logRR.psi","logSE", "RFlevel")

  for(j in 1:length(unique(res$level))){
    temp<-res[res$level==unique(res$level)[j],]
  fit<-rma(yi=logRR.psi, vi=logRR.var, data=temp, method=method, measure="RR")
  est<-data.frame(study="Pooled estimate", temp$variable[1],temp$level[1] ,fit$b, fit$se, RFlevel=NA)
  colnames(est)<-c("study","variable","level","logRR.psi","logSE","RFlevel")
  RMAest<-rbind(RMAest, est)
  }
  
RMAest$RR<-exp(RMAest$logRR)
RMAest$RR.CI1<-exp(RMAest$logRR - 1.96 * RMAest$logSE)
RMAest$RR.CI2<-exp(RMAest$logRR + 1.96 * RMAest$logSE)

#rename SE column to var just to allow binding to original results:
#not used again in calculations
colnames(RMAest)[5]<-"logRR.var"
    return(RMAest)
}


#Fit RE and FE models for continious outcome
meta_fun_diff <- function(res, method, reference=NULL){
  
  #ID reference level
  levels<-unique(res$level)
  res<-res[!is.na(res$ATE),]
  if(is.null(reference)){
    reference<- levels[!(levels %in% unique(res$level))]
  }
  require(metafor)
  RMAest<-data.frame(study="Pooled estimate", res$variable[1],reference ,b=NA, se=NA, RFlevel=NA)
  colnames(RMAest)<-c("study","variable","level","ATE","SE", "RFlevel")

  for(j in 1:length(unique(res$level))){
    temp<-res[res$level==unique(res$level)[j],]
  fit<-rma(yi=ATE, vi=var, data=temp, method=method, measure="GEN")
  est<-data.frame(study="Pooled estimate", temp$variable[1],temp$level[1] ,fit$b, fit$se, RFlevel=NA)
  colnames(est)<-c("study","variable","level","ATE","SE","RFlevel")
  RMAest<-rbind(RMAest, est)
  }
  
RMAest$CI1<-(RMAest$ATE - 1.96 * RMAest$SE)
RMAest$CI2<-(RMAest$ATE + 1.96 * RMAest$SE)

#rename SE column to var just to allow binding to original results:
#not used again in calculations
colnames(RMAest)[5]<-"var"
RMAest$var<-RMAest$var^2
    return(RMAest)
}






#----------------------------------------------------
# Clean results dfs function
#----------------------------------------------------
#d=stunt_unadj[[i]]
cleandf <- function(d, RF_levels, meta_method="REML", measure="RR"){
  
  #remove grant identifier
  d$study<- gsub("^k.*?-" , "", d$STUDYID)
  
  #Create harmonized risk factor level
  if(all(RF_levels %in% c("Q1","Q2","Q3","Q4"))){
    d <- d %>% group_by(study, COUNTRY) %>% mutate(RFlevel=level, level=paste0("Q",row_number())) 
  }else{
    d$RFlevel<-d$level
  }
  
  
  #create labels for country-specific cohorts
  d$country <- as.character(d$COUNTRY)
  d$country[d$country=="TANZANIA, UNITED REPUBLIC OF"] <- "TANZANIA"
  d$study <- paste0(d$study, " ", d$country)
  
  if(measure=="RR"){
  #drop studies with <5 cases in either reference or comparison level
  d <- d %>% filter(is.na(b) | !(b<5 | d<5))
  
  d_meantab <- d %>% subset(.,select=c(study,meanLevel:meanY, mean.CI1, mean.CI2))
  colnames(d_meantab) <- c("Study","Level", "Number in quartile", "Cumulative incidence of wasting", "95% CI-lower", "95% CI-upper")

  d <- d %>% subset(.,select=c(study, variable, level, logRR.psi, logRR.var,
                                               RR, RR.CI1, RR.CI2, RFlevel)) 

      
  #Estimate and merge in pooled RR's
  d_RE<-meta_fun(d, method=meta_method)
  d<-rbind(d_RE,d) %>% subset(., select=-c(logRR.psi,  logRR.var))

  # Clean dataframes for plotting
  #Sort levels/studies for plot
  suppressWarnings(ranked.est <- d %>% group_by(study) %>% summarize(max.est = max(RR, na.rm=T)) %>% arrange(max.est))
  }else{
    #drop studies with <5 cases in either reference or comparison level
  d <- d %>% filter(is.na(refN) | refN>4)
  
  d_meantab <- d %>% subset(.,select=c(study,meanLevel:meanY, mean.CI1, mean.CI2))
  colnames(d_meantab) <- c("Study","Level", "Number in quartile", "Mean", "95% CI-lower", "95% CI-upper")

  d <- d %>% subset(.,select=c(study, variable, level, ATE, var, CI1, CI2, RFlevel)) 

      
  #Estimate and merge in pooled RR's
  d_RE<-meta_fun_diff(d, method=meta_method)
  d<-rbind(d_RE,d) 

  # Clean dataframes for plotting
  #Sort levels/studies for plot
  suppressWarnings(ranked.est <- d %>% group_by(study) %>% summarize(max.est = max(ATE, na.rm=T)) %>% arrange(max.est))
  }
  
  #Order by max estimate
  ranked.est$order <- 1:nrow(ranked.est)
  d <- left_join(d, ranked.est, by="study")
  d$order[d$study=="Pooled estimate"] <- 0
  
  #Drop studies with no estimated RR and order by size of RR and arrange the dataframe
  d$level <- factor(d$level, levels=RF_levels)
  d <- d %>% filter(!(max.est=="-Inf")) %>% arrange(order, level)
  d$x<-1:nrow(d)
  d$study <- factor(d$study , as.character(unique(d$study )))

  return(d)
}








#----------------------------------------------------
# RR Plot function
#----------------------------------------------------
forestplot_fun <- function(d, reflevel, title, units="", measure="RR", levels=NULL, free_Y=T){
  
  #change names of risk factor levels
  if(!is.null(levels)){
    d$level <- factor(as.numeric(d$level))
    levels(d$level) <- levels
  }
  
  #plot
n <- nrow(d)

plot_df <- d %>% filter(as.character(level)!=as.character(reflevel))
plot_df$size<-ifelse(plot_df$study=="Pooled estimate",1,0.5)

if(measure=="RR"){
  plot_df <- plot_df %>% filter(!is.na(RR))
}else{
  plot_df <- plot_df %>% filter(!is.na(ATE))
}

forestplot<-ggplot(data=plot_df) + 
  labs(title = d$variable, x = "Cohort", y = paste0("Risk ratio (reference = ",reflevel," ", units,")")) +
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) 

  if(length(levels)>2){
    if(free_Y){
      forestplot<-forestplot + facet_wrap(~level, ncol=1, scales = "free_y") 
    }else{
      forestplot<-forestplot + facet_wrap(~level, ncol=1) 
    }
  }



if(measure=="RR"){
  forestplot<-forestplot +
    geom_pointrange( mapping=aes(x=study, y=RR, ymin=RR.CI1, ymax=RR.CI2, colour=study, size=size)) +
    scale_y_continuous(breaks=c(0.125,0.25,0.5,1,2,4), trans='log10')  +
    geom_hline(yintercept = 1) +
    geom_vline(xintercept=1.5, color="grey20", linetype=2) 
}else{
  forestplot<-forestplot +
    geom_pointrange( mapping=aes(x=study, y=ATE, ymin=CI1, ymax=CI2, colour=study, size=size)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept=1.5, color="grey20", linetype=2) 
}


  #Facet by study and drop overall if risk factor is month
if(grepl("month", plot_df$variable[1])){
forestplot<- ggplot(data=plot_df[plot_df$study!="Pooled estimate",]) + 
  labs(title = d$variable, x = "Month", y = paste0("Risk ratio (reference = ",reflevel," ", units,")")) +
  coord_flip() +
  scale_fill_manual(values=cbPalette[-1]) +
  scale_colour_manual(values=cbPalette[-1]) +
   facet_wrap(~study) 
  if(measure=="RR"){
    forestplot<-forestplot +
      geom_pointrange( mapping=aes(x=level, y=RR, ymin=RR.CI1, ymax=RR.CI2, colour=study, size=size)) +
      scale_y_continuous(breaks=c(0.125,0.25,0.5,1,2,4), trans='log10')  +
      geom_hline(yintercept = 1) +
      geom_vline(xintercept=1.5, color="grey20", linetype=2) 
  }else{
    forestplot<-forestplot +
      geom_pointrange( mapping=aes(x=level, y=ATE, ymin=CI1, ymax=CI2, colour=study, size=size)) +
      geom_hline(yintercept = 0) +
      geom_vline(xintercept=1.5, color="grey20", linetype=2) 
}
}


#Add in formatting universal to all plot types
forestplot<-forestplot + 
  coord_flip() +
  scale_size_continuous(range = c(0.5, 1))+
  ggtitle(title) +
    theme(panel.border = element_blank(), 
    strip.background = element_blank(),
    legend.position="none")


return(forestplot)
}




#----------------------------------------------------
# Descriptive Epi plot function
#----------------------------------------------------


desc_epi_metaplot <- function(d, 
                              stat="Wasting\nincidence\nrate",
                              ylabel="Wasting incidence rate per 1000 days",
                              title="Wasting incidence rate",
                              xlabel="Child age stratification"){
if(!is.null(stat)){
  p <- ggplot(d[d$statistic==stat,]) +
              geom_point(aes(x=strata, y=Mean, fill=stratacol, color=stratacol), size = 4) +
              geom_linerange(aes(x=strata, ymin = Lower.95.CI, ymax = Upper.95.CI, color=stratacol), 
                 alpha=0.5, size = 3) 
}else{
  p <- ggplot(d) + 
              geom_point(aes(x=statistic, y=Mean, fill=stratacol, color=stratacol), size = 4) +
              geom_linerange(aes(x=statistic, ymin = Lower.95.CI, ymax = Upper.95.CI, color=stratacol), 
                 alpha=0.5, size = 3) 
}
  
p <- p + 
  facet_wrap(~country_cohort, nrow = 2)  + 
  scale_fill_manual(values=cbPalette) +
  scale_colour_manual(values=cbPalette) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.background = element_blank(),
        legend.position="none",
        strip.text.x = element_text(size=8),
        axis.text.x = element_text(size=8)) +
  ylab(ylabel)+
  ggtitle(title) + 
  xlab(xlabel)

 return(p)
}





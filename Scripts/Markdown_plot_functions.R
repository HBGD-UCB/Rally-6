

#hbgdki pallet
tableau10 <- c("#1F77B4","#FF7F0E","#2CA02C","#D62728",
  "#9467BD","#8C564B","#E377C2","#7F7F7F","#BCBD22","#17BECF")



#----------------------------------------------
#Plot functions
#----------------------------------------------

#GAM curve
WHZ_curve<-function(df){
  theme_set(theme_bw())
  
  # grab a color blind friendly palette
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cols <- cbPalette[c(1,3,7)]
  
  p <- gmsn_age<-ggplot(df, aes(x = AGEDAYS)) +
    geom_smooth(aes(y=WHZ), color=tableau10[1]) +
    geom_jitter(aes(y=WHZ,x=AGEDAYS), height = 0.2, width=0.2,  alpha = 0.1, size=0.5)+
    labs(y = "WHZ",
         x = "Child Age (Days)",
         title = "") +
    theme(strip.background = element_blank())
  
  return(p)
}


#GAM curve
HAZ_curve<-function(df){
  theme_set(theme_bw())
  
  # grab a color blind friendly palette
  cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cols <- cbPalette[c(1,3,7)]
  
  p <- gmsn_age<-ggplot(df, aes(x = AGEDAYS)) +
    geom_smooth(aes(y=WHZ), color=tableau10[2]) +
    geom_jitter(aes(y=WHZ,x=AGEDAYS), height = 0.2, width=0.2,  alpha = 0.1, size=0.5)+
    labs(y = "WHZ",
         x = "Child Age (Days)",
         title = "") +
    theme(strip.background = element_blank())
  
  return(p)
}


#spaghetti plot
spaghetti<-function(df, Zscore="WHZ"){
  theme_set(theme_bw())
  
    df <- df %>% group_by(SUBJID) %>%
    mutate(rand=runif(n = 1)) %
    
    ord<- df %>% slice(1) %>% ungroup() %>% arrange(desc(rand))
  cutoff_10<-ord$rand[10]
  
  
  set.seed(12345)
  d <- df %>% group_by(SUBJID) %>%
    mutate(alpha=  ifelse(rand >= cutoff_10 , 1, 0.1))
  
  if(Zscore=="WHZ"){
  p <- ggplot(d, aes(x=AGEDAYS, y=WHZ)) + 
    geom_line() + guides(colour=FALSE) + xlab("Child Age (Days)") +
    ylab("WHZ") + aes(alpha=alpha, group=factor(SUBJID)) + guides(alpha=FALSE)
  }
  
  if(Zscore=="HAZ"){
  p <- ggplot(d, aes(x=AGEDAYS, y=WHZ)) + 
    geom_line() + guides(colour=FALSE) + xlab("Child Age (Days)") +
    ylab("HAZ") + aes(alpha=alpha, group=factor(SUBJID)) + guides(alpha=FALSE)
  }
  return(p)
}


#Heatmap
heatmap<-function(d){
  theme_set(theme_bw())
  
  #generate whz categories
  d$wastcat <- 0
  d$wastcat[d$WHZ<(-2)] <- -1
  d$wastcat[d$WHZ<(-3)] <- -2
  
  #Sum an ad-hoc total wasting score to rank children by  
  d <- d %>% group_by(SUBJID) %>% mutate(wastscore=sum(wastcat, na.rm=T)) %>% ungroup()
  table(d$wastscore)
  
  #Add level for missingness
  d$wastcat[is.na(d$WHZ)] <- 1
  d$wastcat <- as.factor(d$wastcat)
  d$SUBJID <- as.factor(d$SUBJID)
  table(d$wastcat)
  
  #make ordered childnum by amount of wasting
  d <- d %>%
    arrange(wastscore) %>%  
    mutate(SUBJID = factor(SUBJID, unique(SUBJID))) 
  d$childnum <- as.numeric(d$SUBJID)
  head(as.data.frame(d))
  
  #Create a child age in months variable
  d$agemonths<-floor(d$AGEDAYS/30)
  
  
  if(sum(d$wastcat== -2) > 0){
    levels(d$wastcat)<-c("Severely wasted", "Wasted", "Not wasted")
    cbPalette <- c("#56B4E9", "#E69F00", "#c1c1c1")
  }else{
    levels(d$wastcat)<-c( "Wasted", "Not wasted")
    cbPalette <- c( "#E69F00", "#c1c1c1")   
  }
  
  p <- ggplot(d,aes(x=agemonths,y=childnum)) + 
    geom_tile(aes(fill=wastcat)) + scale_fill_manual(values=cbPalette) +
    xlab("Child Age (Months)") + ylab("Child number")
  ggtitle("Heatmap of wasting and severe wasting episodes") +
    theme(legend.title=element_blank())
  
  return(p)
}





#Incidence/prevalence figures
means_plot<-function(df){
  
  #df<-df[df$strata!=">24 months",]
  df$statistic <- factor(df$statistic, levels = unique(df$statistic))
  
  theme_set(theme_bw())
  cbPalette <- c( "#56B4E9" , "#999999" , "#999999", "#999999" , "#999999", "#999999" ) #, #f7a809, "#56B4E9",  "#E69F00",)
  p <- ggplot(df, aes(`Child age stratification`)) + 
    geom_point(aes(x=strata, y=Mean, fill=strata, color=strata), size = 4) +
    geom_linerange(aes(x=strata, ymin = Lower.95.CI, ymax = Upper.95.CI, color=strata), 
                   alpha=0.5, size = 3) +
    facet_wrap(~statistic, scales = "free")  + 
    scale_fill_manual(values=cbPalette) +
    scale_colour_manual(values=cbPalette) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    theme(strip.background = element_blank(),
          legend.position="none")
  
  return(p)
}







# 

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



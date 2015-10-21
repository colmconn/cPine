rm(list=ls())
graphics.off()

##library(gmodels)
library(gdata)
library(ggplot2)
library(reshape)

## http://wiki.stdout.org/rcookbook/Graphs/Plotting%20means%20and%20error%20bars%20(ggplot2)/
## Summarizes data.
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
    require(plyr)

    ## New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }
    
    ## This is does the summary; it's not easy to understand...
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun= function(xx, col, na.rm) {
                     c(N    = length2(xx[,col], na.rm=na.rm),
                       mean = mean (xx[,col], na.rm=na.rm),
                       median = median (xx[,col], na.rm=na.rm),
                       mad = mad (xx[,col], na.rm=na.rm),                       
                       sd   = sd   (xx[,col], na.rm=na.rm),
                       min  = min  (xx[,col], na.rm=na.rm),
                       max  = max  (xx[,col], na.rm=na.rm),                       
                       nacount  = sum  (is.na((xx[,col])))
                       )
                   },
                   measurevar,
                   na.rm
                   )
   
    ## Rename the "mean" column    
    datac <- rename(datac, c("mean"=measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  ## Calculate standard error of the mean
    
    ## Confidence interval multiplier for standard error
    ## Calculate t-statistic for confidence interval: 
    ## e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
  }

if ( Sys.info()["sysname"] == "Darwin" ) {
   ## Group.data.dir="/Volumes/opt/mriAnalyses/MDD/data/Group.data"
  ## Group.results.dir="/Volumes/opt/mriAnalyses/MDD/data/Group.results"  
  ## scriptsDir="/Volumes/opt/mriAnalyses/MDD/script"

  Admin.data.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/admin"
  Group.data.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.data"    
  Group.results.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.results"  
  scriptsDir="/Volumes/PROMISEPEGASUS/yangdata/cPine/script"
} else if ( Sys.info()["sysname"] == "Linux" ) {
  Sys.setenv(AFNI_R_DIR="/abin")
  ## Group.data.dir="/mnt/nfs/yangdata/restingstate/data/Group.data"
  ## Group.results.dir="/mnt/nfs/yangdata/restingstate/data/Group.results"
  ## scriptsDir="/mnt/nfs/yangdata/restingstate/scripts"
}

## this file stores the order of the subjects in each of the following BRIK files
#subjectOrderFilename=file.path(Group.data.dir, paste("subjectOrder", groups, seed, "REML.csv", sep="."))
## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
##demographicsFilename=file.path(Group.data.dir, "demographics2.csv")

demographicsFilename=file.path(Admin.data.dir, "data_entry31Jan2013.csv")

## the file containins the hits, misses, false alarms etc
summaryTabFilename=file.path(Group.data.dir, "summary.tab")

cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=TRUE, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!"))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))
demographics$ID=as.factor(gsub("_[ABCD]", "", as.character(demographics$ID), fixed=FALSE))

cat("*** Reading", summaryTabFilename, "\n")
summaryTab=read.table(summaryTabFilename, header=TRUE, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!"))
cat(sprintf("*** Read Pine task summary data for %s unique subjects\n",  length(unique(summaryTab$Subject))))

summaryTab$Subject=gsub("_[ABCD]", "", as.character(summaryTab$Subject))
summaryTab$Subject[summaryTab$Subject=="300"]="169/300"
summaryTab$Subject=as.factor(summaryTab$Subject)

summaryTab$Group=demographics[match(summaryTab$Subject, demographics$ID), "Grp"]
summaryTab$Group=drop.levels(summaryTab$Group)
 
summaryTab.melted=melt(summaryTab, id.vars=c("Subject", "Group"),
  measure.vars=c("Hits", "Misses", "CorrectRejections", "FalseAlarms", "Omissions"), variable_name="measure")    

## graph=graph=ggplot(roistats.summary, aes(x=Group, y=value, fill=Group)) +
##   theme_bw(base_size =  my.base.size) +
##   geom_bar(position=position_dodge()) +
##   geom_errorbar(aes(ymin=value-se, ymax=value+se), width=.2, position=position_dodge(.9)) +
##   facet_wrap(~cluster, scales="free") +
##   scale_fill_brewer(palette="Set1") +
##   labs(y=ylabel) +
##   opts(title = graphTitle,
##        legend.position="none",
##        ##panel.grid.major = theme_blank(),
##        panel.grid.minor = theme_blank(),
##        axis.title.x=theme_blank(),
##        axis.title.x = theme_text(size=my.base.size, vjust=0),
##        axis.title.y = theme_text(size=my.base.size, vjust=0.4, angle =  90),
##        plot.title=theme_text(size=my.base.size*1.2, vjust=1))


graph=ggplot(summaryTab.melted, aes(x=value, fill=measure) ) +
  geom_bar(position="fill") +
  facet_wrap(~Group) 

print(graph)

summarySE(summaryTab.melted, measure="value", groupvars=c("Group", "measure"))



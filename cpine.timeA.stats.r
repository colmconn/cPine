rm(list=ls())
graphics.off()

##library(gmodels)
library(gdata)
library(ggplot2)
library(compute.es)

stack <- function(){ 
   it <- list() 
   res <- list( 
               push=function(x){ 
                 it[[length(it)+1]] <<- x 
               }, 
               pop=function(){ 
                 val <- it[[length(it)]] 
                 it <<- it[-length(it)] 
                 return(val) 
               }, 
               value=function(){ 
                 return(it) 
               } 
               ) 
   class(res) <- "stack" 
   res 
} 
print.stack <- function(x,...){ 
   print(x$value()) 
} 
push <- function(stack,obj){ 
   stack$push(obj) 
} 
pop <- function(stack){ 
   stack$pop() 
}

make.significance.indications <- function(pValues, which.pValues=c(1)) {

  Signif=symnum(pValues, corr=FALSE, na=FALSE, cutpoints = c(0,  .001,.01, .05, .1, 1),
    symbols   =  c("***", "**", "*", ".", " "))
  f=format(Signif)

  ## only return the first one as we're only interested in marking significant group effects
  return(f[which.pValues])
}

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

## Reads the seed file and does the (crude) equivalent of BAS variable
## substitution
readSeedsFile <- function (inSeedsFile) {
  cat("*** Reading seed from", inSeedsFile, "\n")
  table=scan(inSeedsFile, what=character())
  table=gsub("$scriptsDir", scriptsDir, table, fixed=TRUE)

  return (table)
}

## extracts the seed name from a file path name pointing to a NIfTI
## file containing the seed
getSeedName <- function(inSeedPath){
  name=basename(inSeedPath)
  return(gsub("\\.nii.*", "", name))
}


makeTableString <- function(inGroup, inMean, inSe, inMin, inMax, inNaCount, inMissingData=TRUE) {
  ##  st=paste(round(inMean, 1), " / ", round(inSe, 1),    

  st=paste(round(inMean, 1), " ± ", round(inSe, 1),
    " (", round(inMin, 1), "-", round(inMax, 1), ")", ifelse(inMissingData & inNaCount > 0, paste(" [", inNaCount, "]", sep=""), ""), sep="")
  return(st)
}

groups="mddAndCtrl"
contrast="allRemVsAllNotrem"

## this file stores the order of the subjects in each of the following BRIK files
subjectOrderFilename=file.path(Group.data.dir, paste("subjectOrder", groups, contrast, "REML.csv", sep="."))

## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
##demographicsFilename=file.path(Group.data.dir, "demographics2.csv")

demographicsFilename=file.path(Admin.data.dir, "data_entry31Jan2013.csv")

cat("*** Reading", subjectOrderFilename, "\n")
subjectOrder=read.csv(subjectOrderFilename, header=T)
subjectOrder$subject=as.factor(gsub("_A", "", as.character(subjectOrder$subject), fixed=TRUE))

colnames(subjectOrder)=c("subject")

subjectOrder$subject=as.factor(gsub("_A", "", as.character(subjectOrder$subject), fixed=TRUE))

##cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(subjectOrder$subject))))

cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!"))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))
demographics$ID=as.factor(gsub("_[ABCD]", "", as.character(demographics$ID), fixed=FALSE))


##nodata=data.frame(subject=c(107, 108, 124, 146, 148, 153, 303, 308, 341, 347, 357, 117, 161, "169/300", 313, 316, 340, 345, 344))

selectedColumns=c(
  "Grp", "Gender", "DOB", "MRI", "Race", "Hand", "SES", "Tanner1", "Tanner2", "TannerAvg",
  "CGI.CGAS",                 "CDRS.raw",                  "CDRS.tscore",               "SFMQ.total",              
  "MADRS.raw",                "WASI.PERF",                 "WASI.Full.4",               "PSWQ", 			"CoRum",
  "RADS.DM",                  "RADS.AN",                   "RADS.NS",                   "RADS.SC",
  "RADS.DM.Tscore",           "RADS.AN.Tscore",            "RADS.NS.Tscore",            "RADS.SC.Tscore",  "RADS.Total.Tscore",
  "BDI.II",                   "CDI",                       "MASC.total",                "MASC.tscore",           "RSQ.fixed",
  "SCARED.som.panic",         "SCARED.gen.anx",            "SCARED.seper.anx",         
  "SCARED.soc.phobia",        "SCARED.school.phobia",      "BIS",                       "BAS.drive",                
  "BAS.funseek",              "BAS.reward.responsiveness", "PsycAche.total")

if (length(setdiff(selectedColumns, colnames(demographics))) != 0) {
  cat ("*** The following column name is causing problmes", setdiff(selectedColumns, colnames(demographics)), "\n")
  stop("There is a mismatch between the selected columns and the columns in the demographics data frame. Stopping.")
}

mgd=cbind(subjectOrder, demographics[match(subjectOrder$subject, demographics$ID), selectedColumns])
  
##mgd.dti=cbind(dti.subjectOrder, demographics[match(dti.subjectOrder$subject, demographics$ID), selectedColumns])

##nodata.mgd=cbind(nodata, demographics[match(nodata$subject, demographics$ID), c("Grp", "Gender")])

nonParametricRegexp="Tanner|SES|PSWQ|CGI.CGAS|CoRum|RSQ|Hand"

##stop()

psychMeasures=list(
  list(variable="age.in.years", name="Age at time of scan (years)"),
  list(variable="Hand",         name="Edinburgh Handedness Inventory"),
  list(variable="SES",          name="Hollingshead Socioeconomic Score"),
  list(variable="TannerAvg",    name="Tanner Score"),
  ##list(variable="WASI.PERF",    name="Wechsler Abbreviated Scale of Intelligence (Performance)"),
  list(variable="WASI.Full.4",  name="Wechsler Abbreviated Scale of Intelligence"),
  
  list(variable="CGI.CGAS",     name="Children's Global Assessment Scale"),
  list(variable="PSWQ",         name="Penn State Worry Questionnaire"),
  list(variable="CoRum",        name="Corumination Questionnaire"),
  list(variable="RSQ.fixed",    name="Ruminative Responses Styles Questionnaire"),
  
  list(variable="CDRS.tscore",  name="Children's Depression Rating Scale (Standardized)"),

  ## list(variable="RADS.DM",      name="Reynolds Adolescent Depression Scale Dysphoric Mood"),
  ## list(variable="RADS.AN",      name="Reynolds Adolescent Depression Scale Anhedonia/Negative Affect"),
  ## list(variable="RADS.NS",      name="Reynolds Adolescent Depression Scale Negative Self-evaluation"),
  ## list(variable="RADS.SC",      name="Reynolds Adolescent Depression Scale Somatic Complaints"),    

  list(variable="RADS.DM.Tscore",      name="Reynolds Adolescent Depression Scale Dysphoric Mood (Standardized)"),
  list(variable="RADS.AN.Tscore",      name="Reynolds Adolescent Depression Scale Anhedonia/Negative Affect (Standardized)"),
  list(variable="RADS.NS.Tscore",      name="Reynolds Adolescent Depression Scale Negative Self-evaluation (Standardized)"),
  list(variable="RADS.SC.Tscore",      name="Reynolds Adolescent Depression Scale Somatic Complaints (Standardized)"),
  list(variable="RADS.Total.Tscore",      name="Reynolds Adolescent Depression Scale Total (Standardized)"),  

  
  ## list(variable="SFMQ.total", name="SFMQ"),
  list(variable="MADRS.raw",    name="Montgomery-Asberg Depression Rating Scale"),
  list(variable="BDI.II",       name="Beck Depression Inventory II"),
  list(variable="CDI",          name="Children's Depression Inventory"),
  list(variable="MASC.total", name="Multidimensional Anxiety Scale for Children"),
  list(variable="MASC.tscore", name="Multidimensional Anxiety Scale for Children (Standardized)")  
  ## list(variable="SCARED.som.panic", name="SCARED Som. Panic"),
  ## list(variable="SCARED.gen.anx", name="SCARED Ganeralized Anxiety"),
  ## list(variable="SCARED.seper.anx", name="SCARED Seperation Anxiety"),
  ## list(variable="SCARED.soc.phobia", name="SCARED Social Phobia"),
  ## list(variable="BIS", name="Barrett's Impulsivity Scale"),
  ## list(variable="BAS.drive", name="BAS.drive"),
  ## list(variable="BAS.funseek", name="BAS.funseek"),
  ## list(variable="BAS.reward.responsiveness", name="BAS Reward Responsiveness"),
  ## list(variable="PsycAche.total", name="PsycAche.total")
  )

analyse <- function(inData, inGenderGraphName=NULL, inRaceGraphName=NULL, inRunGenderAnalysis=FALSE) {
  ## this complicated looking regexp stuff cleans up the years in DOB
  ## and MRI with 4 digits to be just 2 digits
  inData$DOB=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", inData$DOB)
  inData$MRI=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", inData$MRI)
  
  ## inData$DOB=as.Date(inData$DOB, "%m/%d/%y")
  ## inData$MRI=as.Date(inData$MRI, "%m/%d/%y")
  ## inData$age.in.days=difftime(inData$MRI, inData$DOB, units="days")
  ## inData$age.in.days=as.numeric(inData$age.in.days)
  
  ## inData$age.in.weeks=difftime(inData$MRI, inData$DOB, units="weeks")
  ## inData$age.in.weeks=as.numeric(inData$age.in.weeks)
  ## inData$age.in.years=(inData$age.in.weeks)/52
  
  inData$Grp=drop.levels(inData$Grp)
  inData$Gender=drop.levels(inData$Gender)
  ## inData$Race=drop.levels(inData$Race)
  inData$Race=as.factor(tolower(drop.levels(inData$Race)))
  
  cat("\n*** Gender\n")
  gender.table=table(inData[, c("Gender", "Grp")])
  ##ethnicity.table=table(inData$ethnicity, inData$Group)
  print(addmargins(gender.table))
  print(prop.test(gender.table))
  
  cat("\n*** Race\n")
  race.table=table(inData[, c("Race", "Grp")])
  ##ethnicity.table=table(inData$ethnicity, inData$Group)
  print(addmargins(race.table))
  print(prop.test(race.table))
  
  cat("*** Now performing tests on psychMeasures\n")
  mystack <- stack()
  header="Characteristic,Control,MDD,DF,Stat.,pValue,EffectSize,Signif."
  for (i in 1: length(psychMeasures)) {
    variable=psychMeasures[[i]]$variable
    name=psychMeasures[[i]]$name
    
    cat("################################################################################\n");
    cat("Summary for ", name, "variable = ", variable, "\n")

    if (any(inData[, variable] < 0, na.rm=TRUE )) {
      cat ("****************************************************************************************************\n")
      cat (sprintf("*** The following subjects have data less than zero for %s (%s)\n", name, variable))
      print(as.vector(inData[inData[, variable] < 0, "subject"]))
      print(as.vector(inData[inData[, variable] < 0, "Grp"]))
      cat ("*** Now setting these values to NA\n")
      inData[inData[, variable] < 0 & !is.na(inData[, variable]), variable]=NA
      inData[inData[, variable] < 0 & !is.na(inData[, variable]), variable]=NA      
      cat ("****************************************************************************************************\n")      
    }
    
    sm.df=summarySE(inData, measure=variable, groupvars=c("Grp"), na.rm=TRUE)
    print(sm.df)
    
    ##print(sm.df)
    ctrl.string=""
    mdd.string=""
    test=""
    if (any(grep(nonParametricRegexp, variable))) {
      ctrl.string=makeTableString(sm.df[2, 1], inMean=sm.df[2, "median"],  sm.df[2, "mad"], sm.df[2, "min"], sm.df[2, "max"], sm.df[2, "nacount"], inMissingData=TRUE)
      mdd.string=makeTableString(sm.df[1, 1], inMean=sm.df[1, "median"],  sm.df[1, "mad"], sm.df[1, "min"], sm.df[1, "max"], sm.df[1, "nacount"], inMissingData=TRUE)
      
      name=paste(name, "†", sep="")
    } else {
      ctrl.string=makeTableString(sm.df[2, 1], sm.df[2, variable],  sm.df[2, "se"], sm.df[2, "min"], sm.df[2, "max"], sm.df[2, "nacount"], inMissingData=TRUE)
      mdd.string=makeTableString(sm.df[1, 1], sm.df[1, variable],  sm.df[1, "se"], sm.df[1, "min"], sm.df[1, "max"], sm.df[1, "nacount"], inMissingData=TRUE)
    }
    if (any(is.na(inData[, variable]))) {
      cat (sprintf("*** The following subjects have missing data for %s (%s)\n", name, variable))
      print(as.vector(inData[is.na(inData[, variable]), "subject"]))
      print(as.vector(inData[is.na(inData[, variable]), "Grp"]))      
    }

    ##cat("*** Analyzing ", variable, "\n")
    if (sm.df[1, "N"] > 3 && sm.df[2, "N"] > 3) {
      if (any(grep(nonParametricRegexp, variable))) {
        ##cat("Control length: ", length(inData[inData$Grp=="NCL", variable]), "MDD length: ", length(inData[inData$Grp=="MDD", variable]), "\n")
        if( inherits(test <- try(wilcox.test(inData[inData$Grp=="NCL", variable],
                                             inData[inData$Grp=="MDD", variable]),
                                 silent=FALSE),
                     "try-error") ) {
          test <- 0
        }
      } else {
        if( inherits(test <- try(t.test(inData[inData$Grp=="NCL", variable],
                                        inData[inData$Grp=="MDD", variable]),
                                 silent=FALSE),
                     "try-error") ) {
          test <- 0
        }
      }
      print(test)
      if (is.list (test) ) {
        if  ( ! any(grep(nonParametricRegexp, variable))) {
          ##cat("a=" , length(inData[inData$Grp=="NCL", variable]), "\n")
          ##cat("b=" , length(inData[inData$Grp=="MDD", variable]), "\n")
          my.es=tes(test$statistic, length(inData[inData$Grp=="NCL", variable]), length(inData[inData$Grp=="MDD", variable]))
          ## g.t will give you Hedge's g change it to d.t to get cohen's d
          effectSize=round(my.es$MeanDifference["g.t"], 5)
        } else {
          ## This is the propabilty of superiority see Erceg-Hurn
          ## DM, Mirosevich VM. Modern robust statistical methods: an
          ## easy way to maximize the accuracy and power of your
          ## research. Am Psychol. 2008;63(7):591–601.
          effectSize=round(test$statistic/(length(inData[inData$Grp=="NCL", variable]) * length(inData[inData$Grp=="MDD", variable])), 5)
        }
      }
    } else {
      cat ("*** Insufficient opservations\n")
    } ## end of if (sm.df[1, "N"] > 3 && sm.df[2, "N"] > 3) {
    st=paste(name, ctrl.string, mdd.string,
      ifelse(is.list(test), ifelse(is.null(test$parameter), "NA", round(test$parameter, 2)), ""),
      ifelse(is.list(test), round(test$statistic, 2), ""),
      ifelse(is.list(test), round(test$p.value, 5), ""),
      effectSize,
      ifelse(is.list(test), make.significance.indications(test$p.value), ""), sep=",")
    push(mystack, st)
  }

  l=mystack$value()
  cat(header, "\n")
  for (i in 1:length(l)) {
    cat (l[[i]], "\n")
  }
  ##stop()

##   if ( inRunGenderAnalysis ) {
## #########################################################################################################################################################################

##     for ( group in c("NCL", "MDD" )) {
##       cat ("\n\n####################################################################################################\n")
##       cat ("####################################################################################################\n")
##       cat("### Now performing tests on Males vs Females ", ifelse(group=="NCL", "Normal Controls", "MDD"), " psychMeasures\n")
##       cat ("####################################################################################################\n")
##       cat ("####################################################################################################\n")

##       ## Just pick out the normal controls
##       aData=inData[inData$Grp==group, ]
##       aData$Race=drop.levels(aData$Race)
      
##       cat("\n*** Gender\n")
##       gender.table=table(aData[, "Gender"])
##       ##ethnicity.table=table(aData$ethnicity, aData$Gender)
##       print(addmargins(gender.table))
##       print(prop.test(gender.table))
      
##       cat("\n*** Race\n")
##       race.table=table(aData[, c("Race", "Gender")])
##       ##ethnicity.table=table(aData$ethnicity, aData$Gender)
##       print(addmargins(race.table))
##       print(prop.test(race.table))
      
##       mystack <- stack()
##       header="Characteristic,Male,Female,DF,Stat.,Signif."
##       for (i in 1: length(psychMeasures)) {
##         variable=psychMeasures[[i]]$variable
##         name=psychMeasures[[i]]$name
        
##         ##cat("################################################################################\n");
##         ##cat("Summary for ", name, "variable = ", variable, "\n")
##         sm.df=summarySE(aData, measure=variable, groupvars=c("Gender"), na.rm=TRUE)
##         ##print(sm.df)
##         male.string=""
##         female.string=""
##         test=""
##         if (any(grep(nonParametricRegexp, variable))) {
##           male.string=makeTableString(sm.df[2, 1], inMean=sm.df[2, "median"],  sm.df[2, "mad"], sm.df[2, "min"], sm.df[2, "max"], sm.df[2, "nacount"], inMissingData=FALSE)
##           female.string=makeTableString(sm.df[1, 1], inMean=sm.df[1, "median"],  sm.df[1, "mad"], sm.df[1, "min"], sm.df[1, "max"], sm.df[1, "nacount"], inMissingData=FALSE)
          
##           name=paste(name, "†", sep="")
##         } else {
##           male.string=makeTableString(sm.df[2, 1], sm.df[2, variable],  sm.df[2, "se"], sm.df[2, "min"], sm.df[2, "max"], sm.df[2, "nacount"], inMissingData=FALSE)
##           female.string=makeTableString(sm.df[1, 1], sm.df[1, variable],  sm.df[1, "se"], sm.df[1, "min"], sm.df[1, "max"], sm.df[1, "nacount"], inMissingData=FALSE)
##         }
##         ##cat("*** Analyzing ", variable, "\n")
##         if (sm.df[1, "N"] > 3 && sm.df[2, "N"] > 3) {
##           if (any(grep(nonParametricRegexp, variable))) {
##             ##cat("Male length: ", length(aData[aData$Gender=="M", variable]), "Female length: ", length(aData[aData$Gender=="F", variable]), "\n")
##             if( inherits(test <- try(wilcox.test(aData[aData$Gender=="M", variable],
##                                                  aData[aData$Gender=="F", variable]),
##                                      silent=FALSE),
##                          "try-error") ) {
##               test <- 0
##             }
##           } else {
##             if( inherits(test <- try(t.test(aData[aData$Gender=="M", variable],
##                                             aData[aData$Gender=="F", variable]),
##                                      silent=FALSE),
##                          "try-error") ) {
##               test <- 0
##             }
##           }
##           print(test)
##         } else {
##           cat ("*** Insufficient opservations\n")
##         }
##         st=paste(name, male.string, female.string,
##           ifelse(is.list(test), ifelse(is.null(test$parameter), "NA", round(test$parameter, 2)), ""),
##           ifelse(is.list(test), round(test$statistic, 2), ""),
##           ifelse(is.list(test), round(test$p.value, 5), ""),
##           ifelse(is.list(test), make.significance.indications(test$p.value), ""), sep=",")
##         push(mystack, st)
##       }
      
##       l=mystack$value()
##       cat(header, "\n")
##       for (i in 1:length(l)) {
##         cat (l[[i]], "\n")
##       }
##     } ## end of for group in c("NCL", "MDD") {
##   } ## end of  if ( 1 == 0 ) {
  
  ## range=seq(min(floor(inData$age.in.years)), max(ceiling(inData$age.in.years)), by=1)
  ## my.base.size=18  
  ## graph=
  ##   ggplot(inData, aes(x=age.in.years)) +
  ##     geom_histogram(aes(y = ..count.., fill=Gender),width=0.4, position = position_dodge(width=0.5), breaks=range, drop=FALSE) +
  ##       ##theme_bw(base_size =  my.base.size) +
  ##       scale_fill_brewer(palette="Set1") +
  ##         labs(x="Age (years)") +
  ##           opts(title="Distribution of genders in controls",
  ##                ##panel.grid.major = theme_blank(),
  ##                panel.grid.minor = theme_blank(),
  ##                ##axis.title.x=theme_blank()##,
  ##                axis.title.x = theme_text(size=my.base.size, vjust=0),
  ##                axis.title.y = theme_text(size=my.base.size, vjust=0.4, angle =  90),
  ##                plot.title=theme_text(size=my.base.size*1.2, vjust=1)
  ##                )
  ## print(graph)
  ## ggsave(inGenderGraphName, graph)


  ## graph=
  ##   ggplot(inData, aes(x=Race)) +
  ##     geom_histogram(aes(y = ..count.., fill=Gender),width=0.4, position = position_dodge(width=0.5), drop=FALSE) +
  ##       ##theme_bw(base_size =  my.base.size) +
  ##       scale_fill_brewer(palette="Set1") +
  ##         labs(x="Race") +
  ##           opts(title="Distribution of Race in controls",
  ##                ##panel.grid.major = theme_blank(),
  ##                panel.grid.minor = theme_blank(),
  ##                axis.title.x=theme_blank(),
  ##                axis.title.x = theme_text(size=my.base.size, vjust=0),
  ##                axis.title.y = theme_text(size=my.base.size, vjust=0.4, angle =  90),
  ##                axis.text.x  = theme_text(angle=90, size=my.base.size*0.5),
  ##                plot.title=theme_text(size=my.base.size*1.2, vjust=1)
  ##                )
  ## print(graph)
  ## ggsave(inRaceGraphName, graph)
}

if ( 1 == 1 ) { 

  cat("####################################################################################################\n")
  cat("### Resting State Statistics\n")
  
  mgd$DOB=as.Date(mgd$DOB, "%m/%d/%Y")
  mgd$MRI=as.Date(mgd$MRI, "%m/%d/%Y")
  mgd$age.in.days=difftime(mgd$MRI, mgd$DOB, units="days")
  mgd$age.in.days=as.numeric(mgd$age.in.days)
  
  mgd$age.in.weeks=difftime(mgd$MRI, mgd$DOB, units="weeks")
  mgd$age.in.weeks=as.numeric(mgd$age.in.weeks)
  mgd$age.in.years=(mgd$age.in.weeks)/52
  
  ##sink("../data/Group.results/restingStateStatictics28Nov.txt")
  ## analyse(mgd, "../data/Group.results/restingStateDistributionOfGenders.png", "../data/Group.results/restingStateDistributionOfRaces.png")
  
  analyse(mgd)
  
  
  ##sink()
  
  graphics.off()
  ## cat("\n*** Age\n")
  ## print(summary(aov(age.in.days ~ Grp * Gender, data=mgd)))
  ## print(t.test(subset(mgd, Grp=="NCL", "age.in.days")[[1]], subset(mgd, Grp=="MDD", "age.in.days")[[1]]))
  
  ## cat("\n*** CDRS.CGI\n")
  ## print(summary(aov(CDRS.CGI ~ Grp * Gender, data=mgd)))
  ## print(t.test(subset(mgd, Grp=="NCL", "CDRS.CGI")[[1]], subset(mgd, Grp=="MDD", "CDRS.CGI")[[1]]))
  
  ## cat("\n*** CDRS.raw\n")
  ## print(summary(aov(CDRS.raw ~ Grp * Gender, data=mgd)))
  ## print(t.test(subset(mgd, Grp=="NCL", "CDRS.raw")[[1]], subset(mgd, Grp=="MDD", "CDRS.raw")[[1]]))
  
  ## cat("\n*** MADRS.raw\n")
  ## print(summary(aov(CDRS.raw ~ Grp * Gender, data=mgd)))
  ## print(t.test(subset(mgd, Grp=="NCL", "MADRS.raw")[[1]], subset(mgd, Grp=="MDD", "MADRS.raw")[[1]]))
  
  ## cat("\n*** BDI.II\n")
  ## print(summary(aov(CDRS.raw ~ Grp * Gender, data=mgd)))
  ## print(t.test(subset(mgd, Grp=="NCL", "BDI.II")[[1]], subset(mgd, Grp=="MDD", "BDI.II")[[1]]))
  
}

rm(list=ls())
graphics.off()

##library(gmodels)
library(gdata)
library(ggplot2)
library(lubridate)
library(compute.es)
library(reshape)
library(orddom)

source("scoreMasc.r")

########################################################################################################################################################################################################
### START OF FUNCTIONS #################################################################################################################################################################################
########################################################################################################################################################################################################

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



analyse <- function(inData, inGenderGraphName=NULL, inRaceGraphName=NULL, inRunGenderAnalysis=FALSE) {

    mystack <- stack()
    header="Characteristic,Control,MDD,Stat.,pValue,Effect Size (95% CI),Signif."
    
    inData$Grp=drop.levels(inData$Grp)
    inData$Gender=drop.levels(inData$Gender)
    ## inData$Race=drop.levels(inData$Race)
    inData$Race=as.factor(tolower(drop.levels(inData$Race)))
    
    cat("\n*** Gender\n")
    gender.table=table(inData[, c("Gender", "Grp")])
    gender.test=prop.test(gender.table)
    gender.table=addmargins(gender.table)
    print(gender.table)
    print(gender.test)

    ##ethnicity.table=table(inData$ethnicity, inData$Group)
    ## print(gender.table)
    
    csvLine=paste("Number of participants (n)", gender.table["Sum", "NCL"], gender.table["Sum", "MDD"], ",,,,,", sep=',')
    push(mystack, csvLine)
    csvLine=paste("Gender (M / F)", paste(gender.table["M", "NCL"], gender.table["F", "NCL"], sep="/"),
        paste(gender.table["M", "MDD"], gender.table["F", "MDD"], sep="/"),
        sprintf("Chi(%0.2f) = %0.2f",
                round(gender.test$parameter, 2), round(gender.test$statistic, 2)),
        round(gender.test$p.value, 2), "", make.significance.indications(gender.test$p.value), sep=",")
    push(mystack, csvLine)
    
    ## cat("\n*** Race\n")
    ## race.table=table(inData[, c("Race", "Grp")])
    ## ##ethnicity.table=table(inData$ethnicity, inData$Group)
    ## print(addmargins(race.table))
    ## print(prop.test(race.table))
    
    cat("*** Now performing tests on psychMeasures\n")
    
    for (i in 1: length(psychMeasures)) {
        variable=psychMeasures[[i]]$variable
        name=psychMeasures[[i]]$name
        
        ##cat("################################################################################\n");
        ##cat("Summary for ", name, "variable = ", variable, "\n")
        if (is.factor(inData[, variable])) {
            message (sprintf("%s, is a factor. Skipping\n", variable))
        } else {
            
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
            ##print(sm.df)
            
            ##print(sm.df)
            ctrl.string=""
            mdd.string=""
            test=""
            if (doNonParametricTestsOnVariable(variable)) {
                ctrl.string=makeTableString(sm.df[2, 1], inMean=sm.df[2, "median"],  sm.df[2, "mad"], sm.df[2, "min"], sm.df[2, "max"], sm.df[2, "nacount"], inMissingData=TRUE)
                mdd.string=makeTableString(sm.df[1, 1], inMean=sm.df[1, "median"],  sm.df[1, "mad"], sm.df[1, "min"], sm.df[1, "max"], sm.df[1, "nacount"], inMissingData=TRUE)
                
                ##name=paste(name, "†", sep="")
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
                if (doNonParametricTestsOnVariable(variable)) {
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
                ## print(test)
            } else {
                cat ("*** Insufficient opservations\n")
            } ## end of if (sm.df[1, "N"] > 3 && sm.df[2, "N"] > 3) {

            var.statistic=""
            var.df=""
            var.pvalue=""
            var.parameter=""
            var.significance=""
            var.effect.size=""
            
            if (is.list(test)) {
                var.statistic=round(test$statistic, 2)
                var.df=ifelse(is.numeric(test$parameter), round(test$parameter, 2), "")
                var.pvalue=round(test$p.value, 2)
                var.parameter=ifelse(is.null(test$parameter), "NA", round(test$parameter, 2))
                var.significance=make.significance.indications(test$p.value)

                if (doNonParametricTestsOnVariable(variable)) {
                    ## compute Probability of Superiority here
                    if (var.pvalue < 0.1) {

                        if ( 1 == 1 ) {
                            orddom.ps   =dmes     (na.omit(inData[inData$Grp=="NCL", variable]), na.omit(inData[inData$Grp=="MDD", variable]))
                            orddom.ps.ci=dmes.boot(na.omit(inData[inData$Grp=="NCL", variable]), na.omit(inData[inData$Grp=="MDD", variable]), theta.es="PSc")

                            var.effect.size=round(orddom.ps$PSc, 3)
                            ## upper bound on the 95% Confidence interval for the effect size
                            es.ci.lb=round(orddom.ps.ci$theta.bci.lo, 2)
                            ## upper bound on the 95% Confidence interval for the effect size                        
                            es.ci.ub=round(orddom.ps.ci$theta.bci.up, 2)
                        }
                        else {
                            var.effect.size=round(var.statistic / ( length(inData[inData$Grp=="NCL", variable]) *  length(inData[inData$Grp=="MDD", variable]) ), 2)
                            es.ci.ub=0
                            es.ci.lb=0
                        }
                        
                        st = sprintf("%s,%s,%s,W = %0.0f,%0.2f,PS=%s (%s; %s),%s",
                            name, paste(ctrl.string, "†", sep=""), paste(mdd.string, "†", sep=""),
                            round(var.statistic, 2),
                            round(var.pvalue, 2),  var.effect.size, es.ci.lb , es.ci.ub, var.significance)
                    } else {
                        st = sprintf("%s,%s,%s,W = %0.0f,%0.2f,,%s",
                            name, ctrl.string, mdd.string,
                            round(var.statistic, 2),
                            round(var.pvalue, 2), var.significance )
                    }
                } else {
                    if (var.pvalue < 0.1) {
                        es=tes(var.statistic, length(inData[inData$Grp=="NCL", variable]), length(inData[inData$Grp=="MDD", variable]))
                        ## upper bound on the 95% Confidence interval for the effect size
                        es.ci.lb=round(es$l.g, 2)
                        ## upper bound on the 95% Confidence interval for the effect size                        
                        es.ci.ub=round(es$u.g, 2)
                        var.effect.size=round(es$g, 2)
                        st = sprintf("%s,%s,%s,t(%0.2f) = %0.2f,%0.2f,g=%s (%0.2f; %0.2f),%s",
                            name, ctrl.string, mdd.string,
                            round(var.parameter, 2), round(var.statistic, 2),
                            round(var.pvalue, 2),    var.effect.size, es.ci.lb, es.ci.ub, var.significance )
                    } else {
                        st = sprintf("%s,%s,%s,t(%0.2f) = %0.2f,%0.2f,,%s",
                            name, ctrl.string, mdd.string,
                            round(var.parameter, 2), round(var.statistic, 2),
                            round(var.pvalue, 2),    var.significance )
                    }
                }
            } ## end of if (is.list(test)) {

            ## st=paste(name, ctrl.string, mdd.string,
            ##     var.parameter, var.statistic, var.pvalue, var.effect.size, var.significance, sep=",")
            push(mystack, st)
        } ## end of if (is.factor(inData[, variable])) { ... } else { ... }
    }

    cat("################################################################################\n");
    cat("Summary statistics table\n")  
    l=mystack$value()
    cat(header, "\n")
    for (i in 1:length(l)) {
        cat (l[[i]], "\n")
    }
    ##stop()
}

doNonParametricTestsOnVariable <- function (inVariableName) {

    nonParametricRegexp="Tanner|SES|PSWQ|CGI.CGAS|CoRum|RSQ|Hand"

    if ( inVariableName == "NS" ) {
        return (TRUE)
    } else if (inVariableName == "HA" ) {
        return (TRUE)
    } else if (inVariableName == "RD" ) {
        return (TRUE)
    } else if (inVariableName == "P" ) {
        return (TRUE)
    } else if (inVariableName == "SD" ) {
        return (TRUE)
    } else if (inVariableName == "C" ) {
        return (TRUE)
    } else if (inVariableName == "ST" ) {
        return (TRUE)
    } else if (any(grep(nonParametricRegexp, inVariableName)) ) {
        return (TRUE)
    } else {
        return (FALSE)
    }
}

##########################################################################################################################################################################
### END OF FUNCTIONS #####################################################################################################################################################
##########################################################################################################################################################################


if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    cat(paste("Sorry can't set data directories for this computer\n"))
}

Admin.data.dir=normalizePath(file.path(root.dir, "sanDiego/cPine/data/admin"))
Regressor.data.dir=normalizePath(file.path(root.dir, "sanDiego/cPine/data/regressors"))
Config.data.dir=normalizePath(file.path(root.dir, "sanDiego/cPine/data/config"))
Group.data.dir=normalizePath(file.path(root.dir, "sanDiego/cPine/data/Group.data"))
Group.results.dir=normalizePath(file.path(root.dir, "sanDiego/cPine/data/Group.results"))
scriptsDir=normalizePath(file.path(root.dir, "sanDiego/cPine/scripts"))

## Admin.data.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/admin"
## Regressor.data.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/regressors"
## Config.data.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/config"
## Group.data.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.data"
## Group.results.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.results"  
## scriptsDir="/Volumes/PROMISEPEGASUS/yangdata/cPine/script"

cleanSubjectsFilename=file.path(Config.data.dir, "clean.mdd.list.txt")
demographicsFilename=file.path(Admin.data.dir, "data_entry_current_021113.csv")
original.control.subjectList.filename=file.path(Config.data.dir, "control.subjectList.txt")
original.mdd.subjectList.filename=file.path(Config.data.dir, "mdd.subjectList.txt")
regressor.summary.filename=file.path(Regressor.data.dir, "summary.tab")

contrastName="allEmotiveVsNeutral"

ctrl.subjectOrder.filename=file.path(Group.data.dir, paste("subjectOrder.ctrlOnly", contrastName, "REML.csv", sep="."))
mdd.subjectOrder.filename=file.path(Group.data.dir, paste("subjectOrder.mddOnly",contrastName, "REML.csv", sep="."))

### Original subject lists
cat ("########## Reading the original subject lists\n")

cat("Reading original CONTROL subject list: ", original.control.subjectList.filename, "\n")
original.control.subjectList=read.table(original.control.subjectList.filename, header=FALSE)
colnames(original.control.subjectList)=c("subject")
cat(sprintf("Read %d subjects from the original CONTROL subject list\n", dim(original.control.subjectList)[1]))

cat("Reading original MDD subject list: ", original.mdd.subjectList.filename, "\n")
original.mdd.subjectList=read.table(original.mdd.subjectList.filename, header=FALSE)
colnames(original.mdd.subjectList)=c("subject")
cat(sprintf("Read %d subjects from the original MDD subject list\n", dim(original.mdd.subjectList)[1]))


### Subjects actually in the bucket files
cat ("########## Reading the bucklet lists\n")

cat("Reading list of CONTROLs subject in the bucket file: ", ctrl.subjectOrder.filename, "\n")
ctrl.subjectList=read.table(ctrl.subjectOrder.filename, header=TRUE)
colnames(ctrl.subjectList)=c("subject")
cat(sprintf("Read %d subjects from the CONTROL subject bucket list\n", dim(ctrl.subjectList)[1]))

cat("Reading list of MDDs subject in the bucket file: ", mdd.subjectOrder.filename, "\n")
mdd.subjectList=read.table(mdd.subjectOrder.filename, header=TRUE)
colnames(mdd.subjectList)=c("subject")
cat(sprintf("Read %d subjects from the MDD subject bucket list\n", dim(mdd.subjectList)[1]))


### Subjects that were not in the bucket files, either due to error or missing regressors etc.
cat ("########## The following subjects were in the original subjects list but were automatically excluded from bucket generation due to missing data\n")

nodata.ctrl.filename=file.path(Group.data.dir, paste("nodata.ctrlOnly", contrastName, "REML.txt", sep="."))
cat("Reading missing subject information file for CONTROLs: ", nodata.ctrl.filename, "\n")
if (file.info(nodata.ctrl.filename)$size > 0 ) {
    nodata.ctrl=read.table(nodata.ctrl.filename)
    colnames(nodata.ctrl)=c("subject")
} else {
    nodata.ctrl=data.frame("subject"=c())
}
## nodata.ctrl$subject=as.factor(gsub("_A", "", as.character(nodata.ctrl$subject), fixed=TRUE))
cat(sprintf("Read missing subject information for %s unique CONTROL subjects\n",  length(unique(nodata.ctrl$subject))))
cat("There was no data for the following CTRL subjects:\n")
cat(as.vector(nodata.ctrl$subject))
cat("\n")

nodata.mdd.filename=file.path(Group.data.dir, paste("nodata.mddOnly", contrastName, "REML.txt", sep="."))
cat("Reading missing subject information file for MDDs: ", nodata.mdd.filename, "\n")
if (file.info(nodata.mdd.filename)$size > 0 ) {
    nodata.mdd=read.table(nodata.mdd.filename)
    colnames(nodata.mdd)=c("subject")
} else {
    nodata.mdd=data.frame("subject"=c())
}
##nodata.mdd$subject=as.factor(gsub("_A", "", as.character(nodata.mdd$subject), fixed=TRUE))
cat(sprintf("Read missing subject information for %s unique MDD subjects\n",  length(unique(nodata.mdd$subject))))
cat("There was no data for the following MDD subjects:\n")
cat(as.vector(nodata.mdd$subject))
cat("\n")


### Subjects that were not in the bucket files, either due to too much motion.
cat ("########## The following subjects were in the original subjects list but were automatically excluded from bucket generation due to much motion\n")
doNotAnalyze.ctrl.filename=file.path(Group.data.dir, paste("doNotAnalyse.ctrlOnly", contrastName, "REML.txt", sep="."))
cat("Reading missing subject information file for CONTROLs: ", doNotAnalyze.ctrl.filename, "\n")
if (file.info(doNotAnalyze.ctrl.filename)$size > 0 ) {
    doNotAnalyze.ctrl=read.table(doNotAnalyze.ctrl.filename)
    colnames(doNotAnalyze.ctrl)=c("subject")
} else {
    doNotAnalyze.ctrl=data.frame("subject"=c())
}
## doNotAnalyze.ctrl$subject=as.factor(gsub("_A", "", as.character(doNotAnalyze.ctrl$subject), fixed=TRUE))
cat(sprintf("Read missing subject information for %s unique CONTROL subjects\n",  length(unique(doNotAnalyze.ctrl$subject))))
cat("There was no data for the following CTRL subjects:\n")
cat(as.vector(doNotAnalyze.ctrl$subject))
cat("\n")

doNotAnalyze.mdd.filename=file.path(Group.data.dir, paste("doNotAnalyse.mddOnly", contrastName, "REML.txt", sep="."))
cat("Reading missing subject information file for MDDs: ", doNotAnalyze.mdd.filename, "\n")
if (file.info(doNotAnalyze.mdd.filename)$size > 0 ) {
    doNotAnalyze.mdd=read.table(doNotAnalyze.mdd.filename)
    colnames(doNotAnalyze.mdd)=c("subject")
} else {
    doNotAnalyze.mdd=data.frame("subject"=c())
}
##doNotAnalyze.mdd$subject=as.factor(gsub("_A", "", as.character(doNotAnalyze.mdd$subject), fixed=TRUE))
cat(sprintf("Read missing subject information for %s unique MDD subjects\n",  length(unique(doNotAnalyze.mdd$subject))))
cat("There was no data for the following MDD subjects:\n")
cat(as.vector(doNotAnalyze.mdd$subject))
cat("\n")


## Create the subjectOrder variable used throughout the rest of the script
## the clean subject list has no timepoint suffix
subjectOrder=data.frame("subject"=c(as.vector(mdd.subjectList$subject), as.vector(ctrl.subjectList$subject)))
colnames(subjectOrder)=c("subject")
subjectOrder$subject=as.vector(subjectOrder$subject)
## due to some sort of screw up 169/300_A is the same subject so fix it here.
badSubjectIds=which (subjectOrder$subject=="300_A")
if (length(badSubjectIds) > 0 ) {
    subjectOrder[badSubjectIds, "subject"]="169/300_A"
}
subjectOrder$subject=as.factor(gsub("_A", "", as.character(subjectOrder$subject), fixed=TRUE))
cat(sprintf("*** The subjectOrder data frame has %s unique subjects\n",  length(unique(subjectOrder$subject))))

##stop("Stooooooooooooooooooooooooooooooooooooooooooooooooooping")

###The clean list
cat("Reading the clean subjects list from: ", cleanSubjectsFilename, "\n")
cleanSubjects=read.csv(cleanSubjectsFilename, header=FALSE)
colnames(cleanSubjects)=c("subject")
cat(sprintf("Read %d subjects from the clean subject  list\n", dim(cleanSubjects)[1]))

## due to some sort of screw up 169/300_A is the same subject so fix it here.
mdd.subjectList$subject=as.vector(mdd.subjectList$subject)

badSubjectIds=which (mdd.subjectList$subject=="300_A")
if (length(badSubjectIds) > 0 ) {
    mdd.subjectList[badSubjectIds, "subject"]="169/300_A"
}
##mdd.subjectList[which (mdd.subjectList$subject=="300_A"), "subject"]="169/300_A"

mdd.subjectList$subject=as.factor(gsub("_A", "", as.character(mdd.subjectList$subject), fixed=TRUE))

notOnCleanList=setdiff(as.vector(mdd.subjectList$subject), as.vector(cleanSubjects$subject))
cat(sprintf("The following %d subjects were not on the clean list\n", length(notOnCleanList)))
print(notOnCleanList)

###stop("Stooooooooooooooooooooooooooooooooooooooooooooooooooping")

## colnames(mdd.subjectOrder)=c("subject")
## colnames(ctrl.subjectOrder)=c("subject")
## mdd.subjectOrder$subject=as.factor(gsub("_A", "", as.character(mdd.subjectOrder$subject), fixed=TRUE))
## ctrl.subjectOrder$subject=as.factor(gsub("_A", "", as.character(ctrl.subjectOrder$subject), fixed=TRUE))

## nodata.ctrl.filename=file.path(Group.data.dir, "nodata.ctrlOnly.allRemVsAllNotrem.REML.txt")
## cat("*** Reading missing subject information file for controls: ", nodata.ctrl.filename, "\n")
## nodata.ctrl=read.table(nodata.ctrl.filename)
## colnames(nodata.ctrl)=c("subject")
## nodata.ctrl$subject=as.factor(gsub("_A", "", as.character(nodata.ctrl$subject), fixed=TRUE))
## cat(sprintf("*** Read missing subject information for %s unique subjects\n",  length(unique(nodata.ctrl$subject))))
## cat(sprintf("Before removing the CTRL subjects for which there is no data, there were %d subjects\n", length(unique(ctrl.subjectOrder$subject))))
## ctrl.subjectOrder=data.frame(ctrl.subjectOrder$subject[ ! ctrl.subjectOrder$subject %in% nodata.ctrl$subject])
## colnames(ctrl.subjectOrder)=c("subject")
## ctrl.subjectOrder$subject=drop.levels(ctrl.subjectOrder$subject)
## cat(sprintf("After removing the CTRL subjects for which there is no data, there are %d subjects\n", length(unique(ctrl.subjectOrder$subject))))
## cat("There was no data for the following CTRL subjects:\n")
## cat(as.vector(nodata.ctrl$subject))
## cat("\n")

## nodata.mdd.filename=file.path(Group.data.dir, "nodata.mddOnly.allRemVsAllNotrem.REML.txt")
## cat("*** Reading missing subject information file for MDDs: ", nodata.mdd.filename, "\n")
## nodata.mdd=read.table(nodata.mdd.filename)
## colnames(nodata.mdd)=c("subject")
## nodata.mdd$subject=as.factor(gsub("_A", "", as.character(nodata.mdd$subject), fixed=TRUE))
## cat(sprintf("*** Read missing subject information for %s unique subjects\n",  length(unique(nodata.mdd$subject))))
## cat(sprintf("Before removing the MDD subjects for which there is no data, there were %d subjects\n", length(unique(mdd.subjectOrder$subject))))
## mdd.subjectOrder=data.frame(mdd.subjectOrder$subject[ ! mdd.subjectOrder$subject %in% nodata.mdd$subject])
## colnames(mdd.subjectOrder)=c("subject")
## mdd.subjectOrder$subject=drop.levels(mdd.subjectOrder$subject)
## cat(sprintf("After removing the MDD subjects for which there is no data, there are %d subjects\n", length(unique(mdd.subjectOrder$subject))))    
## cat("There was no data for the following MDD subjects:\n")
## cat(as.vector(nodata.mdd$subject))
## cat("\n")

## subjectOrder=data.frame("subject"=c(intersect(mdd.subjectOrder$subject, cleanSubjects[cleanSubjects$GOOD=="TRUE", "ID"]), as.vector(ctrl.subjectOrder$subject)))

## colnames(subjectOrder)=c("subject")
## subjectOrder$subject=as.factor(gsub("_A", "", as.character(subjectOrder$subject), fixed=TRUE))
## cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(subjectOrder$subject))))

cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT"))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))
demographics$ID=as.factor(gsub("_[ABCD]", "", as.character(demographics$ID), fixed=FALSE))

cat("*** Reading", regressor.summary.filename, "\n")
regressorSummary=read.table(regressor.summary.filename, header=TRUE, sep="\t")
cat(sprintf("*** Read regressor summary data for %s unique subjects\n",  length(unique(regressorSummary$Subject))))
regressorSummary$Subject=as.vector(regressorSummary$Subject)

## drop all timepoints except for the A (baseline)
regressorSummary=regressorSummary[! grepl ("_[BCDEF]", regressorSummary$Subject), ]

badSubjectIds=which (regressorSummary$Subject=="300_A")
if (length(badSubjectIds) > 0 ) {
    regressorSummary[badSubjectIds, "Subject"]="169/300_A"
}

##regressorSummary[which (regressorSummary$Subject=="300_A"), "Subject"]="169/300_A"
regressorSummary$Subject=as.factor(gsub("_[ABCD]", "", as.character(regressorSummary$Subject), fixed=FALSE))


##nodata=data.frame(subject=c(107, 108, 124, 146, 148, 153, 303, 308, 341, 347, 357, 117, 161, "169/300", 313, 316, 340, 345, 344))

selectedColumns=c(
    "Grp", "Gender", "DOB", "MRI", "Race", "Hand", "SES", "Tanner1", "Tanner2", "TannerAvg",
    "CGI.CGAS",                 "CDRS.raw",                  "CDRS.tscore",               "SFMQ.total",              
    "MADRS.raw",                "WASI.PERF",                 "WASI.Full.4",               "PSWQ", 			"CoRum",
    "RADS.DM",                  "RADS.AN",                   "RADS.NS",                   "RADS.SC",
    "RADS.DM.Tscore",           "RADS.AN.Tscore",            "RADS.NS.Tscore",            "RADS.SC.Tscore",  ##"RADS.Total.Tscore",
    "BDI.II",                   "CDI",                       "MASC.total",                "MASC.tscore",           "RSQ.fixed",
    "SCARED.som.panic",         "SCARED.gen.anx",            "SCARED.seper.anx",         
    "SCARED.soc.phobia",        "SCARED.school.phobia",      "BIS",                       "BAS.drive",                
    "BAS.funseek",              "BAS.reward.responsiveness", "PsycAche.total", "NS",	"HA",	"RD",	"P",	"SD",	"C",	"ST")

if (length(setdiff(selectedColumns, colnames(demographics))) != 0) {
    message ("*** The following column name(s) is causing problmes: ", setdiff(selectedColumns, colnames(demographics)), "\n")
    stop("There is a mismatch between the selected columns and the columns in the demographics data frame. Stopping.")
    
}

####################################################################################################
### Set up the data frames to work out how many were excluded for motion and missing data
####################################################################################################

original.subjectList=rbind(original.control.subjectList, original.mdd.subjectList)
original.subjectList$subject=gsub("_A", "", as.character(original.subjectList$subject), fixed=TRUE)
## due to some sort of screw up 169/300_A is the same subject so fix it here.
original.subjectList[which (original.subjectList$subject=="300"), ]="169/300"
original.subjectList.mgd=cbind(original.subjectList, demographics[match(original.subjectList$subject, demographics$ID), c("Grp", "Gender")])
original.subjectList.mgd$subject=as.factor(original.subjectList.mgd$subject)
original.subjectList.mgd$Grp=drop.levels(original.subjectList.mgd$Grp)
original.subjectList.mgd$Gender=drop.levels(original.subjectList.mgd$Gender)

nodata.subjectList=rbind(nodata.ctrl, nodata.mdd)
nodata.subjectList$subject=gsub("_A", "", as.character(nodata.subjectList$subject), fixed=TRUE)
## due to some sort of screw up 169/300_A is the same subject so fix it here.
if (length(which (nodata.subjectList$subject=="300")) > 0 ) {
    nodata.subjectList[which (nodata.subjectList$subject=="300"), ]="169/300"
}
nodata.subjectList.mgd=cbind(nodata.subjectList, demographics[match(nodata.subjectList$subject, demographics$ID), c("Grp", "Gender")])
nodata.subjectList.mgd$subject=as.factor(nodata.subjectList.mgd$subject)
nodata.subjectList.mgd$Grp=drop.levels(nodata.subjectList.mgd$Grp)
nodata.subjectList.mgd$Gender=drop.levels(nodata.subjectList.mgd$Gender)


doNotAnalyze.subjectList=rbind(doNotAnalyze.ctrl, doNotAnalyze.mdd)
doNotAnalyze.subjectList=rbind(doNotAnalyze.ctrl, doNotAnalyze.mdd)
doNotAnalyze.subjectList$subject=gsub("_A", "", as.character(doNotAnalyze.subjectList$subject), fixed=TRUE)
## due to some sort of screw up 169/300_A is the same subject so fix it here.
if (length(which (doNotAnalyze.subjectList$subject=="300")) > 0 ) {
    doNotAnalyze.subjectList[which (doNotAnalyze.subjectList$subject=="300"), ]="169/300"
}
doNotAnalyze.subjectList.mgd=cbind(doNotAnalyze.subjectList, demographics[match(doNotAnalyze.subjectList$subject, demographics$ID), c("Grp", "Gender")])
doNotAnalyze.subjectList.mgd$subject=as.factor(doNotAnalyze.subjectList.mgd$subject)
doNotAnalyze.subjectList.mgd$Grp=drop.levels(doNotAnalyze.subjectList.mgd$Grp)
doNotAnalyze.subjectList.mgd$Gender=drop.levels(doNotAnalyze.subjectList.mgd$Gender)

cat("\n*** Table of individuals originally included in analysis\n")
original.gender.table=table(original.subjectList.mgd[, c("Gender", "Grp")])
##original.gender.test=prop.test(original.gender.table)
original.gender.table=addmargins(original.gender.table)
print(original.gender.table)

cat("\n*** Table of individuals dropped because of missing data\n")
nodata.gender.table=table(nodata.subjectList.mgd[, c("Gender", "Grp")])
##nodata.gender.test=prop.test(nodata.gender.table)
nodata.gender.table=addmargins(nodata.gender.table)
print(nodata.gender.table)

cat("\n*** Table of individuals dropped because of excessive motion\n")
doNotAnalyze.gender.table=table(doNotAnalyze.subjectList.mgd[, c("Gender", "Grp")])
##doNotAnalyze.gender.test=prop.test(doNotAnalyze.gender.table)
doNotAnalyze.gender.table=addmargins(doNotAnalyze.gender.table)
print(doNotAnalyze.gender.table)

####################################################################################################
### Set up the data frames for statistical analysis
####################################################################################################

mgd=cbind(subjectOrder,
    demographics    [match(subjectOrder$subject, demographics$ID), selectedColumns],
    regressorSummary[match(subjectOrder$subject, regressorSummary$Subject), -1])
mgd$Grp=drop.levels(mgd$Grp)
mgd$SCARED.som.panic=as.numeric(as.vector(mgd$SCARED.som.panic))
mgd$MASC.tscore=as.numeric(mgd$MASC.tscore)

### now fix the age

## this complicated looking regexp stuff cleans up the years in DOB
## and MRI with 4 digits to be just 2 digits
mgd$DOB=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", mgd$DOB)
mgd$MRI=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", mgd$MRI)

## inData$DOB=as.Date(inData$DOB, "%m/%d/%y")
## inData$MRI=as.Date(inData$MRI, "%m/%d/%y")
## inData$age.in.days=difftime(inData$MRI, inData$DOB, units="days")
## inData$age.in.days=as.numeric(inData$age.in.days)

## inData$age.in.weeks=difftime(inData$MRI, inData$DOB, units="weeks")
## inData$age.in.weeks=as.numeric(inData$age.in.weeks)
## inData$age.in.years=(inData$age.in.weeks)/52

mgd$DOB=as.Date(mgd$DOB, "%m/%d/%y")
mgd$MRI=as.Date(mgd$MRI, "%m/%d/%y")
mgd$age.in.days=difftime(mgd$MRI, mgd$DOB, units="days")
mgd$age.in.days=as.numeric(mgd$age.in.days)

mgd$age.in.weeks=difftime(mgd$MRI, mgd$DOB, units="weeks")
mgd$age.in.weeks=as.numeric(mgd$age.in.weeks)
mgd$age.in.years=(mgd$age.in.weeks)/52

### Now do the MASC scoring

mgd.dim=dim(mgd)
for (r in seq(1, mgd.dim[1]) ) {
    ## cat("##################################################\n")
    subjectNumber=mgd[r, "subject"]
    gender=mgd[r, "Gender"]
    age=round(mgd[r, "age.in.years"], 0)
    old.masc.tscore=mgd[r, "MASC.tscore"]

    ## cat(sprintf("r=%d subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f Old MASC tscore=%0.0f,\n", r, subjectNumber, gender, age, mgd[r, "MASC.total"], old.masc.tscore))
    
    new.masc.tscore=scoreMasc(gender, age, mgd[r, "MASC.total"])
    if (is.na(new.masc.tscore) ) {
        warning(sprintf ("Couldn't set a MASC tscore for subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f", subjectNumber, gender, age, mgd[r, "MASC.total"]))
    }
    
    mgd[r, "MASC.tscore"]=new.masc.tscore
    
    ## cat (sprintf("Old MASC tscore=%0.0f, new MASC tscore=%0.0f\n", old.masc.tscore, new.masc.tscore))
}

##stop("stopppppppppping")



## stop("no more ")

psychMeasures=list(
    list(variable="age.in.years", name="Age at time of scan (years)"),
    ## list(variable="Hand",         name="Edinburgh Handedness Inventory"),
    list(variable="SES",          name="Hollingshead Socioeconomic Score"),
    list(variable="TannerAvg",    name="Tanner Score"),
    ##list(variable="WASI.PERF",    name="Wechsler Abbreviated Scale of Intelligence (Performance)"),
    list(variable="WASI.Full.4",  name="Wechsler Abbreviated Scale of Intelligence"),
    
    list(variable="CGI.CGAS",     name="Children's Global Assessment Scale"),
    ## list(variable="PSWQ",         name="Penn State Worry Questionnaire"),
    ## list(variable="CoRum",        name="Corumination Questionnaire"),
    ## list(variable="RSQ.fixed",    name="Ruminative Responses Styles Questionnaire"),
    
    list(variable="CDRS.tscore",  name="Children's Depression Rating Scale (Standardized)"),
    
    ## list(variable="RADS.DM",      name="Reynolds Adolescent Depression Scale Dysphoric Mood"),
    ## list(variable="RADS.AN",      name="Reynolds Adolescent Depression Scale Anhedonia/Negative Affect"),
    ## list(variable="RADS.NS",      name="Reynolds Adolescent Depression Scale Negative Self-evaluation"),
    ## list(variable="RADS.SC",      name="Reynolds Adolescent Depression Scale Somatic Complaints"),    
    
    ## list(variable="RADS.DM.Tscore",      name="Reynolds Adolescent Depression Scale Dysphoric Mood (Standardized)"),
    ## list(variable="RADS.AN.Tscore",      name="Reynolds Adolescent Depression Scale Anhedonia/Negative Affect (Standardized)"),
    ## list(variable="RADS.NS.Tscore",      name="Reynolds Adolescent Depression Scale Negative Self-evaluation (Standardized)"),
    ## list(variable="RADS.SC.Tscore",      name="Reynolds Adolescent Depression Scale Somatic Complaints (Standardized)"),
    ##list(variable="RADS.Total.Tscore",      name="Reynolds Adolescent Depression Scale Total (Standardized)"),  
    
    
    ## list(variable="SFMQ.total", name="SFMQ"),
    ##list(variable="MADRS.raw",    name="Montgomery-Asberg Depression Rating Scale"),
    list(variable="BDI.II",       name="Beck Depression Inventory II"),
    list(variable="CDI",          name="Children's Depression Inventory"),
    ## list(variable="MASC.total", name="Multidimensional Anxiety Scale for Children"),
    list(variable="MASC.tscore", name="Multidimensional Anxiety Scale for Children (Standardized)")
    
    ## list(variable="SCARED.som.panic", name="SCARED Som. Panic"),
    ## list(variable="SCARED.gen.anx", name="SCARED Ganeralized Anxiety"),
    ## list(variable="SCARED.seper.anx", name="SCARED Seperation Anxiety"),
    ## list(variable="SCARED.soc.phobia", name="SCARED Social Phobia"),
    ## list(variable="BIS", name="Behavioral Inhibition System"),
    ## list(variable="BAS.drive", name="Behavioral Approach System: Drive"),
    ## list(variable="BAS.funseek", name="Behavioral Approach System: Fun Seeking"),
    ## list(variable="BAS.reward.responsiveness", name="Behavioral Approach System: Reward Responsiveness"),


    ## ## TCI related variables
    ## list(variable="NS", name="Novelty Seeking"),
    ## list(variable="HA", name="Harm Avoidance"),
    ## list(variable="RD", name="Reward Dependence"),
    ## list(variable="P",  name="Persistence"),
    ## list(variable="SD", name="Self-Directedness"),
    ## list(variable="C",  name="Cooperativeness"),
    ## list(variable="ST", name="Self-Transcendence")
    
    ## list(variable="PsycAche.total", name="PsycAche.total")
    )
stop()
cat("####################################################################################################\n")
cat("### Demographic and Psychiatric Statistics\n")

##sink("../data/Group.results/restingStateStatictics28Nov.txt")
## analyse(mgd, "../data/Group.results/restingStateDistributionOfGenders.png", "../data/Group.results/restingStateDistributionOfRaces.png")

##stop("check out the orddom package\n")

##analyse(mgd)

cat("####################################################################################################\n")
cat("### Behavioral Statistics\n")
## there were 50 stimuli in the Pine task (I think)
nStimuli=48

mystack <- stack()

countVariables=
    c("TotalRecallHits", "RecallFearfulHits", "RecallHappyHits", "RecallNeutralHits", "RecallSadHits",
      "TotalRecallMisses", "TotalRecallCorrectRejections", "TotalRecallFalseAlarms", "TotalRecallOmissions")
proportionVariables=
    c("TotalRecallHitsProportion", "RecallFearfulHitsProportion", "RecallHappyHitsProportion",
      "RecallNeutralHitsProportion", "RecallSadHitsProportion", "TotalRecallMissesProportion",
      "TotalRecallCorrectRejectionsProportion", "TotalRecallFalseAlarmsProportion", "TotalRecallOmissionsProportion")

behavioralData = mgd[, c("subject", "Grp", proportionVariables)]
proportionVariablesNumbers = grep ("Proportion", colnames(behavioralData), fixed=TRUE)

## now check for any subjects were the proportions are all zero, these
## should be dropped

## select those rows where all of the proportion variables are 0, these should be dropped
droppedSubjects = apply(behavioralData[, proportionVariablesNumbers] == 0, 1, all)
if (sum(droppedSubjects) > 0 ) {
    cat ("####################################################################################################\n")
    cat ("########## WARNING ### WARNING ### WARNING ### WARNING ### WARNING ### WARNING ### WARNING #########\n")
    cat ("####################################################################################################\n")    
    
    selectedSubjects = ! droppedSubjects
    cat("The following subjects were dropped because they had all 0 proportions\n")
    print(as.character(behavioralData[droppedSubjects, "subject"]))
    print(as.character(behavioralData[droppedSubjects, "Grp"]))
    cat ("####################################################################################################\n")
    cat ("########## WARNING ### WARNING ### WARNING ### WARNING ### WARNING ### WARNING ### WARNING #########\n")
    cat ("####################################################################################################\n")    
}
## now check for behavioral 
## stop()
## qnorm can't handle zeros (0) or 1 (one) so replace them with 0.01
## and 0.99, respectively
for ( variable in c("TotalRecallHitsProportion", "RecallFearfulHitsProportion", "RecallHappyHitsProportion",
                    "RecallNeutralHitsProportion", "RecallSadHitsProportion") ) { 
    if (any(behavioralData[, variable]==0) ) {
        behavioralData[behavioralData[, variable]==0, variable] = 0.01
    }
    if (any(behavioralData[, variable]==1) ) {
        behavioralData[behavioralData[, variable]==1, variable] = 0.99
    }
    qnorm.variable=paste(variable, "qnorm", sep=".")
    behavioralData[, qnorm.variable] = qnorm(behavioralData[, variable] )
}
if (any(behavioralData$TotalRecallFalseAlarmsProportion==0) ) {
    behavioralData[behavioralData$TotalRecallFalseAlarmsProportion==0, "TotalRecallFalseAlarmsProportion"] = 0.01
}

behavioralData$TotalRecallFalseAlarms.qnorm=qnorm(behavioralData$TotalRecallFalseAlarms)


summarizeByGroup <- function (inVariable) {

    sm.df = summarySE(behavioralData, measure=inVariable, groupvars=c("Grp"), na.rm=TRUE)
    ctrl.string = paste(
        round(sm.df[sm.df$Grp=="NCL", inVariable], 2), ## mean
        " ± ",
        round(sm.df[sm.df$Grp=="NCL", "sd"],     2), ## standard deviation
        sep="")
    mdd.string = paste(
        round(sm.df[sm.df$Grp=="MDD", inVariable], 2), ## mean
        " ± ",
        round(sm.df[sm.df$Grp=="MDD", "sd"],     2), ## standard deviation
        sep="")
    ## the sub/gsub calls remove "Recall", extraneous white space, and
    ## replace each uppercase letter with itself prepended by a space
    text = paste (sub("^[ ]", "", sub("[ ]?Recall ", " ", gsub('([a-z])([A-Z])', "\\1 \\2", inVariable, fixed=FALSE) )),
        ctrl.string,
        mdd.string, sep=",")

    return(list("sm.df" = sm.df, "text"=text))

}

## push(mystack, "Absolute Counts")
## push(mystack, "Characteristic,Control,MDD")
## for ( variable in countVariables ) {
##     st=summarizeByGroup(variable)
##     push(mystack, st[["st"]])
## }
push(mystack, "Proportions")
tbl=table(behavioralData[, "Grp"])
push(mystack, sprintf("Characteristic,Control (n=%d),MDD (n=%d)", tbl["NCL"], tbl["MDD"]))
for ( variable in proportionVariables ) {
    st=summarizeByGroup(variable)

    ## only compute the dPrime measures for Totsal Recall hits,
    ## fearful, happy, neutral and sad hits
    if (any(grep("Fearful|Happy|Neutral|Sad|TotalRecallHits", variable))) {
        ## stats go here
        qnorm.variable  = paste(variable, "qnorm",  sep=".")
        dPrime.variable = paste(variable, "dPrime", sep=".")
        
        ##cat(sprintf("qnorm.variable = %s , dPrime.variable = %s\n", qnorm.variable, dPrime.variable))
        behavioralData[, dPrime.variable] = behavioralData[, qnorm.variable] - behavioralData$TotalRecallFalseAlarms.qnorm
    }
    ## behavioralData$deePrime=behavioralData$TotalRecallHitsProportion.qnorm - behavioralData$TotalRecallFalseAlarms.qnorm
    
    push(mystack, st[["text"]])
}

dPrime.variables = colnames(behavioralData)[grep("dPrime", colnames(behavioralData))]
dPrimeColumnNumbers = grep ("dPrime$", colnames(behavioralData), fixed=FALSE)

## list of the subjects being dropped

##################################################
## Subject selection
##################################################
##
## only keep subjects who performed well on the task, that is where
## dPrime for the various stimuli types is > 0. This would mean that
## they got more correctly remembered faces than false alarms (i.e.,
## faces indicated as remembered but were actually novel at the
## recall phase of the experiment
##
dropSubjectsWithDPrimtLessThanZero=FALSE
if (dropSubjectsWithDPrimtLessThanZero) {
    selectedSubjects = apply(behavioralData[, dPrimeColumnNumbers ] > 0, 1, all)
    droppedSubjects = !selectedSubjects
    noDroppedSubjects = sum(droppedSubjects)
    if ( noDroppedSubjects > 0 ) { 
        cat ("####################################################################################################\n")
        cat ("########## WARNING ### WARNING ### WARNING ### WARNING ### WARNING ### WARNING ### WARNING #########\n")
        cat ("####################################################################################################\n")    
        
        cat(sprintf("The following %d subjects were dropped because they had d' values < 0\n", length(as.character(behavioralData[droppedSubjects, "subject"]))))
        print(as.character(behavioralData[droppedSubjects, "subject"]))
        print(as.character(behavioralData[droppedSubjects, "Grp"]))
        cat("Number of dropped subjects per group: \n")
        print(table(behavioralData[droppedSubjects, "Grp"]))
        behavioralData = behavioralData[selectedSubjects, ]
        
        cat ("####################################################################################################\n")
        cat ("########## WARNING ### WARNING ### WARNING ### WARNING ### WARNING ### WARNING ### WARNING #########\n")
        cat ("####################################################################################################\n")    
    }
}
##
##################################################

push(mystack, "d'")
tbl=table(behavioralData[, "Grp"])
push(mystack, sprintf("Characteristic,Control (n=%d),MDD (n=%d),Statistic,p value,Significance", tbl["NCL"], tbl["MDD"]))

for ( variable in dPrime.variables ) {

    sm.df = summarySE(behavioralData, measure=variable, groupvars=c("Grp"), na.rm=TRUE)
    ctrl.string = paste(
        round(sm.df[sm.df$Grp=="NCL", variable], 1), ## mean
        " ± ",
        round(sm.df[sm.df$Grp=="NCL", "sd"],     1), ## standard deviation
        sep="")
    mdd.string = paste(
        round(sm.df[sm.df$Grp=="MDD", variable], 1), ## mean
        " ± ",
        round(sm.df[sm.df$Grp=="MDD", "sd"],     1), ## standard deviation
        sep="")
    ## the sub/gsub calls remove "Recall", extraneous white space, and
    ## replace each uppercase letter with itself prepended by a space
    text = paste (sub("^[ ]", "", sub("[ ]?Recall ", " ", gsub('([a-z])([A-Z])', "\\1 \\2", variable, fixed=FALSE) )),
        ctrl.string,
        mdd.string, sep=",")
    
    test = 0
    if( inherits(test <- try(t.test(behavioralData[behavioralData$Grp=="MDD", variable],
                                    behavioralData[behavioralData$Grp=="NCL", variable]),
                             silent=FALSE),
                 "try-error") ) {
        cat("Got an error, setting test to 0\n")
        test = 0
    }
    stat = ""
    if (is.list(test)) {
        text = sprintf("%s,t(%0.2f) = %0.2f,%0.2f,%s",
            text, round(test$parameter, 2), round(test$statistic, 2),
            round(test$p.value, 2), make.significance.indications(test$p.value) )
    }
    
    push(mystack, text)
}

l=mystack$value()
for (i in 1:length(l)) {
    cat (l[[i]], "\n")
}

##stop()

my.base.size=18   
my_theme=
    theme_bw(base_size =  my.base.size) +
    theme(
        legend.position="bottom",
        ##panel.grid.major = element_blank(),
        ##panel.grid.minor = element_blank(),
        ##axis.title.x=element_blank(),
        axis.title.x = element_text(size=my.base.size, vjust=0),
        axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
        plot.title=element_text(size=my.base.size*1.2, vjust=1))

## proportionColumnNumbers = grep ("Proportion$", colnames(behavioralData), fixed=FALSE)
## for ( ty in  proportionColumnNumbers ) {
##     for ( group in levels(behavioralData$Grp) ) {
##         cat ("####################################################################################################\n")
##         boxout=boxplot.stats(behavioralData[behavioralData$Grp==group, ty])$out
##         outname=as.character(behavioralData[behavioralData$Grp==group, "subject"])
##         outliers=outname[behavioralData[behavioralData$Grp==group, ty] %in% boxout]
##         cat (sprintf("The following subjects are outliers for %s in the %s group\n", colnames(behavioralData)[ty], group))
##         if (length(outliers) == 0) {
##             cat("No outliers\n")
##         } else {
##             print(outliers)
##             print(boxout)
##         }
##     }
## }

                                        # now create a graph of the proportions of hits, omissions, false alarms, and correct rejections

## melted.behavioral=melt(behavioralData[, c(1, 2, proportionColumnNumbers)], id.vars=c("subject", "Grp"), measure.vars=colnames(behavioralData)[proportionColumnNumbers], variable_name="deeElements")
## graph=ggplot(melted.behavioral, aes(x=Grp, y=value, fill=deeElements)) +
##     geom_boxplot() +
##     scale_fill_brewer(palette="Set1",
##                       name="Proportions: ", 
##                       labels=gsub("Total ", "", gsub('([a-z])([A-Z])', "\\1 \\2",
##                           sub("Recall([A-Za-z]*)Proportion", "\\1", colnames(behavioralData)[grep ("Proportion$", colnames(behavioralData), fixed=FALSE)])))) +
##     labs(y="Proportion") +
##     my_theme
## print(graph)


##dPrimeColumnNumbers = grep ("dPrime$", colnames(behavioralData), fixed=FALSE)
outliersList=c()
for ( ty in dPrimeColumnNumbers ) {
    ##cat ("####################################################################################################\n")
    bs = boxplot.stats(behavioralData[, ty])
    iqr = diff(bs$stats[c(2, 4)])
    coef=1.5
    outlierIndices <- if (!is.na(iqr)) {
        behavioralData[, ty] < (bs$stats[2L] - (coef * iqr)) | behavioralData[, ty] > (bs$stats[4L] + (coef * iqr))
    }
    noOutliers = sum(outlierIndices)
    
    cat (sprintf("### %s boxplot.stats says there are %d outliers, we think there are %d outliers\n", colnames(behavioralData)[ty], length(bs$out), noOutliers))
    outliers=as.character(behavioralData[outlierIndices, "subject"])
    
    if (length(outliers) == 0) {
        outliersList = c(outliersList, colnames(behavioralData)[ty])
    } else {
        cat ("####################################################################################################\n")
        cat (sprintf("The following subjects are outliers for %s\n", colnames(behavioralData)[ty]))
        cat ( sprintf ("Outlier cutoffs: lower bound = %.3f upper bound = %.3f\n", (bs$stats[2L] - (coef * iqr)), (bs$stats[4L] + (coef * iqr))))
        print (outliers)
        print (as.character(behavioralData[outlierIndices, "Grp"]))
        print (behavioralData[outlierIndices, ty])
    }
}

if ( length(outliersList)  > 0 ) {
    cat ("####################################################################################################\n")
    cat( sprintf ("No outliers were found for the following variables\n"))
    print (outliersList)
}

fixSubjectOrderTable <- function (inSubjectOrderTable) {
    inSubjectOrderTable$subject=gsub("_A", "", as.character(inSubjectOrderTable$subject), fixed=TRUE)
    inSubjectOrderTable[inSubjectOrderTable$subject=="300", "subject"] = "169/300"
    inSubjectOrderTable$subject=as.factor(inSubjectOrderTable$subject)

    return(inSubjectOrderTable)
}

roistats.filename=file.path(Group.results.dir, "roiStats.fwhm4.2.pine.mddAndCtrl.happyVsSad.group.F-value.txt")
roistats=read.table(roistats.filename, header=TRUE)
## now only keep the 3rd and 4th ROIs. These are respectively the
## Fusiform and insula.
roistats=roistats[, c("Mean_3", "Mean_4")]

subjectOrder.filename=file.path(Group.data.dir, "subjectOrder.mddAndCtrl.happyVsSad.REML.csv")
subjectOrder=read.table(subjectOrder.filename, header=TRUE)
subjectOrder=fixSubjectOrderTable(subjectOrder)

roistats=cbind(subjectOrder, roistats)

behavioralData=cbind(behavioralData,
    roistats    [match(behavioralData$subject, roistats$subject), c("Mean_3", "Mean_4")],
    demographics[match(behavioralData$subject, demographics$ID), c("RADS.SC.Tscore")])
colnames(behavioralData)[dim(behavioralData)[2]]="RADS.SC.Tscore"

behavioralData.forCorrelations=behavioralData[, c("subject", "Grp", "RecallHappyHitsProportion.dPrime", "RecallSadHitsProportion.dPrime", "RADS.SC.Tscore", "Mean_3", "Mean_4")]
behavioralData.forCorrelations=rename(behavioralData.forCorrelations, c("Mean_3"="L.Insula", "Mean_4"="L.Fusiform.Gyrus"))
behavioralData.forCorrelations.mddOnly=subset(behavioralData.forCorrelations, Grp=="MDD")

for (roi in c("L.Insula", "L.Fusiform.Gyrus")) {
    for (dp in c("RecallHappyHitsProportion.dPrime", "RecallSadHitsProportion.dPrime", "RADS.SC.Tscore") ) {
        cat("### Correlation of", dp, "with contrast of happy vs sad activity in the", roi, "\n")
        ct=cor.test(behavioralData.forCorrelations.mddOnly[, dp], behavioralData.forCorrelations.mddOnly[, roi], method="spearman")
        print(ct)
        graph=ggplot(behavioralData.forCorrelations.mddOnly, aes_string(x=dp, y=roi)) +
            geom_point() +
                labs(x=dp, y=sprintf("Happy - Sad in %s", roi)) +
                    my_theme
        ggsave(file.path(Group.results.dir, sprintf("%s-%s.scatterBoxplot.pdf", dp, roi)), graph)
    }
}

## rads.lfusiform.graph=ggplot(behavioralData.forCorrelations.mddOnly, aes_string(x="RADS.SC.Tscore", y="L.Fusiform.Gyrus")) +
##     geom_point() +
##     labs(x="RADS Somatic Complaints t-score", y="Happy - Sad in L Fusiform Gyrus") +
##     my_theme


## now create a graph of the d' values for hits (total, happy, fearful, neutral, and sad)
melted.behavioral=melt(behavioralData[, c(1, 2, dPrimeColumnNumbers)], id.vars=c("subject", "Grp"), measure.vars=colnames(behavioralData)[dPrimeColumnNumbers], variable_name="deeElements")
graph=ggplot(melted.behavioral, aes(x=Grp, y=value, fill=deeElements)) +
    geom_boxplot() +
    scale_fill_brewer(palette="Set1",
                      name="d': ", 
                      labels=gsub(" Hits", "",
                          gsub('([a-z])([A-Z])', "\\1 \\2",
                               sub("Recall([A-Za-z]*)Proportion.dPrime", "\\1", colnames(behavioralData)[grep ("dPrime$", colnames(behavioralData), fixed=FALSE)]))
                          )
                      ) +
    labs(y="d'") +
    ggtitle("Task Performance Boxplot") +
    my_theme
##quartz()
##print(graph)
ggsave(file.path(Group.results.dir, "taskPerformanceBoxplot.pdf"))


cat("### Between-group false alarm rate\n")
false.alarm.ttest =
    t.test(behavioralData[behavioralData$Grp=="MDD", "TotalRecallFalseAlarmsProportion"],
           behavioralData[behavioralData$Grp=="NCL", "TotalRecallFalseAlarmsProportion"]) 
print(false.alarm.ttest)

total.recall.ttest=t.test(melted.behavioral[melted.behavioral$deeElements=="TotalRecallHitsProportion.dPrime", "value"])
cat("### Total Recall hits\n")
print(total.recall.ttest)

between.group.total.recall.hits.ttest=
    t.test(melted.behavioral[melted.behavioral$deeElements=="TotalRecallHitsProportion.dPrime" & melted.behavioral$Grp=="MDD", "value"],
           melted.behavioral[melted.behavioral$deeElements=="TotalRecallHitsProportion.dPrime" & melted.behavioral$Grp=="NCL", "value"])
cat("### Between Group: Total Recall hits\n")
print(between.group.total.recall.hits.ttest)

melted.behavioral.emotionsOnly = melted.behavioral[melted.behavioral$deeElements!="TotalRecallHitsProportion.dPrime", ]
melted.behavioral.emotionsOnly$Grp = drop.levels(melted.behavioral.emotionsOnly$Grp)

emotion.recall.aov = aov(value ~ Grp*deeElements, data=melted.behavioral.emotionsOnly)
cat("### Emotion Recall ANOVA\n")
print(summary(emotion.recall.aov))
print(pairwise.t.test(melted.behavioral.emotionsOnly$value, melted.behavioral.emotionsOnly$deeElements))
cat("### Means for the deePrime values\n")
print(aggregate(value ~ deeElements, data=melted.behavioral.emotionsOnly, mean))


## dPrimeColumnNumbers = grep ("dPrime$", colnames(behavioralData), fixed=FALSE)
## outliersList=list()
## for ( ty in dPrimeColumnNumbers ) {
##     for ( group in levels(behavioralData$Grp) ) {
##         cat ("####################################################################################################\n")
##         bs = boxplot.stats(behavioralData[behavioralData$Grp==group, ty])
##         iqr = diff(bs$stats[c(2, 4)])
##         coef=1.5
##         outlierIndices <- if (!is.na(iqr)) {
##             behavioralData[behavioralData$Grp == "MDD", ty] < (bs$stats[2L] - (coef * iqr)) | behavioralData[behavioralData$Grp == "MDD", ty] > (bs$stats[4L] + (coef * iqr))
##         }
##         ## cat ( sprintf ("lower = %.3f upper= %.3f\n", (bs$stats[2L] - (coef * iqr)), (bs$stats[4L] + (coef * iqr))))
##         ## print(outlierIndices)
##         ## print(behavioralData[behavioralData$Grp == "MDD", ty])
##         ## print(behavioralData[behavioralData$Grp == "MDD", "subject"])
##         ## cat("bs$out:")
##         ## print(bs$out)
##         noOutliers = sum(outlierIndices)

##         cat (sprintf("### %s %s boxplot.stats says there are %d outliers, we think there are %d outliers\n", colnames(behavioralData)[ty], group, length(bs$out), noOutliers))
##         ## outlierIndices = behavioralData[behavioralData$Grp == "MDD", ty] < bs$stat[1] | behavioralData[behavioralData$Grp == "MDD", ty] > bs$stat[5]
##         a = behavioralData[behavioralData$Grp == "MDD", c(1, ty)]
##         ## print(a)
##         ## print (a[outlierIndices, "subject"])
##         outliers=as.character(a[outlierIndices, "subject"])

##         if (length(outliers) == 0) {
##             ## cat("No outliers\n")
##             outliersList[[group]] = c(outliersList[[group]], colnames(behavioralData)[ty])
##         } else {
##             cat ("####################################################################################################\n")
##             cat (sprintf("The following subjects are outliers for %s in the %s group\n", colnames(behavioralData)[ty], group))            
##             print (outliers)
##             ##print (boxout)
##         }
##     }
## }
## for (group in names(outliersList) ) {
##     if ( length(outliersList[[group]])  > 0 ) {
##         cat ("####################################################################################################\n")
##         cat( sprintf ("No outliers were found for the following variables in the %s group\n", group))
##         print (outliersList[[group]])
##     }
## }

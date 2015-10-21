rm(list=ls())
graphics.off()

library(gdata)
library(reshape)
library(ggplot2)
library(robustbase)
library(MASS)
##source('pcor.R')

source("scoreMasc.r")

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

capwords <- function(s, strict = FALSE) {
    cap <- function(s) paste(toupper(substring(s,1,1)),
                             {s <- substring(s,2); if(strict) tolower(s) else s},
                             sep = "", collapse = " " )
    sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

make.significance.indications <- function(pValues, which.pValues=c(1)) {

    Signif=symnum(pValues, corr=FALSE, na=FALSE, cutpoints = c(0,  .001,.01, .05, .1, 1),
        symbols   =  c("***", "**", "*", ".", " "))
    f=format(Signif)

    ## only return the first one as we're only interested in marking significant group effects
    return(f[which.pValues])
}

substituteShortLabels <- function(inLevel) {
    returnSubstitutedLabel = gsub("[0-9]+ ", "", gsub("Inf", "Inferior", gsub("Sup", "Superior", gsub("Gy", "Gyrus",
        gsub("^R", "Right", gsub("^L", "Left", inLevel, fixed=FALSE), fixed=FALSE), fixed=TRUE), fixed=TRUE), fixed=TRUE))
    
    return (returnSubstitutedLabel)
}

stderror <- function(x) sd(x)/sqrt(length(x))

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
                       c( N    = length2(xx[,col], na.rm=na.rm),
                         mean = mean   (xx[,col], na.rm=na.rm),
                         sd   = sd     (xx[,col], na.rm=na.rm)
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

makePublicationTable <- function(inClusterWhereAmI, inClusters, inRoistats, inCom=TRUE) {
    hemisphere=gsub("[^RL]", "", substr(inClusterWhereAmI, 1, 1))
    ##print(hemisphere)
    if ( inCom ) {
        locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "CM RL", "CM AP", "CM IS")], 0))
    } else {
        locations=cbind(gsub("^[RL] ", "", inClusterWhereAmI), hemisphere, round(inClusters[, c("Volume", "MI RL", "MI AP", "MI IS")], 0))
    }
    
    ##print(locations)
    ##agg=aggregate(inRoistats[, grep("Mean", colnames(inRoistats))], list(inRoistats$Group), mean)
    ##print(agg)
    ##mns=round(t(agg[, -1]), 2)
    ##print(mns)
    ##colnames(mns)=levels(agg[,1])
    ##pubTable=cbind(locations, mns)
    pubTable=locations  
    ##print(pubTable)
    colnames(pubTable)=c("Structure", colnames(pubTable)[-1])
    rownames(pubTable)=NULL
    ##print(pubTable)
    return(pubTable)
}

loadPctChangeRoiStats <- function (inContrast) {
    ## stores the names of the percentage change file names to be read    
    fnames=c()
    ## stores the names of the stimuli that went into computing the
    ## contrast. this is used to add a column to the roistats
    ## data.frame that identifies each row as comming from a
    ## particular %cs score for the corresponding stimulus
    stims=c()
    for ( stimulus in stimuli[[inContrast]] )  {
        roistats.filename=
            file.path(Group.results.dir, sprintf("roiStats.pctChange.%s.stimulus.fwhm%0.1f.%s.%s.%s.%s.contrast.%s.txt", stimulus, usedFwhm, task, groups, contrast, fLabel, groups))
        
        if ( file.exists (roistats.filename) ) {
            fnames=c(fnames, roistats.filename)
            stims=c(stims, stimulus)
        } else {
            cat(sprintf("*** No such file : %s\n", roistats.filename))
        }
    }
    roistats=do.call(rbind, lapply(fnames,
        function(X) {
            data.frame(read.table(X, header=T, sep=""))
        }
        ))
    
    roistats$stimulus=as.factor(rep(capwords(stims), each=dim(roistats)[1]/length(stims)))
    roistats=roistats[, -1]
    
    return(roistats)
}

makeGraphTitle <- function ( inContrast, inSeedCount, inSeed ) {
    
    return( 
        sprintf("%s: %d. %s Seed", 
                paste (capwords(unlist(strsplit(contrast, "Vs", fixed=TRUE))), collapse=" Vs. "),
                inSeedCount,
                inSeed)
        )
}

getSeedRoiId<- function (inSeedName) {
    return( as.integer(sub("roi([0-9]).*", "\\1", inSeedName, perl=TRUE)))
}

fixDates <- function (inData) {
    ## this complicated looking regexp stuff cleans up the years in DOB
    ## and MRI with 4 digits to be just 2 digits
    ## month day year
    inData$DOB=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", inData$DOB)
    inData$MRI=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", inData$MRI)
    
    ## now convert to year/month/day
    inData$DOB=sub("([0-9]{1,2})/([0-9]{1,2})/([0-9]{2})", "\\3/\\1/\\2", inData$DOB)
    inData$MRI=sub("([0-9]{1,2})/([0-9]{1,2})/([0-9]{2})", "\\3/\\1/\\2", inData$MRI)
    
    inData$DOB=as.Date(inData$DOB, "%y/%m/%d")
    inData$MRI=as.Date(inData$MRI, "%y/%m/%d")

    return(inData)
}

computeAge <- function(inData) {
    
    age.in.weeks=difftime(inData$MRI, inData$DOB, units="weeks")
    age.in.weeks=as.numeric(age.in.weeks)

    inData$age.in.years=age.in.weeks/52

    return(inData)
}



os=system('uname -s', intern=T)
if (os == "Darwin") {
    admin.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/admin/")
    Config.data.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/config"
    Ppi.seeds.data.dir=file.path(Config.data.dir, "ppi_seeds")
    
    Group.data.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.data/")
    Group.results.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.results")
    Group.results.ppi.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.results/ppi")    
    Subject.data=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/")
} else {
    admin.dir=file.path("/data/sanDiego/cPine/data/admin/")
    Config.data.dir="/data/sanDiego/cPine/data/config"
    Ppi.seeds.data.dir=file.path(Config.data.dir, "ppi_seeds")
    
    Group.data.dir=file.path("/data/sanDiego/cPine/data/Group.data/")
    Group.results.dir=file.path("/data/sanDiego/cPine/data/Group.results")
    Group.results.ppi.dir=file.path("/data/sanDiego/cPine/data/Group.results/ppi")    
    Subject.data=file.path("/data/sanDiego/cPine/data/")

    ## Group.data.dir=file.path("/mnt/nfs/yangdata/restingstate/data/Group.data/")
    ## Group.results.ppi.dir=file.path("/mnt/nfs/yangdata/restingstate/data/Group.results/")
    ## Subject.data=file.path("/mnt/nfs/yangdata/restingstate/data/")
}

clust.header = c("Volume", "CM RL", "CM AP", "CM IS", "minRL",
    "maxRL", "minAP", "maxAP", "minIS", "maxIS", "Mean", "SEM", "Max Int",
    "MI RL", "MI AP", "MI IS")

## the vector of stimuli that were presented to the subject
## stimuli=c("happy", "fearful", "sad", "neutral", "happyRem", "happyNotrem", "fearfulRem", "fearfulNotrem", "sadRem", "sadNotrem", "neutralRem", "neutralNotrem")
stimuli=list(
    "fearfulVsHappy"            = c("fearful",    "happy"),
    "fearfulVsNeutral"          = c("fearful",    "Neutral"),
    "fearfulVsSad"              = c("fearful",    "sad"),
    "happyVsNeutral"            = c("happy",      "neutral"),
    "happyVsSad"                = c("happy",      "sad"),
    "neutralVsSad"              = c("neutral",    "sad"),
    "allEmotiveVsNeutral"       = c("fearful",    "happy", "sad", "neutral"),
    "happyRemVsHappyNotrem"     = c("happyRem",   "happyNotrem"),
    "fearfulRemVsFearfulNotrem" = c("fearfulRem", "fearfulNotrem"),
    "neutralRemVsNeutralNotrem" = c("neutralRem", "neutralNotrem"),
    "sadRemVsSadNotrem"         = c("sadRem",     "sadNotRem"),
    "allRemVsAllNotrem"         = c("happyRem",   "happyNotrem", "fearfulRem", "fearfulNotrem", "sadRem", "sadNotrem", "neutralRem", "neutralNotrem")
    )

##contrasts=c("fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad", "allEmotiveVsNeutral", "happyRemVsHappyNotrem",
##    "fearfulRemVsFearfulNotrem", "neutralRemVsNeutralNotrem", "sadRemVsSadNotrem", "allRemVsAllNotrem")

## contrasts=c("fearfulVsHappyExtractedRightSgAcc", "fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad")
## ##a.priori.rois=c("HarvardOxford-sub-maxprob-thr25-3mm-left-amygdala", "HarvardOxford-sub-maxprob-thr25-3mm-right-amygdala", "sgacc.left.3mm", "sgacc.right.3mm")

## contrasts=c("happyVsNeutral", "happyVsSad", "neutralVsSad")

contrasts=c("happyVsSad")

## this list is used to pick out only those ROIs from the
## between-group analysis that should have their connected clusters
## correlated against the psych measures

contrastToRoiNumbersList=list(
    "happyVsNeutral" = c(3),
##    "happyVsSad"     = c(3, 4),
    "happyVsSad"     = c(4),    
    "neutralVsSad"   = c(3)
    )
allSeeds=FALSE
contrastsToSuffices=list()
for ( contrast in contrasts ) {
    seedListFile=file.path(Ppi.seeds.data.dir, paste(contrast, "seeds.txt", sep="."))
    numberOfSeeds=strsplit(system(paste("wc", "-l", seedListFile), intern=TRUE), "[ ]+", fixed=FALSE)[[1]][2]

    contrastsToSuffices[[contrast]] = c(
                           if (allSeeds) {
                               sapply(seq(1, numberOfSeeds),
                                      function(ii) {
                                          sprintf("roi%d.seed.%s.%s", ii, contrast, "ROIinteraction")
                                      })
                           } else {
                               sapply(contrastToRoiNumbersList[[contrast]],
                                      function(ii) {
                                          sprintf("roi%d.seed.%s.%s", ii, contrast, "ROIinteraction")
                                      })
                           },
                           if (exists("a.priori.rois") ) {
                               apply(expand.grid(contrast, a.priori.rois), 1,
                                     function(x) {
                                         return(sprintf("%s.seed.%s.ROIinteraction", x[2], x[1]))
                                     })
                           }
                           )
}

##stop("Check the contrastsToSuffices list\n")

#contrasts=c("fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad")
## contrasts=c("fearfulVsSad", "happyVsSad", "neutralVsSad")

## a.priori.roi="apriori.rois.withVmpfc.3mm"

## contrastsToSuffices=list()
## for (contrast in contrasts) {
##     contrastsToSuffices[[contrast]] =
##         apply(expand.grid(contrast, a.priori.roi), 1,
##               function(x) {
##                   seedListFile=file.path(Ppi.seeds.data.dir, paste(x[1], x[2], "seeds.txt", sep="."))
##                   numberOfSeeds=strsplit(system(paste("wc", "-l", seedListFile), intern=TRUE), "[ ]+", fixed=FALSE)[[1]][2]
##                   unlist(lapply(seq.int(1, numberOfSeeds),
##                                 function (y) {
##                                     sprintf("roi%d.seed.%s%s.%s", y, x[1], x[2], "ROIinteraction")
##                                 } ))
##               } )
## }
##print(contrastsToSuffices)

##stop("Check suffixes list\n")



graphTitles=list(
    "fearfulVsHappyExtractedRightSgAcc" = "Fearful vs. Happy (sgACC)",
    "fearfulVsHappy" = "Fearful vs. Happy",
    "fearfulVsNeutral" = "Fearful vs. Neutral",
    "fearfulVsSad" = "Fearful vs. Sad",
    "happyVsNeutral" = "Happy vs. Neutral",
    "happyVsSad" = "Happy vs. Sad",
    "neutralVsSad" = "Neutral vs. Sad",
    "allEmotiveVsNeutral" = "( Fearful, Happy, & Sad ) vs. Neutral",
    "happyRemVsHappyNotrem" = "Happy Rememebred vs. Happy Not Remembered",
    "fearfulRemVsFearfulNotrem" = "Fearful Remembered vs. Fearful Not Remembered",
    "neutralRemVsNeutralNotrem" = "Neutral Remembered vs. Neutral Not Remembered",
    "sadRemVsSadNotrem" = "Sad Remembered vs. Sad Not Remembered",
    "allRemVsAllNotrem" = "All Remembered vs. All Not Remembered")

run.all=TRUE
run.male=FALSE
run.female=FALSE
####################################################################################################

groups="mddAndCtrl"
fLabel="group.F-value"
usedFwhm=4.2
task="pine"

## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
demographicsFilename=file.path(admin.dir, "data_entry_current_021113.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "#N/A", "#VALUE", "#VALUE!", "n/a"))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))

ksads.filename=file.path(admin.dir, "KSADSDiagnoses_tiffany.csv")
cat("Reading KSADS data from: ", ksads.filename, "\n")
ksads=read.csv(ksads.filename, header=TRUE, skip=1)
cat(sprintf("Read %d subjects from the KSADS file\n", dim(ksads)[1]))

regressionVariables=list(
    ##list(variable="CGI.CGAS",       name="Children's Global Assessment Scale"),
    ##list(variable="PSWQ",           name="Penn State Worry Questionnaire"),
    ##list(variable="CoRum",          name="Corumination Questionnaire"),

    ##list(variable="RSQ.fixed",      name="Ruminative Responses Questionnaire"),

    ##list(variable="CoRumII",        name="Corumination Questionnaire II"),

    list(variable="CDRS.tscore",    name="Children's Depression Rating Scale (Standardized)"),

    ## list(variable="RADS.DM.Tscore", name="Reynolds Adolescent Depression Scale Dysphoric Mood (Standardized)"),
    list(variable="RADS.AN.Tscore", name="Reynolds Adolescent Depression Scale Anhedonia/Negative Affect (Standardized)"),
    ## list(variable="RADS.NS.Tscore", name="Reynolds Adolescent Depression Scale Negative Self-evaluation (Standardized)"),
    list(variable="RADS.SC.Tscore", name="Reynolds Adolescent Depression Scale Somatic Complaints (Standardized)"),
    ## list(variable="MADRS.raw",      name="Montgomery-Asberg Depression Scale"),

    list(variable="BDI.II",         name="Beck Depression Inventory II"),
    ##list(variable="CDI",            name="Children's Depression Inventory"),

    list(variable="MASC.tscore",    name="Multidimensional Anxiety Scale for Children (Standardized)"),
    ## list(variable="BIS", name="Behavioral Inhibition System"),
    ## list(variable="BAS.drive", name="Behavioral Approach System Drive"),
    ## list(variable="BAS.funseek", name="Behavioral Approach System Fun Seeking"),
    ## list(variable="BAS.reward.responsiveness", name="Behavioral Approach System Reward Responsiveness"),

    ## ## TCI related variables
    ## list(variable="NS", name="Novelty Seeking"),
    ## list(variable="HA", name="Harm Avoidance"),
    ## list(variable="RD", name="Reward Dependence"),
    ## list(variable="P",  name="Persistence"),
    ## list(variable="SD", name="Self-Directedness"),
    ## list(variable="C",  name="Cooperativeness"),
    ## list(variable="ST", name="Self-Transcendence")

    list(variable="First.Onset", name="KSADS: Age at First Onset"),
    list(variable="Duration",    name="KSADS: Duration")
    
    )
## ## extract the list of variable names from the regressionVariables
## list rvs=regression variables select DOB, MRI, and MASC.total
## because we will need these later when we are calculating the MASC
## scores
rvs=c("DOB", "MRI", "MASC.total", unlist(regressionVariables)[grep ("variable", names(unlist(regressionVariables)))])

correlationCsvTableFilename=file.path(Group.results.ppi.dir, "correlationTable.csv")
if (file.exists(correlationCsvTableFilename)) {
    file.remove(correlationCsvTableFilename)
}
cat("Writing correlation table to", correlationCsvTableFilename, "\n")

correlationTableFilename=file.path(Group.results.ppi.dir, "correlationTable.txt")
if (file.exists(correlationTableFilename)) {
    file.remove(correlationTableFilename)
}
cat("Writing correlation to", correlationTableFilename, "\n")

my.base.size=12   
my_theme=
    theme_bw(base_size =  my.base.size) +
    theme(
        ## legend.position="none",
        ## panel.grid.major = element_line(color="white"),
        ## panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.x = element_text(size=my.base.size, vjust=0),
        axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
        plot.title=element_text(size=my.base.size*1.2, vjust=1))


## ajp_theme=
##     theme(base_size =  my.base.size,
##           ##legend.position="none",
##           panel.grid.major = element_blank(),
##           panel.grid.minor = element_blank(),
##           axis.title.x=element_blank(),
##           axis.title.x = element_text(size=my.base.size, vjust=0),
##           axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
##           plot.title=element_text(size=my.base.size*1.2, vjust=1))

mystack <- stack()
push(mystack, "Contrast,gender,seed,seed location,cluster,RL,AP,IS,RegressionVariable,S,DoF,pValue,Rho,Significance")

contrastCount=1

sink(correlationTableFilename)
for ( contrast in contrasts ) {

    cat("####################################################################################################\n")
    cat(sprintf("*** Graphing ROIs for the %s contrast\n", contrast))

    seed.clusterLocationsFilename=file.path(Group.results.dir, sprintf("clusterLocations.fwhm%0.1f.%s.%s.%s.%s.csv", usedFwhm, task, groups, contrast, fLabel))
    ##seed.clusterLocationsFilename=file.path(Group.results.dir, sprintf("clusterLocations.fwhm%0.1f.%s.%s.%s.%s.%s.csv", usedFwhm, task, groups, contrast, fLabel, a.priori.roi))
    cat("*** Reading seed cluster locations from", seed.clusterLocationsFilename, "\n")
    ## the gsub here chews up multiple consequtive spaces and replaces them with a single space
    seed.clusterWhereAmI=gsub(" +", " ", scan(file=seed.clusterLocationsFilename, what='character', sep=','))

    seedCount=1
    for ( seed in contrastsToSuffices[[contrast]] ) {

        seed.roi.id = getSeedRoiId(seed)
        
        cat(sprintf("*** Running correlations for seed: %s\n", seed))
        
        roistats.filename=file.path(Group.results.ppi.dir, sprintf("roiStats.fwhm%0.1f.%s.%s.%s.%s.txt", usedFwhm, task, groups, seed, fLabel))
        if(file.exists(roistats.filename)) {
            
            ## roistats contains the avergae from the seed in each ROI,
            ## you do not what to graph this
            roistats=read.table(roistats.filename, header=T, sep="")
            ## dump the first column as it's only the file name
            roistats=roistats[, -1]
            clusterCount=length(grep("Mean", colnames(roistats)))
            if (clusterCount > 0 ) {
                cat(sprintf("*** %d ** clusters found in %s\n", clusterCount, roistats.filename))
                
### Most of the following code up the the first long row of # is book-keeping to get the data frame in order
                
                clustersFilename=sprintf("clust.fwhm%0.1f.%s.%s.%s.%s.txt", usedFwhm, task, groups, seed, fLabel)
                cat("*** Reading", file.path(Group.results.ppi.dir, clustersFilename), "\n")
                clusters=read.table(file.path(Group.results.ppi.dir, clustersFilename))
                colnames(clusters) = clust.header
                ## this command contains the locations, as text, of the clusters and is the output of a perl script
                
                clusterLocationsFilename=file.path(Group.results.ppi.dir, sprintf("clusterLocations.fwhm%0.1f.%s.%s.%s.%s.csv", usedFwhm, task, groups, seed, fLabel))
                cat("*** Reading cluster locations from", clusterLocationsFilename, "\n")
                ## the gsub here chews up multiple consequtive spaces and replaces them with a single space
                clusterWhereAmI=gsub(" +", " ", scan(file=clusterLocationsFilename, what='character', sep=','))
                
                
                ## this file stores the order of the subjects in each of the following BRIK files
                subjectOrderFilename=file.path(Group.data.dir, paste("subjectOrder", groups, seed, "REML.z-score", "csv", sep="."))
                cat("*** Reading", subjectOrderFilename, "\n")
                subjectOrder=read.csv(subjectOrderFilename, header=T)
                subjectOrder$subject=gsub("_A", "", as.character(subjectOrder$subject), fixed=TRUE)
                subjectOrder[subjectOrder$subject=="300", "subject"] = "169/300"
                subjectOrder$subject=as.factor(subjectOrder$subject)
                cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(subjectOrder$subject))))
                
                ##roistats$subject=as.factor(as.vector(subjectOrder$subject))
                ##roistats$Group=demographics[match(roistats$subject, demographics$ID), c("Grp")]
                ##roistats$Group=drop.levels(roistats$Group)
                
                mgd=cbind(subjectOrder,
                    demographics[match(subjectOrder$subject, demographics$ID), c("Grp", "Gender", rvs[grep ("First\\.Onset|Duration", rvs, invert=TRUE)])]
                    )
                ## add this column with the subject ID concatenated with the timepoint
                ## so we can pull the correct data from the SLES data frame. This is
                ## done because the SLES data frame contains data for all time points.
                mgd$subjectWithTimePoint=as.factor(paste(mgd$subject, "A", sep="_"))

                mgd=cbind(mgd,
                    ksads[match(mgd$subjectWithTimePoint, ksads$Subject.Number), c("First.Onset", "Duration")])
                ## remove the now unnecessary subject with timepoint column
                mgd$subjectWithTimePoint=NULL

                mgd$Grp=drop.levels(mgd$Grp) 
                mgd=cbind(roistats, mgd)

                mgd=fixDates(mgd)
                mgd=computeAge(mgd)

                ##stop("check\n")
                
                if ( "MASC.tscore" %in% colnames(mgd) ) {
                    mgd$MASC.tscore=as.numeric(mgd$MASC.tscore)
                
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
                }
                cat("Dropping all subjects but MDD from the mgd data frame\n")
                mgd=mgd[mgd$Grp=="MDD", ]

                ## check if there is a RSQ.fixed column and make sure
                ## it is a vector not a factor. The latter will be the
                ## case if there are NAs in the demographics tables
                ## when it is read in. NAs will cause problems for
                ## cor.test, which required vector arguments not
                ## factors.

                if ( "RSQ.fixed" %in% colnames(mgd) ) {
                    mgd$RSQ.fixed = as.numeric(as.vector(mgd$RSQ.fixed))
                }

                ##stop("check mgd\n")
                
                melted.roistats=melt(mgd, id.vars=c("subject", "Grp", "Gender", rvs), measure.vars=paste("Mean_", seq(1, clusterCount), sep=""), variable_name="cluster")
                
                melted.roistats$cluster=factor(melted.roistats$cluster,
                    levels=c(paste("Mean_", seq(1, clusterCount), sep="")),
                    labels=paste(1:length(clusterWhereAmI), clusterWhereAmI))
                
                ## ## cat ("####################################################################################################\n")
                ## ## cat ("### roistats START\n")
                ## ## print(roistats)
                ## ## cat ("### roistats END\n")
                ## ## cat ("####################################################################################################\n")
                ## if (seed == "allRemVsAllNotrem") {
                ##     ## roistats$memory=ifelse(grepl("Notrem", roistats$stimulus), "Forgotten", "Remembered")
                ##     ## melted.roistats=melt(roistats, id.vars=c("subject", "Group", "memory"), measure.vars=paste("Mean_", seq(1, clusterCount), sep=""), variable_name="cluster")
                ## } else if (seed == "allEmotiveVsNeutral") {
                ##     ## roistats$emotion=ifelse(grepl("Neutral", roistats$stimulus), "Neutral", "Emotive")              
                ##     ## melted.roistats=melt(roistats, id.vars=c("subject", "Group", "emotion"), measure.vars=paste("Mean_", seq(1, clusterCount), sep=""), variable_name="cluster")
                ## } else {
                ##     melted.roistats=melt(roistats, id.vars=c("subject", "Group"), measure.vars=paste("Mean_", seq(1, clusterCount), sep=""), variable_name="cluster")
                ## }
                
                ## melted.roistats$cluster=factor(melted.roistats$cluster,
                ##     levels=c(paste("Mean_", seq(1, clusterCount), sep="")),
                ##     labels=paste(seq(1, clusterCount), clusterWhereAmI))
                ## ##labels=paste(1:length(clusterWhereAmI), clusterWhereAmI))
                
                ## stop("Check the data before correlations\n")
                
                for ( j in 1:length(regressionVariables ) ) {
                    ## clCounter=clusterCounter
                    clCounter=1
                    for ( level in levels(melted.roistats$cluster) ) {
                        cat("####################################################################################################\n")
                        cat("*** Seed & Contrast: ", seed, "\n")
                        cat("*** Running correlations for ", level, "on", regressionVariables[[j]]$name, "\n")
                        
                        mdl.df.all=   melted.roistats[melted.roistats$cluster %in% level, ]
                        mdl.df.male=  melted.roistats[melted.roistats$Gender=="M"   & melted.roistats$cluster %in% level, ]
                        mdl.df.female=melted.roistats[melted.roistats$Gender=="F" & melted.roistats$cluster %in% level, ]

                        ct.all=cor.test(mdl.df.all[, regressionVariables[[j]]$variable], mdl.df.all$value, method="spearman", na.action="na.omit")
                        cat("### BOTH GENDERS COMBINED\n")
                        print(ct.all)

                        ct.male=cor.test(mdl.df.male[, regressionVariables[[j]]$variable], mdl.df.male$value, method="spearman", na.action="na.omit")
                        cat("### MALE ONLY\n")
                        print(ct.male)

                        ct.female=cor.test(mdl.df.female[, regressionVariables[[j]]$variable], mdl.df.female$value, method="spearman", na.action="na.omit")
                        cat("### FEMALE ONLY\n")
                        print(ct.female)

                        
                        ##csvLine=sprintf("%s,%s,%s,%.2f,%.2f,%.5f,%.2f,%s",
                        ##  contrast, makeGraphTitle(level), regressionVariables[[j]]$name, ct$statistic, ifelse(is.null(ct$parameter), NA, ct$parameter),
                        ##  ct$p.value, ct$estimate, make.significance.indications(ct$p.value))
                        if ( run.all ) {
                            csvLine=sprintf("%s,%s,%s,%s,%s,%s,%s, %.2f,%.2f,%.5f,%.2f,%s",
                                ## contrast, "all", makeGraphTitle(level), paste(clusters[clCounter, c("CM RL", "CM AP", "CM IS")], collapse=","), regressionVariables[[j]]$name, ct.all$statistic,
                                contrast, "all", seed, seed.clusterWhereAmI[seed.roi.id], level, paste(clusters[clCounter, c("CM RL", "CM AP", "CM IS")], collapse=","), regressionVariables[[j]]$name, ct.all$statistic,                            
                                ifelse(is.null(ct.all$parameter), NA, ct.all$parameter),
                                ct.all$p.value, ct.all$estimate, make.significance.indications(ct.all$p.value))
                            cat(csvLine, "\n")
                            push(mystack, csvLine)
                        }

                        if ( run.male ) {
                            csvLine=sprintf("%s,%s,%s,%s,%s,%s,%s, %.2f,%.2f,%.5f,%.2f,%s",
                                ## contrast, "male", makeGraphTitle(level), paste(clusters[clCounter, c("CM RL", "CM AP", "CM IS")], collapse=","), regressionVariables[[j]]$name, ct.male$statistic,
                                contrast, "male", seed, seed.clusterWhereAmI[seedCount], level, paste(clusters[clCounter, c("CM RL", "CM AP", "CM IS")], collapse=","), regressionVariables[[j]]$name, ct.male$statistic,                            
                                ifelse(is.null(ct.male$parameter), NA, ct.male$parameter),
                                ct.male$p.value, ct.male$estimate, make.significance.indications(ct.male$p.value))
                            cat(csvLine, "\n")
                            push(mystack, csvLine)
                        }

                        if ( run.female ) { 
                            csvLine=sprintf("%s,%s,%s,%s,%s,%s,%s %.2f,%.2f,%.5f,%.2f,%s",
                                ##contrast, "female", makeGraphTitle(level), paste(clusters[clCounter, c("CM RL", "CM AP", "CM IS")], collapse=","), regressionVariables[[j]]$name, ct.female$statistic,
                                contrast, "female", seed, seed.clusterWhereAmI[seedCount], level, paste(clusters[clCounter, c("CM RL", "CM AP", "CM IS")], collapse=","), regressionVariables[[j]]$name, ct.female$statistic,                            
                                ifelse(is.null(ct.female$parameter), NA, ct.female$parameter),
                                ct.female$p.value, ct.female$estimate, make.significance.indications(ct.female$p.value))
                            cat(csvLine, "\n")
                            push(mystack, csvLine)
                        }
                        
                        clCounter=clCounter+1
                        cat("clCounter is now", clCounter, "\n")
                    } ## end of for ( level in levels(melted.roistats$cluster) )

                    ##stop("Stopping")
                } ## end of for ( level in levels(roistats.summary$cluster) )
                
                
            } ## end of if (clusterCount > 0 )
        } else {
            cat("### ", roistats.filename, "does not exist. Skipping\n")
        } ## end of if(file.exists(roistats.filename))
        seedCount=seedCount+1
    } ##     for ( seed in contrastsToSuffices[[contrast]]
    contrastCount=contrastCount + 1
} ## end of for (contrast in contrasts)

sink()

sink(file=correlationCsvTableFilename)
l=mystack$value()
for (i in 1:length(l)) {
    cat (l[[i]], "\n")
}
cat("\n")
sink()


## now remove all the NAs, if there are any
## outlier.matrix=outlier.matrix[rowSums(is.na(outlier.matrix))!=3, ]
## outlier.df=as.data.frame(outlier.matrix)
## colnames(outlier.df)=c("contrast", "cluster", "subject")

## outlier.df.count=
##   ddply(outlier.df, .variables=c("contrast", "subject"),
##   .fun=
##         function(xx, col) {
##           N=nrow(xx)
##         }
##         )
## ## Rename the "V1" column    
## outlier.df.count <- rename(outlier.df.count, c("V1"="count"))
## old.colnames=colnames(outlier.df.count)
## outlier.df.count=cbind(outlier.df.count,
##   mgd[match(outlier.df.count$subject, mgd$subject), c("Group")]
##   )
## colnames(outlier.df.count) = c(old.colnames, "Group")

## outlier.df.count.by.group=ddply(outlier.df.count, .variables=c("subject", "Group"),
##   .fun=
##   function(xx, col) {
##     N=nrow(xx)
##   }
##   )
## ## Rename the "V1" column    
## outlier.df.count.by.group <- rename(outlier.df.count.by.group, c("V1"="count"))

## ## now try to count the number of volumes censored for each subject in the table
## countCensoredVolumes <- function (inSubjectNumber) {

##   tbl=scan(file.path(Subject.data, paste(inSubjectNumber, "_A", sep=""), "functional", paste("censor_", inSubjectNumber, "_A_combined_2.1D", sep="")), what=integer())
##   censoredCount=length(tbl)-sum(tbl)

##   return(censoredCount)
## }

## ##apply(outlier.df.count.by.group$subject, 
## outlier.df.count.by.group$motionAndOutlierCensorCount=unlist(lapply(outlier.df.count.by.group$subject, countCensoredVolumes))

## #now print it out
## print(with(outlier.df.count.by.group, outlier.df.count.by.group[order(count, Group, decreasing = TRUE), ]))
## print(with(outlier.df.count.by.group, outlier.df.count.by.group[order(motionAndOutlierCensorCount, count, decreasing = TRUE), ]))

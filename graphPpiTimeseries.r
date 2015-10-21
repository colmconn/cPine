rm(list=ls())
graphics.off()

library(gdata)
library(reshape)
library(ggplot2)
library(robustbase)
library(MASS)
##source('pcor.R')

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

makeCleanedSeedTimeseriesFilename <- function( inSubjectNumber, inInfix) {
    return(file.path(Subject.data, inSubjectNumber, "functional", "ppi", paste(inSubjectNumber, task, inInfix , "contrast", "cleaned", "seedTimeseries", "1D", sep=".")))
}

makeCleanedClusterTimeseriesFilename <- function(inSubjectNumber, inInfix) {
    return(file.path(Subject.data, inSubjectNumber, "functional", "ppi", paste(inSubjectNumber, task, inInfix , "contrast", "cleaned", "clusterTimeseries", "1D", sep=".")))
}

makeAllInOnePpiGraphFilename <- function(inSubjectNumber, inInfix, inFileType="pdf") {
    return(file.path(Subject.data, inSubjectNumber, "functional", "ppi", paste(inSubjectNumber, task, inInfix , "contrast", "ppi", "allInOne", "connectedClusters", inFileType, sep=".")))
}

makeIndividualClusterPpiGraphFilename <- function(inSubjectNumber, inInfix, inClusterName, inFileType="pdf") {
    return(file.path(Subject.data, inSubjectNumber, "functional", "ppi", paste(inSubjectNumber, task, inInfix , "contrast", "ppi", gsub(" ", ".", inClusterName, fixed=TRUE) , "connectedClusters", inFileType, sep=".")))
}

makeRoiInteractionSuffix <- function(inInfix) {
    return(paste(inInfix , "ROIinteraction", sep="."))
}

readSubjectOrderFile <- function (inSeed) {
    ## this file stores the order of the subjects in each of the following BRIK files
    subjectOrderFilename=file.path(Group.data.dir, paste("subjectOrder", groups, makeRoiInteractionSuffix(inSeed), "REML.z-score", "csv", sep="."))
    cat("*** Reading", subjectOrderFilename, "\n")
    subjectOrder=read.csv(subjectOrderFilename, header=T)
    cat(sprintf("*** Read subject order data for %s unique subjects\n",  length(unique(subjectOrder$subject))))

    return(levels(subjectOrder$subject))
}


readDemograhicsFile <- function() {
    ## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
    demographicsFilename=file.path(admin.dir, "data_entry_current_021113.csv")
    cat("*** Reading", demographicsFilename, "\n")
    demographics=read.csv(demographicsFilename, header=T)
    cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))

    return(demographics)
}

readSeedLocations <- function (inContrast) { 
    seed.clusterLocationsFilename=file.path(Group.results.dir, sprintf("clusterLocations.fwhm%0.1f.%s.%s.%s.%s.csv", usedFwhm, task, groups, inContrast, fLabel))
    cat("*** Reading seed cluster locations from", seed.clusterLocationsFilename, "\n")
    ## the gsub here chews up multiple consequtive spaces and replaces them with a single space
    seed.clusterWhereAmI=gsub(" +", " ", scan(file=seed.clusterLocationsFilename, what='character', sep=','))

    return (seed.clusterWhereAmI)
}

getStudyGroupForSubject <- function(inDemographics, inSubject) {
    subject=gsub("_A", "", inSubject, fixed=TRUE)
    if (subject == "300" ) {
        subject = "169/300" 
    }
    
    group=demographics[match(subject, inDemographics$ID), c("Grp")]

    if ( ! group %in% c("NCL", "MDD") ) {
        stop(sprintf("Got group %s for subject %s. It should have been one of NCL, or MDD. Stopping.\n", group, subject))
    }

    return(group)
}

makeHtmlFiles <- function(inListOfFiles) {
    for (contrast in names(inListOfFiles) ) {
        for ( seed in names(inListOfFiles[[contrast]]) ) {

            htmlFilename=file.path(Group.results.ppi.dir, paste("ppi", contrast, seed, "html", sep="."))
            cat(sprintf("Making html file: %s\n", htmlFilename))
            cat("<html><body>\n", file=htmlFilename)

            ##cat(sprintf("Length of inListOfFiles[[contrast]][[seed]] is %d\n",  length(inListOfFiles[[contrast]][[seed]] ) ))
            if ( length(inListOfFiles[[contrast]][[seed]] ) > 0 ) { 
                for ( ii in seq.int(1, length(inListOfFiles[[contrast]][[seed]] ) ) ) {
                    cat( sprintf("<img src=\"%s\" />\n", inListOfFiles[[contrast]][[seed]][[ii]]), file=htmlFilename, append=TRUE)
                } ## end of for (subject in inListOfFiles[[contrast]][[seed]] ) {
            } else {
                cat(sprintf("<p>No clusters were functionally connected to %s in the %s contrast.</p>\n", seed, contrast), file=htmlFilename, append=TRUE) 
            }
            cat("</html></body>" , file=htmlFilename, append=TRUE)
        } ## end of for ( seed in names(inListOfFiles[[contrast]]) ) {
    } ## end of for (contrast in names(inListOfFiles) ) {
}

####################################################################################################
### Variables definitions are below this comment
####################################################################################################

os=system('uname -s', intern=T)
if (os == "Darwin") {
    admin.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/admin/")
    Config.data.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/config"
    Ppi.seeds.data.dir=file.path(Config.data.dir, "ppi_seeds")

    regressors.data.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/regressors"
    
    Group.data.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.data/")
    Group.results.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.results")
    Group.results.ppi.dir=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.results/ppi")    
    Subject.data=file.path("/Volumes/PROMISEPEGASUS/yangdata/cPine/data/")
} else {
    stop("Don't know where your directories are!")
}

## the vector of stimuli that were presented to the subject
contrastsToStimuliMap=list(
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

contrasts=c("fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad")
a.priori.rois=c("HarvardOxford-sub-maxprob-thr25-3mm-left-amygdala", "HarvardOxford-sub-maxprob-thr25-3mm-right-amygdala") ##, "anteriorInsula.3mm", "posteriorInsula.3mm", "sgacc.left.3mm", "sgacc.right.3mm")
##a.priori.rois=c("sgacc.left.3mm", "sgacc.right.3mm")

groups="mddAndCtrl"
fLabel="group.F-value"
usedFwhm=4.2
task="pine"

my.base.size=12   
my_theme=
    theme_bw(base_size =  my.base.size) +
    theme(
        ##legend.position="none",
                                        #panel.grid.major = element_line(color="white"),
        ##panel.grid.minor = element_blank(),
        ##axis.title.x=element_blank(),
        axis.title.x = element_text(size=my.base.size, vjust=0),
        axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
        plot.title=element_text(size=my.base.size*1.2, vjust=1))

####################################################################################################
### main code below this comment
####################################################################################################

contrastsToInfixes=list()
for ( contrast in contrasts ) {
    seedListFile=file.path(Ppi.seeds.data.dir, paste(contrast, "seeds.txt", sep="."))
    numberOfSeeds=as.integer(strsplit(system(paste("wc", "-l", seedListFile), intern=TRUE), "[ ]+", fixed=FALSE)[[1]][2])
    
    contrastsToInfixes[[contrast]]=list()
    
    contrastsToInfixes[[contrast]][["wholeBrainRois"]] = 
        sapply(seq(1, numberOfSeeds),
               function(ii) {
                   sprintf("roi%d.seed.%s", ii, contrast)
               })
    
    
    if (exists("a.priori.rois") ) {
        contrastsToInfixes[[contrast]][["aprioriRois"]] =
            apply(expand.grid(contrast, a.priori.rois), 1,
                  function(x) {
                      return(sprintf("%s.seed.%s", x[2], x[1]))
                  })
    } ## end of if (exists("a.priori.rois") ) {

} ## end of for ( contrast in contrasts ) {

## read the demographics file so we can get the group of the subject
demographics=readDemograhicsFile()
## listOfModels = list()
listOfFiles = list()

for ( contrast in contrasts) {

    cat("####################################################################################################\n")
    cat(sprintf("*** Graphing ROIs for the %s contrast\n", contrast))

    seed.clusterWhereAmI=readSeedLocations(contrast)
    ## listOfModels[[contrast]]=list()
    listOfFiles[[contrast]]=list()    
    
    seedCount=1
    ##for ( seed in contrastsToInfixes[[contrast]][["wholeBrainRois"]] ) {
    for ( seed in contrastsToInfixes[[contrast]][["aprioriRois"]] ) {        
        
        cat(sprintf("*** Graphing seed: %s\n", seed))
        subjects = readSubjectOrderFile(seed)

        ## listOfModels[[contrast]][[seed]]=list()
        listOfFiles[[contrast]][[seed]]=list()

        for (subjectNumber in subjects) {
            
            cat(sprintf("*** Graphing connected regions for subject: %s\n", subjectNumber))
            ## listOfModels[[contrast]][[seed]][[subjectNumber]]=list()
            
            ## read the seed timeseries
            cleanedSeedTimeseriesFilename = makeCleanedSeedTimeseriesFilename(subjectNumber, seed)
            if (file.exists (cleanedSeedTimeseriesFilename)) {
                cat(sprintf("*** Reading cleaned seed timeseries from: %s\n", cleanedSeedTimeseriesFilename))
                cleanedSeedTimeseries = read.table(cleanedSeedTimeseriesFilename, header=FALSE)
                colnames(cleanedSeedTimeseries) = seed.clusterWhereAmI[seedCount]

                ## now read the locations (as human readable text names) of the clusters connected to the seed as determined by the PPI analysis
                clusterLocationsFilename=file.path(Group.results.ppi.dir, sprintf("clusterLocations.fwhm%0.1f.%s.%s.%s.%s.csv", usedFwhm, task, groups, makeRoiInteractionSuffix(seed), fLabel))
                if ( file.exists(clusterLocationsFilename)) {
                    cat("*** Reading cluster locations from", clusterLocationsFilename, "\n")
                    ## the gsub here chews up multiple consequtive spaces and replaces them with a single space
                    clusterWhereAmI=gsub(" +", " ", scan(file=clusterLocationsFilename, what='character', sep=','))
                    
                    ## now read the timeseries of all the clusters to which the
                    ## seed is connected according to the PPI analysis
                    cleanedClusterTimeseriesFilename  = makeCleanedClusterTimeseriesFilename(subjectNumber, seed)
                    cat(sprintf("*** Reading cleaned cluster timeseries from: %s\n", cleanedClusterTimeseriesFilename))
                    cleanedClusterTimeseries = read.table(cleanedClusterTimeseriesFilename, header=FALSE)

                    functionallyConnectedClusterCount = dim(cleanedClusterTimeseries)[2]
                    
                    colnames(cleanedClusterTimeseries) = paste("Mean", seq.int(1, functionallyConnectedClusterCount), sep= "_") 

                    ## a 3D array to hold the largest and smallest
                    ## values for each subject, this is used to scale
                    ## the group level PPI graphs at the end of this
                    ## file
                    ## if ( is.null( listOfModels[[contrast]][[seed]][["maximaAndMinima"]] ) ) {
                        
                    ##     listOfModels[[contrast]][[seed]][["maximaAndMinima"]] = array(NA, dim=c(2, functionallyConnectedClusterCount, length(subjects)))
                    ##     dimnames(listOfModels[[contrast]][[seed]][["maximaAndMinima"]])[[1]] = c("min", "max")
                    ##     dimnames(listOfModels[[contrast]][[seed]][["maximaAndMinima"]])[[2]] = colnames(cleanedClusterTimeseries)
                    ##     dimnames(listOfModels[[contrast]][[seed]][["maximaAndMinima"]])[[3]] = subjects
                    ## }
                    
                    ## now we need to read the contrast regressor used in the PPI
                    contrastRegressorFilename = file.path(regressors.data.dir, paste(subjectNumber, contrast, "contrast.ppi.binary.1D", sep="."))
                    cat(sprintf("*** Reading %s contrast regressor from: %s\n", contrast, contrastRegressorFilename))
                    contrastRegressor = read.table( contrastRegressorFilename, header=FALSE)

                    stimuli = capwords(contrastsToStimuliMap[[contrast]])
                    contrastRegressorAsFactor = factor(contrastRegressor[, 1], levels=unique(contrastRegressor[ ,1]), labels=c("Other", stimuli))

                    stimulus.df = data.frame(
                        contrastRegressorAsFactor,
                        cleanedSeedTimeseries,
                        cleanedClusterTimeseries)
                    stimulus.df = subset(stimulus.df, contrastRegressorAsFactor != "Other")
                    ## drop the unused Other level 
                    stimulus.df$contrastRegressorAsFactor = drop.levels(stimulus.df$contrastRegressorAsFactor)
                    
                    ## now we need to demean the timeseries seperately
                    ## for each of the stimulus types

                    for (stimulus in levels(stimulus.df$contrastRegressorAsFactor) ) {
                        ## 1. the seed which is always the second column in the data frame
                        stimulus.df[stimulus.df$contrastRegressorAsFactor == stimulus, 2] =
                            scale(stimulus.df[stimulus.df$contrastRegressorAsFactor == stimulus, 2], center=TRUE, scale=FALSE)
                        ## 2. the cluster time series, which are always 3...end of the data frame
                        stimulus.df[stimulus.df$contrastRegressorAsFactor == stimulus, 3:dim(stimulus.df)[2]] =
                            scale(stimulus.df[stimulus.df$contrastRegressorAsFactor == stimulus, 3:dim(stimulus.df)[2]], center=TRUE, scale=FALSE)
                    } ## end of for (stimulus in levels(stimulus.df$contrastRegressorAsFactor) 

                    ## ## now figure out what the minima and maxima are
                    ## if ( is.null(dim(stimulus.df[, colnames(cleanedClusterTimeseries)])) ) {
                    ##     ## handle the case where there is only one
                    ##     ## cluster connected to the seed. This
                    ##     ## produces a vector when the single column is
                    ##     ## selected from stimulus.df. the apply
                    ##     ## function can't handle dimensionless objects
                    ##     listOfModels[[contrast]][[seed]][["maximaAndMinima"]][ , , subjectNumber] = 
                    ##         as.matrix(c(min=min(stimulus.df[, colnames(cleanedClusterTimeseries)]), max=max(stimulus.df[, colnames(cleanedClusterTimeseries)])), nrows=2)
                        
                    ## } else {
                    ##     listOfModels[[contrast]][[seed]][["maximaAndMinima"]][ , , subjectNumber] = 
                    ##         apply(stimulus.df[, colnames(cleanedClusterTimeseries)], 2, function(x) { c(min=min(x), max=max(x)) } )
                    ## }
                    
                    ## now melt the stimulus.df to facilitate graphing
                    melted.stimulus.df = melt(stimulus.df, id.vars=colnames(stimulus.df)[1:2], measure.vars=colnames(stimulus.df)[3:dim(stimulus.df)[2]], variable_name="cluster")

                    ## now fix the names of the names of the clusters so
                    ## they are not Mean_1 but the cluster name according
                    ## to the T&T atlas
                    melted.stimulus.df$cluster=factor( melted.stimulus.df$cluster,
                        levels=c(paste("Mean", seq.int(1, functionallyConnectedClusterCount), sep="_")),
                        labels=paste(seq.int(1, functionallyConnectedClusterCount), clusterWhereAmI))

                    ## the seed timeseries is always in column 2
                    x.axis=seed
                    y.axis="value"
                    group="contrastRegressorAsFactor"
                    xlabel=paste(substituteShortLabels(gsub(".", " ", x.axis, fixed=TRUE)), "Seed")
                    ylabel="Functionally Connected Cluster Activity"

                    subjectStudyGroup=getStudyGroupForSubject(demographics, subjectNumber)

                    cat(sprintf("Got %s for study group for %s subject\n", as.character(subjectStudyGroup), subjectNumber))
                    if ( ! subjectStudyGroup %in% c("NCL", "MDD") ) {
                        stop("Got an invalid group name (%s) for subject %s\n", subjectStudyGroup, subjectNumber)
                    }
                    ## graphTitle=sprintf("Psychophysiological Interaction\n%s vs. %s for %s Subject %s",
                    ##     stimuli[1], stimuli[2], ifelse(subjectStudyGroup == "NCL", "Control", subjectStudyGroup), subjectNumber)
                    graphTitle=sprintf("Psychophysiological Interaction\n%s vs. %s for %s Subject %s",
                        stimuli[1], stimuli[2], subjectStudyGroup, subjectNumber)


                    ## make linear models of the fit to each of the
                    ## stimulus types and save them so we can plot
                    ## them on one graph later
                    ## fixedFormula = as.formula(paste(y.axis, "~", x.axis, "+", group))
                    ## for  ( clst in levels(melted.stimulus.df$cluster) ) {
                    ##     mdl = lm(fixedFormula, data=subset(melted.stimulus.df, cluster==clst))
                    ##     listOfModels[[contrast]][[seed]][[subjectNumber]][[clst]] = mdl
                    ## }


                    graph=ggplot(melted.stimulus.df, aes_string(x=x.axis, y=y.axis, shape=group, color=group)) +
                        geom_point() + 
                            facet_wrap(~cluster, scales="free_y") +
                                labs(x=xlabel, y=ylabel) +
                                    ggtitle(graphTitle) +
                                        scale_color_brewer(name="Stimulus:", palette="Set1") +
                                            scale_shape_discrete(name="Stimulus:") +
                                                geom_smooth(method=lm, se=FALSE, fullrange=T) + 
                                                    my_theme + theme(legend.position="bottom")
                    
                    
                    ##for (fileType in c("pdf", "png")) {
                    for (fileType in c("png")) {                        
                        imageFilename=makeAllInOnePpiGraphFilename(subjectNumber, seed, inFileType=fileType)
                        if ( fileType == "png" ) {
                            listOfFiles[[contrast]][[seed]][[ length(listOfFiles[[contrast]][[seed]]) + 1 ]] <- imageFilename
                        }
                        cat("*** Saving all-in-one PPI graph to", imageFilename, "\n")
                        ggsave(imageFilename, graph)
                        ##print(graph)
                    } ## end of for (fileType in c("pdf, png")) {
                    
                    ## for ( level in levels(melted.stimulus.df$cluster) ) {

                    ##     graph=ggplot(melted.stimulus.df[melted.stimulus.df$cluster == level, ], aes_string(x=x.axis, y=y.axis, shape=group, color=group)) +
                    ##         geom_point() + 
                    ##             labs(x=xlabel, y=ylabel) +
                    ##                 ggtitle(paste(graphTitle, substituteShortLabels(level), sep="\n")) +
                    ##                     scale_color_brewer(name="Stimulus:", palette="Set1") +
                    ##                         scale_shape_discrete(name="Stimulus:") +
                    ##                             geom_smooth(method=lm, se=FALSE, fullrange=T) + 
                    ##                                 my_theme + theme(legend.position="bottom")
                    

                    ##     for (fileType in c("pdf", "png")) {
                    ##         imageFilename=makeIndividualClusterPpiGraphFilename(subjectNumber, seed, level, inFileType=fileType)
                    ##         listOfFiles[[contrast]][[seed]][[ length(listOfFiles[[contrast]][[seed]]) + 1 ]] <- imageFilename
                    ##         cat("*** Saving individual PPI graph to", imageFilename, "\n")
                    ##         ggsave(imageFilename, graph)
                    ##     } ## end of for (fileType in c("pdf, png")) {   

                    ## } ## end of for ( level in levels(melted.stimulus.df$cluster) ) {

                } else {
                    cat(sprintf("*** No such file: %s\n*** Skipping seed because no clusters were connected to the seed\n", clusterLocationsFilename))
                }

            } ## end of if (file.exists (cleanedSeedTimeseriesFilename)) {
            else {
                cat(sprintf("*** No such file: %s\n", cleanedSeedTimeseriesFilename))
            }
        } ## end of for (subjectNumber in subjectOrder$subject) { 

        seedCount = seedCount + 1
    } ## end of for ( seed in contrastsToSuffices[[contrast]][["wholeBrainRois"]] ) {
} ## end of for ( contrast in contrasts) { 

makeHtmlFiles(listOfFiles)

## ## numberOfNewTimepoints = 50

## for (contrast in names(listOfModels) ) {
##     seedCount = 1
##     for ( seed in names(listOfModels[[contrast]]) ) {

##         ## if you want to figure out how it works try the commented following at the end of this file
##         extrema=matrix(c(
##             apply(t(listOfModels[[contrast]][[seed]][["maximaAndMinima"]][1, , ]), 2, min),
##             apply(t(listOfModels[[contrast]][[seed]][["maximaAndMinima"]][2, , ]), 2, max)),
##             nrow=2, byrow=TRUE)
##         colnames(extrema)= paste("Mean", seq.int(1, dim(extrema)[2]), sep="_")
##         rownames(extrema) = c("min", "max")

##         subjects = names(listOfModels[[contrast]][[seed]])[grep ("[0-9]{3}_[A]", names(listOfModels[[contrast]][[seed]]))]
##         for (subject in subjects ) {

##             ## this code computes the min and max across the columns (second
##             ## dimension, cleanedClusterTimeseries) and subjects (3rd dimension)
##             ## read the seed timeseries
##             cleanedSeedTimeseries = makeCleanedSeedTimeseriesFilename(subjectNumber, seed)
##             if (file.exists (cleanedSeedTimeseries)) {
##                 cat(sprintf("*** Reading cleaned seed timeseries from: %s\n", cleanedSeedTimeseries))
##                 cleanedSeedTimeseries = read.table(cleanedSeedTimeseries, header=FALSE)
##                 colnames(cleanedSeedTimeseries) = seed.clusterWhereAmI[seedCount]


##                 ## now read the locations (as human readable text names) of the clusters connected to the seed as determined by the PPI analysis
##                 clusterLocationsFilename=file.path(Group.results.ppi.dir, sprintf("clusterLocations.fwhm%0.1f.%s.%s.%s.%s.csv", usedFwhm, task, groups, makeRoiInteractionSuffix(seed), fLabel))
##                 if ( file.exists(clusterLocationsFilename)) {
##                     cat("*** Reading cluster locations from", clusterLocationsFilename, "\n")
##                     ## the gsub here chews up multiple consequtive spaces and replaces them with a single space
##                     clusterWhereAmI=gsub(" +", " ", scan(file=clusterLocationsFilename, what='character', sep=','))

                    
##                     cl = "Mean_1"
##                     clIndex = as.integer(strsplit(cl, "_", fixed=TRUE)[[1]][2])
##                     newGrid = expand.grid(
##                         seq(min(cleanedSeedTimeseries), max(cleanedSeedTimeseries), numberOfNewTimepoints),
##                         seq(extrema["min", cl], extrema["max", cl], numberOfNewTimepoints),
##                         capwords(contrastsToStimuliMap[[contrast]])
##                         )
##                     colnames(newGrid) = c(
##                     newGrid$subject = rep(subject, dim(newGrid)[1])
##                     newGrid$group   = rep(getStudyGroupForSubject(demographics, subjectNumber), dim(newGrid)[1])
                    
##                     newGrid$ts = stats::predict(listOfModels[[contrast]][[seed]][[subject]][[clIndex]], newdata=newGrid)

##                     if (exists("oldGrid") ) {
##                         oldGrid = rbind(oldGrid, newGrid)
##                     } else {
##                         oldGrid = newGrid
##                     }
##                 } ## end of if ( file.exists(clusterLocationsFilename)) {
##             } ## end of if (file.exists (cleanedSeedTimeseries)) {
##         } ## end of for (subject in listOfModels[[contrast]][[seed]] ) {
##     } ## end of for ( seed in names(listOfModels[[contrast]]) ) {
## } ## end of for (contrast in names(listOfModels) ) {


## ##################################################################################################
## code to figure out how the computation of the min and max work
## > a[, , 1]=1:6
## > a[, , 2]=(1:6) + 6
## > a[, , 3]=(1:6) -24
## > a
## , , 1

##      [,1] [,2] [,3]
## [1,]    1    3    5
## [2,]    2    4    6

## , , 2

##      [,1] [,2] [,3]
## [1,]    7    9   11
## [2,]    8   10   12

## , , 3

##      [,1] [,2] [,3]
## [1,]  -23  -21  -19
## [2,]  -22  -20  -18

## > a[, 1, ]
##      [,1] [,2] [,3]
## [1,]    1    7  -23
## [2,]    2    8  -22
## > a[1, , ]
##      [,1] [,2] [,3]
## [1,]    1    7  -23
## [2,]    3    9  -21
## [3,]    5   11  -19
## > t(a[1, , ])
##      [,1] [,2] [,3]
## [1,]    1    3    5
## [2,]    7    9   11
## [3,]  -23  -21  -19
## > apply(t(a[1, , ]), 2, min)
## [1] -23 -21 -19
## > a
## , , 1

##      [,1] [,2] [,3]
## [1,]    1    3    5
## [2,]    2    4    6

## , , 2

##      [,1] [,2] [,3]
## [1,]    7    9   11
## [2,]    8   10   12

## , , 3

##      [,1] [,2] [,3]
## [1,]  -23  -21  -19
## [2,]  -22  -20  -18

## > apply(t(a[1, , ]), 2, max)
## [1]  7  9 11
## > matrix(c(apply(t(a[1, , ]), 2, min), apply(t(a[1, , ]), 2, max)), nrow=2, byrow=TRUE)
##      [,1] [,2] [,3]
## [1,]  -23  -21  -19
## [2,]    7    9   11
## ##################################################################################################

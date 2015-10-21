rm(list=ls())
graphics.off()

library(gdata)
library(reshape)
library(ggplot2)
library(robustbase)
library(MASS)
library(boot)

source("scoreMasc.r")

ncpus=1
if ( Sys.info()["sysname"] == "Darwin" ) {
    ncpus=as.integer(strsplit(system("sysctl hw.ncpu", intern=T), ' ')[[1]][2])
    cat(paste("Found" , ncpus, ifelse(ncpus == 1, "cpu", "cpus"), "\n"))
} else if ( Sys.info()["sysname"] == "Linux" ) {
    ncpus=as.integer(system("cat /proc/cpuinfo |grep -c processor", intern=T))
    cat(paste("Found" , ncpus, ifelse(ncpus == 1, "cpu", "cpus"), "\n"))
} else {
    cat(paste("Sorry can't determine the number of CPUs in this architecture defaulting to", ncpus, "\n"))
}

##source('pcor.R')

## colorblind friendly paletts 
                                        # The palette with grey:
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

                                        # The palette with black:
cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


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
    returnSubstitutedLabel = gsub("[0-9]+ ", "", gsub("Inf", "Inferior", gsub("Sup", "Superior", gsub("Gy", "Gyrus", gsub("R", "Right", gsub("L", "Left", inLevel, fixed=TRUE), fixed=TRUE), fixed=TRUE), fixed=TRUE), fixed=TRUE))
    
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

makeGraphTitle <- function(inLevel) {

    graphTitle= gsub("^[0-9]+[ ]*", "", gsub("Inf", "Inferior", gsub("Gy", "Gyrus", gsub("Sup", "Superior", inLevel))))

    return(graphTitle)
}

clust.header = c("Volume", "CM RL", "CM AP", "CM IS", "minRL",
    "maxRL", "minAP", "maxAP", "minIS", "maxIS", "Mean", "SEM", "Max Int",
    "MI RL", "MI AP", "MI IS")

if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    cat(paste("Sorry can't set data directories for this computer\n"))
}

Admin.data.dir=normalizePath(file.path(root.dir, "sanDiego/cPine/data/admin"))
Group.data.dir=normalizePath(file.path(root.dir, "sanDiego/cPine/data/Group.data"))
Group.results.dir=normalizePath(file.path(root.dir, "sanDiego/cPine/data/Group.results"))

####################################################################################################
## These variables control what barcharts are created

createIndividualRegressionImageFiles=FALSE
createSingleRegressionPdfFile=TRUE

runRegressions=TRUE

####################################################################################################
groups="mddAndCtrl"
fLabel="group.F-value"
usedFwhm=4.2
task="pine"

## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
##demographicsFilename=file.path(Group.data.dir, "demographics2.csv")


demographicsFilename=file.path(Admin.data.dir, "data_entry_current_021113.csv")

cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "#N/A", "#VALUE", "#VALUE!", "n/a"))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))

ksads.filename=file.path(Admin.data.dir, "KSADS.age.of.onset.pine.csv")
cat("Reading KSADS data from: ", ksads.filename, "\n")
ksads=read.csv(ksads.filename, header=TRUE)
ksads$Subject.Number=gsub("A", "", ksads$Subject.Number, fixed=TRUE)
ksads[ksads$Subject.Number=="300", "Subject.Number"] = "169/300"
cat(sprintf("Read %d subjects from the KSADS file\n", dim(ksads)[1]))

## these are the ACC ROIs
##seeds=c(paste("acc", "i", 1:9, "L", sep=""), paste("acc", "i", 1:9, "R", sep=""), paste("acc", "s", 1:7, "L", sep=""), paste("acc", "s", 1:7, "R", sep=""))
seeds=c(paste("acc", "i", 8:9, "L", sep=""), paste("acc", "i", 8:9, "R", sep=""))
##, paste("acc", "s", 1:7, "L", sep=""), paste("acc", "s", 1:7, "R", sep="")) 

regressionVariables=list(
    ##list(variable="CGI.CGAS",       name="Children's Global Assessment Scale"),
    ##list(variable="PSWQ",           name="Penn State Worry Questionnaire"),
    ##list(variable="CoRum",          name="Corumination Questionnaire")#,
    ##list(variable="RSQ.fixed",      name="Ruminative Responses Questionnaire"),
    ##list(variable="CoRumII",        name="Corumination Questionnaire II"),
    list(variable="CDRS.tscore",    name="Children's Depression Rating Scale (Standardized)"),
    ## list(variable="RADS.DM.Tscore", name="Reynolds Adolescent Depression Scale Dysphoric Mood (Standardized)"),
    ## list(variable="RADS.AN.Tscore", name="Reynolds Adolescent Depression Scale Anhedonia/Negative Affect (Standardized)"),
    ## list(variable="RADS.NS.Tscore", name="Reynolds Adolescent Depression Scale Negative Self-evaluation (Standardized)"),
    ## list(variable="RADS.SC.Tscore", name="Reynolds Adolescent Depression Scale Somatic Complaints (Standardized)"),
    ## list(variable="MADRS.raw",      name="Montgomery-Asberg Depression Scale"),
    list(variable="BDI.II",         name="Beck Depression Inventory II"),
    ##list(variable="CDI",            name="Children's Depression Inventory"),
    ##list(variable="MASC.tscore",    name="Multidimensional Anxiety Scale for Children (Standardized)")

    list(variable="age.of.onset", name="KSADS: Age of onset"),
    list(variable="number.of.episodes", name="KSADS: Number of episodes")
    
    )
## ## extract the list of variable names from the regressionVariables
## list rvs=regression variables select DOB, MRI, and MASC.total
## because we will need these later when we are calculating the MASC
## scores
rvs=c("DOB", "MRI", "MASC.total", unlist(regressionVariables)[grep ("variable", names(unlist(regressionVariables)))])
## ## select only the columns we want to perform regressions on from the hnrcData data frame
##m=match(unlist(regressionVariables)[grep ("variable", names(unlist(regressionVariables)))], colnames(demographics))
## ## remove NAs caused by the BDI and POMS not coming from the hnrcData frame
##m=m[!is.na(m)]
## m now contains only the column numbers of the neuropsych variables we'll regress the %cs against

if (runRegressions) {
    regressionsDirectory=file.path(Group.results.dir) #, "regressions")
    if ( ! file.exists(regressionsDirectory) ) {
        dir.create(regressionsDirectory, recursive=TRUE)
    }
    regressionsFilename= file.path(regressionsDirectory, paste("regressions.contrast", "txt", sep="."))
    if (file.exists(regressionsFilename) ) {
        file.remove(regressionsFilename)
    }
    coRumFilename= file.path(regressionsDirectory, paste("coRum.contrast", "txt", sep="."))
    if (file.exists(coRumFilename) ) {
        file.remove(coRumFilename)
    }
}

##contrasts=c("fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad", "allEmotiveVsNeutral", "happyRemVsHappyNotrem", "fearfulRemVsFearfulNotrem", "neutralRemVsNeutralNotrem", "sadRemVsSadNotrem", "allRemVsAllNotrem")
##contrasts=c("fearfulVsHappyExtractedRightSgAcc") ##,
contrasts=c("fearfulVsHappy", "fearfulVsNeutral", "fearfulVsSad", "happyVsNeutral", "happyVsSad", "neutralVsSad") ##, "allEmotiveVsNeutral") ##, "happyRemVsHappyNotrem", "fearfulRemVsFearfulNotrem", "neutralRemVsNeutralNotrem", "sadRemVsSadNotrem", "allRemVsAllNotrem")

mystack <- stack()
header="Contrast,gender,cluster,RL,AP,IS,RegressionVariable,S,DoF,pValue,Rho,Significance"

## this should be set to the number of subjects you have in the bucket
## files and hence the subjectOrder file.
expectedNumberOfSubjects=67

##allPvalues=vector("numeric", length=0)
seperateGenders=FALSE
useMask=FALSE
if ( useMask ) {
    mask="apriori.rois.withVmpfc.3mm"
}

for (contrast in contrasts ) {

    cat("####################################################################################################\n")
    cat(sprintf("*** Graphing ROIs for the %s contrast\n", contrast))
    ##roiStats.fwhm4.2.pine.ctrlOnly.neutralVsSad.group.F-value.txt

    if (useMask) {
        roistats.filename=file.path(Group.results.dir, sprintf("roiStats.fwhm%0.1f.%s.%s.%s.%s.%s.txt", usedFwhm, task, groups, contrast, fLabel, mask))
    } else {
        roistats.filename=file.path(Group.results.dir, sprintf("roiStats.fwhm%0.1f.%s.%s.%s.%s.txt", usedFwhm, task, groups, contrast, fLabel))
    }
    
    if(file.exists(roistats.filename)) {
        roistats=read.table(roistats.filename, header=T, sep="")
        if (dim(roistats)[1] != expectedNumberOfSubjects) {
            stop(sprintf("Check the roistats variable. It has dim[%d, %d] but am expecting %d subjects\n", dim(roistats)[1], dim(roistats)[2], expectedNumberOfSubjects))
        }
        ## dump the first column as it's only the file name
        ##roistats=roistats[, -1]
        clusterCount=length(grep("Mean", colnames(roistats)))
        if (clusterCount > 0 ) {
            cat(sprintf("*** %d ** clusters found in %s\n", clusterCount, roistats.filename))
            
### Most of the following code up the the first long row of # is bookkeeping to get the data frame in order

            if (useMask) {
                clustersFilename=sprintf("clust.fwhm%0.1f.%s.%s.%s.%s.%s.txt", usedFwhm, task, groups, contrast, fLabel, mask)
            } else { 
                clustersFilename=sprintf("clust.fwhm%0.1f.%s.%s.%s.%s.txt", usedFwhm, task, groups, contrast, fLabel)
            }
            cat("*** Reading", file.path(Group.results.dir, clustersFilename), "\n")
            clusters=read.table(file.path(Group.results.dir, clustersFilename))
            colnames(clusters) = clust.header
            ## this command contains the locations, as text, of the clusters and is the output of a perl script

            if (useMask) {
                clusterLocationsFilename=file.path(Group.results.dir, sprintf("clusterLocations.fwhm%0.1f.%s.%s.%s.%s.%s.csv", usedFwhm, task, groups, contrast, fLabel, mask))
            } else {
                clusterLocationsFilename=file.path(Group.results.dir, sprintf("clusterLocations.fwhm%0.1f.%s.%s.%s.%s.csv", usedFwhm, task, groups, contrast, fLabel))
            }
            cat("*** Reading cluster locations from", clusterLocationsFilename, "\n")
            ## the gsub here chews up multiple consequtive spaces and replaces them with a single space
            clusterWhereAmI=gsub(" +", " ", scan(file=clusterLocationsFilename, what='character', sep=','))
            
            ## this file stores the order of the subjects in each of the following BRIK files
            subjectOrderFilename=file.path(Group.data.dir, paste("subjectOrder", groups, contrast, "REML.csv", sep="."))
            cat("*** Reading", subjectOrderFilename, "\n")
            subjectOrder=read.csv(subjectOrderFilename, header=T)
            subjectOrder$subject=as.vector(subjectOrder$subject)
            ## due to some sort of screw up 169/300_A is the same subject so fix it here.
            badSubjectIds=which (subjectOrder$subject=="300_A")
            if (length(badSubjectIds) > 0 ) {
                subjectOrder[badSubjectIds, "subject"]="169/300_A"
            }
            subjectOrder$subject=as.factor(gsub("_A", "", as.character(subjectOrder$subject), fixed=TRUE))
            cat(sprintf("*** The subjectOrder data frame has %s unique subjects\n",  length(unique(subjectOrder$subject))))
            
            mgd=cbind(subjectOrder,
                demographics[match(subjectOrder$subject, demographics$ID), c("Grp", "Gender", rvs[!grepl("age.of.onset|number.of.episodes", rvs)] )],
                ksads[match(subjectOrder$subject, ksads$Subject.Number), c("Subject.Number", "age.of.onset", "number.of.episodes")]
                )
            ##colnames(mgd)=c("subject", "Group")
            ## drop the unused levels from the Group factor
            mgd$Grp=drop.levels(mgd$Grp) 
            mgd=cbind(roistats, mgd)

            ## this complicated looking regexp stuff cleans up the years in DOB
            ## and MRI with 4 digits to be just 2 digits
            mgd$DOB=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", mgd$DOB)
            mgd$MRI=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", mgd$MRI)
            mgd$DOB=as.Date(mgd$DOB, "%m/%d/%y")
            mgd$MRI=as.Date(mgd$MRI, "%m/%d/%y")
            mgd$age.in.days=difftime(mgd$MRI, mgd$DOB, units="days")
            mgd$age.in.days=as.numeric(mgd$age.in.days)
            
            mgd$age.in.weeks=difftime(mgd$MRI, mgd$DOB, units="weeks")
            mgd$age.in.weeks=as.numeric(mgd$age.in.weeks)
            mgd$age.in.years=(mgd$age.in.weeks)/52

            if ("MASC.tscore" %in% colnames(mgd) ) {
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
            ## now drop the normal controls
            cat("Dropping all subjects but MDD from the mgd data frame\n")
            mgd=mgd[mgd$Grp=="MDD", ]
            mgd[, c("File", "Sub.brick")]=list(NULL)
            
            melted.roistats=melt(mgd, id.vars=c("subject", "Grp", "Gender", rvs), measure.vars=paste("Mean_", seq(1, clusterCount), sep=""), variable_name="cluster")
            
            melted.roistats$cluster=factor(melted.roistats$cluster,
                levels=c(paste("Mean_", seq(1, clusterCount), sep="")),
                labels=paste(1:length(clusterWhereAmI), clusterWhereAmI))

            ##stop("Check mgd data frame\n")

            cat("*** Now running regression models. Check", regressionsFilename, "for the output\n")
            
            ##regressionsCsvFilename= file.path(regressionsDirectory, paste("regressions", model, "csv", sep="."))
            ##stop()]

            sink(regressionsFilename, append=TRUE)
            for ( j in 1:length(regressionVariables ) ) {
                ## clCounter=clusterCounter
                clCounter=1
                for ( level in levels(melted.roistats$cluster) ) {
                    cat("####################################################################################################\n")
                    cat("*** Contrast: ", contrast, "\n")
                    cat("*** Running regression model for ", level, "on", regressionVariables[[j]]$name, "\n")
                    
                    ##csvLine=paste(level, group, phase, glt, regressionVariables[[i]]$name, sep=",")
                    regressionFormula=as.formula(paste("value",  "~", regressionVariables[[j]]$variable, sep=" "))
                    cat("*** Regression formula: ")
                    print(regressionFormula)
                    mdl.df.all=   melted.roistats[melted.roistats$cluster %in% level, ]

                    if (seperateGenders) {
                        mdl.df.male=  melted.roistats[melted.roistats$Gender=="M"   & melted.roistats$cluster %in% level, ]
                        mdl.df.female=melted.roistats[melted.roistats$Gender=="F" & melted.roistats$cluster %in% level, ]
                    }
                    
                    ct.all=cor.test(mdl.df.all[, regressionVariables[[j]]$variable], mdl.df.all$value, method="spearman")
                    cat("### BOTH GENDERS COMBINED\n")
                    print(ct.all)
                    if (seperateGenders) {
                        ct.male=cor.test(mdl.df.male[, regressionVariables[[j]]$variable], mdl.df.male$value, method="spearman")
                        cat("### MALE ONLY\n")
                        print(ct.male)

                        ct.female=cor.test(mdl.df.female[, regressionVariables[[j]]$variable], mdl.df.female$value, method="spearman")
                        cat("### FEMALE ONLY\n")
                        print(ct.female)
                    }
                    mdl.df.rlm=mdl.df.all[, c("value", regressionVariables[[j]]$variable, "Gender")]
                    regression.formula=as.formula(paste("value", "~", regressionVariables[[j]]$variable, sep=" "))
                    lm.model.nogender=lm(regression.formula, data=mdl.df.rlm)
                    print(summary(lm.model.nogender))

                    ##allPvalues = c(allPvalues, ct.all$p.value)

                    regression.formula=as.formula(paste("value", "~", regressionVariables[[j]]$variable, "+", "Gender", sep=" "))
                    lm.model.withgender=lm(regression.formula, data=mdl.df.rlm)
                    print(summary(lm.model.withgender))
                    
                    ##csvLine=sprintf("%s,%s,%s,%.2f,%.2f,%.5f,%.2f,%s",
                    ##  contrast, makeGraphTitle(level), regressionVariables[[j]]$name, ct$statistic, ifelse(is.null(ct$parameter), NA, ct$parameter),
                    ##  ct$p.value, ct$estimate, make.significance.indications(ct$p.value))
                    csvLine=sprintf("%s,%s,%s,%s,%s, %.2f,%.2f,%.5f,%.2f,%s",
                        contrast, "all", makeGraphTitle(level), paste(clusters[clCounter, c("CM RL", "CM AP", "CM IS")], collapse=","), regressionVariables[[j]]$name, ct.all$statistic,
                        ifelse(is.null(ct.all$parameter), NA, ct.all$parameter),
                        ct.all$p.value, ct.all$estimate, make.significance.indications(ct.all$p.value))
                    cat(csvLine, "\n")
                    push(mystack, csvLine)

                    if (seperateGenders) {
                        csvLine=sprintf("%s,%s,%s,%s,%s, %.2f,%.2f,%.5f,%.2f,%s",
                            contrast, "male", makeGraphTitle(level), paste(clusters[clCounter, c("CM RL", "CM AP", "CM IS")], collapse=","), regressionVariables[[j]]$name, ct.male$statistic,
                            ifelse(is.null(ct.male$parameter), NA, ct.male$parameter),
                            ct.male$p.value, ct.male$estimate, make.significance.indications(ct.male$p.value))
                        cat(csvLine, "\n")
                        push(mystack, csvLine)
                        
                        csvLine=sprintf("%s,%s,%s,%s,%s, %.2f,%.2f,%.5f,%.2f,%s",
                            contrast, "female", makeGraphTitle(level), paste(clusters[clCounter, c("CM RL", "CM AP", "CM IS")], collapse=","), regressionVariables[[j]]$name, ct.female$statistic,
                            ifelse(is.null(ct.female$parameter), NA, ct.female$parameter),
                            ct.female$p.value, ct.female$estimate, make.significance.indications(ct.female$p.value))
                        cat(csvLine, "\n")
                        push(mystack, csvLine)
                    }  
                    ##stop()
                    
                    ## if (ct.all$p.value < 0.05 ) {
                    ##   ## make scatter plots for the ROIs that are significant
                    ##   imageDirectory=file.path(Group.results.dir, contrast)
                    ##   if ( ! file.exists(imageDirectory) ) {
                    ##     dir.create(imageDirectory)
                    ##   }
                    ##   imageFilename=file.path(imageDirectory, sprintf("%s.fwhm%0.1f.%s.%s.and.%s.pdf", gsub(" +", ".", level),  usedFwhm, task, contrast, regressionVariables[[j]]$variable))
                    ##   message(paste("*** Creating", imageFilename, "\n"))
                    
                    ##   ylabel="RSFC"
                    ##   graphTitle=makeGraphTitle(level)
                    ##   ##sprintf("%s", toupper(contrast))
                    ##   my.base.size=18

                    ##   graph=ggplot(mdl.df.all, aes(x=eval(parse(text=regressionVariables[[j]]$variable)), y=value)) +
                    ##     theme_bw(base_size =  my.base.size) +
                    ##       geom_point() +
                    ##         geom_smooth(method="rlm", color="black") +
                    ##         scale_fill_brewer(palette="Set1") +
                    ##           labs(title = graphTitle, y=ylabel, x=regressionVariables[[j]]$name) +
                    ##             theme(legend.position="none",
                    ##                   ##panel.grid.major = element_blank(),
                    ##                   panel.grid.minor = element_blank(),
                    ##                   axis.title.x = element_text(size=my.base.size, vjust=0),
                    ##                   axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
                    ##                   plot.title=element_text(size=my.base.size*1.2, vjust=1))
                    ##   ggsave(imageFilename)
                    ## } ## end of if (ct.all$p.value < 0.05 )
                    
                    
                    clCounter=clCounter+1
                    cat("clCounter is now", clCounter, "\n")
                } ## end of for ( level in levels(melted.roistats$cluster) )

                ##stop("Stopping")
            } ## end of for ( level in levels(roistats.summary$cluster) )
            
            sink()
            ##stop()
            
        } ## end of if (clusterCount > 0 )
    } else {
        cat("### ", roistats.filename, "does not exist. Skipping\n")
    } ## end of if(file.exists(roistats.filename))
} ## end of for (contrast in contrasts)

regressionsCsvFilename=file.path(regressionsDirectory, paste("regressions.contrast", "csv", sep="."))
if (file.exists(regressionsCsvFilename) ) {
    file.remove(regressionsCsvFilename)
}

cat("*** Check", regressionsCsvFilename, "for the CSV output\n")
l=mystack$value()
cat(header, "\n", file=regressionsCsvFilename)
for (i in 1:length(l)) {
    cat (l[[i]], "\n", file=regressionsCsvFilename, append=TRUE)
}

graphics.off()

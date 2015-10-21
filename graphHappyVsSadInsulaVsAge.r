rm(list=ls())
graphics.off()

library(reshape)
library(ggplot2)

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

fixSubjectOrderTable <- function (inSubjectOrderTable) {
    inSubjectOrderTable$subject=gsub("_A", "", as.character(inSubjectOrderTable$subject), fixed=TRUE)
    inSubjectOrderTable[inSubjectOrderTable$subject=="300", "subject"] = "169/300"
    inSubjectOrderTable$subject=as.factor(inSubjectOrderTable$subject)

    return(inSubjectOrderTable)
}

fixDates <- function (inRoistats) {
    ## this complicated looking regexp stuff cleans up the years in DOB
    ## and MRI with 4 digits to be just 2 digits
    ## month day year
    inRoistats$DOB=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", inRoistats$DOB)
    inRoistats$MRI=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", inRoistats$MRI)
    
    ## now convert to year/month/day
    inRoistats$DOB=sub("([0-9]{1,2})/([0-9]{1,2})/([0-9]{2})", "\\3/\\1/\\2", inRoistats$DOB)
    inRoistats$MRI=sub("([0-9]{1,2})/([0-9]{1,2})/([0-9]{2})", "\\3/\\1/\\2", inRoistats$MRI)
    
    inRoistats$DOB=as.Date(inRoistats$DOB, "%y/%m/%d")
    inRoistats$MRI=as.Date(inRoistats$MRI, "%y/%m/%d")

    return(inRoistats)
}

computeAge <- function(inRoistats) {
    
    age.in.weeks=difftime(inRoistats$MRI, inRoistats$DOB, units="weeks")
    age.in.weeks=as.numeric(age.in.weeks)

    inRoistats$age.in.years=age.in.weeks/52

    return(inRoistats)
}

makeGraphTitle <-function(inRoi) {

    return(sub(".VS", " Ventral Striatum", inRoi, fixed=TRUE))

}

make.significance.indications <- function(pValues, which.pValues=c(1)) {

  Signif=symnum(pValues, corr=FALSE, na=FALSE, cutpoints = c(0,  .001,.01, .05, .1, 1),
    symbols   =  c("***", "**", "*", ".", " "))
  f=format(Signif)

  ## only return the first one as we're only interested in marking significant group effects
  return(f[which.pValues])
}

####################################################################################################
### Code
####################################################################################################

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
Group.results.ppi.dir=normalizePath(file.path(root.dir, "sanDiego/cPine/data/Group.results/ppi"))
conditions=c("happy", "sad")

## These filenames are from the PPI analysis looking at the
## connectivity of the cluster that contains the AIC identified from
## the happy vs sad contrast on the whole-brain task-related results

subjectOrder.ppi.filename=file.path(Group.data.dir, "subjectOrder.mddAndCtrl.roi4.seed.happyVsSad.ROIinteraction.REML.z-score.csv")
roistats.ppi.filename=file.path(Group.results.ppi.dir, "roiStats.fwhm4.2.pine.mddAndCtrl.roi4.seed.happyVsSad.ROIinteraction.group.F-value.txt")
cluster.locations.ppi.filename=file.path(Group.results.ppi.dir, "clusterLocations.fwhm4.2.pine.mddAndCtrl.roi4.seed.happyVsSad.ROIinteraction.group.F-value.csv")

### END OF PPI FILENAMES

## These filenames are from the PPI analysis looking at the
## connectivity of the cluster that contains the AIC identified from
## the happy vs sad contrast on the whole-brain task-related results

subjectOrderFiles=c("subjectOrder.mddAndCtrl.happy.REML.csv", "subjectOrder.mddAndCtrl.sad.REML.csv")
subjectOrderFiles=sapply(subjectOrderFiles, function(x) { file.path(Group.data.dir, x) })
names(subjectOrderFiles)=NULL

## these roistats are from the between-group whole brain task results
roistats.files=c("roiStats.pctChange.happy.stimulus.fwhm4.2.pine.mddAndCtrl.happyVsSad.group.F-value.contrast.mddAndCtrl.txt", "roiStats.pctChange.sad.stimulus.fwhm4.2.pine.mddAndCtrl.happyVsSad.group.F-value.contrast.mddAndCtrl.txt")
roistats.files=sapply(roistats.files, function(x) { file.path(Group.results.dir, x) })
names(roistats.files)=NULL

cluster.locations.filename=file.path(Group.results.dir, "clusterLocations.fwhm4.2.pine.mddAndCtrl.happyVsSad.group.F-value.csv")

### END OF WHOLE-BRAIN TASK-RELATED FILENAMES

demographicsFilename=file.path(Admin.data.dir, "0-data_entry_current_10152013.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT"))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))
demographics$ID=as.factor(gsub("_[ABCD]", "", as.character(demographics$ID), fixed=FALSE))

## this file stores all of the SLES scores
slesFilename=file.path(Admin.data.dir, "SLES_20140820.csv")
cat("*** Reading", slesFilename, "\n")
sles=read.csv(slesFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT", "."))
cat(sprintf("*** Read SLES data for %s unique subjects\n",  length(unique(sles$subject))))
sles=sles[, seq(1, grep ("X.80_comment", colnames(sles)))]

ksads.filename=file.path(Admin.data.dir, "KSADS.age.of.onset.pine.csv")
cat("Reading KSADS data from: ", ksads.filename, "\n")
ksads=read.csv(ksads.filename, header=TRUE)
ksads$Subject.Number=gsub("A", "", ksads$Subject.Number, fixed=TRUE)
ksads[ksads$Subject.Number=="300", "Subject.Number"] = "169/300"
cat(sprintf("Read %d subjects from the KSADS file\n", dim(ksads)[1]))

subjectOrder=do.call(rbind, lapply(subjectOrderFiles,
    function(X) {
        data.frame(read.table(X, header=T, sep=""))
    }
    ))
subjectOrder=fixSubjectOrderTable(subjectOrder)

roistats=do.call(rbind, lapply(roistats.files,
    function(X) {
        data.frame(read.table(X, header=T, sep=""))
    }
    ))
roistats=cbind(subjectOrder, roistats)
cluster.locations=gsub(" +", " ", scan(file=cluster.locations.filename, what='character', sep=','))

## Now pick only the cluster with the insula componant
roistats=roistats[, c("subject", "Mean_4")]
roistats$stimulus=as.factor(rep(c("happy", "sad"), each=length(unique(roistats$subject))))
##roistats=roistats[, -which(colnames(roistats) %in% c("File", "Sub.brick"))]

roistats=cbind(roistats,
    demographics[match(roistats$subject, demographics$ID), c("Grp", "Gender", "DOB", "MRI", "BDI.II")]
    )
roistats=computeAge(fixDates(roistats))

## add this column with the subject ID concatenated with the timepoint
## so we can pull the correct data from the SLES data frame. This is
## done because the SLES data frame contains data for all time points.
roistats$subjectWithTimePoint=as.factor(paste(roistats$subject, "A", sep="_"))
roistats=cbind(
    roistats,
    sles [match(roistats$subjectWithTimePoint, sles$subject),
          c("num.stressful.events", "num.severe.events", "total.sum.stress", "num.1.little.to.no.effect", "num.2.some.effect", "num.3.moderate.effect", "num.4.great.effect")],
    ksads[match(roistats$subject, ksads$Subject.Number),
          c("number.of.episodes", "age.of.onset")]
    )
## remove the now unnecessary subject with timepoint column
roistats$subjectWithTimePoint=NULL

roistats$num.stressful.events.sqrt=sqrt(roistats$num.stressful.events)
roistats$num.severe.events.sqrt   =sqrt(roistats$num.severe.events)
roistats$total.sum.stress.sqrt    =sqrt(roistats$total.sum.stress)

melted.roistats=melt(roistats, id.vars=c("subject", "stimulus",
                                   "Grp", "Gender", "DOB", "MRI", "age.in.years", ## demographic data
                                   "BDI.II",
                                   "num.stressful.events", "num.severe.events", "total.sum.stress", ## SLES data 
                                   "number.of.episodes", "age.of.onset" ## KSADS data
                                   ),
    measure.vars=c("Mean_4"), variable_name="cluster")

melted.roistats$cluster=factor(melted.roistats$cluster,
    levels="Mean_4",
    labels=cluster.locations[4])


## now add all of the demographics and psychiatric variables ot the
## PPI-derived roistats df

subjectOrder.ppi=read.csv(subjectOrder.ppi.filename, header=TRUE)
subjectOrder.ppi=fixSubjectOrderTable(subjectOrder.ppi)

roistats.ppi=read.table(roistats.ppi.filename, header=TRUE)
roistats.ppi=cbind(subjectOrder.ppi, roistats.ppi)
cluster.locations.ppi=gsub(" +", " ", scan(file=cluster.locations.ppi.filename, what='character', sep=','))

ppi.cluster.numbers=c(1, 3, 8)
ppi.clusters.colnames=paste("Mean", ppi.cluster.numbers, sep="_")

## Now pick only the insula and Fusiform gyrus components only
roistats.ppi=roistats.ppi[, c("subject", ppi.clusters.colnames)]

roistats.ppi=cbind(
    roistats.ppi,
    demographics[match(roistats.ppi$subject, demographics$ID), c("Grp", "Gender", "DOB", "MRI", "BDI.II")]
    )
roistats.ppi=computeAge(fixDates(roistats.ppi))

## add this column with the subject ID concatenated with the timepoint
## so we can pull the correct data from the SLES data frame. This is
## done because the SLES data frame contains data for all time points.
roistats.ppi$subjectWithTimePoint=as.factor(paste(roistats.ppi$subject, "A", sep="_"))
roistats.ppi=cbind(
    roistats.ppi,
    sles [match(roistats.ppi$subjectWithTimePoint, sles$subject),
          c("num.stressful.events", "num.severe.events", "total.sum.stress", "num.1.little.to.no.effect", "num.2.some.effect", "num.3.moderate.effect", "num.4.great.effect")],
    ksads[match(roistats.ppi$subject, ksads$Subject.Number),
          c("number.of.episodes", "age.of.onset")]
    )
## remove the now unnecessary subject with timepoint column
roistats.ppi$subjectWithTimePoint=NULL

roistats.ppi$num.stressful.events.sqrt=sqrt(roistats.ppi$num.stressful.events)
roistats.ppi$num.severe.events.sqrt   =sqrt(roistats.ppi$num.severe.events)
roistats.ppi$total.sum.stress.sqrt    =sqrt(roistats.ppi$total.sum.stress)

melted.roistats.ppi=melt(roistats.ppi, id.vars=c("subject",
                                           "Grp", "Gender", "DOB", "MRI", "age.in.years", ## demographic data
                                           "BDI.II",
                                           "num.stressful.events", "num.severe.events", "total.sum.stress", ## SLES data 
                                           "number.of.episodes", "age.of.onset" ## KSADS data
                                           ),
    measure.vars=ppi.clusters.colnames, variable_name="cluster")

melted.roistats.ppi$cluster=factor(melted.roistats.ppi$cluster,
    levels=ppi.clusters.colnames,
    labels=cluster.locations.ppi[ppi.cluster.numbers])

##stop()
my.base.size=14
my_theme=
    theme_bw(base_size =  my.base.size) +
    theme(
        ##legend.position="none",
        legend.position="bottom",        
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),

        ##remove the panel border
        ##panel.border = element_blank(),

        ## add back the axis lines
        axis.line=element_line(colour = "grey50"),
        
        ##axis.title.x=element_blank(),
        axis.title.x = element_text(size=my.base.size, vjust=0),
        axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
        plot.title=element_text(size=my.base.size*1.2, vjust=1))


## Now do the graphing and analysis of the between-group whole-brain task-related results
                                        #if ( 1 == 0 ) {
for ( variable in c("age.in.years", "BDI.II", "num.stressful.events", "num.severe.events", "total.sum.stress", "number.of.episodes", "age.of.onset") ) {

    if ( variable %in% c("number.of.episodes", "age.of.onset") ) {
        roistats.ss=subset(roistats, Grp=="MDD")
        melted.ss=subset(melted.roistats, Grp=="MDD")
    } else {
        roistats.ss=roistats
        melted.ss=melted.roistats
    }
    
    graph=ggplot(melted.ss, aes_string(x=variable, y="value", color="stimulus")) +
        geom_jitter() +
            stat_smooth(method=lm) +
                scale_color_brewer(name="Stimulus", palette="Set1") +
                    facet_grid(.~Grp, scales = "free_x") +
                        labs(y="Percentage Change", title="Left Insula") +
                            my_theme + theme(legend.position="right")
    if (variable == "age.in.years") {
        graph = graph + labs(x="Age (years)")
        model.formula=as.formula(paste("Mean_4", "~", variable, "+", "stimulus", "+", "Grp"))
    } else if ( variable == "num.stressful.events") {
        graph = graph + labs(x="SLES: Number of stressful events")
        model.formula=as.formula(paste("Mean_4", "~", variable, ".sqrt +", "stimulus", "+", "Grp", sep=""))        
    } else if ( variable == "num.severe.events") {
        graph = graph + labs(x="SLES: Number of severe events")
        model.formula=as.formula(paste("Mean_4", "~", variable, ".sqrt +", "stimulus", "+", "Grp", sep=""))        
    } else if ( variable == "total.sum.stress") {
        graph = graph + labs(x="SLES: Total Sum Stress")
        model.formula=as.formula(paste("Mean_4", "~", variable, ".sqrt +", "stimulus", "+", "Grp", sep=""))

    } else if ( variable == "number.of.episodes") {
        graph = graph + labs(x="KSADS: Number of Episodes")
        model.formula=as.formula(paste("Mean_4", "~", variable, " +", "stimulus", sep=""))
    } else if ( variable == "age.of.onset") {
        graph = graph + labs(x="KSADS: Age of First onset")
        model.formula=as.formula(paste("Mean_4", "~", variable, " +", "stimulus", sep=""))
    } else if ( variable == "Duration") {
        graph = graph + labs(x="KSADS: Duration")
        model.formula=as.formula(paste("Mean_4", "~", variable, " +", "stimulus", sep=""))
    }
            
    ## dev.new(); print(graph)
    ggsave(file.path(Group.results.dir, sprintf("%sVsInsulaActivation.happyVsSad.pdf", variable)), graph)
    
    model=lm(model.formula, data=roistats.ss, na.action="na.omit")
    print(anova(model))
    
}
                                        #}
## Now do the graphing and analysis of the PPI-based functional
## connectivity of the regions connected with the insular blob
## identified in teh happyVsSad contrast on the between-group
## whole-brain task-related results

correlations.stack <- stack()
for ( cluster.number in ppi.cluster.numbers ) {

    roi=cluster.locations.ppi[cluster.number]
    
    for ( variable in c("BDI.II", "num.stressful.events", "num.severe.events", "total.sum.stress", "number.of.episodes", "age.of.onset") ) {
##    for ( variable in c("number.of.episodes", "age.of.onset") ) {        

        cat ("####################################################################################################\n")
        cat(sprintf("Correlation of the functional connectivity between the L Claustrum/Insula and %s with %s\n", roi, variable))
        
        roistats.ss=subset(roistats.ppi, Grp=="MDD")
        melted.ss  =subset(melted.roistats.ppi, Grp=="MDD" & cluster==roi)

        ## dump those subjects with ages of onset less than 6 since
        ## they are extreme outliers
#        cat("Rejecting ", sum(roistats.ss$First.Onset <= 9, na.rm=TRUE), " subjects with an age of onset <= 9\n")
#        roistats.ss=subset(roistats.ss, First.Onset > 9)
#        melted.ss  =subset(melted.ss, First.Onset > 9)
        
        graph=ggplot(melted.ss, aes_string(x=variable, y="value", color="Grp")) +
            geom_jitter() +
                stat_smooth(method=lm) +
                    scale_color_brewer(name="Group", palette="Set1") +
                        labs(y="Functional Connectivity (Z-score)", title=sprintf("Claustrum/Insula <-> %s", roi)) +
                            my_theme + theme(legend.position="bottom")

        v=variable
        y=paste("Mean", cluster.number, sep="_")
        if (variable == "age.in.years") {
            graph = graph + labs(x="Age (years)")
            model.formula=as.formula(paste(y, "~", variable))
        } else if ( variable == "BDI.II") {
            graph = graph + labs(x="BDI-II")
            model.formula=as.formula(paste(y, "~", variable, sep=""))        
        } else if ( variable == "num.stressful.events") {
            graph = graph + labs(x="SLES: Number of stressful events")
            model.formula=as.formula(paste(y, "~", variable, ".sqrt ", sep=""))
            v=paste(variable, "sqrt", sep=".")
        } else if ( variable == "num.severe.events") {
            graph = graph + labs(x="SLES: Number of severe events")
            model.formula=as.formula(paste(y, "~", variable, ".sqrt ", sep=""))
            v=paste(variable, "sqrt", sep=".")            
        } else if ( variable == "total.sum.stress") {
            graph = graph + labs(x="SLES: Total Sum Stress")
            model.formula=as.formula(paste(y, "~", variable, ".sqrt ", sep=""))
            v=paste(variable, "sqrt", sep=".")

        } else if ( variable == "number.of.episodes") {
            graph = graph + labs(x="KSADS: Number of episodes")
            model.formula=as.formula(paste(y, "~", variable, sep=""))
        } else if ( variable == "age.of.onset") {
            graph = graph + labs(x="KSADS: Age of First onset")
            model.formula=as.formula(paste(y, "~", variable, sep=""))
        } else if ( variable == "Duration") {
            graph = graph + labs(x="KSADS: Duration")
            model.formula=as.formula(paste(y, "~", variable, sep=""))
        }
        
        dev.new(); print(graph)
        image.filename=file.path(Group.results.dir, sprintf("%sVsInsulaTo%s.connectivity.pdf", variable, gsub(" ", ".", roi, fixed=TRUE)))
        cat("Writing ", image.filename, "\n")
        ggsave(image.filename, graph)
        
        model=lm(model.formula, data=roistats.ss, na.action="na.omit")
        model.summary=summary(model)
        print(model.summary)
        print(anova(model))

        print(cor.test(roistats.ss[, variable], roistats.ss[, y], method="spearman"))
        
        model.estimate=coef(model.summary)[v, "Estimate"]
        model.stderr=coef(model.summary)[v, "Std. Error"]
        model.tvalue=round(coef(model.summary)[v, "t value"], 3)
        model.df=model.summary$df[2]
        model.nparameters=model.summary$df[1]
        n=sum(! is.na(roistats.ss[, variable]))
        nacount=sum(is.na(roistats.ss[, variable]))
        model.pvalue=round(coef(model.summary)[v, "Pr(>|t|)"], 4)
        model.significance=make.significance.indications(model.pvalue)
        csvLine=paste(roi, variable, model.estimate, model.stderr, model.tvalue, n, model.nparameters, model.df, nacount, model.pvalue, model.significance, sep=",")
        push (correlations.stack, csvLine)
    }
}
cat("roi,variable,estimate,stderr,tvalue,n,nparameters,df,nacount,pvalue,significance\n")
l=correlations.stack$value()
for (i in 1:length(l)) {
    cat (l[[i]], "\n")
}
cat("\n")
## model=lm(Mean_4 ~ age.in.years + stimulus + Grp, data=roistats)
## print(anova(model))

## ggsave(file.path(Group.results.dir, "ageVsInsulaActivation.happyVsSad.pdf"), graph)

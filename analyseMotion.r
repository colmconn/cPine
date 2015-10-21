rm(list=ls())
graphics.off()
library(nlme)
library(ggplot2)
library(gdata)

dataRoot="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/"
motionTableHeader=c("roll", "pitch", "yaw", "dS",  "dL",  "dP")

Group.results.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/Group.results"

admin.data.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/admin"
config.dir="/Volumes/PROMISEPEGASUS/yangdata/cPine/data/config"  

groups="mddAndCtrl"

## this file stores the order of the subjects in each of the following BRIK files
#subjectOrderFilename=file.path(Group.data.dir, paste("subjectOrder", groups, seed, "REML.csv", sep="."))
## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice


#cat("*** Reading", subjectOrderFilename, "\n")
#subjectOrder=read.csv(subjectOrderFilename, header=T)
##subjectOrder$subject=as.factor(gsub("_A", "", as.character(subjectOrder$subject), fixed=TRUE))

demographicsFilename=file.path(admin.data.dir, "data_entry_current_021113.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a"))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))
demographics$ID=as.factor(gsub("_[ABCD]", "", as.character(demographics$ID), fixed=FALSE))
demographics$Grp=drop.levels(demographics$Grp)

##stop()

mddSubjects=read.table(file.path(config.dir, "mdd.subjectList.txt"), header=FALSE)
controlSubjects=read.table(file.path(config.dir, "control.subjectList.txt"), header=FALSE)
bothGroups=rbind(mddSubjects, controlSubjects)
allSubjects=data.frame(
    "subjects" = bothGroups,
    "Group"    = demographics[match(gsub("_A", "", bothGroups$V1, fixed=TRUE), demographics$ID), "Grp"]
    )
colnames(allSubjects)=c("subjects", "Group")
allSubjects[allSubjects$subject=="300_A", "Group"]="MDD"
allSubjects$subjects=drop.levels(allSubjects$subjects)

nSubjects=dim(allSubjects)[1]

censor.matrix=matrix(0, nrow=nSubjects, ncol=1)
censor.matrix.rownames=c()

subjectCount=1
for (subject in allSubjects$subject) {
    censorFile=file.path(dataRoot, subject, "functional", paste("censor", subject, "combined_2.1D", sep="_"))
    if (file.exists(censorFile)) {
        cat("*** Attempting to read", censorFile, "\n")
        censorTable=read.table(censorFile, col.names=c("censor"), header=FALSE)
        numberOfCensoredVolumes=sum(1-censorTable[,1])
        totalNumberOfVolumes=length(censorTable[,1])
        censor.matrix[subjectCount, 1] = numberOfCensoredVolumes
        censor.matrix.rownames=c(censor.matrix.rownames, subject)
        
    } else {
        cat("*** No such file", censorFile, "\n")
    }
    subjectCount=subjectCount+1
}
colnames(censor.matrix)="numberOfCensoredVolumes"
censor.matrix.rownames=gsub("_A", "", censor.matrix.rownames, fixed=TRUE)
rownames(censor.matrix)=censor.matrix.rownames

censor.df=as.data.frame(censor.matrix)
censor.df$subject=as.factor(rownames(censor.matrix))
censor.df$Group=demographics[match(censor.matrix.rownames, demographics$ID), "Grp"]
censor.df$Gender=demographics[match(censor.matrix.rownames, demographics$ID), "Gender"]

censor.df$Group=drop.levels(censor.df$Group)
censor.df$drop=ifelse(censor.df$numberOfCensoredVolumes > 0.2*totalNumberOfVolumes, "drop", "keep")
censor.df$drop=as.factor(censor.df$drop)

censor.df[censor.df$subject==300, "Group"]="MDD"
censor.df[censor.df$subject==300, "Gender"]="F"
censor.complete.df=censor.df[censor.df$drop == "keep" , ]

cat("Number of subjects per group\n")

print(table(censor.complete.df$Group))

cat("The following subjects should be dropped:\n")
print(rownames(censor.df)[censor.df$drop=="drop"])
print(as.vector(censor.df[censor.df$drop=="drop", "Group"]))

cat ("*** Only the subjects that can be included in the final analysis\n")
cat (" *** NOTE ***\nThis table does not account for subjects dropped because of missing data.\nTo get the final N per group for your analysis you should use analyseGroupStatsOnly.r and use the tables there\n*** NOTE ***\n")
print(table(censor.df$Group, censor.df$drop))


cat("\n*** NOTE ***\nStatistics refer only to those subjects who are not dropped due to excessive motion and outlier count\n*** NOTE ***\n")

censor.test=wilcox.test(censor.complete.df[censor.complete.df$Group=="MDD", "numberOfCensoredVolumes"], censor.complete.df[censor.complete.df$Group=="NCL", "numberOfCensoredVolumes"])
print(censor.test)

my.base.size=14
my_theme =
    theme_bw(base_size =  my.base.size) +
    theme(##title = "Proportion of risky choices",
          ##legend.position="none",
          ##panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          ##axis.title=element_blank(),
          axis.title.x = element_text(size=my.base.size, vjust=0),
          axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
          plot.title=element_text(size=my.base.size*1.2, vjust=1))

quartz()
censor.boxplot=ggplot(censor.df, aes(x=Group, numberOfCensoredVolumes, fill=Group)) +
    geom_boxplot() +
    xlab("Group") +
    ylab("Number of volumes censored") +
    ggtitle("Volumes censored in the Pine analysis") +
    scale_fill_brewer(palette="Set1") +
    my_theme

print(censor.boxplot)

quartz()
censor.barplot=ggplot(censor.df, aes(x=subject, fill=Group, y=numberOfCensoredVolumes)) +
    geom_bar(stat="identity") +
    scale_fill_brewer(palette="Set1") +
    geom_hline(yintercept=0.2*totalNumberOfVolumes, color="black") +
    annotate("text", x=10, y=(0.2*totalNumberOfVolumes) + 4, label=sprintf("20%% (n=%d)", round(0.2*totalNumberOfVolumes, 0))) +
    my_theme +
    labs(x="Subject", y="Number of volumes censored") +
    ggtitle("Volumes censored in the Pine analysis") +
    theme(axis.text.x=element_text(angle=-90))

print(censor.barplot)


## we need to ammend the list of subjects to drop 121/NCL and 137/MDD
## because these were excluded due to missing data. The subjects to be
## excluded for this reason can be determined using
## analyseGroupStats.r and when that script has finished run
## print(nodata.subjectList.mgd[, c("Gender", "Grp", "subject")]) to
## get the list of excluded subjects

##stop()
##allSubjects =  allSubjects[ - which (allSubjects$subjects %in% c("121_A", "137_A")), ]

censor.df=censor.df[-which (censor.df$subject %in% c(121, 137)), ]
censor.df=censor.df[censor.df$drop=="keep", ]

nSubjects=dim(censor.df)[1]
##allSubjects$subjects=data.frame(subject = drop.levels(allSubjects$subjects))

cat("Number of subjects per group\n")

print(table(censor.df$Group))

cat ("Motion Analyses\n")
for ( func in c("min", "mean", "max")) { 

    motion.matrix=matrix(0, nrow=nSubjects, ncol=6)
    motion.matrix.rownames=c()
    
    subjectCount=1
    for (subject in censor.df$subject) {
        motionFile=file.path(dataRoot, paste(subject, "A", sep="_"), "functional", "motion_demean.1D")
        if (file.exists(motionFile)) {
            ##cat("*** Attempting to read", motionFile, "\n")
            motionTable=read.table(motionFile, header=FALSE)
            
            ## change mean here to max to get the maximum movement excursion
            ## instead of average
            motionExcursion=apply(motionTable, 2, func)
            motion.matrix[subjectCount, ] = motionExcursion
            motion.matrix.rownames=c(motion.matrix.rownames, subject)
            
        } else {
            cat("*** No such file", motionFile, "\n")
        }
        subjectCount=subjectCount+1
    }
    
    colnames(motion.matrix)=motionTableHeader
    motion.matrix.rownames=gsub("_A", "", motion.matrix.rownames, fixed=TRUE)
    rownames(motion.matrix)=motion.matrix.rownames
    
    motion.df=as.data.frame(motion.matrix)
    motion.df$subject=as.factor(rownames(motion.matrix))
    motion.df$Group=demographics[match(motion.matrix.rownames, demographics$ID), "Grp"]
    ##motion.df$BDI.II=demographics[match(motion.matrix.rownames, demographics$ID), "BDI.II"]
    motion.df$Gender=demographics[match(motion.matrix.rownames, demographics$ID), "Gender"]
    ##motion.df$inFinalAnalysis=ifelse(is.na(match(motion.df$subject, gsub("_A", "", subjectsWithTooFewTimePoints$subject, fixed=TRUE))), "included", "excluded")
    ##motion.df$inFinalAnalysis=as.factor(motion.df$inFinalAnalysis)
    
    ## now fix up the 300 subject cos somehow it's got 2 IDs.
    motion.df[motion.df$subject==300, "Group"]="MDD"
    motion.df[motion.df$subject==300, "Gender"]="F"
    ##includedInFinalAnalysis.df=motion.df[motion.df$inFinalAnalysis=="included", ]
    motion.complete.df=motion.df
    
    ## only the subjects included in the final analysis
    cat (sprintf("*** Only the subjects included in the final analysis: %s motion\n", func))
    ## Y=motion.matrix[motion.df$inFinalAnalysis=="included", 1:6]
    ## model.lme=lme(Y ~ Group, random = ~ 1 | subject, data=includedInFinalAnalysis.df)
    ## model.lme.anova=anova(model.lme)
    ## print(model.lme.anova)
    
    model2.lme=lme(cbind(roll, pitch, yaw, dS,  dL, dP) ~ Group, random = ~ 1 | subject, data=motion.complete.df)
    model2.lme.anova=anova(model2.lme)
    print(model2.lme.anova)
    
    ## cat ("*** All subjects\n")
    ## ##motion.df=motion.df[, c("subject", "Group", "Gender", "inFinalAnalysis")]
    ## model.lme=lme(cbind(roll, pitch, yaw, dS,  dL, dP) ~ Group*inFinalAnalysis, random = ~ 1 | subject, data=motion.df) 
    ## model.lme.anova=anova(model.lme)
    ## print(model.lme.anova)
}

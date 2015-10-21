rm(list=ls())
graphics.off()

library(gdata)
library(reshape)
library(ggplot2)

source("scoreMasc.r")

capwords <- function(s, strict = FALSE) {
  cap <- function(s) paste(toupper(substring(s,1,1)),
                           {s <- substring(s,2); if(strict) tolower(s) else s},
                           sep = "", collapse = " " )
  sapply(strsplit(s, split = " "), cap, USE.NAMES = !is.null(names(s)))
}

readSubjectLists <- function( inSubjects ) {
    subjects=do.call(rbind, lapply(subjectListFilenames,
        function(X) {
            data.frame(read.table(X, header=F, sep=""))
        }
        ))

    return (as.vector(subjects[, 1]))
}

buildRoiStatsFileNames <- function (inSubjects, inStimuli ) {

    subjectsByStimuli = expand.grid(inSubjects, inStimuli)

    filenames=apply(subjectsByStimuli, 1,
        function(x) {
            sprintf("%s/%s/functional/roistats.%s.%s.amygdala.mask.txt", data.root, x[1], x[1], x[2])
        }
        )
    return(filenames)
}

readRoiStatsFiles <- function(inRoiStatsFiles, inStimuli) {

    roistats=do.call(rbind, lapply(inRoiStatsFiles,
        function(X) {
            data.frame(read.table(X, header=T, sep=""))
        }
        ))
    
    roistats$stimulus=as.factor(rep(capwords(inStimuli), each=dim(roistats)[1]/length(inStimuli)))
    ##roistats=roistats[, -1]
    
    return(roistats)
}

groups="mddAndCtrl"
fLabel="group.F-value"
usedFwhm=4.2
task="pine"

stimuli=c("fearful", "happy", "neutral", "sad")

root="/Volumes/PROMISEPEGASUS/yangdata/cPine"
data.root=file.path(root, "data")
admin.dir=file.path(data.root, "admin/")
config.dir=file.path(data.root, "config/")    
Group.data.dir=file.path(data.root, "Group.data/")
Group.results.dir=file.path(data.root, "Group.results/")

subjectListFilenames=c(file.path(config.dir, "control.subjectList.txt"), file.path(config.dir, "mdd.subjectList.txt"))
subjectList=readSubjectLists(subjectListFilenames)

## this file stores all of the demographics of interest, such as ID, group, and drug(s) of choice
demographicsFilename=file.path(admin.dir, "data_entry_current_021113.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename,  na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!"), header=T)
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))

### now fix the age

## this complicated looking regexp stuff cleans up the years in DOB
## and MRI with 4 digits to be just 2 digits
demographics$DOB=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", demographics$DOB)
demographics$MRI=sub("([0-9]{1,2})/([0-9]{1,2})/[0-9]{2}([0-9]{2})", "\\1/\\2/\\3", demographics$MRI)

## inData$DOB=as.Date(inData$DOB, "%m/%d/%y")
## inData$MRI=as.Date(inData$MRI, "%m/%d/%y")
## inData$age.in.days=difftime(inData$MRI, inData$DOB, units="days")
## inData$age.in.days=as.numeric(inData$age.in.days)

## inData$age.in.weeks=difftime(inData$MRI, inData$DOB, units="weeks")
## inData$age.in.weeks=as.numeric(inData$age.in.weeks)
## inData$age.in.years=(inData$age.in.weeks)/52

demographics$DOB=as.Date(demographics$DOB, "%m/%d/%y")
demographics$MRI=as.Date(demographics$MRI, "%m/%d/%y")
demographics$age.in.days=difftime(demographics$MRI, demographics$DOB, units="days")
demographics$age.in.days=as.numeric(demographics$age.in.days)

demographics$age.in.weeks=difftime(demographics$MRI, demographics$DOB, units="weeks")
demographics$age.in.weeks=as.numeric(demographics$age.in.weeks)
demographics$age.in.years=(demographics$age.in.weeks)/52

my.base.size=12   
my_theme=
    theme_bw(base_size =  my.base.size) +
    theme(
        ##legend.position="none",
        ##panel.grid.major = element_blank(),
        ##panel.grid.minor = element_blank(),
        ## axis.title.x=element_blank(),
        axis.title.x = element_text(size=my.base.size, vjust=0),
        axis.title.y = element_text(size=my.base.size, vjust=0.4, angle =  90),
        plot.title=element_text(size=my.base.size*1.2, vjust=1))

pd <- position_dodge(.2) # move them .1 to the left and right

roistatsFilenames=buildRoiStatsFileNames(subjectList, stimuli)
roistatsFileExists=file.exists(roistatsFilenames)
numberOfNonexistantRoistatsFiles=sum(!roistatsFileExists)

if (numberOfNonexistantRoistatsFiles > 0 ) {
    cat(sprintf("*** The following %d files do not exist and will be removed from the ROI stats file list\n", numberOfNonexistantRoistatsFiles))
    print(roistatsFilenames[! roistatsFileExists])
}
roistatsFilenames=roistatsFilenames[roistatsFileExists]
roistats=readRoiStatsFiles(roistatsFilenames, stimuli)
## change the file name to be the subject name
roistats$File=as.vector(roistats$File)
roistats$File=sub("(([0-9]*)_[ABCD])(.*)", "\\2", roistats$File)
## rename the File column to subject
colnames(roistats)=c("subject", colnames(roistats)[-1])
## now fix up the mess with 169/300 having two subject IDs
roistats[roistats$subject=="300", "subject"] = "169/300"
## convert back to a factor
roistats$subject=as.factor(roistats$subject)

roistats$Group=demographics[match(roistats$subject, demographics$ID), c("Grp")]
roistats$Group=drop.levels(roistats$Group)

roistats$Sub.brick = NULL

demographics=demographics[match(unique(roistats$subject), demographics$ID), ]

### Now do the MASC scoring

demographics.dim=dim(demographics)
for (r in seq(1, demographics.dim[1]) ) {
  ## cat("##################################################\n")
  subjectNumber=demographics[r, "ID"]
  gender=demographics[r, "Gender"]
  age=round(demographics[r, "age.in.years"], 0)
  old.masc.tscore=demographics[r, "MASC.tscore"]

  ## cat(sprintf("r=%d subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f Old MASC tscore=%0.0f,\n", r, subjectNumber, gender, age, demographics[r, "MASC.total"], old.masc.tscore))
  
  new.masc.tscore=scoreMasc(gender, age, demographics[r, "MASC.total"])
  if (is.na(new.masc.tscore) ) {
    warning(sprintf ("Couldn't set a MASC tscore for subjectNumber=%s gender=%s age=%0.0f MASC Raw Score=%0.0f", subjectNumber, gender, age, demographics[r, "MASC.total"]))
  }
    
  demographics[r, "MASC.tscore"]=new.masc.tscore
  
  ## cat (sprintf("Old MASC tscore=%0.0f, new MASC tscore=%0.0f\n", old.masc.tscore, new.masc.tscore))
}

roistats = cbind(roistats, demographics[match(roistats$subject, demographics$ID), c("MASC.tscore", "CDRS.tscore")])
#colnames(roistats)=c(colnames(roistats)[-length(colnames(roistats))], "MASC.tscore")
## some of the subjects are missing MASC values so drop the incomplete
## cases from the data frame
roistats=roistats[complete.cases(roistats), ]

clusterWhereAmI=c("L Amygdala", "R Amygdala")

clusterCount=length(grep("Mean", colnames(roistats)))

melted.roistats=melt(roistats,
    id.vars=c("subject", "Group", "stimulus", "MASC.tscore", "CDRS.tscore"),
    measure.vars=paste("Mean_", seq(1, clusterCount), sep=""),
    variable_name="cluster")

melted.roistats$cluster=factor(melted.roistats$cluster,
    levels=c(paste("Mean_", seq(1, clusterCount), sep="")),
    labels=paste(seq(1, clusterCount), clusterWhereAmI))

graph.masc.tscore=ggplot(melted.roistats, aes(x=MASC.tscore, y=value, shape=Group, color=Group)) +
    geom_point() +
    stat_smooth(method="rlm") +
    facet_grid(stimulus~cluster, scales="free_y") +
    scale_color_brewer(palette="Set1") +
    ylab("Percentage Signal Change") + xlab("MASC t score") +
    my_theme + theme(legend.position="bottom")

print(graph.masc.tscore)

## graph.cdrsr.tscore=ggplot(melted.roistats, aes(x=CDRS.tscore, y=value, shape=Group, color=Group)) +
##     geom_point() +
##     stat_smooth(method="rlm") +
##     facet_wrap(~cluster, scales="free_y") +
##     scale_color_brewer(palette="Set1") +
##     ylab("Percentage Signal Change") + xlab("CDRS-R t score") +
##     my_theme + theme(legend.position="bottom")

## quartz()
## print(graph.cdrsr.tscore)

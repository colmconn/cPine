rm(list=ls())

fixSubjectOrderTable <- function (inSubjectOrderTable) {
    inSubjectOrderTable$subject=gsub("_A", "", as.character(inSubjectOrderTable$subject), fixed=TRUE)
    inSubjectOrderTable[inSubjectOrderTable$subject=="300", "subject"] = "169/300"
    inSubjectOrderTable$subject=as.factor(inSubjectOrderTable$subject)

    return(inSubjectOrderTable)
}

if ( Sys.info()["sysname"] == "Darwin" ) {
    root.dir="/Volumes/data"
} else if ( Sys.info()["sysname"] == "Linux" ) {
    root.dir="/data"
} else {
    cat(paste("Sorry can't set data directories for this computer\n"))
}

Admin.data.dir=normalizePath(file.path(root.dir, "sanDiego/cPine/data/admin"))
Group.data.dir=normalizePath(file.path(root.dir, "sanDiego/cPine/data/Group.data"))

subjectOrder.filename=file.path(Group.data.dir, "subjectOrder.mddOnly.fearful.REML.csv")
subjectOrder=read.csv(subjectOrder.filename, header=TRUE)
subjectOrder=fixSubjectOrderTable(subjectOrder)

demographicsFilename=file.path(Admin.data.dir, "0-data_entry_current_10152013.csv")
cat("*** Reading", demographicsFilename, "\n")
demographics=read.csv(demographicsFilename, header=T, na.strings = c("NA", "<NA>", "#N/A", "#VALUE", "#VALUE!", "n/a", "N/A", "#DIV/0!", "IGNORE THIS SUBJECT"))
cat(sprintf("*** Read demographics data for %s unique subjects\n",  length(unique(demographics$ID))))
demographics$ID=as.factor(gsub("_[ABCD]", "", as.character(demographics$ID), fixed=FALSE))


ksads.filename=file.path(Admin.data.dir, "KSADS.age.of.onset.pine.csv")
cat("Reading KSADS data from: ", ksads.filename, "\n")
ksads=read.csv(ksads.filename, header=TRUE)
ksads$Subject.Number=gsub("A", "_A", ksads$Subject.Number, fixed=TRUE)
##ksads[ksads$Subject.Number=="300", "Subject.Number"] = "169/300"
cat(sprintf("Read %d subjects from the KSADS file\n", dim(ksads)[1]))


ksads.filename=file.path(Admin.data.dir, "KSADSDiagnoses_tiffany.csv")
cat("Reading KSADS data from: ", ksads.filename, "\n")
ksads.tiff=read.csv(ksads.filename, header=TRUE, skip=1)
cat(sprintf("Read %d subjects from the KSADS file\n", dim(ksads.tiff)[1]))

binded.ksads=
    cbind(
        ksads[ , c("Subject.Number", "age.of.onset", "number.of.episodes")],
        ksads.tiff[match(ksads$Subject.Number, ksads.tiff$Subject.Number), c("First.Onset", "X..Episodes")])
colnames(binded.ksads)=c("Subject.Number", "Eva.age.of.onset", "Eva.number.of.episodes", "Tiff.First.Onset", "Tiff.X..Episodes")

binded.ksads$Subject.Number=gsub("_A", "", binded.ksads$Subject.Number, fixed=TRUE)
binded.ksads[binded.ksads$Subject.Number=="300", "Subject.Number"] = "169/300"

binded.ksads=
    cbind(
        demographics[match(binded.ksads$Subject.Number, demographics$ID), c("Grp", "Gender")],
        binded.ksads)

mdd.binded.ksads=subset(binded.ksads, Grp=="MDD")

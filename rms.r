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

makeTableString <- function(inGroup, inMean, inSe, inMin, inMax, inNaCount, inMissingData=TRUE) {
  ##  st=paste(round(inMean, 1), " / ", round(inSe, 1),    

  st=paste(round(inMean, 1), " ± ", round(inSe, 1),
    " (", round(inMin, 1), "-", round(inMax, 1), ")", ifelse(inMissingData & inNaCount > 0, paste(" [", inNaCount, "]", sep=""), ""), sep="")
  return(st)
}


runModelStatistics <- function (inIndividualModels, inCreateBoxplots=FALSE, inAlpha=0.05) {

    corrected.p.value=inAlpha / ( length(levels(inIndividualModels$Group)) * length(levels(inIndividualModels$cluster)) )
    cat (sprintf("*** Bonferroni correct p value is %0.4f\n", round(corrected.p.value, 4)))

    parameters=c("delay", "k", "rise", "decay", "timeOfExtrema", "extrema")
    mystack <- stack()
    header="Parameter,Cluster,Stimulus,Control,MDD,DF,Stat.,pValue,Signif."

    for ( param in parameters ) {
        individualModels.summarySE=summarySE(individualModels, measurevar=param, groupvars=c("stimulus", "cluster", "Group"), na.rm=TRUE)
        for ( cl in levels(inIndividualModels$cluster) ) {
            for ( stim in levels(inIndividualModels$stimulus )) {

                ctrl.string=""
                mdd.string=""
                test=""
                sm.df=subset(individualModels.summarySE, stimulus==stim & cluster==cl)
                ctrl.string=makeTableString(sm.df[2, 1], inMean=sm.df[2, "median"],  sm.df[2, "mad"], sm.df[2, "min"], sm.df[2, "max"], sm.df[2, "nacount"], inMissingData=TRUE)
                mdd.string =makeTableString(sm.df[1, 1], inMean=sm.df[1, "median"],  sm.df[1, "mad"], sm.df[1, "min"], sm.df[1, "max"], sm.df[1, "nacount"], inMissingData=TRUE)
                
                cat ("####################################################################################################\n")
                cat (sprintf("Wilcox test for cluster: %s stimulus: %s parameter: %s\n", cl, stim, param ))
                
                test <- wilcox.test(inIndividualModels[inIndividualModels$Group=="NCL" & inIndividualModels$cluster==cl, param],
                                    inIndividualModels[inIndividualModels$Group=="MDD" & inIndividualModels$cluster==cl, param])
                
                if( inherits(test <- try(wilcox.test(inIndividualModels[inIndividualModels$Group=="NCL" & inIndividualModels$cluster==cl, param],
                                                     inIndividualModels[inIndividualModels$Group=="MDD" & inIndividualModels$cluster==cl, param]),
                                         silent=FALSE),
                             "try-error") ) {
                    test <- 0
                }
                
                var.statistic=""
                var.pvalue=""
                var.parameter=""
                var.significance=""
                
                if (is.list(test)) {
                    var.statistic=round(test$statistic, 2)
                    var.parameter=ifelse(is.null(test$parameter), "NA", round(test$parameter, 2))
                    var.pvalue=round(test$p.value, 2)
                    var.significance=make.significance.indications(test$p.value)
                } ## end of if (is.list(test)) {
                print (test)
                print(sm.df)
                st=paste(param, cl, stim, ctrl.string, mdd.string,
                    var.parameter, var.statistic, var.pvalue, var.significance, sep=",")
                push(mystack, st)
            } ## end of for ( stim inlevels(melted.individual.models$stimulus )) {
        } ## end of for ( cl in levels(melted.individual.models$cluster) ) {
    } ## end of for ( param in parameters ) {
    
    cat("################################################################################\n");
    cat("Summary statistics table\n")  
    l=mystack$value()
    cat(header, "\n")
    for (i in 1:length(l)) {
        cat (l[[i]], "\n")
    }
} ## end of runModelStatistics <- function (inNndividualModels) {

graphParameterBoxplots <- function (inIndividualModels, inContrast) {

    imageDirectory=file.path(Group.results.dir, inContrast)
    if ( ! file.exists(imageDirectory) ) {
        dir.create(imageDirectory)
    }

    parameters=c("delay", "k", "rise", "decay", "timeOfExtrema")
    melted.individual.models=melt(inIndividualModels,
        id.vars=c("subject", "Group", "cluster", "stimulus"),
        measure.vars=parameters,
        variable_name="parameter")

    for (cl in levels(melted.individual.models$cluster) ) {

        df.cl=subset ( melted.individual.models, cluster %in% cl )
        
        graph=ggplot(df.cl, aes(x=parameter, y=value)) +
            geom_boxplot() +
                scale_y_log10() +
                facet_grid ( stimulus ~ Group ) +
                    ggtitle (substituteShortLabels(cl)) +
                        xlab("Parameter") +
                            ylab("Value") +
                                my_theme + theme(legend.position="bottom")
        quartz()
        print(graph)
        stop()
        imageFilename=file.path(imageDirectory, sprintf("boxplot.irf.%s.%s.%s.pdf", gsub(" +", ".", level),  task, inContrast, fLabel))
        cat(paste("*** Creating individual boxplot image named", imageFilename, "\n"))
        ggsave(imageFilename, graph)

    } ## end of for (cl in levels(melted.individual.models$cluster) ) {
} ## end of graphParameterBoxplots <- function (inIndividualModels) {

parameters=c("delay", "k", "rise", "decay", "timeOfExtrema")
melted.individual.models=melt(individualModels,
        id.vars=c("subject", "Group", "cluster", "stimulus"),
        measure.vars=parameters,
        variable_name="parameter")

groupAverageRoistats=ddply(melted.individual.models, .(stimulus, cluster, Group, parameter), summarize,
    "min"     = min(value, na.rm=TRUE),
    "1st Qu." = quantile(value, 0.25, na.rm=TRUE),
    "Median"  = median(value, na.rm=TRUE),
    "Mean"    = mean(value, na.rm=TRUE),
    "3rd Qu." = quantile(value, 0.75, na.rm=TRUE),
    "Max."    = max(value, na.rm=TRUE),
    "NA's"    = sum(is.na(value))
    )

ddply (s107, .(subject, stimulus, Group),
       .fun=colwise(
           .fun=function (xx) {
               diff(xx, 1)
           },
           c(paste("Mean_", seq(1, 3), sep=""))
           )
       )

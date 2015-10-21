rm(list=ls())
graphics.off()
library(ggplot2)
library(plyr)
library(minpack.lm)

gammaVar <- function (x, inDelay, inK, inRise, inDecay, inIntercept=0) {
    retVal=rep(0.0,length(x))
    ## which indices of x have positive differences with inDelay, no
    ## negative logs
    idx=which (! x-inDelay < 0 )
    if ( inDecay >= 0 ) {
        retVal[idx]=inK * exp( log(x[idx]-inDelay) * inRise ) * exp(-(x[idx]-inDelay)/inDecay) - inIntercept
    }
    return(retVal)
}

## gammaVar <- function (x, inDelay, inK, inRise, inDecay) {
##     retVal=rep(0.0,length(x))
##     ## which indices of x have positive differences with inDelay, no
##     ## negative lags allowed
##     idx=which (! x-inDelay < 0 )
##     if ( inDecay >= 0 ) {
##         retVal[idx]=inK * exp( log(x[idx]-inDelay) * inRise ) * exp(-(x[idx]-inDelay)/inDecay)
##     }
##     return(retVal)
## }

##ctrl=nls.control(maxiter = 2048, warnOnly=TRUE)
ctrl=nls.lm.control(maxiter = 1e5)

## time=0:14
## b=8.6
## c=0.547
## b=7
## c=0.4
## hrf=time^b * exp(-time/c) + rnorm(length(time), 0, 0.1)
## df=data.frame(time=time, hrf=hrf)
## p=ggplot(data=df, aes(x=time, y=hrf) ) + geom_line() + geom_point() + ggtitle("Simpler Cohen HRF")
## quartz()
## print(p)


## cat("Trying to use nls to fit the simpler Cohen model\n")
## nlsModel=try(nls(hrf ~ time ^ b * exp(-time/c),
##     start=list(b=7, c=0.4),
##     data=df,
##     control=ctrl))

## if ( ! inherits(nlsModel,  "try-error") ) {
##     nt=seq(0, 14, by=0.1)
##     fit=predict(nlsModel, list(time=nt))
##     predicted.df=data.frame(time=nt, hrf=fit)
##     p = p + geom_line(data=predicted.df, color="blue")
##     print(p)
## }

## print(summary(nlsModel))

## ## cat(":Trying to use nlmrt to perform the fit\n")
## ## nlmrtModel=try(nlxb(hrf ~ k* time ^ b * exp(-time/c),
## ##     start=list(b=7, c=0.4, k=1), data=df, trace=TRUE))

## ## m.2 <- nls(y ~ rhs(x, intercept, power), data = ds, start = list(intercept = 0,
## ## + power = 2), trace = T)


## cat("Trying to use nls to fit the more complex HRF\n")
## nlsComplexModel=try(
##     nls(hrf ~ gammaVar(time, delay, k, rise, decay),
##         start=list(delay=1, k=2, rise=2, decay=1.5), data=df,
##         lower=c(0, -Inf, -Inf, 0),
##         algorithm="port",
##         control=ctrl
##         )
##     )
## if ( ! inherits(nlsComplexModel,  "try-error") ) {
##     nt=seq(0, 14, by=0.1)
##     fit=predict(nlsComplexModel, list(time=nt))
##     predicted.df=data.frame(time=nt, hrf=fit)
##     p = p + geom_line(data=predicted.df, color="red")
##     print(p)
## }

## print(summary(nlsComplexModel))

cat ("Now trying with real data, fingers crossed\n")
roistats=structure(list(subject = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 4L, 4L, 4L, 4L, 4L, 4L, 4L, 4L), .Label = c("107", "108", 
"109", "111", "112", "114", "122", "123", "127", "132", "134", 
"135", "138", "139", "141", "142", "143", "144", "145", "147", 
"148", "151", "152", "153", "155", "157", "159", "160", "162", 
"163", "165", "169/300", "303", "304", "306", "307", "310", "312", 
"317", "318", "319", "320", "321", "322", "323", "324", "325", 
"328", "330", "331", "332", "336", "337", "339", "343", "344", 
"346", "347", "348", "349", "354", "356", "357", "358", "359", 
"360", "361", "362", "363", "367", "373", "376"), class = "factor"), 
    Sub.brick = c(0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 0L, 1L, 2L, 
    3L, 4L, 5L, 6L, 7L, 0L, 1L, 2L, 3L, 4L, 5L, 6L, 7L, 0L, 1L, 
    2L, 3L, 4L, 5L, 6L, 7L), Mean_4 = c(-33.367175, -26.732865, 
    -4.206342, 50.478874, -19.778826, -4.578424, 39.251347, -36.009782, 
    20.861298, -14.821766, -17.755583, -9.803989, 4.717233, -10.829425, 
    3.709272, 6.902878, -7.893618, -32.548031, -10.36734, 29.95981, 
    -3.530568, 18.217021, 21.265448, 21.878167, 12.078371, 1.575327, 
    8.084844, -0.829687, -22.570601, -44.547768, 6.730416, -6.574744
    ), stimulus = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
    1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
    2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L), .Label = c("Fearful", 
    "Sad"), class = "factor"), Group = structure(c(2L, 2L, 2L, 
    2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 
    2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("MDD", 
    "NCL"), class = "factor")), .Names = c("subject", "time", 
"hrf", "stimulus", "Group"), row.names = c(1L, 2L, 3L, 4L, 
5L, 6L, 7L, 8L, 289L, 290L, 291L, 292L, 293L, 294L, 295L, 296L, 
577L, 578L, 579L, 580L, 581L, 582L, 583L, 584L, 865L, 866L, 867L, 
868L, 869L, 870L, 871L, 872L), class = "data.frame")


## The input to 3dDeconvolve is scalled to a mean of 10,000 so the
## results of 3dDeconvolce are too big by a factor of 100, divide by
## 100 here to correct for this
roistats$hrf=(roistats$hrf/100)
## multiply by 2 to convert from TRs to seconds
roistats$time=roistats$time * 2

nonLinearModel <- function ( inDf ) {
    try(
        nlsLM(hrf ~ gammaVar(time, delay, k, rise, decay),
              start=list(delay=1, k=1, rise=2, decay=1.5),
              data=inDf,
              ##lower=c(0, -Inf, -Inf, 0),
              ##upper=c(Inf, 1, Inf, Inf),
              ## algorithm="port",
              control=ctrl
              )
        )
    ## try(
    ##     nlsLM(hrf ~ gammaVar(time, 0, k, rise, decay, intercept),
    ##         start=list( k=1, rise=2, decay=1.5, intercept=0),
    ##         data=inDf,
    ##           ##lower=c(0, -Inf, -Inf, 0),
    ##         ##upper=c(Inf, 1, Inf, Inf),
    ##         ## algorithm="port",
    ##         control=ctrl
    ##         )
    ##     )
}


## apply the nonLinearModel function for each combination of subject and stimlus
complexHrfModels=dlply( roistats, .(subject, stimulus), nonLinearModel)

##stop("stopping")

## code for plotting below here

subjectsAndStimuli=strsplit(names(complexHrfModels), ".", fixed=TRUE)
subjectsToGroupMap=unique(roistats[, c("subject", "Group")])
nt=seq(0, max(roistats$time), by=0.1)
numberOfNewTimepoints=length(nt)

## lapply(complexHrfModels, function(z) { attributes(deparse(substitute(z)))$names  } )

## each element of the complexHrfModels list is an object of class nls
## iterate over the connents of this list and generate a predicted fit
## of the model to new time data. Thisblock of code returns a data
## frame organized just like roistats, that it is it contains the
## following columns: subject, time, hrf (the fitted data) stimulus
## type (fearful, sad, etc) and the group of the subject

fitModels <- function(ii) {
    subject=subjectsAndStimuli[[ii]][1]
    stimulus=subjectsAndStimuli[[ii]][2]
    group=as.vector(subjectsToGroupMap[subjectsToGroupMap$subject==subject, "Group"])

    cat(sprintf("subject: %s stimulus: %s\n", subject, stimulus))
    
    if ( ! inherits (complexHrfModels[[ii]], "try-error") ) {
        fit=predict(complexHrfModels[[ii]], list(time=nt))
    } else {
        warning(sprintf("Got a try-error for %s %s\n", subject, stimulus))
        fit=rep(0, numberOfNewTimepoints)
    }
    
    ## construct the data frame to return
    predicted.df=data.frame(
        subject=rep(subject, numberOfNewTimepoints),
        time=nt,
        hrf=fit,
        stimulus=rep(stimulus, numberOfNewTimepoints),
        group=rep(group, numberOfNewTimepoints)
        )
    return(predicted.df)
}

fittedModels=do.call(rbind, lapply(seq_along(complexHrfModels), fitModels))


cat("Plotting the actual data\n")
allData.plot=ggplot(roistats, aes(x=time, y=hrf, color=subject)) +
    geom_point() +
    geom_line(linetype=3, size=1.2) +
    geom_line(data=fittedModels) +
    ggtitle("Real Data") +
    xlab("Time (s)") +
    ylab("IRF") +
    facet_wrap(~stimulus, scales="free_y")
quartz()
print(allData.plot) 

lapply(seq_along(complexHrfModels),
       function(ii) {

           subject=subjectsAndStimuli[[ii]][1]
           stimulus=subjectsAndStimuli[[ii]][2]
           group=as.vector(subjectsToGroupMap[subjectsToGroupMap$subject==subject, "Group"])
           cat("##################################################\n")
           cat(sprintf("*** Model summary for %s in the %s group for %s stimulus\n", subject, group, stimulus))
           print(complexHrfModels[[ii]])
       }
       )



#roistats$hrf=roistats$hrf/max(roistats$hrf)

## cat("Trying to use nls to fit the simpler Cohen model to the read data\n")
## roistats.nlsModel=try(nls(hrf ~ time ^ b * exp(-time/c),
##     start=list(b=7, c=0.4),
##     data=roistats,
##     control=ctrl))

##print(summary(roistats.nlsModel))
## print("Plotting the actual data\n")
## roistats.p=ggplot(data=roistats, aes(x=time, y=hrf) ) + geom_line() + geom_point() + ggtitle("Real Data")
## quartz()
## print(roistats.p)

## if ( ! inherits(roistats.nlsModel,  "try-error") ) {
##     nt=seq(0, 14, by=0.1)
##     fit=predict(roistats.nlsModel, list(time=nt))
##     predicted.df=data.frame(time=nt, hrf=fit)
##     roistats.p = roistats.p + geom_line(data=predicted.df, color="blue")
##     print("Plotting the simpler Cohen gamma model\n")
##     print(roistats.p)
## }

## cat("Trying to use nls to fit the more complex HRF\n")
## roistats.nlsComplexModel=try(
##     nls(hrf ~ gammaVar(time, 0, k, rise, decay, intercept),
##         start=list(k=max(roistats$hrf), rise=2, decay=1.5, intercept = 0), data=df,
##         ##lower=c(0, -Inf, -Inf, 0),
##         algorithm="port",
##         control=ctrl
##         )
##     )
## if ( ! inherits(roistats.nlsComplexModel,  "try-error") ) {
##     nt=seq(0, 14, by=0.1)
##     fit=predict(roistats.nlsComplexModel, list(time=nt))
##     predicted.df=data.frame(time=nt, hrf=fit)
##     roistats.p = roistats.p + geom_line(data=predicted.df, color="red")
##     print("Plotting the more complex gamma model\n")    
##     print(roistats.p)
## }

## print(summary(roistats.nlsComplexModel))






## k=2
## delay=1
## rise=2
## decay=1.5
## complexHrf=gammaVar(time, delay, k, rise, decay) + rnorm(length(time), 0, 0.1)
## complex.df=data.frame(time=time, hrf=complexHrf)
## pc=ggplot(data=complex.df, aes(x=time, y=hrf) ) + geom_line() + geom_point() + ggtitle("Complex HRF")
## quartz()
## print(pc)



## version
## packageVersion("nlmrt")

rm(list=ls())
##graphics.off()
library(ggplot2)

##well behaved data
## roistats=structure(list(subject = structure(c(4L, 4L, 4L, 4L, 4L, 4L, 
## 4L, 4L), .Label = c("107", "108", "109", "111", "112", "114", 
## "122", "123", "127", "132", "134", "135", "138", "139", "141", 
## "142", "143", "144", "145", "147", "148", "151", "152", "153", 
## "155", "157", "159", "160", "162", "163", "165", "169/300", "303", 
## "304", "306", "307", "310", "312", "317", "318", "319", "320", 
## "321", "322", "323", "324", "325", "328", "330", "331", "332", 
## "336", "337", "339", "343", "344", "346", "347", "348", "349", 
## "354", "356", "357", "358", "359", "360", "361", "362", "363", 
## "367", "373", "376"), class = "factor"), time = c(0, 2, 4, 6, 
## 8, 10, 12, 14), hrf = c(0.20861298, -0.14821766, -0.17755583, 
## -0.09803989, 0.04717233, -0.10829425, 0.03709272, 0.06902878), 
##     stimulus = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("Fearful", 
##     "Sad"), class = "factor"), Group = structure(c(1L, 1L, 1L, 
##     1L, 1L, 1L, 1L, 1L), .Label = c("MDD", "NCL"), class = "factor")), .Names = c("subject", 
## "time", "hrf", "stimulus", "Group"), row.names = 289:296, class = "data.frame")

## not so well behaved data

## roistats=structure(list(subject = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 
## 1L, 1L), .Label = c("107", "108", "109", "111", "112", "114", 
## "122", "123", "127", "132", "134", "135", "138", "139", "141", 
## "142", "143", "144", "145", "147", "148", "151", "152", "153", 
## "155", "157", "159", "160", "162", "163", "165", "169/300", "303", 
## "304", "306", "307", "310", "312", "317", "318", "319", "320", 
## "321", "322", "323", "324", "325", "328", "330", "331", "332", 
## "336", "337", "339", "343", "344", "346", "347", "348", "349", 
## "354", "356", "357", "358", "359", "360", "361", "362", "363", 
## "367", "373", "376"), class = "factor"), time = c(0, 2, 4, 6, 
## 8, 10, 12, 14), hrf = c(-0.33367175, -0.26732865, -0.04206342, 
## 0.50478874, -0.19778826, -0.04578424, 0.39251347, -0.36009782
## ), stimulus = structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), .Label = c("Fearful", 
## "Sad"), class = "factor"), Group = structure(c(2L, 2L, 2L, 2L, 
## 2L, 2L, 2L, 2L), .Label = c("MDD", "NCL"), class = "factor")), .Names = c("subject", 
## "time", "hrf", "stimulus", "Group"), row.names = c(NA, 8L), class = "data.frame")

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


gammaVar <- function (x, inDelay, inK, inRise, inDecay) {
    retVal=rep(0.0,length(x))
    ## which indices of x have positive differences with inDelay, no
    ## negative logs
    idx=which (! x-inDelay < 0 )
    if ( inDecay >= 0 ) {
        retVal[idx]=inK * exp( log(x[idx]-inDelay) * inRise ) * exp(-(x[idx]-inDelay)/inDecay)
    }
    return(retVal)
}

gammaVarSSE <- function (par, x, y) {
    delay=par[1]
    k=par[2]
    rise=par[3]
    decay=par[4]

    yhat=gammaVar(x, delay, k, rise, decay)

    sse=sum((y - yhat)^2)

    return(sse)
}

nonLinearModel <- function ( inDf ) {
    ##    try(
    model=optim(c(1, 1, 2, 1.5), gammaVarSSE,
        method="L-BFGS-B",
        lower=c(-Inf, -Inf, -Inf, 0),
        hessian=TRUE, x=inDf$time, y=inDf$hrf)
    ## model=powell(c(1, 1, 2, 1.5), gammaVarSSE,
    ##     x=inDf$time, y=inDf$hrf)
    ##        )
}

## apply the nonLinearModel function for each combination of subject and stimlus
complexHrfModels=dlply( roistats, .(subject, stimulus), nonLinearModel)

subjectsAndStimuli=strsplit(names(complexHrfModels), ".", fixed=TRUE)
subjectsToGroupMap=unique(roistats[, c("subject", "Group")])
nt=seq(0, max(roistats$time), by=0.1)
numberOfNewTimepoints=length(nt)

fitModels <- function(ii) {
    subject=subjectsAndStimuli[[ii]][1]
    stimulus=subjectsAndStimuli[[ii]][2]
    group=as.vector(subjectsToGroupMap[subjectsToGroupMap$subject==subject, "Group"])

    cat(sprintf("subject: %s stimulus: %s\n", subject, stimulus))

    model=complexHrfModels[[ii]]

    if ( ! inherits (model, "try-error") ) {
        fit=gammaVar(nt, model$par[1], model$par[2], model$par[3], model$par[4])
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

roistats$hrf2=roistats$hrf+2
ddply(roistats, .(time, stimulus, Group),
      .fun=colwise(
          .fun=function (xx) {
              c(mean=mean(xx))
          },
          c("hrf", "hrf2"))
      )

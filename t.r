findExtrema <- function ( inDf ) {
    retVal = c(time=NA, value=NA)
    tp=turnpoints(inDf$value)
    if ( length(tp$tppos) > 1 ) {
        print(tp$tppos)
        stop("Can't handle more than one turnpoint in the model HFR")
    } else {
        if (is.na(tp$tppos)) {
            retVal = c(time=NA, value=NA)
        } else {
            retVal = c(time=inDf[tp$tppos, "time"], value=inDf[tp$tppos, "value"])
        }
    }
    return(retVal)
}


for ( sj in levels(fittedModels$subject) ) {
    for ( stim in levels(fittedModels$stimulus ) ) {
        for ( cl in  levels(fittedModels$cluster ) ) {
            cat(sprintf("subject = %s, stim = %s, cl = %s\n", sj, stim, cl))

            tp=turnpoints(fittedModels[ fittedModels$subject==sj &  fittedModels$stimulus==stim & fittedModels$cluster==cl, 'value'])
            if ( length(tp$tppos) > 1 ) {
                print(tp$tppos)
                stop("Can't handle more than one turnpoint in the model HFR")
            } else {
                if (is.na(tp$tppos)) {
                    retVal = c(time=NA, value=NA)
                } else {
                    retVal = c(time=fittedModels[tp$tppos, "time"], value=fittedModels[tp$tppos, "value"])                    
                }
            }
            print(retVal)
            ##findExtrema(fittedModels[fittedModels$subject==sj & fittedModels$stimulus==stim & fittedModels$cluster==cl, ] )
        }
    }
}

peaks=ddply(fittedModels, groupVars, findExtrema )

deltaT=0.1
x = seq(0, 10, by=deltaT)

f1=0.1
f2=1.2
h1=1
h2=2
h3=4
h4=4
k=1


    ## t1=[0:sampling_rate:h1];
    ## y1=f1/2*cos(2*pi*1/(2*h1).*t1)-f1/2;            #initial dip
    
    ## t2=[0:sampling_rate:h2-sampling_rate];
    ## y2=-(1+f1)/2*cos(2*pi*1/(2*h2).*t2)+1/2-f1/2;   #rise to peak
    
    ## t3=[0:sampling_rate:h3-sampling_rate];
    ## y3=f2/2*cos(2*pi*1/(2*h3).*t3)+1-f2/2;          #drop and undershoot
    
    ## t4=[0:sampling_rate:h4-sampling_rate];
    ## y4=-(f2-1)/2*cos(2*pi*1/(2*h4).*t4)-(f2-1)/2;   #recovery to baseline
    
    ## t_total=[t1 t2+h1 t3+h1+h2 t4+h1+h2+h3];
    ## y_total=[y1 y2 y3 y4];

    

t1=seq(0, h1, by=deltaT)
y1=f1/2*cos(2*pi*1/(2*h1)*t1)-f1/2            #initial dip

t2=seq(0, h2-deltaT, by=deltaT)
y2=-(1+f1)/2*cos(2*pi*1/(2*h2)*t2)+1/2-f1/2   #rise to peak

t3=seq(0, h3-deltaT, by=deltaT)
y3=f2/2*cos(2*pi*1/(2*h3)*t3)+1-f2/2          #drop and undershoot

t4=seq(0, h4-deltaT, by=deltaT)
y4=-(f2-1)/2*cos(2*pi*1/(2*h4)*t4)-(f2-1)/2   #recovery to baseline

t_total=c(t1, t2+h1, t3+h1+h2, t4+h1+h2+h3);
y_total=k * c(y1, y2, y3, y4)

df=data.frame(
    time = t_total,
    y=y_total)

graph=ggplot(data=df, aes(x=time, y=y)) +
    geom_line() +
    geom_vline(xintercept=h1,                color="red",   linetype=2) +
    geom_vline(xintercept=h1 + h2,           color="green", linetype=2) +
    geom_vline(xintercept=h1 + h2 + h3,      color="blue",  linetype=2) +
    geom_vline(xintercept=h1 + h2 + h3 + h4, color="green", linetype=2) +
    ggtitle("Cosine 4 HRF") + xlab("Time") + ylab("HRF")
print(graph)

stop()


cosine4 <- function (x, f1, f2, h1, h2, h3, h4, k, deltaT=0.1  ) {
    t1=seq(0, h1, by=deltaT)
    y1=f1/2*cos(2*pi*1/(2*h1)*t1)-f1/2            #initial dip
    
    t2=seq(0, h2-deltaT, by=deltaT)
    y2=-(1+f1)/2*cos(2*pi*1/(2*h2)*t2)+1/2-f1/2   #rise to peak
    
    t3=seq(0, h3-deltaT, by=deltaT)
    y3=f2/2*cos(2*pi*1/(2*h3)*t3)+1-f2/2          #drop and undershoot
    
    t4=seq(0, h4-deltaT, by=deltaT)
    y4=-(f2-1)/2*cos(2*pi*1/(2*h4)*t4)-(f2-1)/2   #recovery to baseline
    
    ## t_total=c(t1 t2+h1 t3+h1+h2 t4+h1+h2+h3);
    y_total=k * c(y1, y2, y3, y4)[1:length(x)]

}



## par a vector of parameters
## x a vector of times
## y a vector of values of the real data at each of the time points in x

cosine4SSE <- function (par, x, y) {
    ##           x,  f1     f2      h1      h2      h3      h4      k
    yhat=cosine4(x, par[1], par[2], par[3], par[4], par[5], par[6], par[7])

    sse=sum((y - yhat)^2)

    return(sse)
}


optimizeModelParameters <- function ( inDf ) {
    if ( ! all(c("time", "value") %in% colnames(inDf) ) ) {
        stop ("*** optimizeModelParameters: The column names of inDf do not contain (time, value)\n")
    }
    
    model=try(
        optim(c(1, 1, 2, 1.5), consien4SSE,
              method="L-BFGS-B",
              hessian=TRUE, x=inDf$time, y=inDf$value)
        ## model=powell(c(1, 1, 2, 1.5), gammaVarSSE,
        ##     x=inDf$time, y=inDf$hrf)
        ##        )
        )
    return(model)
}

df=data.frame(
    time = seq(0, 10, by=0.1),
    y=cosine4(time, 0.1, 0.3, 1, 6, 6, 4, 1))

graph=ggplot(data=df, aes(x=time, y=y)) +
    geom_line() +
    ggtitle("Cosine 4 HRF") + xlab("TIme") + ylab("HRF")

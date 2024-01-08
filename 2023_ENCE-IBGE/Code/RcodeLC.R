##---------- The Lee-Carter method -----------------
##
##  Modelling and forecasting mortality with the 
##  Lee-Carter method
##
##  Code prepared for the Short course for the 
##  DAAD-CAPES(PROBAL) research collaboration
## 
##  Author: Ugofilippo Basellini
##  Date: 12 December 2023
## 
##--------------------------------------------------

## cleaning the workspace
rm(list=ls(all=TRUE))

## load useful packages
library(demography)
library(graphics)
library(fields)
library(viridis)

## loading life-table function
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("LifetableMX.R")

## loading England & Wales data
load("GBRTENW.Rdata")

## dimensions of matrices
m <- length(a)
n <- length(t)
n1 <- length(t1)
nF <- length(tF)

## plotting observed log-mortality
par(mar=c(3,3,1,1))
matplot(a,LMX,t="l",col=viridis(n1),lty=1,xlab="",ylab="")
mtext("ages",side=1,line=1.5);mtext("log-mortality",side=2,line=2)
rany <- range(LMX,finite=T)
for(i in 1:n1){
  plot(a, LMX[,i], col=viridis(n1)[i], pch=16, ylim=rany, main=t1[i],xlab="",ylab="")
  mtext("ages",side=1,line=1.5);mtext("log-mortality",side=2,line=2)
  Sys.sleep(0.1)
}

## observed life expectancy
e0obs <- apply(exp(LMX),2,e0.mx,x=a,sex="F")
plot(t1,e0obs,pch=16,col="blue",xlab="",ylab="")
mtext("years",side=1,line=1.5);mtext("life exp",side=2,line=2)


##---- step 1: derive ALPHA -----

## compute the mean of observed rates
Alpha <- apply(LMX,1,mean)

## plot
matplot(a,LMX,t="l",lty=1,lwd=0.8,col="grey70",xlab="",ylab="")
mtext("ages",side=1,line=1.5);mtext("log-mortality",side=2,line=2)
lines(a,Alpha,col="darkgreen",lwd=2)

##---- step 2: derive BETA and KAPPA -----

## derive matrix of residuals ("centred" mortality rates)
LMXres <- LMX - Alpha

## plotting the matrix
image.plot(t1,a,t(LMXres),axes="F")
axis(1);axis(2,las=2);box()
mtext("ages",side=2,line=2);mtext("year",side=1,line=2)
dev.off()

## performing the SVD
LCsvd <- svd(LMXres,nu=1,nv=1)

## extract first left- and right-singular vectors of svd
Beta <- c(LCsvd$u)
Kappa1 <- c(LCsvd$v)

## plotting singular vectors
par(mar=c(3,3,1,1))
plot(a,Beta)
plot(t1,Kappa1)

## including the constraints
## constraint 1 
sum.Beta <- sum(Beta)
Beta <- Beta/sum.Beta
sum(Beta)

## constraint 2
Kappa1 <- LCsvd$d[1]*Kappa1*sum.Beta
sum(Kappa1)

## plotting 
plot(a,Beta)
plot(t1,Kappa1)


## second step: adjust kappa to replicate number of deaths in each year
koptim <- function(par,alpha,beta,sum.dx,Exp){
  kappa <- par[1]
  lmx.lc <- alpha+beta*kappa
  z.lc <- exp(lmx.lc)*Exp
  sum.z.lc <- sum(z.lc)
  diff.lc <- abs(sum.dx-sum.z.lc)
  return(diff.lc)
}

Kappa <- numeric(n1)
for (i in 1:n1){
  KappaSecStep <- optimize(f=koptim,interval=c(-100,100),alpha=Alpha,beta=Beta,sum.dx=sum(Y[,i]),
                           Exp=E[,i])
  Kappa[i] <- KappaSecStep$minimum
}
plot(t1,Kappa)

## fitted log-mortality
Ones <- matrix(1,n1)
ETAlc <- Alpha%*%t(Ones) + Beta%*%t(Kappa)

## plotting
matplot(a,ETAlc,t="l",lty=1,lwd=0.8,col=viridis(n1,option = "A"),xlab="",ylab="")
mtext("ages",side=1,line=1.5);mtext("log-mortality",side=2,line=2)

rany <- range(LMX,finite=T)
i <- 1
for(i in 1:n1){
  plot(a, LMX[,i], col=viridis(n1)[i], pch=1, ylim=rany, main=t1[i],xlab="",ylab="")
  mtext("ages",side=1,line=1.5);mtext("log-mortality",side=2,line=2)
  lines(a,ETAlc[,i],col=viridis(n1)[i],lwd=2)
  Sys.sleep(0.1)
}

## fitted life expectancy
e0lc <- apply(exp(ETAlc),2,e0.mx,x=a,sex="F")
plot(t1,e0obs,pch=1,col="blue",xlab="",ylab="",lwd=2)
mtext("years",side=1,line=1.5);mtext("life exp",side=2,line=2)
points(t1,e0lc,pch=4,col="orange",lwd=2)

##---- step 3: forecasting -----

## time series object for Kappa
Kts <- ts(c(Kappa), start = t1[1])

## random walk model with drift (ARIMA (0,1,0) + drift)
modK <- Arima(Kts, order=c(0,1,0), include.drift=TRUE)

## Kappa forecasts
predK <- forecast(modK,h=nF)

## plotting forecasts
plot(predK)

plot(t1,Kappa,t="l",lwd=2,ylim=range(Kappa,predK$lower),xlim=range(t))
lines(tF,predK$mean,col=2,lwd=2)
matlines(tF,predK$upper,col=c(3,4),lwd=2)
matlines(tF,predK$lower,col=c(3,4),lwd=2)

## simulation matrix for future kappas
nS <- 200
SIMk <- matrix(NA,nrow=nF,ncol=nS)

## simulate future K
i <- 2
for(i in 1:nS){
  ## generate simulation with bootsrapping
  kappa.sim <- simulate(modK, nsim=nF,
                        future=TRUE, bootstrap=TRUE)
  # plot(t1,Kappa,ylim=range(Kappa,kappa.sim),xlim=range(t))
  # lines(tF,kappa.sim)
  
  ## derive the bootsrap values
  SIMk[,i] <- kappa.sim
}

## plotting simulations 
plot(t1,Kappa,ylim=range(Kappa,SIMk),xlim=range(t))
matlines(tF,SIMk,col="grey70")

## define PI level
lev <- 95
lev.p <- lev/100

## compute median and PIs from simulations
Kappa.low <- apply(SIMk,1,quantile,prob=(1-lev.p)/2)
Kappa.upp <- apply(SIMk,1,quantile,prob=1-(1-lev.p)/2)
Kappa.med <- apply(SIMk,1,median)

## compare analytical PIs with those of simulations
plot(t1,Kappa,ylim=range(Kappa,SIMk),xlim=range(t))
## median, low and up from RWD
lines(tF,predK$lower[,2],col=3,lwd=2)
lines(tF,predK$upper[,2],col=3,lwd=2)
lines(tF,predK$mean,col=3,lwd=2)
## median, low and up from simulations
lines(tF,Kappa.low,col=4,lwd=2)
lines(tF,Kappa.upp,col=4,lwd=2)
lines(tF,Kappa.med,col=4,lwd=2)

## for each simulation, compute forecast rates and e0

## object to store results
ETAlc.fore <- array(NA,dim = c(m,nF,nS))
e0lc.fore <- matrix(NA,nF,nS)

## loop over simulations
i <- 1
for (i in 1:nS){
  Ones <- matrix(1,nF)
  ETAlc.fore[,,i] <- Alpha%*%t(Ones) + Beta%*%t(SIMk[,i]) 
  e0lc.fore[,i] <- apply(exp(ETAlc.fore[,,i]),2,e0.mx,x=a,sex="F")
}

## derive PIs for e0
e0lc.fore.low <- apply(e0lc.fore,1,quantile,prob=(1-lev.p)/2)
e0lc.fore.upp <- apply(e0lc.fore,1,quantile,prob=1-(1-lev.p)/2)
e0lc.fore.med <- apply(e0lc.fore,1,median)


## derive PIs for ETA
ETAlc.fore.low <- ETAlc.fore.upp <- 
  ETAlc.fore.med <- matrix(NA,m,nF)
i <- 1
for (i in 1:m){
  ETAlc.fore.low[i,] <- apply(ETAlc.fore[i,,],1,quantile,prob=(1-lev.p)/2)
  ETAlc.fore.upp[i,] <- apply(ETAlc.fore[i,,],1,quantile,prob=1-(1-lev.p)/2)
  ETAlc.fore.med[i,] <- apply(ETAlc.fore[i,,],1,median)
}

## plotting forecast rates
my.age <- 30
ylim <- range(LMX[which(a==my.age),],
              ETAlc.fore.low[which(a==my.age),])
plot(t1,LMX[which(a==my.age),],t="n",xlab = "",ylab="",
     ylim=ylim,xlim=range(t),axes=F)
axis(1);axis(2,las=2)
grid();box()
title(paste("age",my.age))
mtext("year",side=1,line=2)
mtext("log-mortality",side=2,las=3,line=1.8)
points(t1,LMX[which(a==my.age),],pch=16,col="grey40",cex=0.75)
lines(t1,ETAlc[which(a==my.age),],col="blue",lwd=2)
matlines(tF,ETAlc.fore[which(a==my.age),,],col="grey70")
lines(tF,ETAlc.fore.med[which(a==my.age),],col="darkblue",lwd=2)
lines(tF,ETAlc.fore.low[which(a==my.age),],col="darkblue",lwd=2,lty=2)
lines(tF,ETAlc.fore.upp[which(a==my.age),],col="darkblue",lwd=2,lty=2)


## plotting life expectancy
ylim <- range(e0obs,e0lc,e0lc.fore)
plot(t1,e0obs,t="n",xlab = "",ylab="",
     ylim=ylim,xlim=range(t),axes=F)
axis(1);axis(2,las=2);grid();box()
title("Life expectancy at birth")
mtext("year",side=1,line=2.2)
mtext(expression(e[0]),side=2,las=1,line=1.9)
points(t1,e0obs,pch=16,col="grey40",cex=0.75)
lines(t1,e0lc,col="blue",lwd=2)
matlines(tF,e0lc.fore,col="grey70")
lines(tF,e0lc.fore.med,col="darkblue",lwd=2)
lines(tF,e0lc.fore.upp,col="darkblue",lwd=2,lty=2)
lines(tF,e0lc.fore.low,col="darkblue",lwd=2,lty=2)

## END
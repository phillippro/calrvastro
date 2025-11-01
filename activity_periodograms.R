Pact <- c()
act.names <- c()
act.data <- list()
logBFs <- c()
plotf <- TRUE
####Sindex
if(Sindex){
    acts <- cbind(tab[iss,1],ss,dss)
    act.data <- c(act.data,list(acts))
    act.name <- 'Sindex'
    periodogram.type <- 'BFP'
    source('activity.R')
    act.names <- c(act.names,act.name)
    Pact <- rbind(Pact,per$Popt[1:Nact.max])
    logBFs <-c(logBFs,per$logBF[which.max(per$power)])
}

####Halpha
if(Halpha){
    acts <- cbind(tab[ihh,1],hh,dhh)
    act.data <- c(act.data,list(acts))
    act.name <- 'Halpha'
    periodogram.type <- 'BFP'
    source('activity.R')
    act.names <- c(act.names,act.name)
    Pact <- rbind(Pact,per$Popt[1:Nact.max])
    logBFs <-c(logBFs,per$logBF[which.max(per$power)])
}

####photon count
if(PhotonCount){
    acts <- cbind(tab[,1],tab[,6],sqrt(tab[,6]))
    act.data <- c(act.data,list(acts))
    act.name <- 'Photon Count'
    periodogram.type <- 'BFP'
    source('activity.R')
    act.names <- c(act.names,act.name)
    Pact <- rbind(Pact,per$Popt[1:Nact.max])
    logBFs <-c(logBFs,per$logBF[which.max(per$power)])
}

###window function
acts <- cbind(tab[,1],rep(1,nrow(tab)),rnorm(nrow(tab),0.1,1e-4))
act.data <- c(act.data,list(acts))
act.name <- 'Window Function'
periodogram.type <- 'LS'
source('activity.R')
act.names <- c(act.names,act.name)
Pact <- rbind(Pact,per$ps[1:Nact.max])
logBFs <-c(logBFs,0)

names(act.data) <-rownames(Pact) <- act.names
colnames(Pact) <- paste(1:Nact.max,'signal')
kappa <- 1

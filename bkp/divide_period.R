pmin <- 1.1#d; min(diff(trv))
#pmin <- 125#d; min(diff(trv))
if(grepl('NA',id)){
    pmax <- 8e4#the upper limit of period in component B
}else{
    pmax <- 2*(tmax-tmin)#adjustable
}
logP.range <- log(pmax)-log(pmin)
cat('pmax=',pmax,'\n')
dlogP <- logP.range/Nbin.per
Pmin <- pmin*exp((nbin.per-1)*dlogP)
Pmax <- pmin*exp(nbin.per*dlogP)
if(!exists('plows') | !exists('Plows')){
    Plows <- plows <- Pmin
    Pups <- pups <- Pmax
}else{
    tmp <- period.division(plows,pups,Pmax,Pmin)
    Pups <- tmp$pup
    Plows <- tmp$plow
}
logPmin = log(Pmin)
logPmax = log(Pmax)
nvmin = 1/Pmax
nvmax = 1/Pmin
######
logPmins <- log(Plows)
logPmaxs <- log(Pups)
nvmins <- 1/Plows
nvmaxs <- 1/Pups
ind0 <- sample(length(Plows))[1]
cat('period.optimize=',period.optimize,'\n')
calcP <- TRUE
if(exists('estimation')){
    if(estimation){
        calcP <- FALSE
    }
}
if(period.optimize){
    ofac <- max(Nbin.per/2,1)
    if(!exists('res.all')){
        res.all <- res
    }
    res.glst <-  glst(t=trv.all,y=res.all,err=eRV.all,ofac=1,fmax=1/Pmin,fmin=1/Pmax)
    Pini <- res.glst$P[ind[which.max(res.glst$power[ind])]]
    Kini <- res.glst$Kopt
}else{
    Pini <- runif(1,Plows[ind0],Pups[ind0])
}
#Pini <- 131
#Kini <- 4.26
cat('Pini=',Pini,'\n')
if(FALSE){
    cat('logPmins=',logPmins,'\n')
    cat('logPmaxs=',logPmaxs,'\n')
    cat('Plows=',Plows,'\n')
    cat('Pups=',Pups,'\n')
    cat('plows=',plows,'\n')
    cat('pups=',pups,'\n')
    cat('Pmax=',Pmax,'\n')
    cat('Pmin=',Pmin,'\n')
}
if(period.par=='nv'){
    per.ini <- 1/Pini
    per.max <- nvmax
    per.min <- nvmin
}else if(period.par=='P'){
    per.ini <- Pini
    per.max <- Pmax
    per.min <- Pmin
}else{
    per.ini <- log(Pini)
    per.max <- logPmax
    per.min <- logPmin
}

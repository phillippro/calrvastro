source('mcmc_func.R')
options(scipen=9)
yr2d <- 365.25
show.digit <- function(x,Ndig,type='format'){
    if(type=='format'){
        format(round(x,Ndig),nsmall=Ndig)
    }else if(type=='signif'){
        signif(x,Ndig)
    }
}
#    load(paste0('phase',version,'.Robj'))
#    load(paste0('phase',version,'.Robj'))
load('phase1.Robj')
out.phase1 <- out.phase
load('phase2.Robj')
out.phase2 <- out.phase
out.phase <- c(out.phase1,out.phase2)

fs <- names(out.phase)
#pars <- read.table('../../data/UVES/uves_mass.txt',header=TRUE)
planets <- read.csv('../../data/code/SigClassify/ranking_complex1.csv')
targets <-  gsub('.+/|_.+','',fs)
n0 <- gsub(' ','',as.character(planets[,'Name']))
n1 <- gsub(' ','',as.character(planets[,'ID']))
n2 <- gsub(' ','',as.character(planets[,'StarKnown']))
inds <- c()
for(s in targets){
    s1 <- gsub('GJ','GL',s)
    ind <- which(n0==s | n1==s | n2==s | n0==s1 | n1==s1 | n2==s1)
    if(length(ind)==0){
        if(s=='GJ9349') ind <- which(n0=='HIP54532')
        if(s=='GJ251') ind <- which(n0=='HD265866')
        if(s=='GJ880') ind <- which(n0=='HD216899')
        if(s=='GJ4070') ind <- which(n0=='HIP91699')
        if(s=='GJ2056') ind <- which(n0=='HIP34785')
        if(s=='GJ9066') ind <- which(n0=='GL83.1')
        if(s=='GJ480') ind <- which(n0=='HIP61706')
        if(s=='GJ3072') ind <- which(n0=='HIP4845')
    }
#    cat(s,'\n')
#    cat(s1,'\n')
#    cat('ind=',ind,'\n')
    inds <- c(inds,ind)
}
#stop()
tab <- planets[inds,]
#inds <- sort(as.numeric(gsub('GJ','',n2)),index.return=TRUE)$ix
#targets <- target[inds]

Ns <- 1e5
#value.type <- 'mq'#map+quantile
#value.type <- 'ms'#mean+sd
value.type <- c('ms','mq')#mean+sd and map+quantile
type <- 'format'
xup <- 'x99per'
xlow <- 'x1per'

table.out <- c()
stars <- c()
Ndig0 <- Ndig <- 2
pp <-c('b','c','d','e','f','g','h','i')
ind.out <- c()
for(k5 in 1:length(fs)){
    f <- fs[k5]
    star <- target <- targets[k5]
    Nsig <- out.phase[[f]]$Nsig
    par.stat <- out.phase[[f]]$par.stat[[paste0('sig',Nsig)]]
    par.kep <- out.phase[[f]]$par.opt
    Mstar <- as.numeric(as.character(tab[k5,'Ms']))
    eMstar <- Mstar*0.1
#    eMstar <-as.numeric(as.character(pars[ind2,'e_StellarMass']))
    ind <-grep('per',colnames(par.stat))
    inds <- sort(par.stat['mean',ind],index.return=TRUE)$ix
    for(i in 1:length(inds)){
        for(vt in value.type){
            par.kep0 <- par.kep
            par.stat0 <- par.stat
#            ind.out <- c(ind.out,ind2)
            j <- inds[i]
            ind.kep <- (j-1)*5+(1:5)
            if(vt=='ms'){
#                ns <- gsub('HD','HD ',target)
                ns <- gsub('HD','HD ',star)
                ns <- gsub('HIP','HIP ',ns)
                ns <- gsub('GL','GL ',ns)
                ns <- gsub('GJ','GJ ',ns)
                tmp <-paste(ns,pp[i])
                stars <- c(stars,tmp)
            }else{
                tmp <-''
            }
            for(k in ind.kep){
                if(k%%5==1){
                    ps <- exp(rnorm(Ns,par.stat['mean',k],par.stat['sd',k]))
                    val <- data.distr(ps,plot=FALSE)
                    if(length(val)>nrow(par.stat)){
                        par.stat[,k] <- val[-1]
                    }else{
                        par.stat[,k] <- val
                    }
                    popt <- par.kep[k] <- exp(par.kep[k])
                    if(par.stat['sd',k]<1e-3){
                        Ndig <- max(3,ceiling(log10(1/par.stat['sd',k])))
                    }else{
                        Ndig <- 3
                    }
                }else if(k%%5==4 | k%%5==0){
                    Ndig <- 0
                }else{
                    Ndig <- Ndig0
                }
#cat('k=',k,'\n')
#cat('Ndig=',Ndig,'\n\n')
                if(k%%5==2){
                    ks <- rnorm(Ns,par.stat['mean',k],par.stat['sd',k])
                    kopt <- par.kep[k]
                }
                if(k%%5==3){
                    es <-rnorm(5*Ns,par.stat['mean',k],par.stat['sd',k])
                    es <- es[es>0 & es<1][1:Ns]
                    eopt <- par.kep[k]
                }
                if(k%%5==4|k%%5==0){
                    par.stat[,k] <-par.stat[,k]*180/pi
                    par.kep[k] <- par.kep[k]*180/pi
                }
                if(vt=='ms'){
                    str <- paste0('$',show.digit(par.stat['mean',k],Ndig,type),'\\pm',show.digit(par.stat['sd',k],Ndig,type),'$')
                }else{
                                        #                    str <- paste0('$',show.digit(par.kep[k],Ndig,type),'[',show.digit(par.stat[xlow,k],Ndig,type),',',show.digit(par.stat[xup,k],Ndig,type),']$')
                    dpar.lower <- -par.stat[xlow,k]+par.kep[k]
                    dpar.upper <-par.stat[xup,k]-par.kep[k]
                    str <- paste0('$',show.digit(par.kep[k],Ndig,type),'_{-',show.digit(dpar.lower,Ndig,type),'}^{+',show.digit(dpar.upper,Ndig,type),'}$')
                }
                if(k%%5==0) str <- paste0(str,'&\\\\')
                tmp <- c(tmp,str)
            }
###calculate semi-major axis
###calculate mass
            ms <-rnorm(5*Ns,Mstar,eMstar)
            es <- es[ms>0][1:Ns]
            ma.opt <- K2msini(kopt,popt,eopt,Mstar)
            ma <- K2msini(ks,ps,es,ms)
            msini <- ma$me
            a <-ma$a
            msini.stat <- data.distr(msini,plot=FALSE)
            a.stat <- data.distr(a,plot=FALSE)
            if(vt=='ms'){
                str.m <- paste0('$',show.digit(msini.stat['mean'],Nsig),'\\pm',show.digit(msini.stat['sd'],Nsig),'$')
                str.a <- paste0('$',show.digit(a.stat['mean'],Ndig0+1),'\\pm',show.digit(a.stat['sd'],Ndig0+1),'$')
            }else{
#                str.m <-  paste0('$',show.digit(ma.opt$me,Ndig,type),'[',show.digit(msini.stat[xlow],Ndig,type),',',show.digit(msini.stat[xup],Ndig,type),']$')
                dm.lower <- -msini.stat[xlow]+ma.opt$me
                dm.upper <- msini.stat[xup]-ma.opt$me
                str.m <- paste0('$',show.digit(ma.opt$me,Ndig0,type),'_{-',show.digit(dm.lower,Ndig0,type),'}^{+',show.digit(dm.upper,Ndig0,type),'}$')
#                str.a <-  paste0('$',show.digit(ma.opt$a,Ndig0+1,type),'[',show.digit(a.stat['x10per'],Ndig0+1,type),',',show.digit(a.stat[xup],Ndig0+1,type),']$')
                da.lower <- -a.stat[xlow]+ma.opt$a
                da.upper <- a.stat[xup]-ma.opt$a
                str.a <- paste0('$',show.digit(ma.opt$a,Ndig0+1,type),'_{-',show.digit(da.lower,Ndig0+1,type),'}^{+',show.digit(da.upper,Ndig0+1,type),'}$')
            }
            tmp <- c(tmp[1],str.m,str.a,tmp[-1])
#            cat(tmp,'\n')
            table.out <- rbind(table.out,tmp)
            par.stat <- par.stat0
            par.kep <- par.kep0
        }
    }
}
fout <- 'table_planet.txt'
cat(fout,'\n')
write.table(table.out,file=fout,quote=FALSE,sep='&',row.names=FALSE,col.names=FALSE)


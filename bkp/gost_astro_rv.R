library(magicaxis)
library(graphics)
library(RColorBrewer)
source('mcmc_func.R')
options(scipen = 0)
library(paletteer)
library("plotrix")
source('OrbitFunction.R')
source('timing_function.R')
source('general_function.R')
source('sofa_function.R')
set.seed(9999)
args <- commandArgs(trailingOnly=TRUE)
if(length(args)>0){
    star <- as.character(args[1])
    cal <- as.numeric(args[-1])
}else{
#    star <- 'HD22049'
#   star <- 'GJ317'
#   star <- 'HD164604'
   star <- 'UCAC4569-026385'
#   star <- 'HD209100'
#   star <- 'HD222237'
#    cal <- c(2025,1,1)
    cal <- c(2024,1,1)
}
cal <- cbind(t(cal))
comptype <- 'companion'
Nmc <- Nsamp <- 0
#Nmc <- Nsamp <- 1000
m22 <- read.table('../data/combined/mamajek22.txt',header=TRUE)
ms <- as.numeric(m22[,'Msun'])
mg <- as.numeric(m22[,'M_G'])
ind <- which(!is.na(ms) & !is.na(mg))
mrl.m22 <- approxfun(ms[ind],mg[ind])
mlow.m22 <- min(ms[ind])
mup.m22 <- max(ms[ind])
version <- 1
set.seed(100)
ff <- read.table('1companion_files.txt')[,1]
#ff <- read.table('bh_companion.txt')[,1]
stars <- gsub('\\/.+','',ff)
if(star!=''){
    ind <- which(stars==star)
}else{
                                        #   ind <- 1:length(stars)
    ind <- 1:10
}
if(length(ind)==0){
    stop('No file is found for this target!')
}
fpdf <- ff[ind]
fobj <- paste0('results/',gsub('.pdf','.Robj',fpdf))
###modifiy input file names
                                        #Nsig <- length(fs)
planets <- read.csv('../data/code/SigClassify/ranking_complex2.csv')
n0 <- gsub(' ','',as.character(planets[,'Name']))
n1 <- gsub(' ','',as.character(planets[,'ID']))
n2 <- gsub(' ','',as.character(planets[,'StarKnown']))
###collect all data needed for the plot
cat('load ',fobj,'\n')
if(!exists('dP')) dP <- 0
gdr1plx <- 'raw'
#if(!exists('out')) load(fobj)
load(fobj)
if(!any(names(out)=='cov.astro')) out$cov.astro <- cov.astro
if(target=='UCAC4569-026385'){
#if(FALSE){
#    tadd <- c(2460371.00370,2460371.01898,2460371.03425)%%24e5-tmin%%24e5
    tadd <- c(2460371.00370,2460371.01898,2460371.03425)
    i0 <- which(out$lamost$RV[,3]==390)[1]
    rv0 <- 1630-out$lamost$RV[i0,2]#the RV0 subtracted from the raw data
    rvadd <- c(-27.31,-27.31,-27.31)*1e3-rv0
    ervadd <- c(470,520,520)
    t <- c(out$lamost$RV[,1],tadd)
    rv <- c(out$lamost$RV[,2],rvadd)
    erv <- c(out$lamost$RV[,3],ervadd)
    out$lamost$RV <- cbind(t,rv,erv)
    out$all <- data.frame(trv.all=t,RV.all=rv,eRV.all=erv,ins.all=rep('lamost',length(t)))
}

#load(fobj)
#load(fobj)
if(!exists('eta0')) eta0 <- rep(NA,10)
source('mcmc_func.R')
if(is.null(out$comp.epoch)) out$comp.epoch <- 1
if(is.null(out$cats)) out$cats <- out$gdrs
if(is.null(out$cat.ind)) out$cat.ind <- out$gdr.ind
if(is.null(out$ins.rv)) out$ins.rv <- out$ins
if(max(tsim)-min(tsim)<max(Popt)){
    tsim <- min(tmin,out$astrometry[,1])+seq(0,max(Popt),length.out=10000)
}
ins <- out$ins.rv
nqp <- list()
for(k in ins){
    nqp[[k]] <- out[[k]]$noise$nqp
}
ell.col <- c('grey','orange')
ell.pch <- 20
rel.pch <- 15
cc <- brewer.pal.info[brewer.pal.info$category=='qual',]
cc1 <- unlist(mapply(brewer.pal,cc$maxcolors,rownames(cc)))
cc2 <- cc1[c(1:3,16,5:7,9:12,14:15)]
alpha <- 80
size <- 1
DT <- 1
if(!is.null(out$data.epoch)){
    Npanel <- 4
    margin <- 1.2
    oma <- c(2,6,0,1)
#    oma <- c(2,7,0,1)
#    mgp <- c(2,1,0)
    mgp <- c(2,0.5,0)
}else{
    Npanel <- 3
    margin <- 1.2
    oma <- c(2,6,0,1)
    mgp <- c(2,1,0)
#    oma <- c(4,6,3,1)
}
mm <- c(1.2*margin,3*margin,2*margin,0)
pch.size <- 1
line.size <- 2
fit.col <- tcol('red',90)
fit.opt <- tcol(c('red','red'),20)
Nstar <- length(stars)

if(grepl('\\+_',fobj)){
    imaging <- TRUE
}else{
    imaging <- FALSE
}
dd <- './'
if(star=='HD209100'| star=='HD22049'){
 dd <- '/Users/ffeng/Documents/projects/dwarfs/paper/nearest_jupiters/arxiv/'
}
if(star=='UCAC4569-026385') dd <- '/Users/ffeng/Documents/projects/dwarfs/paper/BH/'
if(star=='HD222237') dd <- '/Users/ffeng/Documents/projects/dwarfs/paper/HD222237/'

#fpdf <- paste0(comptype,'_plot/',star,'_N',Nmc,'_fit.pdf')
fpdf <- paste0(dd,star,'_N',Nmc,'_fit_v2.pdf')
cat('output:\n',fpdf,'\n')
pdf(fpdf,4*Npanel,4)
ss0 <- split.screen( figs = c(1,Npanel))
par(oma=oma,cex=size,cex.lab=1.2*size,cex.axis=size,mgp=mgp,cex.main=1.2*size)
ss1 <- split.screen(rbind(c(0,1,0.3,1),c(0,1,0,0.3)),screen=ss0[1])
par(mar=c(0,0,0,0),oma=oma)
Nw0 <- Nset <- length(out$ins.rv)
pchs <- rep(20,Nw0)
if(Nw0>10){
    pchs[1:round(Nw0/2)] <- 15
}
###split the RV plot into two sub-panel
instr <- ins
###simulation for ylim
jitter <- TRUE
parallax <- FALSE
if(!exists('bases')) bases <- rep('natural',5)
astrometry <- 8
kep.type <- 'pure'
yr2d <- 365.25
period.par <- 'logP'
offset <- TRUE
time.unit <- 1
Np <- length(ind)
prior.type <- 'mt'
pars <- par.opt
trv.all <- tall <- out$all$trv.all
if(star=='UCAC4569-026385'){
    out$tiall <- data.frame( bjd=tall,  instr='lamost', type='rv')
    out$ind.all$rv$lamost <- c(out$ind.all$rv$lamost,length(out$ind.all$rv$lamost)+(1:3))
}
ind <- grep('per',names(pars))
Popt <- exp(par.opt[ind])
indI <- grep('Inc',names(pars))
par.opt[indI] <- par.opt[indI]%%pi
t3 <- seq(min(tall),max(tall),length.out=1e4)
if(sd(out$all$RV.all)<1e3){
    rv.unit <- 'm/s'
}else{
    rv.unit <- 'km/s'
}
rvs <- RV.kepler(pars.kep=par.opt,tt=t3,kep.only=TRUE,nqp=nqp)$rv
ind <- sort(out$all$trv.all,index.return=TRUE)$ix
erv.all <- unlist(sapply(1:Nset,function(k) out[[ins[k]]]$RV[,3]))[ind]
RV2.all <- unlist(sapply(out$ins.rv,function(ins) loglikelihood(par.opt,bases=bases,prediction=TRUE,nqp=nqp)$res[[ins]]$res1+RV.kepler(par.opt,kep.only=TRUE,bases=bases,nqp=nqp)$rv[[ins]]))

#bins <- wtb(trv.all-tmin,RV2.all,erv.all,dt=100)
bins <- wtb(out$all$trv.all-tmin,RV2.all,out$all$eRV.all,dt=DT)
res.bins <- wtb(out$all$trv.all-tmin,out$res.comb.all[[paste0('sig',Nsig)]],out$all$eRV.all,dt=DT)
####binning for individual data set
res.list <-  loglikelihood(par.opt,bases=bases,prediction=TRUE,nqp=nqp)$res
resbin <- rvbin <- list()
RV2.all <- unlist(sapply(instr,function(i) res.list[[i]]$res1+RV.kepler(par.opt,kep.only=TRUE,bases=bases,nqp=nqp)$rv[[i]]))
if(rv.unit=='km/s'){
    rvs <- rvs/1e3
    RV2.all <- RV2.all/1e3
}
for(k in 1:Nset){
#    kepres <- res.list[[ins[k]]]$res2
    res <- wtb(out[[ins[k]]]$RV[,1]%%24e5-tmin%%24e5,res.list[[ins[k]]]$res2,out[[ins[k]]]$RV[,3],dt=DT)
    resbin[[k]] <- res
    rv <- RV.kepler(par.opt,tt=res[,1]+tmin,kep.only=TRUE,bases=rep('natural',10),nqp=nqp)$rv
    res[,2] <- rv+res[,2]
    rvbin[[k]] <- res
}
ylim <- range(RV2.all,rvs)
###first plot
screen(ss1[1])
par(mar=c(0,mm[-1]))
for(k in 1:Nset){
    t <- out[[ins[k]]]$RV[,1]%%24e5-tmin%%24e5
    rv <- res.list[[ins[k]]]$res2+RV.kepler(par.opt,kep.only=TRUE,bases=bases,nqp=nqp)$rv[[ins[k]]]
    erv <- out[[ins[k]]]$RV[,3]
    if(rv.unit=='km/s'){
        rv <- rv/1e3
        erv <- erv/1e3
    }
#    c1 <- cc2[k]
    c1 <- 'grey'
    if(k==1){
        plot(t,rv,xlab=paste0('BJD-',round(tmin%%24e5+24e5,3)),ylab='',col=tcol(c1,alpha),pch=pchs[k],xlim=c(min(tall%%24e5-tmin%%24e5),max(tall%%24e5-tmin%%24e5)),ylim=ylim,xaxt='n',yaxt='n',cex=1.2*pch.size,main='RV')#paste0(star,' RV')
#        magaxis(2,cex=size,majorn=3)
        magaxis(2,cex=size)
    }else{
        points(t,rv,pch=pchs[k],col=tcol(c1,alpha),cex=pch.size)
    }
#    arrows(t,rv-erv,t,rv+erv,length=0.05,angle=90,code=3,col=tcol(c1,alpha))
###binning
#    bins <- wtb(t,rv,erv,dt=DT)
    bins <- rvbin[[k]]
    if(rv.unit=='km/s'){
        bins[,2:3] <- bins[,2:3]*1e-3
    }
    if(star=='UCAC4569-026385'){
        points(t,rv,pch=pchs[k],col=cc2[k],cex=pch.size)
        arrows(t,rv-erv,t,rv+erv,length=0.05,angle=90,code=3,col=cc2[k])
#rvs <- RV.kepler(pars.kep=par.opt,tt=tt,kep.only=TRUE,nqp=nqp)$rv
#if(rv.unit=='km/s'){
#    rvs <- rvs/1e3
#}
#lines(tt%%24e5-tmin%%24e5,rvs,col=fit.opt[1],lwd=line.size)$rv
    }else{
        points(bins[,1],bins[,2],col=cc2[k],pch=pchs[k])
        arrows(bins[,1],bins[,2]-bins[,3],bins[,1],bins[,2]+bins[,3],length=0.05,angle=90,code=3,col=cc2[k])
    }
}
#points(bins[,1],bins[,2],col='black',pch=20)
#arrows(bins[,1],bins[,2]-bins[,3],bins[,1],bins[,2]+bins[,3],length=0.05,angle=90,code=3,col='black')

mtext(side=2,text=as.expression(paste0('RV [',rv.unit,']')),line=1.5,xpd=NA,cex=size)
ind <- which(!grepl('HARPS',instr))
if(length(ind)>0) instr[ind] <- toupper(instr[ind])
ind <- grep('UCLES',instr)
if(length(ind)>0) instr[ind] <- 'AAT'
ind <- grep('SOPHIE2',instr)
if(length(ind)>0) instr[ind] <- 'SOPHIEpost'
ind <- grep('SOPHIE1',instr)
if(length(ind)>0) instr[ind] <- 'SOPHIEpre'
ind <- grep('SOPHIEPRE',instr)
if(length(ind)>0) instr[ind] <- 'SOPHIEpre'
ind <- grep('SOPHIEPOST',instr)
if(length(ind)>0) instr[ind] <- 'SOPHIEpost'
ind <- grep('HARPSPRE',instr)
if(length(ind)>0) instr[ind] <- 'HARPSpre'
ind <- grep('HARPSPOST',instr)
if(length(ind)>0) instr[ind] <- 'HARPSpost'
ind <- grep('Subaru',instr)
if(length(ind)>0) instr[ind] <- 'HDS'
ind <- grep('MCD',instr)
if(length(ind)>0) instr[ind] <- 'HRS'
ind <- grep('SUBARU',instr)
if(length(ind)>0) instr[ind] <- 'HDS'
ind <- grep('C07|C14|C98',instr)
if(length(ind)>0) instr[ind] <- gsub('C','COR',instr[ind])
ind <- grep('CORALIE07|CORALIE14|CORALIE98',instr)
if(length(ind)>0) instr[ind] <- gsub('CORALIE','COR',instr[ind])
instr[instr=='L20'] <- 'AL19'
if(exists('star.name')){
    sn <- gsub('GL','GJ',star.name)
}else{
    sn <- gsub('GL','GJ',star)
}
if(!imaging){
    inset <- c(-0.6,0)
}else{
#    inset <- c(-0.82,-0.1)
    inset <- c(-0.7,-0.1)
}
                                        #        legend('topleft',legend=c(instr,'Hipparcos','Gaia','Best fit'),col=c(cols[1:Nw0],ell.col,fit.opt,'black'),lty=c(rep(NA,Nw0+2),1,NA),pch=c(rep(20,Nw0),ell.pch,NA,3),xpd=NA,inset=inset,bty='n',title=sn,lwd=c(rep(NA,Nw0+2),line.size,NA),pt.cex=c(rep(1,Nw0),rep(1,2),NA,1),y.intersp=1.1)
font.size <- 1.2
size <- 0.8
if(star=="UCAC4569-026385"){
    sn <- 'G3425'
}

if(!imaging){
                                        #            legend('topleft',legend=c(instr,'Hipparcos IAD','Gaia GOST','Best fit'),col=c(cols[1:Nw0],ell.col,fit.opt[1],'black'),lty=c(rep(NA,Nw0+2),1),pch=c(rep(20,Nw0),rep(ell.pch,2),NA,3),xpd=NA,inset=inset,bty='n',lwd=c(rep(NA,Nw0+2),line.size,NA),pt.cex=c(rep(font.size,Nw0),rep(font.size,2),NA,font.size),y.intersp=1.1,cex=size)#title=sn
    instr.show <- instr
    instr.show[instr.show=='APF'] <- 'APFp'
    instr.show[instr.show=='APFJ'] <- 'APFh'
    instr.show[instr.show=='KECK'] <- 'HIRESp'
    instr.show[instr.show=='KECKJ'] <- 'HIRESh'
    legend('topleft',legend=c(instr.show,'Best fit'),col=c(cc2[1:Nw0],fit.opt[1],'black'),lty=c(rep(NA,Nw0),1),pch=c(pchs,NA),xpd=NA,inset=inset,bty='n',lwd=c(rep(NA,Nw0),line.size),pt.cex=c(rep(font.size,Nw0),NA),y.intersp=1.1,cex=size,title=sn)
}else{
    ins0 <- unique(unlist(sapply(1:Nsig,function(n) names(out$rel[[n]]))))
    inss <- toupper(ins0)
    if(sn=='HD42581' & any(inss=='ASTROMETRY')) inss[inss=='ASTROMETRY'] <- 'MB21'
    if(sn=='HD182488' & any(inss=='ASTROMETRY')) inss[inss=='ASTROMETRY'] <- 'BB18'
    if(sn=='HD4113' & any(inss=='SPHERE')) inss[inss=='SPHERE'] <- 'AC18'
    if(sn=='GJ3677' & any(inss=='E10')) inss[inss=='E10'] <- 'DE10'
    if(sn=='HD13724' & any(inss=='R20')) inss[inss=='R20'] <- 'ER20'
    if(sn=='HD7449' & any(inss=='R16')) inss[inss=='R16'] <- 'TR16'
    if(sn=='GJ494' & any(inss=='M19')) inss[inss=='M19'] <- 'AM19'
    if(sn=='HD39060' & any(inss=='G20')) inss[inss=='G20'] <- 'SL21'
    if(sn=='HD39060' & any(inss=='HCI')) inss[inss=='HCI'] <- 'AL20'
    Nins <- length(inss)
    col.rel <- cc[1:Nins]
    names(col.rel) <- ins0
                                        #            Nrel <- length(out$rel)
    Nrel <- 1
#    legend('topleft',legend=c(instr,'Hipparcos IAD','Gaia GOST',inss,rep('Best fit',Nrel)),col=c(cc2[1:Nw0],ell.col,col.rel,fit.opt[1:Nrel],'black'),lty=c(rep(NA,Nw0+2+Nins),rep(1,Nrel),NA),pch=c(pchs,rep(ell.pch,2),rep(rel.pch,Nins),rep(NA,Nrel),3),xpd=NA,inset=inset,bty='n',lwd=c(rep(NA,Nw0+2+Nins),rep(line.size,Nrel),NA),pt.cex=c(rep(1,Nw0+2+Nins),rep(NA,Nrel),1),y.intersp=1.1,title=sn)
    legend('topleft',legend=c(instr,'Hipparcos IAD','Gaia GOST',inss,rep('Best fit',Nrel)),col=c(cc2[1:Nw0],ell.col,col.rel),lty=c(rep(NA,Nw0+2+Nins),rep(1,Nrel)),pch=c(pchs,rep(ell.pch,2),rep(rel.pch,Nins),rep(NA,Nrel)),xpd=NA,inset=inset,bty='n',lwd=c(rep(NA,Nw0+2+Nins),rep(line.size,Nrel)),pt.cex=c(rep(1,Nw0+2+Nins),rep(NA,Nrel)),y.intersp=1.1,title=sn)
}
tspan <- max(tall)-min(tall)
tt <- seq(min(tall),max(tall),length.out=1e4)
                                        #    if(Nmc>0){
if(FALSE){
    t1 <- min(min(tall),mean(tall)-max(Popt)/2)
    t2 <- max(max(tall),mean(tall)+max(Popt)/2)
    tt <- seq(t1,t2,length.out=1e4)
}
##
if(Nmc>0){
    mc <- out$mcmc.opt[[paste0('sig',Nsig)]]
    ind.sample <- inds <- sample(1:nrow(mc),Nsamp)
    Npar <- length(par.opt)
    for(j3 in inds){
        rvs <- RV.kepler(pars.kep=mc[j3,1:Npar],tt=tt,kep.only=TRUE,nqp=nqp)$rv
        if(rv.unit=='km/s'){
            rvs <- rvs/1e3
        }
        lines(tt%%24e5-tmin%%24e5,rvs,col=fit.col,lwd=2)
    }
}
rvs <- RV.kepler(pars.kep=par.opt,tt=tt,kep.only=TRUE,nqp=nqp)$rv
if(rv.unit=='km/s'){
    rvs <- rvs/1e3
}
lines(tt%%24e5-tmin%%24e5,rvs,col=fit.opt[1],lwd=line.size)$rv
#stop()
###residual plot
screen(ss1[2])
par(mar=c(mm[1:2],0,mm[4]))
rv.planet <- RV.kepler(pars.kep=par.opt,tt=trv.all,kep.only=TRUE,nqp=nqp)$rv#plus offset
#drv <- unlist(sapply(out$ins.rv,function(i) res.list[[i]]$res2))
drv <- unlist(sapply(out$ins.rv,function(i) res.list[[i]]$res1))
ylim <- range(-3*sd(drv),3*sd(drv))
xlim <- range(tall)%%24e5-tmin%%24e5
#ylim <- range(unlist(sapply(1:Nset,function(j2) resbin[[j2]])))
for(j2 in 1:Nset){
    x <- out[[ins[j2]]]$RV[,1]%%24e5-tmin%%24e5
    y <- res.list[[ins[j2]]]$res2
#    y <- res.list[[ins[j2]]]$res1
    ey <- out[[ins[j2]]]$RV[,3]
    if(star=="UCAC4569-026385"){
        y <- y/1e3
        ey <- ey/1e3
        ylim <- ylim/1e3
    }
    t <- out[[ins[j2]]]$RV[,3]
#    c1 <- cc2[j2]
    c1 <- 'grey'
    if(j2==1){
#        plot(x,y,xlab='',ylab='',ylim=ylim,xlim=xlim,col=tcol(cc2[j2],alpha),pch=pchs[j2],cex=size,xaxt='n',yaxt='n')
        plot(x,y,xlab='',ylab='',ylim=ylim,xlim=xlim,col=tcol(c1,alpha),pch=pchs[j2],cex=size,xaxt='n',yaxt='n')
        magaxis(1:2,cex=size,majorn=3)
    }else{
#        points(x,y,col=tcol(cc2[j2],alpha),pch=20,cex=1)
        points(x,y,col=tcol(c1,alpha),pch=pchs[j2],cex=1)
    }
#    arrows(x,y-ey,x,y+ey,length=0.05,angle=90,code=3,col=tcol(c1,alpha))
###binning
#    res.bins <- wtb(x,y,ey,dt=DT)
    res.bins <- resbin[[j2]]
    if(star=="UCAC4569-026385") res.bins[,2:3] <- res.bins[,2:3]/1e3#m/s->km/s
    if(star=='UCAC4569-026385'){
        points(x,y,col=cc2[j2],pch=pchs[j2])
        arrows(x,y-ey,x,y+ey,length=0.05,angle=90,code=3,col=cc2[j2])
    }else{
        points(res.bins[,1],res.bins[,2],col=cc2[j2],pch=pchs[j2])
        arrows(res.bins[,1],res.bins[,2]-res.bins[,3],res.bins[,1],res.bins[,2]+res.bins[,3],length=0.05,angle=90,code=3,col=cc2[j2])
    }
}
#points(res.bins[,1],res.bins[,2],col='black',pch=20)
#arrows(res.bins[,1],res.bins[,2]-res.bins[,3],res.bins[,1],res.bins[,2]+res.bins[,3],length=0.05,angle=90,code=3,col='black')
                                        #    par(mgp=c(3,1,0),xpd=NA)
mtext(side=1,text=paste0('JD-',round(tmin%%24e5+24e5,3)),xpd=NA,line=2,cex=1.2*size)
mtext(side=2,text=as.expression(paste0('O-C [',rv.unit,']')),line=1.5,xpd=NA,cex=1.2*size)
abline(h=0,lty=2,lwd=2)
####astrometry
#dev.off()
#stop()
for(j5 in 1:Nsig){
#   source('paper_ellipse_plot.R')
    source('paper_ellipse_plot_general.R')
}
close.screen(all = TRUE)
dev.off()




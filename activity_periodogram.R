library(magicaxis)
source('periodograms.R')
source('periodoframe.R')
#load('../output/LHS1140/keppure_priormt_poly11_Ndata144_quantifyTRUE_1per1_Nw1_LHS1140_ind0_2planet_ARMA01_Nsamp10000000_tem1_acc0.43_pretem1P92.1d24.6d_negLmax398_optpar.Robj')
#load('../output/LHS1140/keppure_priormt_poly11_Ndata144_quantifyTRUE_1per1_Nw1_LHS1140_ind0_2planet_GPqp_Nsamp4000000_tem1_acc0.67_pretem1P93.1d24.5d_negLmax400_optpar.Robj')
load('../output/LHS1140/keppure_priormt_poly11_Ndata144_quantifyTRUE_1per1_Nw1_LHS1140_ind0_3planet_ARMA00_Nsamp10000000_tem1_acc0.26_pretem1P91.4d24.6d3.8d_negLmax390_optpar.Robj')
Popt <- exp(par.stat['opt',grepl('per',colnames(par.stat))])
Prot <- 131

#tab <- read.table('../data/aperture/LHS1140/LHS1140.dat',header=TRUE)
tab <- read.table('../data/aperture/LHS1140/LHS1140_HARPS.dat',header=TRUE,check.names=FALSE)
etab <- read.table('../data/aperture/LHS1140/LHS1140_HARPS_error.dat',header=TRUE,check.names=FALSE)
t <- tab[,1]
RV <- tab[,2]
eRV <- tab[,3]
proxies <- tab[,3+1:7]
eproxies <- etab[,1:7]
#Sindex <- tab[,'S-index']
#Halpha <- tab[,'Halpha']
#FWHM <- tab[,'FWHM']
#BIS <- tab[,'BIS']
#names <- c('RV',colnames(proxies))
names <- c(colnames(proxies),'Window Function')
names[7] <- 'Differential RV'
#per.type <- 'glst'
per.type <- 'glst'
cat('per.type=',per.type,'\n')
Nd <- 1
lwd <- 3
FAP <- FALSE
ofac <- 10
indicators <- c(1,1)
fmax <- 1/1.2
#fmax <- 1/10
#fmax <- 1/70
Dt <- max(t)-min(t)
#Dt <- 200
dy <- 1.1
###############loop
#for(k1 in rev(1:Nd)){
######
cat('Periodogram for RVs\n')
separate <- FALSE
Nma <- 0
if(separate){
if(per.type=='glst'){
    rv.per <- glst(t=t,y=RV,err=eRV,ofac=ofac,fmax=fmax,fmin=1/Dt)
    cat('Periodogram for Sindex\n')
    Sindex.per <- glst(t=t,y=Sindex,err=rnorm(length(t),1,0.01),ofac=ofac,fmax=fmax,fmin=1/Dt)
#    cat('Periodogram for BIS\n')
    BIS.per <- glst(t=t,y=BIS,err=rnorm(length(t),1,0.01),ofac=ofac,fmax=fmax,fmin=1/Dt)
    cat('Periodogram for FWHM\n')
    FWHM.per <- glst(t=t,y=FWHM,err=rnorm(length(t),1,0.01),ofac=ofac,fmax=fmax,fmin=1/Dt)
    Halpha.per <- glst(t=t,y=Halpha,err=rnorm(length(t),1,0.01),ofac=ofac,fmax=fmax,fmin=1/Dt)
    if(FALSE){
        cat('Periodogram for window function\n')
        win.per <- lsp(times=t,x=rep(1,length(t)),ofac=ofac,fmax=fmax,fmin=1/Dt)
    }
    Nma <- 0
}else{
    Inds <- 0
    Indices <- NA
    Nma <- 1
    rv.per <- BFP(t=t,y=RV,dy=eRV,ofac=ofac,fmax=fmax,fmin=1/Dt,Indices=Indices,progress=FALSE,Inds=Inds,Nma=Nma)
    cat('Periodogram for Sindex\n')
    Sindex.per <- BFP(t=t,y=scale(Sindex),dy=rep(0.1,length(t)),ofac=ofac,fmax=fmax,fmin=1/Dt,progress=FALSE,Inds=Inds,Nma=Nma)
    cat('Periodogram for BIS\n')
    BIS.per <- BFP(t=t,y=scale(BIS),dy=rep(0.1,length(t)),ofac=ofac,fmax=fmax,fmin=1/Dt,progress=FALSE,Inds=Inds,Nma=Nma)
    cat('Periodogram for FWHM\n')
    FWHM.per <- BFP(t=t,y=scale(FWHM),dy=rep(0.1,length(t)),ofac=ofac,fmax=fmax,fmin=1/Dt,progress=FALSE,Inds=Inds,Nma=Nma)
    Halpha.per <- BFP(t=t,y=scale(Halpha),dy=rep(0.1,length(t)),ofac=ofac,fmax=fmax,fmin=1/Dt,progress=FALSE,Inds=Inds,Nma=Nma)
}
}
#####RVs
#paste0('GLST; p value:',format(rv.glst$pvalue,digit=3))
pdf.file <- paste0('activity_compare_',per.type,'_Pmax',round(Dt),'_Pmin',round(1/fmax),'_ofac',ofac,'_Nma',Nma,'.pdf')
cat('output pdf\n',pdf.file,'\n')
pos.leg <- 'topright'
size <- 2
pdf(pdf.file,16,8)
par(mfrow=c(2,4),mar=c(0,0,0,0),oma=c(7,7,1,1),cex.axis=size,cex.lab=size)
if(FALSE){
power <- rv.per$power
plot(rv.per$P,power,xlab='Period[day]',ylab='Power of RV',log='x',type='l',xaxt='n',yaxt='n',ylim=c(min(power),dy*max(power)))
abline(h=rv.per$sig.level[1],lty=2)
abline(h=rv.per$sig.level[2],lty=2)
abline(h=rv.per$sig.level[3],lty=2)
abline(v=Popt,col='red',lty=3,lwd=lwd)
abline(v=Prot,col='steelblue')
magaxis(side=2)    
#abline(v=c(13.9,35.4,94,168,640),col='blue',lty=3,lwd=lwd)
legend(pos.leg,legend='RV',bty='n',cex=size)
}
####indexes


####plots for proxies
###Sindex
for(j in 1:(ncol(proxies)+1)){
    Nma <- 0
    Inds <- 0
    if(j<=ncol(proxies)){
    if(per.type=='BFP'){
        tmp <- BFP(t=t,y=scale(proxies[,j]),dy=eproxies[,j],ofac=ofac,fmax=fmax,fmin=1/Dt,progress=FALSE,Inds=Inds,Nma=Nma)
    }else{
        if(all(eproxies[,j]==0)){
            err <- rep(0.1,length(proxies[,j]))
        }else{
            err <- eproxies[,j]/sd(proxies[,j])
        }
        tmp <- glst(t=t,y=scale(proxies[,j]),err=err,ofac=ofac,fmax=fmax,fmin=1/Dt)
    }
    plot(tmp$P,tmp$power,xlab='Period[day]',ylab='Power',log='x',type='l',xaxt='n',yaxt='n',ylim=c(min(tmp$power),dy*max(tmp$power,tmp$sig.level[1])))
    if(j==1){
        mtext('Power',side=2,line=5,outer=TRUE,cex=1.5)
    }else if(j==5){
        mtext('Period [day]',side=1,line=5,outer=TRUE,cex=1.5)
    }
}else{
    tmp <- lsp(times=t,x=rep(1,length(t)),ofac=ofac,to=fmax,from=1/Dt)
    plot(tmp$P,tmp$power,xlab='Period[day]',ylab='Power',log='x',type='l',main=tit,xaxt='n',yaxt='n',ylim=c(min(tmp$power),max(tmp$power)+0.1*(max(tmp$power)-min(tmp$power))))
}
    if(FAP){
        abline(h=tmp$sig.level[1],lty=2)
        abline(h=tmp$sig.level[2],lty=2)
        abline(h=tmp$sig.level[3],lty=2)
    }
#    if(j==1|j==5){
    if(FALSE){
        magaxis(side=2,usepar=TRUE)    
    }
    if(j>4){
        magaxis(side=1,tcl=-1,mgp=c(3,2,0))  
    }
    abline(v=Popt,col='red',lty=3,lwd=lwd)
    abline(v=Prot,col='steelblue')
    legend(pos.leg,legend=names[j],bty='n',cex=size)
}


###window function
if(FALSE){
    plot(win.per$P,win.per$power,xlab='Period[day]',ylab='Power of window function',log='x',type='l',main=tit,xaxt='n',yaxt='n',ylim=c(min(win.per$power),max(win.per$power)+0.1*(max(win.per$power)-min(win.per$power))))
    if(FAP){
        abline(h=win.per$sig.level[1],lty=2)
        abline(h=win.per$sig.level[2],lty=2)
        abline(h=win.per$sig.level[3],lty=2)
    }
    abline(v=Popt,col='red',lty=3,lwd=lwd)
abline(v=Prot,col='steelblue')
    legend(pos.leg,legend='Window Function',bty='n',cex=size)
}
dev.off()

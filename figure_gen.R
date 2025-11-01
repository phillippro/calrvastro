####make plots; note that the plot_func.R can only run once since the mcmc.out would be thinning again if the code is run twice
##load data
##optional
#load(paste0(fname,'.Robj'))
#fname <- substr(fname,26,nchar(fname))
###
##########posterior over parameters
#global parameters
source('mcmc_func.R',local=TRUE)
source('periodograms.R')
source('data2figure.R')
if(Nw<=7){
    cols <- c('black','blue','green','purple','cyan','brown','pink')#colors used to show different data sets
}else{
    cols <- rainbow(Nw)#colors used to show different data sets
}
####naming the output file
plot.figure <- TRUE
if(!exists('emax')){
    emax <- 1
}

if(!any(fix.par=='')){
    fix.name <- paste0('fix',fix.par,collapse='and')
}else{
    fix.name <- 'nofix'
}
if(!exists('data.type')){
   data.type <- 'normal'
}
if(exists('ep.all')){
   Nep.all <- length(ep.all)
}else{
   Nep.all <- 0
}
fname <-  paste0('kep',kep.type,'_prior',prior.type,'_poly',Npoly,Npoly.sub,'_Ndata',Ndata,'_quantify',quantify,'_',Nbin.per,'per',nbin.per,'_Nw',Nw,'_',ids[1],'_',Np,'planet_',noise.models[1],'_Nsamp',Niter,'_tem',format(tem,digit=3),'_Npar',length(par.opt),'_acc',format(acceptance*100,digit=2))
if(!exists('tem0')){
    tem0 <- ''
}
if(length(Popt)>4){
    fname<- paste0(fname,'_pretem',format(tem0,digit=3),'P',paste0(round(c(Popt[1:3],Popt[length(Popt)]),digit=1),collapse='d'),'d_negLmax',format(abs(max(loglike.out)),digits=3))
}else{
    fname<- paste0(fname,'_pretem',format(tem0,digit=3),'P',paste0(round(Popt,digit=1),collapse='d'),'d_negLmax',format(-1*max(loglike.out),digits=3))
}
folder.name <- '/car-data/ffeng/dwarfs/output/'
if(!file.exists(folder.name)){
    folder.name <- '../output/'
}
folder.name <- paste0(folder.name,star.name,'/')
cat('target=',target,'\n')
cat('ID=',ID,'\n')
#folder.name <- paste0(folder.name,ID,'/')
if(!file.exists(folder.name)){
    system(paste('mkdir',folder.name))
}
fname <- paste0(folder.name,fname)
#####################################
###make graphs
#####################################
if(plot.figure){
pdf.name <- paste(fname,'.pdf',sep='')
pdf(pdf.name,width=16,height=16)
par(mfrow=c(4,4),mar=c(5,5,5,1),cex.axis=1.5,cex.lab=1.5)
####statistics of parameter inference
par.stat <- c()
####
################################################################
###partI: plot of the data, periodograms and residuals
################################################################    
#################################################
####phase-folded data and model prediction
#################################################
RV2.all <- c()
for(k1 in 1:Nw){
    var <- names(par.data[[ids[k1]]])
    for(k in 1:length(var)){
        assign(var[k],par.data[[ids[k1]]][[var[k]]])
    }
    source('data_periodogram.R')
    RV2.all <- c(RV2.all,RV2)
}
########################################################
###PartII: plots for combined data sets
########################################################
if(Nw>1){
####RV data
trv.all <- c()
RV.all <- c()
eRV.all <- c()
res.all <- c()
for(i3 in 1:Nw){
    var <- names(par.data[[ids[i3]]])
    for(k in 1:length(var)){
        assign(var[k],par.data[[ids[i3]]][[var[k]]])
    }
    trv.all <- c(trv.all,trv)
    RV.all <- c(RV.all,RV)
    eRV.all <- c(eRV.all,eRV)
    res.all <- c(res.all,res)
    if(i3==1){
        plot(trv,RV2,xlab=expression((t-t[0])*'[days]'),ylab='RV [m/s]',ylim=range(RV2.all),main='Combined data set',col=cols[i3],pch=20,cex=0.4,xlim=c(tmin,tmax))
    }else{
        points(trv,RV2,pch=20,cex=0.4,col=cols[i3])
    }
    arrows(trv,RV2-eRV,trv,RV2+eRV,length=0.05,angle=90,code=3,col=cols[i3])
}
lines(ts,RV.kepler(pars.kep=par.opt,tt=ts,kep.only=TRUE)[[1]],col="red")
legend('topright',legend=paste0('RMS=',format(sd(RV[!is.na(RV)]),digit=3),'m/s'),bty='n')
####phase folded
data.combine <- TRUE
if(Np>0){
    source('phase_fold.R')
}
####combined data periodogram
ind <- sort(trv.all,index.return=TRUE)$ix
trv.all <- trv.all[ind]
RV.all <- RV.all[ind]
eRV.all <- eRV.all[ind]
rv.glst <- glst(t=trv.all,y=RV.all,err=eRV.all,ofac=ofac)
plot(rv.glst$P,rv.glst$logp,xlab='Period[day]',ylab='Power of Combined RV',log='x',type='l',main='')
#paste0('GLST; p value:',format(rv.glst$pvalue,digit=3))
abline(h=rv.glst$sig.level[1],lty=2)
abline(h=rv.glst$sig.level[2],lty=3)
abline(h=rv.glst$sig.level[3],lty=4)
abline(v=Popt,lty=2)
####combined window function
win.bgls <- bgls(t=trv.all,y=rep(mean(RV.all),length(trv.all)),err=eRV.all, ofac=ofac,fmax=fmax)
plot(win.bgls$P,win.bgls$logp,xlab='Period[day]',ylab='Power of window function',log='x',type='l',main=tit)
abline(v=Popt,lty=2)
##############################################################################
####residual
for(i3 in 1:Nw){
    var <- names(par.data[[ids[i3]]])
    for(k in 1:length(var)){
        assign(var[k],par.data[[ids[i3]]][[var[k]]])
    }
    if(i3==1){
#ylim=c(rmin,rmax),
        plot(trv,res,xlab=expression((t-t[0])*'[days]'),ylab='RV(data-model) [m/s]',main='Residual of Combined Data Set',col=cols[i3],pch=20,cex=0.4,xlim=c(tmin,tmax),ylim=range(RV2.all))
    }else{
        points(trv,res,pch=20,cex=0.4,col=cols[i3])
    }
    arrows(trv,res-eRV,trv,res+eRV,length=0.05,angle=90,code=3,col=cols[i3])
}
####residual periodogram
res.all.glst <- glst(t=trv.all,y=res.all,err=eRV.all,ofac=ofac)
plot(res.all.glst$P,res.all.glst$power,xlab='Period[day]',ylab='Power of Combined Residuals',log='x',type='l',main='')
ind <- which.max(res.all.glst$power)
if(length(ind)>0){
    pmax <- res.all.glst$P[ind]
    power.max <- res.all.glst$power[ind]
    par(xpd=TRUE)
    text(pmax,1.1*power.max,pos=3,labels=format(pmax,digit=3),col='red',cex=1.0)
    arrows(pmax,1.1*power.max,pmax,1.05*power.max,col='red',length=0.05)
    par(xpd=FALSE)
}
abline(h=res.all.glst$sig.level[1],lty=2)
abline(h=res.all.glst$sig.level[2],lty=3)
abline(h=res.all.glst$sig.level[3],lty=4)
}
###############################################################################
###save data; put it here because I want to save the res.all.glst periodograms
################################################################################
obj.name <- paste(fname,'.Robj',sep='')
cat('Robj name: ',obj.name,'\n')
#save(list=c('mcmc.out',list(ls(all=TRUE))[[1]]),file=obj.name)
remove.vars <- c()
if(exists('Np1')){
    remove.vars <- c(remove.vars,'Np1')
}
if(exists('fs')){
    remove.vars <- c(remove.vars,'fs')
}
if(exists('newpar0')){
    remove.vars <- c(remove.vars,'newpar0')
}
if(exists('estimation')){
    remove.vars <- c(remove.vars,'estimation')
}
if(exists('Nbin.per1')){
    remove.vars <- c(remove.vars,'Nbin.per1','nbin.per1','res.glst.all')
}
save(list=setdiff(ls(), remove.vars),file=obj.name)
#######################################
###part III: trace plot; likelihood and posterior distribution
#######################################
#posteriors and likelihoods
plot(index.mcmc,post.out1,type='l',xlab='Iterations',ylab='Logarithm Posterior')
plot(index.mcmc,loglike.out1,type='l',xlab='Iterations',ylab='Logarithm Likelihood')
label.expressions <- plot.labels(Np=Np,noise.model=noise.model,p=p,q=q)
#tracing parameters
ind.tmp <- 0
opt.par.arr <- array(data=NA,dim=c(10,Npar+2,Np))#at most 10 peaks per planet model
for(j in 1:Npar){
    if(Np>0 & any(j==ind.per)){
        tit <- paste(which(j==ind.per),'planet',sep='')
        #trace
        nam <- names(startvalue)[j]
	if(!fixP & any(ind.per==j) & !any(fix.par==nam)){
        plot(index.mcmc,P.out[,ind.tmp+1],type='l',xlab='Iterations',ylab=label.expressions[j+2*ind.tmp],log='y',main=tit)
        plot(index.mcmc,P.out[,ind.tmp+1],type='l',xlab='Iterations',ylab=label.expressions[j+1+2*ind.tmp],main=tit)
        plot(index.mcmc,nv.out[,ind.tmp+1],type='l',xlab='Iterations',ylab=label.expressions[j+2+2*ind.tmp],main=tit)
        ylim.like <- c(mean(like.bin[,j]),max(like.bin[,j]))
        ylim.post <- c(mean(post.bin[,j]),max(post.bin[,j]))
#        cat('head(like.bin[,j])=',head(like.bin[,j]),'\n')
        #likelihood
        plot(P.bin[,ind.tmp+1],like.bin[,j],type='l',ylab='Logarithm likelihood',xlab=label.expressions[j+2*ind.tmp],log='x',ylim=ylim.like,main=tit)
        plot(P.bin[,ind.tmp+1],like.bin[,j],type='l',ylab='Logarithm likelihood',xlab=label.expressions[j+2*ind.tmp+1],ylim=ylim.like,main=tit)
        plot(nv.bin[,ind.tmp+1],like.bin[,j],type='l',ylab='Logarithm likelihood',xlab=label.expressions[j+2*ind.tmp+2],ylim=ylim.like,main=tit)
        #posterior
        plot(P.bin[,ind.tmp+1],post.bin[,j],type='l',xlab=label.expressions[j+2*ind.tmp],ylab='Logarithm posterior',log='x',ylim=ylim.post,main=tit)
#some thresholds
###################
###highlight the peaks
        Pmax.global <- P.bin[which.max(post.bin[,j]),ind.tmp+1]
        ind <- which(post.bin[,j]>(max(post.bin[,j])+log(0.1)))
        P.peaks <- P.bin[ind,ind.tmp+1]
        per.peaks <- per.bin[ind,ind.tmp+1]
        post.peaks <- post.bin[ind,j]
        like.peaks <- like.bin[ind,j]     
        if((max(P.bin)-min(P.bin))<(0.3*Pmax.global) & mean(diff(P.peaks))<5 & length(P.peaks)>10){
            Nc2p <- Nc2p+1
        }
        for(m in 1:length(post.peaks)){
            if(m<=10){
                index <- which.min(abs(post.out-post.peaks[m]))[1]
                opt.par.arr[m,,ind.tmp+1] <- c(mcmc.out[index,],post.peaks[m],like.peaks[m])
            }
        }
        da <- 0.1
        abline(v=Pmax.global,col='red',lty=3,lwd=3)
        legend('topright',legend=format(Pmax.global,digit=3),bty='n',col='red')
#        par(xpd=TRUE)
#        text(P.peaks,rep(da*diff(ylim)+max(post.bin[,j]),length(ind)),pos=3,labels=format(P.peaks,digit=3),col='red',cex=1.0)
#        par(xpd=FALSE)
#        arrows(P.peaks,rep(da*diff(ylim)+max(post.bin[,j]),length(ind)),P.peaks,rep(0.5*da*diff(ylim)+max(post.bin[,j]),length(ind)),col='red',length=da)
###############
        plot(P.bin[,ind.tmp+1],post.bin[,j],type='l',ylab='Logarithm Posterior',xlab=label.expressions[j+2*ind.tmp+1],ylim=ylim.post,main=tit)
        plot(nv.bin[,ind.tmp+1],post.bin[,j],type='l',ylab='Logarithm Posterior',xlab=label.expressions[j+2*ind.tmp+2],ylim=ylim.post,main=tit)
###histograms
        data.distr(logP[,ind.tmp+1],'log(Period[day])','Frequency',main=tit)
        p1 <- data.distr(P[,ind.tmp+1],label.expressions[j+2*ind.tmp+1],'Frequency',main=tit)
        par.stat <- cbind(par.stat,p1)
        data.distr(nv[,ind.tmp+1],label.expressions[j+2*ind.tmp+2],'Frequency',main=tit)
	}else{
            p1 <- c(exp(mcmc.out[1,j]),rep(0,10))
            par.stat <- cbind(par.stat,p1)
        }
        ind.tmp <- ind.tmp+1
    }else{
        nam <- names(startvalue)[j]
        if(!any(fix.par==nam)){
            if(j>(Np*Nkeppar)){
                                        #            nap <- as.integer(unlist(strsplit(nam,''))[nchar(nam)])
                tit.data <- ''
            }else{
                tit.data <- ''
            }
            lab <- label.expressions[j+2*ind.tmp]
####Note: beta=log(tau)
            if(grepl('beta',names(par.opt)[j])){
                str <- unlist(strsplit(names(par.opt)[j],''))
                i2 <- as.integer(str[length(str)])
                plot(index.mcmc,exp(mcmc.out1[,j]),type='l',xlab='Iterations',ylab=expression(tau*'[day]'),main=tit.data,log='y')
                                        #posterior
                betas <- grep(paste0('beta',i2),names(par.opt))
                ind <- which(betas==j)
                plot(tau.bin[,ind,i2],post.bin[,j],type='l',ylab='Logarithm Posterior',xlab=expression(tau*'[day]'),main=tit.data,log='x')
                                        #histograms
                p2 <- data.distr(mcmc.out[,j],expression('log('*tau*'[day])'),'Frequency',startvalue[j],main=tit.data)
            }else{
                                        #trace
                plot(index.mcmc,mcmc.out1[,j],type='l',xlab='Iterations',ylab=lab,main=tit.data)
                                        #posterior
                plot(par.bin[,j],post.bin[,j],type='l',ylab='Logarithm Posterior',xlab=lab,main=tit.data)
                                        #histograms
                p2 <- data.distr(mcmc.out[,j],lab,'Frequency',startvalue[j],main=tit.data)
            }
            par.stat <- cbind(par.stat,p2)
        }else{
        p1 <- c(mcmc.out[1,j],rep(0,10))
        par.stat <- cbind(par.stat,p1)
    }
    }
}
dm <- par.opt.post-par.stat[1,]
dp <- par.stat[2,]-par.opt.post
par.stat <- rbind(par.opt.post,dm,dp,par.stat)
colnames(par.stat) <- names(par.opt)
rownames(par.stat) <- c('opt','dminus','dplus','x1per','x99per','x10per','x90per','xminus.1sig','xplus.1sig','mode','mean','sd','skewness','kurtosis')
dev.off()
#################################
###save data
##################################
####figures for papers and presentations
#source('paper_plot.R')
###an initial criteria to judge whether the chain converge to the global maxima
cat('conv=',conv,'\n')
if(conv==TRUE){
    fname1<- paste0(fname,'-c2p',Nc2p)
}else{
    fname1 <- fname
}
tmp <- paste0(fname1,'.pdf')
file.rename(pdf.name,tmp)
pdf.name <- tmp
cat('output pdf:\n')
cat(pdf.name,'\n')
par.file <- paste0(fname,'_optpar.Robj')
save(conv,Np,Nc2p,opt.par.arr,par.stat,par.opt,par.opt.like,par.opt.post,file=par.file)
cat('output par file:',par.file,'\n')
#obj2.name <- paste(fname,'_plot.Robj',sep='')
#cat('Robj name: ',obj2.name,'\n')
#save(list=ls(all=TRUE),file=obj2.name)
}

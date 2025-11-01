##This code is to compare the imaged positions and predicted positions of a companion based on RV+Hipparcos+Gaia solutions.library(RColorBrewer)
library(RColorBrewer)
source('mcmc_func.R')
nepoch <- 2
astro.type <- 'hg3'
fall <- 'results/HD4747/HD4747_fix1_relativityFALSE_Niter2100000_Ncores16_ofac2_Nset3_hg3_transit0_P11987_acc10_sinI_lnlmax-300.Robj'
fhg2 <- 'results/HD4747/HD4747_fix1_relativityFALSE_Niter2200000_Ncores16_ofac2_Nset3_astro_transit0_P11999_acc2.4_sinI_lnlmax-300.Robj'
fhgca2 <- 'results/HD4747/HD4747_fix1_relativityFALSE_Niter2000000_Ncores8_ofac2_Nset3_hg2--_transit0_P12091_acc6.9_sinI_lnlmax-301.Robj'
fhgca3 <- 'results/HD4747/HD4747_fix1_relativityFALSE_Niter2100000_Ncores16_ofac2_Nset3_hgca3--_transit0_P12067_acc16_sinI_lnlmax-296.Robj'
fhg3 <- 'results/HD4747/HD4747_fix1_relativityFALSE_Niter1200000_Ncores16_ofac2_Nset3_hg3_transit0_P12118_acc4.9_sinI_lnlmax-292.Robj'
solutions <- c('RV+HEDR3+IMG','RV+HDR2','RV+HGCA2','RV+HEDR3','RV+HGCA3')
if(!exists('out')){
    load(fall)
    par.opt0 <- par.opt
}
if(!exists('full')) load(fall,envir=full <- new.env())
if(!exists('hg2')) load(fhg2,envir=hg2 <- new.env())
if(!exists('hg3')) load(fhg3,envir=hg3 <- new.env())
if(!exists('hgca2')) load(fhgca2,envir=hgca2 <- new.env())
if(!exists('hgca3')) load(fhgca3,envir=hgca3 <- new.env())
data.astrometry <- out$astrometry
if(out$Nrv>0){
    Ein  <-  RV.kepler(pars.kep=par.opt,bases=bases)$E
}else{
    Ein <- NULL
}

####simulation
nc <- 1
source('../../pexo/code/timing_function.R')
source('../../pexo/code/general_function.R')
source('../../pexo/code/sofa_function.R')
                                        #    ind <- sort(Popt,index.return=TRUE)$ix
trel <- out$trel
                                        #    Dt <- Popt[j]-(max(trel)-min(trel))
Dt <- max(0,(max(trel)-min(trel))-Popt)
#tsim <- seq(min(trel)-Dt/2,max(trel)+Dt/2,length.out=1e3)
tsim <- seq(min(trel)-Dt,max(trel)+Dt,length.out=1e4)


yrs <- time_Jd2yr(cbind(tsim,0))
                                        #    tts <- seq(max(floor(min(yrs)),1990),min(ceiling(max(yrs)),2030),length.out=10)
tts <- seq(max(floor(min(yrs)),1990),min(ceiling(max(yrs)),2030),length.out=5)
tts <- unique(round(tts))
                                        #    tts <- tts[tts>1990 & tts<2030]
                                        #    tts <- tts[tts>2010 & tts<2030]
                                        #    tts <- c(tts,2022.1)
ii <- unlist(lapply(tts, function(tt) which.min(abs(yrs-tt))))
jjs <- trel[ii]

cc <- c('black','blue','orange','steelblue','green','darkgreen',rainbow(10))
cols.model <- brewer.pal(5, "Set1")
fpdf <- 'HD4747_solution_comparison.pdf'
cat(fpdf,'\n')
pdf(fpdf,6,6)
size <- 1.2
par(mar=c(5,5,4,1),cex.lab=size)
Nsamp <- 0
for(i1 in 1:5){
    for(j1 in 0:Nsamp){
        if(j1==0){
            if(i1==1){
                par.opt <- par.opt0
            }
            if(i1==2){
                par.opt <- hg2$par.opt
            }
            if(i1==3){
                par.opt <- hgca2$par.opt
            }
            if(i1==4){
                par.opt <- hg3$par.opt
            }
            if(i1==4){
                par.opt <- hgca3$par.opt
            }
            Npar <- length(par.opt)
        }else{
            if(i1==1){
                mc <- full$out$mcmc.opt$sig1
            }
            if(i1==2){
                mc <- hg2$out$mcmc.opt$sig1
            }
            if(i1==3){
                mc <- hgca2$out$mcmc.opt$sig1
            }
            if(i1==4){
                mc <- hg3$out$mcmc.opt$sig1
            }
            if(i1==5){
                mc <- hgca3$out$mcmc.opt$sig1
            }
            k1 <- sample(1:nrow(mc),1)
            par.opt <- mc[k1,1:Npar]
        }
        col.model <- cols.model[i1]
        planet <- astrometry.rel(par.opt,tsim=tsim)
        range.raS <- range(unlist(planet$rel$raB))
        range.decS <- range(unlist(planet$rel$decB))
        j <- 1
        n <- paste0('p',j)
        inss <- names(out$rel[[n]])
        inss <- inss[inss!='tot']
        draS <- planet$rel$raB[[n]]
        ddecS <-planet$rel$decB[[n]]
        if(j1==0){
            range.raP <- range(unlist(lapply(names(out$rel),function(n) unlist(lapply(names(out$rel[[n]]),function(i) out$rel[[n]][[i]][,'dra'])))))
            range.decP <- range(unlist(lapply(names(out$rel),function(n) unlist(lapply(names(out$rel[[n]]),function(i) out$rel[[n]][[i]][,'ddec'])))))

#########relative astrometry part

            tPs <- draPs <- ddecPs <- c()
            for(i in inss){
                draPs <- c(draPs,out$rel[[n]][[i]][,'dra'])
                ddecPs <- c(ddecPs,out$rel[[n]][[i]][,'ddec'])
                tPs <- c(tPs,out$rel[[n]][[i]][,1])
            }
        }


####pdf plots
                                        #    xlim <- range(range.raP,range.raS)
                                        #    ylim <- range(range.decP,range.decS)
        xlim <- c(-150,250)
        ylim <- c(-700,-490)
        if(Nsamp>0){
            alpha <- 0.2
            beta <- 1
            xlim <- xlim+c(-alpha*abs(diff(xlim)),alpha*abs(diff(xlim)))
            ylim <- ylim+c(-beta*abs(diff(ylim)),beta*abs(diff(ylim)))
        }
        if(i1==1 & j1==0){
            plot(draPs,ddecPs,xlab=expression(Delta*alpha*'* [mas]'),ylab=expression(Delta*delta*' [mas]'),xlim=rev(xlim),ylim=ylim,col='white')
        }
        if(j1==0){
            for(i in inss){
                tP <- out$rel[[n]][[i]][,1]
                draP <- out$rel[[n]][[i]][,'dra']
                ddecP <- out$rel[[n]][[i]][,'ddec']
                edraP <- out$rel[[n]][[i]][,'edra']
                eddecP <- out$rel[[n]][[i]][,'eddec']
                points(draP,ddecP,col='black',pch=20,cex=0.8)
                points(0,0,pch='+',col='black',cex=2)
                if(target=='HD42581'){
                    t0 <- sum(time_Yr2jd(2022))
                    t1 <- sum(time_Yr2jd(2024))
                    jj <- which(tsim>t0&tsim<t1)
                    lines(draS[jj],ddecS[jj],col='purple',lwd=4,lty=1)
                    rho <- sqrt(draS[jj]^2+ddecS[jj]^2)/1000#as
                    cat('rho=',rho,'as\n')
                }
                pos <- rep(3,length(ii))
                pos[which(draS[ii]>mean(draS[ii]))] <- 4
####show years
                if(TRUE){
#                if(FALSE){
                    points(draS[ii],ddecS[ii],pch='x',col=col.model,cex=0.5)
                    if(i1==1 | i1==4) pos <- rep(1,length(ii))
                    text(draS[ii],ddecS[ii],labels=tts,pos=pos,xpd=NA,col=col.model,cex=0.6)
                }
                inds <- c()
                for(k in 1:length(draP)){
###                inds <- c(inds,which.min(abs(out$rel[[n]]$tot[k,1]-tsim)))
                    inds <- c(inds,which.min(abs(tP[k]-tsim)))
                }
                points(draS[inds],ddecS[inds],pch=20,col=col.model)
                arrows(draS[inds],ddecS[inds],draP,ddecP,length=0.001,angle=0,code=1,col=tcol(col.model,50))
                arrows(draP-edraP,ddecP,draP+edraP,ddecP,length=0.01,angle=90,code=3)
                arrows(draP,ddecP-eddecP,draP,ddecP+eddecP,length=0.01,angle=90,code=3)
                nc <- nc+1
            }
            lab.ins <- inss
        }
###add model prediction
        if(j1==0){
            lines(draS,ddecS,col=col.model,lty=1)
        }else{
            lines(draS,ddecS,col=tcol(col.model,80),lty=1)
        }
    }
}
legend('topright',inset=c(-0.1,0),legend=solutions,col=cols.model,bty='n',pch=NA,lty=1,xpd=NA,cex=0.8)
dev.off()

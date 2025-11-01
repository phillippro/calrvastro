source('mcmc_func.R')
library(magicaxis)
library(RColorBrewer)
#if(!exists('out')){
#    load('results/HD209100/HD209100_calFALSE_coplanarFALSE_resonanceFALSE_relativityFALSE_Niter2200000_Ncores16_ofac2_Nset5_hg123_transit0_P15663_Esd1_acc1.5_rvcbarycentric_lnlmax-8150.Robj')
#    load('results/HD22049/HD22049_calFALSE_coplanarFALSE_resonanceFALSE_relativityFALSE_Niter2100000_Ncores16_ofac2_Nset13_hg123_transit0_P2695_Esd1_acc3_rvcbarycentric_lnlmax-4267.Robj')
#    load('results/HD222237/HD222237_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter120000_sppriorFALSE_dr1FALSE_230617_Nset4_hg123_Nsig1_P22120_Esd1_astro5TRUE_acc0.93_rvcbarycentric_lnpmax-707_lnlmax-679.Robj')
#}
#load('results/HD222237/HD222237_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter1100000_sppriorFALSE_dr1FALSE_230617_Nset4_hg123_Nsig1_P21646_Esd1_astro5TRUE_acc0.81_rvcbarycentric_lnpmax-705_lnlmax-677.Robj')
load('results/UCAC4569-026385/UCAC4569-026385_calFALSE_coplanarFALSE_resonanceFALSE_staFALSE_relativityFALSE_Niter2200000_dr1FALSE_230617_Nset1_hg123_Nsig1_P879_Esd1_astro5TRUE_acc0.29_rvcbarycentric_priorf3_lnpmax-225_lnlmax-200.Robj')
Np <- Nsig
mc0 <- mc <- out$mcmc.opt[[paste0('sig',Nsig)]]

getLevel <- function(kk,prob=0.95){
    dx <- diff(kk$x[1:2])
    dy <- diff(kk$y[1:2])
    sz <- sort(kk$z)
    c1 <- cumsum(sz) * dx * dy
    approx(c1, sz, xout = 1 - prob)$y
}
#sigmas <- c(0.682,0.954, 0.997)
#sigmas <- c(0.682,0.954)
#sigmas <- 0.682
label.expressions <- plot.labels.simple(names(par.opt))
Npar <- 7
labs <- nams <- label.expressions[1:Npar]
if(target=="UCAC4569-026385") labs <- gsub(' \\[b\\]','',labs)

indP <- grep('per',names(par.opt))
popt <- exp(par.opt[indP])#days
ps <- exp(mc0[,indP])#day
###convert to years
if(all(popt>1e3)){
    popt <- popt/365.25
    labs[indP] <- gsub('day','yr',labs[indP])
    ps <- ps/365.25
}
Nsamp <- nrow(mc)
#Nsamp <- min(nrow(mc),5e5)
#Nsamp <- 1e4

#par.opt[indP] <- popt
mc[,indP] <- ps
dc <- matrix(mc[,ncol(mc)],ncol=1)
colnames(dc) <- colnames(mc)[ncol(mc)]
nlevel <- 3
my.cols <- rev(brewer.pal(nlevel, "RdYlBu"))
#sigmas <- 1-(1-pnorm(1:5))*2#1 to 5 sigmas
sigmas <- 1-(1-pnorm(1))*2#1 to 5 sigmas
fpair <- paste0(fname,'_sample.pdf')
fpair <- paste0('~/Desktop/',target,'_post_pair_N',round(Nsamp/1e3),'k.pdf')
#contourf <- FALSE
contourf <- TRUE
pdf(fpair,width=16,height=16)
#Npar <- length(par.opt)
                                        #Nbasic <- 5*Np+10
size <- 1
par.names <- colnames(mc)
par(mfcol=c(Npar,Npar),mar=c(0,0,0,0),oma=c(4,4,2,1),cex.axis=0.8,cex.lab=size,cex=size)
nsamp <- 1e4
inds <- sample(nrow(mc),Nsamp)
par.name <- names(par.opt)
ii <- sample(inds,nsamp)
for(j in 1:Npar){
    for(k in 1:Npar){
        x <- mc[inds,j]
        y <- mc[inds,k]
#        if(grepl('lnP',labs[j])){
#            labs[j] <- "P[b]*'[day]'"
#            x <- exp(x)
#        }
#        if(grepl('lnP',labs[k])){
#            labs[k] <- "P[b]*'[day]'"
#            y <- exp(y)
#        }
        xlab <- labs[j]
        ylab <- labs[k]
        if(target=="UCAC4569-026385" & grepl('^K',xlab)){
            x <- x/1e3
            xlab <- gsub('\\[m','\\[km',xlab)
        }
        if(target=="UCAC4569-026385" & grepl('^K',ylab)){
            y <- y/1e3
            ylab <- gsub('\\[m','\\[km',ylab)
        }
        if(k>j){
            if(contourf)   kk <- kde2d(x,y, n=10)
            my.cols <- 'black'
            if(contourf) levels <- getLevel(kk,prob=sigmas)
                                        #            levels <- pretty(zlim, n=1)#68,95, and 99% contours
#            if(length(inds)<=1e3 & FALSE){
            if(TRUE){
                xaxt='n'
                yaxt='n'
                if(j==1) yaxt <- 's'
                if(k==Npar) xaxt <- 's'
                plot(x[ii],y[ii],xaxt='n',yaxt='n',cex=0.5,pch=20,col='grey')
#                plot(x[jj],y[jj],xaxt='n',yaxt='n',cex=0.5,pch=20,col='grey')
                side <- c()
                if(xaxt=='s'){
                    side <- c(side,1)
                }
                if(yaxt=='s'){
                    side <- c(side,2)
                }
                if(length(side)>0){
                    magaxis(side)
                    magaxis((1:4)[-side],label=FALSE)
                }else{
                    magaxis(1:4,label=FALSE)
                }
               if(contourf) contour(kk, drawlabels=FALSE, nlevels=length(levels), levels=levels,col=my.cols,add=TRUE,axes=FALSE)
            }else{
#                contourLines(kk, drawlabels=FALSE, nlevels=length(levels), levels=levels,col=my.cols,xaxt='n',yaxt='n')
                if(contourf) contourLines(kk$x,kk$y,kk$z, nlevels=length(levels),levels=levels)
            }
                                        #            legend('topleft',inset=c(-0.2,-0.1),legend=paste0(format(cor(x,y),digit=2)),bty='n',text.col='red',cex=size)
        }else if(k==j){
            hist(x,main='',axes=FALSE,breaks=20,col = "gray",border=FALSE)
            magaxis(c(1,2,4),labels=FALSE)
#            magaxis(3,major=3)
            magaxis(3)
            tmp <- data.distr(x,plotf=FALSE)
            abline(v=c(tmp[c('xlow','xup')]),lwd=2,lty=1,col='blue')
                                        #            abline(v=c(tmp['median']),lwd=2,lty=2)
            abline(v=tmp['med'],lwd=2,lty=1,col='red')
                                        #            ind <- which.max(post.out)
                                        #            abline(v=mc[ind,j],lwd=2,lty=2)
        }else{
            plot(1, type="n", axes=F, xlab="", ylab="",main='')
        }

        if(j==1 & k>=j){
            if(j>1 | k>1) mtext(side=2,text=parse(text=ylab),line=2,cex=1.3*size)
        }
        if(k==Npar& k>=j){
#            if(j<Npar | k!=j) mtext(side=1,text=parse(text=xlab),line=2,cex=1.3*size)
            mtext(side=1,text=parse(text=xlab),line=2,cex=1.3*size)
        }
#            stop()
#        cat('k=',k,';j=',j,'\n')
    }

}
dev.off()
cat('output MCMC posterior samples:\n',fpair,'\n')

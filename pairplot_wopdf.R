library(RColorBrewer)
Np <- Nsig
mc <- out$mcmc.opt[[paste0('sig',Nsig)]]

getLevel <- function(kk,prob=0.95){
    dx <- diff(kk$x[1:2])
    dy <- diff(kk$y[1:2])
    sz <- sort(kk$z)
    c1 <- cumsum(sz) * dx * dy
    approx(c1, sz, xout = 1 - prob)$y
}
#sigmas <- c(0.682,0.954, 0.997)
sigmas <- c(0.682,0.954)
label.expressions <- plot.labels.simple(names(par.opt))
Npar <- 7
labs <- nams <- label.expressions[1:Npar]

indP <- grep('per',names(par.opt))
popt <- exp(par.opt[indP])
ps <- exp(mc[indP,])#day
###convert to years
if(all(popt>1e3)){
    popt <- popt/365.25
    labs[indP] <- gsub('day','yr',labs[indP])
    ps <- ps/365.25
}
par.opt[indP] <- popt
mc[indP,] <- ps
dc <- matrix(mc[,ncol(mc)],ncol=1)
colnames(dc) <- colnames(mc)[ncol(mc)]
nlevel <- 3
my.cols <- rev(brewer.pal(nlevel, "RdYlBu"))
#sigmas <- 1-(1-pnorm(1:5))*2#1 to 5 sigmas
sigmas <- 1-(1-pnorm(1:2))*2#1 to 5 sigmas
size <- 1
par.names <- colnames(mc)
par(mfcol=c(Npar,Npar),mar=c(0,0,0,0),oma=c(4,4,2,1),cex.axis=0.8,cex.lab=size,cex=size)
#Nsamp <- nrow(mc)
Nsamp <- 1e5
#Nsamp <- 1e4
nsamp <- 1e3
inds <- sample(nrow(mc),Nsamp)
par.name <- names(par.opt)
ii <- sample(1:length(inds),nsamp)
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
        if(k>j){
            kk <- kde2d(x,y, n=10)
            my.cols <- 'black'
            levels <- getLevel(kk,prob=sigmas)
                                        #            levels <- pretty(zlim, n=1)#68,95, and 99% contours
#            if(length(inds)<=1e3 & FALSE){
            if(TRUE){
                xaxt='n'
                yaxt='n'
                if(j==1) yaxt <- 's'
                if(k==Npar) xaxt <- 's'
                plot(x[ii],y[ii],xaxt='n',yaxt='n',cex=0.5,pch=20,col='grey')
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
                contour(kk, drawlabels=FALSE, nlevels=length(levels), levels=levels,col=my.cols,add=TRUE,axes=FALSE)
            }else{
#                contourLines(kk, drawlabels=FALSE, nlevels=length(levels), levels=levels,col=my.cols,xaxt='n',yaxt='n')
                contourLines(kk$x,kk$y,kk$z, nlevels=length(levels),levels=levels)
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
    }
}
cat('output MCMC posterior samples:\n',fpair,'\n')

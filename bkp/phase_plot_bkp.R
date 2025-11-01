Ps <- pars <- Pshows <- par.shows <- c()
if(any(Ndet==length(noise.types)) | TRUE){
    ind <- which(Ndet==length(noise.types))
    ysim.cum <- 0
    for(j in 1:Nmax){
        ind.ma <- grep('MA',noise.types)
        ind.show <- 1
        if(length(ind.ma)>0){
            ind.show <- ind.ma
        }
        Pshow <- Psig[j,ind.show]
        par.show <- opt.par[[ind.show]][j,]
        par.shows <- rbind(par.shows,par.show)
        Pshows <- c(Pshows,Pshow)
        y <- ws[[ind.show]][,j]
        if(j>1){
            ycum <- ws[[ind.show]][,j]+rowSums(sigs[[ind.show]][,1:(j-1),drop=FALSE])
        }else{
            ycum <- y
        }
        tt2 <- tt%%Pshow
        tsim2 <- tsim%%Pshow
        ##        if(any(j==ind)){
        plot(tt2,y,xlab=paste0('orbital phase [days]'),ylab='RV (data-model) [m/s]')#pre-whitened data
        legend('topright',legend=paste0('P=',round(Pshow,2),'d'),bty='n',col='red')
        try(arrows(tt2,y+tab0[,3],tt2,y-tab0[,3],length=0.05,angle=90,code=3,col='grey'),TRUE)
        ysim <- par.show['A']*cos(2*pi/Pshow*tsim)+par.show['B']*sin(2*pi/Pshow*tsim)
        points(tsim2,ysim,col='red',pch=20,cex=0.2)
        ysim.cum <- ysim.cum+ysim
        ##        }
        if(j==Nmax & j!=1){
            plot(trv,ycum,xlab=paste0('JD'),ylab='RV (data-model) [m/s]',main='all signals')
            try(arrows(trv,ycum+tab0[,3],trv,ycum-tab0[,3],length=0.05,angle=90,code=3,col='grey'),TRUE)
            ysim <- par.show['A']*cos(2*pi/Pshow*tsim)+par.show['B']*sin(2*pi/Pshow*tsim)
            lines(tsim+min(trv),ysim.cum,col='red')
        }
    }
}



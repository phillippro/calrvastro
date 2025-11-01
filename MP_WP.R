###This routine is to calculate the time-varying and wavelength-varying periodgram
Nbin0 <- 10
dTs <- sort(diff(trv.all),decreasing=TRUE)
wps <- mps <- list()
Pkep <- out$Popt[[paste0('sig',Nsig)]]
#par.opt <- out$par.stat[[paste0('sig',Nsig)]]['xopt',]
t0 <- combine.list(lapply(out$ins.rv,function(i) out[[i]]$RV[,1]))
index <- sort(t0,index.return=TRUE)$ix

###calculate residual
tmp <- cal.residual(par.opt,bases=bases)
rvmn.list <- tmp$res.noise# rv-noise
rvma.list <- tmp$res.all#rv-signal-noise
y0 <- lapply(out$ins.rv,function(k) rvmn.list[[k]])#sort
y0 <- combine.list(y0)#combine list
rvmn.comb <- y0[index]#

for(j3 in 1:Nsig){
    y1 <- list()
    par1 <- rep(0,length(par.opt))
    if(any(grepl('Omega',names(par.opt)))) Npar <- 7
    par1[(j3-1)*Npar+1:Npar] <-par.opt[(j3-1)*Npar+1:Npar]
    names(par1) <- names(par.opt)
    if(any(names(par.opt)=='Mstar')) par1['Mstar'] <- par.opt['Mstar']
    rv.kep <- cal.residual(par1,bases=bases)$rv.sig# residual
#    rv.kep <- RV.kepler(pars.kep=par1,kep.only=TRUE)
    for(i in out$ins.rv){
        y1[[i]] <-rvma.list[[i]]+rv.kep[[i]]
    }
    y0 <- lapply(out$ins.rv,function(k) y1[[k]])#sort
    y0 <- combine.list(y0)#combine list
    y <- y0[index]#form a residual combined set
    if((tspan>10) & (Ndata>20) & (!any(dTs-median(dTs)>max(365,tspan/2))) & (Pkep[j3]<2*tspan)){
        if(sum(dTs[1:2])>tspan/2){
            Nbin <- 2
        }else{
            Nbin <- Nbin0
        }
        Nma <- 0#Nma.opt
        Nar <- 0

###check number of RVs in each bin
        fdt <- rev(seq(0.1,1,by=0.1))
        for(j in 1:length(fdt)){
            Dt <- fdt[j]*tspan
            cw <- check.window(trv.all,Dt,Nbin)
            if(all(cw>length(trv.all)/Nbin0)) break()
        }
        if(Nbin==2 | Dt==tspan) Dt <- tspan/2
###
        if(!exists('indices')){
            indices <- NULL
        }
        mp.res <- try(MP.norm(t=trv.all,y=y,dy=eRV.all,Dt=Dt,nbin=Nbin,ofac=ofac,fmin=fmin,fmax=fmax,per.type='BFP',Nma=Nma,Indices=indices),TRUE)
        if(class(mp.res)=='try-error'){
            mp.res <- try(MP.norm(t=trv.all,y=y,dy=eRV.all,Dt=Dt,nbin=Nbin,ofac=min(2,10*ofac),fmin=fmin,fmax=fmax,per.type='BFP',Nma=Nma,Indices=indices),TRUE)
        }
	if(class(mp.res)!='try-error'){
        tt <- trv.all-tmin
        xlim <- range(tt)
###plot
        xx <- mp.res$tmid-tmin
        yy <- mp.res$P
#        power <- (mp.res$powers-min(mp.res$powers))/length(xx)
        power <- (mp.res$powers-min(mp.res$powers))/(max(mp.res$powers)-min(mp.res$powers))
        zz <- t(power)
        mps[[paste0('sig',j3)]] <- list(xx=xx,yy=yy,zz=zz,y=y,power=mp.res$powers)
	}else{
        mps[[paste0('sig',j3)]] <- NULL
	}
    }
##calculate WP
    inds <- grep('9AP',names(out))
    if(length(inds)>0){
        for(ind in inds){
            dy <- y <- zz <- c()
            data <- out[[ind]]
            nn <- names(out)[ind]
            Nap <- (ncol(data)-1)/2
            xx <- 1:Nap

            ##subtract other signals
            if(Nsig>1){
                par.sig <- par.opt
                par.sig[(j3-1)*5+1:5] <- 0
                names(par.sig) <- names(par.opt)
                rv.sim <- RV.kepler(pars.kep=par.sig,tt=data[,1]%%2400000,kep.only=TRUE,bases=bases)
            }else{
                rv.sim <- rep(0,nrow(data))
            }

            for(j in 1:Nap){
                rv <- data[,2*(j-1)+2]-rv.sim
                erv <- data[,2*(j-1)+3]
                per <- BFP(data[,1],rv,erv,Nma=0,Nar=0,Indices=NULL,ofac=ofac,model.type='man',fmin=fmin,fmax=fmax,quantify=TRUE,progress=FALSE,GP=FALSE,noise.only=FALSE,Nsamp=1,gp.par=FALSE,renew=TRUE)
                y <- cbind(y,rv)
                dy <- cbind(dy,data[,2*(j-1)+3])
                p <- rev(per$power)
                                        #            power <- (p-min(p))/sum(1/erv^2)
                                        #            power <- (p-min(p))/nrow(data)
                power <- (p-min(p))/(max(p)-min(p))
                zz <- rbind(zz,as.numeric(power))#relative power
            }
            t <- data[,1]
            yy <- rev(per$P)
            wps[[paste0('sig',j3)]][[nn]] <- list(xx=xx,yy=yy,zz=zz,t=t,y=y,dy=dy)
        }
    }
}
out[['MP']] <- mps
out[['WP']] <- wps

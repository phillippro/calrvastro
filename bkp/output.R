###plot setting up
if(!file.exists(paste0('results/',target))){
    system(paste0('mkdir results/',target))
}
fname <- paste0('results/',target,'/',target,'_',ins,'_',paste(noise.types,collapse='_'),'_Esd',Esd,'_ofac',ofac,'_xi',xi,'_N',Niter0,'_new')
fpdf <- paste0(fname,'.pdf')
fobj <- paste0(fname,'.Robj')
pdf(fpdf,16,16)
par(mfrow=c(4,4))

if(nrow(tab)>2){
###plot data and BFP periodograms
    cat('periodograms and likelihoods for data and activity indicators!\n')
    source('data_periodogram.R')
    if(!noMCMC){
####model comparison
        source('model_selection.R')

###differential RVs diagnostics (later)
        Nsig <- 0
        mcmc.all <- res.all <- stat.all <- par.all <- list()#for different planet models
        logBF3 <- logBF5 <- logP.all <- logL.all <- c()
        bf3 <- bf5 <- 10
        Popt.all <- list()
####MCMC constrain of plausible signals using the optimal noise model
        source('setting_up.R',local=TRUE)
        RV.res <- RV
        for(j3 in 1:(Nmax+1)){
            cat('\n\nsignal=',j3-1,'\n')
            cat('\nrun mcmc!\n')
            source('run_mcmc.R')
            if(j3>1){
                bf3 <- logL.max -logL.all[length(logL.all)]-1.5*log(Ndata)
                bf5 <- logL.max -logL.all[length(logL.all)]-2.5*log(Ndata)
            }
            cat('bf3=',bf3,'\n')
            cat('bf5=',bf5,'\n')
            if(j3<=1 | bf3>=5){
                stat.all[[j3]] <- par.stat
                par.all[[j3]] <- par.opt
                if(j3>1){
                    if(period.par=='P'){
                        P <- par.opt[grep('per',names(par.opt))]
                    }else{
                        P <- exp(par.opt[grep('per',names(par.opt))])
                    }
                }else{
                    P <- NA
                }
                cat('mcmc.all j3=',j3,'\n')
                Popt.all[[j3]] <- P
                mcmc.all[[j3]] <- mcmc.out
                if(j3>1){
                    logBF3 <- c(logBF3,bf3)
                    logBF5 <- c(logBF5,bf5)
                }
                logP.all <- c(logP.all,logP.max)
                logL.all <- c(logL.all,logL.max)

###
                RV.kep <- RV.kepler(pars.kep=par.opt,kep.only=TRUE,Np.kep=j3-1)$rv[[1]]
                                        #    RV.kep <- RV.kepler(pars.kep=par.opt,Np.kep=j3-1)[[1]]
                RV.res <- RV-RV.kep

#######find the next signal through rednoise residual periodograms
                cat('\nFind ',j3,' signal!\n')
                per <- BFP(tt,RV.res,tab[,3],Nma=1,Nar=0,Indices=NULL,ofac=ofac,model.type='man',fmin=fmin0,fmax=fmax0,quantify=TRUE,progress=FALSE)
                Popt <- per$P[which.max(per$power)]
                Kopt <- min(as.numeric(sqrt(per$pars[which.max(per$power),'A']^2+per$pars[which.max(per$power),'B']^2)),Kmax)
                res.all[[j3]] <- RV.res
                if(j3>1){
                    Nsig <- Nsig+1
                }
            }else{
                cat('Failed to pass log(BF)>5\n')
                Npar <- Npar-Nkeppar
                break()
            }
        }
        for(j in 1:Npar){
            tmp <- data.distr(mcmc.all[[Nsig+1]][,j],xlab=colnames(mcmc.all[[Nsig+1]])[j],ylab='Freq.',plotf=TRUE)
        }

        if(Nsig>0){
### If there is any statistically significant signals, show the phase curve and MP; otherwise, don't show anything.
            index <- Nsig+1
            par.kep <- par.all[[index]]
            Pkep <- Popt.all[[index]]

####phase curve and report parameters and their uncertainties in pdf
            source('phase_plot.R')

                                        #    pdf('test.pdf',16,16)
            alpha <- 0.0
            size <- 2
            gammas <- c()
            par(cex.axis=size,cex.lab=size,cex=size)
###calculate moving periodogram
            cat('\ncalculate MP!\n')
            for(j3 in 1:Nsig){
                RV.res <- res.all[[j3]]
                if(tspan>10 & Ndata>20 & !any(dTs-median(dTs)>max(365,tspan/2))){
                    if(sum(dTs[1:2])>tspan/2){
                        Nbin <- 2
                    }else{
                        Nbin <- Nbin0
                    }
####check the data property and whether it is suitable for MP
                    source('MPcal.R')
                }else{
                    par(mfrow=c(2,2))
                    gamma <- NA
                }
                gammas <- c(gammas,gamma)
            }

###show the parameters
            source('show_table.R')
####finish plotting
        }
    }
    time.end <- proc.time()
    duration <- (time.end-time.start)[3]/3600#h
    cat('total time consumed:',round(duration,2),'h\n')
}else{
    Nsig <- 0
    plot.new()
}
dev.off()

##final output
cat('output Robj:\n',fobj,'\n')
cat('output pdf:\n',fpdf,'\n')
#save(solutions,out,Popt.all,res.all,opt.par,Psig,Ndet,Nact,durs,file=fobj)
save(list = ls(all=T),file=fobj)


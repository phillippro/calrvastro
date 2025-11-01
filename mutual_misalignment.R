source('../mcmc_func.R')
ds  <- list.dirs('.')
ds <- ds[ds!='.']
ds  <- gsub('\\./','',ds)
stars <- c()
pars <- t(rep(NA,81))
i3 <- which(ds=='BD-114672')
#for(k in 1:length(ds)){
for(k in 1:10){
#for(k in i3){
    d <- ds[k]
    fs <- list.files(d,pattern='hg123.+Robj',full.name=TRUE)
    if(length(fs)>0){
    lls <- gsub('.+lnlmax|\\.Robj','',fs)
    ind <- which.max(as.numeric(lls))
    per <- gsub('.+_P|_acc.+','',fs[ind])
    if(grepl('d',per)){
        tmp <- gsub('\\/.+','',fs[ind])
	cat(tmp,';',k,'/',length(ds),'\n')
        load(fs[ind],env=e0 <- new.env())
        par.stat <- e0$out$par.stat[[paste0('sig',e0$Nsig)]]
	cn <- colnames(par.stat)
	ii <- grep('Inc',cn)
	jj <- grep('per',cn)
	kk <- grep('Omega',cn)
	ll <- grep('^e',cn)
	nn <- grep('K',cn)
	cat(length(ii),'planets\n')
	for(i in 1:length(ii)){
	inc.opt <- par.stat['med',ii[i]]
	inc.minus <- par.stat['xminus.1sig',ii[i]]
	inc.plus <- par.stat['xplus.1sig',ii[i]]
	per.opt <- par.stat['med',jj[i]]
	per.minus <- par.stat['xminus.1sig',jj[i]]
	per.plus <- par.stat['xplus.1sig',jj[i]]
	Omega.opt <- par.stat['med',kk[i]]  
	Omega.minus <- par.stat['xminus.1sig',kk[i]]
	Omega.plus <- par.stat['xplus.1sig',kk[i]]
	e.opt <- par.stat['med',ll[i]]  
        e.minus <- par.stat['xminus.1sig',ll[i]]
	e.plus <- par.stat['xplus.1sig',ll[i]]
	K.opt <- par.stat['med',nn[i]]  
        K.minus <- par.stat['xminus.1sig',nn[i]]
	K.plus <- par.stat['xplus.1sig',nn[i]]
	Ks <- e0$out$mcmc.opt[[paste0('sig',e0$Nsig)]][,nn[i]]
	es <- e0$out$mcmc.opt[[paste0('sig',e0$Nsig)]][,ll[i]]
	incs <- e0$out$mcmc.opt[[paste0('sig',e0$Nsig)]][,ii[i]]
	ps <- exp(e0$out$mcmc.opt[[paste0('sig',e0$Nsig)]][,jj[i]])
	Ms <- rnorm(length(ps)*2,e0$out$Mstar,e0$out$eMstar)
	Ms <- Ms[Ms>0][1:length(ps)]
        msini <- K2msini.full(Ks,ps,es,Ms)
        mps <- msini$mj/sin(incs)#Mjup
	mm <- quantile(mps,c(0.16,0.5,0.84))
	mp.opt <- mm[2]
	mp.minus <- mm[1]
	mp.plus <- mm[3]
	inc1 <- e0$out$mcmc.opt[[paste0('sig',e0$Nsig)]][,ii[i]]
	Omega1 <- e0$out$mcmc.opt[[paste0('sig',e0$Nsig)]][,kk[i]]
	psis <- c()
	for(j in 1:length(ii)){
	    if(j!=i){
               inc2 <- e0$out$mcmc.opt[[paste0('sig',e0$Nsig)]][,ii[j]]
	       Omega2 <- e0$out$mcmc.opt[[paste0('sig',e0$Nsig)]][,kk[j]]
               psi <- acos(cos(inc1)*cos(inc2)+cos(Omega2-Omega1)*sin(inc1)*sin(inc2))*180/pi#mutual inclination
	       psis <- rbind(psis,quantile(psi,c(0.16,0.5,0.84)))
            }	    
	}
	psi.opt <- paste(psis[,2],collapse='d')
	psi.minus <- paste(psis[,1],collapse='d')
	psi.plus <- paste(psis[,3],collapse='d')
        tmp <- c(tmp,c(per.opt,per.minus,per.plus,inc.opt,inc.minus,inc.plus,Omega.opt,Omega.minus,Omega.plus,psi.opt,psi.minus,psi.plus,e.opt,e.minus,e.plus,K.opt,K.minus,K.plus,mp.opt,mp.minus,mp.plus))
        }
	if(length(tmp)<ncol(pars)) tmp  <- c(tmp,rep(NA,ncol(pars)-length(tmp)))
	pars <- data.frame(rbind(pars,tmp))
    }
    }   
}
colnames(pars) <- c('name',as.vector(outer(c('Per.opt','Per.lower','Per.upper','Inc.opt','Inc.lower','Inc.upper','Omega.opt','Omega.lower','Omega.upper','psi.opt','psi.lower','psi.upper','e.opt','e.minus','e.plus','K.opt','K.minus','K.plus','mp.opt','mp.minus','mp.pllus'),1:4,paste0)))
fout <- 'misaligned_starpar.txt'
cat(fout,'\n')
write.table(pars[-1,],file=fout,quote=FALSE,row.names=FALSE)

if(!dir.exists('/malibu')){
    dir1 <- 'results/'
}else{
    dir1 <- '/malibu/ffeng/astro_output/'
}
dir <- paste0(dir1,target,'/')
if(!file.exists(dir)){
    system(paste('mkdir',dir))
}
#if(!file.exists('results')) system('mkdir results')
#dir <- paste0('results/',target,'/')
#if(!file.exists(dir)) system(paste('mkdir',dir))
mcmc <- out$mcmc.opt[[paste0('sig',Nsig)]]
acceptance <- 100*(1-mean(duplicated(mcmc)))
Npar <- ncol(mcmc)-2
lnlmax <- max(mcmc[,ncol(mcmc)])
lnpmax <- max(mcmc[,ncol(mcmc)-1])
#par.opt <- mcmc[which.max(mcmc[,ncol(mcmc)-1]),1:Npar]
if(!any(names(out)=='relativity')) out$relativity <- FALSE
if(!exists('calibrate')) calibrate <- FALSE
if(!exists('rvc.type')) rvc.type <- FALSE
if(!exists('astro5')) astro5 <- TRUE
fname <-  paste0(dir,target,'_cal',calibrate,'_coplanar',coplanar,'_resonance',resonance,'_sta',stability,'_relativity',out$relativity,'_Niter',Niter,'_jitter',jitter,'_240715_Nset',Nset,'_Nrv',out$Nrv,'_',astro.type,'_Nsig',Nsig,'_P',paste0(round(Popt),collapse='d'),'_Esd',Esd,'_astro5',astro5,'_acc',format(acceptance,digit=2),'_rvc',rvc.type,'_priorf',priorf,'_lnpmax',format(lnpmax,digit=2),'_lnlmax',format(lnlmax,digit=2))
#save(list=ls(all=TRUE),file=paste0(fname,'.Robj'))
fobj <- paste0(fname,'.Robj')
cat(fobj,'\n')
if(!exists('save.memory')) save.memory <- TRUE
if(!exists('Nchain')) Nchain <- 1
if(!exists('Tmin')) Tmin <- 1000
#save(list=c('out','Nset','target','Nmax','Niter0','Ncores','ofac','fmin','fmax','Nsamp','Esd','fname','trv.all','RV.all','eRV.all','Nsig','Pmin','Pmax','tmin','tmax','ind.transit','prior.type','period.par','time.unit','Popt','par.opt','sig.type','noise.types','Nkeppar','basis','Rhat','conv','ind.chain','chain.type','Niter','tspan','Nchain','save.memory','Nchain','bases','basis','basis2','offset','astrometry','Nastro','cov.astro','Nmin','Nmax','nepoch','astro.type','photovary','tsim','rel.type','Tmin','eta0','coplanar','resonance','bc','fpar','rvc.type','astro5','priorf','use.gdr1','gdr1.epoch','gdr2.epoch'),file=fobj)
save(list=ls(all=TRUE),file=fobj)
